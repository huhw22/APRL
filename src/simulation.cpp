#include "simulation.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "particle_initializer.h"
#include "runtime_control.h"
#include "runtime_util.h"

namespace fel
{
  namespace
  {
    const int TRANSFER_UPPER_COUNT = 710;
    const int TRANSFER_UPPER_DATA = 711;
    const int TRANSFER_LOWER_COUNT = 712;
    const int TRANSFER_LOWER_DATA = 713;

    struct ParticlePacket
    {
      Double position[3];
      Double properVelocity[3];
      Double charge;
      Double mass;
      Double weight;
      std::uint64_t id;
      std::uint64_t sourceId;
    };

    struct TransferPacket
    {
      ParticlePacket particle;
      Double segmentStart[3];
    };

    ParticlePacket packParticle(const RelativisticParticleSI& particle)
    {
      ParticlePacket packet = {};
      for (unsigned int component = 0; component < 3; ++component)
        {
          packet.position[component] = particle.position[component];
          packet.properVelocity[component] =
            particle.properVelocity[component];
        }
      packet.charge = particle.charge;
      packet.mass = particle.mass;
      packet.weight = particle.weight;
      packet.id = particle.id;
      packet.sourceId = particle.sourceId;
      return packet;
    }

    RelativisticParticleSI unpackParticle(const ParticlePacket& packet)
    {
      RelativisticParticleSI particle;
      for (unsigned int component = 0; component < 3; ++component)
        {
          particle.position[component] = packet.position[component];
          particle.properVelocity[component] =
            packet.properVelocity[component];
        }
      particle.charge = packet.charge;
      particle.mass = packet.mass;
      particle.weight = packet.weight;
      particle.id = packet.id;
      particle.sourceId = packet.sourceId;
      return particle;
    }

    std::size_t exactCells(Double length, Double resolution,
                           const char* axis)
    {
      if (!(length > 0.0) || !(resolution > 0.0) ||
          !std::isfinite(length) || !std::isfinite(resolution))
        throw std::invalid_argument("E/B mesh lengths and resolutions must be positive");
      const Double requested = length / resolution;
      const std::size_t cells = static_cast<std::size_t>(
        std::llround(requested));
      const Double error = std::abs(static_cast<Double>(cells) - requested);
      if (cells < 3 || error > 1.0e-9 * std::max(1.0, requested))
        throw std::invalid_argument(std::string("E/B mesh on axis ") + axis +
          " must contain an integer number of at least three cells");
      return cells;
    }

    int checkedBytes(std::size_t records, std::size_t recordBytes)
    {
      if (records > static_cast<std::size_t>(INT_MAX) / recordBytes)
        throw std::overflow_error("MPI particle message exceeds int byte count");
      return static_cast<int>(records * recordBytes);
    }

    void exchangeTransfers(const std::vector<TransferPacket>& send,
                           int destination, int source,
                           int countTag, int dataTag,
                           MPI_Comm communicator,
                           std::vector<TransferPacket>& receive)
    {
      if (send.size() > static_cast<std::size_t>(INT_MAX))
        throw std::overflow_error("Too many particles cross one MPI interface");
      const int sendCount = static_cast<int>(send.size());
      int receiveCount = 0;
      MPI_Sendrecv(&sendCount, 1, MPI_INT, destination, countTag,
                   &receiveCount, 1, MPI_INT, source, countTag,
                   communicator, MPI_STATUS_IGNORE);
      if (receiveCount < 0)
        throw std::runtime_error("Negative MPI particle transfer count");
      receive.resize(static_cast<std::size_t>(receiveCount));
      const int sendBytes = checkedBytes(send.size(), sizeof(TransferPacket));
      const int receiveBytes = checkedBytes(receive.size(),
                                            sizeof(TransferPacket));
      MPI_Sendrecv(send.empty() ? NULL : &send[0], sendBytes, MPI_BYTE,
                   destination, dataTag,
                   receive.empty() ? NULL : &receive[0], receiveBytes,
                   MPI_BYTE, source, dataTag, communicator,
                   MPI_STATUS_IGNORE);
    }

    bool clipSegmentToBox(const FieldVector<Double>& start,
                          const FieldVector<Double>& end,
                          const Double lower[3], const Double upper[3],
                          FieldVector<Double>& clipped,
                          Double& exitFraction)
    {
      bool outside = false;
      exitFraction = 1.0;
      for (unsigned int axis = 0; axis < 3; ++axis)
        {
          if (end[axis] < lower[axis])
            {
              outside = true;
              exitFraction = std::min(exitFraction,
                (lower[axis] - start[axis]) /
                (end[axis] - start[axis]));
            }
          else if (end[axis] > upper[axis])
            {
              outside = true;
              exitFraction = std::min(exitFraction,
                (upper[axis] - start[axis]) /
                (end[axis] - start[axis]));
            }
        }
      if (!outside) return false;
      if (!(exitFraction >= 0.0) || !(exitFraction <= 1.0) ||
          !std::isfinite(exitFraction))
        throw std::runtime_error("Cannot clip escaped particle trajectory");
      clipped = start;
      for (unsigned int axis = 0; axis < 3; ++axis)
        {
          clipped[axis] += exitFraction * (end[axis] - start[axis]);
          clipped[axis] = std::max(lower[axis],
                                   std::min(upper[axis], clipped[axis]));
        }
      return true;
    }
  }

  Simulation::Simulation(const SimulationConfig& config,
                         MPI_Comm communicator)
    : config_(config), communicator_(communicator), rank_(0), size_(1),
      globalGeometry_(), localGeometry_(), localZOffset_(0),
      globalOriginBox_(0.0), localOriginBox_(0.0), frame_(), sources_(),
      fields_(), halo_(), incident_(), particles_(), detectors_(),
      trajectoryWriter_(),
      trajectoryRhythmSI_(0.0), nextTrajectorySampleTime_(0.0),
      trajectorySamplesSinceFlush_(0), timeBoxSI_(0.0),
      totalTimeBoxSI_(0.0), step_(0), interrupted_(false),
      configuredStopReached_(false), lostParticles_(0), stopReason_()
  {
    if (communicator_ == MPI_COMM_NULL)
      throw std::invalid_argument("E/B solver communicator cannot be null");
    MPI_Comm_rank(communicator_, &rank_);
    MPI_Comm_size(communicator_, &size_);
  }

  void Simulation::solve()
  {
    initialize();
    if (config_.trajectory.enabled)
      sampleTrajectory();
    if (detectors_)
      detectors_->sampleFieldPlanes(*fields_, timeBoxSI_);

    while (timeBoxSI_ < totalTimeBoxSI_)
      {
        if ((step_ % config_.runtime.stopCheckIntervalSteps == 0) &&
            synchronizedStopRequested())
          break;

        fields_->clearCurrent();
        fields_->advanceMagnetic();
        if (incident_)
          incident_->correctAfterMagneticUpdate(
            *fields_, sources_, localOriginBox_, timeBoxSI_, frame_);

        pushDepositAndMigrate();
        halo_->advanceElectric(*fields_);
        if (incident_)
          incident_->correctAfterElectricUpdate(
            *fields_, sources_, localOriginBox_, timeBoxSI_, frame_);

        timeBoxSI_ += globalGeometry_.dt;
        ++step_;
        if (config_.trajectory.enabled &&
            timeBoxSI_ >= nextTrajectorySampleTime_)
          sampleTrajectory();
        if (detectors_)
          detectors_->sampleFieldPlanes(*fields_, timeBoxSI_);
        if (configuredStopReached())
          {
            configuredStopReached_ = true;
            break;
          }
      }

    if (!interrupted_) synchronizedStopRequested();
    const bool completed = !interrupted_ && configuredStopReached_;
    finalizeTrajectoryOutput(completed);
    if (detectors_) detectors_->close(completed);
    if (configuredStopReached_ && rank_ == 0)
      logRoot(communicator_, "Configured stop reached: " + stopReason_);
    if (interrupted_ && rank_ == 0)
      logRoot(communicator_,
        "Direct E/B run stopped after a complete field step; committed trajectory and detector outputs remain readable.");
    if (!interrupted_ && !configuredStopReached_ &&
        timeBoxSI_ >= totalTimeBoxSI_)
      throw std::runtime_error(
        "mesh.duration was exhausted before the configured physical stop; trajectory and detector outputs are readable but marked incomplete");
  }

  void Simulation::initialize()
  {
    if (!(config_.mesh.boostGamma >= 1.0) ||
        !std::isfinite(config_.mesh.boostGamma))
      throw std::invalid_argument(
        "Direct E/B solver requires boost_gamma >= 1");

    initializeGeometry();
    initializeParticles();
    initializeSources();

    fields_.reset(new EBFieldGrid(localGeometry_));
    halo_.reset(new EBZSlabHaloExchange(communicator_));
    halo_->installPhysicalBoundaryMask(*fields_);

    if (sources_.maxwellIncidentWaveCount() > 0)
      {
        const EBTFSFRegion region(
          1, globalGeometry_.nx - 1,
          1, globalGeometry_.ny - 1,
          1, globalGeometry_.nz - 1);
        incident_.reset(new EBMaxwellIncidentInjector(
          localGeometry_, region, localZOffset_, globalGeometry_.nz));
        incident_->initialize(*fields_, sources_, localOriginBox_,
                              timeBoxSI_, frame_);
      }

    initializeTrajectoryOutput();
    initializeDetectorOutput();
    if (rank_ == 0)
      {
        logRoot(communicator_, config_.runtime.interactive() ?
          "Runtime strategy: interactive local-test mode; SIGINT/SIGTERM stops at a complete field step." :
          "Runtime strategy: HPC throughput mode; no signal polling or periodic durability flushes.");
        logRoot(communicator_,
          "Direct SI E/B simulation active; no A/phi state is allocated.");
        logRoot(communicator_,
          "WARNING: the initial Gauss-consistent particle field and CPML are not implemented; this is not yet a final radiation-production solver.");
        logRoot(communicator_,
          "WARNING: particle subcycling is not implemented; the E/B field step must resolve every prescribed device field.");
      }
  }

  void Simulation::initializeGeometry()
  {
    const Double dx = config_.mesh.resolution[0];
    const Double dy = config_.mesh.resolution[1];
    const Double dz = config_.mesh.resolution[2];
    const std::size_t nx = exactCells(config_.mesh.lengths[0], dx, "x");
    const std::size_t ny = exactCells(config_.mesh.lengths[1], dy, "y");
    const std::size_t nz = exactCells(config_.mesh.lengths[2], dz, "z");
    if (nz < static_cast<std::size_t>(size_))
      throw std::invalid_argument(
        "Direct E/B grid requires at least one z cell per MPI rank");

    const Double inverseSpacingSquared =
      1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz);
    const Double dt = 0.95 /
      (SI::c * std::sqrt(inverseSpacingSquared));
    globalGeometry_ = EBGridGeometry(nx, ny, nz, dx, dy, dz, dt);
    const Slab local = slabForRank(rank_);
    localZOffset_ = local.offset;
    localGeometry_ = EBGridGeometry(nx, ny, local.cells,
                                    dx, dy, dz, dt);

    for (unsigned int axis = 0; axis < 3; ++axis)
      globalOriginBox_[axis] =
        config_.mesh.center[axis] -
        0.5 * static_cast<Double>(axis == 0 ? nx : (axis == 1 ? ny : nz)) *
        (axis == 0 ? dx : (axis == 1 ? dy : dz));
    localOriginBox_ = globalOriginBox_;
    localOriginBox_[2] += static_cast<Double>(localZOffset_) * dz;
    totalTimeBoxSI_ = config_.mesh.duration;
    if (!(totalTimeBoxSI_ > 0.0) || !std::isfinite(totalTimeBoxSI_))
      throw std::invalid_argument(
        "Direct E/B total-time must be a positive box-frame duration");
  }

  void Simulation::initializeParticles()
  {
    const Double firstPhysicalEntrance =
      firstBeamlinePhysicalEntranceLab();
    const Double firstInteractionEntrance =
      firstBeamlineInteractionEntranceLab();
    const Double boxReferenceZ = config_.mesh.center[2];
    frame_.setOriginsFromGamma(config_.mesh.boostGamma, SI::c, 0.0,
      config_.reference.initialCenterZ,
      boxReferenceZ);

    particles_ = ParticleInitializer::create(config_, communicator_);
    const unsigned long long localCount =
      static_cast<unsigned long long>(particles_.size());
    unsigned long long idOffset = 0;
    MPI_Exscan(&localCount, &idOffset, 1, MPI_UNSIGNED_LONG_LONG,
               MPI_SUM, communicator_);
    if (rank_ == 0) idOffset = 0;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      particles_[index].id = idOffset + index + 1;

    SIBunchPlacement placement;
    placement.firstInteractionEntranceLab = firstInteractionEntrance;
    placement.referencePositionLab = config_.reference.initialCenterZ;
    /* A box-frame z cell spans gamma*dz at fixed box time in the lab.  Keep
     * that as the minimum recommended clearance beyond the relative head. */
    placement.recommendationMarginLab =
      config_.mesh.boostGamma * globalGeometry_.dz;

    const SIBunchPlacementReport placementReport =
      SIBunchPreprocessor::placeRelativeLabSnapshot(
        particles_, placement, communicator_);
    const SIBunchBoostReport localBoostReport =
      SIBunchPreprocessor::boostLabSnapshotToBoxTimeZero(
        particles_, frame_);
    SIBunchBoostReport boostReport;
    const Double localEarliest = particles_.empty() ?
      std::numeric_limits<Double>::infinity() :
      localBoostReport.earliestLabEventTime;
    const Double localLatest = particles_.empty() ?
      -std::numeric_limits<Double>::infinity() :
      localBoostReport.latestLabEventTime;
    MPI_Allreduce(&localEarliest, &boostReport.earliestLabEventTime,
                  1, MPI_DOUBLE, MPI_MIN, communicator_);
    MPI_Allreduce(&localLatest, &boostReport.latestLabEventTime,
                  1, MPI_DOUBLE, MPI_MAX, communicator_);
    MPI_Allreduce(&localBoostReport.maximumAbsoluteDriftTime,
                  &boostReport.maximumAbsoluteDriftTime,
                  1, MPI_DOUBLE, MPI_MAX, communicator_);

    Double localEventHead = -std::numeric_limits<Double>::infinity();
    for (std::size_t index = 0; index < particles_.size(); ++index)
      {
        Double eventTimeLab = 0.0;
        Double eventZLab = 0.0;
        frame_.boxToLab(0.0, particles_[index].position[2],
                        eventTimeLab, eventZLab);
        localEventHead = std::max(localEventHead, eventZLab);
      }
    Double eventHead = 0.0;
    MPI_Allreduce(&localEventHead, &eventHead, 1, MPI_DOUBLE, MPI_MAX,
                  communicator_);
    const Double entranceTolerance = 64.0 *
      std::numeric_limits<Double>::epsilon() *
      std::max(1.0, std::max(std::abs(firstInteractionEntrance),
                             std::abs(eventHead)));
    if (eventHead >= firstInteractionEntrance - entranceTolerance)
      {
        const Double transformedHeadOffset = eventHead -
          config_.reference.initialCenterZ;
        const Double recommendedCenter = firstInteractionEntrance -
          transformedHeadOffset - placement.recommendationMarginLab;
        std::ostringstream message;
        message << "Initial Lorentz transform places the bunch front at or "
          "inside the first element interaction region: interaction_start_z="
          << firstInteractionEntrance
          << " m, transformed_front_z=" << eventHead
          << " m. Set beam.reference.initial_center_z to at most "
          << recommendedCenter << " m (" << recommendedCenter /
               config_.inputUnits.length
          << " in the configured length unit). The physical first-element "
          "entrance remains z=" << firstPhysicalEntrance << " m.";
        throw std::runtime_error(message.str());
      }
    redistributeParticles();
    validateParticlesInsideGlobalBox();

    if (rank_ == 0)
      {
        std::ostringstream placementMessage;
        placementMessage << "E/B bunch placement: particles="
          << placementReport.particles << ", reference z [m]="
          << placementReport.referencePositionLab
          << ", relative head z [m]=" << placementReport.relativeHeadLab
          << ", physical first entrance [m]=" << firstPhysicalEntrance
          << ", interaction start [m]=" << firstInteractionEntrance
          << ", transformed interaction gap [m]="
          << firstInteractionEntrance - eventHead;
        logRoot(communicator_, placementMessage.str());
        std::ostringstream boostMessage;
        boostMessage << "Free-drift Lorentz events [s]: "
          << boostReport.earliestLabEventTime << " to "
          << boostReport.latestLabEventTime;
        logRoot(communicator_, boostMessage.str());
      }
  }

  void Simulation::initializeSources()
  {
    for (std::size_t i = 0; i < config_.waves.size(); ++i)
      sources_.addWave(config_.waves[i].source, config_.waves[i].evolution);
    for (std::size_t i = 0; i < config_.magnets.size(); ++i)
      sources_.addPrescribedLabMagnet(config_.magnets[i]);
  }

  void Simulation::initializeTrajectoryOutput()
  {
    if (!config_.trajectory.enabled) return;
    trajectoryRhythmSI_ = config_.trajectory.rhythm;
    if (!(trajectoryRhythmSI_ > 0.0))
      throw std::invalid_argument(
        "Direct E/B trajectory rhythm must be positive");
    const std::string baseName = joinPath(
      config_.trajectory.directory,
      config_.trajectory.basename.empty() ?
        "trajectory" : config_.trajectory.basename);
    createDirectories(config_.trajectory.directory, communicator_);
    std::ostringstream filename;
    filename << baseName << "-rank-" << std::setfill('0')
             << std::setw(5) << rank_ << ".h5";
    trajectoryWriter_.open(filename.str(), rank_, size_,
      config_.trajectory.bufferRecords, config_.trajectory.compression,
      config_.runtime.interactive());
    nextTrajectorySampleTime_ = 0.0;
    trajectorySamplesSinceFlush_ = 0;
  }

  void Simulation::initializeDetectorOutput()
  {
    if (!config_.detectors.enabled()) return;
    detectors_.reset(new LabDetectorManager(
      config_.detectors, globalGeometry_,
      globalOriginBox_, localOriginBox_, frame_, communicator_));
    if (rank_ == 0)
      {
        std::ostringstream message;
        message << "Laboratory detector planes active: fields="
                << config_.detectors.fieldPlanes.size()
                << ", particles="
                << config_.detectors.particlePlanes.size()
                << "; detector HDF5 is written only by MPI rank zero.";
        logRoot(communicator_, message.str());
      }
  }

  void Simulation::sampleTrajectory()
  {
    if (!trajectoryWriter_.isOpen()) return;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      {
        const RelativisticParticleSI& particle = particles_[index];
        TrajectoryRecord record;
        record.particleId = particle.id;
        record.sourceId = particle.sourceId;
        frame_.boxToLab(timeBoxSI_, particle.position[2],
                        record.time, record.position[2]);
        record.position[0] = particle.position[0];
        record.position[1] = particle.position[1];
        const Double gammaBox =
          BoostFrameTransform::gammaFromProperVelocity(
            particle.properVelocity);
        record.properVelocity[0] = particle.properVelocity[0];
        record.properVelocity[1] = particle.properVelocity[1];
        record.properVelocity[2] = frame_.gamma() *
          (particle.properVelocity[2] + frame_.beta() * gammaBox);
        record.charge = particle.charge;
        record.weight = particle.weight;
        trajectoryWriter_.append(record);
      }
    if (config_.runtime.interactive())
      {
        ++trajectorySamplesSinceFlush_;
        if (trajectorySamplesSinceFlush_ >=
            config_.trajectory.flushEverySamples)
          {
            trajectoryWriter_.flush();
            trajectorySamplesSinceFlush_ = 0;
          }
      }
    const Double tolerance = 32.0 *
      std::numeric_limits<Double>::epsilon() *
      std::max(1.0, std::abs(timeBoxSI_));
    do
      nextTrajectorySampleTime_ += trajectoryRhythmSI_;
    while (nextTrajectorySampleTime_ <= timeBoxSI_ + tolerance);
  }

  void Simulation::finalizeTrajectoryOutput(bool completed)
  {
    if (trajectoryWriter_.isOpen()) trajectoryWriter_.close(completed);
  }

  Simulation::Slab Simulation::slabForRank(int rank) const
  {
    const std::size_t base = globalGeometry_.nz /
      static_cast<std::size_t>(size_);
    const std::size_t remainder = globalGeometry_.nz %
      static_cast<std::size_t>(size_);
    Slab slab;
    slab.cells = base +
      (static_cast<std::size_t>(rank) < remainder ? 1 : 0);
    slab.offset = static_cast<std::size_t>(rank) * base +
      std::min(static_cast<std::size_t>(rank), remainder);
    return slab;
  }

  int Simulation::ownerRank(Double boxZ) const
  {
    const Double coordinate =
      (boxZ - globalOriginBox_[2]) / globalGeometry_.dz;
    const Double tolerance = 64.0 *
      std::numeric_limits<Double>::epsilon() *
      std::max(1.0, std::abs(coordinate));
    if (coordinate < -tolerance ||
        coordinate > static_cast<Double>(globalGeometry_.nz) + tolerance)
      return -1;
    std::size_t cell = coordinate >= static_cast<Double>(globalGeometry_.nz) ?
      globalGeometry_.nz - 1 :
      static_cast<std::size_t>(std::max(0.0, std::floor(coordinate)));
    const std::size_t base = globalGeometry_.nz /
      static_cast<std::size_t>(size_);
    const std::size_t remainder = globalGeometry_.nz %
      static_cast<std::size_t>(size_);
    const std::size_t longCells = (base + 1) * remainder;
    if (cell < longCells)
      return static_cast<int>(cell / (base + 1));
    return static_cast<int>(remainder + (cell - longCells) / base);
  }

  void Simulation::redistributeParticles()
  {
    std::vector<int> sendCounts(static_cast<std::size_t>(size_), 0);
    std::vector<int> destinations(particles_.size(), -1);
    int localInvalid = 0;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      {
        destinations[index] = ownerRank(particles_[index].position[2]);
        if (destinations[index] < 0)
          localInvalid = 1;
        else
          ++sendCounts[static_cast<std::size_t>(destinations[index])];
      }
    int globalInvalid = 0;
    MPI_Allreduce(&localInvalid, &globalInvalid, 1, MPI_INT, MPI_MAX,
                  communicator_);
    if (globalInvalid)
      throw std::out_of_range(
        "Boosted bunch lies outside the direct E/B longitudinal box");

    std::vector<int> receiveCounts(static_cast<std::size_t>(size_), 0);
    MPI_Alltoall(&sendCounts[0], 1, MPI_INT,
                 &receiveCounts[0], 1, MPI_INT, communicator_);
    std::vector<int> sendDisplacements(static_cast<std::size_t>(size_), 0);
    std::vector<int> receiveDisplacements(static_cast<std::size_t>(size_), 0);
    for (int process = 1; process < size_; ++process)
      {
        sendDisplacements[process] = sendDisplacements[process - 1] +
          sendCounts[process - 1];
        receiveDisplacements[process] = receiveDisplacements[process - 1] +
          receiveCounts[process - 1];
      }
    const int totalReceive = receiveDisplacements[size_ - 1] +
      receiveCounts[size_ - 1];
    std::vector<ParticlePacket> sendPackets(particles_.size());
    std::vector<int> cursor(sendDisplacements);
    for (std::size_t index = 0; index < particles_.size(); ++index)
      sendPackets[static_cast<std::size_t>(
        cursor[static_cast<std::size_t>(destinations[index])]++)] =
        packParticle(particles_[index]);
    std::vector<ParticlePacket> receivePackets(
      static_cast<std::size_t>(totalReceive));

    std::vector<int> sendBytes(sendCounts);
    std::vector<int> receiveBytes(receiveCounts);
    std::vector<int> sendByteDisplacements(sendDisplacements);
    std::vector<int> receiveByteDisplacements(receiveDisplacements);
    for (int process = 0; process < size_; ++process)
      {
        sendBytes[process] = checkedBytes(
          static_cast<std::size_t>(sendCounts[process]),
          sizeof(ParticlePacket));
        receiveBytes[process] = checkedBytes(
          static_cast<std::size_t>(receiveCounts[process]),
          sizeof(ParticlePacket));
        sendByteDisplacements[process] = checkedBytes(
          static_cast<std::size_t>(sendDisplacements[process]),
          sizeof(ParticlePacket));
        receiveByteDisplacements[process] = checkedBytes(
          static_cast<std::size_t>(receiveDisplacements[process]),
          sizeof(ParticlePacket));
      }
    MPI_Alltoallv(sendPackets.empty() ? NULL : &sendPackets[0],
                  &sendBytes[0], &sendByteDisplacements[0], MPI_BYTE,
                  receivePackets.empty() ? NULL : &receivePackets[0],
                  &receiveBytes[0], &receiveByteDisplacements[0], MPI_BYTE,
                  communicator_);
    particles_.clear();
    particles_.reserve(receivePackets.size());
    for (std::size_t index = 0; index < receivePackets.size(); ++index)
      particles_.push_back(unpackParticle(receivePackets[index]));
  }

  void Simulation::pushDepositAndMigrate()
  {
    ChargeConservingCurrentDepositor depositor(*fields_, localOriginBox_);
    std::vector<RelativisticParticleSI> retained;
    retained.reserve(particles_.size());
    std::vector<TransferPacket> sendLower;
    std::vector<TransferPacket> sendUpper;
    const Double lowerZ = localOriginBox_[2];
    const Double upperZ = lowerZ +
      static_cast<Double>(localGeometry_.nz) * localGeometry_.dz;
    const Double globalLower[3] = {
      globalOriginBox_[0], globalOriginBox_[1], globalOriginBox_[2]
    };
    const Double globalUpper[3] = {
      globalOriginBox_[0] + static_cast<Double>(globalGeometry_.nx) * globalGeometry_.dx,
      globalOriginBox_[1] + static_cast<Double>(globalGeometry_.ny) * globalGeometry_.dy,
      globalOriginBox_[2] + static_cast<Double>(globalGeometry_.nz) * globalGeometry_.dz
    };
    unsigned long long localLost = 0;

    for (std::size_t index = 0; index < particles_.size(); ++index)
      {
        RelativisticParticleSI particle = particles_[index];
        const FieldVector<Double> start(particle.position);
        RelativisticBorisPusher::pushFromGridAndPrescribedLab(
          particle, *fields_, localOriginBox_, sources_, frame_,
          timeBoxSI_, globalGeometry_.dt);
        if (detectors_)
          detectors_->captureParticleStep(
            particles_[index], particle, timeBoxSI_,
            timeBoxSI_ + globalGeometry_.dt);

        const bool crossLower = rank_ > 0 &&
          particle.position[2] < lowerZ;
        const bool crossUpper = rank_ + 1 < size_ &&
          particle.position[2] >= upperZ;
        const bool crossesInterface = crossLower || crossUpper;
        Double interfaceFraction = 1.0;
        if (crossesInterface)
          {
            const Double boundary = crossLower ? lowerZ : upperZ;
            interfaceFraction = (boundary - start[2]) /
              (particle.position[2] - start[2]);
            if (!(interfaceFraction >= 0.0 && interfaceFraction <= 1.0))
              throw std::runtime_error(
                "Invalid MPI-interface particle crossing");
          }

        FieldVector<Double> clipped(0.0);
        Double exitFraction = 1.0;
        const bool particleOutside = clipSegmentToBox(
          start, particle.position, globalLower, globalUpper,
          clipped, exitFraction);
        const Double crossingTolerance = 64.0 *
          std::numeric_limits<Double>::epsilon();
        if (particleOutside &&
            (!crossesInterface ||
             exitFraction <= interfaceFraction + crossingTolerance))
          {
            depositor.depositSegment(start, clipped, particle.charge);
            ++localLost;
            continue;
          }

        if (crossLower || crossUpper)
          {
            FieldVector<Double> interfacePosition(start);
            for (unsigned int axis = 0; axis < 3; ++axis)
              interfacePosition[axis] += interfaceFraction *
                (particle.position[axis] - start[axis]);
            depositor.depositSegment(start, interfacePosition,
                                     particle.charge);
            TransferPacket transfer = {};
            transfer.particle = packParticle(particle);
            for (unsigned int axis = 0; axis < 3; ++axis)
              transfer.segmentStart[axis] = interfacePosition[axis];
            (crossLower ? sendLower : sendUpper).push_back(transfer);
          }
        else
          {
            depositor.depositSegment(start, particle.position,
                                     particle.charge);
            retained.push_back(particle);
          }
      }

    std::vector<TransferPacket> receiveLower;
    std::vector<TransferPacket> receiveUpper;
    exchangeTransfers(sendUpper, halo_->upperRank(), halo_->lowerRank(),
      TRANSFER_UPPER_COUNT, TRANSFER_UPPER_DATA,
      communicator_, receiveLower);
    exchangeTransfers(sendLower, halo_->lowerRank(), halo_->upperRank(),
      TRANSFER_LOWER_COUNT, TRANSFER_LOWER_DATA,
      communicator_, receiveUpper);

    const std::vector<TransferPacket>* received[2] = {
      &receiveLower, &receiveUpper
    };
    for (unsigned int side = 0; side < 2; ++side)
      for (std::size_t index = 0; index < received[side]->size(); ++index)
        {
          const TransferPacket& transfer = (*received[side])[index];
          const RelativisticParticleSI particle =
            unpackParticle(transfer.particle);
          FieldVector<Double> segmentStart(0.0);
          for (unsigned int axis = 0; axis < 3; ++axis)
            segmentStart[axis] = transfer.segmentStart[axis];
          FieldVector<Double> clipped(0.0);
          Double exitFraction = 1.0;
          if (clipSegmentToBox(segmentStart, particle.position,
                               globalLower, globalUpper,
                               clipped, exitFraction))
            {
              depositor.depositSegment(segmentStart, clipped,
                                       particle.charge);
              ++localLost;
            }
          else
            {
              depositor.depositSegment(segmentStart, particle.position,
                                       particle.charge);
              retained.push_back(particle);
            }
        }
    particles_.swap(retained);
    unsigned long long globalLost = 0;
    MPI_Allreduce(&localLost, &globalLost, 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, communicator_);
    lostParticles_ += globalLost;
    if (detectors_) detectors_->collectParticleCrossings();
  }

  bool Simulation::synchronizedStopRequested()
  {
    if (!config_.runtime.interactive()) return false;
    int local = RuntimeControl::stopRequested() ? 1 : 0;
    int global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, communicator_);
    interrupted_ = interrupted_ || global != 0;
    return interrupted_;
  }

  bool Simulation::configuredStopReached()
  {
    const unsigned long long localCount =
      static_cast<unsigned long long>(particles_.size());
    unsigned long long globalCount = 0;
    MPI_Allreduce(&localCount, &globalCount, 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, communicator_);

    if (config_.stop.mode == StopMode::ReferenceCenterZ)
      {
        const Double referenceZ = frame_.labZFromBoxZT(
          config_.mesh.center[2], timeBoxSI_);
        const Double tolerance = 64.0 *
          std::numeric_limits<Double>::epsilon() *
          std::max(1.0, std::max(std::abs(referenceZ),
                                 std::abs(config_.stop.referenceZ)));
        if (referenceZ < config_.stop.referenceZ - tolerance) return false;
        std::ostringstream message;
        message << "boost reference centre reached lab z=" << referenceZ
                << " m (target " << config_.stop.referenceZ << " m)";
        if (lostParticles_ > 0)
          message << "; excluded escaped particles=" << lostParticles_;
        stopReason_ = message.str();
        return true;
      }

    const Double interactionExit = lastBeamlineInteractionExitLab();
    int localNotPast = 0;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      {
        Double eventTimeLab = 0.0;
        Double eventZLab = 0.0;
        frame_.boxToLab(timeBoxSI_, particles_[index].position[2],
                        eventTimeLab, eventZLab);
        const Double tolerance = 64.0 *
          std::numeric_limits<Double>::epsilon() *
          std::max(1.0, std::max(std::abs(eventZLab),
                                 std::abs(interactionExit)));
        if (eventZLab <= interactionExit + tolerance)
          localNotPast = 1;
      }
    int globalNotPast = 0;
    MPI_Allreduce(&localNotPast, &globalNotPast, 1, MPI_INT, MPI_MAX,
                  communicator_);
    if (globalNotPast) return false;

    std::ostringstream message;
    if (globalCount == 0)
      message << "no valid particles remain before the last interaction "
              << "boundary at z=" << interactionExit << " m";
    else
      message << "all " << globalCount
              << " valid particles passed the last interaction boundary at z="
              << interactionExit << " m";
    if (lostParticles_ > 0)
      message << "; excluded escaped particles=" << lostParticles_;
    stopReason_ = message.str();
    return true;
  }

  Double Simulation::firstBeamlinePhysicalEntranceLab() const
  {
    Double entrance = std::numeric_limits<Double>::infinity();
    for (std::size_t element = 0;
         element < config_.beamlineElements.size(); ++element)
      entrance = std::min(entrance,
        config_.beamlineElements[element].physicalEntrance);
    return entrance;
  }

  Double Simulation::firstBeamlineInteractionEntranceLab() const
  {
    Double entrance = std::numeric_limits<Double>::infinity();
    for (std::size_t element = 0;
         element < config_.beamlineElements.size(); ++element)
      entrance = std::min(entrance,
        config_.beamlineElements[element].interactionEntrance);
    return entrance;
  }

  Double Simulation::lastBeamlineInteractionExitLab() const
  {
    Double exit = -std::numeric_limits<Double>::infinity();
    for (std::size_t element = 0;
         element < config_.beamlineElements.size(); ++element)
      exit = std::max(exit,
        config_.beamlineElements[element].interactionExit);
    return exit;
  }

  void Simulation::validateParticlesInsideGlobalBox() const
  {
    const Double upper[3] = {
      globalOriginBox_[0] + static_cast<Double>(globalGeometry_.nx) * globalGeometry_.dx,
      globalOriginBox_[1] + static_cast<Double>(globalGeometry_.ny) * globalGeometry_.dy,
      globalOriginBox_[2] + static_cast<Double>(globalGeometry_.nz) * globalGeometry_.dz
    };
    int localInvalid = 0;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      for (unsigned int axis = 0; axis < 3; ++axis)
        if (particles_[index].position[axis] < globalOriginBox_[axis] ||
            particles_[index].position[axis] > upper[axis])
          localInvalid = 1;
    int globalInvalid = 0;
    MPI_Allreduce(&localInvalid, &globalInvalid, 1, MPI_INT, MPI_MAX,
                  communicator_);
    if (globalInvalid)
      throw std::out_of_range(
        "Boosted bunch does not fit inside the configured direct E/B box");
  }
}
