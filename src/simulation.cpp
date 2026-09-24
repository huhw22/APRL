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
    const int CARRIER_UPPER_COUNT = 714;
    const int CARRIER_UPPER_DATA = 715;
    const int CARRIER_LOWER_COUNT = 716;
    const int CARRIER_LOWER_DATA = 717;
    const int RETIREMENT_UPPER_COUNT = 718;
    const int RETIREMENT_UPPER_DATA = 719;
    const int RETIREMENT_LOWER_COUNT = 720;
    const int RETIREMENT_LOWER_DATA = 721;
    const Double PLANCK_EV_SECOND = 4.135667696e-15;

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
      Double segmentStartFraction;
    };

    struct CarrierTransferPacket
    {
      Double position[3];
      Double velocity[3];
      Double charge;
      Double currentWeight;
      Double segmentStart[3];
      Double segmentStartFraction;
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

    CarrierTransferPacket packCarrier(
        const ParticleCPMLCarrier& carrier,
        const FieldVector<Double>& segmentStart,
        Double segmentStartFraction)
    {
      CarrierTransferPacket packet = {};
      for (unsigned int component = 0; component < 3; ++component)
        {
          packet.position[component] = carrier.position[component];
          packet.velocity[component] = carrier.velocity[component];
          packet.segmentStart[component] = segmentStart[component];
        }
      packet.charge = carrier.charge;
      packet.currentWeight = carrier.currentWeight;
      packet.segmentStartFraction = segmentStartFraction;
      return packet;
    }

    ParticleCPMLCarrier unpackCarrier(
        const CarrierTransferPacket& packet)
    {
      ParticleCPMLCarrier carrier;
      for (unsigned int component = 0; component < 3; ++component)
        {
          carrier.position[component] = packet.position[component];
          carrier.velocity[component] = packet.velocity[component];
        }
      carrier.charge = packet.charge;
      carrier.currentWeight = packet.currentWeight;
      return carrier;
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

    std::string transverseGridRecommendation(
        const char* axis, Double length, Double dz, Double inputLengthUnit)
    {
      const Double cellsReal = length / dz;
      const Double tolerance = 64.0 *
        std::numeric_limits<Double>::epsilon() *
        std::max(1.0, std::abs(cellsReal));
      const std::size_t maximumCells = cellsReal >= 1.0 ?
        static_cast<std::size_t>(std::floor(cellsReal + tolerance)) : 0;
      std::ostringstream message;
      message << axis << " spacing must be >= dz=" << dz << " m";
      if (maximumCells >= 3)
        {
          const Double recommended = length /
            static_cast<Double>(maximumCells);
          message << "; for the configured " << axis
                  << " length use at most " << maximumCells
                  << " cells, for example resolution="
                  << recommended << " m ("
                  << recommended / inputLengthUnit
                  << " in the configured length unit)";
        }
      else
        message << "; this transverse length cannot contain three valid cells, so increase it to at least "
                << 3.0 * dz << " m";
      return message.str();
    }

    int checkedBytes(std::size_t records, std::size_t recordBytes)
    {
      if (records > static_cast<std::size_t>(INT_MAX) / recordBytes)
        throw std::overflow_error("MPI particle message exceeds int byte count");
      return static_cast<int>(records * recordBytes);
    }

    template<typename Packet>
    void exchangePackets(const std::vector<Packet>& send,
                         int destination, int source,
                         int countTag, int dataTag,
                         MPI_Comm communicator,
                         std::vector<Packet>& receive)
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
      const int sendBytes = checkedBytes(send.size(), sizeof(Packet));
      const int receiveBytes = checkedBytes(receive.size(),
                                            sizeof(Packet));
      MPI_Sendrecv(send.empty() ? NULL : &send[0], sendBytes, MPI_BYTE,
                   destination, dataTag,
                   receive.empty() ? NULL : &receive[0], receiveBytes,
                   MPI_BYTE, source, dataTag, communicator,
                   MPI_STATUS_IGNORE);
    }

  }

  Simulation::Simulation(const SimulationConfig& config,
                         MPI_Comm communicator)
    : config_(config), detectorConfig_(config.detectors),
      beamlineElements_(config.beamlineElements),
      communicator_(communicator), rank_(0), size_(1),
      globalGeometry_(), localGeometry_(), localZOffset_(0),
      globalOriginBox_(0.0), localOriginBox_(0.0), frame_(), sources_(),
      fields_(), halo_(), incident_(), particles_(), particleCPML_(),
      pmlCarriers_(), retirementCarriers_(), particleBoundary_(), detectors_(),
      trajectoryWriter_(),
      trajectoryRhythmSI_(0.0), nextTrajectorySampleTime_(0.0),
      trajectorySamplesSinceFlush_(0), timeBoxSI_(0.0),
      totalTimeBoxSI_(0.0), step_(0), interrupted_(false),
      configuredStopReached_(false), lostParticles_(0),
      retirementEntryCount_(0), retirementEntryCharge_(0.0),
      retirementExitCount_(0), retirementExitResidualCharge_(0.0),
      peakRetirementCarriers_(0), stopReason_()
  {
    if (communicator_ == MPI_COMM_NULL)
      throw std::invalid_argument("E/B solver communicator cannot be null");
    MPI_Comm_rank(communicator_, &rank_);
    MPI_Comm_size(communicator_, &size_);
    for (std::size_t face = 0; face < 6; ++face)
      {
        cpmlEntryCount_[face] = 0;
        cpmlEntryCharge_[face] = 0.0;
        directOuterCount_[face] = 0;
        directOuterCharge_[face] = 0.0;
        carrierOuterCount_[face] = 0;
        carrierOuterCharge_[face] = 0.0;
      }
    peakPmlCarriers_ = 0;
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
        particleBoundary_->beginStep();
        halo_->advanceMagnetic(*fields_);
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
    reportParticleBoundaryLosses();
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
    initializeSources();

    if (sources_.maxwellIncidentWaveCount() > 0 &&
        config_.boundary.type == EBBoundaryType::Cpml)
      throw std::runtime_error(
        "CPML is active, but TF/SF seed injection has not yet been placed inside the CPML interior. This stage intentionally validates the no-seed absorbing boundary first; refusing an overlapping boundary configuration.");
    if (config_.mesh.fieldSolver == EBMaxwellSolver::CowanZ &&
        sources_.maxwellIncidentWaveCount() > 0)
      throw std::runtime_error(
        "Cowan-z field advance is active, but the present TF/SF incident-wave correction is still Yee-specific. Select mesh.field_solver: yee for this run until the generalized Cowan TF/SF stencil is connected; refusing to inject a numerically inconsistent seed field.");

    particleCPML_.reset(new ParticleCPMLRegion(
      globalGeometry_, globalOriginBox_,
      config_.boundary.type == EBBoundaryType::Cpml ?
        config_.boundary.cpml : EBCPMLParameters()));
    initializeParticles();
    fields_.reset(new EBFieldGrid(localGeometry_));
    halo_.reset(new EBZSlabHaloExchange(communicator_));
    particleBoundary_.reset(new ParticleOpenBoundary(
      localGeometry_, localOriginBox_, globalGeometry_, globalOriginBox_));
    if (config_.boundary.type == EBBoundaryType::Cpml)
      fields_->setBoundary(std::unique_ptr<EBBoundaryOperator>(
        new EBConvolutionalPML(
          localGeometry_, globalGeometry_.nz, localZOffset_,
          config_.boundary.cpml,
          halo_->lowerRank() == MPI_PROC_NULL,
          halo_->upperRank() == MPI_PROC_NULL)));
    else
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
    const unsigned long long localBoundaryBytes =
      static_cast<unsigned long long>(
        fields_->memoryFootprint().boundaryBytes);
    unsigned long long totalBoundaryBytes = 0;
    unsigned long long maximumBoundaryBytes = 0;
    MPI_Reduce(&localBoundaryBytes, &totalBoundaryBytes, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator_);
    MPI_Reduce(&localBoundaryBytes, &maximumBoundaryBytes, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, communicator_);
    if (rank_ == 0)
      {
        logRoot(communicator_, config_.runtime.interactive() ?
          "Runtime strategy: interactive local-test mode; SIGINT/SIGTERM stops at a complete field step." :
          "Runtime strategy: HPC throughput mode; no signal polling or periodic durability flushes.");
        logRoot(communicator_,
          "Direct SI E/B simulation active; no A/phi state is allocated.");
        if (config_.mesh.fieldSolver == EBMaxwellSolver::CowanZ)
          {
            const EBCowanCoefficients& coefficient =
              fields_->cowanCoefficients();
            const Double sx = SI::c * globalGeometry_.dt /
                              globalGeometry_.dx;
            const Double sy = SI::c * globalGeometry_.dt /
                              globalGeometry_.dy;
            const Double sz = SI::c * globalGeometry_.dt /
                              globalGeometry_.dz;
            const Double samples = 16.0;
            std::ostringstream solverMessage;
            solverMessage << std::setprecision(10)
              << "Maxwell solver: Cowan-z controlled-dispersion FDTD; "
              << "axis Courant S=(" << sx << ", " << sy << ", "
              << sz << "), squared aspect r=(" << coefficient.ratio[0]
              << ", " << coefficient.ratio[1] << ", 1), where "
              << "r=((dz/dx)^2, (dz/dy)^2, 1).";
            logRoot(communicator_, solverMessage.str());
            std::ostringstream dispersionMessage;
            dispersionMessage << std::setprecision(10)
              << "Transverse vacuum dispersion at " << samples
              << " cells/wavelength: x vp/c="
              << EBFieldGrid::axisPhaseVelocityRatio(sx, samples)
              << ", vg/c="
              << EBFieldGrid::axisGroupVelocityRatio(sx, samples)
              << "; y vp/c="
              << EBFieldGrid::axisPhaseVelocityRatio(sy, samples)
              << ", vg/c="
              << EBFieldGrid::axisGroupVelocityRatio(sy, samples)
              << ". The z-axis values are exactly 1 for resolved vacuum modes.";
            logRoot(communicator_, dispersionMessage.str());
            std::ostringstream coefficientMessage;
            coefficientMessage << std::setprecision(10)
              << "Cowan smoothing coefficients: alpha=("
              << coefficient.alpha[0] << ", "
              << coefficient.alpha[1] << ", "
              << coefficient.alpha[2] << "), beta=("
              << coefficient.beta[0] << ", "
              << coefficient.beta[1] << ", "
              << coefficient.beta[2] << "), delta_xy="
              << coefficient.deltaXY << ", delta_yz="
              << coefficient.deltaYZ << ", delta_zx="
              << coefficient.deltaZX << ".";
            logRoot(communicator_, coefficientMessage.str());
            logRoot(communicator_,
              "Cowan-z stability guard passed: dz is the smallest spacing, the coefficient factorization is valid, and c*dt=dz. This is the analytic Courant boundary, so problem-scale convergence testing is still required.");
          }
        else
          logRoot(communicator_,
            "Maxwell solver: standard Yee FDTD regression mode.");
        if (config_.boundary.type == EBBoundaryType::Cpml)
          {
            std::ostringstream boundaryMessage;
            boundaryMessage << std::setprecision(8)
              << "Boundary: unsplit CFS-CPML, cells=("
              << config_.boundary.cpml.cells[0] << ", "
              << config_.boundary.cpml.cells[1] << ", "
              << config_.boundary.cpml.cells[2] << "), order="
              << config_.boundary.cpml.polynomialOrder
              << ", target reflection="
              << config_.boundary.cpml.targetReflection
              << ", kappa_max=" << config_.boundary.cpml.kappaMax
              << ", alpha_fraction="
              << config_.boundary.cpml.alphaFraction
              << "; auxiliary memory total="
              << static_cast<Double>(totalBoundaryBytes) /
                 (1024.0 * 1024.0)
              << " MiB, maximum rank="
              << static_cast<Double>(maximumBoundaryBytes) /
                 (1024.0 * 1024.0) << " MiB.";
            logRoot(communicator_, boundaryMessage.str());
          }
        else
          logRoot(communicator_,
            "Boundary: PEC regression mode; no absorbing-layer state allocated.");
        logRoot(communicator_,
          "Outer particle fallback: charge-conserving open absorption remains active on non-CPML faces and for residual carrier cleanup; production mode retains only face totals.");
        if (particleCPML_->enabled())
          logRoot(communicator_,
            "CPML particle policy: physical trajectories and detector participation end at the inner CPML surface; lightweight ballistic carriers continue with per-particle conductivity-matched current damping and no diagnostic output.");
        if (config_.particleRetirement.enabled)
          {
            std::ostringstream retirementMessage;
            retirementMessage << std::setprecision(10)
              << "EXPERIMENTAL particle retirement: lab z=["
              << config_.particleRetirement.entranceZ << ", "
              << config_.particleRetirement.exitZ()
              << "] m; physical particle push/output ends at entry and a "
                 "ballistic current carrier is tapered with a C2 quintic "
                 "profile. Net charge removal is not continuity exact; a "
                 "matched zero-radiation baseline is mandatory.";
            logRoot(communicator_, retirementMessage.str());
          }
        logRoot(communicator_,
          "WARNING: the initial Gauss-consistent particle field is not implemented; this is not yet a final radiation-production solver.");
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
    if (nz < 2 * static_cast<std::size_t>(size_))
      throw std::invalid_argument(
        "Direct E/B grid requires at least two z cells per MPI rank");

    if (config_.boundary.type == EBBoundaryType::Cpml)
      {
        const std::size_t globalCells[3] = {nx, ny, nz};
        for (unsigned int axis = 0; axis < 3; ++axis)
          if (config_.boundary.cpml.cells[axis] > 0 &&
              2 * config_.boundary.cpml.cells[axis] >= globalCells[axis])
            {
              std::ostringstream message;
              message << "CPML on axis " << "xyz"[axis]
                      << " uses " << config_.boundary.cpml.cells[axis]
                      << " cells per face, but the global mesh has only "
                      << globalCells[axis]
                      << " cells; opposite layers must leave a non-PML interior";
              throw std::invalid_argument(message.str());
            }
        const std::size_t baseCells = nz /
          static_cast<std::size_t>(size_);
        const std::size_t remainder = nz %
          static_cast<std::size_t>(size_);
        const std::size_t lowerEndCells = baseCells +
          (remainder > 0 ? 1 : 0);
        const std::size_t upperEndCells = baseCells +
          (static_cast<std::size_t>(size_ - 1) < remainder ? 1 : 0);
        if (config_.boundary.cpml.cells[2] > lowerEndCells ||
            config_.boundary.cpml.cells[2] > upperEndCells)
          {
            std::ostringstream message;
            message << "The current z-slab CPML implementation keeps each z absorbing layer on one endpoint MPI rank. cpml cells z="
                    << config_.boundary.cpml.cells[2]
                    << ", endpoint slab cells=(" << lowerEndCells
                    << ", " << upperEndCells
                    << "); reduce MPI ranks, increase nz, or reduce the z CPML thickness.";
            throw std::invalid_argument(message.str());
          }
      }

    if (config_.mesh.fieldSolver == EBMaxwellSolver::CowanZ)
      {
        const Double tolerance = 64.0 *
          std::numeric_limits<Double>::epsilon() *
          std::max(dx, std::max(dy, dz));
        if (dx + tolerance < dz || dy + tolerance < dz)
          {
            std::ostringstream message;
            message << "Cowan-z requires z to have the smallest grid spacing so that c*dt=dz is stable and dispersion-free. ";
            if (dx + tolerance < dz)
              message << transverseGridRecommendation(
                "x", config_.mesh.lengths[0], dz,
                config_.inputUnits.length) << ". ";
            if (dy + tolerance < dz)
              message << transverseGridRecommendation(
                "y", config_.mesh.lengths[1], dz,
                config_.inputUnits.length) << ".";
            throw std::invalid_argument(message.str());
          }
      }

    Double dt = 0.0;
    if (config_.mesh.fieldSolver == EBMaxwellSolver::CowanZ)
      dt = dz / SI::c;
    else
      {
        const Double inverseSpacingSquared =
          1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz);
        dt = 0.95 / (SI::c * std::sqrt(inverseSpacingSquared));
      }
    globalGeometry_ = EBGridGeometry(
      nx, ny, nz, dx, dy, dz, dt, config_.mesh.fieldSolver);
    const Slab local = slabForRank(rank_);
    localZOffset_ = local.offset;
    localGeometry_ = EBGridGeometry(nx, ny, local.cells,
                                    dx, dy, dz, dt,
                                    config_.mesh.fieldSolver);

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

    initializeFieldDetectorRegions();
    validateBeamlineExclusionRules();
    const Double firstPhysicalEntrance =
      firstBeamlinePhysicalEntranceLab();
    const Double firstInteractionEntrance =
      firstBeamlineInteractionEntranceLab();

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

  void Simulation::initializeFieldDetectorRegions()
  {
    const unsigned long long localParticleCount =
      static_cast<unsigned long long>(particles_.size());
    unsigned long long globalParticleCount = 0;
    MPI_Allreduce(&localParticleCount, &globalParticleCount, 1,
                  MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator_);

    Double localMaximumGamma = 1.0;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      localMaximumGamma = std::max(localMaximumGamma,
        BoostFrameTransform::gammaFromProperVelocity(
          particles_[index].properVelocity));
    Double maximumGamma = 1.0;
    MPI_Allreduce(&localMaximumGamma, &maximumGamma, 1, MPI_DOUBLE,
                  MPI_MAX, communicator_);
    if (!(maximumGamma >= 1.0) || !std::isfinite(maximumGamma))
      throw std::runtime_error(
        "Cannot derive a finite laboratory gamma for field-detector references");

    /* size_x and size_y are full mesh widths.  Using their complete diagonal
     * is deliberately more conservative than the centre-to-corner radius and
     * also leaves room for a modest off-axis bunch envelope. */
    const Double rho = std::hypot(config_.mesh.lengths[0],
                                  config_.mesh.lengths[1]);
    if (!(rho > 0.0) || !std::isfinite(rho))
      throw std::runtime_error(
        "Field-detector reference radius must be finite and positive");

    for (std::size_t detector = 0;
         detector < detectorConfig_.fieldPlanes.size(); ++detector)
      {
        FieldDetectorPlaneConfig& plane =
          detectorConfig_.fieldPlanes[detector];
        plane.referenceRho = rho;
        plane.referenceGamma = maximumGamma;
        plane.referenceDistance = plane.particleBackgroundReference ?
          maximumGamma * rho : 0.0;
        plane.referenceEntranceZ = plane.z - plane.referenceDistance;
        if (plane.retirementFrequencyProtection)
          {
            if (!config_.particleRetirement.enabled)
              throw std::logic_error(
                "Frequency-protected detector has no particle retirement");
            const Double inverseGamma2 =
              1.0 / (maximumGamma * maximumGamma);
            const Double betaSquared = 1.0 - inverseGamma2;
            if (!(betaSquared > 0.0))
              throw std::runtime_error(
                "Retirement frequency protection requires moving relativistic particles");
            const Double beta = std::sqrt(betaSquared);
            /* Stable evaluation of 1/beta-1 avoids subtracting two nearly
             * equal numbers at the high gamma values this guard targets. */
            const Double delayFactor =
              inverseGamma2 / (beta * (1.0 + beta));
            plane.retirementRequiredLength =
              static_cast<Double>(plane.retirementMinimumCycles) *
              PLANCK_EV_SECOND * SI::c /
              (plane.retirementMinimumPhotonEnergyEV * delayFactor);
            plane.retirementObservedCycles =
              config_.particleRetirement.length * delayFactor *
              plane.retirementMinimumPhotonEnergyEV /
              (SI::c * PLANCK_EV_SECOND);
            plane.retirementCausalDistance = maximumGamma * rho;
            plane.retirementCausalEntranceZ =
              plane.z - plane.retirementCausalDistance;
            const Double marginLab =
              config_.mesh.boostGamma * globalGeometry_.dz;
            const Double lengthTolerance = 128.0 *
              std::numeric_limits<Double>::epsilon() *
              std::max(1.0, plane.retirementRequiredLength);
            if (config_.particleRetirement.length + lengthTolerance <
                plane.retirementRequiredLength)
              {
                std::ostringstream message;
                message << std::setprecision(10)
                        << "Field detector '" << plane.name
                        << "' retirement frequency protection requires length "
                        << plane.retirementRequiredLength << " m ("
                        << plane.retirementRequiredLength /
                             config_.inputUnits.length
                        << " in the configured length unit) for "
                        << plane.retirementMinimumCycles << " cycles at "
                        << plane.retirementMinimumPhotonEnergyEV
                        << " eV using maximum lab gamma=" << maximumGamma
                        << "; configured length is "
                        << config_.particleRetirement.length << " m.";
                throw std::runtime_error(message.str());
              }
            if (config_.particleRetirement.exitZ() + marginLab >
                plane.retirementCausalEntranceZ)
              {
                const Double recommended = config_.particleRetirement.exitZ() +
                  plane.retirementCausalDistance + marginLab;
                std::ostringstream message;
                message << std::setprecision(10)
                        << "Field detector '" << plane.name
                        << "' causal guard starts at z="
                        << plane.retirementCausalEntranceZ
                        << " m before the retirement exit clearance. Move "
                           "detector z to at least " << recommended << " m ("
                        << recommended / config_.inputUnits.length
                        << " in the configured length unit).";
                throw std::runtime_error(message.str());
              }
          }
        if (plane.referenceValidation &&
            globalParticleCount > static_cast<unsigned long long>(
              plane.referenceValidationMaximumParticles))
          {
            std::ostringstream message;
            message << "Field detector '" << plane.name
                    << "' two-plane validation is limited to "
                    << plane.referenceValidationMaximumParticles
                    << " particles, but the input contains "
                    << globalParticleCount
                    << ". Disable particle_background.validation for the "
                       "production run or raise maximum_particles explicitly.";
            throw std::runtime_error(message.str());
          }
      }

    std::size_t fieldIndex = 0;
    for (std::size_t element = 0;
         element < beamlineElements_.size(); ++element)
      if (beamlineElements_[element].role ==
          BeamlineElementRole::FieldDetectorPlane)
        {
          if (fieldIndex >= detectorConfig_.fieldPlanes.size())
            throw std::logic_error(
              "Field-detector beamline extent count is inconsistent");
          const FieldDetectorPlaneConfig& plane =
            detectorConfig_.fieldPlanes[fieldIndex++];
          if (plane.particleBackgroundReference)
            beamlineElements_[element].interactionEntrance =
              plane.referenceEntranceZ;
          else if (plane.retirementFrequencyProtection)
            beamlineElements_[element].interactionEntrance =
              config_.particleRetirement.entranceZ;
          else
            beamlineElements_[element].interactionEntrance = plane.z;
          beamlineElements_[element].interactionExit = plane.z;
        }
    if (fieldIndex != detectorConfig_.fieldPlanes.size())
      throw std::logic_error(
        "Field-detector configuration count is inconsistent");

    if (rank_ == 0)
      for (std::size_t detector = 0;
           detector < detectorConfig_.fieldPlanes.size(); ++detector)
        {
          const FieldDetectorPlaneConfig& plane =
            detectorConfig_.fieldPlanes[detector];
          std::ostringstream message;
          message << std::setprecision(10) << "Field detector '"
                  << plane.name << "': ";
          if (plane.particleBackgroundReference)
            {
              message
                << "diagnostic-only ballistic reference region lab z=["
                << plane.referenceEntranceZ << ", " << plane.z
                << "] m, length=" << plane.referenceDistance
                << " m, rho_guard=" << plane.referenceRho
                << " m, gamma_guard=" << plane.referenceGamma
                << ". Physical particle and Maxwell evolution are unchanged.";
              if (plane.referenceValidation)
                message << " Two-plane straight-line validation is enabled for "
                        << globalParticleCount << " particles (limit "
                        << plane.referenceValidationMaximumParticles << ").";
            }
          else if (plane.retirementFrequencyProtection)
            message
              << "frequency-protected retirement region lab z=["
              << config_.particleRetirement.entranceZ << ", " << plane.z
              << "] m; configured/required taper length="
              << config_.particleRetirement.length << "/"
              << plane.retirementRequiredLength << " m, protected band starts "
              << "at " << plane.retirementMinimumPhotonEnergyEV
              << " eV, requested/actual cycles="
              << plane.retirementMinimumCycles << "/"
              << plane.retirementObservedCycles
              << ", causal guard lab z=[" << plane.retirementCausalEntranceZ
              << ", " << plane.z << "] m, gamma_guard="
              << plane.referenceGamma << ".";
          else
            message << "raw laboratory field plane; no particle-background "
                       "diagnostic or retirement-frequency guard.";
          logRoot(communicator_, message.str());
        }
  }

  void Simulation::validateBeamlineExclusionRules() const
  {
    const Double marginLab = config_.mesh.boostGamma * globalGeometry_.dz;
    const auto overlaps = [](Double firstEntrance, Double firstExit,
                             Double secondEntrance, Double secondExit)
      {
        const Double scale = std::max(1.0,
          std::max(std::abs(firstEntrance),
            std::max(std::abs(firstExit),
              std::max(std::abs(secondEntrance), std::abs(secondExit)))));
        const Double tolerance = 128.0 *
          std::numeric_limits<Double>::epsilon() * scale;
        return std::max(firstEntrance, secondEntrance) <
               std::min(firstExit, secondExit) - tolerance;
      };

    /* Magnetic interaction regions are mutually exclusive.  Touching compact
     * support boundaries are allowed; overlapping fringe supports are not. */
    for (std::size_t first = 0; first < config_.magnets.size(); ++first)
      for (std::size_t second = first + 1;
           second < config_.magnets.size(); ++second)
        if (overlaps(config_.magnets[first].interactionEntranceLab(),
                     config_.magnets[first].interactionExitLab(),
                     config_.magnets[second].interactionEntranceLab(),
                     config_.magnets[second].interactionExitLab()))
          {
            std::ostringstream message;
            message << "Magnetic element interaction regions overlap: magnet "
                    << first << " lab z=["
                    << config_.magnets[first].interactionEntranceLab() << ", "
                    << config_.magnets[first].interactionExitLab()
                    << "] m and magnet " << second << " lab z=["
                    << config_.magnets[second].interactionEntranceLab() << ", "
                    << config_.magnets[second].interactionExitLab()
                    << "] m. Separate the physical devices or reduce their "
                       "configured fringe support.";
            throw std::runtime_error(message.str());
          }

    Double lastMagneticExit = -std::numeric_limits<Double>::infinity();
    for (std::size_t magnet = 0; magnet < config_.magnets.size(); ++magnet)
      lastMagneticExit = std::max(lastMagneticExit,
        config_.magnets[magnet].interactionExitLab());

    for (std::size_t detector = 0;
         detector < detectorConfig_.fieldPlanes.size(); ++detector)
      {
        const FieldDetectorPlaneConfig& plane =
          detectorConfig_.fieldPlanes[detector];
        if (!plane.particleBackgroundReference &&
            !plane.retirementFrequencyProtection)
          continue;
        const Double regionEntrance = plane.particleBackgroundReference ?
          plane.referenceEntranceZ : config_.particleRetirement.entranceZ;
        for (std::size_t magnet = 0;
             magnet < config_.magnets.size(); ++magnet)
          if (overlaps(regionEntrance, plane.z,
                       config_.magnets[magnet].interactionEntranceLab(),
                       config_.magnets[magnet].interactionExitLab()))
            {
              const Double recommended = plane.particleBackgroundReference ?
                lastMagneticExit + plane.referenceDistance + marginLab :
                lastMagneticExit + marginLab;
              std::ostringstream message;
              message << std::setprecision(10)
                      << "Field detector '" << plane.name
                      << "' has a "
                      << (plane.particleBackgroundReference ?
                          "diagnostic ballistic-reference" :
                          "frequency-protected retirement")
                      << " region lab z=[" << regionEntrance << ", "
                      << plane.z
                      << "] m that overlaps magnetic element " << magnet
                      << " interaction region lab z=["
                      << config_.magnets[magnet].interactionEntranceLab()
                      << ", "
                      << config_.magnets[magnet].interactionExitLab()
                      << "] m. "
                      << (plane.particleBackgroundReference ?
                          "Move detector z" : "Move retirement entrance_z")
                      << " to at least " << recommended << " m ("
                      << recommended / config_.inputUnits.length
                      << " in the configured length unit). Field-detector "
                         "regions may overlap each other; particle detector "
                         "planes do not participate in exclusion checks.";
              throw std::runtime_error(message.str());
            }
      }
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
    if (!detectorConfig_.enabled()) return;
    detectors_.reset(new LabDetectorManager(
      detectorConfig_, globalGeometry_,
      globalOriginBox_, localOriginBox_, frame_, communicator_));
    if (rank_ == 0)
      {
        std::ostringstream message;
        message << "Laboratory detector planes active: fields="
                << detectorConfig_.fieldPlanes.size()
                << ", particles="
                << detectorConfig_.particlePlanes.size()
                << "; detector HDF5 is written only by MPI rank zero.";
        logRoot(communicator_, message.str());
      }
  }

  void Simulation::sampleTrajectory()
  {
    if (!trajectoryWriter_.isOpen()) return;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      appendTrajectoryEvent(
        particles_[index], timeBoxSI_, TrajectoryEvent::Sample,
        ParticleBoundaryFace::LowerX);
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
      std::max(trajectoryRhythmSI_,
        std::max(std::abs(timeBoxSI_),
                 std::abs(nextTrajectorySampleTime_)));
    do
      nextTrajectorySampleTime_ += trajectoryRhythmSI_;
    while (nextTrajectorySampleTime_ <= timeBoxSI_ + tolerance);
  }

  void Simulation::appendTrajectoryEvent(
      const RelativisticParticleSI& particle,
      Double timeBox,
      TrajectoryEvent event,
      ParticleBoundaryFace face)
  {
    if (!trajectoryWriter_.isOpen()) return;
    TrajectoryRecord record;
    record.particleId = particle.id;
    record.sourceId = particle.sourceId;
    frame_.boxToLab(timeBox, particle.position[2],
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
    record.event = static_cast<std::uint8_t>(event);
    record.boundaryFace =
      (event == TrajectoryEvent::CpmlEntry ||
       event == TrajectoryEvent::DomainExit) ?
      static_cast<std::int8_t>(face) : -1;
    trajectoryWriter_.append(record);
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
    std::vector<ParticleCPMLCarrier> retainedCarriers;
    retainedCarriers.reserve(pmlCarriers_.size());
    std::vector<ParticleCPMLCarrier> retainedRetirementCarriers;
    retainedRetirementCarriers.reserve(retirementCarriers_.size());
    std::vector<TransferPacket> sendLower;
    std::vector<TransferPacket> sendUpper;
    std::vector<CarrierTransferPacket> sendCarrierLower;
    std::vector<CarrierTransferPacket> sendCarrierUpper;
    std::vector<CarrierTransferPacket> sendRetirementLower;
    std::vector<CarrierTransferPacket> sendRetirementUpper;
    const Double lowerZ = localOriginBox_[2];
    const Double upperZ = lowerZ +
      static_cast<Double>(localGeometry_.nz) * localGeometry_.dz;
    const Double crossingTolerance = 64.0 *
      std::numeric_limits<Double>::epsilon();
    const Double carrierCutoff = 64.0 *
      std::numeric_limits<Double>::epsilon();
    unsigned long long localLost = 0;

    const auto interfaceCrossing = [&](const FieldVector<Double>& start,
                                       const FieldVector<Double>& end,
                                       bool& crossLower,
                                       bool& crossUpper)
    {
      crossLower = rank_ > 0 && end[2] < lowerZ;
      crossUpper = rank_ + 1 < size_ && end[2] >= upperZ;
      if (!crossLower && !crossUpper) return 1.0;
      const Double boundary = crossLower ? lowerZ : upperZ;
      const Double fraction = (boundary - start[2]) /
        (end[2] - start[2]);
      if (!(fraction >= 0.0 && fraction <= 1.0) ||
          !std::isfinite(fraction))
        throw std::runtime_error(
          "Invalid MPI-interface particle crossing");
      return fraction;
    };

    const auto labPlaneCrossing = [&](
        const FieldVector<Double>& start,
        const FieldVector<Double>& end,
        Double startFraction, Double planeZ,
        ParticleBoundaryHit& hit)
    {
      if (!config_.particleRetirement.enabled) return false;
      const Double startTime = timeBoxSI_ +
        startFraction * globalGeometry_.dt;
      const Double endTime = timeBoxSI_ + globalGeometry_.dt;
      const Double startZ = frame_.labZFromBoxZT(start[2], startTime);
      const Double endZ = frame_.labZFromBoxZT(end[2], endTime);
      const Double scale = std::max(1.0,
        std::max(std::abs(startZ), std::max(std::abs(endZ),
                                            std::abs(planeZ))));
      const Double tolerance = 128.0 *
        std::numeric_limits<Double>::epsilon() * scale;
      if (!(endZ > startZ + tolerance) ||
          startZ >= planeZ - tolerance || endZ < planeZ - tolerance)
        return false;
      hit.fraction = std::max(0.0, std::min(1.0,
        (planeZ - startZ) / (endZ - startZ)));
      /* Retirement planes are not domain faces, but initialize the shared
       * hit record so passing it through generic event selection never reads
       * an indeterminate enum value. Retirement output stores face=-1. */
      hit.face = ParticleBoundaryFace::UpperZ;
      for (unsigned int axis = 0; axis < 3; ++axis)
        hit.position[axis] = start[axis] + hit.fraction *
          (end[axis] - start[axis]);
      return true;
    };

    const auto retirementPrimitive = [](Double coordinate)
    {
      const Double s = std::max(0.0, std::min(1.0, coordinate));
      const Double s2 = s * s;
      const Double s4 = s2 * s2;
      return s - 2.5 * s4 + 3.0 * s4 * s - s4 * s2;
    };

    const auto retirementWeight = [&](Double labZ)
    {
      const Double s = (labZ - config_.particleRetirement.entranceZ) /
        config_.particleRetirement.length;
      if (s <= 0.0) return 1.0;
      if (s >= 1.0) return 0.0;
      const Double s2 = s * s;
      const Double s3 = s2 * s;
      return 1.0 - 10.0 * s3 + 15.0 * s3 * s - 6.0 * s3 * s2;
    };

    const auto retirementAverageWeight = [&](Double firstLabZ,
                                               Double secondLabZ)
    {
      const Double first = (firstLabZ -
        config_.particleRetirement.entranceZ) /
        config_.particleRetirement.length;
      const Double second = (secondLabZ -
        config_.particleRetirement.entranceZ) /
        config_.particleRetirement.length;
      if (std::abs(second - first) <=
          64.0 * std::numeric_limits<Double>::epsilon())
        return retirementWeight(0.5 * (firstLabZ + secondLabZ));
      return (retirementPrimitive(second) -
              retirementPrimitive(first)) / (second - first);
    };

    const auto processCarrierSegment =
      [&](ParticleCPMLCarrier carrier,
          const FieldVector<Double>& segmentStart,
          Double segmentStartFraction,
          bool allowTransfer)
    {
      if (!(segmentStartFraction >= 0.0 &&
            segmentStartFraction <= 1.0) ||
          !std::isfinite(segmentStartFraction))
        throw std::runtime_error(
          "Invalid CPML-carrier segment time fraction");
      if (!(carrier.currentWeight >= 0.0) ||
          !std::isfinite(carrier.currentWeight))
        throw std::runtime_error("Invalid CPML carrier current weight");
      if (carrier.currentWeight <= carrierCutoff) return;

      bool crossLower = false;
      bool crossUpper = false;
      const Double interfaceFraction = interfaceCrossing(
        segmentStart, carrier.position, crossLower, crossUpper);
      const bool crossesInterface = crossLower || crossUpper;
      ParticleBoundaryHit outerHit;
      const bool exitsOuter = particleBoundary_->firstExit(
        segmentStart, carrier.position, outerHit);

      Double localEndFraction = 1.0;
      FieldVector<Double> segmentEnd(carrier.position);
      enum CarrierEnd { RetainCarrier, TransferCarrier, RemoveCarrier };
      CarrierEnd endAction = RetainCarrier;
      if (exitsOuter &&
          (!crossesInterface ||
           outerHit.fraction <= interfaceFraction + crossingTolerance))
        {
          localEndFraction = outerHit.fraction;
          segmentEnd = outerHit.position;
          endAction = RemoveCarrier;
        }
      else if (crossesInterface)
        {
          if (!allowTransfer)
            throw std::runtime_error(
              "A CPML carrier crossed more than one MPI interface in one field step");
          localEndFraction = interfaceFraction;
          for (unsigned int axis = 0; axis < 3; ++axis)
            segmentEnd[axis] = segmentStart[axis] +
              localEndFraction *
              (carrier.position[axis] - segmentStart[axis]);
          endAction = TransferCarrier;
        }

      const Double globalEndFraction = segmentStartFraction +
        localEndFraction * (1.0 - segmentStartFraction);
      const Double duration = (globalEndFraction - segmentStartFraction) *
        globalGeometry_.dt;
      FieldVector<Double> midpoint(0.0);
      for (unsigned int axis = 0; axis < 3; ++axis)
        midpoint[axis] = 0.5 *
          (segmentStart[axis] + segmentEnd[axis]);
      const Double midpointDamping =
        particleCPML_->currentDampingFactor(
          segmentStart, midpoint, 0.5 * duration);
      const Double endDamping = particleCPML_->currentDampingFactor(
        segmentStart, segmentEnd, duration);
      const Double depositedWeight =
        carrier.currentWeight * midpointDamping;
      const Double endWeight = carrier.currentWeight * endDamping;
      depositor.depositSegment(segmentStart, segmentEnd,
                               carrier.charge * depositedWeight);

      if (endAction == RemoveCarrier)
        {
          particleBoundary_->depositOutgoingFlux(
            outerHit, carrier.charge * endWeight);
          const std::size_t face =
            static_cast<std::size_t>(outerHit.face);
          ++carrierOuterCount_[face];
          carrierOuterCharge_[face] += carrier.charge * endWeight;
          return;
        }
      if (endWeight <= carrierCutoff) return;

      carrier.currentWeight = endWeight;
      if (endAction == TransferCarrier)
        {
          const CarrierTransferPacket packet = packCarrier(
            carrier, segmentEnd, globalEndFraction);
          (crossLower ? sendCarrierLower : sendCarrierUpper).push_back(
            packet);
        }
      else
        {
          carrier.position = segmentEnd;
          retainedCarriers.push_back(carrier);
        }
    };

    const auto processRetirementSegment = [&](
        ParticleCPMLCarrier carrier,
        const FieldVector<Double>& segmentStart,
        Double segmentStartFraction, bool allowTransfer)
    {
      if (!(segmentStartFraction >= 0.0 &&
            segmentStartFraction <= 1.0) ||
          !std::isfinite(segmentStartFraction))
        throw std::runtime_error(
          "Invalid retirement-carrier segment time fraction");

      bool crossLower = false;
      bool crossUpper = false;
      const Double interfaceFraction = interfaceCrossing(
        segmentStart, carrier.position, crossLower, crossUpper);
      const bool crossesInterface = crossLower || crossUpper;
      ParticleBoundaryHit cpmlHit;
      ParticleBoundaryHit outerHit;
      ParticleBoundaryHit retirementExitHit;
      const bool entersCpml = particleCPML_->firstEntry(
        segmentStart, carrier.position, cpmlHit);
      const bool exitsOuter = particleBoundary_->firstExit(
        segmentStart, carrier.position, outerHit);
      const bool reachesRetirementExit = labPlaneCrossing(
        segmentStart, carrier.position, segmentStartFraction,
        config_.particleRetirement.exitZ(), retirementExitHit);

      enum RetirementEnd { RetainRetirement, TransferRetirement,
                           EnterCpmlFromRetirement, ExitOuterFromRetirement,
                           CompleteRetirement };
      RetirementEnd endAction = RetainRetirement;
      Double localEndFraction = 1.0;
      FieldVector<Double> segmentEnd(carrier.position);
      ParticleBoundaryHit selectedHit;
      const auto select = [&](Double fraction, RetirementEnd action,
                              const ParticleBoundaryHit& hit)
      {
        if (fraction <= localEndFraction + crossingTolerance)
          {
            localEndFraction = fraction;
            endAction = action;
            selectedHit = hit;
            segmentEnd = hit.position;
          }
      };
      if (crossesInterface)
        {
          if (!allowTransfer)
            throw std::runtime_error(
              "A retirement carrier crossed more than one MPI interface in one field step");
          ParticleBoundaryHit hit;
          hit.fraction = interfaceFraction;
          for (unsigned int axis = 0; axis < 3; ++axis)
            hit.position[axis] = segmentStart[axis] +
              interfaceFraction *
              (carrier.position[axis] - segmentStart[axis]);
          select(interfaceFraction, TransferRetirement, hit);
        }
      if (exitsOuter)
        select(outerHit.fraction, ExitOuterFromRetirement, outerHit);
      if (entersCpml)
        select(cpmlHit.fraction, EnterCpmlFromRetirement, cpmlHit);
      if (reachesRetirementExit)
        select(retirementExitHit.fraction, CompleteRetirement,
               retirementExitHit);

      const Double globalEndFraction = segmentStartFraction +
        localEndFraction * (1.0 - segmentStartFraction);
      const Double firstTime = timeBoxSI_ +
        segmentStartFraction * globalGeometry_.dt;
      const Double secondTime = timeBoxSI_ +
        globalEndFraction * globalGeometry_.dt;
      const Double firstLabZ = frame_.labZFromBoxZT(
        segmentStart[2], firstTime);
      const Double secondLabZ = frame_.labZFromBoxZT(
        segmentEnd[2], secondTime);
      const Double depositedWeight = std::max(0.0, std::min(1.0,
        retirementAverageWeight(firstLabZ, secondLabZ)));
      const Double endWeight = retirementWeight(secondLabZ);
      depositor.depositSegment(segmentStart, segmentEnd,
                               carrier.charge * depositedWeight);

      if (endAction == CompleteRetirement)
        {
          ++retirementExitCount_;
          retirementExitResidualCharge_ += carrier.charge * endWeight;
          return;
        }
      if (endAction == ExitOuterFromRetirement)
        {
          particleBoundary_->depositOutgoingFlux(
            selectedHit, carrier.charge * endWeight);
          const std::size_t face =
            static_cast<std::size_t>(selectedHit.face);
          ++carrierOuterCount_[face];
          carrierOuterCharge_[face] += carrier.charge * endWeight;
          return;
        }
      carrier.currentWeight = endWeight;
      if (endAction == EnterCpmlFromRetirement)
        {
          processCarrierSegment(carrier, selectedHit.position,
                                globalEndFraction, allowTransfer);
        }
      else if (endAction == TransferRetirement)
        {
          const CarrierTransferPacket packet = packCarrier(
            carrier, segmentEnd, globalEndFraction);
          (crossLower ? sendRetirementLower :
                        sendRetirementUpper).push_back(packet);
        }
      else
        {
          carrier.position = segmentEnd;
          retainedRetirementCarriers.push_back(carrier);
        }
    };

    for (std::size_t index = 0; index < pmlCarriers_.size(); ++index)
      {
        ParticleCPMLCarrier carrier = pmlCarriers_[index];
        const FieldVector<Double> start(carrier.position);
        carrier.position.pmv(globalGeometry_.dt, carrier.velocity);
        processCarrierSegment(carrier, start, 0.0, true);
      }
    for (std::size_t index = 0;
         index < retirementCarriers_.size(); ++index)
      {
        ParticleCPMLCarrier carrier = retirementCarriers_[index];
        const FieldVector<Double> start(carrier.position);
        carrier.position.pmv(globalGeometry_.dt, carrier.velocity);
        processRetirementSegment(carrier, start, 0.0, true);
      }

    const auto processActiveSegment =
      [&](RelativisticParticleSI particle,
          const FieldVector<Double>& segmentStart,
          Double segmentStartFraction,
          bool allowTransfer)
    {
      bool crossLower = false;
      bool crossUpper = false;
      const Double interfaceFraction = interfaceCrossing(
        segmentStart, particle.position, crossLower, crossUpper);
      const bool crossesInterface = crossLower || crossUpper;
      ParticleBoundaryHit cpmlHit;
      ParticleBoundaryHit outerHit;
      ParticleBoundaryHit retirementHit;
      const bool entersCpml = particleCPML_->firstEntry(
        segmentStart, particle.position, cpmlHit);
      const bool exitsOuter = particleBoundary_->firstExit(
        segmentStart, particle.position, outerHit);
      const bool entersRetirement = labPlaneCrossing(
        segmentStart, particle.position, segmentStartFraction,
        config_.particleRetirement.entranceZ, retirementHit);

      enum ActiveEnd { RetainActive, TransferActive,
                       EnterCpml, ExitOuter, EnterRetirement };
      ActiveEnd endAction = RetainActive;
      Double localEndFraction = 1.0;
      FieldVector<Double> segmentEnd(particle.position);
      ParticleBoundaryHit selectedHit;
      if (entersRetirement &&
          (!entersCpml || retirementHit.fraction <=
            cpmlHit.fraction + crossingTolerance) &&
          (!exitsOuter || retirementHit.fraction <=
            outerHit.fraction + crossingTolerance) &&
          (!crossesInterface || retirementHit.fraction <=
            interfaceFraction + crossingTolerance))
        {
          endAction = EnterRetirement;
          localEndFraction = retirementHit.fraction;
          segmentEnd = retirementHit.position;
          selectedHit = retirementHit;
        }
      else if (entersCpml &&
          (!entersRetirement || cpmlHit.fraction <=
            retirementHit.fraction + crossingTolerance) &&
          (!exitsOuter ||
           cpmlHit.fraction <= outerHit.fraction + crossingTolerance) &&
          (!crossesInterface ||
           cpmlHit.fraction <= interfaceFraction + crossingTolerance))
        {
          endAction = EnterCpml;
          localEndFraction = cpmlHit.fraction;
          segmentEnd = cpmlHit.position;
          selectedHit = cpmlHit;
        }
      else if (exitsOuter &&
               (!entersRetirement || outerHit.fraction <=
                 retirementHit.fraction + crossingTolerance) &&
               (!crossesInterface ||
                outerHit.fraction <=
                  interfaceFraction + crossingTolerance))
        {
          endAction = ExitOuter;
          localEndFraction = outerHit.fraction;
          segmentEnd = outerHit.position;
          selectedHit = outerHit;
        }
      else if (crossesInterface)
        {
          if (!allowTransfer)
            throw std::runtime_error(
              "An active particle crossed more than one MPI interface in one field step");
          endAction = TransferActive;
          localEndFraction = interfaceFraction;
          for (unsigned int axis = 0; axis < 3; ++axis)
            segmentEnd[axis] = segmentStart[axis] +
              localEndFraction *
              (particle.position[axis] - segmentStart[axis]);
        }

      depositor.depositSegment(segmentStart, segmentEnd, particle.charge);
      const Double globalEndFraction = segmentStartFraction +
        localEndFraction * (1.0 - segmentStartFraction);
      if (endAction == EnterCpml)
        {
          const std::size_t face =
            static_cast<std::size_t>(selectedHit.face);
          ++cpmlEntryCount_[face];
          cpmlEntryCharge_[face] += particle.charge;
          ++localLost;

          ParticleCPMLCarrier carrier;
          carrier.position = particle.position;
          const Double gamma =
            BoostFrameTransform::gammaFromProperVelocity(
              particle.properVelocity);
          carrier.velocity = particle.properVelocity;
          carrier.velocity *= SI::c / gamma;
          carrier.charge = particle.charge;
          carrier.currentWeight = 1.0;
          processCarrierSegment(
            carrier, selectedHit.position, globalEndFraction,
            allowTransfer);
        }
      else if (endAction == EnterRetirement)
        {
          ++retirementEntryCount_;
          retirementEntryCharge_ += particle.charge;
          ParticleCPMLCarrier carrier;
          carrier.position = particle.position;
          const Double gamma =
            BoostFrameTransform::gammaFromProperVelocity(
              particle.properVelocity);
          carrier.velocity = particle.properVelocity;
          carrier.velocity *= SI::c / gamma;
          carrier.charge = particle.charge;
          carrier.currentWeight = 1.0;
          processRetirementSegment(
            carrier, selectedHit.position, globalEndFraction,
            allowTransfer);
        }
      else if (endAction == ExitOuter)
        {
          particleBoundary_->depositOutgoingFlux(
            selectedHit, particle.charge);
          const std::size_t face =
            static_cast<std::size_t>(selectedHit.face);
          ++directOuterCount_[face];
          directOuterCharge_[face] += particle.charge;
          ++localLost;
        }
      else if (endAction == TransferActive)
        {
          TransferPacket transfer = {};
          transfer.particle = packParticle(particle);
          for (unsigned int axis = 0; axis < 3; ++axis)
            transfer.segmentStart[axis] = segmentEnd[axis];
          transfer.segmentStartFraction = globalEndFraction;
          (crossLower ? sendLower : sendUpper).push_back(transfer);
        }
      else
        retained.push_back(particle);
    };

    for (std::size_t index = 0; index < particles_.size(); ++index)
      {
        RelativisticParticleSI particle = particles_[index];
        const FieldVector<Double> start(particle.position);
        RelativisticBorisPusher::pushFromGridAndPrescribedLab(
          particle, *fields_, localOriginBox_, sources_, frame_,
          timeBoxSI_, globalGeometry_.dt);

        ParticleBoundaryHit diagnosticHit;
        ParticleBoundaryHit cpmlHit;
        ParticleBoundaryHit outerHit;
        ParticleBoundaryHit retirementHit;
        const bool entersCpml = particleCPML_->firstEntry(
          start, particle.position, cpmlHit);
        const bool exitsOuter = particleBoundary_->firstExit(
          start, particle.position, outerHit);
        const bool entersRetirement = labPlaneCrossing(
          start, particle.position, 0.0,
          config_.particleRetirement.entranceZ, retirementHit);
        bool hasTerminalEvent = false;
        TrajectoryEvent terminalEvent = TrajectoryEvent::DomainExit;
        Double terminalFraction = 1.0;
        const auto selectTerminal = [&](
            bool present, const ParticleBoundaryHit& hit,
            TrajectoryEvent event)
        {
          if (present &&
              hit.fraction <= terminalFraction + crossingTolerance)
            {
              hasTerminalEvent = true;
              terminalFraction = hit.fraction;
              diagnosticHit = hit;
              terminalEvent = event;
            }
        };
        selectTerminal(exitsOuter, outerHit,
                       TrajectoryEvent::DomainExit);
        selectTerminal(entersCpml, cpmlHit,
                       TrajectoryEvent::CpmlEntry);
        selectTerminal(entersRetirement, retirementHit,
                       TrajectoryEvent::RetirementEntry);

        RelativisticParticleSI diagnosticEnd(particle);
        Double diagnosticTime = timeBoxSI_ + globalGeometry_.dt;
        if (hasTerminalEvent)
          {
            diagnosticEnd.position = diagnosticHit.position;
            diagnosticTime = timeBoxSI_ +
              diagnosticHit.fraction * globalGeometry_.dt;
          }
        if (detectors_)
          detectors_->captureParticleStep(
            particles_[index], diagnosticEnd,
            timeBoxSI_, diagnosticTime);
        if (hasTerminalEvent && config_.trajectory.enabled)
          {
            appendTrajectoryEvent(
              diagnosticEnd, diagnosticTime,
              terminalEvent, diagnosticHit.face);
          }
        processActiveSegment(particle, start, 0.0, true);
      }

    std::vector<TransferPacket> receiveLower;
    std::vector<TransferPacket> receiveUpper;
    exchangePackets(sendUpper, halo_->upperRank(), halo_->lowerRank(),
      TRANSFER_UPPER_COUNT, TRANSFER_UPPER_DATA,
      communicator_, receiveLower);
    exchangePackets(sendLower, halo_->lowerRank(), halo_->upperRank(),
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
          processActiveSegment(
            particle, segmentStart, transfer.segmentStartFraction, false);
        }

    std::vector<CarrierTransferPacket> receiveCarrierLower;
    std::vector<CarrierTransferPacket> receiveCarrierUpper;
    exchangePackets(sendCarrierUpper,
      halo_->upperRank(), halo_->lowerRank(),
      CARRIER_UPPER_COUNT, CARRIER_UPPER_DATA,
      communicator_, receiveCarrierLower);
    exchangePackets(sendCarrierLower,
      halo_->lowerRank(), halo_->upperRank(),
      CARRIER_LOWER_COUNT, CARRIER_LOWER_DATA,
      communicator_, receiveCarrierUpper);
    const std::vector<CarrierTransferPacket>* receivedCarriers[2] = {
      &receiveCarrierLower, &receiveCarrierUpper
    };
    for (unsigned int side = 0; side < 2; ++side)
      for (std::size_t index = 0;
           index < receivedCarriers[side]->size(); ++index)
        {
          const CarrierTransferPacket& transfer =
            (*receivedCarriers[side])[index];
          ParticleCPMLCarrier carrier = unpackCarrier(transfer);
          FieldVector<Double> segmentStart(0.0);
          for (unsigned int axis = 0; axis < 3; ++axis)
            segmentStart[axis] = transfer.segmentStart[axis];
          processCarrierSegment(
            carrier, segmentStart, transfer.segmentStartFraction, false);
        }

    std::vector<CarrierTransferPacket> receiveRetirementLower;
    std::vector<CarrierTransferPacket> receiveRetirementUpper;
    exchangePackets(sendRetirementUpper,
      halo_->upperRank(), halo_->lowerRank(),
      RETIREMENT_UPPER_COUNT, RETIREMENT_UPPER_DATA,
      communicator_, receiveRetirementLower);
    exchangePackets(sendRetirementLower,
      halo_->lowerRank(), halo_->upperRank(),
      RETIREMENT_LOWER_COUNT, RETIREMENT_LOWER_DATA,
      communicator_, receiveRetirementUpper);
    const std::vector<CarrierTransferPacket>* receivedRetirement[2] = {
      &receiveRetirementLower, &receiveRetirementUpper
    };
    for (unsigned int side = 0; side < 2; ++side)
      for (std::size_t index = 0;
           index < receivedRetirement[side]->size(); ++index)
        {
          const CarrierTransferPacket& transfer =
            (*receivedRetirement[side])[index];
          ParticleCPMLCarrier carrier = unpackCarrier(transfer);
          FieldVector<Double> segmentStart(0.0);
          for (unsigned int axis = 0; axis < 3; ++axis)
            segmentStart[axis] = transfer.segmentStart[axis];
          processRetirementSegment(
            carrier, segmentStart, transfer.segmentStartFraction, false);
        }
    particles_.swap(retained);
    pmlCarriers_.swap(retainedCarriers);
    retirementCarriers_.swap(retainedRetirementCarriers);
    peakPmlCarriers_ = std::max(peakPmlCarriers_, pmlCarriers_.size());
    peakRetirementCarriers_ = std::max(
      peakRetirementCarriers_, retirementCarriers_.size());
    unsigned long long globalLost = 0;
    MPI_Allreduce(&localLost, &globalLost, 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, communicator_);
    lostParticles_ += globalLost;
    if (detectors_) detectors_->collectParticleCrossings();
  }

  void Simulation::reportParticleBoundaryLosses() const
  {
    unsigned long long localCounts[18] = {};
    Double localCharges[18] = {};
    for (std::size_t face = 0; face < 6; ++face)
      {
        localCounts[face] = cpmlEntryCount_[face];
        localCharges[face] = cpmlEntryCharge_[face];
        localCounts[6 + face] = directOuterCount_[face];
        localCharges[6 + face] = directOuterCharge_[face];
        localCounts[12 + face] = carrierOuterCount_[face];
        localCharges[12 + face] = carrierOuterCharge_[face];
      }
    unsigned long long globalCounts[18] = {};
    Double globalCharges[18] = {};
    MPI_Reduce(localCounts, globalCounts, 18, MPI_UNSIGNED_LONG_LONG,
               MPI_SUM, 0, communicator_);
    MPI_Reduce(localCharges, globalCharges, 18, MPI_DOUBLE,
               MPI_SUM, 0, communicator_);
    const unsigned long long localPeak =
      static_cast<unsigned long long>(peakPmlCarriers_);
    unsigned long long maximumRankPeak = 0;
    MPI_Reduce(&localPeak, &maximumRankPeak, 1, MPI_UNSIGNED_LONG_LONG,
               MPI_MAX, 0, communicator_);
    const unsigned long long localRetirement[3] = {
      retirementEntryCount_, retirementExitCount_,
      static_cast<unsigned long long>(peakRetirementCarriers_)
    };
    unsigned long long globalRetirement[2] = {};
    unsigned long long maximumRetirementPeak = 0;
    MPI_Reduce(localRetirement, globalRetirement, 2,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator_);
    MPI_Reduce(&localRetirement[2], &maximumRetirementPeak, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, communicator_);
    const Double localRetirementCharge[2] = {
      retirementEntryCharge_, retirementExitResidualCharge_
    };
    Double globalRetirementCharge[2] = {};
    MPI_Reduce(localRetirementCharge, globalRetirementCharge, 2,
               MPI_DOUBLE, MPI_SUM, 0, communicator_);
    if (rank_ != 0) return;

    const char* labels[3] = {
      "CPML physical-domain exits",
      "Direct outer-domain exits",
      "Residual CPML-carrier cleanup"
    };
    for (std::size_t category = 0; category < 3; ++category)
      {
        unsigned long long total = 0;
        for (std::size_t face = 0; face < 6; ++face)
          total += globalCounts[6 * category + face];
        if (total == 0) continue;
        std::ostringstream message;
        message << std::setprecision(10) << labels[category] << ":";
        for (std::size_t face = 0; face < 6; ++face)
          {
            const ParticleBoundaryFace boundaryFace =
              static_cast<ParticleBoundaryFace>(face);
            const std::size_t offset = 6 * category + face;
            message << " " << particleBoundaryFaceName(boundaryFace)
                    << "=" << globalCounts[offset]
                    << " (Q=" << globalCharges[offset] << " C)";
            if (face + 1 < 6) message << ",";
          }
        logRoot(communicator_, message.str());
      }
    if (maximumRankPeak > 0)
      {
        std::ostringstream message;
        message << "Peak lightweight CPML carriers on any MPI rank="
                << maximumRankPeak << " ("
                << maximumRankPeak * sizeof(ParticleCPMLCarrier)
                << " bytes of live carrier records).";
        logRoot(communicator_, message.str());
      }
    if (globalRetirement[0] > 0 || globalRetirement[1] > 0)
      {
        std::ostringstream message;
        message << std::setprecision(10)
                << "Experimental particle retirement: entries="
                << globalRetirement[0] << " (Q="
                << globalRetirementCharge[0] << " C), completed="
                << globalRetirement[1] << ", residual Q at taper exit="
                << globalRetirementCharge[1]
                << " C, maximum live carriers on any rank="
                << maximumRetirementPeak << " ("
                << maximumRetirementPeak * sizeof(ParticleCPMLCarrier)
                << " bytes).";
        logRoot(communicator_, message.str());
      }
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
         element < beamlineElements_.size(); ++element)
      entrance = std::min(entrance,
        beamlineElements_[element].physicalEntrance);
    return entrance;
  }

  Double Simulation::firstBeamlineInteractionEntranceLab() const
  {
    Double entrance = std::numeric_limits<Double>::infinity();
    for (std::size_t element = 0;
         element < beamlineElements_.size(); ++element)
      entrance = std::min(entrance,
        beamlineElements_[element].interactionEntrance);
    return entrance;
  }

  Double Simulation::lastBeamlineInteractionExitLab() const
  {
    Double exit = -std::numeric_limits<Double>::infinity();
    for (std::size_t element = 0;
         element < beamlineElements_.size(); ++element)
      exit = std::max(exit,
        beamlineElements_[element].interactionExit);
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
    int localInCpml = 0;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      for (unsigned int axis = 0; axis < 3; ++axis)
        if (particles_[index].position[axis] < globalOriginBox_[axis] ||
            particles_[index].position[axis] > upper[axis])
          localInvalid = 1;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      if (!particleCPML_->containsPhysical(particles_[index].position))
        localInCpml = 1;
    int globalInvalid = 0;
    int globalInCpml = 0;
    MPI_Allreduce(&localInvalid, &globalInvalid, 1, MPI_INT, MPI_MAX,
                  communicator_);
    MPI_Allreduce(&localInCpml, &globalInCpml, 1, MPI_INT, MPI_MAX,
                  communicator_);
    if (globalInvalid)
      throw std::out_of_range(
        "Boosted bunch does not fit inside the configured direct E/B box");
    if (globalInCpml)
      throw std::out_of_range(
        "Boosted bunch overlaps the non-physical CPML particle region; enlarge the mesh or reduce the CPML thickness");
  }
}
