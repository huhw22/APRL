#include "simulation.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>

#include "particle_initializer.h"
#include "runtime_control.h"
#include "runtime_util.h"

namespace aprl
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
      runMetadata_(),
      globalGeometry_(), localGeometry_(), localZOffset_(0),
      globalOriginBox_(0.0), localOriginBox_(0.0),
      referenceCenterBoxZ_(0.0), frame_(), sources_(),
      fields_(), halo_(), incident_(), particles_(), particleCPML_(),
      pmlCarriers_(), retirementCarriers_(), particleBoundary_(), detectors_(),
      trajectoryWriter_(), energyLedgerWriter_(), lastEnergyLedgerRecord_(),
      trajectoryRhythmSI_(0.0), nextTrajectorySampleTime_(0.0),
      trajectorySamplesSinceFlush_(0), particleSubsteps_(1),
      energyLedgerPrescribedWork_(0.0L),
      energyLedgerRemovedKinetic_(0.0L),
      energyLedgerRemovedTotal_(0.0L),
      energyLedgerInitialKinetic_(0.0L),
      energyLedgerInitialField_(0.0L),
      energyLedgerRemovedMacroparticles_(0),
      energyLedgerLastSampleStep_(0), energyLedgerHasReference_(false),
      timeBoxSI_(0.0),
      totalTimeBoxSI_(0.0), step_(0), interrupted_(false),
      configuredStopReached_(false), loopWallStart_(0.0),
      estimatedStepSeconds_(0.0), estimatedMaximumSteps_(0),
      modeledLocalPeakBytes_(0), lostParticles_(0),
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
        energyLedgerPreviousPower_[face] = 0.0L;
        energyLedgerOutwardEnergy_[face] = 0.0L;
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
    if (config_.energyLedger.enabled)
      sampleEnergyLedger(true);
    loopWallStart_ = MPI_Wtime();

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

        if (config_.energyLedger.enabled)
          advanceEnergyLedgerStep();

        timeBoxSI_ += globalGeometry_.dt;
        ++step_;
        if (config_.trajectory.enabled &&
            timeBoxSI_ >= nextTrajectorySampleTime_)
          sampleTrajectory();
        if (detectors_)
          detectors_->sampleFieldPlanes(*fields_, timeBoxSI_);
        if (config_.energyLedger.enabled)
          sampleEnergyLedger();
        if (config_.runtime.resourceReport &&
            config_.runtime.resourceProgressIntervalSteps > 0 &&
            step_ % config_.runtime.resourceProgressIntervalSteps == 0)
          reportResourceProgress("progress");
        if (configuredStopReached())
          {
            configuredStopReached_ = true;
            break;
          }
      }

    if (!interrupted_) synchronizedStopRequested();
    const bool completed = !interrupted_ && configuredStopReached_;
    if (config_.energyLedger.enabled)
      sampleEnergyLedger(true);
    finalizeTrajectoryOutput(completed);
    if (detectors_) detectors_->close(completed);
    finalizeEnergyLedger(completed);
    reportParticleBoundaryLosses();
    reportResourceProgress("final");
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

    preflightOutputPaths();
    runMetadata_ = initializeRunMetadata(config_, communicator_);
    initializeGeometry();
    validateAndReportRadiationResolution();
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
    reportDurationEstimate();
    fields_.reset(new EBFieldGrid(localGeometry_));
    halo_.reset(new EBZSlabHaloExchange(communicator_));
    particleBoundary_.reset(new ParticleOpenBoundary(
      localGeometry_, localOriginBox_, globalGeometry_, globalOriginBox_,
      localZOffset_));
    if (config_.boundary.type == EBBoundaryType::Cpml)
      fields_->setBoundary(std::unique_ptr<EBBoundaryOperator>(
        new EBConvolutionalPML(
          localGeometry_, globalGeometry_.nz, localZOffset_,
          config_.boundary.cpml,
          halo_->lowerRank() == MPI_PROC_NULL,
          halo_->upperRank() == MPI_PROC_NULL)));
    else
      halo_->installPhysicalBoundaryMask(*fields_);
    calibrateResourceEstimate();
    fields_->clearFields();
    initializeParticleSelfField();

    if (sources_.maxwellIncidentWaveCount() > 0)
      {
        const EBTFSFRegion region(
          1, globalGeometry_.nx - 1,
          1, globalGeometry_.ny - 1,
          1, globalGeometry_.nz - 1);
        incident_.reset(new EBMaxwellIncidentInjector(
          localGeometry_, region, localZOffset_, globalGeometry_.nz));
        incident_->initialize(*fields_, sources_, localOriginBox_,
                              timeBoxSI_, frame_, false);
      }

    initializeTrajectoryOutput();
    initializeDetectorOutput();
    initializeEnergyLedger();
    reportResourceEstimate();
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
        const std::size_t slabBase = globalGeometry_.nz /
          static_cast<std::size_t>(size_);
        const std::size_t slabRemainder = globalGeometry_.nz %
          static_cast<std::size_t>(size_);
        std::ostringstream meshMessage;
        meshMessage << std::setprecision(10)
          << "Mesh contract: cells=(" << globalGeometry_.nx << ", "
          << globalGeometry_.ny << ", " << globalGeometry_.nz
          << "), cell_size_m=(" << globalGeometry_.dx << ", "
          << globalGeometry_.dy << ", " << globalGeometry_.dz
          << "), derived_extent_m=(" << config_.mesh.extent[0] << ", "
          << config_.mesh.extent[1] << ", " << config_.mesh.extent[2]
          << "); z slabs use integer offsets with cells/rank in ["
          << slabBase << ", "
          << slabBase + (slabRemainder > 0 ? 1 : 0)
          << "] and remainder=" << slabRemainder << ".";
        logRoot(communicator_, meshMessage.str());
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
        if (!config_.initialSelfField.enabled)
          logRoot(communicator_,
            "WARNING: initial_self_field is disabled; the Maxwell state does not satisfy Gauss's law for the input bunch and startup radiation can contaminate the result.");
      }
  }

  void Simulation::preflightOutputPaths() const
  {
    if (config_.output.overwrite) return;
    std::vector<std::string> paths;
    if (rank_ == 0)
      {
        paths.push_back(config_.output.manifest);
        for (std::size_t i = 0; i < config_.detectors.fieldPlanes.size(); ++i)
          {
            const FieldDetectorPlaneConfig& plane =
              config_.detectors.fieldPlanes[i];
            paths.push_back(joinPath(config_.detectors.directory,
              plane.name + ".h5"));
            if (plane.particleBackgroundReference)
              paths.push_back(joinPath(config_.detectors.directory,
                plane.name + "-ballistic-reference.h5"));
          }
        for (std::size_t i = 0;
             i < config_.detectors.particlePlanes.size(); ++i)
          paths.push_back(joinPath(config_.detectors.directory,
            config_.detectors.particlePlanes[i].name + ".h5"));
        if (config_.energyLedger.enabled)
          paths.push_back(joinPath(config_.energyLedger.directory,
            config_.energyLedger.filename));
      }
    if (config_.trajectory.enabled)
      {
        std::ostringstream name;
        name << (config_.trajectory.basename.empty() ?
          "trajectory" : config_.trajectory.basename)
             << "-rank-" << std::setfill('0') << std::setw(5)
             << rank_ << ".h5";
        paths.push_back(joinPath(config_.trajectory.directory, name.str()));
      }

    std::string conflict;
    for (std::size_t i = 0; i < paths.size(); ++i)
      {
        struct stat status;
        if (::stat(paths[i].c_str(), &status) == 0)
          {
            conflict = paths[i];
            break;
          }
      }
    int localConflict = conflict.empty() ? 0 : 1;
    int anyConflict = 0;
    MPI_Allreduce(&localConflict, &anyConflict, 1, MPI_INT, MPI_MAX,
                  communicator_);
    if (anyConflict)
      {
        if (!conflict.empty())
          throw std::runtime_error("Refusing to overwrite existing output: " +
            conflict +
            "; set output.overwrite: true only for an intentional rerun");
        throw std::runtime_error(
          "Another MPI rank found an existing output; refusing the run. "
          "Set output.overwrite: true only for an intentional rerun");
      }
  }

  void Simulation::validateAndReportRadiationResolution() const
  {
    const RadiationResolutionConfig& target = config_.radiationResolution;
    if (!target.enabled)
      {
        if (rank_ == 0)
          logRoot(communicator_,
            "WARNING: radiation_resolution is disabled; no target-band "
            "Nyquist check is being made. This does not affect generic "
            "particle/field runs, but production radiation cards should "
            "declare maximum_photon_energy_eV.");
        return;
      }

    const Double photonEnergyLab = target.maximumPhotonEnergyEV;
    const Double gamma = config_.mesh.boostGamma;
    const Double rapidityFactor = gamma +
      BoostFrameTransform::gammaBetaFromGamma(gamma);
    const Double onAxisDoppler = 1.0 / rapidityFactor;
    const Double periodLab = PLANCK_EV_SECOND / photonEnergyLab;
    const Double wavelengthLab = SI::c * periodLab;
    const Double periodBox = periodLab / onAxisDoppler;
    const Double wavelengthBox = wavelengthLab / onAxisDoppler;
    const Double gridPoints = wavelengthBox / globalGeometry_.dz;
    const Double maxwellSamples = periodBox / globalGeometry_.dt;
    const Double nyquistTolerance = 128.0 *
      std::numeric_limits<Double>::epsilon();

    if (!(gridPoints > 2.0 * (1.0 + nyquistTolerance)) ||
        !(maxwellSamples > 2.0 * (1.0 + nyquistTolerance)))
      {
        std::ostringstream message;
        message << std::setprecision(10)
          << "Radiation target at " << photonEnergyLab
          << " eV is below the strict on-axis +z Nyquist resolution: "
          << "lab wavelength=" << wavelengthLab
          << " m, boosted wavelength=" << wavelengthBox
          << " m, z points/wavelength=" << gridPoints
          << ", Maxwell samples/boosted cycle=" << maxwellSamples
          << ". Both must exceed 2. Refine mesh.cell_size z below "
          << 0.5 * wavelengthBox << " m and/or reduce the Maxwell dt below "
          << 0.5 * periodBox
          << " s. Cowan-z uses dt=dz/c, so refining dz satisfies both. "
             "This guard only covers a forward on-axis mode.";
        throw std::runtime_error(message.str());
      }

    if (rank_ == 0)
      {
        std::ostringstream message;
        message << std::setprecision(10)
          << "Radiation-resolution target (forward on-axis lab mode): "
          << "maximum photon energy=" << photonEnergyLab
          << " eV, lab wavelength=" << wavelengthLab
          << " m, boosted photon energy="
          << photonEnergyLab * onAxisDoppler
          << " eV, boosted wavelength=" << wavelengthBox
          << " m; z points/wavelength=" << gridPoints
          << ", Maxwell/current samples per boosted cycle="
          << maxwellSamples << ". Strict Nyquist passed.";
        logRoot(communicator_, message.str());
        if (gridPoints < target.warningGridPointsPerWavelength ||
            maxwellSamples < target.warningMaxwellSamplesPerCycle)
          {
            std::ostringstream warning;
            warning << std::setprecision(8)
              << "WARNING: target-band grid/current sampling is below a "
                 "configured quality level (z=" << gridPoints << " versus "
              << target.warningGridPointsPerWavelength
              << " points/wavelength, time=" << maxwellSamples << " versus "
              << target.warningMaxwellSamplesPerCycle
              << " Maxwell samples/cycle"
              << "). This is not a hard validity limit: Cowan-z has exact "
                 "resolved on-axis vacuum phase velocity, but source "
                 "deposition, interpolation and amplitude still require a "
                 "problem-specific refinement scan.";
            logRoot(communicator_, warning.str());
          }
      }

    const Double detectorFieldStepLab = globalGeometry_.dt / gamma;
    for (std::size_t index = 0;
         index < config_.detectors.fieldPlanes.size(); ++index)
      {
        const FieldDetectorPlaneConfig& plane =
          config_.detectors.fieldPlanes[index];
        const Double requestedSteps = plane.rhythm / detectorFieldStepLab;
        const Double roundingTolerance = 128.0 *
          std::numeric_limits<Double>::epsilon() *
          std::max(1.0, std::abs(requestedSteps));
        const Double fieldSteps = std::max(
          1.0, std::ceil(requestedSteps - roundingTolerance));
        const Double maximumRealizedGap =
          fieldSteps * detectorFieldStepLab;
        const Double detectorSamples = periodLab / maximumRealizedGap;
        if (!(detectorSamples > 2.0 * (1.0 + nyquistTolerance)))
          {
            std::ostringstream message;
            message << std::setprecision(10)
              << "Field detector '" << plane.name << "' cannot Nyquist-sample "
              << photonEnergyLab << " eV: requested rhythm=" << plane.rhythm
              << " s, one Maxwell step at the fixed lab plane="
              << detectorFieldStepLab
              << " s, conservative realized gap=" << maximumRealizedGap
              << " s, samples/cycle=" << detectorSamples
              << ". More than 2 are required.";
            if (!(detectorFieldStepLab < 0.5 * periodLab))
              message << " Refine the Maxwell step below "
                      << 0.5 * gamma * periodLab
                      << " s in the boosted frame";
            else
              {
                const Double maximumSteps = std::max(1.0, std::floor(
                  0.5 * periodLab / detectorFieldStepLab -
                  nyquistTolerance));
                message << " Set detector rhythm to at most "
                        << maximumSteps * detectorFieldStepLab
                        << " s (an integer number of fixed-plane Maxwell "
                           "increments below half a target cycle)";
              }
            message << ". No higher-order temporal reconstruction is assumed.";
            throw std::runtime_error(message.str());
          }
        if (rank_ == 0)
          {
            std::ostringstream message;
            message << std::setprecision(10)
              << "Field detector '" << plane.name
              << "' target-band cadence: requested rhythm=" << plane.rhythm
              << " s, fixed-plane Maxwell increment="
              << detectorFieldStepLab
              << " s, conservative realized gap=" << maximumRealizedGap
              << " s, samples/cycle=" << detectorSamples
              << ". Strict Nyquist passed.";
            logRoot(communicator_, message.str());
            if (detectorSamples <
                target.warningDetectorSamplesPerCycle)
              {
                std::ostringstream warning;
                warning << std::setprecision(8)
                  << "WARNING: field detector '" << plane.name
                  << "' has fewer than "
                  << target.warningDetectorSamplesPerCycle
                  << " samples per target cycle. Spectral amplitude and "
                     "phase should be checked by reducing both detector "
                     "rhythm and, when necessary, the Maxwell step.";
                logRoot(communicator_, warning.str());
              }
          }
      }

    if (rank_ == 0)
      logRoot(communicator_,
        "Radiation-resolution scope: no optical-wavelength condition is "
        "imposed on dx/dy for an on-axis paraxial carrier. Transverse "
        "envelope, angular aperture, macro-particle noise, CPML reflection "
        "and spectral-window convergence remain problem-dependent checks.");
  }

  void Simulation::initializeParticleSelfField()
  {
    if (!config_.initialSelfField.enabled) return;
    std::size_t staticFieldGuardCells[3] = {0, 0, 0};
    if (config_.boundary.type == EBBoundaryType::Cpml)
      for (unsigned int axis = 0; axis < 3; ++axis)
        staticFieldGuardCells[axis] = config_.boundary.cpml.cells[axis];
    const EBGaussInitializationReport report =
      EBGaussFieldInitializer::initialize(
        *fields_, globalGeometry_, localZOffset_, localOriginBox_, particles_,
        config_.initialSelfField.model,
        staticFieldGuardCells,
        config_.initialSelfField.relativeTolerance,
        config_.initialSelfField.maximumIterations, communicator_);
    if (rank_ == 0)
      {
        std::ostringstream message;
        message << std::setprecision(10)
          << "Initial particle self-field: "
          << (config_.initialSelfField.model ==
                InitialSelfFieldModel::RelativisticPoisson ?
                "relativistic-Poisson rigid-beam" :
                "legacy electrostatic-Poisson")
          << " construction "
          << "converged in " << report.iterations
          << " CG iterations; post-check relative L2 residual="
          << report.relativeResidual
          << ", maximum absolute residual="
          << report.maximumGaussResidual << " V/m^2, total charge="
          << report.totalCharge << " C, mean beta_box_z="
          << report.meanBetaZ << ", mean beta_box_transverse="
          << report.meanBetaTransverse << ", initial electric energy="
          << report.initialElectricEnergy << " J, initial magnetic energy="
          << report.initialMagneticEnergy
          << " J, static zero-potential guard cells=("
          << report.staticFieldGuardCells[0] << ", "
          << report.staticFieldGuardCells[1] << ", "
          << report.staticFieldGuardCells[2]
          << "), bunch RMS cells=("
          << report.rmsPositionCells[0] << ", "
          << report.rmsPositionCells[1] << ", "
          << report.rmsPositionCells[2] << "), centre-to-static-boundary "
             "padding/RMS=("
          << report.centrePaddingRms[0] << ", "
          << report.centrePaddingRms[1] << ", "
          << report.centrePaddingRms[2]
          << "), temporary scalar memory/rank<="
          << static_cast<Double>(report.temporaryBytes) /
             (1024.0 * 1024.0) << " MiB. The potential is now released; "
             "only SI E/B remains. The reported static zero-potential "
             "boundary still requires padding convergence. ";
        if (config_.initialSelfField.model ==
            InitialSelfFieldModel::RelativisticPoisson)
          message << "The relativistic model represents a rigid "
                     "common-velocity bunch; velocity-spread and "
                     "kinetic-equilibrium errors require separate "
                     "convergence tests.";
        else
          message << "The legacy model deliberately omits the moving-bunch "
                     "magnetic field and is retained only for regression.";
        logRoot(communicator_, message.str());
        for (unsigned int axis = 0; axis < 3; ++axis)
          {
            if (report.rmsPositionCells[axis] < 2.0)
              {
                std::ostringstream warning;
                warning << "WARNING: initial bunch RMS on axis " << axis
                        << " is only " << report.rmsPositionCells[axis]
                        << " cells. A zero-width or under-resolved CIC "
                           "distribution can create macro-particle self-field "
                           "relaxation. Resolve each finite bunch RMS with at "
                           "least two cells (preferably four or more), or "
                           "treat the run as a sheet/point regression rather "
                           "than a collective-field result.";
                logRoot(communicator_, warning.str());
              }
            if (report.centrePaddingRms[axis] > 0.0 &&
                report.centrePaddingRms[axis] < 4.0)
              {
                std::ostringstream warning;
                warning << "WARNING: initial static-field boundary on axis "
                        << axis << " is only "
                        << report.centrePaddingRms[axis]
                        << " bunch RMS from the centroid. Enlarge mesh cells "
                           "at fixed cell_size and require the entrance "
                           "energy result to converge with padding.";
                logRoot(communicator_, warning.str());
              }
          }
      }
  }

  void Simulation::calibrateResourceEstimate()
  {
    const Double stepsReal = std::ceil(totalTimeBoxSI_ /
                                       globalGeometry_.dt);
    if (!(stepsReal > 0.0) || !std::isfinite(stepsReal) ||
        stepsReal > static_cast<Double>(
          std::numeric_limits<std::size_t>::max()))
      throw std::overflow_error(
        "Configured duration produces an invalid field-step count");
    estimatedMaximumSteps_ = static_cast<std::size_t>(stepsReal);
    if (!config_.runtime.resourceReport) return;

    if (sources_.maxwellIncidentWaveCount() > 0)
      {
        estimatedStepSeconds_ = 0.0;
        if (rank_ == 0)
          logRoot(communicator_,
            "[resource] phase=pre-run-calibration status=skipped reason=incident-wave-state");
        return;
      }

    fields_->clearCurrent();
    MPI_Barrier(communicator_);
    const Double fieldStart = MPI_Wtime();
    for (unsigned int repetition = 0;
         repetition < config_.runtime.resourceCalibrationSteps;
         ++repetition)
      {
        halo_->advanceMagnetic(*fields_);
        halo_->advanceElectric(*fields_);
        fields_->clearCurrent();
      }
    MPI_Barrier(communicator_);
    const Double localFieldStep =
      (MPI_Wtime() - fieldStart) /
      static_cast<Double>(config_.runtime.resourceCalibrationSteps);
    Double fieldStep = 0.0;
    MPI_Allreduce(&localFieldStep, &fieldStep, 1, MPI_DOUBLE,
                  MPI_MAX, communicator_);

    const std::size_t maximumSamples = 4096;
    const std::size_t sampleCount =
      std::min(maximumSamples, particles_.size());
    const std::size_t stride = sampleCount > 0 ?
      (particles_.size() + sampleCount - 1) / sampleCount : 1;
    ChargeConservingCurrentDepositor depositor(*fields_,
                                               localOriginBox_);
    std::size_t completedSamples = 0;
    MPI_Barrier(communicator_);
    const Double particleStart = MPI_Wtime();
    for (unsigned int repetition = 0;
         repetition < config_.runtime.resourceCalibrationSteps;
         ++repetition)
      for (std::size_t index = 0;
           index < particles_.size() && completedSamples <
             sampleCount * static_cast<std::size_t>(
               repetition + 1);
           index += stride)
        {
          RelativisticParticleSI particle(particles_[index]);
          const FieldVector<Double> start(particle.position);
          RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
            particle, *fields_, localOriginBox_, sources_, frame_,
            timeBoxSI_, globalGeometry_.dt, particleSubsteps_);
          bool local = true;
          for (unsigned int axis = 0; axis < 3; ++axis)
            {
              const Double upper = localOriginBox_[axis] +
                static_cast<Double>(axis == 0 ? localGeometry_.nx :
                  (axis == 1 ? localGeometry_.ny : localGeometry_.nz)) *
                (axis == 0 ? localGeometry_.dx :
                  (axis == 1 ? localGeometry_.dy : localGeometry_.dz));
              local = local &&
                particle.position[axis] >= localOriginBox_[axis] &&
                particle.position[axis] <= upper;
            }
          if (local)
            depositor.depositSegment(start, particle.position,
                                     particle.charge);
          ++completedSamples;
        }
    const Double localParticleElapsed = MPI_Wtime() - particleStart;
    fields_->clearCurrent();
    const Double localParticleStep = completedSamples > 0 ?
      localParticleElapsed /
        static_cast<Double>(completedSamples) *
        static_cast<Double>(particles_.size()) : 0.0;
    Double particleStep = 0.0;
    MPI_Allreduce(&localParticleStep, &particleStep, 1, MPI_DOUBLE,
                  MPI_MAX, communicator_);

    estimatedStepSeconds_ =
      (fieldStep + particleStep) * config_.runtime.timeSafetyFactor;
    if (!(estimatedStepSeconds_ > 0.0) ||
        !std::isfinite(estimatedStepSeconds_))
      estimatedStepSeconds_ = 0.0;
  }

  void Simulation::reportResourceEstimate()
  {
    if (!config_.runtime.resourceReport) return;

    const EBMemoryFootprint field = fields_->memoryFootprint();
    long double modeled = static_cast<long double>(field.totalBytes()) +
      static_cast<long double>(halo_->memoryBytes()) +
      static_cast<long double>(particleBoundary_->memoryBytes()) +
      static_cast<long double>(trajectoryWriter_.memoryBytes()) +
      static_cast<long double>(energyLedgerWriter_.memoryBytes()) +
      static_cast<long double>(particles_.capacity()) *
        sizeof(RelativisticParticleSI) +
      static_cast<long double>(pmlCarriers_.capacity()) *
        sizeof(ParticleCPMLCarrier) +
      static_cast<long double>(retirementCarriers_.capacity()) *
        sizeof(ParticleCPMLCarrier);

    /* During migration the old particle vector, retained vector and transfer
     * packets coexist.  This bounded estimate deliberately covers three
     * additional live records per initially local particle. */
    modeled += 3.0L * static_cast<long double>(particles_.size()) *
      static_cast<long double>(std::max(
        sizeof(RelativisticParticleSI), sizeof(TransferPacket)));

    const long double fieldPointBytes = 6.0L * sizeof(Double);
    const long double planePoints =
      static_cast<long double>(globalGeometry_.nx) *
      static_cast<long double>(globalGeometry_.ny);
    if (!detectorConfig_.fieldPlanes.empty())
      modeled += planePoints * fieldPointBytes;
    if (rank_ == 0)
      {
        const long double particlePlaneRecordBytes = 112.0L;
        const long double validationEntryBytes = 192.0L;
        for (std::size_t detector = 0;
             detector < detectorConfig_.fieldPlanes.size(); ++detector)
          {
            const FieldDetectorPlaneConfig& plane =
              detectorConfig_.fieldPlanes[detector];
            modeled += static_cast<long double>(plane.bufferSamples) *
              (sizeof(Double) + planePoints * fieldPointBytes);
            if (plane.particleBackgroundReference)
              modeled += static_cast<long double>(
                plane.referenceBufferRecords) *
                particlePlaneRecordBytes;
            if (plane.referenceValidation)
              modeled += static_cast<long double>(
                plane.referenceValidationMaximumParticles) *
                validationEntryBytes;
          }
        for (std::size_t detector = 0;
             detector < detectorConfig_.particlePlanes.size(); ++detector)
          modeled += static_cast<long double>(
            detectorConfig_.particlePlanes[detector].bufferRecords) *
            particlePlaneRecordBytes;
      }

    /* The matrix-free initial Poisson solve owns four vertex slabs only
     * during initialization. Include that transient allocation in the
     * conservative peak even though detector buffers are opened later. */
    if (config_.initialSelfField.enabled)
      modeled += 4.0L * sizeof(Double) *
        static_cast<long double>(localGeometry_.nx + 1) *
        static_cast<long double>(localGeometry_.ny + 1) *
        static_cast<long double>(localGeometry_.nz + 2);

    modeled *= config_.runtime.memorySafetyFactor;
    const long double byteMaximum =
      static_cast<long double>(
        std::numeric_limits<std::uint64_t>::max());
    modeledLocalPeakBytes_ = static_cast<std::uint64_t>(
      std::min(modeled, byteMaximum));
    const std::uint64_t resident = currentResidentBytes();
    const long double residentWithMargin =
      static_cast<long double>(resident) *
      config_.runtime.memorySafetyFactor;
    if (residentWithMargin >
        static_cast<long double>(modeledLocalPeakBytes_))
      modeledLocalPeakBytes_ = static_cast<std::uint64_t>(
        std::min(residentWithMargin, byteMaximum));

    const unsigned long long localModeled =
      static_cast<unsigned long long>(modeledLocalPeakBytes_);
    const unsigned long long localResident =
      static_cast<unsigned long long>(resident);
    unsigned long long totalModeled = 0;
    unsigned long long maximumModeled = 0;
    unsigned long long totalResident = 0;
    unsigned long long maximumResident = 0;
    MPI_Reduce(&localModeled, &totalModeled, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator_);
    MPI_Reduce(&localModeled, &maximumModeled, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, communicator_);
    MPI_Reduce(&localResident, &totalResident, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator_);
    MPI_Reduce(&localResident, &maximumResident, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, communicator_);

    const unsigned long long localParticles =
      static_cast<unsigned long long>(particles_.size());
    unsigned long long globalParticles = 0;
    MPI_Reduce(&localParticles, &globalParticles, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator_);

    long double rawFieldOutputBytes = 0.0L;
    for (std::size_t detector = 0;
         detector < detectorConfig_.fieldPlanes.size(); ++detector)
      {
        const Double labDuration =
          totalTimeBoxSI_ / frame_.gamma();
        const long double samples = std::ceil(
          labDuration / detectorConfig_.fieldPlanes[detector].rhythm) +
          2.0L;
        rawFieldOutputBytes += samples * planePoints *
          fieldPointBytes;
      }
    long double rawTrajectoryBytes = 0.0L;
    if (config_.trajectory.enabled && trajectoryRhythmSI_ > 0.0)
      {
        const long double samples = std::ceil(
          totalTimeBoxSI_ / trajectoryRhythmSI_) + 2.0L;
        rawTrajectoryBytes = samples *
          static_cast<long double>(globalParticles) *
          sizeof(TrajectoryRecord);
      }

    if (rank_ == 0)
      {
        const Double mib = 1024.0 * 1024.0;
        const Double gib = 1024.0 * mib;
        const Double estimatedWall =
          estimatedStepSeconds_ > 0.0 ?
          estimatedStepSeconds_ *
            static_cast<Double>(estimatedMaximumSteps_) : 0.0;
        std::ostringstream work;
        work << std::setprecision(10)
          << "[resource] phase=pre-run ranks=" << size_
          << " grid_cells=" << globalGeometry_.nx << "x"
          << globalGeometry_.ny << "x" << globalGeometry_.nz
          << " particles=" << globalParticles
          << " dt_s=" << globalGeometry_.dt
          << " steps_duration_upper=" << estimatedMaximumSteps_
          << " calibrated_step_s=" << estimatedStepSeconds_
          << " estimated_wall_s=" << estimatedWall
          << " time_safety_factor=" << config_.runtime.timeSafetyFactor;
        logRoot(communicator_, work.str());

        std::ostringstream memory;
        memory << std::setprecision(10)
          << "[resource] phase=pre-run estimated_peak_rank_max_mib="
          << static_cast<Double>(maximumModeled) / mib
          << " estimated_peak_total_mib="
          << static_cast<Double>(totalModeled) / mib
          << " current_rss_rank_max_mib="
          << static_cast<Double>(maximumResident) / mib
          << " current_rss_total_mib="
          << static_cast<Double>(totalResident) / mib
          << " memory_safety_factor="
          << config_.runtime.memorySafetyFactor
          << " uncompressed_field_output_gib="
          << static_cast<Double>(rawFieldOutputBytes) / gib
          << " uncompressed_trajectory_output_gib="
          << static_cast<Double>(rawTrajectoryBytes) / gib;
        logRoot(communicator_, memory.str());
        logRoot(communicator_,
          "[resource] note=wall estimate is a zero-field kernel microbenchmark plus initial particle push/deposition estimate; detector I/O, filesystem contention, particle migration and early physical stopping can change it. Calibrate safety factors per machine and scale. Lines are root-only stdout and remain intact under sbatch redirection.");
      }
  }

  void Simulation::reportResourceProgress(const char* phase)
  {
    if (!config_.runtime.resourceReport) return;
    const unsigned long long localCurrent =
      static_cast<unsigned long long>(currentResidentBytes());
    const unsigned long long localPeak = std::max(
      localCurrent,
      static_cast<unsigned long long>(peakResidentBytes()));
    unsigned long long totalCurrent = 0;
    unsigned long long maximumCurrent = 0;
    unsigned long long totalPeak = 0;
    unsigned long long maximumPeak = 0;
    MPI_Reduce(&localCurrent, &totalCurrent, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator_);
    MPI_Reduce(&localCurrent, &maximumCurrent, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, communicator_);
    MPI_Reduce(&localPeak, &totalPeak, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator_);
    MPI_Reduce(&localPeak, &maximumPeak, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, communicator_);
    if (rank_ != 0) return;

    const Double elapsed = loopWallStart_ > 0.0 ?
      MPI_Wtime() - loopWallStart_ : 0.0;
    const Double secondsPerStep = step_ > 0 ?
      elapsed / static_cast<Double>(step_) : 0.0;
    const std::size_t remainingSteps =
      estimatedMaximumSteps_ > step_ ?
      estimatedMaximumSteps_ - step_ : 0;
    const Double projectedRemaining =
      secondsPerStep * static_cast<Double>(remainingSteps);
    const Double mib = 1024.0 * 1024.0;
    std::ostringstream message;
    message << std::setprecision(10)
      << "[resource] phase=" << (phase ? phase : "unknown")
      << " step=" << step_
      << " elapsed_s=" << elapsed
      << " measured_step_s=" << secondsPerStep
      << " duration_upper_remaining_s=" << projectedRemaining
      << " current_rss_rank_max_mib="
      << static_cast<Double>(maximumCurrent) / mib
      << " current_rss_total_mib="
      << static_cast<Double>(totalCurrent) / mib
      << " peak_rss_rank_max_mib="
      << static_cast<Double>(maximumPeak) / mib
      << " peak_rss_total_mib="
      << static_cast<Double>(totalPeak) / mib;
    logRoot(communicator_, message.str());
  }


  void Simulation::initializeGeometry()
  {
    const Double dx = config_.mesh.cellSize[0];
    const Double dy = config_.mesh.cellSize[1];
    const Double dz = config_.mesh.cellSize[2];
    const std::size_t nx = config_.mesh.cells[0];
    const std::size_t ny = config_.mesh.cells[1];
    const std::size_t nz = config_.mesh.cells[2];
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
            message << std::setprecision(10)
              << "Cowan-z requires z to have the smallest cell_size so "
                 "c*dt=dz is stable and dispersion-free.";
            if (dx + tolerance < dz)
              message << " Set mesh.cell_size x to at least " << dz
                      << " m (" << dz / config_.inputUnits.length
                      << " in the configured length unit).";
            if (dy + tolerance < dz)
              message << " Set mesh.cell_size y to at least " << dz
                      << " m (" << dz / config_.inputUnits.length
                      << " in the configured length unit).";
            message << " Keep mesh.cells as explicit integers; the physical "
                       "extent will be recomputed by multiplication.";
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
      globalOriginBox_[axis] = std::fma(
        -0.5, config_.mesh.extent[axis], config_.mesh.center[axis]);
    localOriginBox_ = globalOriginBox_;
    localOriginBox_[2] = std::fma(
      static_cast<Double>(localZOffset_), dz, globalOriginBox_[2]);
    totalTimeBoxSI_ = config_.mesh.duration;
    if (!(totalTimeBoxSI_ > 0.0) || !std::isfinite(totalTimeBoxSI_))
      throw std::invalid_argument(
        "Direct E/B total-time must be a positive box-frame duration");
  }

  void Simulation::initializeParticles()
  {
    ParticleInitializationReport inputReport;
    particles_ = ParticleInitializer::create(
      config_, communicator_, inputReport);
    const unsigned long long localCount =
      static_cast<unsigned long long>(particles_.size());
    unsigned long long idOffset = 0;
    MPI_Exscan(&localCount, &idOffset, 1, MPI_UNSIGNED_LONG_LONG,
               MPI_SUM, communicator_);
    if (rank_ == 0) idOffset = 0;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      particles_[index].id = idOffset + index + 1;


    unsigned long long globalCount = 0;
    MPI_Allreduce(&localCount, &globalCount, 1,
                  MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator_);
    long double localRepresentedElectrons = 0.0L;
    long double localMinimumMacroElectrons =
      std::numeric_limits<long double>::infinity();
    long double localMaximumMacroElectrons = 0.0L;
    Double localMinimumInputWeight =
      std::numeric_limits<Double>::infinity();
    Double localMaximumInputWeight = 0.0;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      {
        const RelativisticParticleSI& particle = particles_[index];
        const long double represented = std::abs(
          static_cast<long double>(particle.charge) /
          static_cast<long double>(SI::elementaryCharge));
        if (!(represented > 0.0L) ||
            !std::isfinite(static_cast<Double>(represented)) ||
            !(particle.mass > 0.0) || !std::isfinite(particle.mass) ||
            !(particle.weight > 0.0) || !std::isfinite(particle.weight))
          throw std::runtime_error(
            "Particle input produced an invalid macro charge, mass or weight");
        const Double chargeToMass = particle.charge / particle.mass;
        const Double electronChargeToMass =
          -SI::elementaryCharge / SI::electronMass;
        if (std::abs(chargeToMass / electronChargeToMass - 1.0) >
            256.0 * std::numeric_limits<Double>::epsilon())
          throw std::runtime_error(
            "Particle macro charge and mass do not preserve the electron charge-to-mass ratio");
        localRepresentedElectrons += represented;
        localMinimumMacroElectrons = std::min(
          localMinimumMacroElectrons, represented);
        localMaximumMacroElectrons = std::max(
          localMaximumMacroElectrons, represented);
        localMinimumInputWeight = std::min(
          localMinimumInputWeight, particle.weight);
        localMaximumInputWeight = std::max(
          localMaximumInputWeight, particle.weight);
      }
    long double globalRepresentedElectrons = 0.0L;
    long double globalMinimumMacroElectrons = 0.0L;
    long double globalMaximumMacroElectrons = 0.0L;
    Double globalMinimumInputWeight = 0.0;
    Double globalMaximumInputWeight = 0.0;
    MPI_Allreduce(&localRepresentedElectrons,
                  &globalRepresentedElectrons, 1,
                  MPI_LONG_DOUBLE, MPI_SUM, communicator_);
    MPI_Allreduce(&localMinimumMacroElectrons,
                  &globalMinimumMacroElectrons, 1,
                  MPI_LONG_DOUBLE, MPI_MIN, communicator_);
    MPI_Allreduce(&localMaximumMacroElectrons,
                  &globalMaximumMacroElectrons, 1,
                  MPI_LONG_DOUBLE, MPI_MAX, communicator_);
    MPI_Allreduce(&localMinimumInputWeight,
                  &globalMinimumInputWeight, 1,
                  MPI_DOUBLE, MPI_MIN, communicator_);
    MPI_Allreduce(&localMaximumInputWeight,
                  &globalMaximumInputWeight, 1,
                  MPI_DOUBLE, MPI_MAX, communicator_);
    const long double electronMismatch = std::abs(
      globalRepresentedElectrons -
      static_cast<long double>(config_.beam.electrons)) /
      static_cast<long double>(config_.beam.electrons);
    if (electronMismatch > 1.0e-12L)
      throw std::runtime_error(
        "Normalized macro-particle charges do not reproduce beam.input.electrons to 1e-12 relative accuracy");
    if (rank_ == 0)
      {
        std::ostringstream message;
        message << std::setprecision(10)
                << "Macro-particle normalization: records="
                << globalCount << ", represented electrons="
                << static_cast<Double>(globalRepresentedElectrons)
                << ", electrons per macro range=["
                << static_cast<Double>(globalMinimumMacroElectrons)
                << ", "
                << static_cast<Double>(globalMaximumMacroElectrons)
                << "], input relative-weight range=["
                << globalMinimumInputWeight << ", "
                << globalMaximumInputWeight << "].";
        logRoot(communicator_, message.str());
      }

    SILabPlaneProjectionReport planeProjection;
    if (inputReport.laboratoryPlaneCoordinates)
      {
        if (!config_.reference.inputPlaneZSet)
          throw std::invalid_argument(
            "Particle HDF5 format version 3 or 4 contains fixed-lab-plane records; beam.reference.input_plane_z is required in the shared laboratory coordinate system");
        if (inputReport.laboratoryPlaneTimeCoordinates)
          planeProjection =
            SIBunchPreprocessor::projectLabPlaneEventsToSnapshot(
              particles_, config_.reference.inputPlaneZ,
              config_.reference.initialCenterZ,
              config_.beam.positionOffset[2], communicator_);
        else
          planeProjection = SIBunchPreprocessor::projectLabPlaneToSnapshot(
            particles_, config_.reference.inputPlaneZ,
            config_.reference.initialCenterZ, communicator_);
      }
    else if (config_.reference.inputPlaneZSet && rank_ == 0)
      logRoot(communicator_,
        "WARNING: beam.reference.input_plane_z is ignored because this particle input is a legacy common-time snapshot rather than an HDF5-v3/v4 lab-plane record set.");

    initializeFieldDetectorRegions();
    validateBeamlineExclusionRules();
    const bool hasMagneticElements = !config_.magnets.empty();
    const Double firstPhysicalEntrance = hasMagneticElements ?
      firstMagneticPhysicalEntranceLab() :
      std::numeric_limits<Double>::quiet_NaN();
    const Double firstInteractionEntrance = hasMagneticElements ?
      firstMagneticInteractionEntranceLab() :
      std::numeric_limits<Double>::quiet_NaN();

    SIBunchPlacement placement;
    placement.hasFirstInteractionEntrance = hasMagneticElements;
    placement.firstInteractionEntranceLab = firstInteractionEntrance;
    placement.referencePositionLab = config_.reference.initialCenterZ;
    /* A box-frame z cell spans gamma*dz at fixed box time in the lab.  Keep
     * that as the minimum recommended clearance beyond the relative head. */
    placement.recommendationMarginLab =
      config_.mesh.boostGamma * globalGeometry_.dz;

    const SIBunchPlacementReport placementReport =
      SIBunchPreprocessor::placeRelativeLabSnapshot(
        particles_, placement, communicator_);

    Double physicalHeadLimit = std::numeric_limits<Double>::infinity();
    Double placementTolerance = 0.0;
    if (hasMagneticElements)
      {
        physicalHeadLimit = firstInteractionEntrance -
          placement.recommendationMarginLab;
        const Double placementScale = std::max(1.0,
          std::max(std::abs(firstInteractionEntrance),
                   std::abs(placementReport.headAfterLab)));
        placementTolerance = 128.0 *
          std::numeric_limits<Double>::epsilon() * placementScale;
        if (placementReport.headAfterLab >=
            physicalHeadLimit - placementTolerance)
          {
            const Double maximumCenter = physicalHeadLimit -
              placementReport.relativeHeadLab;
            std::ostringstream message;
            message << std::setprecision(10)
              << "The reconstructed common-time laboratory bunch does not fit "
                 "before the first magnetic interaction region: bunch_head_z="
              << placementReport.headAfterLab
              << " m, required_head_limit=" << physicalHeadLimit
              << " m, magnetic_interaction_start_z="
              << firstInteractionEntrance << " m. Set "
                 "beam.reference.initial_center_z to at most "
              << maximumCenter << " m (" << maximumCenter /
                   config_.inputUnits.length
              << " in the configured length unit)";
            if (inputReport.laboratoryPlaneCoordinates)
              {
                const Double minimumCenter =
                  planeProjection.recommendedMinimumReferencePosition;
                message << "; forward projection from input_plane_z requires it "
                           "to be at least " << minimumCenter << " m ("
                        << minimumCenter / config_.inputUnits.length << ")";
                if (minimumCenter >= maximumCenter - placementTolerance)
                  message << ". No admissible centre interval remains: move the "
                             "Elegant plane upstream, move/redefine the magnetic "
                             "interaction region downstream, reduce the bunch "
                             "longitudinal span, or refine dz to reduce the "
                             "one-cell safety margin";
              }
            message << ".";
            throw std::runtime_error(message.str());
          }
      }

    /* Anchor t_box=0 to the downstream bunch-front event at the reconstructed
     * lab snapshot.  The remaining particles are synchronized relative to
     * this event.  Their potentially large lab-time span is not interpreted
     * as a required physical drift upstream of the input plane. */
    frame_.setOriginsFromGamma(config_.mesh.boostGamma, SI::c, 0.0,
      placementReport.headAfterLab, 0.0);
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

    MPI_Allreduce(&localBoostReport.minimumLabGamma,
                  &boostReport.minimumLabGamma,
                  1, MPI_DOUBLE, MPI_MIN, communicator_);
    MPI_Allreduce(&localBoostReport.maximumLabGamma,
                  &boostReport.maximumLabGamma,
                  1, MPI_DOUBLE, MPI_MAX, communicator_);
    MPI_Allreduce(
      &localBoostReport.maximumRelativeMomentumRoundTripError,
      &boostReport.maximumRelativeMomentumRoundTripError,
      1, MPI_DOUBLE, MPI_MAX, communicator_);
    MPI_Allreduce(
      &localBoostReport.maximumRelativeGammaRoundTripError,
      &boostReport.maximumRelativeGammaRoundTripError,
      1, MPI_DOUBLE, MPI_MAX, communicator_);

    const Double roundTripLimit = 1.0e-10;
    if (boostReport.maximumRelativeMomentumRoundTripError >
          roundTripLimit ||
        boostReport.maximumRelativeGammaRoundTripError > roundTripLimit)
      {
        std::ostringstream message;
        message << std::setprecision(10)
                << "Lorentz boost round-trip error exceeds "
                << roundTripLimit << ": relative momentum error="
                << boostReport.maximumRelativeMomentumRoundTripError
                << ", relative gamma error="
                << boostReport.maximumRelativeGammaRoundTripError
                << ". Reduce mesh.boost_gamma, remove nonphysical extreme "
                   "input momenta, or choose a frame whose particle "
                   "momenta remain representable in double precision.";
        throw std::runtime_error(message.str());
      }

    Double localBoostedMinimum = std::numeric_limits<Double>::infinity();
    Double localBoostedMaximum = -std::numeric_limits<Double>::infinity();
    for (std::size_t index = 0; index < particles_.size(); ++index)
      {
        localBoostedMinimum = std::min(
          localBoostedMinimum, particles_[index].position[2]);
        localBoostedMaximum = std::max(
          localBoostedMaximum, particles_[index].position[2]);
      }
    Double boostedMinimum = 0.0;
    Double boostedMaximum = 0.0;
    MPI_Allreduce(&localBoostedMinimum, &boostedMinimum, 1,
                  MPI_DOUBLE, MPI_MIN, communicator_);
    MPI_Allreduce(&localBoostedMaximum, &boostedMaximum, 1,
                  MPI_DOUBLE, MPI_MAX, communicator_);
    const Double boxTranslation = config_.mesh.center[2] -
      0.5 * (boostedMinimum + boostedMaximum);
    for (std::size_t index = 0; index < particles_.size(); ++index)
      particles_[index].position[2] += boxTranslation;
    boostedMinimum += boxTranslation;
    boostedMaximum += boxTranslation;
    frame_.setOriginsFromGamma(config_.mesh.boostGamma, SI::c, 0.0,
      placementReport.headAfterLab, boxTranslation);
    referenceCenterBoxZ_ = frame_.boxZFromLabZT(
      config_.reference.initialCenterZ, 0.0);

    Double localMaximumAdvance = 0.0;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      {
        const Double gammaBox =
          BoostFrameTransform::gammaFromProperVelocity(
            particles_[index].properVelocity);
        const Double betaBoxZ =
          particles_[index].properVelocity[2] / gammaBox;
        localMaximumAdvance = std::max(localMaximumAdvance,
          frame_.gamma() * SI::c * globalGeometry_.dt *
          (frame_.beta() + betaBoxZ));
      }
    Double maximumAdvance = 0.0;
    MPI_Allreduce(&localMaximumAdvance, &maximumAdvance, 1, MPI_DOUBLE,
                  MPI_MAX, communicator_);
    particleSubsteps_ = 1;
    Double minimumUndulatorPeriod =
      std::numeric_limits<Double>::infinity();
    for (std::size_t magnet = 0; magnet < config_.magnets.size(); ++magnet)
      if (config_.magnets[magnet].type ==
          SIMagnetType::PlanarUndulator)
        minimumUndulatorPeriod = std::min(
          minimumUndulatorPeriod, config_.magnets[magnet].period);
    if (std::isfinite(minimumUndulatorPeriod))
      {
        if (!(maximumAdvance > 0.0) || !std::isfinite(maximumAdvance))
          throw std::runtime_error(
            "Cannot derive a positive lab advance per Maxwell step for undulator sampling");
        const Double fieldSteps =
          minimumUndulatorPeriod / maximumAdvance;
        const Double requestedSteps = static_cast<Double>(
          config_.mesh.particleStepsPerUndulatorPeriod);
        const Double substepRatio = requestedSteps / fieldSteps;
        const Double roundingTolerance = 128.0 *
          std::numeric_limits<Double>::epsilon() *
          std::max(1.0, std::abs(substepRatio));
        const Double requiredSubstepsReal = std::max(
          1.0, std::ceil(substepRatio - roundingTolerance));
        if (!std::isfinite(requiredSubstepsReal) ||
            requiredSubstepsReal > static_cast<Double>(
              std::numeric_limits<unsigned int>::max()))
          throw std::runtime_error(
            "Required particle substep count is outside the supported unsigned-int range; refine the Maxwell time step");
        const unsigned int requiredSubsteps =
          static_cast<unsigned int>(requiredSubstepsReal);
        if (requiredSubsteps > config_.mesh.maximumParticleSubsteps)
          {
            const Double maximumDt = globalGeometry_.dt * fieldSteps *
              static_cast<Double>(config_.mesh.maximumParticleSubsteps) /
              requestedSteps;
            const Double maximumDz = SI::c * maximumDt;
            std::ostringstream message;
            message << std::setprecision(10)
                    << "Undulator integration requires "
                    << requiredSubsteps
                    << " particle substeps per Maxwell step, above "
                       "mesh.maximum_particle_substeps="
                    << config_.mesh.maximumParticleSubsteps
                    << ". The shortest period has " << fieldSteps
                    << " Maxwell steps and requests "
                    << config_.mesh.particleStepsPerUndulatorPeriod
                    << " particle samples. Increase "
                       "mesh.maximum_particle_substeps to at least "
                    << requiredSubsteps
                    << ", or reduce the Maxwell time step to at most "
                    << maximumDt << " s";
            if (config_.mesh.fieldSolver == EBMaxwellSolver::CowanZ)
              message << " by setting mesh.cell_size z to at most "
                      << maximumDz << " m (" << maximumDz /
                           config_.inputUnits.length
                      << " in the configured length unit); keep dx and dy "
                         "not smaller than the new dz.";
            else
              message << " by refining the Yee mesh; its CFL time step is "
                         "derived from all three spacings.";
            throw std::runtime_error(message.str());
          }
        particleSubsteps_ = requiredSubsteps;
        if (rank_ == 0)
          {
            std::ostringstream message;
            message << std::setprecision(10)
                    << "Particle subcycling: Maxwell steps per shortest "
                       "undulator period=" << fieldSteps
                    << ", Boris substeps per Maxwell step="
                    << particleSubsteps_
                    << ", realized prescribed-device samples per period="
                    << fieldSteps * static_cast<Double>(particleSubsteps_)
                    << ", requested="
                    << config_.mesh.particleStepsPerUndulatorPeriod
                    << ", shortest period=" << minimumUndulatorPeriod
                    << " m. Grid E/B, charge-conserving current deposition, "
                       "and detector cadence remain on the Maxwell step; "
                       "subcycling does not relax radiation-band field sampling.";
            logRoot(communicator_, message.str());
          }
      }
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
    if (hasMagneticElements)
      {
        const Double entranceTolerance = 64.0 *
          std::numeric_limits<Double>::epsilon() *
          std::max(1.0, std::max(std::abs(firstInteractionEntrance),
                                 std::abs(eventHead)));
        if (eventHead >= physicalHeadLimit - entranceTolerance)
          {
            const Double transformedHeadOffset = eventHead -
              config_.reference.initialCenterZ;
            const Double recommendedCenter = physicalHeadLimit -
              transformedHeadOffset;
            std::ostringstream message;
            message << std::setprecision(10)
              << "Head-anchored Lorentz synchronization places a particle too "
                 "close to the first magnetic interaction region: "
                 "interaction_start_z="
              << firstInteractionEntrance
              << " m, transformed_front_z=" << eventHead
              << " m. Set beam.reference.initial_center_z to at most "
              << recommendedCenter << " m (" << recommendedCenter /
                   config_.inputUnits.length
              << " in the configured length unit). The physical first-magnet "
              "entrance remains z=" << firstPhysicalEntrance << " m.";
            throw std::runtime_error(message.str());
          }
      }
    redistributeParticles();
    validateParticlesInsideGlobalBox();

    if (rank_ == 0)
      {
        if (inputReport.laboratoryPlaneCoordinates)
          {
            std::ostringstream projectionMessage;
            projectionMessage << std::setprecision(10)
              << "Elegant lab-plane reconstruction (HDF5 v"
              << inputReport.fileFormatVersion << "): input_plane_z="
              << planeProjection.inputPlaneLab
              << " m, snapshot_reference_z="
              << planeProjection.referencePositionLab
              << " m, per-particle forward distance range=["
              << planeProjection.minimumForwardDistance << ", "
              << planeProjection.maximumForwardDistance
              << "] m. ";
            if (inputReport.laboratoryPlaneTimeCoordinates)
              projectionMessage
                << "The original arrival-time events were synchronized with "
                   "their individual velocities; mean beta_z="
                << planeProjection.meanLongitudinalBeta
                << ", common time offset="
                << planeProjection.referenceTimeOffsetLab << " s.";
            else
              projectionMessage
                << "Transverse coordinates were advanced with ux/uz and "
                   "uy/uz; longitudinal offsets retain their input timing "
                   "definition.";
            logRoot(communicator_, projectionMessage.str());
          }
        std::ostringstream placementMessage;
        placementMessage << std::setprecision(10)
          << "E/B bunch placement: particles="
          << placementReport.particles << ", reference z [m]="
          << placementReport.referencePositionLab
          << ", relative z range [m]=["
          << placementReport.relativeTailLab << ", "
          << placementReport.relativeHeadLab << "]";
        if (hasMagneticElements)
          placementMessage
            << ", physical first magnetic entrance [m]="
            << firstPhysicalEntrance
            << ", magnetic interaction start [m]="
            << firstInteractionEntrance
            << ", transformed interaction gap [m]="
            << firstInteractionEntrance - eventHead
            << ", required one-cell lab margin [m]="
            << placement.recommendationMarginLab << ".";
        else
          placementMessage
            << ". No magnetic element is configured; detector-only or "
               "element-free propagation has no magnetic entrance "
               "placement constraint.";
        logRoot(communicator_, placementMessage.str());
        std::ostringstream boostMessage;
        boostMessage << std::setprecision(10)
          << "Head-anchored Lorentz synchronization: anchor_lab_z="
          << placementReport.headAfterLab << " m, virtual lab event times [s]=["
          << boostReport.earliestLabEventTime << " to "
          << boostReport.latestLabEventTime << "], maximum synchronization "
             "span=" << boostReport.maximumAbsoluteDriftTime << " s ("
          << SI::c * boostReport.maximumAbsoluteDriftTime
          << " light-metres; diagnostic only, not a required physical drift), "
             "centred boosted z range [m]=[" << boostedMinimum << ", "
          << boostedMaximum << "]"
          << "; lab gamma range=[" << boostReport.minimumLabGamma
          << ", " << boostReport.maximumLabGamma
          << "], maximum relative boost round-trip errors: momentum="
          << boostReport.maximumRelativeMomentumRoundTripError
          << ", gamma=" << boostReport.maximumRelativeGammaRoundTripError
          << ". Individual-electron lab-energy roundoff scale at gamma_max="
          << std::numeric_limits<Double>::epsilon() *
             boostReport.maximumLabGamma * SI::electronMass *
             SI::c * SI::c / SI::elementaryCharge
          << " eV; energy-loss analysis must subtract per-particle gamma "
             "with compensated or long-double accumulation, not two rounded "
             "total beam energies.";
        logRoot(communicator_, boostMessage.str());
        if (boostReport.maximumRelativeMomentumRoundTripError > 1.0e-13 ||
            boostReport.maximumRelativeGammaRoundTripError > 1.0e-13)
          logRoot(communicator_,
            "WARNING: Lorentz round-trip error is above 1e-13. Consider a less aggressive boost_gamma and verify exported gamma differences before interpreting small particle-energy losses.");
      }
  }

  void Simulation::reportDurationEstimate() const
  {
    Double estimatedBoxTime = 0.0;
    bool reachable = true;
    const char* estimateModel = 0;

    if (config_.stop.mode == StopMode::ReferenceCenterZ)
      {
        estimateModel = "boost-reference-worldline";
        const Double currentZ = frame_.labZFromBoxZT(
          referenceCenterBoxZ_, 0.0);
        const Double distance = config_.stop.referenceZ - currentZ;
        const Double labAdvancePerBoxSecond =
          frame_.gammaBeta() * SI::c;
        if (distance > 0.0)
          {
            reachable = labAdvancePerBoxSecond > 0.0 &&
                        std::isfinite(labAdvancePerBoxSecond);
            if (reachable)
              estimatedBoxTime = distance / labAdvancePerBoxSecond;
          }
      }
    else
      {
        estimateModel = "initial-ballistic-particles";
        const Double interactionExit = lastBeamlineInteractionExitLab();
        Double localMaximumTime = 0.0;
        int localUnreachable = 0;
        for (std::size_t index = 0; index < particles_.size(); ++index)
          {
            Double eventTimeLab = 0.0;
            Double eventZLab = 0.0;
            frame_.boxToLab(0.0, particles_[index].position[2],
                            eventTimeLab, eventZLab);
            const Double remaining = interactionExit - eventZLab;
            if (remaining <= 0.0) continue;
            const Double gammaBox =
              BoostFrameTransform::gammaFromProperVelocity(
                particles_[index].properVelocity);
            const Double betaBoxZ =
              particles_[index].properVelocity[2] / gammaBox;
            const Double labAdvancePerBoxSecond =
              frame_.gamma() * SI::c * (frame_.beta() + betaBoxZ);
            if (!(labAdvancePerBoxSecond > 0.0) ||
                !std::isfinite(labAdvancePerBoxSecond))
              localUnreachable = 1;
            else
              localMaximumTime = std::max(localMaximumTime,
                remaining / labAdvancePerBoxSecond);
          }
        int globalUnreachable = 0;
        MPI_Allreduce(&localUnreachable, &globalUnreachable, 1, MPI_INT,
                      MPI_MAX, communicator_);
        MPI_Allreduce(&localMaximumTime, &estimatedBoxTime, 1, MPI_DOUBLE,
                      MPI_MAX, communicator_);
        reachable = globalUnreachable == 0;
      }

    if (!reachable || !std::isfinite(estimatedBoxTime))
      {
        if (rank_ == 0)
          {
            std::ostringstream warning;
            warning << "[duration-estimate] stop_mode="
                    << (config_.stop.mode == StopMode::ReferenceCenterZ ?
                        "reference-center-z" : "after-last-element")
                    << " model=" << estimateModel
                    << " status=unreachable-from-initial-kinematics "
                       "configured_box_s=" << totalTimeBoxSI_
                    << ". The configured mesh.duration remains the hard "
                       "runtime guard; choose a moving boost reference or "
                       "verify the initial longitudinal particle velocities.";
            logRoot(communicator_, warning.str());
          }
        return;
      }

    const Double estimatedStepsReal = std::max(1.0,
      std::ceil(estimatedBoxTime / globalGeometry_.dt));
    const Double suggestedDuration =
      estimatedStepsReal * globalGeometry_.dt;
    if (rank_ != 0) return;

    std::ostringstream message;
    message << std::setprecision(10)
            << "[duration-estimate] stop_mode="
            << (config_.stop.mode == StopMode::ReferenceCenterZ ?
                "reference-center-z" : "after-last-element")
            << " model=" << estimateModel
            << " estimated_box_s=" << estimatedBoxTime
            << " estimated_steps=" << estimatedStepsReal
            << " minimum_step_rounded_box_s=" << suggestedDuration
            << " configured_box_s=" << totalTimeBoxSI_
            << " headroom="
            << (suggestedDuration > 0.0 ?
                totalTimeBoxSI_ / suggestedDuration : 0.0) << ".";
    logRoot(communicator_, message.str());

    const Double comparisonTolerance = 64.0 *
      std::numeric_limits<Double>::epsilon() *
      std::max(totalTimeBoxSI_, suggestedDuration);
    if (totalTimeBoxSI_ + comparisonTolerance < suggestedDuration)
      {
        std::ostringstream warning;
        warning << std::setprecision(10)
                << "WARNING: mesh.duration is shorter than the startup "
                   "stop-time estimate. Set mesh.duration to at least "
                << suggestedDuration / config_.inputUnits.time
                << " in the configured time unit ("
                << suggestedDuration << " s), then add operational margin. ";
        if (config_.stop.mode == StopMode::AfterLastElement)
          warning << "This estimate propagates the initial particle states "
                     "ballistically; magnetic dynamics can change the actual "
                     "stop time.";
        else
          warning << "This estimate follows the same inertial boost-reference "
                     "worldline used by the stop predicate.";
        logRoot(communicator_, warning.str());
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
    const Double rho = std::hypot(config_.mesh.extent[0],
                                  config_.mesh.extent[1]);
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
      config_.runtime.interactive(), config_.output.overwrite,
      &runMetadata_);
    nextTrajectorySampleTime_ = 0.0;
    trajectorySamplesSinceFlush_ = 0;
  }

  void Simulation::initializeDetectorOutput()
  {
    if (!detectorConfig_.enabled()) return;
    detectors_.reset(new LabDetectorManager(
      detectorConfig_, globalGeometry_,
      globalOriginBox_, localOriginBox_, frame_, communicator_, NULL,
      config_.output.overwrite, &runMetadata_));
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

  long double Simulation::localInteriorFieldEnergy() const
  {
    std::size_t layer[3] = {0, 0, 0};
    if (config_.boundary.type == EBBoundaryType::Cpml)
      for (unsigned int axis = 0; axis < 3; ++axis)
        layer[axis] = config_.boundary.cpml.cells[axis];

    const std::size_t iBegin = layer[0];
    const std::size_t iEnd = globalGeometry_.nx - layer[0];
    const std::size_t jBegin = layer[1];
    const std::size_t jEnd = globalGeometry_.ny - layer[1];
    const std::size_t globalKBegin = layer[2];
    const std::size_t globalKEnd = globalGeometry_.nz - layer[2];
    const std::size_t localGlobalEnd = localZOffset_ + localGeometry_.nz;
    const std::size_t ownedKBegin = std::max(globalKBegin, localZOffset_);
    const std::size_t ownedKEnd = std::min(globalKEnd, localGlobalEnd);
    if (iBegin >= iEnd || jBegin >= jEnd || ownedKBegin >= ownedKEnd)
      return 0.0L;

    const long double volume =
      static_cast<long double>(globalGeometry_.dx) *
      static_cast<long double>(globalGeometry_.dy) *
      static_cast<long double>(globalGeometry_.dz);
    long double energy = 0.0L;
    for (std::size_t globalK = ownedKBegin;
         globalK < ownedKEnd; ++globalK)
      {
        const std::size_t k = globalK - localZOffset_;
        for (std::size_t j = jBegin; j < jEnd; ++j)
          for (std::size_t i = iBegin; i < iEnd; ++i)
            energy += static_cast<long double>(
              fields_->radiationSampleCell(i, j, k).energyDensity) * volume;
      }
    return energy;
  }

  void Simulation::localInteriorBoundaryPower(long double power[6]) const
  {
    for (unsigned int face = 0; face < 6; ++face) power[face] = 0.0L;
    if (config_.boundary.type != EBBoundaryType::Cpml) return;

    const std::size_t cx = config_.boundary.cpml.cells[0];
    const std::size_t cy = config_.boundary.cpml.cells[1];
    const std::size_t cz = config_.boundary.cpml.cells[2];
    const std::size_t iBegin = cx;
    const std::size_t iEnd = globalGeometry_.nx - cx;
    const std::size_t jBegin = cy;
    const std::size_t jEnd = globalGeometry_.ny - cy;
    const std::size_t globalKBegin = cz;
    const std::size_t globalKEnd = globalGeometry_.nz - cz;
    const std::size_t localGlobalEnd = localZOffset_ + localGeometry_.nz;
    const std::size_t ownedKBegin = std::max(globalKBegin, localZOffset_);
    const std::size_t ownedKEnd = std::min(globalKEnd, localGlobalEnd);
    if (iBegin >= iEnd || jBegin >= jEnd || globalKBegin >= globalKEnd)
      return;

    if (ownedKBegin < ownedKEnd && cx > 0)
      {
        const long double area =
          static_cast<long double>(globalGeometry_.dy) *
          static_cast<long double>(globalGeometry_.dz);
        const std::size_t lowerI = iBegin;
        const std::size_t upperI = iEnd - 1;
        for (std::size_t globalK = ownedKBegin;
             globalK < ownedKEnd; ++globalK)
          {
            const std::size_t k = globalK - localZOffset_;
            for (std::size_t j = jBegin; j < jEnd; ++j)
              {
                power[0] -= static_cast<long double>(
                  fields_->radiationSampleCell(lowerI, j, k).poynting[0]) *
                  area;
                power[1] += static_cast<long double>(
                  fields_->radiationSampleCell(upperI, j, k).poynting[0]) *
                  area;
              }
          }
      }

    if (ownedKBegin < ownedKEnd && cy > 0)
      {
        const long double area =
          static_cast<long double>(globalGeometry_.dx) *
          static_cast<long double>(globalGeometry_.dz);
        const std::size_t lowerJ = jBegin;
        const std::size_t upperJ = jEnd - 1;
        for (std::size_t globalK = ownedKBegin;
             globalK < ownedKEnd; ++globalK)
          {
            const std::size_t k = globalK - localZOffset_;
            for (std::size_t i = iBegin; i < iEnd; ++i)
              {
                power[2] -= static_cast<long double>(
                  fields_->radiationSampleCell(i, lowerJ, k).poynting[1]) *
                  area;
                power[3] += static_cast<long double>(
                  fields_->radiationSampleCell(i, upperJ, k).poynting[1]) *
                  area;
              }
          }
      }

    if (cz > 0)
      {
        const long double area =
          static_cast<long double>(globalGeometry_.dx) *
          static_cast<long double>(globalGeometry_.dy);
        const std::size_t lowerGlobalK = globalKBegin;
        const std::size_t upperGlobalK = globalKEnd - 1;
        if (lowerGlobalK >= localZOffset_ &&
            lowerGlobalK < localGlobalEnd)
          {
            const std::size_t k = lowerGlobalK - localZOffset_;
            for (std::size_t j = jBegin; j < jEnd; ++j)
              for (std::size_t i = iBegin; i < iEnd; ++i)
                power[4] -= static_cast<long double>(
                  fields_->radiationSampleCell(i, j, k).poynting[2]) * area;
          }
        if (upperGlobalK >= localZOffset_ &&
            upperGlobalK < localGlobalEnd)
          {
            const std::size_t k = upperGlobalK - localZOffset_;
            for (std::size_t j = jBegin; j < jEnd; ++j)
              for (std::size_t i = iBegin; i < iEnd; ++i)
                power[5] += static_cast<long double>(
                  fields_->radiationSampleCell(i, j, k).poynting[2]) * area;
          }
      }
  }

  void Simulation::initializeEnergyLedger()
  {
    if (!config_.energyLedger.enabled) return;
    createDirectories(config_.energyLedger.directory, communicator_);
    if (rank_ == 0)
      energyLedgerWriter_.open(joinPath(config_.energyLedger.directory,
          config_.energyLedger.filename), size_,
        config_.energyLedger.bufferRecords,
        config_.energyLedger.compression,
        config_.runtime.interactive(), config_.output.overwrite,
        &runMetadata_);
    localInteriorBoundaryPower(energyLedgerPreviousPower_);
    if (rank_ == 0)
      {
        std::ostringstream message;
        message << "Energy ledger active: root-only HDF5='"
                << joinPath(config_.energyLedger.directory,
                            config_.energyLedger.filename)
                << "', volume/spread sample interval="
                << config_.energyLedger.sampleIntervalSteps
                << " field steps. Boundary power is accumulated every step; "
                   "analytic-device work is accumulated along Boris substeps.";
        logRoot(communicator_, message.str());
        if (sources_.maxwellIncidentWaveCount() > 0)
          logRoot(communicator_,
            "WARNING: energy-ledger closure_valid is false because incident-wave boundary work is not yet included.");
      }
  }

  void Simulation::advanceEnergyLedgerStep()
  {
    long double currentPower[6] = {};
    localInteriorBoundaryPower(currentPower);
    const long double halfDt = 0.5L *
      static_cast<long double>(globalGeometry_.dt);
    for (unsigned int face = 0; face < 6; ++face)
      {
        energyLedgerOutwardEnergy_[face] += halfDt *
          (energyLedgerPreviousPower_[face] + currentPower[face]);
        energyLedgerPreviousPower_[face] = currentPower[face];
      }
  }

  void Simulation::sampleEnergyLedger(bool force)
  {
    if (!config_.energyLedger.enabled) return;
    if (energyLedgerHasReference_ && step_ == energyLedgerLastSampleStep_)
      return;
    if (!force &&
        step_ % config_.energyLedger.sampleIntervalSteps != 0)
      return;

    long double local[16] = {};
    Double localMinimumGamma = std::numeric_limits<Double>::infinity();
    Double localMaximumGamma =
      -std::numeric_limits<Double>::infinity();
    for (std::size_t index = 0; index < particles_.size(); ++index)
      {
        const RelativisticParticleSI& particle = particles_[index];
        const Double gammaBox =
          BoostFrameTransform::gammaFromProperVelocity(
            particle.properVelocity);
        const long double mass =
          static_cast<long double>(particle.mass);
        const long double c = static_cast<long double>(SI::c);
        const long double rest = mass * c * c;
        const long double u2 = static_cast<long double>(
          particle.properVelocity.norm2());
        local[0] += rest * u2 /
          static_cast<long double>(gammaBox + 1.0);
        local[1] += rest * static_cast<long double>(gammaBox);
        const long double represented = mass /
          static_cast<long double>(SI::electronMass);
        local[2] += represented;
        FieldVector<Double> properVelocityLab(0.0);
        frame_.properVelocityBoxToLab(
          particle.properVelocity, properVelocityLab);
        const Double gammaLab =
          BoostFrameTransform::gammaFromProperVelocity(properVelocityLab);
        local[3] += represented * static_cast<long double>(gammaLab);
        local[4] += represented * static_cast<long double>(gammaBox);
        local[5] += represented * static_cast<long double>(
          frame_.labZFromBoxZT(particle.position[2], timeBoxSI_));
        localMinimumGamma = std::min(localMinimumGamma, gammaLab);
        localMaximumGamma = std::max(localMaximumGamma, gammaLab);
      }
    local[6] = localInteriorFieldEnergy();
    local[7] = energyLedgerRemovedKinetic_;
    local[8] = energyLedgerRemovedTotal_;
    local[9] = energyLedgerPrescribedWork_;
    for (unsigned int face = 0; face < 6; ++face)
      local[10 + face] = energyLedgerOutwardEnergy_[face];

    long double global[16] = {};
    MPI_Allreduce(local, global, 16, MPI_LONG_DOUBLE, MPI_SUM,
                  communicator_);
    unsigned long long localCounts[3] = {
      static_cast<unsigned long long>(particles_.size()),
      energyLedgerRemovedMacroparticles_, retirementEntryCount_
    };
    unsigned long long globalCounts[3] = {};
    MPI_Allreduce(localCounts, globalCounts, 3,
      MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator_);
    Double globalMinimumGamma = 0.0;
    Double globalMaximumGamma = 0.0;
    MPI_Allreduce(&localMinimumGamma, &globalMinimumGamma, 1,
                  MPI_DOUBLE, MPI_MIN, communicator_);
    MPI_Allreduce(&localMaximumGamma, &globalMaximumGamma, 1,
                  MPI_DOUBLE, MPI_MAX, communicator_);

    const long double weight = global[2];
    const long double meanGammaLab = weight > 0.0L ?
      global[3] / weight : 0.0L;
    const long double meanGammaBox = weight > 0.0L ?
      global[4] / weight : 0.0L;
    const long double meanZLab = weight > 0.0L ?
      global[5] / weight : 0.0L;
    long double localCentred[4] = {};
    if (weight > 0.0L)
      for (std::size_t index = 0; index < particles_.size(); ++index)
        {
          const RelativisticParticleSI& particle = particles_[index];
          const long double represented =
            static_cast<long double>(particle.mass) /
            static_cast<long double>(SI::electronMass);
          const Double gammaBox =
            BoostFrameTransform::gammaFromProperVelocity(
              particle.properVelocity);
          FieldVector<Double> properVelocityLab(0.0);
          frame_.properVelocityBoxToLab(
            particle.properVelocity, properVelocityLab);
          const Double gammaLab =
            BoostFrameTransform::gammaFromProperVelocity(
              properVelocityLab);
          const long double deltaGamma =
            static_cast<long double>(gammaLab) - meanGammaLab;
          const long double deltaGammaBox =
            static_cast<long double>(gammaBox) - meanGammaBox;
          const long double deltaZ = static_cast<long double>(
            frame_.labZFromBoxZT(particle.position[2], timeBoxSI_)) -
            meanZLab;
          localCentred[0] += represented * deltaGamma * deltaGamma;
          localCentred[1] += represented * deltaGammaBox * deltaGammaBox;
          localCentred[2] += represented * deltaZ * deltaZ;
          localCentred[3] += represented * deltaZ * deltaGamma;
        }
    long double globalCentred[4] = {};
    MPI_Allreduce(localCentred, globalCentred, 4, MPI_LONG_DOUBLE,
                  MPI_SUM, communicator_);

    if (!energyLedgerHasReference_)
      {
        energyLedgerInitialKinetic_ = global[0];
        energyLedgerInitialField_ = global[6];
        energyLedgerHasReference_ = true;
      }
    const long double particleAccounted = global[0] + global[7];
    const long double initialDynamic =
      energyLedgerInitialKinetic_ + energyLedgerInitialField_;
    const long double currentDynamic = particleAccounted + global[6];
    long double outward = 0.0L;
    for (unsigned int face = 0; face < 6; ++face)
      outward += global[10 + face];
    const long double residual = currentDynamic - initialDynamic +
      outward - global[9];
    const long double deltaParticle =
      particleAccounted - energyLedgerInitialKinetic_;
    const long double deltaField = global[6] - energyLedgerInitialField_;
    const long double exchangeScale = std::abs(deltaParticle) +
      std::abs(deltaField) + std::abs(outward) + std::abs(global[9]);
    long double varianceGamma = weight > 0.0L ?
      std::max(0.0L, globalCentred[0] / weight) : 0.0L;
    long double varianceGammaBox = weight > 0.0L ?
      std::max(0.0L, globalCentred[1] / weight) : 0.0L;
    const long double varianceZ = weight > 0.0L ?
      std::max(0.0L, globalCentred[2] / weight) : 0.0L;
    const long double covariance = weight > 0.0L ?
      globalCentred[3] / weight : 0.0L;
    const long double gammaResolution = 64.0L *
      static_cast<long double>(std::numeric_limits<Double>::epsilon()) *
      std::max(1.0L, std::abs(meanGammaLab));
    if (varianceGamma < gammaResolution * gammaResolution)
      varianceGamma = 0.0L;
    const long double gammaBoxResolution = 64.0L *
      static_cast<long double>(std::numeric_limits<Double>::epsilon()) *
      std::max(1.0L, std::abs(meanGammaBox));
    if (varianceGammaBox < gammaBoxResolution * gammaBoxResolution)
      varianceGammaBox = 0.0L;
    const long double chirp = varianceZ > 0.0L && varianceGamma > 0.0L ?
      covariance / varianceZ : 0.0L;
    const long double uncorrelatedVariance =
      varianceZ > 0.0L && varianceGamma > 0.0L ?
      std::max(0.0L, varianceGamma - covariance * covariance /
        varianceZ) : varianceGamma;

    EnergyLedgerRecord record;
    record.step = static_cast<std::uint64_t>(step_);
    record.activeMacroparticles = globalCounts[0];
    record.removedMacroparticles = globalCounts[1];
    record.closureValid =
      sources_.maxwellIncidentWaveCount() == 0 && globalCounts[2] == 0 ? 1 : 0;
    record.timeBox = timeBoxSI_;
    record.activeRepresentedElectrons = static_cast<Double>(global[2]);
    record.particleKineticActive = static_cast<Double>(global[0]);
    record.particleTotalActive = static_cast<Double>(global[1]);
    record.particleKineticRemoved = static_cast<Double>(global[7]);
    record.particleTotalRemoved = static_cast<Double>(global[8]);
    record.fieldEnergyInterior = static_cast<Double>(global[6]);
    for (unsigned int face = 0; face < 6; ++face)
      record.outwardFieldEnergyFace[face] =
        static_cast<Double>(global[10 + face]);
    record.outwardFieldEnergy = static_cast<Double>(outward);
    record.prescribedWork = static_cast<Double>(global[9]);
    record.balanceResidual = static_cast<Double>(residual);
    record.relativeBalanceExchange = exchangeScale > 0.0L ?
      static_cast<Double>(residual / exchangeScale) : 0.0;
    record.relativeBalanceInitial = initialDynamic != 0.0L ?
      static_cast<Double>(residual / std::abs(initialDynamic)) : 0.0;
    const long double kineticDenominator = global[0] + global[6];
    record.fieldFractionKinetic = kineticDenominator != 0.0L ?
      static_cast<Double>(global[6] / kineticDenominator) : 0.0;
    const long double totalDenominator = global[1] + global[6];
    record.fieldFractionIncludingRest = totalDenominator != 0.0L ?
      static_cast<Double>(global[6] / totalDenominator) : 0.0;
    record.meanGammaLab = static_cast<Double>(meanGammaLab);
    record.sigmaGammaLab = static_cast<Double>(std::sqrt(varianceGamma));
    record.relativeEnergySpreadLab = meanGammaLab != 0.0L ?
      static_cast<Double>(std::sqrt(varianceGamma) / meanGammaLab) : 0.0;
    record.minimumGammaLab = globalCounts[0] > 0 ? globalMinimumGamma : 0.0;
    record.maximumGammaLab = globalCounts[0] > 0 ? globalMaximumGamma : 0.0;
    record.linearChirpGammaPerM = static_cast<Double>(chirp);
    record.uncorrelatedSigmaGammaLab =
      static_cast<Double>(std::sqrt(uncorrelatedVariance));
    record.meanGammaBox = static_cast<Double>(meanGammaBox);
    record.sigmaGammaBox =
      static_cast<Double>(std::sqrt(varianceGammaBox));
    lastEnergyLedgerRecord_ = record;
    if (rank_ == 0) energyLedgerWriter_.append(record);
    energyLedgerLastSampleStep_ = step_;
  }

  void Simulation::finalizeEnergyLedger(bool completed)
  {
    if (!config_.energyLedger.enabled) return;
    if (rank_ == 0)
      {
        energyLedgerWriter_.close(completed);
        std::ostringstream message;
        message << std::setprecision(10)
          << "[energy] step=" << lastEnergyLedgerRecord_.step
          << " closure_valid="
          << static_cast<unsigned int>(lastEnergyLedgerRecord_.closureValid)
          << " particle_kinetic_active_J="
          << lastEnergyLedgerRecord_.particleKineticActive
          << " removed_particle_kinetic_J="
          << lastEnergyLedgerRecord_.particleKineticRemoved
          << " field_energy_J="
          << lastEnergyLedgerRecord_.fieldEnergyInterior
          << " outward_field_energy_J="
          << lastEnergyLedgerRecord_.outwardFieldEnergy
          << " prescribed_work_J="
          << lastEnergyLedgerRecord_.prescribedWork
          << " residual_J="
          << lastEnergyLedgerRecord_.balanceResidual
          << " relative_exchange="
          << lastEnergyLedgerRecord_.relativeBalanceExchange
          << " field_fraction_kinetic="
          << lastEnergyLedgerRecord_.fieldFractionKinetic
          << " mean_gamma_lab=" << lastEnergyLedgerRecord_.meanGammaLab
          << " relative_energy_spread_lab="
          << lastEnergyLedgerRecord_.relativeEnergySpreadLab
          << " uncorrelated_sigma_gamma_lab="
          << lastEnergyLedgerRecord_.uncorrelatedSigmaGammaLab
          << " chirp_gamma_per_m="
          << lastEnergyLedgerRecord_.linearChirpGammaPerM << ".";
        logRoot(communicator_, message.str());
        if (lastEnergyLedgerRecord_.closureValid &&
            std::abs(lastEnergyLedgerRecord_.relativeBalanceInitial) >
              config_.energyLedger.warningRelativeTolerance)
          {
            std::ostringstream warning;
            warning << std::setprecision(6)
              << "WARNING: energy-ledger residual is "
              << std::abs(lastEnergyLedgerRecord_.relativeBalanceInitial)
              << " of the initial boosted-frame kinetic-plus-field energy, "
                 "above warning_relative_tolerance="
              << config_.energyLedger.warningRelativeTolerance
              << ". Refine the field grid/time step, increase macro-particle "
                 "count, enlarge the initial-field padding, and check that "
                 "the residual converges before interpreting small radiation "
                 "or energy-spread changes.";
            logRoot(communicator_, warning.str());
          }
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
    FieldVector<Double> properVelocityLab(0.0);
    frame_.properVelocityBoxToLab(particle.properVelocity, properVelocityLab);
    for (unsigned int component = 0; component < 3; ++component)
      record.properVelocity[component] = properVelocityLab[component];

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
    const Double upperZ = std::fma(
      static_cast<Double>(localZOffset_ + localGeometry_.nz),
      globalGeometry_.dz, globalOriginBox_[2]);
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
        Double prescribedWork = 0.0;
        RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
          particle, *fields_, localOriginBox_, sources_, frame_,
          timeBoxSI_, globalGeometry_.dt, particleSubsteps_,
          config_.energyLedger.enabled ? &prescribedWork : 0);

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
        if (config_.energyLedger.enabled)
          energyLedgerPrescribedWork_ +=
            static_cast<long double>(prescribedWork) *
            static_cast<long double>(hasTerminalEvent ?
              terminalFraction : 1.0);

        RelativisticParticleSI diagnosticEnd(particle);
        Double diagnosticTime = timeBoxSI_ + globalGeometry_.dt;
        if (hasTerminalEvent)
          {
            diagnosticEnd.position = diagnosticHit.position;
            diagnosticTime = timeBoxSI_ +
              diagnosticHit.fraction * globalGeometry_.dt;
            if (config_.energyLedger.enabled)
              {
                const RelativisticParticleSI& beginning = particles_[index];
                const Double gammaBeginning =
                  BoostFrameTransform::gammaFromProperVelocity(
                    beginning.properVelocity);
                const Double gammaEnd =
                  BoostFrameTransform::gammaFromProperVelocity(
                    particle.properVelocity);
                const long double mass =
                  static_cast<long double>(particle.mass);
                const long double c = static_cast<long double>(SI::c);
                const long double rest = mass * c * c;
                const long double kineticBeginning = rest *
                  static_cast<long double>(beginning.properVelocity.norm2()) /
                  static_cast<long double>(gammaBeginning + 1.0);
                const long double kineticEnd = rest *
                  static_cast<long double>(particle.properVelocity.norm2()) /
                  static_cast<long double>(gammaEnd + 1.0);
                const long double fraction = static_cast<long double>(
                  diagnosticHit.fraction);
                energyLedgerRemovedKinetic_ += kineticBeginning + fraction *
                  (kineticEnd - kineticBeginning);
                energyLedgerRemovedTotal_ += rest *
                  (static_cast<long double>(gammaBeginning) + fraction *
                   static_cast<long double>(gammaEnd - gammaBeginning));
                ++energyLedgerRemovedMacroparticles_;
              }
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
          referenceCenterBoxZ_, timeBoxSI_);
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

  Double Simulation::firstMagneticPhysicalEntranceLab() const
  {
    Double entrance = std::numeric_limits<Double>::infinity();
    for (std::size_t element = 0; element < config_.magnets.size(); ++element)
      entrance = std::min(entrance,
        config_.magnets[element].physicalEntranceLab());
    return entrance;
  }

  Double Simulation::firstMagneticInteractionEntranceLab() const
  {
    Double entrance = std::numeric_limits<Double>::infinity();
    for (std::size_t element = 0; element < config_.magnets.size(); ++element)
      entrance = std::min(entrance,
        config_.magnets[element].interactionEntranceLab());
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
      std::fma(static_cast<Double>(globalGeometry_.nx), globalGeometry_.dx,
               globalOriginBox_[0]),
      std::fma(static_cast<Double>(globalGeometry_.ny), globalGeometry_.dy,
               globalOriginBox_[1]),
      std::fma(static_cast<Double>(globalGeometry_.nz), globalGeometry_.dz,
               globalOriginBox_[2])
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
