#include "simulation.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cerrno>
#include <exception>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>

#include "particle_initializer.h"
#include "checkpoint_io.h"
#include "compute_threads.h"
#include "runtime_control.h"
#include "runtime_util.h"
#include "shutdown_watchdog.h"
#include "teardown_trace.h"
#include "grid_work_audit.h"
#include "field_stage_audit.h"
#include "discrete_field_ledger.h"
#include "inlet_constraint_probe.h"

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

    /* Clip one straight Maxwell-step orbit to the transverse field aperture.
     * The omitted part is a field-free external particle; its first in-box
     * segment is a prescribed boundary source, not a periodic image. */
    bool clipPhysicalSegment(
        const FieldVector<Double>& start,
        const FieldVector<Double>& end,
        const ParticleCPMLRegion& region,
        Double& first, Double& last)
    {
      first = 0.0;
      last = 1.0;
      for (unsigned int axis = 0; axis < 3; ++axis)
        {
          const Double low = region.physicalLower(axis);
          const Double high = region.physicalUpper(axis);
          const Double delta = end[axis] - start[axis];
          if (delta == 0.0)
            {
              if (start[axis] < low || start[axis] > high) return false;
              continue;
            }
          Double a = (low - start[axis]) / delta;
          Double b = (high - start[axis]) / delta;
          if (a > b) std::swap(a, b);
          first = std::max(first, a);
          last = std::min(last, b);
          if (!(last > first)) return false;
        }
      return last > first;
    }

    struct ParticlePacket
    {
      Double position[3];
      Double properVelocity[3];
      Double charge;
      Double mass;
      Double weight;
      std::uint64_t id;
      std::uint64_t sourceId;
      unsigned char frozen;
      Double frozenBetaZ;
    };

    struct TransferPacket
    {
      ParticlePacket particle;
      Double segmentStart[3];
      Double segmentStartFraction;
      Double activationKink[3];
      Double activationKinkFraction;
      unsigned char hasActivationKink;
      unsigned char pendingActivePush;
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
      packet.frozen = particle.frozen ? 1 : 0;
      packet.frozenBetaZ = particle.frozenBetaZ;
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
      particle.frozen = packet.frozen != 0;
      particle.frozenBetaZ = packet.frozenBetaZ;
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
                         std::vector<Packet>& receive,
                         Double* communicationSeconds = NULL)
    {
      if (send.size() > static_cast<std::size_t>(INT_MAX))
        throw std::overflow_error("Too many particles cross one MPI interface");
      const int sendCount = static_cast<int>(send.size());
      int receiveCount = 0;
      Double start = communicationSeconds ? MPI_Wtime() : 0.0;
      MPI_Sendrecv(&sendCount, 1, MPI_INT, destination, countTag,
                   &receiveCount, 1, MPI_INT, source, countTag,
                   communicator, MPI_STATUS_IGNORE);
      if (communicationSeconds) *communicationSeconds += MPI_Wtime() - start;
      if (receiveCount < 0)
        throw std::runtime_error("Negative MPI particle transfer count");
      receive.resize(static_cast<std::size_t>(receiveCount));
      const int sendBytes = checkedBytes(send.size(), sizeof(Packet));
      const int receiveBytes = checkedBytes(receive.size(),
                                            sizeof(Packet));
      start = communicationSeconds ? MPI_Wtime() : 0.0;
      MPI_Sendrecv(send.empty() ? NULL : &send[0], sendBytes, MPI_BYTE,
                   destination, dataTag,
                   receive.empty() ? NULL : &receive[0], receiveBytes,
                   MPI_BYTE, source, dataTag, communicator,
                   MPI_STATUS_IGNORE);
      if (communicationSeconds) *communicationSeconds += MPI_Wtime() - start;
    }

  }

  /* Reuse storage across field steps. Particle preparation is bounded to a
   * small batch so enabling threads does not duplicate the whole bunch. */
  struct Simulation::ParticleWorkspace
  {
    enum { PreparationBatchCapacity = 2048 };
    struct Prepared
    {
      RelativisticParticleSI particle;
      RelativisticParticleSI activeStart;
      Double prescribedWork;
      BorisWork borisWork;
      ParticleEnergyTrace auditTrace;
      Double activeStartTime;
      bool hasActiveSegment;
      bool pendingActivePush;
      bool captureFrozenSegment;
    };
    std::vector<RelativisticParticleSI> retained;
    std::vector<ParticleCPMLCarrier> retainedCarriers, retainedRetirementCarriers;
    std::vector<TransferPacket> sendLower, sendUpper, receiveLower, receiveUpper;
    std::vector<CarrierTransferPacket> sendCarrierLower, sendCarrierUpper;
    std::vector<CarrierTransferPacket> receiveCarrierLower, receiveCarrierUpper;
    std::vector<CarrierTransferPacket> sendRetirementLower, sendRetirementUpper;
    std::vector<CarrierTransferPacket> receiveRetirementLower, receiveRetirementUpper;
    std::vector<Prepared> prepared;
  };

  Simulation::Simulation(const SimulationConfig& config,
                         MPI_Comm communicator)
    : config_(config), detectorConfig_(config.detectors),
      beamlineElements_(config.beamlineElements),
      communicator_(communicator), rank_(0), size_(1),
      runMetadata_(), longitudinalLayoutReport_(),
      hasLongitudinalLayoutReport_(false),
      globalGeometry_(), localGeometry_(), localZOffset_(0),
      globalOriginBox_(0.0), localOriginBox_(0.0),
      referenceCenterBoxZ_(0.0), frame_(), sources_(),
      fields_(), halo_(), incident_(), particles_(),
      particleWorkspace_(new ParticleWorkspace()), particleCPML_(),
      pmlCarriers_(), retirementCarriers_(), particleBoundary_(), detectors_(),
      trajectoryWriter_(), energyLedgerWriter_(), lastEnergyLedgerRecord_(),
      trajectoryRhythmSI_(0.0), nextTrajectorySampleTime_(0.0),
      trajectorySamplesSinceFlush_(0), particleSubsteps_(1),
      energyLedgerPrescribedWork_(0.0L),
      energyLedgerBoundarySourceWork_(0.0L),
      energyLedgerBoundarySourceAbsWork_(0.0L),
      energyLedgerRemovedKinetic_(0.0L),
      energyLedgerRemovedTotal_(0.0L),
      energyLedgerBorisWork_(), energyLedgerRemovedMoment_(),
      energyLedgerReleasedMoment_(), energyLedgerInitialMoment_(),
      energyLedgerReleasedCharge_(0), energyLedgerPartialWorkSubsteps_(0),
      energyLedgerInitialKinetic_(0.0L),
      energyLedgerInitialField_(0.0L),
      energyLedgerRemovedMacroparticles_(0),
      energyLedgerLastSampleStep_(0), energyLedgerHasReference_(false),
      timeBoxSI_(0.0),
      totalTimeBoxSI_(0.0), step_(0), interrupted_(false),
      configuredStopReached_(false), loopWallStart_(0.0),
      estimatedStepSeconds_(0.0), estimatedMaximumSteps_(0),
      modeledLocalPeakBytes_(0), performanceInterval_(0),
      performanceLastStep_(0), performanceWallStart_(0.0),
      localLostParticles_(0), lostParticles_(0),
      ballisticBoundaryEntries_(0), frozenPlaneActivations_(0),
      crossRankActivationPushes_(0),
      boundaryInjectionWorkSamples_(),
      retirementEntryCount_(0), retirementEntryCharge_(0.0),
      retirementExitCount_(0), retirementExitResidualCharge_(0.0),
    peakRetirementCarriers_(0), stopReason_(),
    checkpointEnabled_(config.checkpoint.enabled()), nextCheckpointStep_(0),
    nextCheckpointTime_(0), retirementChangedOnRestart_(false),
    initialMaximumLabGamma_(1.0), initialGlobalParticleCount_(0),
    checkpointReader_()
  {
    if (communicator_ == MPI_COMM_NULL)
      throw std::invalid_argument("E/B solver communicator cannot be null");
    MPI_Comm_rank(communicator_, &rank_);
    MPI_Comm_size(communicator_, &size_);
    std::fill(performanceSeconds_,
              performanceSeconds_ + PerformanceStageCount, 0.0);
    // This runtime-only switch never changes the simulation or output cadence.
    const char* interval = std::getenv("APRL_PROFILE_INTERVAL");
    unsigned long long requestedInterval = 0;
    int invalidInterval = 0;
    if (interval && *interval)
      {
        char* end = NULL;
        errno = 0;
        requestedInterval = std::strtoull(interval, &end, 10);
        invalidInterval = interval[0] < '0' || interval[0] > '9' ||
          errno == ERANGE || !end || *end != '\0' ||
          requestedInterval > std::numeric_limits<std::size_t>::max();
      }
    int anyInvalidInterval = 0;
    MPI_Allreduce(&invalidInterval, &anyInvalidInterval, 1, MPI_INT, MPI_MAX,
                  communicator_);
    unsigned long long minimumInterval = 0, maximumInterval = 0;
    MPI_Allreduce(&requestedInterval, &minimumInterval, 1,
      MPI_UNSIGNED_LONG_LONG, MPI_MIN, communicator_);
    MPI_Allreduce(&requestedInterval, &maximumInterval, 1,
      MPI_UNSIGNED_LONG_LONG, MPI_MAX, communicator_);
    if (anyInvalidInterval || minimumInterval != maximumInterval)
      throw std::invalid_argument(
        "APRL_PROFILE_INTERVAL must be the same non-negative integer on all ranks");
    performanceInterval_ = static_cast<std::size_t>(requestedInterval);
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

  Simulation::~Simulation() = default;

  void Simulation::solve()
  {
    const char* gridAuditDirectory=std::getenv("APRL_GRID_WORK_AUDIT");
    const char* stageAuditDirectory=std::getenv("APRL_FIELD_STAGE_AUDIT");
    const char* geometryAuditDirectory=std::getenv("APRL_PARTICLE_WORK_GEOMETRY");
    const bool invalidGeometry=geometryAuditDirectory && (!*geometryAuditDirectory ||
       !config_.energyLedger.enabled || !config_.mesh.movingWindow.enabled ||
       !config_.mesh.selfFieldFeedback || config_.mesh.ballisticOutsideTransverse ||
       config_.particleRetirement.enabled ||
       config_.mesh.currentFilter!=EBCurrentFilter::None || checkpointEnabled_ ||
       !config_.checkpoint.restartFrom.empty());
    const int geometryFlags[2]={geometryAuditDirectory?1:0,invalidGeometry?1:0};int globalGeometryFlags[2]={};
    MPI_Allreduce(geometryFlags,globalGeometryFlags,2,MPI_INT,MPI_SUM,communicator_);
    if(globalGeometryFlags[1] || (globalGeometryFlags[0]!=0 && globalGeometryFlags[0]!=size_))
      throw std::invalid_argument("Geometry audit requires identical valid settings on all ranks, moving window, enabled ledger, feedback, and no filter/checkpoint/retirement/ballistic outside");
    if(config_.energyLedger.discreteOperatorAudit &&
       (!config_.energyLedger.enabled||!config_.mesh.movingWindow.enabled||
        config_.mesh.currentFilter!=EBCurrentFilter::None||checkpointEnabled_||
        !config_.checkpoint.restartFrom.empty()||config_.mesh.ballisticOutsideTransverse))
      throw std::invalid_argument("discrete_operator_audit requires moving window, enabled ledger, no filter/checkpoint/ballistic outside transverse");
    if(stageAuditDirectory && (!*stageAuditDirectory || !gridAuditDirectory))
      throw std::invalid_argument("APRL_FIELD_STAGE_AUDIT requires a directory and APRL_GRID_WORK_AUDIT");
    if(gridAuditDirectory && (!*gridAuditDirectory ||
        !config_.mesh.movingWindow.enabled || !config_.energyLedger.enabled ||
        config_.mesh.currentFilter!=EBCurrentFilter::None || checkpointEnabled_ ||
        !config_.checkpoint.restartFrom.empty()))
      throw std::invalid_argument("APRL_GRID_WORK_AUDIT requires a directory, moving window, ledger, unfiltered current and no checkpoint");
    const Double initializationStart = performanceNow();
    initialize();
    if(geometryAuditDirectory) {
      if(incident_)throw std::invalid_argument("Particle geometry audit does not support incident fields");
      particleWorkGeometry_.reset(new ParticleWorkGeometryWriter(geometryAuditDirectory,rank_));
      logRoot(communicator_,"[particle-work-geometry] read-only complete active pushes with private neighbor electric continuation; release/physical-boundary exclusions explicit; no state or deposition changes");
    }
    std::unique_ptr<InletConstraintProbe> inletConstraintProbe;
    if(const char* directory=std::getenv("APRL_INLET_CONSTRAINT_PROBE")) {
      if(!*directory||!config_.mesh.movingWindow.enabled)
        throw std::invalid_argument("Inlet constraint probe requires moving window and a directory");
      if(rank_==0)inletConstraintProbe.reset(new InletConstraintProbe(directory,localGeometry_));
    }
    std::unique_ptr<GridWorkAudit> gridWorkAudit;
    std::unique_ptr<FieldStageAudit> fieldStageAudit;
    std::unique_ptr<DiscreteFieldLedger> discreteLedger;
    if(config_.energyLedger.discreteOperatorAudit) {
      if(incident_)throw std::invalid_argument("discrete_operator_audit does not yet support incident waves");
      discreteLedger.reset(new DiscreteFieldLedger(config_.energyLedger.directory,
        localGeometry_,globalGeometry_.nz,localZOffset_,config_.boundary.cpml.cells[0],
        config_.boundary.cpml.cells[1],communicator_));
      movingWindow_->fields->setDiscreteObserver(&DiscreteFieldLedger::observe,discreteLedger.get());
      logRoot(communicator_,"[discrete-energy-ledger] enabled; independent Cowan/Yee curls, native energy and explicit numerical exchanges; existing physical ledger retained");
    }
    if(gridAuditDirectory) {
      gridWorkAudit.reset(new GridWorkAudit(gridAuditDirectory,rank_,
        config_.boundary.cpml.cells[0],config_.boundary.cpml.cells[1]));
      movingWindow_->fields->setBeforeInletObserver(
        &GridWorkAudit::beforeInlet,gridWorkAudit.get());
      logRoot(communicator_,"[grid-work-audit] enabled; read-only native/cell-centred J dot E, pre-inlet endpoint; no discrete Cowan conservation claim");
    }
    if(stageAuditDirectory) {
      fieldStageAudit.reset(new FieldStageAudit(stageAuditDirectory,rank_,
        globalGeometry_.nz,localZOffset_,config_.boundary.cpml.cells[0],config_.boundary.cpml.cells[1]));
      movingWindow_->fields->setStageObserver(&FieldStageAudit::observe,fieldStageAudit.get());
      logRoot(communicator_,"[field-stage-audit] enabled; read-only operation energy/flux/source trace");
    }
    if (performanceInterval_)
      {
        const Double elapsed = MPI_Wtime() - initializationStart;
        Double maximum = 0.0, sum = 0.0;
        MPI_Reduce(&elapsed, &maximum, 1, MPI_DOUBLE, MPI_MAX, 0, communicator_);
        MPI_Reduce(&elapsed, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, communicator_);
        if (rank_ == 0)
          {
            std::ostringstream message;
            message << std::setprecision(10)
                    << "[performance] phase=initialization mean_s=" << sum / size_
                    << " max_s=" << maximum;
            logRoot(communicator_, message.str());
          }
      }
    if (config_.checkpoint.restartFrom.empty() && config_.trajectory.enabled)
      sampleTrajectory();
    if (config_.checkpoint.restartFrom.empty() && detectors_)
      detectors_->sampleFieldPlanes(*fields_, timeBoxSI_);
    if (config_.checkpoint.restartFrom.empty() && config_.energyLedger.enabled)
      sampleEnergyLedger(true);
    if (checkpointEnabled_) maybeWriteCheckpoint();
    loopWallStart_ = MPI_Wtime();
    performanceWallStart_ = loopWallStart_;
    performanceLastStep_ = step_;

    while (timeBoxSI_ < totalTimeBoxSI_)
      {
        Double stageStart = performanceNow();
        const bool requestedStop =
          (step_ % config_.runtime.stopCheckIntervalSteps == 0) &&
          synchronizedStopRequested();
        performanceSeconds_[ControlStage] += performanceNow() - stageStart;
        if (requestedStop) break;

        if(gridWorkAudit) gridWorkAudit->begin(*fields_,step_);
        if(particleWorkGeometry_)particleWorkGeometry_->beginStep();
        if(discreteLedger)discreteLedger->begin(*fields_,step_);
        if(fieldStageAudit)fieldStageAudit->begin(*fields_,step_,particles_.size(),pmlCarriers_.size());

        if (config_.mesh.movingWindow.enabled)
          {
            stageStart = performanceNow();
            admitMovingWindowParticles();
            performanceSeconds_[PendingAdmissionStage] += performanceNow() - stageStart;
          }
        stageStart = performanceNow();
        Double communicationStart = halo_->communicationSeconds();
        fields_->clearCurrent();
        if(inletConstraintProbe)inletConstraintProbe->begin(*fields_,particles_,localOriginBox_,step_);
        particleBoundary_->beginStep();
        if (config_.mesh.movingWindow.enabled)
          {
            const EBElectricHaloView physicalHalo = beginMovingWindowStep();
            halo_->advanceMagnetic(*fields_, physicalHalo);
          }
        else
          halo_->advanceMagnetic(*fields_);
        if (incident_)
          incident_->correctAfterMagneticUpdate(
            *fields_, sources_, localOriginBox_, timeBoxSI_, frame_);
        if (config_.mesh.movingWindow.enabled)
          {
            const Double sampleCommunicationStart =
              movingWindow_->particleHalo->communicationSeconds();
            movingWindow_->particleHalo->exchangeAndInstall(*fields_);
            performanceSeconds_[FieldCommunicationStage] +=
              movingWindow_->particleHalo->communicationSeconds() -
              sampleCommunicationStart;
          }
        performanceSeconds_[FieldStage] += performanceNow() - stageStart;
        performanceSeconds_[FieldCommunicationStage] +=
          halo_->communicationSeconds() - communicationStart;

        if(particleWorkGeometry_)particleWorkGeometry_->exchangeFields(*fields_,communicator_);
        pushDepositAndMigrate();
        if(particleWorkGeometry_)particleWorkGeometry_->finishStep(step_+1);
        if(config_.mesh.movingWindow.referenceSurfaceCurrent)
          applyMovingTransverseReferenceCurrent();
        if(fieldStageAudit)fieldStageAudit->populations(particles_.size(),pmlCarriers_.size());
        if (config_.mesh.movingWindow.enabled)
          movingWindow_->particleHalo->clear(*fields_);
        stageStart = performanceNow();
        communicationStart = halo_->communicationSeconds();
        halo_->advanceElectric(*fields_, config_.mesh.currentFilter);
        if (incident_)
          incident_->correctAfterElectricUpdate(
            *fields_, sources_, localOriginBox_, timeBoxSI_, frame_);
        if(config_.mesh.movingWindow.referenceSurfaceCurrent)
          finishMovingTransverseReferenceCurrent();
        if(inletConstraintProbe)inletConstraintProbe->finish(*fields_,particles_,localOriginBox_);
        if(gridWorkAudit) gridWorkAudit->complete(*fields_);
        if(discreteLedger)discreteLedger->finishElectric(*fields_);
        if(fieldStageAudit)fieldStageAudit->sample("electric_mpi",*fields_);
        performanceSeconds_[FieldStage] += performanceNow() - stageStart;
        performanceSeconds_[FieldCommunicationStage] +=
          halo_->communicationSeconds() - communicationStart;

        stageStart = performanceNow();
        if (config_.energyLedger.enabled)
          {
            finalizeBoundaryInjectionWorkStep();
            if (movingWindow_) finalizeMovingConstraintWork();
            advanceEnergyLedgerStep();
          }

        timeBoxSI_ += globalGeometry_.dt;
        ++step_;
        if (config_.mesh.movingWindow.enabled)
          {
            const Double movingStart = performanceNow();
            finishMovingWindowStep();
            if(discreteLedger)discreteLedger->afterRemap(*fields_);
            if(fieldStageAudit) {
              fieldStageAudit->populations(particles_.size(),pmlCarriers_.size());
              fieldStageAudit->sample("remap",*fields_);
            }
            performanceSeconds_[MovingWindowStage] += performanceNow() - movingStart;
          }
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
        performanceSeconds_[DiagnosticsStage] += performanceNow() - stageStart;
        stageStart = performanceNow();
        const bool reached = configuredStopReached();
        performanceSeconds_[ControlStage] += performanceNow() - stageStart;
        reportPerformance();
        if (checkpointEnabled_) maybeWriteCheckpoint();
        if (reached)
          {
            configuredStopReached_ = true;
            break;
          }
      }

    if(particleWorkGeometry_)particleWorkGeometry_->close();
    if(inletConstraintProbe)inletConstraintProbe->close();
    if(gridWorkAudit) {
      movingWindow_->fields->setBeforeInletObserver(nullptr,nullptr);
      gridWorkAudit->close();
    }
    if(discreteLedger) {
      movingWindow_->fields->setDiscreteObserver(nullptr,nullptr);
      discreteLedger->close(!interrupted_&&configuredStopReached_);
    }
    if(fieldStageAudit) {
      movingWindow_->fields->setStageObserver(nullptr,nullptr);
      fieldStageAudit->close();
    }

    reportPerformance(true);

    logTeardownTrace("final-stop-sync-enter");
    armShutdownWatchdog("final-stop-sync");
    if (!interrupted_) synchronizedStopRequested();
    disarmShutdownWatchdog();
    logTeardownTrace("final-stop-sync-complete");
    const bool completed = !interrupted_ && configuredStopReached_;
    logTeardownTrace("final-energy-sample-enter");
    if (config_.energyLedger.enabled) sampleEnergyLedger(true, true);
    logTeardownTrace("final-energy-sample-complete");

    logTeardownTrace("trajectory-close-enter");
    finalizeTrajectoryOutput(completed);
    logTeardownTrace("trajectory-close-complete");

    logTeardownTrace("detector-close-enter");
    if (detectors_) detectors_->close(completed);
    logTeardownTrace("detector-close-complete");

    logTeardownTrace("ledger-close-enter");
    finalizeEnergyLedger(completed);
    logTeardownTrace("ledger-close-complete");

    logTeardownTrace("particle-boundary-report-enter");
    armShutdownWatchdog("particle-boundary-report");
    reportParticleBoundaryLosses();
    if (config_.mesh.movingWindow.enabled) finalizeMovingWindow();
    disarmShutdownWatchdog();
    logTeardownTrace("particle-boundary-report-complete");

    logTeardownTrace("final-resource-report-enter");
    armShutdownWatchdog("final-resource-report");
    reportResourceProgress("final");
    disarmShutdownWatchdog();
    logTeardownTrace("final-resource-report-complete");
    if (configuredStopReached_ && rank_ == 0)
      logRoot(communicator_, "Configured stop reached: " + stopReason_);
    if (interrupted_ && rank_ == 0)
      logRoot(communicator_,
        "Direct E/B run stopped after a complete field step; committed trajectory and detector outputs remain readable.");
    if (!interrupted_ && !configuredStopReached_ &&
        timeBoxSI_ >= totalTimeBoxSI_)
      {
        // The reductions above are not an output-close barrier: non-root
        // ranks may return while rank zero is still flushing its writers.
        // Complete every local close before this expected exception can
        // reach main's MPI_Abort. Genuine rank-local failures still abort
        // immediately through the ordinary exception handler.
        logTeardownTrace("duration-exhaustion-readiness-sync-enter");
        armShutdownWatchdog("duration-exhaustion-readiness-sync");
        MPI_Barrier(communicator_);
        disarmShutdownWatchdog();
        logTeardownTrace("duration-exhaustion-readiness-sync-complete");
        throw std::runtime_error(
          "mesh.duration was exhausted before the configured physical stop; trajectory and detector outputs are readable but marked incomplete");
      }
    logTeardownTrace("solve-complete");
  }

  void Simulation::initialize()
  {
    if (!(config_.mesh.boostGamma >= 1.0) ||
        !std::isfinite(config_.mesh.boostGamma))
      throw std::invalid_argument(
        "Direct E/B solver requires boost_gamma >= 1");

    preflightOutputPaths();
    initializeSources();

    if (sources_.maxwellIncidentWaveCount() > 0 &&
        config_.boundary.type == EBBoundaryType::Cpml)
      throw std::runtime_error(
        "CPML is active, but TF/SF seed injection has not yet been placed inside the CPML interior. This stage intentionally validates the no-seed absorbing boundary first; refusing an overlapping boundary configuration.");
    if (config_.mesh.fieldSolver == EBMaxwellSolver::CowanZ &&
        sources_.maxwellIncidentWaveCount() > 0)
      throw std::runtime_error(
        "Cowan-z field advance is active, but the present TF/SF incident-wave correction is still Yee-specific. Select mesh.field_solver: yee for this run until the generalized Cowan TF/SF stencil is connected; refusing to inject a numerically inconsistent seed field.");

    if (config_.checkpoint.restartFrom.empty())
      initializeParticles();
    else
      loadCheckpointGeometry();
    if (config_.mesh.movingWindow.enabled) initializeMovingWindow();
    runMetadata_ = initializeRunMetadata(config_, communicator_,
      hasLongitudinalLayoutReport_ ? &longitudinalLayoutReport_ : NULL);
    validateAndReportRadiationResolution();
    if (config_.checkpoint.restartFrom.empty()) reportDurationEstimate();
    fields_.reset(new EBFieldGrid(localGeometry_));
    halo_.reset(new EBZSlabHaloExchange(communicator_));
    particleBoundary_.reset(new ParticleOpenBoundary(
      localGeometry_, localOriginBox_, globalGeometry_, globalOriginBox_,
      localZOffset_));
    if (config_.mesh.movingWindow.enabled)
      installMovingWindowBoundary();
    else if (config_.boundary.type == EBBoundaryType::Cpml)
      fields_->setBoundary(std::unique_ptr<EBBoundaryOperator>(
        new EBConvolutionalPML(
          localGeometry_, globalGeometry_.nz, localZOffset_,
          config_.boundary.cpml,
          halo_->lowerRank() == MPI_PROC_NULL,
          halo_->upperRank() == MPI_PROC_NULL)));
    else
      halo_->installPhysicalBoundaryMask(*fields_);
    if (config_.checkpoint.restartFrom.empty())
      {
        calibrateResourceEstimate();
        fields_->clearFields();
        if (config_.mesh.movingWindow.enabled)
          initializeMovingWindowFields();
        else
          initializeParticleSelfField();
      }
    else
      restoreCheckpointFields();

    if (sources_.maxwellIncidentWaveCount() > 0)
      {
        const EBTFSFRegion region(
          1, globalGeometry_.nx - 1,
          1, globalGeometry_.ny - 1,
          1, globalGeometry_.nz - 1);
        incident_.reset(new EBMaxwellIncidentInjector(
          localGeometry_, region, localZOffset_, globalGeometry_.nz));
        if (config_.checkpoint.restartFrom.empty())
          incident_->initialize(*fields_, sources_, localOriginBox_,
                                timeBoxSI_, frame_, false);
      }

    initializeTrajectoryOutput();
    initializeDetectorOutput();
    initializeEnergyLedger();
    if (!config_.checkpoint.restartFrom.empty())
      restoreCheckpointDiagnostics();
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
            if (config_.mesh.currentFilter ==
                EBCurrentFilter::CompensatedBilinear)
              logRoot(communicator_,
                "Experimental current filter: compensated-bilinear is active on Jx/Jy/Jz after charge-conserving deposition and MPI current summation. The separable one-axis response is H=1-sin^4(k*cell_size/2), so Nyquist current is removed while the long-wave error begins at fourth order. Each axis pass leaves its two physical endpoint samples unchanged. Initial self-field charge is not filtered, and particle-field gather is not adjoint-filtered; energy-closure and radiation-band convergence must therefore be compared explicitly against current_filter: none.");
            else
              logRoot(communicator_,
                "Current filter: none; deposited J enters Ampere's law unchanged.");
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
        if (config_.mesh.ballisticOutsideTransverse)
          logRoot(communicator_,
            "EXPERIMENTAL boundary injection: particles outside the physical "
            "aperture are retained and drift without grid or prescribed "
            "fields, including through CPML/outside the box; trajectories "
            "are current-deposited only on clipped in-aperture segments. "
            "Incoming charge is an approximate boundary source, so global "
            "Gauss/energy closure is not guaranteed; compare aperture and "
            "initialization scans before interpreting radiation.");
        else
          logRoot(communicator_,
            "Outer particle fallback: charge-conserving open absorption remains active on non-CPML faces and for residual carrier cleanup; production mode retains only face totals.");
        if (particleCPML_->enabled() &&
            !config_.mesh.ballisticOutsideTransverse)
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
        if (config_.reference.frozenLabPlaneInjection)
          {
            std::ostringstream message;
            message << std::setprecision(10)
                    << "EXPERIMENTAL frozen lab-plane injection: "
                       "activation_plane_z="
                    << config_.reference.activationPlaneZ << " m; ";
            if (config_.reference.holdTransverseUntilActivation)
              message << "pre-activation x/y are held at their HDF5 "
                         "input-plane values and only z advances using "
                         "u'_z/sqrt(1+u'_z^2). Full transverse momentum "
                         "starts at the activation plane; the abrupt "
                         "current change can radiate. ";
            else
              message << "HDF5-v4 events are synchronized by full 3-D "
                         "ballistic motion; particles drift in x/y/z "
                         "without fields before the activation plane. ";
            message << "Particles remain Maxwell sources inside the "
                       "physical aperture. The crossing step is split at "
                       "the plane and its remainder uses the regular "
                       "Boris push. Check post-plane orbit and field-energy "
                       "convergence before production use.";
            logRoot(communicator_, message.str());
          }
        if (config_.reference.frozenLabPlaneInjection &&
            config_.mesh.ballisticOutsideTransverse)
          logRoot(communicator_,
            "WARNING: frozen particles outside the transverse physical aperture remain in the particle list but do not deposit Maxwell current there. This retains the finite-window boundary approximation for those particles.");
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
        for (const auto& plane : config_.detectors.energyPlanes)
          paths.push_back(joinPath(config_.detectors.directory,
            plane.name + "-energy.h5"));
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
        if (config_.mesh.currentFilter ==
            EBCurrentFilter::CompensatedBilinear)
          {
            const Double filterGain = 1.0 - std::pow(
              std::sin(3.1415926535897932384626433832795 / gridPoints),
              4.0);
            std::ostringstream filterMessage;
            filterMessage << std::setprecision(10)
              << "Current-filter target-band response (forward on-axis "
                 "deposited source): amplitude gain=" << filterGain
              << ", squared-gain proxy=" << filterGain * filterGain
              << " at " << photonEnergyLab
              << " eV. This analytical source response is not a guarantee "
                 "of nonlinear FEL power or bunching accuracy.";
            logRoot(communicator_, filterMessage.str());
          }
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
    std::vector<RelativisticParticleSI> initialFieldParticles;
    const std::vector<RelativisticParticleSI>* initialSource = &particles_;
    const bool includeCpml =
      config_.initialSelfField.includeTransverseCpmlParticles;
    if (config_.mesh.ballisticOutsideTransverse)
      {
        /* The legacy solve requires a complete CIC cell between every source
         * and the inner CPML potential boundary. In the opt-in solve, x/y
         * sources may occupy the whole mesh including CPML; its temporary
         * Poisson grid is padded outside the outer mesh surface. The z
         * boundary retains the existing one-cell protection. */
        const Double spacing[3] = {
          globalGeometry_.dx, globalGeometry_.dy, globalGeometry_.dz
        };
        unsigned long long localCpmlIncluded = 0;
        for (std::size_t index = 0; index < particles_.size(); ++index)
          {
            bool safe = true;
            for (unsigned int axis = 0; axis < 3; ++axis)
              {
                const Double lower = particleCPML_->physicalLower(axis);
                const Double upper = particleCPML_->physicalUpper(axis);
                if (includeCpml && axis < 2)
                  safe = safe && particles_[index].position[axis] >=
                    lower - config_.boundary.cpml.cells[axis] * spacing[axis] &&
                    particles_[index].position[axis] <=
                    upper + config_.boundary.cpml.cells[axis] * spacing[axis];
                else
                  safe = safe && particles_[index].position[axis] >
                    lower + spacing[axis] &&
                    particles_[index].position[axis] < upper - spacing[axis];
              }
            if (safe)
              {
                initialFieldParticles.push_back(particles_[index]);
                if (includeCpml &&
                    (particles_[index].position[0] <
                       particleCPML_->physicalLower(0) ||
                     particles_[index].position[0] >
                       particleCPML_->physicalUpper(0) ||
                     particles_[index].position[1] <
                       particleCPML_->physicalLower(1) ||
                     particles_[index].position[1] >
                       particleCPML_->physicalUpper(1)))
                  ++localCpmlIncluded;
              }
          }
        const unsigned long long localIncluded =
          static_cast<unsigned long long>(initialFieldParticles.size());
        unsigned long long globalIncluded = 0;
        MPI_Allreduce(&localIncluded, &globalIncluded, 1,
          MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator_);
        unsigned long long globalCpmlIncluded = 0;
        MPI_Allreduce(&localCpmlIncluded, &globalCpmlIncluded, 1,
          MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator_);
        if (rank_ == 0)
          {
            std::ostringstream message;
            message << "Boundary-injection initial Poisson source: "
                    << globalIncluded << " particles";
            if (includeCpml)
              message << " inside the full transverse mesh, including "
                      << globalCpmlIncluded << " in the transverse CPML; "
                         "all farther-out particles are omitted. The CPML "
                         "evolution fields and memories still start at zero.";
            else
              message << " at least one complete cell inside the physical "
                         "aperture. All other particles are omitted from "
                         "the initial field.";
            logRoot(communicator_, message.str());
          }
        if (globalIncluded == 0)
          {
            if (rank_ == 0)
              logRoot(communicator_,
                "WARNING: no particles qualify for the initial Poisson "
                "source; the self-field starts at zero and subsequent "
                "boundary injection is approximate.");
            return;
          }
        initialSource = &initialFieldParticles;
      }
    std::size_t staticFieldGuardCells[3] = {0, 0, 0};
    if (config_.boundary.type == EBBoundaryType::Cpml)
      for (unsigned int axis = 0; axis < 3; ++axis)
        staticFieldGuardCells[axis] = config_.boundary.cpml.cells[axis];
    const EBGaussInitializationReport report =
      EBGaussFieldInitializer::initialize(
        *fields_, globalGeometry_, localZOffset_, localOriginBox_,
        *initialSource,
        config_.initialSelfField.model,
        config_.initialSelfField.maximumMeanTransverseBeta,
        staticFieldGuardCells,
        config_.initialSelfField.relativeTolerance,
        config_.initialSelfField.maximumIterations, communicator_,
        includeCpml,
        config_.reference.holdTransverseUntilActivation, &frame_);
    if (rank_ == 0)
      {
        std::ostringstream message;
        message << std::setprecision(10)
          << "Initial particle self-field: "
          << (config_.initialSelfField.model ==
                InitialSelfFieldModel::RelativisticPoisson ?
                "relativistic-Poisson rigid-beam" :
              config_.initialSelfField.model ==
                InitialSelfFieldModel::CowanRelativisticPoisson ?
                "Cowan relativistic-Poisson quasi-steady rigid-beam" :
                "legacy electrostatic-Poisson")
          << " construction "
          << "converged in " << report.iterations
          << " CG iterations; post-check relative L2 residual="
          << report.relativeResidual
          << ", maximum absolute residual="
          << report.maximumGaussResidual << " V/m^2, total charge="
          << report.totalCharge << " C, mean beta_box_z="
          << report.meanBetaZ << ", mean beta_box_transverse="
          << report.meanBetaTransverse
          << ", maximum mean beta_box_transverse="
          << config_.initialSelfField.maximumMeanTransverseBeta
          << ", initial electric energy="
          << report.initialElectricEnergy << " J, initial magnetic energy="
          << report.initialMagneticEnergy
          << " J, static zero-potential guard cells=("
          << report.staticFieldGuardCells[0] << ", "
          << report.staticFieldGuardCells[1] << ", "
          << report.staticFieldGuardCells[2]
          << "), transverse Poisson padding cells="
          << report.transversePoissonPaddingCells
          << ", bunch RMS cells=("
          << report.rmsPositionCells[0] << ", "
          << report.rmsPositionCells[1] << ", "
          << report.rmsPositionCells[2] << "), centre-to-static-boundary "
             "padding/RMS=("
          << report.centrePaddingRms[0] << ", "
          << report.centrePaddingRms[1] << ", "
          << report.centrePaddingRms[2]
          << "), temporary scalar memory/rank<="
          << static_cast<Double>(report.temporaryBytes) /
             (1024.0 * 1024.0) << " MiB including transient CIC charge "
             "lattices. The potential is now released; "
             "only SI E/B remains. The reported static zero-potential "
             "boundary still requires padding convergence. ";
        if (config_.initialSelfField.model ==
            InitialSelfFieldModel::RelativisticPoisson)
          message << "The relativistic model represents a rigid "
                     "common-velocity bunch; velocity-spread and "
                     "kinetic-equilibrium errors require separate "
                     "convergence tests.";
        else if (config_.initialSelfField.model ==
                 InitialSelfFieldModel::CowanRelativisticPoisson)
          message << "The Cowan model matches the spatial Cowan operator "
                     "for a common axial velocity; finite-step translation, "
                     "velocity-spread and kinetic-equilibrium errors remain "
                     "approximations. Particle positions, momenta and held "
                     "release conditions are unchanged.";
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
    if (config_.mesh.movingWindow.enabled)
      {
        estimatedStepSeconds_ = 0.0;
        if (rank_ == 0)
          logRoot(communicator_,
            "[resource] phase=pre-run-calibration status=skipped reason=moving-window-state");
        return;
      }
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
        halo_->advanceElectric(*fields_, config_.mesh.currentFilter);
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
          if (config_.mesh.ballisticOutsideTransverse &&
              !particleCPML_->containsPhysical(start))
            RelativisticBorisPusher::pushPosition(
              particle.position, particle.properVelocity,
              globalGeometry_.dt);
          else
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
          if (local &&
              (!config_.mesh.ballisticOutsideTransverse ||
               (particleCPML_->containsPhysical(start) &&
                particleCPML_->containsPhysical(particle.position))))
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
    modeled += static_cast<long double>(ParticleWorkspace::PreparationBatchCapacity) *
      sizeof(ParticleWorkspace::Prepared);
    if (movingWindow_ && movingWindow_->pending)
      modeled += movingWindow_->pending->memoryBytes();
    if (movingWindow_ && movingWindow_->pendingStatistics)
      modeled += movingWindow_->pendingStatistics->storageBytes();
    if (movingWindow_)
      modeled += static_cast<long double>(movingWindow_->constraintWorkSamples.capacity()) *
        sizeof(SimulationMovingWindowState::ConstraintWorkSample);
    if (config_.mesh.ballisticOutsideTransverse &&
        config_.energyLedger.enabled)
      modeled += static_cast<long double>(particles_.size()) *
        sizeof(BoundaryInjectionWorkSample);

    const long double fieldPointBytes = 6.0L * sizeof(Double);
    const long double planePoints =
      static_cast<long double>(globalGeometry_.nx) *
      static_cast<long double>(globalGeometry_.ny);
    if (!detectorConfig_.fieldPlanes.empty() || !detectorConfig_.energyPlanes.empty())
      modeled += planePoints * fieldPointBytes;
    if (config_.energyLedger.enabled)
      modeled += static_cast<long double>(ParticleWorkspace::PreparationBatchCapacity) *
        particleSubsteps_ * sizeof(ParticleEnergyTrace::Segment);
    if (!detectorConfig_.energyPlanes.empty())
      modeled += static_cast<long double>(ParticleWorkspace::PreparationBatchCapacity) *
        particleSubsteps_ * detectorConfig_.energyPlanes.size() * sizeof(EnergyPlaneCrossing);
    if (rank_ == 0)
      for (const auto& plane : detectorConfig_.energyPlanes)
        modeled += static_cast<long double>(plane.bufferSamples) * 1024.L;
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
      static_cast<unsigned long long>(storedParticleCount());
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
    const std::size_t segmentSteps = step_ - config_.checkpoint.restartStep;
    const Double secondsPerStep = segmentSteps > 0 ?
      elapsed / static_cast<Double>(segmentSteps) : 0.0;
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
            message << " The resolved physical extent is always computed "
                       "from mesh.cells times mesh.cell_size; z may be an "
                       "explicit count or the resolved observation-bounded "
                       "automatic count.";
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
        if (config_.mesh.longitudinal.automatic() &&
            config_.initialSelfField.enabled &&
            config_.mesh.longitudinal.selfFieldPaddingRms < 4.0)
          {
            std::ostringstream warning;
            warning << std::setprecision(8)
              << "WARNING: automatic longitudinal initial self-field uses "
                 "self_field_padding_rms="
              << config_.mesh.longitudinal.selfFieldPaddingRms
              << ". The configured Poisson solve still initializes the "
                 "complete resolved physical interval and enforces the "
                 "discrete Gauss constraint, but its zero-potential surface "
                 "is the current inner CPML boundary rather than a guaranteed "
                 "four-RMS boundary. Treat this as a box-clipped finite-domain "
                 "initial field and scan lower_guard_cells/upper_guard_cells "
                 "or compare a larger explicit box before interpreting "
                 "absolute self-field energy.";
            logRoot(communicator_, warning.str());
          }
      }

    SILabPlaneProjectionReport planeProjection;
    if (config_.reference.frozenLabPlaneInjection &&
        !inputReport.laboratoryPlaneTimeCoordinates)
      throw std::invalid_argument(
        "beam.reference.frozen_lab_plane_injection requires HDF5 format 4 fixed-plane arrival-time records");
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
              config_.beam.positionOffset[2], communicator_,
              config_.reference.holdTransverseUntilActivation);
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
      config_.mesh.boostGamma * config_.mesh.cellSize[2];

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
                             "input plane upstream, move/redefine the magnetic "
                             "interaction region downstream, reduce the bunch "
                             "longitudinal span, or refine dz to reduce the "
                             "one-cell safety margin";
              }
            message << ".";
            throw std::runtime_error(message.str());
          }
      }

    /* The legacy mode anchors box time zero at the reconstructed front.
     * Frozen-plane injection instead starts before every input-plane event.
     * Either way, the large synchronization span is not a real beamline
     * drift; the frozen mode keeps the input-plane transverse coordinates. */
    if (config_.mesh.ballisticOutsideTransverse)
      {
        unsigned long long localLabOutside = 0;
        for (std::size_t index = 0; index < particles_.size(); ++index)
          {
            bool outside = false;
            for (unsigned int axis = 0; axis < 2; ++axis)
              {
                const Double lower = config_.mesh.center[axis] -
                  0.5 * config_.mesh.extent[axis] +
                  config_.boundary.cpml.cells[axis] *
                  config_.mesh.cellSize[axis];
                const Double upper = config_.mesh.center[axis] +
                  0.5 * config_.mesh.extent[axis] -
                  config_.boundary.cpml.cells[axis] *
                  config_.mesh.cellSize[axis];
                outside = outside ||
                  particles_[index].position[axis] < lower ||
                  particles_[index].position[axis] > upper;
              }
            if (outside) ++localLabOutside;
          }
        unsigned long long globalLabOutside = 0;
        MPI_Allreduce(&localLabOutside, &globalLabOutside, 1,
          MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator_);
        if (globalLabOutside && rank_ == 0)
          {
            std::ostringstream message;
            message << "WARNING: "
                    << "reconstructed initial lab snapshot has "
                    << globalLabOutside << " particles outside the "
                       "transverse physical aperture; boundary injection "
                       "is enabled, so this is not fatal.";
            logRoot(communicator_, message.str());
          }
      }
    Double initialFrameTimeLab = 0.0;
    Double initialFrameZLab = placementReport.headAfterLab;
    if (config_.reference.frozenLabPlaneInjection)
      {
        Double localEarliestCrossing =
          std::numeric_limits<Double>::infinity();
        int localInvalidHeldVelocity = 0;
        BoostFrameTransform probeFrame;
        if (config_.reference.holdTransverseUntilActivation)
          probeFrame.setOriginsFromGamma(
            config_.mesh.boostGamma, SI::c, 0.0,
            config_.reference.activationPlaneZ, 0.0);
        for (std::size_t index = 0; index < particles_.size(); ++index)
          {
            const RelativisticParticleSI& particle = particles_[index];
            const Double gamma =
              BoostFrameTransform::gammaFromProperVelocity(
                particle.properVelocity);
            const Double betaZ = particle.properVelocity[2] / gamma;
            /* placeRelativeLabSnapshot() has already translated z into
             * the absolute laboratory coordinate system. */
            const Double snapshotZ = particle.position[2];
            if (config_.reference.holdTransverseUntilActivation)
              {
                FieldVector<Double> boxProperVelocity(0.0);
                probeFrame.properVelocityLabToBox(
                  particle.properVelocity, boxProperVelocity);
                const Double betaHeld = probeFrame.heldLongitudinalBetaZ(
                  boxProperVelocity);
                const Double labAdvance = probeFrame.beta() + betaHeld;
                if (!(labAdvance > 0.0) || !std::isfinite(labAdvance))
                  {
                    localInvalidHeldVelocity = 1;
                    continue;
                  }
                const Double inputTimeLab =
                  (config_.reference.inputPlaneZ - snapshotZ) /
                  (betaZ * SI::c);
                const Double inputToActivation =
                  config_.reference.inputPlaneZ -
                  config_.reference.activationPlaneZ;
                const Double crossingOriginTime = inputTimeLab -
                  inputToActivation *
                  (1.0 + probeFrame.beta() * betaHeld) /
                  (SI::c * labAdvance);
                const Double twoCellMargin =
                  2.0 * config_.mesh.cellSize[2] /
                  (probeFrame.gamma() * SI::c * labAdvance);
                localEarliestCrossing = std::min(localEarliestCrossing,
                  crossingOriginTime - twoCellMargin);
              }
            else
              localEarliestCrossing = std::min(localEarliestCrossing,
                (config_.reference.activationPlaneZ - snapshotZ) /
                (betaZ * SI::c));
          }
        int invalidHeldVelocity = 0;
        MPI_Allreduce(&localInvalidHeldVelocity, &invalidHeldVelocity,
          1, MPI_INT, MPI_MAX, communicator_);
        if (invalidHeldVelocity)
          throw std::invalid_argument(
            "Held-transverse pseudo-trajectory must advance in laboratory z for every particle");
        Double earliestCrossing = 0.0;
        MPI_Allreduce(&localEarliestCrossing, &earliestCrossing,
          1, MPI_DOUBLE, MPI_MIN, communicator_);
        /* The frame origin is two nominal z-Courant steps before the first
         * plane crossing. This makes every source frozen at box time zero. */
        initialFrameTimeLab =
          config_.reference.holdTransverseUntilActivation ?
            earliestCrossing :
            earliestCrossing -
              2.0 * config_.mesh.cellSize[2] /
              (config_.mesh.boostGamma * SI::c);
        initialFrameZLab = config_.reference.activationPlaneZ;
      }
    frame_.setOriginsFromGamma(config_.mesh.boostGamma, SI::c,
      initialFrameTimeLab, initialFrameZLab, 0.0);
    const SIBunchBoostReport localBoostReport =
      SIBunchPreprocessor::boostLabSnapshotToBoxTimeZero(
        particles_, frame_, 0.0,
        config_.reference.holdTransverseUntilActivation,
        config_.reference.inputPlaneZ);
    if (config_.reference.frozenLabPlaneInjection)
      {
        Double physicalLower[2];
        Double physicalUpper[2];
        Double meshLower[2];
        Double meshUpper[2];
        for (unsigned int axis = 0; axis < 2; ++axis)
          {
            meshLower[axis] = config_.mesh.center[axis] -
              0.5 * config_.mesh.extent[axis];
            meshUpper[axis] = config_.mesh.center[axis] +
              0.5 * config_.mesh.extent[axis];
            physicalLower[axis] = meshLower[axis] +
              config_.boundary.cpml.cells[axis] *
              config_.mesh.cellSize[axis];
            physicalUpper[axis] = meshUpper[axis] -
              config_.boundary.cpml.cells[axis] *
              config_.mesh.cellSize[axis];
          }
        unsigned long long localOutside[2] = {};
        for (std::size_t index = 0; index < particles_.size(); ++index)
          {
            bool outsidePhysical = false;
            bool outsideMesh = false;
            for (unsigned int axis = 0; axis < 2; ++axis)
              {
                const Double coordinate = particles_[index].position[axis];
                outsidePhysical = outsidePhysical ||
                  coordinate < physicalLower[axis] ||
                  coordinate > physicalUpper[axis];
                outsideMesh = outsideMesh ||
                  coordinate < meshLower[axis] ||
                  coordinate > meshUpper[axis];
              }
            if (outsidePhysical) ++localOutside[0];
            if (outsideMesh) ++localOutside[1];
          }
        unsigned long long globalOutside[2] = {};
        MPI_Allreduce(localOutside, globalOutside, 2,
          MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator_);
        if (globalOutside[0] && rank_ == 0)
          {
            std::ostringstream warning;
            warning << std::setprecision(6)
                    << "WARNING: frozen-mode boosted t=0 transverse "
                       "ballistic positions outside the physical aperture: "
                    << globalOutside[0] << "/" << globalCount << " ("
                    << 100.0 * static_cast<Double>(globalOutside[0]) /
                         static_cast<Double>(globalCount)
                    << "%); outside the full mesh: "
                    << globalOutside[1] << "/" << globalCount << " ("
                    << 100.0 * static_cast<Double>(globalOutside[1]) /
                         static_cast<Double>(globalCount)
                    << "%). Outside particles cannot be Maxwell sources "
                       "there; with ballistic_outside_transverse they remain "
                       "tracked until entry, otherwise normal particle "
                       "boundary rules apply.";
            logRoot(communicator_, warning.str());
          }
        int localAlreadyPast = 0;
        for (std::size_t index = 0; index < particles_.size(); ++index)
          {
            particles_[index].frozen = true;
            if (config_.reference.holdTransverseUntilActivation)
              particles_[index].frozenBetaZ = frame_.heldLongitudinalBetaZ(
                particles_[index].properVelocity);
            if (frame_.labZFromBoxZT(particles_[index].position[2], 0.0) >=
                config_.reference.activationPlaneZ)
              localAlreadyPast = 1;
          }
        int anyAlreadyPast = 0;
        MPI_Allreduce(&localAlreadyPast, &anyAlreadyPast, 1,
          MPI_INT, MPI_MAX, communicator_);
        if (anyAlreadyPast)
          throw std::runtime_error(
            "Frozen lab-plane injection requires every particle upstream of the input plane at box time zero");
      }
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

    Double boxTranslation = config_.mesh.center[2] -
      0.5 * (boostedMinimum + boostedMaximum);
    if (config_.mesh.longitudinal.automatic())
      {
        if (globalCount == 0)
          throw std::runtime_error(
            "Automatic longitudinal sizing requires at least one particle");

        /* Use represented charge, rather than record count or the optional
         * diagnostic input weight, for the self-field padding statistics. */
        long double localChargeWeight = 0.0L;
        long double localWeightedZ = 0.0L;
        for (std::size_t index = 0; index < particles_.size(); ++index)
          {
            const long double weight = std::abs(
              static_cast<long double>(particles_[index].charge));
            localChargeWeight += weight;
            localWeightedZ += weight *
              static_cast<long double>(particles_[index].position[2]);
          }
        long double globalChargeWeight = 0.0L;
        long double globalWeightedZ = 0.0L;
        MPI_Allreduce(&localChargeWeight, &globalChargeWeight, 1,
                      MPI_LONG_DOUBLE, MPI_SUM, communicator_);
        MPI_Allreduce(&localWeightedZ, &globalWeightedZ, 1,
                      MPI_LONG_DOUBLE, MPI_SUM, communicator_);
        if (!(globalChargeWeight > 0.0L))
          throw std::runtime_error(
            "Automatic longitudinal sizing found no represented charge");
        const long double meanWide =
          globalWeightedZ / globalChargeWeight;
        long double localVarianceNumerator = 0.0L;
        for (std::size_t index = 0; index < particles_.size(); ++index)
          {
            const long double weight = std::abs(
              static_cast<long double>(particles_[index].charge));
            const long double offset =
              static_cast<long double>(particles_[index].position[2]) -
              meanWide;
            localVarianceNumerator += weight * offset * offset;
          }
        long double globalVarianceNumerator = 0.0L;
        MPI_Allreduce(&localVarianceNumerator,
                      &globalVarianceNumerator, 1, MPI_LONG_DOUBLE,
                      MPI_SUM, communicator_);
        const Double boostedMean = static_cast<Double>(meanWide);
        const Double boostedRms = std::sqrt(std::max(0.0,
          static_cast<Double>(globalVarianceNumerator /
                              globalChargeWeight)));

        const Double observationTargetLab =
          config_.stop.mode == StopMode::ReferenceCenterZ ?
            config_.stop.referenceZ : lastBeamlineInteractionExitLab();
        /* The lab snapshot was synchronized onto t_box=0.  Its rightmost
         * record is therefore the causal emission event for the tight
         * forward-light bound, even if it occurs at a slightly different
         * lab time/position than the original lab-snapshot head anchor. */
        const Double initialRadiatingHeadLab = frame_.labZFromBoxZT(
          boostedMaximum, 0.0);
        const Double targetDistanceLab =
          observationTargetLab - initialRadiatingHeadLab;
        if (!(targetDistanceLab > 0.0) ||
            !std::isfinite(targetDistanceLab))
          {
            std::ostringstream message;
            message << std::setprecision(10)
                    << "Automatic longitudinal sizing requires its final "
                       "observation/stop plane to be downstream of the "
                       "initial bunch head: head_z="
                    << initialRadiatingHeadLab
                    << " m, target_z=" << observationTargetLab << " m.";
            throw std::runtime_error(message.str());
          }

        Double ballisticStopTime = 0.0;
        int localUnreachable = 0;
        if (config_.stop.mode == StopMode::ReferenceCenterZ)
          {
            const Double rawReferenceBoxZ = frame_.boxZFromLabZT(
              config_.reference.initialCenterZ, 0.0);
            const Double currentReferenceLab = frame_.labZFromBoxZT(
              rawReferenceBoxZ, 0.0);
            const Double remaining = observationTargetLab -
              currentReferenceLab;
            const Double labAdvancePerBoxSecond =
              frame_.gammaBeta() * SI::c;
            if (!(remaining > 0.0) ||
                !(labAdvancePerBoxSecond > 0.0) ||
                !std::isfinite(labAdvancePerBoxSecond))
              localUnreachable = 1;
            else
              ballisticStopTime = remaining / labAdvancePerBoxSecond;
          }
        else
          {
            Double localMaximumTime = 0.0;
            for (std::size_t index = 0; index < particles_.size(); ++index)
              {
                Double eventTimeLab = 0.0;
                Double eventZLab = 0.0;
                frame_.boxToLab(0.0, particles_[index].position[2],
                                eventTimeLab, eventZLab);
                const Double remaining = observationTargetLab - eventZLab;
                if (remaining <= 0.0) continue;
                const Double gammaBox =
                  BoostFrameTransform::gammaFromProperVelocity(
                    particles_[index].properVelocity);
                const Double betaBoxZ =
                  particles_[index].properVelocity[2] / gammaBox;
                const Double labAdvancePerBoxSecond =
                  frame_.gamma() * SI::c *
                  (frame_.beta() + betaBoxZ);
                if (!(labAdvancePerBoxSecond > 0.0) ||
                    !std::isfinite(labAdvancePerBoxSecond))
                  localUnreachable = 1;
                else
                  localMaximumTime = std::max(localMaximumTime,
                    remaining / labAdvancePerBoxSecond);
              }
            MPI_Allreduce(&localMaximumTime, &ballisticStopTime, 1,
                          MPI_DOUBLE, MPI_MAX, communicator_);
          }
        int globalUnreachable = 0;
        MPI_Allreduce(&localUnreachable, &globalUnreachable, 1, MPI_INT,
                      MPI_MAX, communicator_);
        if (globalUnreachable || !(ballisticStopTime > 0.0) ||
            !std::isfinite(ballisticStopTime))
          throw std::runtime_error(
            "Automatic longitudinal sizing cannot reach the configured stop plane from the initial boosted-frame particle kinematics");
        const Double plannedBallisticTime = ballisticStopTime *
          config_.mesh.longitudinal.ballisticTimeSafetyFactor;
        if (!(plannedBallisticTime > 0.0) ||
            !std::isfinite(plannedBallisticTime))
          throw std::runtime_error(
            "Automatic longitudinal ballistic safety time is not finite");

        Double localBallisticMinimum =
          std::numeric_limits<Double>::infinity();
        Double localBallisticMaximum =
          -std::numeric_limits<Double>::infinity();
        int localBallisticOverflow = 0;
        for (std::size_t index = 0; index < particles_.size(); ++index)
          {
            const Double start = particles_[index].position[2];
            const Double gammaBox =
              BoostFrameTransform::gammaFromProperVelocity(
                particles_[index].properVelocity);
            const Double betaBoxZ =
              particles_[index].properVelocity[2] / gammaBox;
            const Double finish = std::fma(
              betaBoxZ * SI::c, plannedBallisticTime, start);
            if (!std::isfinite(finish))
              {
                localBallisticOverflow = 1;
                continue;
              }
            localBallisticMinimum = std::min(localBallisticMinimum,
              std::min(start, finish));
            localBallisticMaximum = std::max(localBallisticMaximum,
              std::max(start, finish));
          }
        int globalBallisticOverflow = 0;
        MPI_Allreduce(&localBallisticOverflow,
                      &globalBallisticOverflow, 1, MPI_INT, MPI_MAX,
                      communicator_);
        if (globalBallisticOverflow)
          throw std::runtime_error(
            "Automatic longitudinal ballistic particle envelope overflowed");
        Double ballisticMinimum = 0.0;
        Double ballisticMaximum = 0.0;
        MPI_Allreduce(&localBallisticMinimum, &ballisticMinimum, 1,
                      MPI_DOUBLE, MPI_MIN, communicator_);
        MPI_Allreduce(&localBallisticMaximum, &ballisticMaximum, 1,
                      MPI_DOUBLE, MPI_MAX, communicator_);

        const Double radiationAdvance = forwardLightIntersectionAdvance(
          targetDistanceLab, config_.mesh.boostGamma);
        const Double radiationUpper = boostedMaximum + radiationAdvance;
        if (!std::isfinite(radiationUpper))
          throw std::runtime_error(
            "Automatic longitudinal forward-radiation envelope overflowed");

        LongitudinalLayoutInput layoutInput;
        layoutInput.initialMinimum = boostedMinimum;
        layoutInput.initialMaximum = boostedMaximum;
        layoutInput.initialMean = boostedMean;
        layoutInput.initialRms = boostedRms;
        layoutInput.ballisticMinimum = ballisticMinimum;
        layoutInput.ballisticMaximum = ballisticMaximum;
        layoutInput.radiationUpper = radiationUpper;
        layoutInput.dz = config_.mesh.cellSize[2];
        layoutInput.cpmlCells = config_.boundary.cpml.cells[2];
        layoutInput.leftGuardCells =
          config_.mesh.longitudinal.lowerGuardCells;
        layoutInput.rightGuardCells =
          config_.mesh.longitudinal.upperGuardCells;
        layoutInput.selfFieldPaddingRms =
          config_.initialSelfField.enabled ?
            config_.mesh.longitudinal.selfFieldPaddingRms : 0.0;
        layoutInput.alignmentCells =
          config_.mesh.longitudinal.alignmentCells;
        layoutInput.meshCenter = config_.mesh.center[2];

        longitudinalLayoutReport_.layout =
          deriveLongitudinalLayout(layoutInput);
        longitudinalLayoutReport_.observationTargetLab =
          observationTargetLab;
        longitudinalLayoutReport_.initialHeadLab =
          initialRadiatingHeadLab;
        longitudinalLayoutReport_.ballisticStopTime = ballisticStopTime;
        longitudinalLayoutReport_.plannedBallisticTime =
          plannedBallisticTime;
        longitudinalLayoutReport_.radiationAdvance = radiationAdvance;
        hasLongitudinalLayoutReport_ = true;
        config_.mesh.cells[2] =
          longitudinalLayoutReport_.layout.totalCells;
        config_.mesh.extent[2] =
          longitudinalLayoutReport_.layout.extent;
        boxTranslation = longitudinalLayoutReport_.layout.translation;
      }

    initializeGeometry();
    particleCPML_.reset(new ParticleCPMLRegion(
      globalGeometry_, globalOriginBox_,
      config_.boundary.type == EBBoundaryType::Cpml ?
        config_.boundary.cpml : EBCPMLParameters()));
    for (std::size_t index = 0; index < particles_.size(); ++index)
      particles_[index].position[2] += boxTranslation;
    boostedMinimum += boxTranslation;
    boostedMaximum += boxTranslation;
    frame_.setOriginsFromGamma(config_.mesh.boostGamma, SI::c,
      initialFrameTimeLab, initialFrameZLab, boxTranslation);
    referenceCenterBoxZ_ = frame_.boxZFromLabZT(
      config_.reference.initialCenterZ, 0.0);

    if (hasLongitudinalLayoutReport_ && rank_ == 0)
      {
        const LongitudinalLayout& layout =
          longitudinalLayoutReport_.layout;
        std::ostringstream message;
        message << std::setprecision(10)
                << "Automatic longitudinal mesh: mode=observation-bounded"
                << ", target_lab_z="
                << longitudinalLayoutReport_.observationTargetLab
                << " m, earliest-head-light advance="
                << longitudinalLayoutReport_.radiationAdvance
                << " m, ballistic_stop_s="
                << longitudinalLayoutReport_.ballisticStopTime
                << ", planned_ballistic_s="
                << longitudinalLayoutReport_.plannedBallisticTime
                << ", resolved_cells=" << layout.totalCells
                << ", extent=" << layout.extent
                << " m, translation=" << layout.translation
                << " m, physical_z=[" << layout.physicalLower << ", "
                << layout.physicalUpper << "] m, initial_z=["
                << layout.initialMinimum << ", "
                << layout.initialMaximum << "] m, ballistic_z=["
                << layout.ballisticMinimum << ", "
                << layout.ballisticMaximum << "] m, radiation_z_max="
                << layout.radiationUpper
                << " m, right_rounding_padding="
                << layout.rightRoundingPadding
                << " m. CPML lies outside the reported physical interval; "
                   "this layout preserves radiation only through the final "
                   "observation event, not the complete final 3-D field.";
        logRoot(communicator_, message.str());
      }

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
              << "HDF5 lab-plane reconstruction (v"
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
          << (config_.reference.frozenLabPlaneInjection ?
                "Frozen-plane Lorentz synchronization: anchor_lab_z=" :
                "Head-anchored Lorentz synchronization: anchor_lab_z=")
          << initialFrameZLab << " m, origin_lab_t="
          << initialFrameTimeLab << " s, virtual lab event times [s]=["
          << boostReport.earliestLabEventTime << " to "
          << boostReport.latestLabEventTime << "], maximum synchronization "
             "span=" << boostReport.maximumAbsoluteDriftTime << " s ("
          << SI::c * boostReport.maximumAbsoluteDriftTime
          << (config_.reference.holdTransverseUntilActivation ?
                " light-metres; held-transverse pseudo-trajectory), " :
              config_.reference.frozenLabPlaneInjection ?
                " light-metres; full 3-D ballistic synchronization), " :
                " light-metres; diagnostic only, not a required physical drift), ")
          << (config_.mesh.longitudinal.automatic() ?
                "automatically placed boosted z range [m]=[" :
                "centred boosted z range [m]=[")
          << boostedMinimum << ", "
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
        const auto considerParticle = [&](const RelativisticParticleSI& particle)
          {
            Double eventTimeLab = 0.0;
            Double eventZLab = 0.0;
            frame_.boxToLab(0.0, particle.position[2],
                            eventTimeLab, eventZLab);
            const Double remaining = interactionExit - eventZLab;
            if (remaining <= 0.0) return;
            const Double gammaBox =
              BoostFrameTransform::gammaFromProperVelocity(
                particle.properVelocity);
            const Double betaBoxZ =
              particle.properVelocity[2] / gammaBox;
            const Double labAdvancePerBoxSecond =
              frame_.gamma() * SI::c * (frame_.beta() + betaBoxZ);
            if (!(labAdvancePerBoxSecond > 0.0) ||
                !std::isfinite(labAdvancePerBoxSecond))
              localUnreachable = 1;
            else
              localMaximumTime = std::max(localMaximumTime,
                remaining / labAdvancePerBoxSecond);
          };
        for (const auto& particle : particles_) considerParticle(particle);
        if (movingWindow_ && movingWindow_->pending)
          for (const auto& particle : movingWindow_->pending->initialParticles())
            considerParticle(particle);
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
    if (config_.checkpoint.restartFrom.empty())
      {
    const unsigned long long localParticleCount =
      static_cast<unsigned long long>(particles_.size());
    MPI_Allreduce(&localParticleCount, &initialGlobalParticleCount_, 1,
                  MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator_);

    Double localMaximumGamma = 1.0;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      localMaximumGamma = std::max(localMaximumGamma,
        BoostFrameTransform::gammaFromProperVelocity(
          particles_[index].properVelocity));
    MPI_Allreduce(&localMaximumGamma, &initialMaximumLabGamma_, 1, MPI_DOUBLE,
                  MPI_MAX, communicator_);
      }
    const Double maximumGamma = initialMaximumLabGamma_;
    const unsigned long long globalParticleCount = initialGlobalParticleCount_;
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
              config_.mesh.boostGamma * config_.mesh.cellSize[2];
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
    const Double marginLab =
      config_.mesh.boostGamma * config_.mesh.cellSize[2];
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
    if (config_.checkpoint.restartFrom.empty())
      {
        nextTrajectorySampleTime_ = 0.0;
        trajectorySamplesSinceFlush_ = 0;
      }
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
                << ", energy=" << detectorConfig_.energyPlanes.size()
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
              fields_->fieldEnergyDensityCell(i, j, k)) * volume;
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
                  fields_->poyntingComponentCell(lowerI, j, k, 0)) *
                  area;
                power[1] += static_cast<long double>(
                  fields_->poyntingComponentCell(upperI, j, k, 0)) *
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
                  fields_->poyntingComponentCell(i, lowerJ, k, 1)) *
                  area;
                power[3] += static_cast<long double>(
                  fields_->poyntingComponentCell(i, upperJ, k, 1)) *
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
                  fields_->poyntingComponentCell(i, j, k, 2)) * area;
          }
        if (upperGlobalK >= localZOffset_ &&
            upperGlobalK < localGlobalEnd)
          {
            const std::size_t k = upperGlobalK - localZOffset_;
            for (std::size_t j = jBegin; j < jEnd; ++j)
              for (std::size_t i = iBegin; i < iEnd; ++i)
                power[5] += static_cast<long double>(
                  fields_->poyntingComponentCell(i, j, k, 2)) * area;
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
        &runMetadata_, config_.mesh.ballisticOutsideTransverse,
        config_.mesh.movingWindow.enabled);
    if (config_.checkpoint.restartFrom.empty() && !movingWindow_)
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
        if (config_.mesh.ballisticOutsideTransverse)
          logRoot(communicator_,
            "WARNING: energy-ledger closure_valid is false in approximate "
            "boundary-injection mode. The ledger records both raw and "
            "boundary-source-corrected residuals; neither is an exact "
            "whole-bunch conservation proof.");
        if (config_.reference.frozenLabPlaneInjection)
          logRoot(communicator_,
            "WARNING: energy-ledger closure_valid is false for frozen lab-plane injection; the pre-plane force gate and initial field model are not an exact closed-system energy exchange. Use field-energy and orbit comparisons as convergence diagnostics, not an exact energy-conservation proof.");
        if (!config_.mesh.selfFieldFeedback)
          logRoot(communicator_,
            "WARNING: self_field_feedback is disabled. Maxwell sources and fields remain active, but grid E/B do not act on particles; the energy ledger is diagnostic only and closure_valid is false.");
        if (config_.mesh.ballisticOutsideTransverse &&
            config_.mesh.currentFilter != EBCurrentFilter::None)
          logRoot(communicator_,
            "WARNING: boundary-source work is estimated from unfiltered "
            "particle displacement while Maxwell uses filtered current; "
            "set mesh.current_filter: none for the clearest energy audit.");
      }
  }

  void Simulation::advanceEnergyLedgerStep()
  {
    if (movingWindow_)
      {
        // Called after MPI electric interfaces have completed, before remap.
        movingWindow_->fields->finishEnergyLedgerStep(*fields_);
        const auto& ledger = movingWindow_->fields->energyLedger();
        for (unsigned int face = 0; face < 6; ++face)
          energyLedgerOutwardEnergy_[face] = ledger.outwardEnergyFace[face];
        return;
      }
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

  void Simulation::finalizeBoundaryInjectionWorkStep()
  {
    if (boundaryInjectionWorkSamples_.empty()) return;
    for (std::size_t index = 0;
         index < boundaryInjectionWorkSamples_.size(); ++index)
      {
        const BoundaryInjectionWorkSample& sample =
          boundaryInjectionWorkSamples_[index];
        FieldVector<Double> newElectric(0.0);
        FieldVector<Double> ignoredMagnetic(0.0);
        try
          {
            fields_->sampleFieldsPosition(sample.midpoint, localOriginBox_,
              newElectric, ignoredMagnetic);
          }
        catch (const std::out_of_range& error)
          {
            throw std::out_of_range(std::string(
              "boundary-injection ledger new-field sample: ") +
              error.what());
          }
        const long double newDot =
          static_cast<long double>(sample.chargeDisplacement[0]) *
            newElectric[0] +
          static_cast<long double>(sample.chargeDisplacement[1]) *
            newElectric[1] +
          static_cast<long double>(sample.chargeDisplacement[2]) *
            newElectric[2];
        const long double entryWork = -0.5L *
          (static_cast<long double>(sample.oldElectricDot) + newDot);
        energyLedgerBoundarySourceWork_ += entryWork;
        energyLedgerBoundarySourceAbsWork_ += std::abs(entryWork);
      }
    boundaryInjectionWorkSamples_.clear();
  }

  void Simulation::sampleEnergyLedger(bool force,
                                      bool guardFinalCollectives)
  {
    if (!config_.energyLedger.enabled) return;
    if (energyLedgerHasReference_ && step_ == energyLedgerLastSampleStep_)
      return;
    if (!force &&
        step_ % config_.energyLedger.sampleIntervalSteps != 0)
      return;

    long double local[46] = {};
    Double localMinimumGamma = std::numeric_limits<Double>::infinity();
    Double localMaximumGamma =
      -std::numeric_limits<Double>::infinity();
    PendingParticleStatisticsSnapshot pending;
    long double pendingMeanLabZ = 0., pendingVarianceZ = 0., pendingGammaZ = 0.;
    if (movingWindow_ && movingWindow_->pendingStatistics)
      {
        const Double pendingStart = performanceNow();
        pending = movingWindow_->pendingStatistics->snapshot();
        if (pending.count)
          {
            const long double n = step_ - movingWindow_->pendingStatistics->initialStep();
            const long double gamma = frame_.gamma();
            const long double meanZ = pending.mean[2] + n * pending.mean[3];
            pendingMeanLabZ = static_cast<long double>(frame_.labZOrigin()) +
              gamma * (meanZ - frame_.boxZOrigin()) +
              static_cast<long double>(frame_.gammaBeta()) * SI::c * timeBoxSI_;
            pendingVarianceZ = gamma * gamma * (pending.centered[2][2] +
              2 * n * pending.centered[2][3] + n * n * pending.centered[3][3]);
            pendingGammaZ = gamma * (pending.centered[0][2] + n * pending.centered[0][3]);
            local[0] = pending.kineticEnergy;
            local[1] = pending.totalEnergy;
            local[2] = pending.representedWeight;
            local[3] = pending.representedWeight * pending.mean[0];
            local[4] = pending.representedWeight * pending.mean[1];
            local[5] = pending.representedWeight * pendingMeanLabZ;
            local[32] = pending.kineticEnergyLab;
            local[33] = pending.cMomentumZBox;
            local[34] = pending.cMomentumZLab;
            localMinimumGamma = pending.minimumGammaLab;
            localMaximumGamma = pending.maximumGammaLab;
          }
        performanceSeconds_[PendingStatisticsStage] += performanceNow() - pendingStart;
      }
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
        const ParticleEnergyMoment auditMoment =
          ParticleEnergyMoment::sample(particle, frame_);
        local[32] += auditMoment.kineticLab;
        local[33] += auditMoment.cPzBox;
        local[34] += auditMoment.cPzLab;
        localMinimumGamma = std::min(localMinimumGamma, gammaLab);
        localMaximumGamma = std::max(localMaximumGamma, gammaLab);
      }
    local[6] = localInteriorFieldEnergy();
    local[7] = energyLedgerRemovedKinetic_;
    local[8] = energyLedgerRemovedTotal_;
    local[9] = energyLedgerPrescribedWork_;
    for (unsigned int face = 0; face < 6; ++face)
      local[10 + face] = energyLedgerOutwardEnergy_[face];
    local[16] = energyLedgerBoundarySourceWork_;

    if (guardFinalCollectives)
      armShutdownWatchdog("final-energy-reductions");
    local[17] = energyLedgerBoundarySourceAbsWork_;
    if (movingWindow_)
      {
        const auto& ledger = movingWindow_->fields->energyLedger();
        local[18] = ledger.insertedEnergyJ;
        local[19] = ledger.croppedEnergyJ;
        local[20] = ledger.remapEnergyJ;
        local[21] = ledger.remapInterfaceCorrectionJ;
        local[22] = ledger.inletOverwriteEnergyJ;
        local[23] = ledger.diagnosticSeconds;
        local[24] = movingWindow_->frozenConstraintWork;
        local[25] = movingWindow_->carrierConstraintWork;
        local[26] = movingWindow_->constraintDiagnosticSeconds;
      }
    local[27] = energyLedgerBorisWork_.gridWork;
    local[28] = energyLedgerBorisWork_.prescribedWork;
    local[29] = energyLedgerBorisWork_.magneticRotationDeltaK;
    local[30] = energyLedgerBorisWork_.totalDeltaK;
    local[31] = energyLedgerBorisWork_.roundoffResidual;
    local[35] = energyLedgerRemovedMoment_.kineticLab;
    local[36] = energyLedgerRemovedMoment_.cPzBox;
    local[37] = energyLedgerRemovedMoment_.cPzLab;
    local[38] = energyLedgerReleasedMoment_.kineticBox;
    local[39] = energyLedgerReleasedMoment_.kineticLab;
    local[40] = energyLedgerReleasedMoment_.cPzBox;
    local[41] = energyLedgerReleasedMoment_.cPzLab;
    local[42] = energyLedgerReleasedCharge_;
    local[43] = energyLedgerPartialWorkSubsteps_;
    local[44] = energyLedgerRemovedMoment_.kineticBox;
    local[45] = energyLedgerBorisWork_.absoluteWork;
    long double global[46] = {};
    MPI_Allreduce(local, global, 46, MPI_LONG_DOUBLE, MPI_SUM,
                  communicator_);
    unsigned long long localCounts[4] = {
      static_cast<unsigned long long>(storedParticleCount()),
      energyLedgerRemovedMacroparticles_, retirementEntryCount_,
      ballisticBoundaryEntries_
    };
    unsigned long long globalCounts[4] = {};
    MPI_Allreduce(localCounts, globalCounts, 4,
      MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator_);
    Double globalMinimumGamma = 0.0;
    Double globalMaximumGamma = 0.0;
    MPI_Allreduce(&localMinimumGamma, &globalMinimumGamma, 1,
                  MPI_DOUBLE, MPI_MIN, communicator_);
    MPI_Allreduce(&localMaximumGamma, &globalMaximumGamma, 1,
                  MPI_DOUBLE, MPI_MAX, communicator_);
    if (guardFinalCollectives) disarmShutdownWatchdog();

    const long double weight = global[2];
    const long double meanGammaLab = weight > 0.0L ?
      global[3] / weight : 0.0L;
    const long double meanGammaBox = weight > 0.0L ?
      global[4] / weight : 0.0L;
    const long double meanZLab = weight > 0.0L ?
      global[5] / weight : 0.0L;
    long double localCentred[4] = {};
    if (pending.count)
      {
        const long double dg = pending.mean[0] - meanGammaLab;
        const long double db = pending.mean[1] - meanGammaBox;
        const long double dz = pendingMeanLabZ - meanZLab;
        const long double w = pending.representedWeight;
        localCentred[0] = pending.centered[0][0] + w * dg * dg;
        localCentred[1] = pending.centered[1][1] + w * db * db;
        localCentred[2] = std::max(0.0L, pendingVarianceZ) + w * dz * dz;
        localCentred[3] = pendingGammaZ + w * dg * dz;
      }
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
    if (guardFinalCollectives)
      armShutdownWatchdog("final-energy-centred-reduction");
    MPI_Allreduce(localCentred, globalCentred, 4, MPI_LONG_DOUBLE,
                  MPI_SUM, communicator_);
    if (guardFinalCollectives) disarmShutdownWatchdog();

    if (!energyLedgerHasReference_)
      {
        energyLedgerInitialKinetic_ = global[0];
        energyLedgerInitialField_ = global[6];
        energyLedgerInitialMoment_.kineticBox = global[0];
        energyLedgerInitialMoment_.kineticLab = global[32];
        energyLedgerInitialMoment_.cPzBox = global[33];
        energyLedgerInitialMoment_.cPzLab = global[34];
        energyLedgerHasReference_ = true;
      }
    const long double particleAccounted = global[0] + global[7];
    const long double initialDynamic =
      energyLedgerInitialKinetic_ + energyLedgerInitialField_;
    const long double currentDynamic = particleAccounted + global[6];
    long double outward = 0.0L;
    for (unsigned int face = 0; face < 6; ++face)
      outward += global[10 + face];
    const long double rawResidual = currentDynamic - initialDynamic +
      outward - global[9];
    const long double residual = rawResidual - global[16] - global[20] - global[22]
      - global[24] - global[25];
    const long double deltaParticle =
      particleAccounted - energyLedgerInitialKinetic_;
    const long double deltaField = global[6] - energyLedgerInitialField_;
    const long double exchangeScale = std::abs(deltaParticle) +
      std::abs(deltaField) + std::abs(outward) + std::abs(global[9]) +
      std::abs(global[16]) + std::abs(global[20]) + std::abs(global[22]) +
      std::abs(global[24]) + std::abs(global[25]);
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
    record.boundaryInjectionEntries = globalCounts[3];
    record.closureValid =
      sources_.maxwellIncidentWaveCount() == 0 && globalCounts[2] == 0 &&
      !config_.mesh.ballisticOutsideTransverse &&
      !config_.reference.frozenLabPlaneInjection &&
      config_.mesh.selfFieldFeedback ? 1 : 0;
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
    record.boundarySourceWork = static_cast<Double>(global[16]);
    record.boundarySourceAbsWork = static_cast<Double>(global[17]);
    record.movingInsertedField = static_cast<Double>(global[18]);
    record.movingCroppedField = static_cast<Double>(global[19]);
    record.movingRemapFieldChange = static_cast<Double>(global[20]);
    record.movingRemapInterfaceCorrection = static_cast<Double>(global[21]);
    record.inletOverwriteWork = static_cast<Double>(global[22]);
    record.longitudinalBoundaryOutward = movingWindow_ ? static_cast<Double>(
      global[14] + global[15] - global[20] - global[22]) : 0.0;
    record.movingEnergyDiagnosticRankSeconds = static_cast<Double>(global[23]);
    record.frozenConstraintWork = static_cast<Double>(global[24]);
    record.carrierConstraintWork = static_cast<Double>(global[25]);
    record.constraintEnergyDiagnosticRankSeconds = static_cast<Double>(global[26]);
    const long double deltaParticleAudit = global[0] + global[44] - energyLedgerInitialKinetic_;
    record.particleKineticBoxRemovedEvent = static_cast<Double>(global[44]);
    record.particleRemovalInterpolationDelta = static_cast<Double>(global[44] - global[7]);
    record.particleGridWork = static_cast<Double>(global[27]);
    record.prescribedBorisWork = static_cast<Double>(global[28]);
    record.magneticRotationDeltaK = static_cast<Double>(global[29]);
    record.borisPushDeltaK = static_cast<Double>(global[30]);
    record.borisRoundoffResidual = static_cast<Double>(global[31]);
    record.borisAbsoluteWork = static_cast<Double>(global[45]);
    record.particleKineticLabActive = static_cast<Double>(global[32]);
    record.particleCPzBoxActive = static_cast<Double>(global[33]);
    record.particleCPzLabActive = static_cast<Double>(global[34]);
    record.particleKineticLabRemoved = static_cast<Double>(global[35]);
    record.particleCPzBoxRemoved = static_cast<Double>(global[36]);
    record.particleCPzLabRemoved = static_cast<Double>(global[37]);
    record.releasedKineticBox = static_cast<Double>(global[38]);
    record.releasedKineticLab = static_cast<Double>(global[39]);
    record.releasedCPzBox = static_cast<Double>(global[40]);
    record.releasedCPzLab = static_cast<Double>(global[41]);
    record.releasedCharge = static_cast<Double>(global[42]);
    record.partialWorkSubsteps = static_cast<Double>(global[43]);
    const long double deltaLab = global[32] + global[35] -
      energyLedgerInitialMoment_.kineticLab;
    const long double deltaCPz = global[33] + global[36] -
      energyLedgerInitialMoment_.cPzBox;
    record.deltaParticleKineticBox = static_cast<Double>(deltaParticleAudit);
    record.deltaParticleKineticLab = static_cast<Double>(deltaLab);
    record.deltaParticleCPzBox = static_cast<Double>(deltaCPz);
    record.deltaFieldEnergy = static_cast<Double>(deltaField);
    record.boostGamma = frame_.gamma();
    record.boostBeta = frame_.beta();
    record.particleLorentzDeltaResidual = static_cast<Double>(deltaLab -
      static_cast<long double>(frame_.gamma()) * deltaParticleAudit -
      static_cast<long double>(frame_.gammaBeta()) * deltaCPz);
    const long double fieldExchange = deltaField + outward - global[16] -
      global[20] - global[22] - global[24] - global[25] + global[27];
    record.fieldExchangeResidual = static_cast<Double>(fieldExchange);
    record.particleWorkBalanceResidual = static_cast<Double>(deltaParticleAudit -
      global[27] - global[28] - global[29]);
    record.balanceResidualBoris = static_cast<Double>(deltaParticleAudit + deltaField +
      outward - global[28] - global[16] - global[20] - global[22] -
      global[24] - global[25]);
    record.balanceResidualRaw = static_cast<Double>(rawResidual);
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
          << " longitudinal_boundary_outward_J="
          << lastEnergyLedgerRecord_.longitudinalBoundaryOutward
          << " moving_remap_field_change_J="
          << lastEnergyLedgerRecord_.movingRemapFieldChange
          << " inlet_overwrite_work_J="
          << lastEnergyLedgerRecord_.inletOverwriteWork
          << " frozen_constraint_work_J="
          << lastEnergyLedgerRecord_.frozenConstraintWork
          << " carrier_constraint_work_J="
          << lastEnergyLedgerRecord_.carrierConstraintWork
          << " boundary_source_work_J="
          << lastEnergyLedgerRecord_.boundarySourceWork
          << " boundary_source_abs_work_J="
          << lastEnergyLedgerRecord_.boundarySourceAbsWork
          << " boundary_source_abs_work_over_field_energy=";
        if (lastEnergyLedgerRecord_.fieldEnergyInterior > 0.0)
          message << lastEnergyLedgerRecord_.boundarySourceAbsWork /
            lastEnergyLedgerRecord_.fieldEnergyInterior;
        else
          message << "undefined_zero_field_energy";
        message
          << " balance_residual_raw_J="
          << lastEnergyLedgerRecord_.balanceResidualRaw
          << " residual_boris_J=" << lastEnergyLedgerRecord_.balanceResidualBoris
          << " particle_work_residual_J=" << lastEnergyLedgerRecord_.particleWorkBalanceResidual
          << " field_exchange_residual_J=" << lastEnergyLedgerRecord_.fieldExchangeResidual
          << " delta_K_lab_J=" << lastEnergyLedgerRecord_.deltaParticleKineticLab
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
        const Double warningRatio =
          config_.mesh.ballisticOutsideTransverse ?
            std::abs(lastEnergyLedgerRecord_.relativeBalanceExchange) :
            std::abs(lastEnergyLedgerRecord_.relativeBalanceInitial);
        if ((lastEnergyLedgerRecord_.closureValid ||
             config_.mesh.ballisticOutsideTransverse) &&
            warningRatio >
              config_.energyLedger.warningRelativeTolerance)
          {
            std::ostringstream warning;
            warning << std::setprecision(6)
              << "WARNING: energy-ledger residual is "
              << warningRatio
              << (config_.mesh.ballisticOutsideTransverse ?
                  " of the summed modeled energy exchange" :
                  " of the initial boosted-frame kinetic-plus-field energy")
              << ", "
                 "above warning_relative_tolerance="
              << config_.energyLedger.warningRelativeTolerance
              << ". Refine the field grid/time step, increase macro-particle "
                 "count, enlarge the initial-field padding, and check that "
                 "the residual converges before interpreting small radiation "
                 "or energy-spread changes.";
            if (config_.mesh.ballisticOutsideTransverse)
              warning << " Boundary injection is approximate; compare the "
                         "raw residual and boundary-source work, and scan "
                         "the transverse aperture.";
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
    if (!trajectoryWriter_.isOpen() || particle.frozen) return;
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
    if (config_.mesh.movingWindow.enabled && boxZ < globalOriginBox_[2])
      return 0;
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
    ParticleWorkspace& workspace = *particleWorkspace_;
    std::vector<RelativisticParticleSI>& retained = workspace.retained;
    retained.clear();
    retained.reserve(particles_.size());
    std::vector<ParticleCPMLCarrier>& retainedCarriers = workspace.retainedCarriers;
    retainedCarriers.clear();
    retainedCarriers.reserve(pmlCarriers_.size());
    std::vector<ParticleCPMLCarrier>& retainedRetirementCarriers =
      workspace.retainedRetirementCarriers;
    retainedRetirementCarriers.clear();
    retainedRetirementCarriers.reserve(retirementCarriers_.size());
    std::vector<TransferPacket>& sendLower = workspace.sendLower;
    std::vector<TransferPacket>& sendUpper = workspace.sendUpper;
    std::vector<CarrierTransferPacket>& sendCarrierLower = workspace.sendCarrierLower;
    std::vector<CarrierTransferPacket>& sendCarrierUpper = workspace.sendCarrierUpper;
    std::vector<CarrierTransferPacket>& sendRetirementLower = workspace.sendRetirementLower;
    std::vector<CarrierTransferPacket>& sendRetirementUpper = workspace.sendRetirementUpper;
    sendLower.clear(); sendUpper.clear();
    sendCarrierLower.clear(); sendCarrierUpper.clear();
    sendRetirementLower.clear(); sendRetirementUpper.clear();
    Double* communicationSeconds = performanceInterval_ ?
      &performanceSeconds_[ParticleCommunicationStage] : NULL;
    Double depositionStart = performanceNow();
    const Double lowerZ = localOriginBox_[2];
    const Double upperZ = std::fma(
      static_cast<Double>(localZOffset_ + localGeometry_.nz),
      globalGeometry_.dz, globalOriginBox_[2]);
    const Double crossingTolerance = 64.0 *
      std::numeric_limits<Double>::epsilon();
    const Double carrierCutoff = 64.0 *
      std::numeric_limits<Double>::epsilon();
    unsigned long long localLost = 0;
    int localAutomaticLongitudinalExit = 0;

    const auto interfaceCrossing = [&](const FieldVector<Double>& start,
                                       const FieldVector<Double>& end,
                                       bool& crossLower,
                                       bool& crossUpper)
    {
      crossLower = rank_ > 0 && end[2] < lowerZ;
      crossUpper = rank_ + 1 < size_ && end[2] >= upperZ;
      if (!crossLower && !crossUpper) return 1.0;
      const Double boundary = crossLower ? lowerZ : upperZ;
      Double fraction = (boundary - start[2]) /
        (end[2] - start[2]);
      const Double endpointTolerance = 128.0 *
        std::numeric_limits<Double>::epsilon() *
        std::max(globalGeometry_.dz,
          std::max(std::abs(globalOriginBox_[2]), std::abs(boundary)));
      if (fraction < 0.0 &&
          std::abs(start[2] - boundary) <= endpointTolerance)
        fraction = 0.0;
      if (fraction > 1.0 &&
          std::abs(end[2] - boundary) <= endpointTolerance)
        fraction = 1.0;
      if (!(fraction >= 0.0 && fraction <= 1.0) ||
          !std::isfinite(fraction))
        {
          std::ostringstream message;
          message << std::setprecision(17)
                  << "Invalid MPI-interface particle crossing: rank="
                  << rank_ << " start_z=" << start[2]
                  << " end_z=" << end[2]
                  << " boundary_z=" << boundary
                  << " fraction=" << fraction;
          throw std::runtime_error(message.str());
        }
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
      if (movingWindow_ && config_.energyLedger.enabled)
        queueMovingConstraintWork(segmentStart, segmentEnd,
          carrier.charge * depositedWeight, true);

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
          const std::size_t face =
            static_cast<std::size_t>(selectedHit.face);
          if (config_.mesh.longitudinal.automatic() && face >= 4)
            localAutomaticLongitudinalExit = 1;
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

    // Store raw geometric hits only. Diagnostic termination and deposition
    // keep their distinct event-selection rules and floating-point ordering.
    struct ActiveSegmentBoundaryHits
    {
      ParticleBoundaryHit cpmlHit, outerHit, retirementHit;
      bool entersCpml, exitsOuter, entersRetirement;
    };

    const auto processActiveSegment =
      [&](RelativisticParticleSI particle,
          FieldVector<Double> segmentStart,
          Double segmentStartFraction,
          bool allowTransfer,
          const FieldVector<Double>* activationKink,
          Double activationKinkFraction,
          bool pendingActivePush,
          const ActiveSegmentBoundaryHits* cachedBoundaryHits)
    {
      if (config_.mesh.ballisticOutsideTransverse)
        {
          const unsigned int legCount = activationKink ? 2 : 1;
          FieldVector<Double> legStart(segmentStart);
          Double legStartFraction = segmentStartFraction;
          bool transferred = false;
          for (unsigned int leg = 0; leg < legCount; ++leg)
            {
              const bool frozenLeg = activationKink && leg == 0;
              const FieldVector<Double>& legTarget = frozenLeg ?
                *activationKink : particle.position;
              const Double legTargetFraction = frozenLeg ?
                activationKinkFraction : 1.0;
              bool crossLower = false;
              bool crossUpper = false;
              const Double interfaceFraction = interfaceCrossing(
                legStart, legTarget, crossLower, crossUpper);
              const bool crossesInterface = crossLower || crossUpper;
              if (crossesInterface && !allowTransfer)
                throw std::runtime_error(
                  "An injected particle crossed more than one MPI interface "
                  "in one field step");
              const Double localEndFraction = crossesInterface ?
                interfaceFraction : 1.0;
              FieldVector<Double> legEnd(legTarget);
              if (crossesInterface)
                for (unsigned int axis = 0; axis < 3; ++axis)
                  legEnd[axis] = legStart[axis] +
                    localEndFraction *
                    (legTarget[axis] - legStart[axis]);
              Double first = 0.0;
              Double last = 1.0;
              if (clipPhysicalSegment(legStart, legEnd,
                                      *particleCPML_, first, last))
                {
                  FieldVector<Double> entry(0.0);
                  FieldVector<Double> exit(0.0);
                  for (unsigned int axis = 0; axis < 3; ++axis)
                    {
                      const Double advance =
                        legEnd[axis] - legStart[axis];
                      entry[axis] = legStart[axis] + first * advance;
                      exit[axis] = legStart[axis] + last * advance;
                    }
                  depositor.depositSegment(entry, exit, particle.charge);
                  if (first > 0.0 &&
                      !particleCPML_->containsPhysical(legStart))
                    {
                      ++ballisticBoundaryEntries_;
                      if (config_.energyLedger.enabled)
                        {
                          BoundaryInjectionWorkSample sample;
                          sample.midpoint = entry;
                          for (unsigned int axis = 0; axis < 3; ++axis)
                            {
                              sample.midpoint[axis] =
                                0.5 * (entry[axis] + exit[axis]);
                              sample.chargeDisplacement[axis] =
                                particle.charge * (exit[axis] - entry[axis]);
                            }
                          FieldVector<Double> oldElectric(0.0);
                          FieldVector<Double> ignoredMagnetic(0.0);
                          try
                            {
                              fields_->sampleFieldsPosition(
                                sample.midpoint, localOriginBox_, oldElectric,
                                ignoredMagnetic);
                            }
                          catch (const std::out_of_range& error)
                            {
                              throw std::out_of_range(std::string(
                                "boundary-injection ledger old-field sample: ") +
                                error.what());
                            }
                          sample.oldElectricDot =
                            sample.chargeDisplacement * oldElectric;
                          boundaryInjectionWorkSamples_.push_back(sample);
                        }
                    }
                }
              if (crossesInterface)
                {
                  TransferPacket transfer = {};
                  transfer.particle = packParticle(particle);
                  transfer.pendingActivePush = pendingActivePush ? 1 : 0;
                  if (pendingActivePush)
                    transfer.activationKinkFraction =
                      activationKinkFraction;
                  for (unsigned int axis = 0; axis < 3; ++axis)
                    transfer.segmentStart[axis] = legEnd[axis];
                  transfer.segmentStartFraction = legStartFraction +
                    localEndFraction *
                    (legTargetFraction - legStartFraction);
                  if (frozenLeg)
                    {
                      transfer.hasActivationKink = 1;
                      transfer.activationKinkFraction =
                        activationKinkFraction;
                      for (unsigned int axis = 0; axis < 3; ++axis)
                        transfer.activationKink[axis] =
                          (*activationKink)[axis];
                    }
                  (crossLower ? sendLower : sendUpper).push_back(transfer);
                  transferred = true;
                  break;
                }
              legStart = legEnd;
              legStartFraction = legTargetFraction;
            }
          if (!transferred)
            retained.push_back(particle);
          if (config_.mesh.longitudinal.automatic() &&
              (particle.position[2] < particleCPML_->physicalLower(2) ||
               particle.position[2] > particleCPML_->physicalUpper(2)))
            localAutomaticLongitudinalExit = 1;
          return;
        }
      if (activationKink)
        {
          bool crossLower = false;
          bool crossUpper = false;
          const Double interfaceFraction = interfaceCrossing(
            segmentStart, *activationKink, crossLower, crossUpper);
          ParticleBoundaryHit cpmlHit;
          ParticleBoundaryHit outerHit;
          ParticleBoundaryHit retirementHit;
          const bool entersCpml = particleCPML_->firstEntry(
            segmentStart, *activationKink, cpmlHit);
          const bool exitsOuter = particleBoundary_->firstExit(
            segmentStart, *activationKink, outerHit);
          const bool entersRetirement = labPlaneCrossing(
                segmentStart, *activationKink,
                segmentStartFraction,
                config_.particleRetirement.entranceZ, retirementHit);
          if (entersCpml || exitsOuter || entersRetirement)
            throw std::runtime_error(
              "Held-transverse activation crossed a particle boundary "
              "during the frozen part of a field step; enable "
              "mesh.ballistic_outside_transverse for boundary injection");
          if (crossLower || crossUpper)
            {
              if (!allowTransfer)
                throw std::runtime_error(
                  "A held-transverse frozen segment crossed more than one "
                  "MPI interface in one field step");
              FieldVector<Double> interfacePosition(0.0);
              for (unsigned int axis = 0; axis < 3; ++axis)
                interfacePosition[axis] = segmentStart[axis] +
                  interfaceFraction *
                  ((*activationKink)[axis] - segmentStart[axis]);
              depositor.depositSegment(segmentStart, interfacePosition,
                                       particle.charge);
              if (movingWindow_ && config_.energyLedger.enabled)
                queueMovingConstraintWork(segmentStart, interfacePosition,
                  particle.charge, false);
              TransferPacket transfer = {};
              transfer.particle = packParticle(particle);
              transfer.pendingActivePush = pendingActivePush ? 1 : 0;
              if (pendingActivePush)
                transfer.activationKinkFraction =
                  activationKinkFraction;
              transfer.segmentStartFraction = segmentStartFraction +
                interfaceFraction *
                (activationKinkFraction - segmentStartFraction);
              transfer.hasActivationKink = 1;
              transfer.activationKinkFraction =
                activationKinkFraction;
              for (unsigned int axis = 0; axis < 3; ++axis)
                {
                  transfer.segmentStart[axis] =
                    interfacePosition[axis];
                  transfer.activationKink[axis] =
                    (*activationKink)[axis];
                }
              (crossLower ? sendLower : sendUpper).push_back(transfer);
              return;
            }
          depositor.depositSegment(segmentStart, *activationKink,
                                   particle.charge);
          if (movingWindow_ && config_.energyLedger.enabled)
            queueMovingConstraintWork(segmentStart, *activationKink,
              particle.charge, false);
          segmentStart = *activationKink;
          segmentStartFraction = activationKinkFraction;
        }
      bool crossLower = false;
      bool crossUpper = false;
      const Double interfaceFraction = interfaceCrossing(
        segmentStart, particle.position, crossLower, crossUpper);
      const bool crossesInterface = crossLower || crossUpper;
      ActiveSegmentBoundaryHits computedBoundaryHits;
      if (!cachedBoundaryHits)
        {
          if (config_.mesh.movingWindow.enabled)
            computedBoundaryHits.entersCpml = firstMovingWindowEntry(
              segmentStart, particle.position, segmentStartFraction,
              computedBoundaryHits.cpmlHit);
          else
            computedBoundaryHits.entersCpml = particleCPML_->firstEntry(
              segmentStart, particle.position, computedBoundaryHits.cpmlHit);
          computedBoundaryHits.exitsOuter = particleBoundary_->firstExit(
            segmentStart, particle.position, computedBoundaryHits.outerHit);
          computedBoundaryHits.entersRetirement = labPlaneCrossing(
            segmentStart, particle.position, segmentStartFraction,
            config_.particleRetirement.entranceZ,
            computedBoundaryHits.retirementHit);
        }
      const ActiveSegmentBoundaryHits& boundaryHits = cachedBoundaryHits ?
        *cachedBoundaryHits : computedBoundaryHits;
      const ParticleBoundaryHit& cpmlHit = boundaryHits.cpmlHit;
      const ParticleBoundaryHit& outerHit = boundaryHits.outerHit;
      const ParticleBoundaryHit& retirementHit = boundaryHits.retirementHit;
      const bool entersCpml = boundaryHits.entersCpml;
      const bool exitsOuter = boundaryHits.exitsOuter;
      const bool entersRetirement = boundaryHits.entersRetirement;

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
      if (movingWindow_ && config_.energyLedger.enabled && particle.frozen)
        queueMovingConstraintWork(segmentStart, segmentEnd,
          particle.charge, false);
      const Double globalEndFraction = segmentStartFraction +
        localEndFraction * (1.0 - segmentStartFraction);
      if (endAction == EnterCpml)
        {
          const std::size_t face =
            static_cast<std::size_t>(selectedHit.face);
          const bool movingExit = config_.mesh.movingWindow.enabled &&
            selectedHit.face == ParticleBoundaryFace::UpperZ;
          if (movingExit)
            recordMovingWindowExit(particle.charge);
          else
            {
              ++cpmlEntryCount_[face];
              cpmlEntryCharge_[face] += particle.charge;
            }
          ++localLost;
          if (config_.mesh.longitudinal.automatic() && face >= 4 && !movingExit)
            localAutomaticLongitudinalExit = 1;

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
          transfer.pendingActivePush = pendingActivePush ? 1 : 0;
          if (pendingActivePush)
            transfer.activationKinkFraction = activationKinkFraction;
          for (unsigned int axis = 0; axis < 3; ++axis)
            transfer.segmentStart[axis] = segmentEnd[axis];
          transfer.segmentStartFraction = globalEndFraction;
          (crossLower ? sendLower : sendUpper).push_back(transfer);
        }
      else
        retained.push_back(particle);
    };

    const auto prepareParticle = [&](std::size_t index,
                                     ParticleWorkspace::Prepared& result)
      {
        result.particle = particles_[index];
        RelativisticParticleSI& particle = result.particle;
        const FieldVector<Double> start(particle.position);
        if (config_.mesh.movingWindow.enabled &&
            start[2] < globalOriginBox_[2] && !particle.frozen)
          throw std::runtime_error(
            "An active particle is upstream of the moving window; "
            "the inlet guard did not contain its release event");
        result.prescribedWork = 0.0;
        const bool transportAudit = detectors_ && detectors_->hasEnergyPlanes();
        const bool auditEnabled = config_.energyLedger.enabled || transportAudit;
        if (auditEnabled) {
          result.borisWork.clear();
          result.auditTrace.reset(frame_, transportAudit ? detectors_.get() : NULL,
                                  config_.energyLedger.enabled);
          if(particleWorkGeometry_)result.auditTrace.geometry.reset(fields_.get(),localOriginBox_,
            config_.boundary.cpml.cells[0],config_.boundary.cpml.cells[1],!particle.frozen,particleWorkGeometry_->halo());
        }
        result.activeStart = particle;
        result.activeStartTime = timeBoxSI_;
        result.hasActiveSegment = !particle.frozen;
        result.pendingActivePush = false;
        result.captureFrozenSegment = false;
        Double& prescribedWork = result.prescribedWork;
        RelativisticParticleSI& activeStart = result.activeStart;
        Double& activeStartTime = result.activeStartTime;
        bool& hasActiveSegment = result.hasActiveSegment;
        bool& pendingActivePush = result.pendingActivePush;
        if (particle.frozen)
          {
            const bool holdTransverse =
              config_.reference.holdTransverseUntilActivation;
            const Double betaHeld = holdTransverse ?
              particle.frozenBetaZ : 0.0;
            const Double startLabZ = frame_.labZFromBoxZT(
              start[2], timeBoxSI_);
            if (holdTransverse)
              particle.position[2] += SI::c * globalGeometry_.dt *
                betaHeld;
            else
              RelativisticBorisPusher::pushPosition(
                particle.position, particle.properVelocity,
                globalGeometry_.dt);
            const Double endLabZ = frame_.labZFromBoxZT(
              particle.position[2], timeBoxSI_ + globalGeometry_.dt);
            if (endLabZ >= config_.reference.activationPlaneZ)
              {
                if (config_.mesh.movingWindow.enabled &&
                    start[2] < globalOriginBox_[2])
                  throw std::runtime_error(
                    "A pending upstream held particle would release before "
                    "a complete in-window step; check the moving inlet guard");
                const Double crossingFraction =
                  startLabZ >= config_.reference.activationPlaneZ ? 0.0 :
                  std::max(0.0, std::min(1.0,
                    (config_.reference.activationPlaneZ - startLabZ) /
                    (endLabZ - startLabZ)));
                particle.position = start;
                if (holdTransverse)
                  particle.position[2] += SI::c *
                    crossingFraction * globalGeometry_.dt * betaHeld;
                else
                  RelativisticBorisPusher::pushPosition(
                    particle.position, particle.properVelocity,
                    crossingFraction * globalGeometry_.dt);
                particle.frozen = false;
                hasActiveSegment = true;
                activeStart = particle;
                activeStartTime = timeBoxSI_ +
                  crossingFraction * globalGeometry_.dt;
                result.captureFrozenSegment = crossingFraction > 0.0;
                const Double remaining =
                  (1.0 - crossingFraction) * globalGeometry_.dt;
                if (remaining > 0.0)
                  {
                    // The activation point can cross a z-slab interface
                    // during the frozen part of this field step.  The
                    // receiving rank must sample its own E/B before the
                    // active remainder is pushed; the sender's interpolation
                    // halo does not extend to every possible release point.
                    const int activationOwner = ownerRank(
                      activeStart.position[2]);
                    if (activationOwner >= 0 && activationOwner != rank_)
                      pendingActivePush = true;
                    else if (config_.mesh.ballisticOutsideTransverse &&
                        !particleCPML_->containsPhysical(
                          activeStart.position))
                      RelativisticBorisPusher::pushPosition(
                        particle.position, particle.properVelocity,
                        remaining);
                    else
                      RelativisticBorisPusher::
                        pushFromGridAndPrescribedLabSubcycled(
                          particle, *fields_, localOriginBox_, sources_,
                          frame_, activeStartTime, remaining,
                          particleSubsteps_,
                          config_.energyLedger.enabled ?
                            &prescribedWork : 0,
                          config_.mesh.selfFieldFeedback,
                          config_.energyLedger.enabled ? &result.borisWork : NULL,
                          auditEnabled ? &result.auditTrace : NULL);
                  }
              }
            else if (holdTransverse)
              result.captureFrozenSegment = true;
          }
        else if (config_.mesh.ballisticOutsideTransverse &&
            !particleCPML_->containsPhysical(start))
          RelativisticBorisPusher::pushPosition(
            particle.position, particle.properVelocity,
            globalGeometry_.dt);
        else
          RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
            particle, *fields_, localOriginBox_, sources_, frame_,
            timeBoxSI_, globalGeometry_.dt, particleSubsteps_,
            config_.energyLedger.enabled ? &prescribedWork : 0,
            config_.mesh.selfFieldFeedback,
            config_.energyLedger.enabled ? &result.borisWork : NULL,
            auditEnabled ? &result.auditTrace : NULL);
        if(particleWorkGeometry_)result.auditTrace.geometry.finish();
      };

    performanceSeconds_[ParticleDepositStage] += performanceNow() - depositionStart;
    // Only preparation reads E/B and modifies independent particle copies.
    // Deposit, boundary accounting, detector events and migration packing keep
    // their original order, including release events crossing MPI interfaces.
    const std::size_t batchCapacity = std::min<std::size_t>(
      ParticleWorkspace::PreparationBatchCapacity, particles_.size());
    workspace.prepared.resize(batchCapacity);
    for (std::size_t batch = 0; batch < particles_.size(); batch += batchCapacity)
      {
        const std::ptrdiff_t count = static_cast<std::ptrdiff_t>(
          std::min(batchCapacity, particles_.size() - batch));
        const Double preparationStart = performanceNow();
        if (computeThreadCount() == 1 || count < 32)
          {
            for (std::ptrdiff_t offset = 0; offset < count; ++offset)
              prepareParticle(batch + static_cast<std::size_t>(offset),
                workspace.prepared[static_cast<std::size_t>(offset)]);
          }
        else
          {
            std::exception_ptr preparationError;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for schedule(static) num_threads(computeThreadCount())
#endif
            for (std::ptrdiff_t offset = 0; offset < count; ++offset)
              {
                try
                  {
                    prepareParticle(batch + static_cast<std::size_t>(offset),
                      workspace.prepared[static_cast<std::size_t>(offset)]);
                  }
                catch (...)
                  {
#ifdef APRL_HAS_OPENMP
#pragma omp critical(aprl_particle_prepare_error)
#endif
                    {
                      if (!preparationError)
                        preparationError = std::current_exception();
                    }
                  }
              }
            if (preparationError) std::rethrow_exception(preparationError);
          }
        performanceSeconds_[ParticlePrepareStage] += performanceNow() - preparationStart;
        depositionStart = performanceNow();
        for (std::ptrdiff_t offset = 0; offset < count; ++offset)
          {
            const std::size_t index = batch + static_cast<std::size_t>(offset);
            ParticleWorkspace::Prepared& result =
              workspace.prepared[static_cast<std::size_t>(offset)];
            if(particleWorkGeometry_)particleWorkGeometry_->add(result.auditTrace.geometry);
            RelativisticParticleSI& particle = result.particle;
            FieldVector<Double> start(particles_[index].position);
            Double segmentStartFraction = 0.0;
            const Double prescribedWork = result.prescribedWork;
            const RelativisticParticleSI& activeStart = result.activeStart;
            const Double activeStartTime = result.activeStartTime;
            const bool hasActiveSegment = result.hasActiveSegment;
            const bool pendingActivePush = result.pendingActivePush;
            if (particles_[index].frozen && hasActiveSegment)
              {
                ++frozenPlaneActivations_;
                if (config_.energyLedger.enabled) {
                  energyLedgerReleasedMoment_.add(
                    ParticleEnergyMoment::sample(activeStart, frame_));
                  energyLedgerReleasedCharge_ += activeStart.charge;
                }
              }
            if (detectors_ && detectors_->hasEnergyPlanes() &&
                particles_[index].frozen && (!hasActiveSegment || activeStartTime > timeBoxSI_))
              detectors_->captureEnergyParticleStep(particles_[index],
                hasActiveSegment ? activeStart : particle, timeBoxSI_,
                hasActiveSegment ? activeStartTime : timeBoxSI_ + globalGeometry_.dt);
            if (detectors_ && result.captureFrozenSegment)
              detectors_->captureParticleStep(particles_[index],
                hasActiveSegment ? activeStart : particle, timeBoxSI_,
                hasActiveSegment ? activeStartTime : timeBoxSI_ + globalGeometry_.dt);

            // The full upstream source list remains on rank zero. Held
            // preparation advances it without field sampling; only its
            // in-window part contributes to this step's mesh current.
            if (config_.mesh.movingWindow.enabled &&
                start[2] < globalOriginBox_[2])
              {
                if (rank_ != 0 || !particles_[index].frozen ||
                    !particle.frozen || hasActiveSegment || pendingActivePush)
                  throw std::runtime_error(
                    "Moving-window upstream admission requires an unreleased held particle on rank zero");
                if (particle.position[2] < globalOriginBox_[2])
                  {
                    retained.push_back(particle);
                    continue;
                  }
                const Double fraction = (globalOriginBox_[2] - start[2]) /
                  (particle.position[2] - start[2]);
                if (!(fraction >= 0.0 && fraction <= 1.0) ||
                    !std::isfinite(fraction))
                  throw std::runtime_error("Invalid moving-window admission fraction");
                for (unsigned axis = 0; axis < 3; ++axis)
                  start[axis] += fraction * (particle.position[axis] - start[axis]);
                start[2] = globalOriginBox_[2];
                segmentStartFraction = fraction;
              }

            const bool hasActivationKink =
              config_.reference.holdTransverseUntilActivation &&
              particles_[index].frozen && hasActiveSegment &&
              activeStartTime > timeBoxSI_ &&
              activeStartTime < timeBoxSI_ + globalGeometry_.dt;
            const Double activationKinkFraction =
              (activeStartTime - timeBoxSI_) / globalGeometry_.dt;
            const FieldVector<Double>* activationKink =
              hasActivationKink ? &activeStart.position : 0;

            if (pendingActivePush)
              {
                if (activationKinkFraction == 0.0)
                  {
                    // A release event at the start of the step has no frozen
                    // displacement to deposit.  Hand it directly to the owner
                    // instead of attempting a zero-length interface crossing.
                    const int destination = ownerRank(particle.position[2]);
                    if (destination != rank_ - 1 && destination != rank_ + 1)
                      throw std::runtime_error(
                        "Cross-rank zero-length activation is not adjacent");
                    TransferPacket transfer = {};
                    transfer.particle = packParticle(particle);
                    for (unsigned int axis = 0; axis < 3; ++axis)
                      transfer.segmentStart[axis] = start[axis];
                    transfer.pendingActivePush = 1;
                    (destination < rank_ ? sendLower : sendUpper).push_back(
                      transfer);
                  }
                else
                  processActiveSegment(particle, start, segmentStartFraction, true,
                    activationKink, activationKinkFraction, true, 0);
                continue;
              }

            if (config_.mesh.ballisticOutsideTransverse)
              {
                if (config_.energyLedger.enabled) {
                  energyLedgerPrescribedWork_ += prescribedWork;
                  energyLedgerBorisWork_.add(result.borisWork);
                }
                if (detectors_ && detectors_->hasEnergyPlanes() && hasActiveSegment) {
                  if (result.auditTrace.observedSubsteps)
                    detectors_->commitEnergyParticleCrossings(result.auditTrace.crossings,
                      std::numeric_limits<Double>::infinity());
                  else
                    detectors_->captureEnergyParticleStep(activeStart,particle,
                      activeStartTime,timeBoxSI_ + globalGeometry_.dt);
                }
                if (detectors_ && hasActiveSegment &&
                    activeStartTime < timeBoxSI_ + globalGeometry_.dt)
                  detectors_->captureParticleStep(
                    activeStart, particle, activeStartTime,
                    timeBoxSI_ + globalGeometry_.dt);
                processActiveSegment(particle, start, segmentStartFraction, true,
                  activationKink, activationKinkFraction, false, 0);
                continue;
              }

            const FieldVector<Double>& diagnosticStart =
              hasActivationKink ? activeStart.position : start;
            const Double diagnosticStartFraction =
              hasActivationKink ? activationKinkFraction : segmentStartFraction;
            ParticleBoundaryHit diagnosticHit;
            ActiveSegmentBoundaryHits diagnosticBoundaryHits;
            if (config_.mesh.movingWindow.enabled)
              diagnosticBoundaryHits.entersCpml = firstMovingWindowEntry(
                diagnosticStart, particle.position, diagnosticStartFraction,
                diagnosticBoundaryHits.cpmlHit);
            else
              diagnosticBoundaryHits.entersCpml = particleCPML_->firstEntry(
                diagnosticStart, particle.position, diagnosticBoundaryHits.cpmlHit);
            diagnosticBoundaryHits.exitsOuter = particleBoundary_->firstExit(
              diagnosticStart, particle.position, diagnosticBoundaryHits.outerHit);
            diagnosticBoundaryHits.entersRetirement = labPlaneCrossing(
              diagnosticStart, particle.position,
              diagnosticStartFraction,
              config_.particleRetirement.entranceZ,
              diagnosticBoundaryHits.retirementHit);
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
            selectTerminal(diagnosticBoundaryHits.exitsOuter,
                           diagnosticBoundaryHits.outerHit,
                           TrajectoryEvent::DomainExit);
            selectTerminal(diagnosticBoundaryHits.entersCpml,
                           diagnosticBoundaryHits.cpmlHit,
                           config_.mesh.movingWindow.enabled &&
                             diagnosticBoundaryHits.cpmlHit.face == ParticleBoundaryFace::UpperZ ?
                             TrajectoryEvent::DomainExit : TrajectoryEvent::CpmlEntry);
            selectTerminal(diagnosticBoundaryHits.entersRetirement,
                           diagnosticBoundaryHits.retirementHit,
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
                  (diagnosticStartFraction + diagnosticHit.fraction *
                    (1.0 - diagnosticStartFraction)) *
                  globalGeometry_.dt;
                if (config_.energyLedger.enabled)
                  {
                    const ParticleEnergyMoment removed = result.auditTrace.at(
                      diagnosticTime, particle, frame_);
                    energyLedgerRemovedMoment_.add(removed);
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
            if (config_.energyLedger.enabled)
              energyLedgerBorisWork_.add(hasTerminalEvent ? result.auditTrace.prefix(
                diagnosticTime, &energyLedgerPartialWorkSubsteps_) : result.borisWork);
            if (detectors_ && detectors_->hasEnergyPlanes() && hasActiveSegment) {
              if (result.auditTrace.observedSubsteps)
                detectors_->commitEnergyParticleCrossings(result.auditTrace.crossings,
                  hasTerminalEvent ? diagnosticTime : std::numeric_limits<Double>::infinity());
              else
                detectors_->captureEnergyParticleStep(activeStart,diagnosticEnd,
                  activeStartTime,diagnosticTime);
            }
            if (detectors_ && hasActiveSegment &&
                activeStartTime < diagnosticTime)
              detectors_->captureParticleStep(
                activeStart, diagnosticEnd,
                activeStartTime, diagnosticTime);
            if (hasTerminalEvent && config_.trajectory.enabled)
              {
                appendTrajectoryEvent(
                  diagnosticEnd, diagnosticTime,
                  terminalEvent, diagnosticHit.face);
              }
            // Without an activation kink, the diagnostic start and time
            // match this local segment, including clipped inlet admission.
            // Kinked, ballistic and received segments retain their own queries.
            processActiveSegment(particle, start, segmentStartFraction, true,
              activationKink, activationKinkFraction, false,
              hasActivationKink ? 0 : &diagnosticBoundaryHits);
          }

        performanceSeconds_[ParticleDepositStage] += performanceNow() - depositionStart;
      }

    const Double migrationStart = performanceNow();
    std::vector<TransferPacket>& receiveLower = workspace.receiveLower;
    std::vector<TransferPacket>& receiveUpper = workspace.receiveUpper;

    exchangePackets(sendUpper, halo_->upperRank(), halo_->lowerRank(),
      TRANSFER_UPPER_COUNT, TRANSFER_UPPER_DATA,
      communicator_, receiveLower, communicationSeconds);
    exchangePackets(sendLower, halo_->lowerRank(), halo_->upperRank(),
      TRANSFER_LOWER_COUNT, TRANSFER_LOWER_DATA,
      communicator_, receiveUpper, communicationSeconds);

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
          FieldVector<Double> activationKink(0.0);
          if (transfer.hasActivationKink)
            for (unsigned int axis = 0; axis < 3; ++axis)
              activationKink[axis] = transfer.activationKink[axis];
          if (!transfer.pendingActivePush)
            {
              processActiveSegment(
                particle, segmentStart, transfer.segmentStartFraction, false,
                transfer.hasActivationKink ? &activationKink : 0,
                transfer.activationKinkFraction, false, 0);
              continue;
            }

          if (ownerRank(particle.position[2]) != rank_)
            throw std::runtime_error(
              "Cross-rank activation packet has no local release point");
          const std::size_t retainedBefore = retained.size();
          // Deposit the remainder of the frozen trajectory on this slab.
          // Its endpoint is the activation event, not a pushed final state.
          if (transfer.activationKinkFraction > 0.0)
            {
              processActiveSegment(
                particle, segmentStart, transfer.segmentStartFraction, false,
                transfer.hasActivationKink ? &activationKink : 0,
                transfer.activationKinkFraction, false, 0);
              if (retained.size() != retainedBefore + 1)
                throw std::runtime_error(
                  "Cross-rank frozen segment did not reach its activation point");
            }
          RelativisticParticleSI activeStart(
            transfer.activationKinkFraction > 0.0 ?
              retained.back() : particle);
          if (transfer.activationKinkFraction > 0.0)
            retained.pop_back();
          ++crossRankActivationPushes_;
          RelativisticParticleSI activeEnd(activeStart);
          const Double activeStartTime = timeBoxSI_ +
            transfer.activationKinkFraction * globalGeometry_.dt;
          const Double remaining =
            (1.0 - transfer.activationKinkFraction) * globalGeometry_.dt;
          Double prescribedWork = 0.0;
          const bool transportAudit = detectors_ && detectors_->hasEnergyPlanes();
          const bool auditEnabled = config_.energyLedger.enabled || transportAudit;
          BorisWork borisWork;
          ParticleEnergyTrace auditTrace;
          if (auditEnabled) auditTrace.reset(frame_,
            transportAudit ? detectors_.get() : NULL, config_.energyLedger.enabled);
          if (remaining > 0.0)
            {
              if (config_.mesh.ballisticOutsideTransverse &&
                  !particleCPML_->containsPhysical(activeStart.position))
                RelativisticBorisPusher::pushPosition(
                  activeEnd.position, activeEnd.properVelocity, remaining);
              else
                RelativisticBorisPusher::
                  pushFromGridAndPrescribedLabSubcycled(
                    activeEnd, *fields_, localOriginBox_, sources_, frame_,
                    activeStartTime, remaining, particleSubsteps_,
                    config_.energyLedger.enabled ? &prescribedWork : 0,
                    config_.mesh.selfFieldFeedback,
                    config_.energyLedger.enabled ? &borisWork : NULL,
                    auditEnabled ? &auditTrace : NULL);
            }
          Double auditStop = timeBoxSI_ + globalGeometry_.dt;
          bool auditTerminal = false;
          if (auditEnabled && !config_.mesh.ballisticOutsideTransverse) {
            ParticleBoundaryHit cpml, outer, retirement;
            const bool entersCpml = config_.mesh.movingWindow.enabled ?
              firstMovingWindowEntry(activeStart.position,activeEnd.position,
                transfer.activationKinkFraction,cpml) :
              particleCPML_->firstEntry(activeStart.position,activeEnd.position,cpml);
            const bool exitsOuter = particleBoundary_->firstExit(
              activeStart.position,activeEnd.position,outer);
            const bool entersRetirement = labPlaneCrossing(activeStart.position,
              activeEnd.position,transfer.activationKinkFraction,
              config_.particleRetirement.entranceZ,retirement);
            Double fraction = 1.;
            const auto selectAudit = [&](bool hit, const ParticleBoundaryHit& crossing) {
              if (hit && crossing.fraction <= fraction + crossingTolerance) {
                auditTerminal=true; fraction=crossing.fraction;
              }
            };
            selectAudit(exitsOuter,outer); selectAudit(entersCpml,cpml);
            selectAudit(entersRetirement,retirement);
            if (auditTerminal) auditStop=activeStartTime+fraction*remaining;
          }
          if (config_.energyLedger.enabled) {
            energyLedgerPrescribedWork_ += prescribedWork;
            energyLedgerBorisWork_.add(auditTerminal ? auditTrace.prefix(
              auditStop,&energyLedgerPartialWorkSubsteps_) : borisWork);
            if (auditTerminal) {
              const auto removed=auditTrace.at(auditStop,activeEnd,frame_);
              energyLedgerRemovedMoment_.add(removed);
              energyLedgerRemovedKinetic_+=removed.kineticBox;
              energyLedgerRemovedTotal_+=removed.totalBox;
              ++energyLedgerRemovedMacroparticles_;
            }
          }
          if (transportAudit && remaining > 0.) {
            if (auditTrace.observedSubsteps)
              detectors_->commitEnergyParticleCrossings(auditTrace.crossings,
                auditTerminal ? auditStop : std::numeric_limits<Double>::infinity());
            else {
              if (auditTerminal) throw std::logic_error("Terminal audit requires substep history");
              detectors_->captureEnergyParticleStep(activeStart,activeEnd,
                activeStartTime,auditStop);
            }
          }
          if (detectors_ && remaining > 0.0)
            detectors_->captureParticleStep(
              activeStart, activeEnd, activeStartTime,
              timeBoxSI_ + globalGeometry_.dt);
          processActiveSegment(
            activeEnd, activeStart.position,
            transfer.activationKinkFraction, false, 0, 0.0, false, 0);
        }

    std::vector<CarrierTransferPacket>& receiveCarrierLower = workspace.receiveCarrierLower;
    std::vector<CarrierTransferPacket>& receiveCarrierUpper = workspace.receiveCarrierUpper;
    exchangePackets(sendCarrierUpper,
      halo_->upperRank(), halo_->lowerRank(),
      CARRIER_UPPER_COUNT, CARRIER_UPPER_DATA,
      communicator_, receiveCarrierLower, communicationSeconds);
    exchangePackets(sendCarrierLower,
      halo_->lowerRank(), halo_->upperRank(),
      CARRIER_LOWER_COUNT, CARRIER_LOWER_DATA,
      communicator_, receiveCarrierUpper, communicationSeconds);
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

    std::vector<CarrierTransferPacket>& receiveRetirementLower = workspace.receiveRetirementLower;
    std::vector<CarrierTransferPacket>& receiveRetirementUpper = workspace.receiveRetirementUpper;
    receiveRetirementLower.clear(); receiveRetirementUpper.clear();
    if (config_.particleRetirement.enabled)
      {
        exchangePackets(sendRetirementUpper,
          halo_->upperRank(), halo_->lowerRank(),
          RETIREMENT_UPPER_COUNT, RETIREMENT_UPPER_DATA,
          communicator_, receiveRetirementLower, communicationSeconds);
        exchangePackets(sendRetirementLower,
          halo_->lowerRank(), halo_->upperRank(),
          RETIREMENT_LOWER_COUNT, RETIREMENT_LOWER_DATA,
          communicator_, receiveRetirementUpper, communicationSeconds);
      }
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
    localLostParticles_ += localLost;
    int globalAutomaticLongitudinalExit = 0;
    if (config_.mesh.longitudinal.automatic())
      {
        const Double communicationStart = performanceNow();
        MPI_Allreduce(&localAutomaticLongitudinalExit,
                      &globalAutomaticLongitudinalExit, 1, MPI_INT, MPI_MAX,
                      communicator_);
        performanceSeconds_[ParticleCommunicationStage] +=
          performanceNow() - communicationStart;
      }
    if (globalAutomaticLongitudinalExit)
      throw std::runtime_error(
        "An active particle or retirement-current carrier entered a longitudinal CPML in automatic observation-bounded mode. The ballistic envelope underestimated the evolved orbit; increase mesh.longitudinal.ballistic_time_safety_factor or the longitudinal guard cells. This run is not physically complete.");
    performanceSeconds_[ParticleMigrationStage] += performanceNow() - migrationStart;
    const Double diagnosticStart = performanceNow();
    if (detectors_) detectors_->collectParticleCrossings();
    performanceSeconds_[DiagnosticsStage] += performanceNow() - diagnosticStart;
  }

  void Simulation::reportParticleBoundaryLosses() const
  {
    if (config_.reference.frozenLabPlaneInjection)
      {
        const unsigned long long localFrozen = static_cast<unsigned long long>(
          std::count_if(particles_.begin(), particles_.end(),
            [](const RelativisticParticleSI& p) { return p.frozen; })) + pendingParticleCount();
        unsigned long long globalFrozen = 0;
        unsigned long long globalActivations = 0;
        unsigned long long globalCrossRankPushes = 0;
        MPI_Reduce(&localFrozen, &globalFrozen, 1,
          MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator_);
        MPI_Reduce(&frozenPlaneActivations_, &globalActivations, 1,
          MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator_);
        MPI_Reduce(&crossRankActivationPushes_, &globalCrossRankPushes, 1,
          MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator_);
        if (rank_ == 0)
          {
            std::ostringstream message;
            message << "Frozen lab-plane injection: activations="
                    << globalActivations << ", cross-rank activation pushes="
                    << globalCrossRankPushes << ", still frozen="
                    << globalFrozen << ".";
            logRoot(communicator_, message.str());
          }
      }
    if (config_.mesh.ballisticOutsideTransverse)
      {
        unsigned long long localOutside = 0;
        for (std::size_t index = 0; index < particles_.size(); ++index)
          if (!particleCPML_->containsPhysical(particles_[index].position))
            ++localOutside;
        unsigned long long globalOutside = 0;
        unsigned long long globalEntries = 0;
        MPI_Reduce(&localOutside, &globalOutside, 1,
          MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator_);
        MPI_Reduce(&ballisticBoundaryEntries_, &globalEntries, 1,
          MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator_);
        if (rank_ == 0)
          {
            std::ostringstream message;
            message << "Experimental transverse boundary injection: "
                    << "physical-aperture entries=" << globalEntries
                    << ", currently outside=" << globalOutside
                    << ". Entries count crossings, not unique particles; "
                       "external particles are retained without grid or "
                       "prescribed fields.";
            logRoot(communicator_, message.str());
          }
      }
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

  Double Simulation::performanceNow() const
  {
    return performanceInterval_ ? MPI_Wtime() : 0.0;
  }

  void Simulation::reportPerformance(bool final)
  {
    if (!performanceInterval_ || step_ == performanceLastStep_ ||
        (!final && step_ % performanceInterval_ != 0)) return;
    const Double localWall = MPI_Wtime() - performanceWallStart_;
    Double wall = 0.0;
    Double sums[PerformanceStageCount] = {};
    struct ValueRank { Double value; int rank; };
    ValueRank localMax[PerformanceStageCount], maxima[PerformanceStageCount];
    for (int stage = 0; stage < PerformanceStageCount; ++stage)
      {
        localMax[stage].value = performanceSeconds_[stage];
        localMax[stage].rank = rank_;
      }
    MPI_Reduce(&localWall, &wall, 1, MPI_DOUBLE, MPI_MAX, 0, communicator_);
    MPI_Reduce(performanceSeconds_, sums, PerformanceStageCount,
               MPI_DOUBLE, MPI_SUM, 0, communicator_);
    MPI_Reduce(localMax, maxima, PerformanceStageCount,
               MPI_DOUBLE_INT, MPI_MAXLOC, 0, communicator_);
    unsigned long long localCounts[4] = {
      static_cast<unsigned long long>(storedParticleCount()), 0,
      static_cast<unsigned long long>(pendingParticleCount()),
      static_cast<unsigned long long>(pendingParticleCount())
    };
    for (std::size_t index = 0; index < particles_.size(); ++index)
      ++localCounts[particles_[index].frozen ? 2 : 1];
    unsigned long long minimumCounts[4] = {}, maximumCounts[4] = {}, sumCounts[4] = {};
    MPI_Reduce(localCounts, minimumCounts, 4, MPI_UNSIGNED_LONG_LONG,
               MPI_MIN, 0, communicator_);
    MPI_Reduce(localCounts, maximumCounts, 4, MPI_UNSIGNED_LONG_LONG,
               MPI_MAX, 0, communicator_);
    MPI_Reduce(localCounts, sumCounts, 4, MPI_UNSIGNED_LONG_LONG,
               MPI_SUM, 0, communicator_);
    if (rank_ == 0)
      {
        static const char* names[PerformanceStageCount] = {
          "field_total", "particle_prepare", "particle_deposit", "particle_migration_total",
          "diagnostics", "control", "field_communication", "particle_communication",
          "pending_admission", "moving_window_total", "pending_statistics"
        };
        for (int stage = 0; stage < PerformanceStageCount; ++stage)
          {
            std::ostringstream message;
            message << std::setprecision(10)
              << "[performance] phase=" << (final ? "final" : "interval")
              << " step=" << step_ << " steps=" << step_ - performanceLastStep_
              << " elapsed_s=" << wall << " stage=" << names[stage]
              << " mean_s=" << sums[stage] / size_
              << " max_s=" << maxima[stage].value
              << " max_rank=" << maxima[stage].rank;
            logRoot(communicator_, message.str());
          }
        std::ostringstream message;
        message << "[performance] phase=workload step=" << step_
          << " ranks=" << size_ << " threads=" << computeThreadCount()
          << " particles_total=" << sumCounts[0]
          << " particles_min=" << minimumCounts[0] << " particles_max=" << maximumCounts[0]
          << " active_total=" << sumCounts[1]
          << " active_min=" << minimumCounts[1] << " active_max=" << maximumCounts[1]
          << " frozen_total=" << sumCounts[2]
          << " frozen_min=" << minimumCounts[2] << " frozen_max=" << maximumCounts[2]
          << " scheduled_pending_total=" << sumCounts[3]
          << " scheduled_pending_max=" << maximumCounts[3];
        logRoot(communicator_, message.str());
      }
    // Exclude reporting and root's log flush from the following interval on
    // every rank. Profiling itself still adds overhead to whole-run wall time.
    MPI_Barrier(communicator_);
    std::fill(performanceSeconds_, performanceSeconds_ + PerformanceStageCount, 0.0);
    performanceLastStep_ = step_;
    performanceWallStart_ = MPI_Wtime();
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
    if (config_.stop.mode == StopMode::ReferenceCenterZ)
      {
        const Double referenceZ = frame_.labZFromBoxZT(
          referenceCenterBoxZ_, timeBoxSI_);
        const Double tolerance = 64.0 *
          std::numeric_limits<Double>::epsilon() *
          std::max(1.0, std::max(std::abs(referenceZ),
                                 std::abs(config_.stop.referenceZ)));
        if (referenceZ < config_.stop.referenceZ - tolerance) return false;
        MPI_Allreduce(&localLostParticles_, &lostParticles_, 1,
          MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator_);
        std::ostringstream message;
        message << "boost reference centre reached lab z=" << referenceZ
                << " m (target " << config_.stop.referenceZ << " m)";
        if (lostParticles_ > 0)
          message << "; excluded escaped particles=" << lostParticles_;
        stopReason_ = message.str();
        return true;
      }

    const Double interactionExit = lastBeamlineInteractionExitLab();
    // Waiting particles are upstream of the inlet and therefore upstream of
    // the requested last interaction. Do not scan their immutable queue.
    int localNotPast = pendingParticleCount() ? 1 : 0;
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
          {
            localNotPast = 1;
            break;
          }
      }
    int globalNotPast = 0;
    MPI_Allreduce(&localNotPast, &globalNotPast, 1, MPI_INT, MPI_MAX,
                  communicator_);
    if (globalNotPast) return false;

    const unsigned long long localCounts[2] = {
      static_cast<unsigned long long>(storedParticleCount()), localLostParticles_
    };
    unsigned long long globalCounts[2] = {};
    MPI_Allreduce(localCounts, globalCounts, 2, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, communicator_);
    const unsigned long long globalCount = globalCounts[0];
    lostParticles_ = globalCounts[1];

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
    unsigned long long localOutsideTransverse = 0;
    for (std::size_t index = 0; index < particles_.size(); ++index)
      {
        for (unsigned int axis = 0; axis < 3; ++axis)
          if (particles_[index].position[axis] < globalOriginBox_[axis] ||
              particles_[index].position[axis] > upper[axis])
            if (!config_.mesh.ballisticOutsideTransverse || axis == 2)
              localInvalid = 1;
        if (config_.mesh.ballisticOutsideTransverse &&
            (particles_[index].position[0] <
               particleCPML_->physicalLower(0) ||
             particles_[index].position[0] >
               particleCPML_->physicalUpper(0) ||
             particles_[index].position[1] <
               particleCPML_->physicalLower(1) ||
             particles_[index].position[1] >
               particleCPML_->physicalUpper(1)))
          ++localOutsideTransverse;
      }
    for (std::size_t index = 0; index < particles_.size(); ++index)
      if (!particleCPML_->containsPhysical(particles_[index].position) &&
          (!config_.mesh.ballisticOutsideTransverse ||
           particles_[index].position[2] <
             particleCPML_->physicalLower(2) ||
           particles_[index].position[2] >
             particleCPML_->physicalUpper(2)))
        localInCpml = 1;
    int globalInvalid = 0;
    int globalInCpml = 0;
    MPI_Allreduce(&localInvalid, &globalInvalid, 1, MPI_INT, MPI_MAX,
                  communicator_);
    MPI_Allreduce(&localInCpml, &globalInCpml, 1, MPI_INT, MPI_MAX,
                  communicator_);
    if (config_.mesh.ballisticOutsideTransverse)
      {
        unsigned long long globalOutsideTransverse = 0;
        MPI_Allreduce(&localOutsideTransverse, &globalOutsideTransverse, 1,
          MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator_);
        if (rank_ == 0)
          {
            std::ostringstream message;
            message << "Experimental transverse boundary injection: "
                    << globalOutsideTransverse
                    << " boosted particles start outside the physical "
                       "x/y aperture. They remain ballistic without grid "
                       "current until entry; ";
            if (config_.initialSelfField.includeTransverseCpmlParticles)
              message << "only those inside the transverse CPML may "
                         "contribute to the one-time initial Poisson field. ";
            else
              message << "all are excluded from the initial Poisson source. ";
            message << "Gauss-law closure across this artificial boundary "
                       "is approximate.";
            logRoot(communicator_, message.str());
          }
      }
    if (globalInvalid)
      throw std::out_of_range(
        "Boosted bunch does not fit inside the configured direct E/B box");
    if (globalInCpml)
      throw std::out_of_range(
        "Boosted bunch overlaps the non-physical CPML particle region; enlarge the mesh or reduce the CPML thickness");
  }
}
