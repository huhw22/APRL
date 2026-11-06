#ifndef DIRECT_EB_SIMULATION_H
#define DIRECT_EB_SIMULATION_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <mpi.h>

#include "config.h"
#include "eb_bunch.h"
#include "eb_deposition.h"
#include "eb_incident.h"
#include "eb_initial_field.h"
#include "eb_mpi.h"
#include "eb_particles.h"
#include "eb_sources.h"
#include "energy_ledger.h"
#include "particle_energy_audit.h"
#include "lab_detectors.h"
#include "longitudinal_mesh.h"
#include "particle_boundary.h"
#include "run_metadata.h"
#include "trajectory_writer.h"
#include "simulation_moving_window.h"
#include "source_coupling.h"

namespace aprl
{
  class CheckpointArchive;
  /* Direct SI E/B application path. It deliberately owns no A/phi state. */
  class Simulation
  {
  public:
    Simulation(const SimulationConfig& config, MPI_Comm communicator);
    ~Simulation();

    void solve();

  private:
    struct ParticleWorkspace;
    enum PerformanceStage
    {
      FieldStage, ParticlePrepareStage, ParticleDepositStage,
      ParticleMigrationStage, DiagnosticsStage, ControlStage,
      FieldCommunicationStage, ParticleCommunicationStage,
      PendingAdmissionStage, MovingWindowStage, PendingStatisticsStage,
      CouplingTotalStage, CouplingCacheStage, CouplingResidualStage,
      PerformanceStageCount
    };
    struct Slab
    {
      std::size_t offset;
      std::size_t cells;
    };

    struct BoundaryInjectionWorkSample
    {
      FieldVector<Double> midpoint;
      FieldVector<Double> chargeDisplacement;
      Double oldElectricDot;
    };

    void initialize();
    void initializeMovingWindow();
    void installMovingWindowBoundary();
    static void completeMovingInletLayer(void*,EBFieldGrid&,std::size_t,std::ptrdiff_t);
    void initializeMovingWindowFields();
    void admitMovingWindowParticles();
    std::size_t pendingParticleCount() const;
    std::size_t storedParticleCount() const;
    EBElectricHaloView beginMovingWindowStep();
    void finishMovingWindowStep();
    void applyMovingTransverseReferenceCurrent();
    void finishMovingTransverseReferenceCurrent();
    void finalizeMovingWindow();
    bool firstMovingWindowEntry(const FieldVector<Double>& start,
        const FieldVector<Double>& end, Double startFraction,
        ParticleBoundaryHit& hit) const;
    void recordMovingWindowExit(Double charge);
    void queueMovingConstraintWork(const FieldVector<Double>& start,
        const FieldVector<Double>& end, Double charge, bool carrier);
    void finalizeMovingConstraintWork();
    static Double movingWindowFieldValue(void* context, unsigned component,
        std::size_t i, std::size_t j, std::ptrdiff_t k, std::size_t step);
    void loadCheckpointGeometry();
    void restoreCheckpointFields();
    void restoreCheckpointDiagnostics();
    void maybeWriteCheckpoint();
    void writeCheckpoint(const std::string& path);
    void transferCheckpointScalars(std::vector<unsigned char>& bytes,
                                   bool reading, bool geometry);
    void validateCheckpointBranch() const;
    void preflightOutputPaths() const;
    void initializeGeometry();
    void validateAndReportRadiationResolution() const;
    void initializeParticles();
    void reportDurationEstimate() const;
    void initializeSources();
    void initializeParticleSelfField();
    void initializeFieldDetectorRegions();
    void validateBeamlineExclusionRules() const;
    void calibrateResourceEstimate();
    void reportResourceEstimate();
    void reportResourceProgress(const char* phase);
    Double performanceNow() const;
    void reportPerformance(bool final = false);
    void initializeTrajectoryOutput();
    void initializeDetectorOutput();
    void initializeEnergyLedger();
    void advanceEnergyLedgerStep();
    void finalizeBoundaryInjectionWorkStep();
    void sampleEnergyLedger(bool force = false,
                            bool guardFinalCollectives = false);
    void finalizeEnergyLedger(bool completed);
    long double localInteriorFieldEnergy() const;
    void localInteriorBoundaryPower(long double power[6]) const;
    void sampleTrajectory();
    void appendTrajectoryEvent(const RelativisticParticleSI& particle,
                               Double timeBox,
                               TrajectoryEvent event,
                               ParticleBoundaryFace face);
    void finalizeTrajectoryOutput(bool completed);

    void redistributeParticles();
    void pushDepositAndMigrate();
    void solveConservativeSourceStep();
    void commitConservativeElectric();
    bool maybeBalanceParticles();
    bool conservativeLongitudinalHandoff(const ParticleBoundaryHit&,
        const FieldVector<Double>& fullEnd) const;
    void pushParticle(RelativisticParticleSI& p,Double time,Double dt,
        Double* prescribed,BorisWork* work,ParticleSubstepObserver* observer) {
      if(sourceCoupling_)sourceCoupling_->push(p,sources_,frame_,time,dt,particleSubsteps_,prescribed,work,observer);
      else RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(p,*fields_,localOriginBox_,sources_,frame_,time,dt,particleSubsteps_,prescribed,config_.mesh.selfFieldFeedback,work,observer);
    }
    void reportParticleBoundaryLosses() const;
    bool synchronizedStopRequested();
    bool configuredStopReached();

    Slab slabForRank(int rank) const;
    int ownerRank(Double boxZ) const;
    Double firstMagneticPhysicalEntranceLab() const;
    Double firstMagneticInteractionEntranceLab() const;
    Double lastBeamlineInteractionExitLab() const;
    void validateParticlesInsideGlobalBox() const;

    SimulationConfig config_;
    DetectorConfig detectorConfig_;
    std::vector<BeamlineElementExtent> beamlineElements_;
    MPI_Comm communicator_;
    int rank_;
    int size_;
    RunMetadata runMetadata_;
    LongitudinalLayoutReport longitudinalLayoutReport_;
    bool hasLongitudinalLayoutReport_;

    EBGridGeometry globalGeometry_;
    EBGridGeometry localGeometry_;
    std::size_t localZOffset_;
    std::vector<std::size_t> slabCuts_;
    std::size_t particleBalanceEvents_ = 0;
    double particleBalanceSeconds_ = 0.;
    FieldVector<Double> globalOriginBox_;
    FieldVector<Double> localOriginBox_;
    Double referenceCenterBoxZ_;
    BoostFrameTransform frame_;
    SIFieldSourceSet sources_;
    std::unique_ptr<EBFieldGrid> fields_;
    std::unique_ptr<SimulationMovingWindowState> movingWindow_;
    std::unique_ptr<EBZSlabHaloExchange> halo_;
    std::unique_ptr<EBMaxwellIncidentInjector> incident_;
    std::vector<RelativisticParticleSI> particles_;
    std::unique_ptr<ParticleWorkspace> particleWorkspace_;
    std::unique_ptr<ParticleCPMLRegion> particleCPML_;
    std::vector<ParticleCPMLCarrier> pmlCarriers_;
    /* The same compact record is reused for experimental retirement-current
     * carriers; the vector identity keeps CPML and retirement damping paths
     * separate. */
    std::vector<ParticleCPMLCarrier> retirementCarriers_;
    std::unique_ptr<ParticleOpenBoundary> particleBoundary_;
    std::unique_ptr<LabDetectorManager> detectors_;
    std::unique_ptr<ParticleWorkGeometryWriter> particleWorkGeometry_;
    std::unique_ptr<ConservativeSourceCoupling> sourceCoupling_;

    TrajectoryWriter trajectoryWriter_;
    EnergyLedgerWriter energyLedgerWriter_;
    EnergyLedgerRecord lastEnergyLedgerRecord_;
    Double trajectoryRhythmSI_;
    Double nextTrajectorySampleTime_;
    unsigned int trajectorySamplesSinceFlush_;
    unsigned int particleSubsteps_;
    long double energyLedgerPreviousPower_[6];
    long double energyLedgerOutwardEnergy_[6];
    long double energyLedgerPrescribedWork_;
    long double energyLedgerBoundarySourceWork_;
    long double energyLedgerBoundarySourceAbsWork_;
    long double energyLedgerRemovedKinetic_;
    long double energyLedgerRemovedTotal_;
    BorisWork energyLedgerBorisWork_;
    ParticleEnergyMoment energyLedgerRemovedMoment_;
    ParticleEnergyMoment energyLedgerReleasedMoment_;
    ParticleEnergyMoment energyLedgerInitialMoment_;
    long double energyLedgerReleasedCharge_;
    unsigned long long energyLedgerPartialWorkSubsteps_;
    long double energyLedgerInitialKinetic_;
    long double energyLedgerInitialField_;
    unsigned long long energyLedgerRemovedMacroparticles_;
    std::size_t energyLedgerLastSampleStep_;
    bool energyLedgerHasReference_;

    Double timeBoxSI_;
    Double totalTimeBoxSI_;
    std::size_t step_;
    bool interrupted_;
    bool configuredStopReached_;
    Double loopWallStart_;
    Double estimatedStepSeconds_;
    std::size_t estimatedMaximumSteps_;
    std::uint64_t modeledLocalPeakBytes_;
    std::size_t performanceInterval_;
    std::size_t performanceLastStep_;
    Double performanceWallStart_;
    Double performanceSeconds_[PerformanceStageCount];
    unsigned long long localLostParticles_;
    unsigned long long lostParticles_;
    unsigned long long ballisticBoundaryEntries_;
    unsigned long long frozenPlaneActivations_;
    unsigned long long crossRankActivationPushes_;
    std::vector<BoundaryInjectionWorkSample> boundaryInjectionWorkSamples_;
    unsigned long long cpmlEntryCount_[6];
    Double cpmlEntryCharge_[6];
    unsigned long long directOuterCount_[6];
    Double directOuterCharge_[6];
    unsigned long long carrierOuterCount_[6];
    Double carrierOuterCharge_[6];
    std::size_t peakPmlCarriers_;
    unsigned long long retirementEntryCount_;
    Double retirementEntryCharge_;
    unsigned long long retirementExitCount_;
    Double retirementExitResidualCharge_;
    std::size_t peakRetirementCarriers_;
    std::string stopReason_;
    // Allocated/opened only for explicitly requested checkpoint/restart I/O.
    bool checkpointEnabled_;
    std::size_t nextCheckpointStep_;
    std::size_t nextCheckpointTime_;
    bool retirementChangedOnRestart_;
    Double initialMaximumLabGamma_;
    unsigned long long initialGlobalParticleCount_;
    std::unique_ptr<CheckpointArchive> checkpointReader_;
  };
}

#endif
