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
#include "lab_detectors.h"
#include "particle_boundary.h"
#include "trajectory_writer.h"

namespace fel
{
  /* Direct SI E/B application path. It deliberately owns no A/phi state. */
  class Simulation
  {
  public:
    Simulation(const SimulationConfig& config, MPI_Comm communicator);

    void solve();

  private:
    struct Slab
    {
      std::size_t offset;
      std::size_t cells;
    };

    void initialize();
    void initializeGeometry();
    void initializeParticles();
    void initializeSources();
    void initializeParticleSelfField();
    void initializeFieldDetectorRegions();
    void validateBeamlineExclusionRules() const;
    void calibrateResourceEstimate();
    void reportResourceEstimate();
    void reportResourceProgress(const char* phase);
    void initializeTrajectoryOutput();
    void initializeDetectorOutput();
    void sampleTrajectory();
    void appendTrajectoryEvent(const RelativisticParticleSI& particle,
                               Double timeBox,
                               TrajectoryEvent event,
                               ParticleBoundaryFace face);
    void finalizeTrajectoryOutput(bool completed);

    void redistributeParticles();
    void pushDepositAndMigrate();
    void reportParticleBoundaryLosses() const;
    bool synchronizedStopRequested();
    bool configuredStopReached();

    Slab slabForRank(int rank) const;
    int ownerRank(Double boxZ) const;
    Double firstBeamlinePhysicalEntranceLab() const;
    Double firstBeamlineInteractionEntranceLab() const;
    Double lastBeamlineInteractionExitLab() const;
    void validateParticlesInsideGlobalBox() const;

    const SimulationConfig& config_;
    DetectorConfig detectorConfig_;
    std::vector<BeamlineElementExtent> beamlineElements_;
    MPI_Comm communicator_;
    int rank_;
    int size_;

    EBGridGeometry globalGeometry_;
    EBGridGeometry localGeometry_;
    std::size_t localZOffset_;
    FieldVector<Double> globalOriginBox_;
    FieldVector<Double> localOriginBox_;
    BoostFrameTransform frame_;
    SIFieldSourceSet sources_;
    std::unique_ptr<EBFieldGrid> fields_;
    std::unique_ptr<EBZSlabHaloExchange> halo_;
    std::unique_ptr<EBMaxwellIncidentInjector> incident_;
    std::vector<RelativisticParticleSI> particles_;
    std::unique_ptr<ParticleCPMLRegion> particleCPML_;
    std::vector<ParticleCPMLCarrier> pmlCarriers_;
    /* The same compact record is reused for experimental retirement-current
     * carriers; the vector identity keeps CPML and retirement damping paths
     * separate. */
    std::vector<ParticleCPMLCarrier> retirementCarriers_;
    std::unique_ptr<ParticleOpenBoundary> particleBoundary_;
    std::unique_ptr<LabDetectorManager> detectors_;

    TrajectoryWriter trajectoryWriter_;
    Double trajectoryRhythmSI_;
    Double nextTrajectorySampleTime_;
    unsigned int trajectorySamplesSinceFlush_;

    Double timeBoxSI_;
    Double totalTimeBoxSI_;
    std::size_t step_;
    bool interrupted_;
    bool configuredStopReached_;
    Double loopWallStart_;
    Double estimatedStepSeconds_;
    std::size_t estimatedMaximumSteps_;
    std::uint64_t modeledLocalPeakBytes_;
    unsigned long long lostParticles_;
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
  };
}

#endif
