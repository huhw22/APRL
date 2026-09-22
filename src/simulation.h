#ifndef DIRECT_EB_SIMULATION_H
#define DIRECT_EB_SIMULATION_H

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <mpi.h>

#include "config.h"
#include "eb_bunch.h"
#include "eb_deposition.h"
#include "eb_incident.h"
#include "eb_mpi.h"
#include "eb_particles.h"
#include "eb_sources.h"
#include "lab_detectors.h"
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
    void initializeTrajectoryOutput();
    void initializeDetectorOutput();
    void sampleTrajectory();
    void finalizeTrajectoryOutput(bool completed);

    void redistributeParticles();
    void pushDepositAndMigrate();
    bool synchronizedStopRequested();
    bool configuredStopReached();

    Slab slabForRank(int rank) const;
    int ownerRank(Double boxZ) const;
    Double firstBeamlinePhysicalEntranceLab() const;
    Double firstBeamlineInteractionEntranceLab() const;
    Double lastBeamlineInteractionExitLab() const;
    void validateParticlesInsideGlobalBox() const;

    const SimulationConfig& config_;
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
    unsigned long long lostParticles_;
    std::string stopReason_;
  };
}

#endif
