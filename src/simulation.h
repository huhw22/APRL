#ifndef DIRECT_EB_SIMULATION_H
#define DIRECT_EB_SIMULATION_H

#include <cstddef>
#include <memory>
#include <vector>

#include <mpi.h>

#include "config.h"
#include "eb_bunch.h"
#include "eb_deposition.h"
#include "eb_incident.h"
#include "eb_mpi.h"
#include "eb_particles.h"
#include "eb_sources.h"
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
    void sampleTrajectory();
    void finalizeTrajectoryOutput(bool completed);

    void redistributeParticles();
    void pushDepositAndMigrate();
    bool synchronizedStopRequested();

    Slab slabForRank(int rank) const;
    int ownerRank(Double boxZ) const;
    Double firstMagneticEntranceLab() const;
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

    TrajectoryWriter trajectoryWriter_;
    Double trajectoryRhythmSI_;
    Double nextTrajectorySampleTime_;
    unsigned int trajectorySamplesSinceFlush_;

    Double timeBoxSI_;
    Double totalTimeBoxSI_;
    std::size_t step_;
    bool interrupted_;
  };
}

#endif
