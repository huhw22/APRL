#ifndef DIRECT_EB_PARTICLE_INITIALIZER_H
#define DIRECT_EB_PARTICLE_INITIALIZER_H

#include <vector>

#include <mpi.h>

#include "config.h"
#include "eb_particles.h"

namespace fel
{
  struct ParticleInitializationReport
  {
    unsigned long long globalRecords;
    int fileFormatVersion;
    bool laboratoryPlaneCoordinates;
    bool laboratoryPlaneTimeCoordinates;

    ParticleInitializationReport();
  };

  class ParticleInitializer
  {
  public:
    static std::vector<RelativisticParticleSI> create(
        const SimulationConfig& config,
        MPI_Comm communicator,
        ParticleInitializationReport& report);
  };
}

#endif
