#ifndef DIRECT_EB_RUN_METADATA_H
#define DIRECT_EB_RUN_METADATA_H

#include <string>

#include <mpi.h>

#include "config.h"

namespace aprl
{
  struct RunMetadata
  {
    std::string runId;
    std::string configurationPath;
    std::string configurationDigest;
    std::string manifestPath;
    std::string sourceRevision;
    std::string buildDescription;
    std::string mpiLibrary;
    std::string inputParticleIdentity;
    int mpiSize;

    RunMetadata();
  };

  /* Creates one identity on rank zero, broadcasts it, and writes the exact
   * input card plus build/input provenance to the root-only manifest. */
  RunMetadata initializeRunMetadata(const SimulationConfig& config,
                                    MPI_Comm communicator);
}

#endif
