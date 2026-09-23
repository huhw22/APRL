#ifndef TRAJECTORY_RADIATION_READER_H
#define TRAJECTORY_RADIATION_READER_H

#include <vector>

#include "mpi.h"

#include "radiation_types.h"

namespace radiation
{
  std::vector<ParticleTrajectory> loadShotDistributed(
      const ShotConfig& shot,
      bool requireComplete,
      std::size_t chunkRecords,
      MPI_Comm communicator,
      TrajectoryLoadStats& globalStats);
}

#endif
