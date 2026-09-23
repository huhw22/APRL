#include <algorithm>
#include <climits>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "hdf5.h"
#include "mpi.h"

#include "configuration.h"
#include "output_writer.h"
#include "radiation_kernel.h"
#include "trajectory_reader.h"

namespace
{
  int checkedCount(std::size_t count)
  {
    if (count > static_cast<std::size_t>(INT_MAX))
      throw std::overflow_error(
        "A radiation block exceeds the MPI reduction count limit");
    return static_cast<int>(count);
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm communicator = MPI_COMM_WORLD;
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &size);
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  try
    {
      if (argc != 2)
        throw std::invalid_argument(
          "Usage: trajectory_radiation <radiation.yaml>");
      const radiation::RadiationConfig config =
        radiation::loadConfiguration(argv[1]);
      const std::vector<double> photonEnergy =
        config.photonEnergyEV.values();
      const std::vector<long double> omega =
        radiation::photonEnergyToOmega(photonEnergy);
      std::unique_ptr<radiation::FarFieldWriter> writer;
      if (rank == 0)
        {
          writer.reset(new radiation::FarFieldWriter(config));
          std::cout << "Trajectory far-field reconstruction: MPI ranks="
                    << size << ", shots=" << config.shots.size()
                    << ", frequencies=" << config.photonEnergyEV.count
                    << ", angles=" << config.thetaX.count << "x"
                    << config.thetaY.count << std::endl;
          std::cout << "Internal-kink endpoint formulation active; "
                    << "artificial first/last trajectory endpoints are suppressed."
                    << std::endl;
        }

      for (std::size_t shot = 0; shot < config.shots.size(); ++shot)
        {
          radiation::TrajectoryLoadStats stats;
          std::vector<radiation::ParticleTrajectory> trajectories =
            radiation::loadShotDistributed(config.shots[shot],
              config.requireComplete, config.readChunkRecords,
              communicator, stats);
          if (rank == 0)
            std::cout << "Shot " << config.shots[shot].name
                      << ": input records=" << stats.inputRecords
                      << ", unique records=" << stats.uniqueRecords
                      << ", particles=" << stats.particles
                      << ", terminal events=" << stats.terminalEvents
                      << ", duplicate timestamps="
                      << stats.duplicateRecords << std::endl;

          unsigned long long localSkippedMaximum = 0;
          for (std::size_t frequencyOffset = 0;
               frequencyOffset < config.photonEnergyEV.count;
               frequencyOffset += config.frequencyBlock)
            for (std::size_t thetaYOffset = 0;
                 thetaYOffset < config.thetaY.count;
                 thetaYOffset += config.thetaYBlock)
              {
                const std::size_t frequencyCount = std::min(
                  config.frequencyBlock,
                  config.photonEnergyEV.count - frequencyOffset);
                const std::size_t thetaYCount = std::min(
                  config.thetaYBlock,
                  config.thetaY.count - thetaYOffset);
                unsigned long long skipped = 0;
                radiation::RadiationBlock local =
                  radiation::calculateRadiationBlock(trajectories,
                    config, omega, frequencyOffset, frequencyCount,
                    thetaYOffset, thetaYCount, skipped);
                localSkippedMaximum = std::max(localSkippedMaximum, skipped);
                radiation::RadiationBlock reduced = local;
                if (rank == 0)
                  reduced.amplitude.assign(local.amplitude.size(), 0.0L);
                MPI_Reduce(local.amplitude.empty() ? NULL :
                    &local.amplitude[0],
                  rank == 0 && !reduced.amplitude.empty() ?
                    &reduced.amplitude[0] : NULL,
                  checkedCount(local.amplitude.size()), MPI_LONG_DOUBLE,
                  MPI_SUM, 0, communicator);
                if (rank == 0)
                  writer->writeShotBlock(shot, reduced, omega);
              }
          unsigned long long globalSkipped = 0;
          MPI_Reduce(&localSkippedMaximum, &globalSkipped, 1,
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, communicator);
          if (rank == 0 && globalSkipped > 0)
            std::cout << "Shot " << config.shots[shot].name
                      << ": skipped particles with fewer than "
                      << config.minimumRecordsPerParticle
                      << " records=" << globalSkipped << std::endl;
        }

      if (rank == 0)
        {
          writer->finalizeEnsemble();
          writer->close();
          std::cout << "Far-field output committed: "
                    << config.outputFile << std::endl;
        }
      MPI_Barrier(communicator);
      MPI_Finalize();
      return EXIT_SUCCESS;
    }
  catch (const std::exception& error)
    {
      std::cerr << "trajectory_radiation rank " << rank
                << " error: " << error.what() << std::endl;
      MPI_Abort(communicator, EXIT_FAILURE);
      return EXIT_FAILURE;
    }
}
