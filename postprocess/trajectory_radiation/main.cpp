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
#include "time_average.h"
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

  bool referencesMatchAngularGrid(
      const radiation::RadiationConfig& config,
      const std::vector<radiation::ObservationDirection>& references)
  {
    const std::vector<double> thetaX = config.thetaX.values();
    const std::vector<double> thetaY = config.thetaY.values();
    if (references.size() != thetaX.size() * thetaY.size()) return false;
    for (std::size_t y = 0; y < thetaY.size(); ++y)
      for (std::size_t x = 0; x < thetaX.size(); ++x)
        {
          const radiation::ObservationDirection& reference =
            references[y * thetaX.size() + x];
          if (reference.thetaX != thetaX[x] ||
              reference.thetaY != thetaY[y]) return false;
        }
    return true;
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
      std::unique_ptr<radiation::TimeAverageAccumulator> timeAverage;
      const std::vector<double> timeWindowCenters =
        config.timeAverage.windowCenters();
      const std::vector<long double> timeAverageOmega =
        radiation::photonEnergyToOmega(
          config.timeAverage.photonEnergyEV);
      std::vector<radiation::ObservationDirection> timeAverageReferences;
      for (std::size_t reference = 0;
           reference < config.timeAverage.referenceAngles.size(); ++reference)
        timeAverageReferences.push_back(radiation::makeObservationDirection(
          config,
          static_cast<double>(
            config.timeAverage.referenceAngles[reference][0]),
          static_cast<double>(
            config.timeAverage.referenceAngles[reference][1])));
      const bool timeReferencesAreAngularGrid =
        referencesMatchAngularGrid(config, timeAverageReferences);
      if (rank == 0)
        {
          writer.reset(new radiation::FarFieldWriter(config));
          if (config.timeAverage.enabled)
            timeAverage.reset(new radiation::TimeAverageAccumulator(config));
          std::cout << "Trajectory far-field reconstruction: MPI ranks="
                    << size << ", shots=" << config.shots.size()
                    << ", frequencies=" << config.photonEnergyEV.count
                    << ", angles=" << config.thetaX.count << "x"
                    << config.thetaY.count << std::endl;
          std::cout << "Internal-kink endpoint formulation active; "
                    << "artificial first/last trajectory endpoints are suppressed."
                    << std::endl;
          if (timeAverage)
            std::cout << "Time averaging: Hann windows="
                      << timeWindowCenters.size()
                      << " per shot, frequencies="
                      << timeAverageOmega.size() << ", references="
                      << timeAverageReferences.size()
                      << ", accumulator="
                      << static_cast<double>(timeAverage->estimatedBytes()) /
                         (1024.0 * 1024.0) << " MiB" << std::endl;
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

          const radiation::ObservationDirection centralObservation =
            radiation::makeObservationDirection(config, 0.0, 0.0);
          const radiation::ObserverTimeRange localObserverRange =
            radiation::internalKnotObserverTimeRange(trajectories,
              centralObservation, config.minimumRecordsPerParticle);
          radiation::ObserverTimeRange observerRange;
          MPI_Reduce(&localObserverRange.minimum, &observerRange.minimum, 1,
            MPI_LONG_DOUBLE, MPI_MIN, 0, communicator);
          MPI_Reduce(&localObserverRange.maximum, &observerRange.maximum, 1,
            MPI_LONG_DOUBLE, MPI_MAX, 0, communicator);
          MPI_Reduce(&localObserverRange.internalKnots,
            &observerRange.internalKnots, 1, MPI_UNSIGNED_LONG_LONG,
            MPI_SUM, 0, communicator);
          if (rank == 0 && observerRange.internalKnots > 0)
            {
              std::cout << "Shot " << config.shots[shot].name
                        << ": central-axis internal-knot u=t-n.r/c range [s]="
                        << static_cast<double>(observerRange.minimum) << " to "
                        << static_cast<double>(observerRange.maximum)
                        << ", knots=" << observerRange.internalKnots
                        << std::endl;
              if (config.timeAverage.enabled &&
                  (config.timeAverage.startTime < observerRange.minimum ||
                   config.timeAverage.endTime > observerRange.maximum))
                std::cout << "Warning: configured time_average interval "
                          << "extends beyond the central-axis internal-knot "
                          << "range for this shot; edge windows may be "
                          << "underfilled." << std::endl;
            }

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

          if (config.timeAverage.enabled)
            for (std::size_t window = 0;
                 window < timeWindowCenters.size(); ++window)
              {
                unsigned long long unusedSkipped = 0;
                std::vector<long double> localReference =
                  radiation::calculateWindowedRadiation(trajectories,
                    config, timeAverageOmega, timeAverageReferences,
                    static_cast<long double>(timeWindowCenters[window]),
                    static_cast<long double>(
                      config.timeAverage.windowDuration), unusedSkipped);
                std::vector<long double> reducedReference;
                if (rank == 0)
                  reducedReference.assign(localReference.size(), 0.0L);
                MPI_Reduce(localReference.empty() ? NULL :
                    &localReference[0],
                  rank == 0 && !reducedReference.empty() ?
                    &reducedReference[0] : NULL,
                  checkedCount(localReference.size()), MPI_LONG_DOUBLE,
                  MPI_SUM, 0, communicator);
                if (rank == 0)
                  timeAverage->beginSample(reducedReference);

                for (std::size_t thetaYOffset = 0;
                     thetaYOffset < config.thetaY.count;
                     thetaYOffset += config.thetaYBlock)
                  {
                    const std::size_t thetaYCount = std::min(
                      config.thetaYBlock,
                      config.thetaY.count - thetaYOffset);
                    const std::vector<radiation::ObservationDirection>
                      directions = radiation::makeObservationDirections(
                        config, thetaYOffset, thetaYCount);
                    std::vector<long double> reducedTarget;
                    if (timeReferencesAreAngularGrid)
                      {
                        if (rank == 0)
                          {
                            const std::size_t localDirections =
                              thetaYCount * config.thetaX.count;
                            reducedTarget.assign(timeAverageOmega.size() *
                              localDirections * 4, 0.0L);
                            const std::size_t globalDirections =
                              config.thetaY.count * config.thetaX.count;
                            for (std::size_t frequency = 0;
                                 frequency < timeAverageOmega.size();
                                 ++frequency)
                              for (std::size_t direction = 0;
                                   direction < localDirections; ++direction)
                                for (std::size_t component = 0;
                                     component < 4; ++component)
                                  reducedTarget[
                                    (frequency * localDirections + direction) *
                                      4 + component] = reducedReference[
                                    (frequency * globalDirections +
                                      thetaYOffset * config.thetaX.count +
                                      direction) * 4 + component];
                          }
                      }
                    else
                      {
                        std::vector<long double> localTarget =
                          radiation::calculateWindowedRadiation(trajectories,
                            config, timeAverageOmega, directions,
                            static_cast<long double>(
                              timeWindowCenters[window]),
                            static_cast<long double>(
                              config.timeAverage.windowDuration),
                            unusedSkipped);
                        if (rank == 0)
                          reducedTarget.assign(localTarget.size(), 0.0L);
                        MPI_Reduce(localTarget.empty() ? NULL :
                            &localTarget[0],
                          rank == 0 && !reducedTarget.empty() ?
                            &reducedTarget[0] : NULL,
                          checkedCount(localTarget.size()), MPI_LONG_DOUBLE,
                          MPI_SUM, 0, communicator);
                      }
                    if (rank == 0)
                      timeAverage->accumulateAngularBlock(reducedTarget,
                        thetaYOffset, thetaYCount);
                  }
                if (rank == 0) timeAverage->endSample();
              }
        }

      if (rank == 0)
        {
          if (timeAverage)
            {
              timeAverage->write(writer->fileHandle());
              std::cout << "Time-average samples committed: "
                        << timeAverage->sampleCount() << std::endl;
            }
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
