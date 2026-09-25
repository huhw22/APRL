#include "particle_initializer.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "boostframe.h"
#include "particle_file.h"

namespace aprl
{
  ParticleInitializationReport::ParticleInitializationReport()
    : globalRecords(0), fileFormatVersion(0),
      laboratoryPlaneCoordinates(false),
      laboratoryPlaneTimeCoordinates(false)
  {}

  namespace
  {
    const Double kPi = 3.141592653589793238462643383279502884;

    std::uint64_t splitmix64(std::uint64_t value)
    {
      value += UINT64_C(0x9e3779b97f4a7c15);
      value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
      value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
      return value ^ (value >> 31);
    }

    Double uniformOpen(std::uint64_t seed, std::uint64_t particle,
                       unsigned int component)
    {
      const std::uint64_t mixed = splitmix64(
        seed ^ splitmix64(particle + UINT64_C(0x100000001b3) * component));
      const std::uint64_t mantissa = (mixed >> 11) | UINT64_C(1);
      const Double value = static_cast<Double>(mantissa) *
        (1.0 / 9007199254740992.0);
      return std::min(value, 1.0 - std::numeric_limits<Double>::epsilon());
    }

    Double gaussian(std::uint64_t seed, std::uint64_t particle,
                    unsigned int pair)
    {
      const Double first = uniformOpen(seed, particle, 2 * pair);
      const Double second = uniformOpen(seed, particle, 2 * pair + 1);
      return std::sqrt(-2.0 * std::log(first)) *
             std::cos(2.0 * kPi * second);
    }

    FieldVector<Double> meanProperVelocity(const BeamInputConfig& config)
    {
      FieldVector<Double> direction(config.direction);
      const Double norm = direction.norm();
      if (!(norm > 0.0) || !std::isfinite(norm))
        throw std::invalid_argument("Generated beam direction must be nonzero");
      direction /= norm;
      FieldVector<Double> result(0.0);
      result.mv(BoostFrameTransform::gammaBetaFromGamma(config.gamma),
                direction);
      return result;
    }

    std::vector<RelativisticParticleSI> generatedGaussian(
        const BeamInputConfig& config, int rank, int size)
    {
      if (config.macroparticles == 0)
        throw std::invalid_argument(
          "generated-gaussian requires macroparticles > 0");
      const std::size_t total = config.macroparticles;
      const std::size_t base = total / static_cast<std::size_t>(size);
      const std::size_t remainder = total % static_cast<std::size_t>(size);
      const std::size_t localCount = base +
        (static_cast<std::size_t>(rank) < remainder ? 1 : 0);
      const std::size_t offset = static_cast<std::size_t>(rank) * base +
        std::min(static_cast<std::size_t>(rank), remainder);
      const Double macroElectrons = config.electrons /
        static_cast<Double>(total);
      const FieldVector<Double> mean = meanProperVelocity(config);

      std::vector<RelativisticParticleSI> particles;
      particles.reserve(localCount);
      for (std::size_t local = 0; local < localCount; ++local)
        {
          const std::uint64_t global =
            static_cast<std::uint64_t>(offset + local);
          RelativisticParticleSI particle;
          for (unsigned int axis = 0; axis < 3; ++axis)
            {
              particle.position[axis] = config.center[axis] +
                config.positionOffset[axis] +
                config.sigmaPosition[axis] *
                gaussian(config.randomSeed, global, axis);
              particle.properVelocity[axis] = mean[axis] +
                config.sigmaProperVelocity[axis] *
                gaussian(config.randomSeed ^ UINT64_C(0xd6e8feb86659fd93),
                         global, axis);
            }
          particle.charge = -SI::elementaryCharge * macroElectrons;
          particle.mass = SI::electronMass * macroElectrons;
          particle.weight = 1.0;
          particle.sourceId = global + 1;
          particles.push_back(particle);
        }
      return particles;
    }
  }

  std::vector<RelativisticParticleSI> ParticleInitializer::create(
      const SimulationConfig& config, MPI_Comm communicator,
      ParticleInitializationReport& report)
  {
    if (communicator == MPI_COMM_NULL)
      throw std::invalid_argument("Particle initializer communicator is null");
    if (config.beam.type == BeamInputType::Hdf5)
      {
        std::vector<RelativisticParticleSI> particles =
          ParticleHdf5File::readDistributed(
          config.beam.file, config.beam.electrons,
          config.beam.positionOffset, communicator,
          report.globalRecords, report.fileFormatVersion);
        report.laboratoryPlaneCoordinates =
          report.fileFormatVersion >=
          ParticleHdf5File::labPlaneOffsetFormatVersion;
        report.laboratoryPlaneTimeCoordinates =
          report.fileFormatVersion >= ParticleHdf5File::formatVersion;
        return particles;
      }

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);
    report.globalRecords = config.beam.macroparticles;
    return generatedGaussian(config.beam, rank, size);
  }
}
