#include "particle_initializer.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <random>
#include <stdexcept>

namespace fel
{
  namespace
  {
    const Double kPi = 3.141592653589793238462643383279502884;

    Double halton(unsigned int dimension, std::size_t index)
    {
      static const unsigned int primes[20] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37, 41, 43, 47, 53, 59, 61, 67, 71
      };
      if (dimension >= 20)
        throw std::out_of_range("Halton dimension exceeds available primes");
      const unsigned int prime = primes[dimension];
      std::size_t value = index + 1;
      Double denominator = static_cast<Double>(prime);
      Double result = 0.0;
      while (value > 0)
        {
          result += static_cast<Double>(value % prime) / denominator;
          value /= prime;
          denominator *= static_cast<Double>(prime);
        }
      return 1.0 - result;
    }

    Double gaussian(Double first, Double second)
    {
      const Double safe = std::max(first,
        std::numeric_limits<Double>::min());
      return std::sqrt(-2.0 * std::log(safe)) *
             std::sin(2.0 * kPi * second);
    }

    FieldVector<Double> meanProperVelocity(
        const BeamDistributionConfig& config)
    {
      if (!(config.gamma >= 1.0) || !std::isfinite(config.gamma))
        throw std::invalid_argument("Beam gamma must be finite and at least one");
      FieldVector<Double> direction(config.direction);
      const Double norm = direction.norm();
      if (!(norm > 0.0) || !std::isfinite(norm))
        throw std::invalid_argument("Beam direction must be nonzero");
      direction /= norm;
      FieldVector<Double> result(0.0);
      result.mv(BoostFrameTransform::gammaBetaFromGamma(config.gamma),
                direction);
      return result;
    }

    RelativisticParticleSI makeParticle(
        const FieldVector<Double>& position,
        const FieldVector<Double>& properVelocity,
        Double macroElectrons,
        std::uint64_t sourceId = 0)
    {
      if (!(macroElectrons > 0.0) || !std::isfinite(macroElectrons))
        throw std::invalid_argument("Macro-particle electron count must be positive");
      RelativisticParticleSI particle;
      particle.position = position;
      particle.properVelocity = properVelocity;
      particle.charge = -SI::elementaryCharge * macroElectrons;
      particle.mass = SI::electronMass * macroElectrons;
      particle.weight = 1.0;
      particle.sourceId = sourceId;
      return particle;
    }

    Double firstUndulatorPeriod(const SimulationConfig& config)
    {
      for (std::size_t i = 0; i < config.magnets.size(); ++i)
        if (config.magnets[i].type == SIMagnetType::PlanarUndulator)
          return config.magnets[i].period;
      return 0.0;
    }

    void appendEllipsoid(const BeamDistributionConfig& config,
                         const SimulationConfig& simulation,
                         int rank, int size,
                         std::vector<RelativisticParticleSI>& output)
    {
      std::size_t count = config.macroparticles;
      const Double undulatorPeriod = firstUndulatorPeriod(simulation);
      const Double resonantWavelength = undulatorPeriod > 0.0 ?
        undulatorPeriod /
          (2.0 * simulation.mesh.boostGamma * simulation.mesh.boostGamma) :
        0.0;
      const std::size_t groupSize = resonantWavelength > 0.0 ? 4 : 1;
      if (count == 0)
        throw std::invalid_argument("Ellipsoid requires macroparticles > 0");
      if (count % groupSize != 0)
        count += groupSize - count % groupSize;

      const Double macroElectrons = config.electrons /
        static_cast<Double>(count);
      const FieldVector<Double> mean = meanProperVelocity(config);
      const std::size_t baseCount = count / groupSize;

      std::vector<Double> randomValues;
      if (config.generator == "random")
        {
          std::mt19937_64 engine(config.randomSeed);
          std::uniform_real_distribution<Double> uniform(
            std::numeric_limits<Double>::epsilon(),
            1.0 - std::numeric_limits<Double>::epsilon());
          randomValues.resize(baseCount * config.positions.size() * 20);
          for (std::size_t i = 0; i < randomValues.size(); ++i)
            randomValues[i] = uniform(engine);
        }
      else if (config.generator != "halton")
        throw std::invalid_argument("Beam generator must be halton or random");

      const auto sample = [&](unsigned int dimension, std::size_t index) {
        return config.generator == "random" ?
          randomValues[index * 20 + dimension] : halton(dimension, index);
      };

      for (std::size_t copy = 0; copy < config.positions.size(); ++copy)
        for (std::size_t base = static_cast<std::size_t>(rank);
             base < baseCount; base += static_cast<std::size_t>(size))
          {
            const std::size_t sequence = copy * baseCount + base;
            FieldVector<Double> displacement(0.0);
            if (config.profile == "uniform")
              for (unsigned int axis = 0; axis < 3; ++axis)
                displacement[axis] = config.sigmaPosition[axis] *
                  (2.0 * sample(axis, sequence) - 1.0);
            else if (config.profile == "gaussian")
              {
                displacement[0] = config.sigmaPosition[0] *
                  gaussian(sample(0, sequence), sample(1, sequence));
                displacement[1] = config.sigmaPosition[1] *
                  gaussian(sample(0, sequence), 0.25 + sample(1, sequence));
                displacement[2] = config.sigmaPosition[2] *
                  gaussian(sample(2, sequence), sample(3, sequence));
              }
            else
              throw std::invalid_argument("Beam profile must be gaussian or uniform");

            if (std::abs(displacement[0]) >= config.transverseCutoff ||
                std::abs(displacement[1]) >= config.transverseCutoff ||
                std::abs(displacement[2]) >= config.longitudinalCutoff)
              continue;

            FieldVector<Double> momentum(mean);
            momentum[0] += config.sigmaProperVelocity[0] *
              gaussian(sample(4, sequence), sample(5, sequence));
            momentum[1] += config.sigmaProperVelocity[1] *
              gaussian(sample(4, sequence), 0.25 + sample(5, sequence));
            momentum[2] += config.sigmaProperVelocity[2] *
              gaussian(sample(6, sequence), sample(7, sequence));

            for (std::size_t member = 0; member < groupSize; ++member)
              {
                FieldVector<Double> position(config.positions[copy]);
                position += displacement;
                if (resonantWavelength > 0.0)
                  {
                    position[2] -= resonantWavelength *
                      static_cast<Double>(member) /
                      static_cast<Double>(groupSize);
                    if (config.shotNoise)
                      {
                        Double rmsBunching = config.bunchingFactor;
                        if (rmsBunching == 0.0)
                          {
                            const Double electronsPerWavelength =
                              config.electrons * resonantWavelength /
                              (2.0 * config.sigmaPosition[2]);
                            if (!(electronsPerWavelength > 0.0))
                              throw std::invalid_argument(
                                "Shot noise requires positive longitudinal spread");
                            rmsBunching = 1.0 /
                              std::sqrt(electronsPerWavelength);
                          }
                        const Double noise = rmsBunching *
                          std::sqrt(-2.0 * std::log(
                          std::max(sample(8, sequence),
                                   std::numeric_limits<Double>::min())));
                        position[2] -= resonantWavelength / kPi * noise *
                          std::sin(2.0 * kPi * position[2] /
                                   resonantWavelength +
                                   2.0 * kPi * sample(9, sequence));
                      }
                    else if (config.bunchingFactor != 0.0)
                      position[2] -= resonantWavelength / kPi *
                        config.bunchingFactor * std::sin(
                          2.0 * kPi * position[2] / resonantWavelength +
                          config.bunchingPhase);
                  }
                output.push_back(makeParticle(position, momentum,
                                              macroElectrons));
              }
          }
    }

    void appendFile(const BeamDistributionConfig& config,
                    const SimulationConfig& simulation,
                    int rank, int size, MPI_Comm communicator,
                    std::vector<RelativisticParticleSI>& output)
    {
      if (config.file.empty())
        throw std::invalid_argument("File beam requires a file path");
      std::ifstream input(config.file.c_str());
      if (!input)
        throw std::runtime_error("Cannot open beam particle file: " + config.file);
      if (config.macroparticles == 0)
        throw std::invalid_argument("File beam requires macroparticles > 0");
      const Double macroElectrons = config.electrons /
        static_cast<Double>(config.macroparticles);
      FieldVector<Double> position(0.0);
      FieldVector<Double> momentum(0.0);
      std::uint64_t line = 0;
      while (input >> position[0] >> position[1] >> position[2]
                   >> momentum[0] >> momentum[1] >> momentum[2])
        {
          ++line;
          if ((line - 1) % static_cast<std::uint64_t>(size) !=
              static_cast<std::uint64_t>(rank))
            continue;
          for (std::size_t copy = 0; copy < config.positions.size(); ++copy)
            {
              FieldVector<Double> shifted(position);
              shifted *= simulation.inputUnits.length;
              shifted += config.positions[copy];
              output.push_back(makeParticle(shifted, momentum,
                                            macroElectrons, line));
            }
        }
      unsigned long long local = static_cast<unsigned long long>(line);
      unsigned long long minimum = 0;
      unsigned long long maximum = 0;
      MPI_Allreduce(&local, &minimum, 1, MPI_UNSIGNED_LONG_LONG,
                    MPI_MIN, communicator);
      MPI_Allreduce(&local, &maximum, 1, MPI_UNSIGNED_LONG_LONG,
                    MPI_MAX, communicator);
      if (minimum != maximum || maximum != config.macroparticles)
        throw std::runtime_error("Beam file record count does not match macroparticles");
    }

    void appendManual(const BeamDistributionConfig& config,
                      int rank, int size,
                      std::vector<RelativisticParticleSI>& output)
    {
      const FieldVector<Double> momentum = meanProperVelocity(config);
      for (std::size_t i = static_cast<std::size_t>(rank);
           i < config.positions.size(); i += static_cast<std::size_t>(size))
        output.push_back(makeParticle(config.positions[i], momentum,
                                      config.electrons));
    }

    void appendCrystal(const BeamDistributionConfig& config,
                       int rank, int size,
                       std::vector<RelativisticParticleSI>& output)
    {
      const std::size_t sites =
        static_cast<std::size_t>(config.latticeCounts[0]) *
        static_cast<std::size_t>(config.latticeCounts[1]) *
        static_cast<std::size_t>(config.latticeCounts[2]);
      if (sites == 0 || config.macroparticles == 0 ||
          config.macroparticles % sites != 0)
        throw std::invalid_argument(
          "3D crystal macroparticles must be a positive multiple of lattice sites");
      const std::size_t perSite = config.macroparticles / sites;
      const Double macroElectrons = config.electrons /
        static_cast<Double>(config.macroparticles);
      const FieldVector<Double> mean = meanProperVelocity(config);

      for (std::size_t copy = 0; copy < config.positions.size(); ++copy)
        for (unsigned int i = 0; i < config.latticeCounts[0]; ++i)
          for (unsigned int j = 0; j < config.latticeCounts[1]; ++j)
            for (unsigned int k = 0; k < config.latticeCounts[2]; ++k)
              for (std::size_t member = 0; member < perSite; ++member)
                {
                  const std::size_t sequence =
                    ((((copy * config.latticeCounts[0] + i) *
                       config.latticeCounts[1] + j) *
                       config.latticeCounts[2] + k) * perSite + member);
                  if (sequence % static_cast<std::size_t>(size) !=
                      static_cast<std::size_t>(rank))
                    continue;
                  FieldVector<Double> position(config.positions[copy]);
                  position[0] += (static_cast<Double>(i) + 1.0 -
                    0.5 * config.latticeCounts[0]) * config.latticeConstants[0];
                  position[1] += (static_cast<Double>(j) + 1.0 -
                    0.5 * config.latticeCounts[1]) * config.latticeConstants[1];
                  position[2] += (static_cast<Double>(k) + 1.0 -
                    0.5 * config.latticeCounts[2]) * config.latticeConstants[2];
                  FieldVector<Double> momentum(mean);
                  for (unsigned int axis = 0; axis < 3; ++axis)
                    {
                      position[axis] += 0.5 * config.sigmaPosition[axis] *
                        gaussian(halton(2 * axis, sequence),
                                 halton(2 * axis + 1, sequence));
                      momentum[axis] += config.sigmaProperVelocity[axis] *
                        gaussian(halton(6 + 2 * axis, sequence),
                                 halton(7 + 2 * axis, sequence));
                    }
                  output.push_back(makeParticle(position, momentum,
                                                macroElectrons));
                }
    }
  }

  std::vector<RelativisticParticleSI> ParticleInitializer::create(
      const SimulationConfig& config, MPI_Comm communicator)
  {
    if (communicator == MPI_COMM_NULL)
      throw std::invalid_argument("Particle initializer communicator is null");
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);
    std::vector<RelativisticParticleSI> result;
    for (std::size_t i = 0; i < config.beam.size(); ++i)
      {
        switch (config.beam[i].type)
          {
          case BeamDistributionType::Ellipsoid:
            appendEllipsoid(config.beam[i], config, rank, size, result);
            break;
          case BeamDistributionType::File:
            appendFile(config.beam[i], config, rank, size,
                       communicator, result);
            break;
          case BeamDistributionType::Manual:
            appendManual(config.beam[i], rank, size, result);
            break;
          case BeamDistributionType::Crystal3D:
            appendCrystal(config.beam[i], rank, size, result);
            break;
          }
      }
    return result;
  }
}
