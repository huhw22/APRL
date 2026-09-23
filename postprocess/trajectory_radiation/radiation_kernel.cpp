#include "radiation_kernel.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace radiation
{
  namespace
  {
    Vec3 recordPosition(const TrajectoryRecord& record)
    {
      return Vec3(record.position[0], record.position[1],
                  record.position[2]);
    }

    Vec3 segmentBeta(const TrajectoryRecord& start,
                     const TrajectoryRecord& end)
    {
      const long double duration =
        static_cast<long double>(end.time) -
        static_cast<long double>(start.time);
      if (!(duration > 0.0L) || !std::isfinite(duration))
        throw std::runtime_error(
          "Particle trajectory time is not strictly increasing");
      const Vec3 beta = (recordPosition(end) - recordPosition(start)) /
        (constants::c * duration);
      const long double beta2 = dot(beta, beta);
      if (!(beta2 < 1.0L + 1.0e-10L) || !std::isfinite(beta2))
        throw std::runtime_error(
          "A trajectory chord has a non-physical speed");
      return beta;
    }

    Vec3 transverseFactor(const Vec3& direction, const Vec3& beta)
    {
      const long double denominator = 1.0L - dot(direction, beta);
      if (!(denominator > 64.0L *
            std::numeric_limits<long double>::epsilon()))
        throw std::runtime_error(
          "Far-field direction is singular with a trajectory segment");
      return (dot(direction, beta) * direction - beta) / denominator;
    }
  }

  std::vector<ObservationDirection> makeObservationDirections(
      const RadiationConfig& config, std::size_t thetaYOffset,
      std::size_t thetaYCount)
  {
    const std::vector<double> thetaX = config.thetaX.values();
    const std::vector<double> thetaY = config.thetaY.values();
    if (thetaYOffset + thetaYCount > thetaY.size())
      throw std::out_of_range("Observation theta-y block is outside the grid");
    std::vector<ObservationDirection> result;
    result.reserve(thetaYCount * thetaX.size());
    for (std::size_t localY = 0; localY < thetaYCount; ++localY)
      for (std::size_t x = 0; x < thetaX.size(); ++x)
        {
          result.push_back(makeObservationDirection(config, thetaX[x],
            thetaY[thetaYOffset + localY]));
        }
    return result;
  }

  ObservationDirection makeObservationDirection(
      const RadiationConfig& config, double thetaX, double thetaY)
  {
    const Vec3 slopeDirection = config.observationAxis +
      std::tan(static_cast<long double>(thetaX)) * config.horizontalAxis +
      std::tan(static_cast<long double>(thetaY)) * config.verticalAxis;
    ObservationDirection observation;
    observation.direction = normalized(
      slopeDirection, "observation direction");
    const Vec3 horizontalProjection = config.horizontalAxis -
      dot(config.horizontalAxis, observation.direction) *
        observation.direction;
    observation.horizontal = normalized(
      horizontalProjection, "far-field horizontal polarization");
    observation.vertical = normalized(cross(
      observation.direction, observation.horizontal),
      "far-field vertical polarization");
    observation.thetaX = thetaX;
    observation.thetaY = thetaY;
    return observation;
  }

  std::vector<long double> photonEnergyToOmega(
      const std::vector<double>& photonEnergyEV)
  {
    std::vector<long double> result(photonEnergyEV.size());
    for (std::size_t index = 0; index < result.size(); ++index)
      result[index] = static_cast<long double>(photonEnergyEV[index]) *
        constants::elementaryCharge / constants::hbar;
    return result;
  }

  ObserverTimeRange internalKnotObserverTimeRange(
      const std::vector<ParticleTrajectory>& trajectories,
      const ObservationDirection& observation,
      std::size_t minimumRecordsPerParticle)
  {
    ObserverTimeRange result;
    for (std::size_t particle = 0; particle < trajectories.size(); ++particle)
      {
        const ParticleTrajectory& trajectory = trajectories[particle];
        if (trajectory.records.size() < minimumRecordsPerParticle) continue;
        for (std::size_t knot = 1;
             knot + 1 < trajectory.records.size(); ++knot)
          {
            const TrajectoryRecord& record = trajectory.records[knot];
            const long double reducedObserverTime =
              static_cast<long double>(record.time) -
              dot(observation.direction, recordPosition(record)) /
                constants::c;
            result.minimum = std::min(result.minimum, reducedObserverTime);
            result.maximum = std::max(result.maximum, reducedObserverTime);
            ++result.internalKnots;
          }
      }
    return result;
  }

  RadiationBlock calculateRadiationBlock(
      const std::vector<ParticleTrajectory>& trajectories,
      const RadiationConfig& config,
      const std::vector<long double>& omega,
      std::size_t frequencyOffset,
      std::size_t frequencyCount,
      std::size_t thetaYOffset,
      std::size_t thetaYCount,
      unsigned long long& skippedShortParticles)
  {
    if (frequencyOffset + frequencyCount > omega.size())
      throw std::out_of_range("Radiation frequency block is outside the grid");
    const std::vector<ObservationDirection> directions =
      makeObservationDirections(config, thetaYOffset, thetaYCount);
    RadiationBlock block;
    block.frequencyOffset = frequencyOffset;
    block.frequencyCount = frequencyCount;
    block.thetaYOffset = thetaYOffset;
    block.thetaYCount = thetaYCount;
    block.thetaXCount = config.thetaX.count;
    const std::size_t scalars = frequencyCount * thetaYCount *
      block.thetaXCount * 2 * 2;
    block.amplitude.assign(scalars, 0.0L);
    skippedShortParticles = 0;

    for (std::size_t particle = 0; particle < trajectories.size(); ++particle)
      {
        const ParticleTrajectory& trajectory = trajectories[particle];
        if (trajectory.records.size() < config.minimumRecordsPerParticle)
          {
            ++skippedShortParticles;
            continue;
          }
        std::vector<Vec3> beta(trajectory.records.size() - 1);
        for (std::size_t segment = 0; segment < beta.size(); ++segment)
          beta[segment] = segmentBeta(trajectory.records[segment],
                                     trajectory.records[segment + 1]);
        const long double charge = trajectory.charge;

        for (std::size_t directionIndex = 0;
             directionIndex < directions.size(); ++directionIndex)
          {
            const ObservationDirection& observation =
              directions[directionIndex];
            const std::size_t localY =
              directionIndex / block.thetaXCount;
            const std::size_t x =
              directionIndex % block.thetaXCount;
            for (std::size_t knot = 1;
                 knot + 1 < trajectory.records.size(); ++knot)
              {
                const Vec3 delta = transverseFactor(
                  observation.direction, beta[knot]) -
                  transverseFactor(observation.direction, beta[knot - 1]);
                const long double component[2] = {
                  dot(delta, observation.horizontal),
                  dot(delta, observation.vertical)
                };
                if (component[0] == 0.0L && component[1] == 0.0L)
                  continue;
                const TrajectoryRecord& record = trajectory.records[knot];
                const Vec3 position = recordPosition(record);
                const long double retardedTime =
                  static_cast<long double>(record.time) -
                  dot(observation.direction, position) / constants::c;
                for (std::size_t localFrequency = 0;
                     localFrequency < frequencyCount; ++localFrequency)
                  {
                    const long double phase = std::remainder(
                      omega[frequencyOffset + localFrequency] * retardedTime,
                      2.0L * constants::pi);
                    const long double cosine = std::cos(phase);
                    const long double sine = std::sin(phase);
                    for (std::size_t polarization = 0;
                         polarization < 2; ++polarization)
                      {
                        const long double magnitude =
                          charge * component[polarization];
                        block.amplitude[block.scalarIndex(
                          localFrequency, localY, x,
                          polarization, 0)] += magnitude * cosine;
                        block.amplitude[block.scalarIndex(
                          localFrequency, localY, x,
                          polarization, 1)] += magnitude * sine;
                      }
                  }
              }
          }
      }
    return block;
  }

  std::vector<long double> calculateWindowedRadiation(
      const std::vector<ParticleTrajectory>& trajectories,
      const RadiationConfig& config,
      const std::vector<long double>& omega,
      const std::vector<ObservationDirection>& directions,
      long double windowCenter,
      long double windowDuration,
      unsigned long long& skippedShortParticles)
  {
    if (!(windowDuration > 0.0L) || !std::isfinite(windowCenter) ||
        !std::isfinite(windowDuration))
      throw std::invalid_argument("Invalid radiation time window");
    std::vector<long double> amplitude(
      omega.size() * directions.size() * 2 * 2, 0.0L);
    skippedShortParticles = 0;
    const long double halfWindow = 0.5L * windowDuration;

    for (std::size_t particle = 0; particle < trajectories.size(); ++particle)
      {
        const ParticleTrajectory& trajectory = trajectories[particle];
        if (trajectory.records.size() < config.minimumRecordsPerParticle)
          {
            ++skippedShortParticles;
            continue;
          }
        std::vector<Vec3> beta(trajectory.records.size() - 1);
        for (std::size_t segment = 0; segment < beta.size(); ++segment)
          beta[segment] = segmentBeta(trajectory.records[segment],
                                     trajectory.records[segment + 1]);
        const long double charge = trajectory.charge;
        for (std::size_t directionIndex = 0;
             directionIndex < directions.size(); ++directionIndex)
          {
            const ObservationDirection& observation =
              directions[directionIndex];
            for (std::size_t knot = 1;
                 knot + 1 < trajectory.records.size(); ++knot)
              {
                const TrajectoryRecord& record = trajectory.records[knot];
                const Vec3 position = recordPosition(record);
                const long double reducedObserverTime =
                  static_cast<long double>(record.time) -
                  dot(observation.direction, position) / constants::c;
                const long double normalizedOffset =
                  (reducedObserverTime - windowCenter) / halfWindow;
                if (!(std::abs(normalizedOffset) < 1.0L)) continue;
                const long double windowWeight = 0.5L * (1.0L +
                  std::cos(constants::pi * normalizedOffset));
                const Vec3 delta = transverseFactor(
                  observation.direction, beta[knot]) -
                  transverseFactor(observation.direction, beta[knot - 1]);
                const long double component[2] = {
                  dot(delta, observation.horizontal),
                  dot(delta, observation.vertical)
                };
                if (component[0] == 0.0L && component[1] == 0.0L)
                  continue;
                for (std::size_t frequency = 0;
                     frequency < omega.size(); ++frequency)
                  {
                    const long double phase = std::remainder(
                      omega[frequency] * reducedObserverTime,
                      2.0L * constants::pi);
                    const long double cosine = std::cos(phase);
                    const long double sine = std::sin(phase);
                    for (std::size_t polarization = 0;
                         polarization < 2; ++polarization)
                      {
                        const long double magnitude = charge *
                          component[polarization] * windowWeight;
                        const std::size_t index =
                          ((frequency * directions.size() + directionIndex) *
                            2 + polarization) * 2;
                        amplitude[index] += magnitude * cosine;
                        amplitude[index + 1] += magnitude * sine;
                      }
                  }
              }
          }
      }
    return amplitude;
  }
}
