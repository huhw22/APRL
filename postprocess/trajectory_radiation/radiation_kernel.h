#ifndef TRAJECTORY_RADIATION_KERNEL_H
#define TRAJECTORY_RADIATION_KERNEL_H

#include <vector>

#include "radiation_types.h"

namespace radiation
{
  std::vector<ObservationDirection> makeObservationDirections(
      const RadiationConfig& config,
      std::size_t thetaYOffset,
      std::size_t thetaYCount);

  std::vector<long double> photonEnergyToOmega(
      const std::vector<double>& photonEnergyEV);

  RadiationBlock calculateRadiationBlock(
      const std::vector<ParticleTrajectory>& trajectories,
      const RadiationConfig& config,
      const std::vector<long double>& omega,
      std::size_t frequencyOffset,
      std::size_t frequencyCount,
      std::size_t thetaYOffset,
      std::size_t thetaYCount,
      unsigned long long& skippedShortParticles);
}

#endif
