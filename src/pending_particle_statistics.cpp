#include "pending_particle_statistics.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace aprl {
namespace {
struct CompensatedSum {
  long double sum = 0., correction = 0.;
  void add(long double value) {
    const long double next = sum + value;
    correction += std::abs(sum) >= std::abs(value) ?
      (sum - next) + value : (value - next) + sum;
    sum = next;
  }
  long double value() const { return sum + correction; }
};
struct Sample {
  long double weight, kinetic, total, value[4];
  long double kineticLab, cMomentumBox, cMomentumLab;
};
Sample sample(const RelativisticParticleSI& particle,
              const BoostFrameTransform& frame, Double heldDistance) {
  if (!particle.frozen || !(particle.mass > 0.) ||
      !std::isfinite(particle.mass) || !std::isfinite(particle.position[2]) ||
      !std::isfinite(particle.frozenBetaZ))
    throw std::invalid_argument("Invalid frozen particle in pending statistics");
  // Keep the ledger's expression order, including its double norm2() and
  // gamma computations before conversion to long double.
  const Double gammaBox = BoostFrameTransform::gammaFromProperVelocity(
    particle.properVelocity);
  FieldVector<Double> properVelocityLab(0.);
  frame.properVelocityBoxToLab(particle.properVelocity, properVelocityLab);
  const Double gammaLab = BoostFrameTransform::gammaFromProperVelocity(
    properVelocityLab);
  const long double mass = particle.mass, c = SI::c;
  const long double rest = mass * c * c;
  const long double u2 = particle.properVelocity.norm2();
  Sample result;
  result.weight = mass / static_cast<long double>(SI::electronMass);
  result.kinetic = rest * u2 / static_cast<long double>(gammaBox + 1.0);
  result.total = rest * static_cast<long double>(gammaBox);
  result.value[0] = gammaLab;
  result.value[1] = gammaBox;
  result.value[2] = particle.position[2];
  result.value[3] = heldDistance * particle.frozenBetaZ;
  // Keep the old energy/moment expressions above unchanged. Additional lab
  // kinetic energy uses a stable long-double norm, matching the new event
  // audit, without another frame transform or a scan at every field step.
  long double labNorm2 = 0.;
  for (unsigned axis = 0; axis < 3; ++axis)
    labNorm2 += static_cast<long double>(properVelocityLab[axis]) * properVelocityLab[axis];
  result.kineticLab = rest * labNorm2 / (std::sqrt(1.L + labNorm2) + 1.L);
  result.cMomentumBox = rest * particle.properVelocity[2];
  result.cMomentumLab = rest * properVelocityLab[2];
  if (!(result.weight > 0.) || !std::isfinite(result.weight) ||
      !std::isfinite(result.kinetic) || !std::isfinite(result.total) ||
      !std::isfinite(result.kineticLab) || !std::isfinite(result.cMomentumBox) ||
      !std::isfinite(result.cMomentumLab))
    throw std::invalid_argument("Nonfinite pending particle energy");
  for (unsigned axis = 0; axis < 4; ++axis)
    if (!std::isfinite(result.value[axis]))
      throw std::invalid_argument("Nonfinite pending particle moment");
  return result;
}
}

PendingParticleStatisticsSnapshot::PendingParticleStatisticsSnapshot()
  : count(0), representedWeight(0.), kineticEnergy(0.), totalEnergy(0.),
    kineticEnergyLab(0.), cMomentumZBox(0.), cMomentumZLab(0.),
    minimumGammaLab(std::numeric_limits<double>::infinity()),
    maximumGammaLab(-std::numeric_limits<double>::infinity()) {
  std::fill(mean, mean + 4, 0.L);
  for (unsigned i = 0; i < 4; ++i)
    std::fill(centered[i], centered[i] + 4, 0.L);
}

PendingParticleStatistics::PendingParticleStatistics(
    const PendingParticleQueue& queue, const BoostFrameTransform& frame,
    Double heldDistancePerBeta)
  : queue_(queue), frame_(frame), heldDistancePerBeta_(heldDistancePerBeta) {
  if (!(heldDistancePerBeta_ > 0.) || !std::isfinite(heldDistancePerBeta_))
    throw std::invalid_argument("Invalid held displacement for pending statistics");
  std::fill(reference_, reference_ + 4, 0.L);
  const std::size_t count = queue_.initialCount();
  if (count) {
    const Sample first = sample(queue_.orderedInitialParticle(0), frame_, heldDistancePerBeta_);
    std::copy(first.value, first.value + 4, reference_);
  }
  const std::size_t blocks = count / BlockSize + (count % BlockSize != 0);
  suffix_.resize(blocks + 1);
  CompensatedSum kinetic, total, kineticLab, cMomentumBox, cMomentumLab;
  for (std::size_t b = blocks; b > 0; --b) {
    const std::size_t begin = (b - 1) * BlockSize;
    const Snapshot part = block(begin, std::min(count, begin + BlockSize));
    suffix_[b - 1] = combine(part, suffix_[b]);
    kinetic.add(part.kineticEnergy);
    total.add(part.totalEnergy);
    suffix_[b - 1].kineticEnergy = kinetic.value();
    suffix_[b - 1].totalEnergy = total.value();
    kineticLab.add(part.kineticEnergyLab);
    cMomentumBox.add(part.cMomentumZBox);
    cMomentumLab.add(part.cMomentumZLab);
    suffix_[b - 1].kineticEnergyLab = kineticLab.value();
    suffix_[b - 1].cMomentumZBox = cMomentumBox.value();
    suffix_[b - 1].cMomentumZLab = cMomentumLab.value();
  }
}

PendingParticleStatistics::Snapshot PendingParticleStatistics::block(
    std::size_t begin, std::size_t end) const {
  Snapshot result;
  CompensatedSum kinetic, total, kineticLab, cMomentumBox, cMomentumLab;
  for (std::size_t k = begin; k < end; ++k) {
    const Sample item = sample(queue_.orderedInitialParticle(k), frame_, heldDistancePerBeta_);
    long double delta[4], relative[4];
    const long double newWeight = result.representedWeight + item.weight;
    for (unsigned axis = 0; axis < 4; ++axis) {
      relative[axis] = item.value[axis] - reference_[axis];
      delta[axis] = relative[axis] - result.mean[axis];
      result.mean[axis] += item.weight / newWeight * delta[axis];
    }
    for (unsigned i = 0; i < 4; ++i)
      for (unsigned j = i; j < 4; ++j) {
        result.centered[i][j] += item.weight * delta[i] *
          (relative[j] - result.mean[j]);
        result.centered[j][i] = result.centered[i][j];
      }
    result.representedWeight = newWeight;
    kinetic.add(item.kinetic);
    total.add(item.total);
    kineticLab.add(item.kineticLab);
    cMomentumBox.add(item.cMomentumBox);
    cMomentumLab.add(item.cMomentumLab);
    result.minimumGammaLab = std::min(result.minimumGammaLab,
      static_cast<double>(item.value[0]));
    result.maximumGammaLab = std::max(result.maximumGammaLab,
      static_cast<double>(item.value[0]));
    ++result.count;
  }
  result.kineticEnergy = kinetic.value();
  result.totalEnergy = total.value();
  result.kineticEnergyLab = kineticLab.value();
  result.cMomentumZBox = cMomentumBox.value();
  result.cMomentumZLab = cMomentumLab.value();
  return result;
}

PendingParticleStatistics::Snapshot PendingParticleStatistics::combine(
    const Snapshot& first, const Snapshot& second) {
  if (!first.count) return second;
  if (!second.count) return first;
  Snapshot result;
  result.count = first.count + second.count;
  result.representedWeight = first.representedWeight + second.representedWeight;
  result.kineticEnergy = first.kineticEnergy + second.kineticEnergy;
  result.totalEnergy = first.totalEnergy + second.totalEnergy;
  result.kineticEnergyLab = first.kineticEnergyLab + second.kineticEnergyLab;
  result.cMomentumZBox = first.cMomentumZBox + second.cMomentumZBox;
  result.cMomentumZLab = first.cMomentumZLab + second.cMomentumZLab;
  const long double fraction = second.representedWeight / result.representedWeight;
  const long double crossWeight = first.representedWeight * fraction;
  long double delta[4];
  for (unsigned i = 0; i < 4; ++i) {
    delta[i] = second.mean[i] - first.mean[i];
    result.mean[i] = first.mean[i] + fraction * delta[i];
  }
  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = i; j < 4; ++j) {
      result.centered[i][j] = first.centered[i][j] + second.centered[i][j] +
        crossWeight * delta[i] * delta[j];
      result.centered[j][i] = result.centered[i][j];
    }
  result.minimumGammaLab = std::min(first.minimumGammaLab, second.minimumGammaLab);
  result.maximumGammaLab = std::max(first.maximumGammaLab, second.maximumGammaLab);
  return result;
}

PendingParticleStatistics::Snapshot PendingParticleStatistics::snapshot(
    std::size_t admittedCount) const {
  if (admittedCount > queue_.initialCount())
    throw std::out_of_range("Pending statistics cursor exceeds input");
  if (admittedCount == queue_.initialCount()) return Snapshot();
  const std::size_t b = admittedCount / BlockSize;
  Snapshot result;
  if (admittedCount % BlockSize == 0) result = suffix_[b];
  else result = combine(block(admittedCount,
    std::min(queue_.initialCount(), (b + 1) * BlockSize)), suffix_[b + 1]);
  for (unsigned i = 0; i < 4; ++i) result.mean[i] += reference_[i];
  return result;
}

PendingParticleStatistics::Snapshot PendingParticleStatistics::snapshot() const {
  return snapshot(queue_.admittedCount());
}
}
