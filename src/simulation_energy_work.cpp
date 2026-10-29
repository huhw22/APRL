#include "simulation.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace aprl {
namespace {
// Include clipping, sampling and queue overhead, including empty clipped
// segments. This timer never controls the particle or field calculation.
struct ConstraintTimer {
  explicit ConstraintTimer(double& sum) : sum(sum), started(MPI_Wtime()) {}
  ~ConstraintTimer() { sum += MPI_Wtime() - started; }
  double& sum;
  double started;
};

long double displacementDot(const FieldVector<Double>& qdx,
                            const FieldVector<Double>& electric) {
  long double result = 0.;
  for (unsigned axis = 0; axis < 3; ++axis)
    result += static_cast<long double>(qdx[axis]) * electric[axis];
  return result;
}
}

void Simulation::queueMovingConstraintWork(const FieldVector<Double>& start,
    const FieldVector<Double>& end, Double charge, bool carrier) {
  if (!movingWindow_ || !config_.energyLedger.enabled) return;
  auto& w = *movingWindow_;
  ConstraintTimer timer(w.constraintDiagnosticSeconds);
  if (!std::isfinite(charge))
    throw std::runtime_error("Non-finite moving constraint-source charge");
  if (charge == 0.) return;

  // The depositing caller already clips to this rank's owned trajectory.
  // Exclude CPML volume from this estimate just as the field-energy ledger
  // does. No prescribed analytic-device field belongs to this source work.
  Double first = 0., last = 1.;
  for (unsigned axis = 0; axis < 3; ++axis) {
    const Double low = particleCPML_->physicalLower(axis);
    const Double high = particleCPML_->physicalUpper(axis);
    const Double delta = end[axis] - start[axis];
    if (!std::isfinite(start[axis]) || !std::isfinite(end[axis]))
      throw std::runtime_error("Non-finite moving constraint-source segment");
    if (delta == 0.) {
      if (start[axis] < low || start[axis] > high) return;
    } else {
      Double a = (low - start[axis]) / delta;
      Double b = (high - start[axis]) / delta;
      if (a > b) std::swap(a, b);
      first = std::max(first, a);
      last = std::min(last, b);
      if (!(last > first)) return;
    }
  }
  SimulationMovingWindowState::ConstraintWorkSample sample;
  sample.carrier = carrier;
  bool displaced = false;
  for (unsigned axis = 0; axis < 3; ++axis) {
    const Double advance = end[axis] - start[axis];
    const Double entry = start[axis] + first * advance;
    const Double exit = start[axis] + last * advance;
    sample.midpoint[axis] = 0.5 * (entry + exit);
    sample.chargeDisplacement[axis] = charge * (exit - entry);
    displaced = displaced || sample.chargeDisplacement[axis] != 0.;
  }
  if (!displaced) return;
  FieldVector<Double> electric(0.), magnetic(0.);
  // Called inside pushDepositAndMigrate, while the old-E interpolation halo
  // is installed. B is ignored: magnetic forces do no path work.
  fields_->sampleFieldsPosition(sample.midpoint, localOriginBox_, electric,
                               magnetic);
  sample.oldElectricDot = displacementDot(sample.chargeDisplacement, electric);
  w.constraintWorkSamples.push_back(sample);
  w.peakConstraintSamples = std::max(w.peakConstraintSamples,
                                    w.constraintWorkSamples.size());
}

void Simulation::finalizeMovingConstraintWork() {
  if (!movingWindow_ || !config_.energyLedger.enabled) return;
  auto& w = *movingWindow_;
  ConstraintTimer timer(w.constraintDiagnosticSeconds);
  // All ranks participate, including ranks with no samples. Reuse the small
  // neighbour-plane buffers; a local-only sample would clamp Ez at MPI faces.
  w.particleHalo->exchangeAndInstall(*fields_);
  try {
    for (const auto& sample : w.constraintWorkSamples) {
      FieldVector<Double> electric(0.), magnetic(0.);
      fields_->sampleFieldsPosition(sample.midpoint, localOriginBox_, electric,
                                   magnetic);
      const long double work = -0.5L * (sample.oldElectricDot +
        displacementDot(sample.chargeDisplacement, electric));
      if (!std::isfinite(work))
        throw std::runtime_error("Non-finite moving constraint-source work");
      if (sample.carrier) {
        w.carrierConstraintWork += work;
        w.carrierConstraintAbsWork += std::abs(work);
        ++w.carrierConstraintSegments;
      } else {
        w.frozenConstraintWork += work;
        w.frozenConstraintAbsWork += std::abs(work);
        ++w.frozenConstraintSegments;
      }
    }
  } catch (...) {
    w.particleHalo->clear(*fields_);
    throw;
  }
  w.particleHalo->clear(*fields_);
  w.constraintWorkSamples.clear();
}
}
