#include "pending_particle_statistics.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using namespace aprl;
typedef PendingParticleStatisticsSnapshot Snapshot;
void require(bool value, const char* message) {
  if (!value) throw std::runtime_error(message);
}
struct Value { long double w, k, total, x[4]; };
Value value(const RelativisticParticleSI& p, const BoostFrameTransform& f, double h) {
  Value v;
  const double gb = BoostFrameTransform::gammaFromProperVelocity(p.properVelocity);
  FieldVector<double> ul(0.); f.properVelocityBoxToLab(p.properVelocity, ul);
  v.w = static_cast<long double>(p.mass) / SI::electronMass;
  const long double rest = static_cast<long double>(p.mass) * SI::c * SI::c;
  v.k = rest * static_cast<long double>(p.properVelocity.norm2()) / static_cast<long double>(gb + 1.0);
  v.total = rest * gb;
  v.x[0] = BoostFrameTransform::gammaFromProperVelocity(ul);
  v.x[1] = gb; v.x[2] = p.position[2]; v.x[3] = h * p.frozenBetaZ;
  return v;
}
// Independent two-pass calculation over the whole remaining cohort.
Snapshot direct(const PendingParticleQueue& q, const BoostFrameTransform& f,
                double h, std::size_t first) {
  Snapshot result;
  if (first == q.initialCount()) return result;
  const Value anchor = value(q.orderedInitialParticle(first), f, h);
  long double displacement[4] = {};
  for (std::size_t i = first; i < q.initialCount(); ++i) {
    const Value v = value(q.orderedInitialParticle(i), f, h);
    result.representedWeight += v.w;
    result.kineticEnergy += v.k; result.totalEnergy += v.total; ++result.count;
    for (unsigned a = 0; a < 4; ++a) displacement[a] += v.w * (v.x[a] - anchor.x[a]);
    result.minimumGammaLab = std::min(result.minimumGammaLab, static_cast<double>(v.x[0]));
    result.maximumGammaLab = std::max(result.maximumGammaLab, static_cast<double>(v.x[0]));
  }
  for (unsigned a = 0; a < 4; ++a) {
    displacement[a] /= result.representedWeight;
    result.mean[a] = anchor.x[a] + displacement[a];
  }
  for (std::size_t i = first; i < q.initialCount(); ++i) {
    const Value v = value(q.orderedInitialParticle(i), f, h);
    for (unsigned a = 0; a < 4; ++a)
      for (unsigned b = 0; b < 4; ++b)
        result.centered[a][b] += v.w * ((v.x[a] - anchor.x[a]) - displacement[a]) *
          ((v.x[b] - anchor.x[b]) - displacement[b]);
  }
  return result;
}
long double maxMomentRelativeError = 0., maximumAxialRoundingDifference = 0.;
unsigned comparisons = 0;
void compare(const Snapshot& got, const Snapshot& expected) {
  require(got.count == expected.count, "Pending count differs");
  require(got.minimumGammaLab == expected.minimumGammaLab &&
          got.maximumGammaLab == expected.maximumGammaLab, "Pending gamma extrema differ");
  const long double eps = std::numeric_limits<long double>::epsilon();
  const long double a[] = {got.representedWeight, got.kineticEnergy, got.totalEnergy};
  const long double b[] = {expected.representedWeight, expected.kineticEnergy, expected.totalEnergy};
  for (unsigned k = 0; k < 3; ++k)
    require(std::abs(a[k] - b[k]) <= 128 * eps * std::abs(b[k]), "Pending energy/weight differs");
  for (unsigned i = 0; i < 4; ++i) {
    const long double scale = std::max(std::abs(expected.mean[i]), expected.count ?
      std::sqrt(std::abs(expected.centered[i][i]) / expected.representedWeight) : 0.L);
    require(std::abs(got.mean[i] - expected.mean[i]) <= 128 * eps * scale,
            "Pending mean differs");
    for (unsigned j = 0; j < 4; ++j) {
      const long double covarianceScale = std::sqrt(
        std::abs(expected.centered[i][i] * expected.centered[j][j]));
      const long double difference = std::abs(got.centered[i][j] - expected.centered[i][j]);
      require(difference <= 5e-15L * covarianceScale, "Pending centered moment differs");
      if (covarianceScale > 0.) maxMomentRelativeError = std::max(maxMomentRelativeError,
        difference / covarianceScale);
    }
  }
  ++comparisons;
}

void cohort(unsigned mode) {
  const double dz = 1.2543262415359034e-6, dt = dz / SI::c, h = SI::c * dt;
  PendingWindowGeometry g; g.dz = dz; g.heldDistancePerBeta = h;
  g.maximumStep = 50000;
  BoostFrameTransform frame; frame.setOriginsFromGamma(929.1305492858543, SI::c,
    -2.0014583498924204e-13, -.05, .0229021590039289);
  std::vector<RelativisticParticleSI> particles;
  for (std::size_t i = 0; i < 1555; ++i) {
    RelativisticParticleSI p;
    p.id = p.sourceId = i + 1;
    p.frozen = true;
    p.weight = mode == 2 ? 1. : 1. + static_cast<double>(i % 7) / 8.;
    p.mass = SI::electronMass * p.weight; p.charge = -SI::elementaryCharge * p.weight;
    const double gamma = mode == 1 ? 1170. : (mode == 2 ?
      1170. + 1e-8 * static_cast<double>(static_cast<int>(i % 17) - 8) :
      1164. + .012 * static_cast<double>(i % 997));
    FieldVector<double> lab(0.); lab[2] = std::sqrt((gamma - 1.) * (gamma + 1.));
    frame.properVelocityLabToBox(lab, p.properVelocity);
    p.frozenBetaZ = frame.heldLongitudinalBetaZ(p.properVelocity);
    p.position[2] = -.001 - .019 * static_cast<double>((i * 89) % 1555) / 1555.;
    particles.push_back(p);
  }
  PendingParticleQueue queue(g); queue.extract(particles);
  require(particles.empty() && queue.initialCount() == 1555, "Fixture unexpectedly admitted particles");
  PendingParticleStatistics stats(queue, frame, h);
  const std::size_t cursors[] = {0, 1, 127, 511, 512, 513, 1023, 1024, 1554, 1555};
  for (std::size_t cursor : cursors) {
    const Snapshot got = stats.snapshot(cursor);
    compare(got, direct(queue, frame, h, cursor));
    if (cursor == queue.initialCount()) continue;
    // Verify the affine diagnostic moment formula against a second complete
    // cohort scan in long double. This does not assert bitwise binary64 z.
    const std::size_t n = 30000;
    const long double meanZ = got.mean[2] + n * got.mean[3];
    const long double centeredZ = got.centered[2][2] + 2.L * n * got.centered[2][3] +
      static_cast<long double>(n) * n * got.centered[3][3];
    long double directMeanZ = 0., directCenteredZ = 0.;
    for (std::size_t i = cursor; i < queue.initialCount(); ++i) {
      const Value v = value(queue.orderedInitialParticle(i), frame, h);
      directMeanZ += v.w * (v.x[2] + n * v.x[3]);
    }
    directMeanZ /= got.representedWeight;
    if (got.count == 1) {
      const Value v = value(queue.orderedInitialParticle(cursor), frame, h);
      directMeanZ = v.x[2] + n * v.x[3];
    }
    for (std::size_t i = cursor; i < queue.initialCount(); ++i) {
      const Value v = value(queue.orderedInitialParticle(i), frame, h);
      const long double delta = v.x[2] + n * v.x[3] - directMeanZ;
      directCenteredZ += v.w * delta * delta;
    }
    require(std::abs(meanZ - directMeanZ) <= 2e-17L *
      std::max(std::abs(directMeanZ), std::sqrt(directCenteredZ / got.representedWeight)),
      "Affine pending mean formula differs");
    require(std::abs(centeredZ - directCenteredZ) <= 5e-15L * directCenteredZ,
      "Affine pending variance formula differs");
  }
  // Explicitly measure the diagnostic's finite rounding departure from the
  // state contract. It must not be fed back as the admission coordinate.
  for (std::size_t i = 0; i < queue.initialCount(); i += 97) {
    const auto& p = queue.orderedInitialParticle(i);
    const double increment = h * p.frozenBetaZ;
    double advanced = p.position[2];
    for (std::size_t n = 0; n < 30000; ++n) advanced += increment;
    const long double analytic = static_cast<long double>(p.position[2]) + 30000.L * increment;
    const long double difference = std::abs(static_cast<long double>(advanced) - analytic);
    maximumAxialRoundingDifference = std::max(maximumAxialRoundingDifference, difference);
    const long double roundingBound = 30000.L * std::numeric_limits<double>::epsilon() *
      (std::abs(p.position[2]) + 30000.L * std::abs(increment));
    require(difference <= roundingBound, "Axial rounding exceeded forward-error bound");
  }
  bool rejected = false;
  try { stats.snapshot(queue.initialCount() + 1); }
  catch (const std::out_of_range&) { rejected = true; }
  require(rejected, "Invalid statistics cursor was not rejected");
}
}

int main() {
  try {
    cohort(0); cohort(1); cohort(2);
    PendingWindowGeometry g; g.dz = 1.; g.heldDistancePerBeta = 1.; g.maximumStep = 1;
    std::vector<RelativisticParticleSI> none; PendingParticleQueue queue(g); queue.extract(none);
    BoostFrameTransform frame; PendingParticleStatistics empty(queue, frame, 1.);
    compare(empty.snapshot(), Snapshot());
    std::cout.precision(17);
    std::cout << "PASS pending_statistics comparisons=" << comparisons
      << " maximum_centered_relative_error=" << static_cast<double>(maxMomentRelativeError)
      << " maximum_diagnostic_axial_rounding_difference_m=" << static_cast<double>(maximumAxialRoundingDifference)
      << " head_scan_bound=" << PendingParticleStatistics::BlockSize - 1 << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FAIL pending_statistics: " << error.what() << '\n'; return 1;
  }
}
