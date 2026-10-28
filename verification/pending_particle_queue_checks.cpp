#include "pending_particle_queue.h"

#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace {
std::uint64_t bits(double value) {
  std::uint64_t result;
  std::memcpy(&result, &value, sizeof result);
  return result;
}
void require(bool ok, const char* message) {
  if (!ok) throw std::runtime_error(message);
}
double reference(double x, double delta, std::size_t steps) {
  // Keep the historical separately rounded recurrence visible to the compiler.
  volatile double position = x;
  for (std::size_t i = 0; i < steps; ++i) position = position + delta;
  return position;
}
void compare(double x, double delta, std::size_t steps) {
  const double expected = reference(x, delta, steps);
  if (!std::isfinite(expected)) return;
  const double actual = aprl::repeatHeldAddition(x, delta, steps);
  if (bits(expected) != bits(actual)) {
    std::cerr << std::hexfloat << "x=" << x << " delta=" << delta
              << " expected=" << expected << " actual=" << actual
              << std::defaultfloat << " steps=" << steps << '\n';
    throw std::runtime_error("Held jump differs from repeated addition");
  }
}

void additionChecks() {
  const double tiny = std::numeric_limits<double>::denorm_min();
  for (double zero : {0., -0.}) for (double delta : {0., -0., tiny, -tiny})
    for (std::size_t count : {0u, 1u, 2u, 3u, 33u}) compare(zero, delta, count);
  for (double x : {1., -1., 2., -2., tiny, -tiny,
                  std::numeric_limits<double>::min(), -std::numeric_limits<double>::min(),
                  0.02289965035, -0.02289965035})
    for (double delta : {0., tiny, -tiny, 2.9e-7, -2.9e-7})
      for (std::size_t count : {0u, 1u, 2u, 17u, 36721u}) compare(x, delta, count);
  // Entering a positive binade from below can produce an odd new mantissa,
  // followed by a different first tie-to-even stride. This catches the subtle
  // failure of multiplying the first observed ULP stride by the whole count.
  for (int exponent = -1000; exponent < 1000; exponent += 23) {
    const double edge = std::ldexp(1., exponent);
    const double ulp = std::ldexp(1., exponent - 52);
    for (double sign : {-1., 1.})
      for (double location : {std::nextafter(edge, 0.), edge,
                              std::nextafter(edge, std::numeric_limits<double>::infinity())})
        for (double scale : {.5, 1.5, 2.5, 3.5, 17.5})
          for (std::size_t count : {1u, 2u, 3u, 257u})
            compare(sign * location, sign * scale * ulp, count);
    for (double scale : {1.1, 1.3, 1.4, 1.6, 2.3, 3.4})
      for (std::size_t count : {8u, 9u, 10u, 11u, 17u})
        compare(-edge - 10 * ulp, scale * ulp, count);
  }
  // Exact cancellation, zero crossing, subnormal/normal and changing-binade
  // spacing; both signs are used because production z can straddle zero.
  for (double delta : {tiny, 3 * tiny, std::ldexp(1., -1022), .125, 2.9e-7})
    for (double scale : {-100., -1.5, -1., -.5, .5, 1., 100.})
      for (std::size_t count : {1u, 2u, 7u, 1001u}) compare(scale * delta, delta, count);
  std::mt19937_64 random(0x5f013214ULL);
  for (unsigned i = 0; i < 3000; ++i) {
    const double x = (static_cast<double>(random() % 2000001) - 1000000.) * 2.5e-8;
    const double delta = (static_cast<double>(random() % 1000001) - 500000.) * 1e-12;
    compare(x, delta, static_cast<std::size_t>(random() % 4097));
  }
  // A huge count must be handled by ULP/binade jumps, not count additions.
  const std::size_t huge = static_cast<std::size_t>(1000000000ULL);
  require(aprl::repeatHeldAddition(1., std::ldexp(1., -52), huge) ==
          1. + static_cast<double>(huge) * std::ldexp(1., -52), "Huge exact jump failed");
  require(aprl::repeatHeldAddition(1., std::ldexp(1., -54), huge) == 1., "Stagnation failed");
  if (std::numeric_limits<std::size_t>::digits > 52)
    require(bits(aprl::repeatHeldAddition(-std::numeric_limits<double>::min(), tiny,
              static_cast<std::size_t>(UINT64_C(1) << 52))) == bits(0.),
            "Huge subnormal cancellation failed");
  const int oldRounding = std::fegetround();
  std::fesetround(FE_UPWARD);
  bool rejected = false;
  try { aprl::repeatHeldAddition(1., .125, 4); } catch (const std::invalid_argument&) { rejected = true; }
  std::fesetround(oldRounding);
  require(rejected, "Unsupported rounding mode was accepted");
}

void queueChecks(std::size_t initialStep) {
  aprl::PendingWindowGeometry geometry;
  geometry.referenceOriginZ = -.25;
  geometry.initialOffset = 23;
  geometry.dz = .01;
  geometry.heldDistancePerBeta = .01;
  geometry.maximumStep = initialStep + 400;
  aprl::PendingParticleQueue queue(geometry, initialStep);
  std::vector<aprl::RelativisticParticleSI> active(513);
  const double inlet = queue.origin(initialStep);
  for (std::size_t i = 0; i < active.size(); ++i) {
    auto& p = active[i];
    p.id = i + 1; p.sourceId = i + 123;
    p.position[0] = .001; p.position[1] = -.002;
    p.position[2] = inlet + .013 - static_cast<double>(i) * .014;
    p.properVelocity[0] = .02; p.properVelocity[1] = -.03; p.properVelocity[2] = .25;
    p.charge = -1.6e-19; p.mass = 9e-31; p.weight = 2.;
    p.frozen = true; p.frozenBetaZ = .23 + static_cast<double>(i % 17) * .0001;
  }
  active[4].frozen = false; // Never optimize an active source.
  active[5].frozenBetaZ = -.1; // Unsupported direction stays in old path.
  active[6].frozenBetaZ = 0.;
  active[7].position[2] = inlet - geometry.heldDistancePerBeta * active[7].frozenBetaZ;
  std::vector<aprl::RelativisticParticleSI> plain = active;
  queue.extract(active);
  std::vector<bool> pending(plain.size(), true);
  for (const auto& p : active) pending[p.id - 1] = false;
  require(queue.size() + active.size() == plain.size(), "Partition lost particles");
  require(queue.initialParticles().size() == queue.initialCount(), "Initial statistics view differs");
  require(queue.initialAdmissionOrder().size() == queue.initialCount(), "Initial order differs");
  for (std::size_t i = 0; i < queue.initialCount(); ++i)
    require(bits(queue.orderedInitialParticle(i).position[2]) ==
            bits(plain[queue.orderedInitialParticle(i).id - 1].position[2]), "Initial state changed");
  std::size_t observerCount = 0;
  for (std::size_t step = initialStep; step <= geometry.maximumStep; ++step) {
    const double origin = std::fma(static_cast<double>(geometry.initialOffset -
      static_cast<std::ptrdiff_t>(step)), geometry.dz, geometry.referenceOriginZ);
    require(bits(origin) == bits(queue.origin(step)), "Origin FMA changed");
    std::vector<std::size_t> expected;
    for (std::size_t i = 0; i < plain.size(); ++i)
      if (pending[i] && plain[i].position[2] + geometry.heldDistancePerBeta * plain[i].frozenBetaZ >= origin)
        expected.push_back(i + 1);
    const std::size_t before = active.size();
    const std::size_t admitted = queue.admit(step, active,
      [&](const aprl::RelativisticParticleSI& initial, const aprl::RelativisticParticleSI& advanced) {
        require(initial.id == advanced.id && initial.frozen && advanced.frozen, "Admission identity changed");
        require(bits(initial.properVelocity[2]) == bits(advanced.properVelocity[2]), "Admission momentum changed");
        ++observerCount;
      });
    require(admitted == expected.size(), "Wrong admission step/count");
    for (std::size_t i = 0; i < admitted; ++i) {
      const auto& p = active[before + i];
      require(p.id == expected[i], "Admission order is not stable within due bucket");
      require(bits(p.position[2]) == bits(plain[p.id - 1].position[2]), "Admitted position differs bitwise");
      pending[p.id - 1] = false;
    }
    std::size_t pendingCount = 0;
    for (std::size_t i = 0; i < plain.size(); ++i) {
      if (pending[i]) ++pendingCount;
      plain[i].position[2] += geometry.heldDistancePerBeta * plain[i].frozenBetaZ;
    }
    require(queue.size() == pendingCount && queue.admittedCount() == observerCount,
            "Queue counts differ");
  }
  require(queue.size() > 0 && queue.nextDueStep() == geometry.maximumStep + 1,
          "Beyond-horizon pending particles were dropped");
  require(queue.memoryBytes() >= queue.initialCount() * sizeof(aprl::RelativisticParticleSI),
          "Queue memory accounting is incomplete");
  bool rejected = false;
  try { queue.admit(initialStep, active); } catch (const std::logic_error&) { rejected = true; }
  require(rejected, "Backward admission accepted");
}
}

int main() {
  try {
    additionChecks();
    queueChecks(0);
    queueChecks(137);
    std::cout << "PASS pending exact held jumps, event scheduling, immutable source states and counts\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FAIL " << error.what() << '\n';
    return 1;
  }
}
