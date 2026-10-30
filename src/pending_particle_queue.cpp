#include "pending_particle_queue.h"

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace aprl {
namespace {
const std::uint64_t signBit = UINT64_C(1) << 63;
const std::uint64_t exponentMask = UINT64_C(0x7ff0000000000000);
const std::uint64_t fractionMask = UINT64_C(0x000fffffffffffff);

std::uint64_t bits(double value) {
  std::uint64_t result;
  std::memcpy(&result, &value, sizeof result);
  return result;
}
double fromBits(std::uint64_t value) {
  double result;
  std::memcpy(&result, &value, sizeof result);
  return result;
}
std::uint64_t ordered(std::uint64_t value) {
  return (value & signBit) ? ~value : value ^ signBit;
}
double fromOrdered(std::uint64_t value) {
  return fromBits((value & signBit) ? value ^ signBit : ~value);
}
bool sameSpacingRange(std::uint64_t a, std::uint64_t b) {
  return (a & (signBit | exponentMask)) == (b & (signBit | exponentMask));
}

// A separately rounded addition is deliberate: contracting x+increment into
// a multiply-add in callers would not reproduce the old held update.
double add(double x, double increment) { return x + increment; }

double repeatPositive(double x, double increment, std::size_t count) {
  while (count) {
    const double y = add(x, increment);
    if (!std::isfinite(y)) throw std::overflow_error("Held position overflow");
    if (count == 1 || y == x) return y;
    const double z = add(y, increment);
    const std::uint64_t bx = bits(x), by = bits(y), bz = bits(z);
    if (std::isfinite(z) && sameSpacingRange(bx, by) && sameSpacingRange(by, bz)) {
      const std::uint64_t ox = ordered(bx), oy = ordered(by), oz = ordered(bz);
      const std::uint64_t stride = oy - ox;
      // Half-ULP ties can have a different first stride after entering a new
      // binade. Two equal strides establish its stable round-to-even phase.
      if (stride && oz - oy == stride) {
        const std::uint64_t lastBits = (bx & signBit) ?
          bx & (signBit | exponentMask) : bx | fractionMask;
        const std::uint64_t available = (ordered(lastBits) - ox) / stride;
        // Leave the boundary-reaching addition literal. At a negative power
        // of two the next binade has half-sized ULPs: rounding a sum just on
        // that finer side cannot be inferred using the old spacing alone.
        const std::uint64_t safeAvailable = available ? available - 1 : 0;
        const std::size_t jump = static_cast<std::size_t>(
          std::min<std::uint64_t>(safeAvailable, count));
        if (jump) {
          x = fromOrdered(ox + static_cast<std::uint64_t>(jump) * stride);
          count -= jump;
          continue;
        }
      }
    }
    x = y;
    --count;
  }
  return x;
}
}

double repeatHeldAddition(double initial, double increment, std::size_t count) {
  static_assert(sizeof(double) == sizeof(std::uint64_t), "Requires binary64");
  static_assert(std::numeric_limits<double>::is_iec559 &&
                std::numeric_limits<double>::digits == 53, "Requires IEEE binary64");
  if (!std::isfinite(initial) || !std::isfinite(increment))
    throw std::invalid_argument("Held position/increment must be finite");
  if (std::fegetround() != FE_TONEAREST)
    throw std::invalid_argument("Held jump requires round-to-nearest/even");
  if (!count) return initial;
  if (increment < 0.) {
    const double result = -repeatPositive(-initial, -increment, count);
    // Exact cancellation in round-to-nearest produces +0 for either direction.
    return result == 0. ? 0. : result;
  }
  return repeatPositive(initial, increment, count);
}

PendingParticleQueue::PendingParticleQueue(const PendingWindowGeometry& geometry,
    std::size_t initialStep) : geometry_(geometry), initialStep_(initialStep),
    lastStep_(initialStep), cursor_(0), nextDue_(initialStep), extracted_(false) {
  if (!(geometry_.dz > 0.) || !(geometry_.heldDistancePerBeta > 0.) ||
      !std::isfinite(geometry_.dz) || !std::isfinite(geometry_.heldDistancePerBeta) ||
      !std::isfinite(geometry_.referenceOriginZ) || initialStep > geometry_.maximumStep ||
      geometry_.maximumStep >= static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) ||
      geometry_.initialOffset < std::numeric_limits<std::ptrdiff_t>::min() +
        static_cast<std::ptrdiff_t>(geometry_.maximumStep))
    throw std::invalid_argument("Invalid pending-window geometry/horizon");
  // Reject unsupported arithmetic before extracting any source particle.
  repeatHeldAddition(0., 0., 0);
  origin(initialStep_);
  origin(geometry_.maximumStep);
}

double PendingParticleQueue::origin(std::size_t step) const {
  if (step > geometry_.maximumStep) throw std::out_of_range("Pending origin exceeds horizon");
  const std::ptrdiff_t offset = geometry_.initialOffset - static_cast<std::ptrdiff_t>(step);
  const double value = std::fma(static_cast<double>(offset), geometry_.dz, geometry_.referenceOriginZ);
  if (!std::isfinite(value)) throw std::overflow_error("Pending window origin overflow");
  return value;
}

bool PendingParticleQueue::participates(const Particle& particle, std::size_t step) const {
  const double increment = geometry_.heldDistancePerBeta * particle.frozenBetaZ;
  const double end = repeatHeldAddition(particle.position[2], increment, step - initialStep_ + 1);
  return end >= origin(step);
}

std::size_t PendingParticleQueue::firstDue(const Particle& particle) const {
  const double increment = geometry_.heldDistancePerBeta * particle.frozenBetaZ;
  const long double distance = static_cast<long double>(origin(initialStep_)) -
    static_cast<long double>(particle.position[2]) - static_cast<long double>(increment);
  const long double estimate = std::ceil(distance /
    (static_cast<long double>(geometry_.dz) + static_cast<long double>(increment)));
  const std::size_t sentinel = geometry_.maximumStep + 1;
  const long double stepsLeft = static_cast<long double>(geometry_.maximumStep - initialStep_);
  const std::size_t guess = estimate > stepsLeft ? sentinel :
    initialStep_ + static_cast<std::size_t>(std::max(1.L, estimate));
  // The analytic guess only narrows the search. Both sides of the admission
  // event are checked with exact old positions and the actual FMA grid origin.
  if (guess <= geometry_.maximumStep && participates(particle, guess) &&
      !participates(particle, guess - 1)) return guess;
  if (guess == sentinel && !participates(particle, geometry_.maximumStep)) return sentinel;
  std::size_t lower = initialStep_ + 1, upper = sentinel;
  while (lower < upper) {
    const std::size_t middle = lower + (upper - lower) / 2;
    if (participates(particle, middle)) upper = middle;
    else lower = middle + 1;
  }
  return lower;
}

void PendingParticleQueue::extract(std::vector<Particle>& active) {
  if (extracted_) throw std::logic_error("Pending queue can only be initialized once");
  const std::size_t buckets = geometry_.maximumStep - initialStep_ + 2;
  offsets_.assign(buckets + 1, 0);
  pending_.reserve(active.size());
  std::vector<std::size_t> due;
  due.reserve(active.size());
  const double inlet = origin(initialStep_);
  std::size_t retained = 0;
  for (std::size_t i = 0; i < active.size(); ++i) {
    const Particle& particle = active[i];
    if (particle.frozen && particle.position[2] < inlet &&
        particle.frozenBetaZ >= 0. && !participates(particle, initialStep_)) {
      const std::size_t event = firstDue(particle);
      due.push_back(event - initialStep_);
      ++offsets_[event - initialStep_ + 1];
      pending_.push_back(particle);
    } else {
      if (retained != i) active[retained] = particle;
      ++retained;
    }
  }
  active.resize(retained);
  for (std::size_t i = 1; i < offsets_.size(); ++i) offsets_[i] += offsets_[i - 1];
  order_.resize(pending_.size());
  std::vector<std::size_t> insertion(offsets_);
  for (std::size_t i = 0; i < due.size(); ++i) order_[insertion[due[i]]++] = i;
  extracted_ = true;
  advanceNextDue();
}

void PendingParticleQueue::advanceNextDue() {
  while (nextDue_ <= geometry_.maximumStep &&
         offsets_[nextDue_ - initialStep_ + 1] <= cursor_) ++nextDue_;
}

std::size_t PendingParticleQueue::admit(std::size_t step, std::vector<Particle>& active,
    const AdmissionObserver& observer) {
  if (!extracted_) throw std::logic_error("Pending queue has not been initialized");
  if (step < lastStep_ || step > geometry_.maximumStep || step > nextDue_)
    throw std::logic_error("Pending admission skipped an event or exceeds its horizon");
  lastStep_ = step;
  const std::size_t end = offsets_[step - initialStep_ + 1];
  const std::size_t admitted = end - cursor_;
  if (admitted > active.max_size() - active.size()) throw std::length_error("Pending admission overflow");
  // Let vector retain its geometric growth policy. Reserving the exact new
  // size at every admission would copy the complete active set every step.
  while (cursor_ < end) {
    const Particle& initial = pending_[order_[cursor_]];
    Particle particle = initial;
    particle.position[2] = repeatHeldAddition(initial.position[2],
      geometry_.heldDistancePerBeta * initial.frozenBetaZ, step - initialStep_);
    active.push_back(particle);
    if (observer) observer(initial, particle);
    ++cursor_;
  }
  advanceNextDue();
  return admitted;
}

std::size_t PendingParticleQueue::size() const { return pending_.size() - cursor_; }
void PendingParticleQueue::visitDue(std::size_t step,
    const std::function<void(const Particle&)>& visitor) const {
  if(!extracted_ || step<lastStep_ || step>geometry_.maximumStep || step>nextDue_)
    throw std::logic_error("Pending inlet view skipped an event or exceeds its horizon");
  const auto end=offsets_[step-initialStep_+1];
  for(std::size_t n=cursor_;n<end;++n) {
    const auto& initial=pending_[order_[n]];Particle particle=initial;
    particle.position[2]=repeatHeldAddition(initial.position[2],
      geometry_.heldDistancePerBeta*initial.frozenBetaZ,step-initialStep_);
    visitor(particle);
  }
}
std::size_t PendingParticleQueue::nextDueStep() const { return nextDue_; }
const std::vector<PendingParticleQueue::Particle>& PendingParticleQueue::initialParticles() const {
  if (!extracted_ || cursor_) throw std::logic_error("Initial pending statistics must precede admission");
  return pending_;
}
const std::vector<std::size_t>& PendingParticleQueue::initialAdmissionOrder() const {
  if (!extracted_ || cursor_) throw std::logic_error("Initial pending order must precede admission");
  return order_;
}
const PendingParticleQueue::Particle& PendingParticleQueue::orderedInitialParticle(std::size_t index) const {
  return pending_.at(order_.at(index));
}
std::size_t PendingParticleQueue::memoryBytes() const {
  return pending_.capacity() * sizeof(Particle) +
    (order_.capacity() + offsets_.capacity()) * sizeof(std::size_t);
}
}
