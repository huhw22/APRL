#ifndef APRL_VERIFICATION_MOVING_WINDOW_REMAP_H
#define APRL_VERIFICATION_MOVING_WINDOW_REMAP_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "eb_cpml.h"

namespace aprl {
namespace verification {

enum class MovingWindowCurrentPolicy {
  Preserve, // Move J overlap and supply its incoming values, just like E/B.
  Reset     // Move E/B, then clear all J; do not call the provider for J.
};

// Allocation-free callback descriptor. Components are Ex,Ey,Ez,Bx,By,Bz,
// Jx,Jy,Jz (0..8); i,j,k refer to each component's native Yee lattice in the
// NEW window. The caller supplies the correct stagger and stored time level.
// An absent callback supplies zero. The callback must not mutate fields or
// CPML histories and should not throw: a throw can leave a partial remap.
struct MovingWindowIncoming {
  using Function = Double (*)(void* context, unsigned int component,
                             std::size_t i, std::size_t j, std::size_t k,
                             Double newOriginZ);
  void* context = nullptr;
  Function value = nullptr;
};

// Byte counts describe logical stored values, not measured memory-bus
// traffic: a memmove counts its payload once, not read+write. Discarded
// values need no separate memory operation. Times include only helper work;
// incomingWallSeconds includes the caller's callback and new-layer writes.
struct MovingWindowShiftReport {
  int signedCells = 0;
  Double newOriginZ = 0.0;
  std::uint64_t fieldMovedBytes = 0;
  std::uint64_t fieldIncomingBytes = 0;
  std::uint64_t fieldDiscardedBytes = 0;
  std::uint64_t currentResetBytes = 0;
  std::uint64_t historyMovedBytes = 0;
  std::uint64_t historyZeroedBytes = 0;
  std::uint64_t historyDiscardedBytes = 0;
  std::uint64_t fieldMoveBlocks = 0;
  std::uint64_t historyMoveBlocks = 0;
  std::uint64_t historyZeroBlocks = 0;
  std::uint64_t incomingValues = 0;
  double fieldWallSeconds = 0.0;
  double historyWallSeconds = 0.0;
  double incomingWallSeconds = 0.0;
  double totalWallSeconds = 0.0;
};

// Read-only storage description for debugging. This is NOT a second oracle:
// acceptance should independently derive the CPML indexing from its contract.
// In each compressed history, activeAxisIndices are stored in ascending order.
// The per-active-index cross section is (k,j) for x derivatives, (k,i) for y
// derivatives, and (j,i) for z derivatives, with the second index unit-stride.
struct MovingWindowHistoryLayout {
  unsigned int fieldComponent = 0; // Same 0..5 E/B numbering as callback.
  unsigned int derivativeAxis = 0;
  std::array<std::size_t, 3> dimensions = {{0, 0, 0}};
  std::vector<std::size_t> activeAxisIndices;
  std::size_t values = 0;
};

struct MovingWindowPlanSummary {
  int signedCells = 0;
  std::uint64_t fieldMovedBytes = 0;    // Preserve policy, all nine arrays.
  std::uint64_t fieldIncomingBytes = 0;
  std::uint64_t historyMovedBytes = 0;
  std::uint64_t historyZeroedBytes = 0;
  std::size_t historyMoveBlocks = 0;
  std::size_t historyZeroBlocks = 0;
};

// Verification-only, single-rank, whole-window integer-z state migration.
// The existing EBConvolutionalPML must have been built with this geometry,
// globalNz=geometry.nz, localZCellOffset=0, both physical z faces=true, and
// exactly these parameters. Its public interface exposes history sizes but
// not coefficients; the caller remains responsible for that construction
// contract. Runtime checks verify grid shapes, boundary type, and all sizes.
//
// New origin = old origin + signedCells*dz; new[k] = old[k+signedCells]
// wherever both arrays contain that world position. Overlap is moved without
// arithmetic or interpolation. CPML histories retain a world-position value
// only where BOTH old and new derivative profiles are active; all other new
// history slots are zero. No history is copied between derivative terms.
//
// Only constructor-approved shifts are accepted. The constructor precomputes
// all compressed-history operations; shift performs no allocation, no full-
// volume lookup, no history resize, and no CPML reconstruction. Normal field
// updates still use the unmodified production Cowan/Yee and CPML formulas.
// The caller owns origin bookkeeping (prefer an integer global cell offset)
// and supplies newOriginZ; the helper does not accumulate floating origins.
//
// IMPORTANT: moving a CPML profile changes its physical location and the
// coefficient history experienced by a point. This migration rule does not
// establish equivalence to a stationary physical absorbing boundary and does
// not supply exterior sources or prove moving-boundary stability.
class MovingWindowRemap {
public:
  MovingWindowRemap(const EBGridGeometry& geometry,
                    const EBCPMLParameters& parameters,
                    const std::vector<int>& supportedShifts = {-3, -1, 1, 3});
  ~MovingWindowRemap();
  MovingWindowRemap(const MovingWindowRemap&) = delete;
  MovingWindowRemap& operator=(const MovingWindowRemap&) = delete;

  // signedCells==0 is always a strict no-op, including with Reset policy.
  // Nonzero shifts must satisfy abs(signedCells)<=geometry.nz. The new origin
  // is metadata/callback input only, and must agree with the caller's offset.
  MovingWindowShiftReport shift(
      EBFieldGrid& fields, int signedCells, Double newOriginZ,
      MovingWindowIncoming incoming = MovingWindowIncoming(),
      MovingWindowCurrentPolicy currentPolicy = MovingWindowCurrentPolicy::Preserve);

  double setupWallSeconds() const;
  // Persistent helper object/plan/layout storage; excludes allocator overhead,
  // the existing fields/CPML, callback state, and temporary constructor storage.
  std::size_t memoryBytes() const;
  std::size_t preparedShiftCount() const;
  MovingWindowPlanSummary planSummary(int signedCells) const;
  const MovingWindowHistoryLayout& historyLayout(std::size_t index) const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace verification
} // namespace aprl

#endif
