#ifndef APRL_PENDING_PARTICLE_STATISTICS_H
#define APRL_PENDING_PARTICLE_STATISTICS_H

#include <cstddef>
#include <vector>
#include "boostframe.h"
#include "pending_particle_queue.h"

namespace aprl {

// All central moments are weighted SUMS, not variances. Coordinates are:
// 0 = gamma_lab, 1 = gamma_box, 2 = box z at queue.initialStep(),
// 3 = the separately rounded held increment (c * dt) * frozenBetaZ.
struct PendingParticleStatisticsSnapshot {
  std::size_t count;
  long double representedWeight, kineticEnergy, totalEnergy;
  long double mean[4], centered[4][4];
  double minimumGammaLab, maximumGammaLab;
  PendingParticleStatisticsSnapshot();
};

// The queue must outlive this object, retaining its initial ordered records.
// Construction visits each pending particle once. Each query visits at most
// BlockSize-1 records in the current partial block and combines one suffix.
// No subtraction of an almost exhausted full-bunch moment is required.
//
// Energy/gamma statistics use the exact existing ledger expressions. Axial
// statistics use z_initial + elapsed_steps * held_increment in long double.
// This is an analytic diagnostic trajectory, not bitwise repeated binary64
// addition. The queue separately reconstructs actual admission coordinates
// exactly; no result from this statistics object drives particle or field
// evolution. Tests and callers must distinguish those two contracts.
class PendingParticleStatistics {
public:
  enum { BlockSize = 512 };
  typedef PendingParticleStatisticsSnapshot Snapshot;

  PendingParticleStatistics(const PendingParticleQueue& queue,
                            const BoostFrameTransform& frame,
                            Double heldDistancePerBeta);
  Snapshot snapshot(std::size_t admittedCount) const;
  Snapshot snapshot() const;
  std::size_t initialStep() const { return queue_.initialStep(); }
  std::size_t storageBytes() const { return suffix_.capacity() * sizeof(Snapshot); }

  static Snapshot combine(const Snapshot& first, const Snapshot& second);

private:
  Snapshot block(std::size_t begin, std::size_t end) const;
  const PendingParticleQueue& queue_;
  BoostFrameTransform frame_;
  Double heldDistancePerBeta_;
  long double reference_[4];
  std::vector<Snapshot> suffix_;
};
}
#endif
