#ifndef APRL_PENDING_PARTICLE_QUEUE_H
#define APRL_PENDING_PARTICLE_QUEUE_H

#include <cstddef>
#include <functional>
#include <vector>
#include "eb_particles.h"

namespace aprl {

// Exactly reproduce count separately rounded binary64 additions in round-to-
// nearest/even mode. Skips constant-spacing ranges by integer ULP arithmetic;
// no analytic z0 + count * increment approximation is used.
double repeatHeldAddition(double initial, double increment, std::size_t count);

struct PendingWindowGeometry {
  double referenceOriginZ = 0.;
  std::ptrdiff_t initialOffset = 0;
  double dz = 0.;
  // Pass SI::c * dt, evaluated in the same order as prepareParticle.
  double heldDistancePerBeta = 0.;
  std::size_t maximumStep = 0;
};

class PendingParticleQueue {
public:
  typedef RelativisticParticleSI Particle;
  typedef std::function<void(const Particle&, const Particle&)> AdmissionObserver;

  explicit PendingParticleQueue(const PendingWindowGeometry& geometry,
                                std::size_t initialStep = 0);

  // One initialization pass. Retains all particles that can participate in the
  // current step; removes only frozen, forward-held, still-upstream particles.
  // Particle values and relative order in the retained vector are unchanged.
  void extract(std::vector<Particle>& active);

  // Invoke before prepareParticle at this full-step index. Only due entries
  // are read. Their z equals the old repeated-addition state at the START of
  // this step, so unchanged prepareParticle handles crossing and deposition.
  // The observer receives (initial_pending_state, admitted_current_state).
  std::size_t admit(std::size_t step, std::vector<Particle>& active,
                    const AdmissionObserver& observer = AdmissionObserver());
  // Read-only view of precisely the next admission states. Used to construct
  // a compatible inlet layer without admitting early or scanning the queue.
  void visitDue(std::size_t step,
                const std::function<void(const Particle&)>& visitor) const;

  std::size_t size() const;
  bool empty() const { return size() == 0; }
  // maximumStep + 1 means no entry is due during the configured run horizon.
  std::size_t nextDueStep() const;
  double origin(std::size_t step) const;

  // For one-time aggregate initialization ONLY, immediately after extract().
  // These are states at initialStep, not lazily updated current positions.
  const std::vector<Particle>& initialParticles() const;
  const std::vector<std::size_t>& initialAdmissionOrder() const;
  const Particle& orderedInitialParticle(std::size_t index) const;
  std::size_t admittedCount() const { return cursor_; }
  std::size_t memoryBytes() const;
  std::size_t initialStep() const { return initialStep_; }
  std::size_t initialCount() const { return pending_.size(); }

private:
  std::size_t firstDue(const Particle& particle) const;
  bool participates(const Particle& particle, std::size_t step) const;
  void advanceNextDue();
  PendingWindowGeometry geometry_;
  std::size_t initialStep_, lastStep_, cursor_, nextDue_;
  bool extracted_;
  std::vector<Particle> pending_;
  std::vector<std::size_t> order_, offsets_;
};
}
#endif
