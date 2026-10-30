#ifndef APRL_SIMULATION_MOVING_WINDOW_H
#define APRL_SIMULATION_MOVING_WINDOW_H
#include <cstdint>
#include <fstream>
#include <memory>
#include <vector>
#include "moving_window_fields.h"
#include "particle_field_halo.h"
#include "pending_particle_queue.h"
#include "pending_particle_statistics.h"
#include "transverse_reference_current.h"
#include "charge_matched_inlet.h"
namespace aprl {
namespace verification { class BlockedCowanReference; }
// No full-domain E/B grid is retained. Particle storage deliberately remains
// the complete input, with upstream waiting particles scheduled separately.
struct SimulationMovingWindowState {
  // Only constrained source segments actually deposited in the current local
  // physical aperture are retained for this field step; pending particles do
  // not create samples. This is a path-work estimate, not a CIC adjoint.
  struct ConstraintWorkSample {
    FieldVector<Double> midpoint;
    FieldVector<Double> chargeDisplacement;
    long double oldElectricDot = 0.;
    bool carrier = false;
  };
  EBGridGeometry referenceGeometry;
  FieldVector<Double> referenceOrigin;
  std::ptrdiff_t offset = 0, initialOffset = 0;
  Double solverBeta = 0., translationBeta = 0., exitLabZ = 0.;
  std::size_t guards[3] = {0,0,0};
#ifdef APRL_HAS_MOVING_WINDOW
  std::unique_ptr<verification::BlockedCowanReference> reference;
  std::unique_ptr<TransverseReferenceCurrent> referenceCurrent;
#endif
  std::unique_ptr<MovingWindowFields> fields;
  std::unique_ptr<EBParticleFieldHalo> particleHalo;
  std::unique_ptr<PendingParticleQueue> pending;
  std::unique_ptr<PendingParticleStatistics> pendingStatistics;
  std::unique_ptr<ChargeMatchedInlet> inletCompletion;
  std::ofstream inletCompletionLedger;
  std::uint64_t admittedParticles = 0;
  std::ofstream ledger;
  unsigned long long exits = 0, croppedCarriers = 0;
  Double exitCharge = 0., croppedCharge = 0.;
  long double insertedFieldEnergy = 0., croppedFieldEnergy = 0.;
  std::vector<ConstraintWorkSample> constraintWorkSamples;
  long double frozenConstraintWork = 0., carrierConstraintWork = 0.;
  long double frozenConstraintAbsWork = 0., carrierConstraintAbsWork = 0.;
  std::uint64_t frozenConstraintSegments = 0, carrierConstraintSegments = 0;
  std::size_t peakConstraintSamples = 0;
  double constraintDiagnosticSeconds = 0.;
  double shiftSeconds = 0., migrationSeconds = 0.;
  std::uint64_t movedBytes = 0, communicatedBytes = 0;
  ~SimulationMovingWindowState();
};
}
#endif
