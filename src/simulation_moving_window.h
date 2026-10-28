#ifndef APRL_SIMULATION_MOVING_WINDOW_H
#define APRL_SIMULATION_MOVING_WINDOW_H
#include <cstdint>
#include <fstream>
#include <memory>
#include "moving_window_fields.h"
#include "particle_field_halo.h"
#include "pending_particle_queue.h"
#include "pending_particle_statistics.h"
namespace aprl {
namespace verification { class BlockedCowanReference; }
// No full-domain E/B grid is retained. Particle storage deliberately remains
// the complete input, with upstream waiting particles scheduled separately.
struct SimulationMovingWindowState {
  EBGridGeometry referenceGeometry;
  FieldVector<Double> referenceOrigin;
  std::ptrdiff_t offset = 0, initialOffset = 0;
  Double solverBeta = 0., translationBeta = 0., exitLabZ = 0.;
  std::size_t guards[3] = {0,0,0};
#ifdef APRL_HAS_MOVING_WINDOW
  std::unique_ptr<verification::BlockedCowanReference> reference;
#endif
  std::unique_ptr<MovingWindowFields> fields;
  std::unique_ptr<EBParticleFieldHalo> particleHalo;
  std::unique_ptr<PendingParticleQueue> pending;
  std::unique_ptr<PendingParticleStatistics> pendingStatistics;
  std::uint64_t admittedParticles = 0;
  std::ofstream ledger;
  unsigned long long exits = 0, croppedCarriers = 0;
  Double exitCharge = 0., croppedCharge = 0.;
  long double insertedFieldEnergy = 0., croppedFieldEnergy = 0.;
  double shiftSeconds = 0., migrationSeconds = 0.;
  std::uint64_t movedBytes = 0, communicatedBytes = 0;
  ~SimulationMovingWindowState();
};
}
#endif
