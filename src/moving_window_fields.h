#ifndef APRL_MOVING_WINDOW_FIELDS_H
#define APRL_MOVING_WINDOW_FIELDS_H

#include <cstddef>
#include <cstdint>
#include <array>
#include <memory>
#include <mpi.h>
#include "eb_cpml.h"

namespace aprl {

// Fixed-world native indices, not indices relative to a rank/window. The
// caller supplies E(n), B(n-1/2); this layer performs no time interpolation.
struct MovingWindowReference {
  using Function = Double (*)(void*, unsigned int, std::size_t, std::size_t,
                             std::ptrdiff_t, std::size_t);
  void* context = nullptr;
  Function value = nullptr;
};

struct MovingWindowFieldShiftReport {
  std::uint64_t fieldMovedBytes = 0;
  std::uint64_t historyMovedBytes = 0;
  std::uint64_t currentClearedBytes = 0;
  std::uint64_t sentBytes = 0;
  std::uint64_t receivedBytes = 0;
  double communicationSeconds = 0.;
  double totalSeconds = 0.;
};

// Rank-local cumulative diagnostics, in joules in the computational frame.
// Use MPI_SUM for energies/counts (remapSteps is a per-rank step count).
// All energies use EBFieldGrid::fieldEnergyDensityCell over the transverse
// physical aperture, excluding xy CPML. Positive overwrite/remap means field
// energy added. Inserted/cropped are positive stored energies, not fluxes.
// These are changes caused by the named boundary operations, not corrections
// fitted to a propagation residual and not an exact Maxwell energy invariant.
struct MovingWindowFieldEnergyLedger {
  bool enabled = false;
  long double inletOverwriteEnergyJ = 0.;
  long double insertedEnergyJ = 0.;
  long double croppedEnergyJ = 0.;
  long double remapEnergyJ = 0.;
  // After summing ranks: remap - (inserted - cropped). This explicitly keeps
  // any contribution from inconsistent duplicated Yee nodes at MPI interfaces.
  long double remapInterfaceCorrectionJ = 0.;
  // Outward Poynting flux through x-,x+,y-,y+,z-,z+, trapezoid-integrated
  // at the same fixed window origin: beginStep and post-CPML/post-MPI,
  // with the measured inlet overwrite contribution removed from the endpoint.
  // This is S dot n only; remap is the separate moving-volume contribution.
  std::array<long double,6> outwardEnergyFace{};
  std::uint64_t overwriteEvents = 0;
  std::uint64_t fluxSteps = 0;
  std::uint64_t remapSteps = 0;
  std::uint64_t diagnosticCellSamples = 0;
  double diagnosticSeconds = 0.;
};

// Fixed z-slab partition, Cowan c*dt=dz, one-cell motion toward -z after each
// complete field step. No z CPML or z PEC; only transverse CPML is installed.
// beginStep -> MPI magnetic -> original particle/J work -> MPI electric ->
// shift is mandatory. The MPI magnetic call must receive physicalElectricHalos.
// The last z face/half-step collar is sacrificial, not a particle sampling area.
// Particle exit, origin updates and fixed-lab diagnostics belong to Simulation.
//
// Shift uses one packed nearest-neighbour exchange, O(transverse area) buffers,
// and local memmove. Node-centred shared planes send old k=localNz-1, NOT the
// duplicated upper node k=localNz. The inserted global-left histories are zero.
// Current normally resets; resetCurrent=false exists for explicit state-remap
// checks and preserves all nine native components without changing arithmetic.
class MovingWindowFields {
public:
  MovingWindowFields(const EBGridGeometry& localGeometry, std::size_t globalNz,
      std::size_t localOffset, MPI_Comm communicator,
      const EBCPMLParameters& parameters, MovingWindowReference reference,
      bool allowNonuniformIntegerSlabs = false);
  ~MovingWindowFields();
  MovingWindowFields(const MovingWindowFields&) = delete;
  MovingWindowFields& operator=(const MovingWindowFields&) = delete;

  void installBoundary(EBFieldGrid& fields);
  // Source coupling: preview the zero-current electric boundary on reversible
  // CPML histories; then undo its field increment and commit once, without a
  // second Maxwell curl. Normal stepping never enters these methods.
  void previewElectricBoundary(bool enabled);
  void undoPreviewElectricBoundary(EBFieldGrid& fields);
  void commitPreparedElectricBoundary(EBFieldGrid& fields);
  // Configure before the first beginStep. Disabled by default: no diagnostic
  // field traversal, extra MPI exchange, or timer is performed when disabled.
  void enableEnergyLedger(bool enabled = true);
  // Optional read-only diagnostic hook after transverse CPML and before the
  // inlet Ex/Ey overwrite. No field/current mutation is permitted by this API.
  using BeforeInletObserver = void (*)(void*, const EBFieldGrid&);
  void setBeforeInletObserver(BeforeInletObserver observer, void* context);
  // Optional completion of the newly inserted global-left layer, before the
  // inserted/remap energies are measured. It must leave retained nodes intact.
  using NewLayerCompletion = void (*)(void*, EBFieldGrid&, std::size_t, std::ptrdiff_t);
  // With centered normal matching, the two Ez half-cells have the same
  // face average as the reference after completion of the interior half-cell.
  void setNewLayerCompletion(NewLayerCompletion completion,void* context,
                             bool centerNormalReference = false);
  using StageObserver = void (*)(void*, const char*, const EBFieldGrid&);
  void setStageObserver(StageObserver observer, void* context);
  using DiscreteObserver = void (*)(void*, const char*, const EBFieldGrid&,
                                   const EBElectricHaloView*);
  void setDiscreteObserver(DiscreteObserver observer,void* context);
  const MovingWindowFieldEnergyLedger& energyLedger() const;
  // Call after MPI electric completion, before shift. The boundary callback
  // precedes MPI interface updates, so cannot provide the full end-step flux.
  // This call is a no-op when disabled; enabled steps require it exactly once.
  void finishEnergyLedgerStep(const EBFieldGrid& fields);
  void beginStep(EBFieldGrid& fields, std::size_t step,
                 std::ptrdiff_t globalWindowOffset);
  EBElectricHaloView physicalElectricHalos(const EBFieldGrid& fields,
                                         std::size_t step) const;
  MovingWindowFieldShiftReport shift(EBFieldGrid& fields,
      std::size_t newStep, std::ptrdiff_t newWindowOffset,
      bool resetCurrent = true);
  // Collective exact ownership change after a completed shift (Idle only).
  // The global native state is copied by integer z index, including all nine
  // E/B/J components and each transverse CPML convolution history. New and
  // old grids must remain alive until this returns; each owns its own boundary.
  // Physical parameters, reference phase, cumulative diagnostics and callbacks
  // are preserved. This operation introduces no field-energy source/sink.
  std::unique_ptr<MovingWindowFields> repartition(const EBFieldGrid& oldFields,
      EBFieldGrid& nextFields, std::size_t newLocalOffset) const;
  std::size_t memoryBytes() const;
  std::uint64_t referenceValues() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
}
#endif
