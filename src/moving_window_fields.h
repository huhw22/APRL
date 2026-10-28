#ifndef APRL_MOVING_WINDOW_FIELDS_H
#define APRL_MOVING_WINDOW_FIELDS_H

#include <cstddef>
#include <cstdint>
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
      const EBCPMLParameters& parameters, MovingWindowReference reference);
  ~MovingWindowFields();
  MovingWindowFields(const MovingWindowFields&) = delete;
  MovingWindowFields& operator=(const MovingWindowFields&) = delete;

  void installBoundary(EBFieldGrid& fields);
  void beginStep(EBFieldGrid& fields, std::size_t step,
                 std::ptrdiff_t globalWindowOffset);
  EBElectricHaloView physicalElectricHalos(const EBFieldGrid& fields,
                                         std::size_t step) const;
  MovingWindowFieldShiftReport shift(EBFieldGrid& fields,
      std::size_t newStep, std::ptrdiff_t newWindowOffset,
      bool resetCurrent = true);
  std::size_t memoryBytes() const;
  std::uint64_t referenceValues() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
}
#endif
