#ifndef APRL_VERIFICATION_HELD_CUT_BOUNDARY_H
#define APRL_VERIFICATION_HELD_CUT_BOUNDARY_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include "eb_cpml.h"
#include "moving_window_remap.h"

namespace aprl { namespace verification {

// A native field query, not a fractional-coordinate interpolation. At state n,
// E is stored at n*dt and B at (n-1/2)*dt. The provider owns time validation.
struct HeldCutReference {
  using Function = Double (*)(void*, unsigned int, std::size_t, std::size_t,
                             std::ptrdiff_t, unsigned int);
  void* context = nullptr;
  Function field = nullptr;
};

// Research-only left source / right causal cut. Only x/y CPML is permitted.
// Both z PEC flags are false. The upper z halo is deliberately explicit and
// may be poisoned for causal tests; it is NOT an absorbing boundary. Right
// end integer faces and a transient half-step collar are outside the trusted
// region. The caller must shift LEFT by one cell per complete step (or shrink
// its trusted region on an unshifted step), and must keep all particle/field
// sampling inside the corresponding trusted native-index range.
//
// beginStep snapshots lower E halos at n; the independent held donor may then
// advance. afterElectricUpdate reads donor Ex/Ey at n+1 at the lower endpoint.
// CPML hooks otherwise use the unmodified TOTAL fields, including the original
// transverse termination. No residual-field CPML, z CPML, or Gauss fix is used.
class HeldCutBoundary : public EBConvolutionalPML {
public:
  HeldCutBoundary(const EBGridGeometry&, const EBCPMLParameters&,
                  const std::size_t& globalOffset, HeldCutReference,
                  Double rightHaloValue = 0.);
  void beginStep(const EBFieldGrid&, unsigned int step);
  EBElectricHaloView electricHalos(const EBFieldGrid&, unsigned int step) const;
  void afterMagneticUpdate(EBFieldGrid&, const EBElectricHaloView&) override;
  void afterElectricUpdate(EBFieldGrid&) override;
  std::size_t memoryBytes() const override;
  std::uint64_t referenceValues() const;
  void setRightHaloValue(Double value);

private:
  enum class Phase { Idle, Magnetic, Electric };
  EBGridGeometry geometry_;
  const std::size_t& offset_;
  std::size_t stepOffset_;
  HeldCutReference reference_;
  Double rightValue_;
  unsigned int step_ = 0;
  Phase phase_ = Phase::Idle;
  std::vector<Double> lower_[3], upper_[3];
  std::uint64_t queries_ = 0;
  void validate(const EBFieldGrid&) const;
  void requirePhase(Phase) const;
};

// Copy native E/B at globalOffset+k. count==0 means all available target
// planes. firstK/count are applied separately to each staggered component.
// Current is not copied: deposition owns J at each physical step.
void copyHeldCutFields(EBFieldGrid& target, const EBFieldGrid& source,
                       std::size_t globalOffset, std::size_t firstK=0,
                       std::size_t count=0);

// Restore matching transverse CPML histories for the same z slab. A remap
// clears the newly created histories; call this with firstK=0,count=1 after
// a left shift to supply the donor's evolving side-PML state at the inlet.
// All retained histories are untouched. Both layouts must have zero z CPML.
void copyHeldCutHistories(EBFieldGrid& target, const MovingWindowRemap& targetLayout,
                          const EBFieldGrid& source, const MovingWindowRemap& sourceLayout,
                          std::size_t globalOffset, std::size_t firstK=0,
                          std::size_t count=0);

}} // namespace aprl::verification
#endif
