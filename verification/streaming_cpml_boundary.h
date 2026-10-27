#ifndef APRL_VERIFICATION_STREAMING_CPML_BOUNDARY_H
#define APRL_VERIFICATION_STREAMING_CPML_BOUNDARY_H

#include <cstddef>
#include <cstdint>

#include "eb_cpml.h"
#include "streaming_reference_boundary.h"

namespace aprl {
namespace verification {

// Verification-only total-field CPML with the existing frozen-reference feed
// at its exterior. The native CPML operates on TOTAL E/B and retains its
// original coefficients, Cowan derivatives, recurrences and zero initial
// histories. No reference is subtracted from fields, currents or histories.
//
// After the native electric update, the existing ReferenceClamp supplies the
// two exterior Ex/Ey planes and reapplies transverse PEC. This deliberately
// replaces the native homogeneous z-face values with prescribed reference
// data: the composite boundary is reference-driven, not an unchanged
// homogeneous CPML termination and not a residual-field CPML.
//
// The base object is a single-rank whole-window CPML (globalNz=local nz,
// localZCellOffset=0, both physical z faces). MovingWindowRemap can access its
// native histories through the existing EBConvolutionalPML base interface.
// WorldPosition history migration remains the caller's responsibility.
//
// Before each complete field step, call beginStep(fields,n), then supply
// electricHalos(fields,n) to fields.advanceMagnetic(). Push/deposit normally
// and call fields.advanceElectric(); only afterward may the caller move
// arrays and update offset. Native E/B time levels are unchanged.
class ReferenceFedCPML : public EBConvolutionalPML {
public:
  ReferenceFedCPML(const EBGridGeometry& geometry,
                   const EBCPMLParameters& parameters,
                   StreamingReferenceCache& reference,
                   const std::size_t& offset,
                   bool zeroFeed = false);
  ~ReferenceFedCPML() override = default;

  void beginStep(const EBFieldGrid& fields, unsigned int step);
  EBElectricHaloView electricHalos(const EBFieldGrid& fields,
                                  unsigned int step);
  void afterMagneticUpdate(EBFieldGrid& fields,
                           const EBElectricHaloView& halo) override;
  void afterElectricUpdate(EBFieldGrid& fields) override;
  std::size_t memoryBytes() const override;
  std::uint64_t boundaryValues() const;
  std::uint64_t haloValues() const;

private:
  StreamingReferenceBoundary feed_;
};

} // namespace verification
} // namespace aprl

#endif
