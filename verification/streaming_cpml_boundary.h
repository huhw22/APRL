#ifndef APRL_VERIFICATION_STREAMING_CPML_BOUNDARY_H
#define APRL_VERIFICATION_STREAMING_CPML_BOUNDARY_H

#include <cstddef>
#include <cstdint>
#include <memory>

#include "eb_cpml.h"
#include "streaming_reference_boundary.h"

namespace aprl {
namespace verification {

// Verification-only total-field CPML with the existing frozen-reference feed
// at its exterior. The native CPML operates on TOTAL E/B and retains its
// original coefficients, Cowan derivatives, recurrences and zero initial
// histories. By default no reference is subtracted from fields, currents or
// histories. The opt-in residualInput mode instead presents E-Eref or B-Bref
// ONLY to the original CPML correction hook, with exact source restoration;
// the preceding native Maxwell advance and the hook's target remain total
// fields. That mode has one set of residual CPML histories, migrated through
// the same base interface. It is limited to z-only CPML and nonzero feed.
//
// After the native electric update, the existing ReferenceClamp supplies the
// two exterior Ex/Ey planes and reapplies transverse PEC. This deliberately
// replaces the native homogeneous z-face values with prescribed reference
// data: the composite boundary is reference-driven, not an unchanged
// homogeneous CPML termination. The default mode is total-field CPML; only
// the explicit residualInput option changes the CPML driving field.
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
                   bool zeroFeed = false,
                   bool residualInput = false,
                   StreamingReferenceCache* diagnosticCpmlReference = nullptr);
  ~ReferenceFedCPML() override;

  void beginStep(const EBFieldGrid& fields, unsigned int step);
  EBElectricHaloView electricHalos(const EBFieldGrid& fields,
                                  unsigned int step);
  void afterMagneticUpdate(EBFieldGrid& fields,
                           const EBElectricHaloView& halo) override;
  void afterElectricUpdate(EBFieldGrid& fields) override;
  void refreshAfterRemap(EBFieldGrid& fields, unsigned int stateStep);
  std::size_t memoryBytes() const override;
  std::uint64_t boundaryValues() const;
  std::uint64_t haloValues() const;
  bool residualInput() const;
  // Includes the one-time backup allocation and its fixed bookkeeping.
  // Already included by memoryBytes(); do not add it to grid memory twice.
  std::size_t residualWorkBytes() const;
  std::uint64_t residualInputValues() const;
  // Backup, reference queries/subtraction, and exact restoration only;
  // excludes the original CPML hook and exterior reference feed.
  double residualInputSeconds() const;

private:
  struct ResidualInputWork;
  StreamingReferenceBoundary feed_;
  std::unique_ptr<ResidualInputWork> residual_;
};

} // namespace verification
} // namespace aprl

#endif
