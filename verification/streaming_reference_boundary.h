#ifndef APRL_VERIFICATION_STREAMING_REFERENCE_BOUNDARY_H
#define APRL_VERIFICATION_STREAMING_REFERENCE_BOUNDARY_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include "eb_field.h"

namespace aprl {
namespace verification {

// Frozen deep copy of one initial Cowan E/B state. By default it has no
// connection to a full-domain evolving simulation. A diagnostic-only callback
// may override queries to replay native fields for error attribution; such a
// replay is not an independent reference supply for a production simulation.
// Component numbers are Ex,Ey,Ez,Bx,By,Bz.
// The cached state is E(0),B(-dt/2); stateStep=n queries E(n),B(n-1/2).
// A B query after its magnetic update therefore uses stateStep=n+1.
// Each component is translated in z by beta*c*stateStep*dt using linear
// interpolation on its own native staggered lattice. This approximation is
// not asserted to reproduce a moving point-particle CIC charge distribution.
class StreamingReferenceCache {
public:
  typedef Double (*DiagnosticQuery)(void* context, unsigned int component,
      std::size_t i, std::size_t j, Double globalK, Double stateStep);
  StreamingReferenceCache(const EBFieldGrid& initial, Double beta);
  Double value(unsigned int component, std::size_t i, std::size_t j,
               Double globalK, Double stateStep) const;
  // Diagnostic replay only. The callback receives the original native query,
  // before frozen translation, and owns its range/time validation. Neither
  // callback context nor its storage is owned or included in memoryBytes().
  // A null callback restores the unchanged frozen-cache query path. Successful
  // queries on either path contribute to the same cumulative sampleCount().
  void setDiagnosticQuery(DiagnosticQuery query, void* context=nullptr);
  const EBGridGeometry& geometry() const;
  Double beta() const;
  std::size_t memoryBytes() const;
  std::uint64_t sampleCount() const;

private:
  struct Component {
    std::size_t nx, ny, nz;
    std::vector<Double> values;
  };
  EBGridGeometry geometry_;
  Double beta_;
  Component components_[6];
  mutable std::uint64_t samples_;
  DiagnosticQuery diagnosticQuery_;
  void* diagnosticContext_;
};

enum class StreamingBoundaryMode {
  ReferenceClamp,
  ResidualMur,
  ZeroFeed
};

// Research boundary for a genuinely smaller z window. Transverse faces keep
// the production PEC convention. The z halo/Ex,Ey end faces use a frozen
// reference, with an optional first-order outgoing Mur condition on F-Fref.
// Mur is an approximation, especially for oblique/near-field content; guard
// sensitivity must be measured. This class is not moving CPML.
//
// Usage, with offset held fixed during each field step:
//   beginStep(fields,n);
//   fields.advanceMagnetic(electricHalos(fields,n));
//   ... push particles and deposit the actual current ...
//   fields.advanceElectric();
// The caller may then shift fields/update offset and fill incoming planes.
// beginStep snapshots the current end/adjacent residuals afresh, so there is
// no unshifted convolution history hidden in this boundary object.
class StreamingReferenceBoundary : public EBBoundaryOperator {
public:
  StreamingReferenceBoundary(const StreamingReferenceCache& reference,
                             const std::size_t& offset,
                             StreamingBoundaryMode mode);
  void beginStep(const EBFieldGrid& fields, unsigned int step);
  EBElectricHaloView electricHalos(const EBFieldGrid& fields,
                                  unsigned int step);
  void afterMagneticUpdate(EBFieldGrid& fields,
                           const EBElectricHaloView& halo) override;
  void afterElectricUpdate(EBFieldGrid& fields) override;
  // Re-establish the already prescribed clamp on the NEW exterior faces
  // after an integer remap, without advancing CPML or any field equation.
  void refreshAfterRemap(EBFieldGrid& fields, unsigned int stateStep);
  std::size_t memoryBytes() const override;
  std::uint64_t boundaryValues() const;
  std::uint64_t haloValues() const;

private:
  void validateGeometry(const EBFieldGrid& fields) const;
  void requirePrepared() const;
  Double referenceValue(unsigned int component, std::size_t i,
                        std::size_t j, Double globalK, Double stateStep) const;

  const StreamingReferenceCache& reference_;
  const std::size_t& offset_;
  StreamingBoundaryMode mode_;
  PerfectElectricConductorBoundary transversePec_;
  unsigned int step_;
  std::size_t stepOffset_;
  bool prepared_;
  std::vector<Double> oldEnd_[2][2];
  std::vector<Double> oldAdjacent_[2][2];
  std::vector<Double> lower_[3], upper_[3];
  std::uint64_t boundaryValues_, haloValues_;
};

}  // namespace verification
}  // namespace aprl
#endif
