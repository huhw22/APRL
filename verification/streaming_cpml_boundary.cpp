#include "streaming_cpml_boundary.h"

namespace aprl {
namespace verification {

ReferenceFedCPML::ReferenceFedCPML(
    const EBGridGeometry& geometry, const EBCPMLParameters& parameters,
    StreamingReferenceCache& reference, const std::size_t& offset,
    bool zeroFeed)
  : EBConvolutionalPML(geometry, geometry.nz, 0, parameters, true, true),
    feed_(reference, offset, zeroFeed ? StreamingBoundaryMode::ZeroFeed :
                                     StreamingBoundaryMode::ReferenceClamp) {}

void ReferenceFedCPML::beginStep(const EBFieldGrid& fields, unsigned int step) {
  feed_.beginStep(fields, step);
}

EBElectricHaloView ReferenceFedCPML::electricHalos(
    const EBFieldGrid& fields, unsigned int step) {
  return feed_.electricHalos(fields, step);
}

void ReferenceFedCPML::afterMagneticUpdate(
    EBFieldGrid& fields, const EBElectricHaloView& halo) {
  EBConvolutionalPML::afterMagneticUpdate(fields, halo);
  feed_.afterMagneticUpdate(fields, halo);
}

void ReferenceFedCPML::afterElectricUpdate(EBFieldGrid& fields) {
  EBConvolutionalPML::afterElectricUpdate(fields);
  feed_.afterElectricUpdate(fields);
}

std::size_t ReferenceFedCPML::memoryBytes() const {
  return EBConvolutionalPML::memoryBytes() + feed_.memoryBytes();
}

std::uint64_t ReferenceFedCPML::boundaryValues() const {
  return feed_.boundaryValues();
}

std::uint64_t ReferenceFedCPML::haloValues() const {
  return feed_.haloValues();
}

} // namespace verification
} // namespace aprl
