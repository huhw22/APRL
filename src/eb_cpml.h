#ifndef DIRECT_EB_CPML_H
#define DIRECT_EB_CPML_H

#include <cstddef>
#include <memory>
#include <vector>

#include "eb_field.h"

namespace aprl
{
  struct EBCPMLParameters
  {
    std::size_t cells[3];
    Double polynomialOrder;
    Double targetReflection;
    Double kappaMax;
    Double alphaFraction;

    EBCPMLParameters();
    bool enabled() const;
    Double maximumConductivityRate(unsigned int axis,
                                   const EBGridGeometry& geometry) const;
  };

  /* Unsplit complex-frequency-shifted convolutional PML.  Only convolution
   * histories whose derivative coordinate lies inside a PML slab are stored;
   * there are no full-domain auxiliary field copies.  For Cowan-z the
   * stretched magnetic derivative is D_i S_i(E), not a fallback Yee D_i(E). */
  class EBConvolutionalPML : public EBBoundaryOperator
  {
  public:
    EBConvolutionalPML(
        const EBGridGeometry& localGeometry,
        std::size_t globalNz,
        std::size_t localZCellOffset,
        const EBCPMLParameters& parameters,
        bool lowerZPhysical,
        bool upperZPhysical);
    virtual ~EBConvolutionalPML();

    virtual void afterMagneticUpdate(
        EBFieldGrid& fields, const EBElectricHaloView& halo);
    virtual void afterElectricUpdate(EBFieldGrid& fields);
    virtual std::size_t memoryBytes() const;

    /* Checkpoint order is magnetic then electric, each with derivative
     * pairs (component, axis): (0,1), (0,2), (1,2), (1,0), (2,0), (2,1).
     * Restore values into the existing history; never resize these arrays.
     * Geometry and CPML parameters must match before restoring a history. */
    std::size_t checkpointHistoryCount() const;
    std::vector<Double>& checkpointHistory(std::size_t index);
    const std::vector<Double>& checkpointHistory(std::size_t index) const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
}

#endif
