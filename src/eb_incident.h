#ifndef DIRECT_EB_INCIDENT_H
#define DIRECT_EB_INCIDENT_H

#include <cstddef>

#include "eb_sources.h"

namespace aprl
{
  /* Node-aligned closed total-field volume.  lower and upper are inclusive
   * primal-node indices and must leave at least one scattered-field cell
   * between the TF/SF surface and every physical grid boundary. */
  struct EBTFSFRegion
  {
    std::size_t lower[3];
    std::size_t upper[3];

    EBTFSFRegion();
    EBTFSFRegion(std::size_t lowerX, std::size_t upperX,
                 std::size_t lowerY, std::size_t upperY,
                 std::size_t lowerZ, std::size_t upperZ);
  };

  /* Closed-surface E/B TF/SF injector.
   *
   * The incident field is present in the total-field volume and absent in the
   * surrounding scattered-field volume.  Surface corrections are applied
   * after the corresponding Yee curl update.  No auxiliary full-domain field
   * arrays or per-surface index lists are allocated. */
  class EBMaxwellIncidentInjector
  {
  public:
    EBMaxwellIncidentInjector(const EBGridGeometry& geometry,
                              const EBTFSFRegion& region);

    /* z-slab view of a global TF/SF volume.  geometry is the local slab,
     * globalZCellOffset is its first global z-cell, and gridOriginBox passed
     * to the sampling methods is the physical origin of that local slab.
     * No surface or volume buffers are allocated. */
    EBMaxwellIncidentInjector(const EBGridGeometry& localGeometry,
                              const EBTFSFRegion& globalRegion,
                              std::size_t globalZCellOffset,
                              std::size_t globalNz);

    const EBTFSFRegion& region() const;
    std::size_t memoryBytes() const;

    /* Initialize E^n and B^(n-1/2) only inside the total-field volume. */
    void initialize(EBFieldGrid& grid,
                    const SIFieldSourceSet& sources,
                    const FieldVector<Double>& gridOriginBox,
                    Double timeEBox,
                    const BoostFrameTransform& frame,
                    bool clearFirst = true) const;

    /* Call immediately after the matching Yee update. */
    void correctAfterMagneticUpdate(
        EBFieldGrid& grid,
        const SIFieldSourceSet& sources,
        const FieldVector<Double>& gridOriginBox,
        Double timeEBox,
        const BoostFrameTransform& frame) const;
    void correctAfterElectricUpdate(
        EBFieldGrid& grid,
        const SIFieldSourceSet& sources,
        const FieldVector<Double>& gridOriginBox,
        Double timeEBox,
        const BoostFrameTransform& frame) const;

    /* Convenience serial/local-slab step. MPI orchestration may call the two
     * corrections separately around its halo exchanges. */
    void advance(EBFieldGrid& grid,
                 const SIFieldSourceSet& sources,
                 const FieldVector<Double>& gridOriginBox,
                 Double timeEBox,
                 const BoostFrameTransform& frame) const;

  private:
    void validate() const;
    void validateGrid(const EBFieldGrid& grid) const;
    FieldVector<Double> point(
        const FieldVector<Double>& origin,
        std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k,
        Double offsetX, Double offsetY, Double offsetZ) const;
    Double incidentElectric(
        const SIFieldSourceSet& sources,
        const FieldVector<Double>& origin,
        std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k,
        Double offsetX, Double offsetY, Double offsetZ,
        unsigned int component, Double time,
        const BoostFrameTransform& frame) const;
    Double incidentMagnetic(
        const SIFieldSourceSet& sources,
        const FieldVector<Double>& origin,
        std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k,
        Double offsetX, Double offsetY, Double offsetZ,
        unsigned int component, Double time,
        const BoostFrameTransform& frame) const;

    EBGridGeometry geometry_;
    EBTFSFRegion region_;
    std::size_t globalZCellOffset_;
    std::size_t globalNz_;
    std::size_t localNodeZ0_;
    std::size_t localNodeZ1_;
    std::size_t localCellZ0_;
    std::size_t localCellZ1_;
    bool hasNodeIntersection_;
    bool hasCellIntersection_;
    bool correctLowerZMagnetic_;
    bool correctUpperZMagnetic_;
    bool correctLowerZElectric_;
    bool correctUpperZElectric_;
  };
}

#endif
