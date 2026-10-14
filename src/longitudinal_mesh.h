#ifndef LONGITUDINAL_MESH_H
#define LONGITUDINAL_MESH_H

#include <cstddef>

namespace aprl
{
  /* Return the boosted-frame distance travelled by a forward light ray
   * before it meets a fixed laboratory plane deltaLab metres downstream of
   * its emission event.  This is
   *
   *   deltaLab / (gamma + sqrt(gamma^2 - 1)),
   *
   * evaluated without forming gamma^2 or subtracting nearly equal values.
   */
  double forwardLightIntersectionAdvance(double deltaLab, double gamma);

  /* Pure longitudinal geometry input.  Every position is expressed in one
   * common, unshifted boosted-frame coordinate system.  The ballistic
   * envelope includes every particle position that must remain physical
   * during the requested run.  radiationUpper is the furthest forward
   * observation event that must precede the +z absorber.
   */
  struct LongitudinalLayoutInput
  {
    double initialMinimum;
    double initialMaximum;
    double initialMean;
    double initialRms;

    double ballisticMinimum;
    double ballisticMaximum;
    double radiationUpper;

    double dz;
    std::size_t cpmlCells;
    std::size_t leftGuardCells;
    std::size_t rightGuardCells;
    double selfFieldPaddingRms;
    std::size_t alignmentCells;
    double meshCenter;

    LongitudinalLayoutInput();
  };

  /* Derived mesh and translated boundaries.  [outerLower, outerUpper] is the
   * complete mesh; [physicalLower, physicalUpper] excludes CPML.  The
   * required interval includes guards.  It is anchored to physicalLower, so
   * cell-ceiling and alignment slack are deliberately placed only on +z.
   */
  struct LongitudinalLayout
  {
    std::size_t totalCells;
    std::size_t physicalCells;
    std::size_t minimumPhysicalCells;
    std::size_t alignmentExtraCells;

    double extent;
    double origin;
    double outerLower;
    double outerUpper;
    double physicalLower;
    double physicalUpper;
    double translation;

    double initialMinimum;
    double initialMaximum;
    double initialMean;
    double ballisticMinimum;
    double ballisticMaximum;
    double radiationUpper;
    double selfFieldLower;
    double selfFieldUpper;
    double contentLower;
    double contentUpper;
    double requiredLower;
    double requiredUpper;
    double rightRoundingPadding;
  };

  /* Runtime provenance attached to an automatically derived layout.  The
   * pure geometry result above deliberately knows nothing about beamline
   * coordinates; this wrapper records the physical events used by the
   * simulation when it resolved the otherwise omitted z cell count. */
  struct LongitudinalLayoutReport
  {
    LongitudinalLayout layout;
    double observationTargetLab;
    double initialHeadLab;
    double ballisticStopTime;
    double plannedBallisticTime;
    double radiationAdvance;
  };

  LongitudinalLayout deriveLongitudinalLayout(
      const LongitudinalLayoutInput& input);
}

#endif
