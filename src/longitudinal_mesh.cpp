#include "longitudinal_mesh.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace aprl
{
  namespace
  {
    void requireFinite(double value, const char* name)
    {
      if (!std::isfinite(value))
        throw std::invalid_argument(std::string(name) + " must be finite");
    }

    std::size_t checkedAdd(std::size_t first, std::size_t second,
                           const char* message)
    {
      if (second > std::numeric_limits<std::size_t>::max() - first)
        throw std::overflow_error(message);
      return first + second;
    }

    std::size_t cellsCovering(double length, double spacing)
    {
      const long double ratio = static_cast<long double>(length) /
                                static_cast<long double>(spacing);
      const long double rounded = std::ceil(ratio);
      if (!(rounded >= 1.0L) ||
          rounded > static_cast<long double>(
            std::numeric_limits<std::size_t>::max()))
        throw std::overflow_error(
          "Longitudinal physical cell count is outside size_t range");
      return static_cast<std::size_t>(rounded);
    }

    std::size_t roundUp(std::size_t value, std::size_t alignment)
    {
      const std::size_t remainder = value % alignment;
      if (remainder == 0) return value;
      return checkedAdd(value, alignment - remainder,
        "Aligned longitudinal cell count overflows size_t");
    }
  }

  double forwardLightIntersectionAdvance(double deltaLab, double gamma)
  {
    requireFinite(deltaLab, "Forward light laboratory distance");
    requireFinite(gamma, "Boost gamma");
    if (deltaLab < 0.0)
      throw std::invalid_argument(
        "Forward light laboratory distance must be nonnegative");
    if (gamma < 1.0)
      throw std::invalid_argument("Boost gamma must be at least one");

    /* gammaBeta/gamma = sqrt((1 - 1/gamma)(1 + 1/gamma)).
     * Dividing before adding the rapidity factors avoids overflow when gamma
     * is close to the largest representable double. */
    const long double gammaWide = static_cast<long double>(gamma);
    const long double inverseGamma = 1.0L / gammaWide;
    const long double beta = std::sqrt(std::max(0.0L,
      (1.0L - inverseGamma) * (1.0L + inverseGamma)));
    const long double advance =
      (static_cast<long double>(deltaLab) / gammaWide) / (1.0L + beta);
    const double result = static_cast<double>(advance);
    if (!std::isfinite(result) || result < 0.0)
      throw std::overflow_error(
        "Forward light intersection advance is not representable");
    return result;
  }

  LongitudinalLayoutInput::LongitudinalLayoutInput()
    : initialMinimum(0.0), initialMaximum(0.0), initialMean(0.0),
      initialRms(0.0), ballisticMinimum(0.0), ballisticMaximum(0.0),
      radiationUpper(0.0), dz(0.0), cpmlCells(0), leftGuardCells(0),
      rightGuardCells(0), selfFieldPaddingRms(0.0), alignmentCells(1),
      meshCenter(0.0)
  {}

  LongitudinalLayout deriveLongitudinalLayout(
      const LongitudinalLayoutInput& input)
  {
    requireFinite(input.initialMinimum, "Initial minimum");
    requireFinite(input.initialMaximum, "Initial maximum");
    requireFinite(input.initialMean, "Initial mean");
    requireFinite(input.initialRms, "Initial RMS");
    requireFinite(input.ballisticMinimum, "Ballistic minimum");
    requireFinite(input.ballisticMaximum, "Ballistic maximum");
    requireFinite(input.radiationUpper, "Radiation upper boundary");
    requireFinite(input.dz, "Longitudinal cell size");
    requireFinite(input.selfFieldPaddingRms, "Self-field RMS padding");
    requireFinite(input.meshCenter, "Mesh centre");

    if (input.initialMinimum > input.initialMaximum)
      throw std::invalid_argument(
        "Initial longitudinal minimum exceeds maximum");
    if (input.initialMean < input.initialMinimum ||
        input.initialMean > input.initialMaximum)
      throw std::invalid_argument(
        "Initial longitudinal mean lies outside the sampled support");
    if (input.ballisticMinimum > input.ballisticMaximum)
      throw std::invalid_argument(
        "Ballistic longitudinal minimum exceeds maximum");
    if (input.initialRms < 0.0)
      throw std::invalid_argument("Initial longitudinal RMS must be nonnegative");
    if (!(input.dz > 0.0))
      throw std::invalid_argument("Longitudinal cell size must be positive");
    if (input.selfFieldPaddingRms < 0.0)
      throw std::invalid_argument(
        "Self-field RMS padding must be nonnegative");
    if (input.alignmentCells == 0)
      throw std::invalid_argument(
        "Longitudinal cell alignment must be positive");

    const double selfFieldDistance =
      input.selfFieldPaddingRms * input.initialRms;
    if (!std::isfinite(selfFieldDistance))
      throw std::overflow_error("Self-field longitudinal padding overflowed");
    const double selfFieldLowerRaw = input.initialMean - selfFieldDistance;
    const double selfFieldUpperRaw = input.initialMean + selfFieldDistance;
    if (!std::isfinite(selfFieldLowerRaw) ||
        !std::isfinite(selfFieldUpperRaw))
      throw std::overflow_error("Self-field longitudinal boundary overflowed");

    const double contentLowerRaw = std::min(
      std::min(input.initialMinimum, input.ballisticMinimum),
      selfFieldLowerRaw);
    const double contentUpperRaw = std::max(
      std::max(std::max(input.initialMaximum, input.ballisticMaximum),
               input.radiationUpper),
      selfFieldUpperRaw);
    const double leftGuard = static_cast<double>(input.leftGuardCells) *
                             input.dz;
    const double rightGuard = static_cast<double>(input.rightGuardCells) *
                              input.dz;
    if (!std::isfinite(leftGuard) || !std::isfinite(rightGuard))
      throw std::overflow_error("Longitudinal guard distance overflowed");
    const double requiredLowerRaw = contentLowerRaw - leftGuard;
    const double requiredUpperRaw = contentUpperRaw + rightGuard;
    const double requiredExtent = requiredUpperRaw - requiredLowerRaw;
    if (!std::isfinite(requiredLowerRaw) ||
        !std::isfinite(requiredUpperRaw) ||
        !(requiredExtent > 0.0) || !std::isfinite(requiredExtent))
      throw std::invalid_argument(
        "Longitudinal content and guards do not form a finite positive interval");

    const std::size_t minimumPhysicalCells =
      cellsCovering(requiredExtent, input.dz);
    const std::size_t twiceCpml = checkedAdd(input.cpmlCells,
      input.cpmlCells, "Longitudinal CPML cell count overflows size_t");
    const std::size_t minimumTotalCells = checkedAdd(
      minimumPhysicalCells, twiceCpml,
      "Longitudinal total cell count overflows size_t");
    const std::size_t totalCells = roundUp(
      minimumTotalCells, input.alignmentCells);
    if (static_cast<long double>(totalCells) >
        std::ldexp(1.0L, std::numeric_limits<double>::digits))
      throw std::overflow_error(
        "Longitudinal total cell count exceeds the exactly representable double index range");
    if (totalCells <= twiceCpml)
      throw std::overflow_error(
        "Aligned longitudinal mesh leaves no physical cells");
    const std::size_t physicalCells = totalCells - twiceCpml;

    const double extent = static_cast<double>(totalCells) * input.dz;
    const double cpmlExtent = static_cast<double>(input.cpmlCells) *
                              input.dz;
    if (!std::isfinite(extent) || !std::isfinite(cpmlExtent))
      throw std::overflow_error("Longitudinal mesh extent overflowed");
    const double origin = input.meshCenter - 0.5 * extent;
    const double outerUpper = origin + extent;
    const double physicalLower = origin + cpmlExtent;
    const double physicalUpper = outerUpper - cpmlExtent;
    const double translation = physicalLower - requiredLowerRaw;
    if (!std::isfinite(origin) || !std::isfinite(outerUpper) ||
        !std::isfinite(physicalLower) || !std::isfinite(physicalUpper) ||
        !std::isfinite(translation))
      throw std::overflow_error(
        "Translated longitudinal mesh boundary overflowed");

    LongitudinalLayout result;
    result.totalCells = totalCells;
    result.physicalCells = physicalCells;
    result.minimumPhysicalCells = minimumPhysicalCells;
    result.alignmentExtraCells = totalCells - minimumTotalCells;
    result.extent = extent;
    result.origin = origin;
    result.outerLower = origin;
    result.outerUpper = outerUpper;
    result.physicalLower = physicalLower;
    result.physicalUpper = physicalUpper;
    result.translation = translation;

    result.initialMinimum = input.initialMinimum + translation;
    result.initialMaximum = input.initialMaximum + translation;
    result.initialMean = input.initialMean + translation;
    result.ballisticMinimum = input.ballisticMinimum + translation;
    result.ballisticMaximum = input.ballisticMaximum + translation;
    result.radiationUpper = input.radiationUpper + translation;
    result.selfFieldLower = selfFieldLowerRaw + translation;
    result.selfFieldUpper = selfFieldUpperRaw + translation;
    result.contentLower = contentLowerRaw + translation;
    result.contentUpper = contentUpperRaw + translation;
    result.requiredLower = requiredLowerRaw + translation;
    result.requiredUpper = requiredUpperRaw + translation;
    result.rightRoundingPadding = physicalUpper - result.requiredUpper;

    const double tolerance = 128.0 *
      std::numeric_limits<double>::epsilon() *
      std::max(1.0, std::max(std::abs(physicalLower),
                             std::abs(physicalUpper)));
    if (std::abs(result.requiredLower - physicalLower) > tolerance ||
        result.requiredUpper > physicalUpper + tolerance ||
        result.rightRoundingPadding < -tolerance)
      throw std::logic_error(
        "Derived longitudinal layout failed its containment invariant");
    return result;
  }
}
