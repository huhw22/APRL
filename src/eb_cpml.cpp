#include "eb_cpml.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace fel
{
  namespace
  {
    struct CPMLPoint
    {
      Double inverseKappaMinusOne;
      Double decay;
      Double drive;
      bool active;

      CPMLPoint()
        : inverseKappaMinusOne(0.0), decay(0.0), drive(0.0),
          active(false)
      {}
    };

    struct CPMLProfile
    {
      unsigned int axis;
      std::vector<CPMLPoint> point;
      std::vector<std::size_t> activeIndices;

      CPMLProfile() : axis(0), point(), activeIndices() {}

      std::size_t memoryBytes() const
      {
        return point.capacity() * sizeof(CPMLPoint) +
          activeIndices.capacity() * sizeof(std::size_t);
      }
    };

    struct CPMLMemoryTerm
    {
      unsigned int axis;
      std::size_t dimensions[3];
      std::vector<std::size_t> baseByAxisIndex;
      std::vector<Double> history;

      CPMLMemoryTerm()
        : axis(0), baseByAxisIndex(), history()
      {
        dimensions[0] = dimensions[1] = dimensions[2] = 0;
      }

      CPMLMemoryTerm(unsigned int derivativeAxis,
                     const std::size_t shape[3],
                     const CPMLProfile& profile)
        : axis(derivativeAxis), baseByAxisIndex(), history()
      {
        for (unsigned int index = 0; index < 3; ++index)
          dimensions[index] = shape[index];
        if (profile.point.size() != dimensions[axis])
          throw std::invalid_argument(
            "CPML profile does not match its field component");
        const std::size_t missing = std::numeric_limits<std::size_t>::max();
        baseByAxisIndex.assign(dimensions[axis], missing);
        const std::size_t crossSection = axis == 0 ?
          dimensions[1] * dimensions[2] :
          (axis == 1 ? dimensions[0] * dimensions[2] :
                       dimensions[0] * dimensions[1]);
        history.assign(profile.activeIndices.size() * crossSection, 0.0);
        for (std::size_t ordinal = 0;
             ordinal < profile.activeIndices.size(); ++ordinal)
          baseByAxisIndex[profile.activeIndices[ordinal]] =
            ordinal * crossSection;
      }

      Double& value(std::size_t i, std::size_t j, std::size_t k)
      {
        const std::size_t axisIndex = axis == 0 ? i : (axis == 1 ? j : k);
        const std::size_t base = baseByAxisIndex[axisIndex];
        if (base == std::numeric_limits<std::size_t>::max())
          throw std::logic_error("Inactive CPML history was requested");
        const std::size_t crossIndex = axis == 0 ?
          k * dimensions[1] + j :
          (axis == 1 ? k * dimensions[0] + i :
                       j * dimensions[0] + i);
        return history[base + crossIndex];
      }

      std::size_t memoryBytes() const
      {
        return baseByAxisIndex.capacity() * sizeof(std::size_t) +
          history.capacity() * sizeof(Double);
      }
    };

    void componentShape(const EBGridGeometry& geometry, bool electric,
                        unsigned int component, std::size_t shape[3])
    {
      if (electric)
        {
          shape[0] = component == 0 ? geometry.nx : geometry.nx + 1;
          shape[1] = component == 1 ? geometry.ny : geometry.ny + 1;
          shape[2] = component == 2 ? geometry.nz : geometry.nz + 1;
        }
      else
        {
          shape[0] = component == 0 ? geometry.nx + 1 : geometry.nx;
          shape[1] = component == 1 ? geometry.ny + 1 : geometry.ny;
          shape[2] = component == 2 ? geometry.nz + 1 : geometry.nz;
        }
    }

    YeeComponent& electricComponent(EBFieldGrid& fields,
                                    unsigned int component)
    {
      if (component == 0) return fields.ex();
      if (component == 1) return fields.ey();
      return fields.ez();
    }

    YeeComponent& magneticComponent(EBFieldGrid& fields,
                                    unsigned int component)
    {
      if (component == 0) return fields.bx();
      if (component == 1) return fields.by();
      return fields.bz();
    }

    Double transformedCorrection(const CPMLPoint& coefficient,
                                 Double rawDifference,
                                 Double& history)
    {
      history = coefficient.decay * history +
                coefficient.drive * rawDifference;
      return coefficient.inverseKappaMinusOne * rawDifference + history;
    }

    template<typename Operation>
    void forActivePoints(const CPMLProfile& profile,
                         const std::size_t begin[3],
                         const std::size_t end[3],
                         Operation operation)
    {
      for (std::size_t ordinal = 0;
           ordinal < profile.activeIndices.size(); ++ordinal)
        {
          const std::size_t coordinate = profile.activeIndices[ordinal];
          if (coordinate < begin[profile.axis] ||
              coordinate >= end[profile.axis])
            continue;
          if (profile.axis == 0)
            {
              const std::size_t i = coordinate;
              for (std::size_t k = begin[2]; k < end[2]; ++k)
                for (std::size_t j = begin[1]; j < end[1]; ++j)
                  operation(i, j, k);
            }
          else if (profile.axis == 1)
            {
              const std::size_t j = coordinate;
              for (std::size_t k = begin[2]; k < end[2]; ++k)
                for (std::size_t i = begin[0]; i < end[0]; ++i)
                  operation(i, j, k);
            }
          else
            {
              const std::size_t k = coordinate;
              for (std::size_t j = begin[1]; j < end[1]; ++j)
                for (std::size_t i = begin[0]; i < end[0]; ++i)
                  operation(i, j, k);
            }
        }
    }

    CPMLProfile makeProfile(unsigned int axis, bool electric,
                            const EBGridGeometry& geometry,
                            std::size_t globalNz,
                            std::size_t localZCellOffset,
                            const EBCPMLParameters& parameters)
    {
      CPMLProfile result;
      result.axis = axis;
      const std::size_t localCells = axis == 0 ? geometry.nx :
        (axis == 1 ? geometry.ny : geometry.nz);
      const std::size_t globalCells = axis == 2 ? globalNz : localCells;
      const std::size_t globalOffset = axis == 2 ? localZCellOffset : 0;
      const std::size_t count = localCells + (electric ? 1 : 0);
      result.point.resize(count);
      const std::size_t thicknessCells = parameters.cells[axis];
      if (thicknessCells == 0) return result;

      const Double thickness = static_cast<Double>(thicknessCells);
      const Double sigmaRateMaximum =
        parameters.maximumConductivityRate(axis, geometry);

      for (std::size_t localIndex = 0; localIndex < count; ++localIndex)
        {
          const Double position = static_cast<Double>(globalOffset +
            localIndex) + (electric ? 0.0 : 0.5);
          Double depth = 0.0;
          if (position < thickness)
            depth = (thickness - position) / thickness;
          const Double upperStart =
            static_cast<Double>(globalCells - thicknessCells);
          if (position > upperStart)
            depth = std::max(depth, (position - upperStart) / thickness);
          if (!(depth > 0.0)) continue;
          depth = std::min(1.0, depth);

          const Double graded = std::pow(depth, parameters.polynomialOrder);
          const Double sigmaRate = sigmaRateMaximum * graded;
          const Double kappa = 1.0 +
            (parameters.kappaMax - 1.0) * graded;
          const Double alphaRate = parameters.alphaFraction *
            sigmaRateMaximum * (1.0 - depth);
          CPMLPoint& coefficient = result.point[localIndex];
          coefficient.inverseKappaMinusOne = 1.0 / kappa - 1.0;
          coefficient.decay = std::exp(
            -(sigmaRate / kappa + alphaRate) * geometry.dt);
          const Double denominator = sigmaRate * kappa +
                                     alphaRate * kappa * kappa;
          coefficient.drive = denominator > 0.0 ?
            sigmaRate * (coefficient.decay - 1.0) / denominator : 0.0;
          coefficient.active = true;
          result.activeIndices.push_back(localIndex);
        }
      return result;
    }
  }

  EBCPMLParameters::EBCPMLParameters()
    : polynomialOrder(3.0), targetReflection(1.0e-8),
      kappaMax(8.0), alphaFraction(0.0)
  {
    cells[0] = cells[1] = cells[2] = 0;
  }

  bool EBCPMLParameters::enabled() const
  {
    return cells[0] > 0 || cells[1] > 0 || cells[2] > 0;
  }

  Double EBCPMLParameters::maximumConductivityRate(
      unsigned int axis, const EBGridGeometry& geometry) const
  {
    if (axis >= 3)
      throw std::out_of_range("CPML conductivity axis must be x, y, or z");
    if (cells[axis] == 0) return 0.0;
    const Double spacing = axis == 0 ? geometry.dx :
      (axis == 1 ? geometry.dy : geometry.dz);
    if (!(spacing > 0.0) || !std::isfinite(spacing) ||
        !(polynomialOrder > 0.0) ||
        !(targetReflection > 0.0 && targetReflection < 1.0))
      throw std::invalid_argument(
        "Cannot evaluate an invalid CPML conductivity profile");
    const Double thickness = static_cast<Double>(cells[axis]) * spacing;
    /* sigma_max/epsilon0 from the polynomial CPML estimate
     * sigma_max=-(m+1)ln(R)/(2 eta0 d), using 1/(eta0 epsilon0)=c. */
    return -(polynomialOrder + 1.0) * SI::c *
      std::log(targetReflection) / (2.0 * thickness);
  }

  struct EBConvolutionalPML::Impl
  {
    EBGridGeometry geometry;
    std::size_t globalNz;
    std::size_t localZCellOffset;
    EBCPMLParameters parameters;
    bool lowerZPhysical;
    bool upperZPhysical;
    CPMLProfile magneticProfile[3];
    CPMLProfile electricProfile[3];
    std::unique_ptr<CPMLMemoryTerm> magneticTerm[3][3];
    std::unique_ptr<CPMLMemoryTerm> electricTerm[3][3];
    PerfectElectricConductorBoundary termination;

    Impl(const EBGridGeometry& geometryValue,
         std::size_t globalNzValue,
         std::size_t localZCellOffsetValue,
         const EBCPMLParameters& parametersValue,
         bool lowerZPhysicalValue,
         bool upperZPhysicalValue)
      : geometry(geometryValue), globalNz(globalNzValue),
        localZCellOffset(localZCellOffsetValue), parameters(parametersValue),
        lowerZPhysical(lowerZPhysicalValue),
        upperZPhysical(upperZPhysicalValue),
        termination(lowerZPhysicalValue, upperZPhysicalValue)
    {
      if (!parameters.enabled())
        throw std::invalid_argument(
          "CPML requires at least one nonzero boundary thickness");
      if (!(parameters.polynomialOrder > 0.0) ||
          !(parameters.targetReflection > 0.0 &&
            parameters.targetReflection < 1.0) ||
          !(parameters.kappaMax >= 1.0) ||
          !(parameters.alphaFraction >= 0.0) ||
          !std::isfinite(parameters.polynomialOrder) ||
          !std::isfinite(parameters.targetReflection) ||
          !std::isfinite(parameters.kappaMax) ||
          !std::isfinite(parameters.alphaFraction))
        throw std::invalid_argument("Invalid CPML grading parameters");
      const std::size_t globalCells[3] = {
        geometry.nx, geometry.ny, globalNz
      };
      for (unsigned int axis = 0; axis < 3; ++axis)
        if (parameters.cells[axis] > 0)
          {
            if (parameters.cells[axis] < 2)
              throw std::invalid_argument(
                "Each enabled CPML direction needs at least two cells");
            if (2 * parameters.cells[axis] >= globalCells[axis])
              throw std::invalid_argument(
                "Opposite CPML slabs must leave a non-PML interior");
          }
      const std::size_t zThickness = parameters.cells[2];
      if (zThickness > 0)
        {
          if (!lowerZPhysical && localZCellOffset < zThickness)
            throw std::invalid_argument(
              "The lower z CPML slab cannot span MPI ranks yet; reduce MPI ranks or CPML cells");
          if (!upperZPhysical && localZCellOffset + geometry.nz >
              globalNz - zThickness)
            throw std::invalid_argument(
              "The upper z CPML slab cannot span MPI ranks yet; reduce MPI ranks or CPML cells");
        }

      for (unsigned int axis = 0; axis < 3; ++axis)
        {
          magneticProfile[axis] = makeProfile(
            axis, false, geometry, globalNz, localZCellOffset, parameters);
          electricProfile[axis] = makeProfile(
            axis, true, geometry, globalNz, localZCellOffset, parameters);
        }

      const unsigned int derivativePairs[6][2] = {
        {0, 1}, {0, 2}, {1, 2}, {1, 0}, {2, 0}, {2, 1}
      };
      for (unsigned int pair = 0; pair < 6; ++pair)
        {
          const unsigned int component = derivativePairs[pair][0];
          const unsigned int derivative = derivativePairs[pair][1];
          std::size_t shape[3];
          componentShape(geometry, false, component, shape);
          magneticTerm[component][derivative].reset(
            new CPMLMemoryTerm(
              derivative, shape, magneticProfile[derivative]));
          componentShape(geometry, true, component, shape);
          electricTerm[component][derivative].reset(
            new CPMLMemoryTerm(
              derivative, shape, electricProfile[derivative]));
        }
    }

    std::size_t memoryBytes() const
    {
      std::size_t result = 0;
      for (unsigned int axis = 0; axis < 3; ++axis)
        {
          result += magneticProfile[axis].memoryBytes();
          result += electricProfile[axis].memoryBytes();
          for (unsigned int derivative = 0; derivative < 3; ++derivative)
            {
              if (magneticTerm[axis][derivative])
                result += magneticTerm[axis][derivative]->memoryBytes();
              if (electricTerm[axis][derivative])
                result += electricTerm[axis][derivative]->memoryBytes();
            }
        }
      return result;
    }
  };

  EBConvolutionalPML::EBConvolutionalPML(
      const EBGridGeometry& localGeometry,
      std::size_t globalNz,
      std::size_t localZCellOffset,
      const EBCPMLParameters& parameters,
      bool lowerZPhysical,
      bool upperZPhysical)
    : impl_(new Impl(localGeometry, globalNz, localZCellOffset,
                     parameters, lowerZPhysical, upperZPhysical))
  {}

  EBConvolutionalPML::~EBConvolutionalPML() {}

  void EBConvolutionalPML::afterMagneticUpdate(
      EBFieldGrid& fields, const EBElectricHaloView& halo)
  {
    const EBGridGeometry& geometry = fields.geometry();
    const Double scale[3] = {
      geometry.dt / geometry.dx,
      geometry.dt / geometry.dy,
      geometry.dt / geometry.dz
    };

    const auto correctTerm = [&](unsigned int targetComponent,
                                 unsigned int sourceComponent,
                                 unsigned int derivativeAxis,
                                 Double sign)
    {
      CPMLProfile& profile = impl_->magneticProfile[derivativeAxis];
      if (profile.activeIndices.empty()) return;
      CPMLMemoryTerm& term =
        *impl_->magneticTerm[targetComponent][derivativeAxis];
      YeeComponent& target = magneticComponent(fields, targetComponent);
      const YeeComponent& source = electricComponent(fields, sourceComponent);
      const std::size_t begin[3] = {0, 0, 0};
      const std::size_t end[3] = {target.nx(), target.ny(), target.nz()};
      forActivePoints(profile, begin, end,
        [&](std::size_t i, std::size_t j, std::size_t k)
        {
          std::ptrdiff_t plus[3] = {
            static_cast<std::ptrdiff_t>(i),
            static_cast<std::ptrdiff_t>(j),
            static_cast<std::ptrdiff_t>(k)
          };
          ++plus[derivativeAxis];
          Double difference = 0.0;
          if (geometry.solver == EBMaxwellSolver::CowanZ)
            difference = fields.smoothedElectric(
              source, sourceComponent, derivativeAxis,
              plus[0], plus[1], plus[2], halo) -
              fields.smoothedElectric(
              source, sourceComponent, derivativeAxis,
              static_cast<std::ptrdiff_t>(i),
              static_cast<std::ptrdiff_t>(j),
              static_cast<std::ptrdiff_t>(k), halo);
          else
            difference = derivativeAxis == 0 ?
              source(i + 1, j, k) - source(i, j, k) :
              (derivativeAxis == 1 ?
                source(i, j + 1, k) - source(i, j, k) :
                source(i, j, k + 1) - source(i, j, k));
          Double& history = term.value(i, j, k);
          const Double correction = transformedCorrection(
            profile.point[derivativeAxis == 0 ? i :
              (derivativeAxis == 1 ? j : k)], difference, history);
          target(i, j, k) += sign * scale[derivativeAxis] * correction;
        });
    };

    correctTerm(0, 2, 1, -1.0);
    correctTerm(0, 1, 2,  1.0);
    correctTerm(1, 0, 2, -1.0);
    correctTerm(1, 2, 0,  1.0);
    correctTerm(2, 1, 0, -1.0);
    correctTerm(2, 0, 1,  1.0);
  }

  void EBConvolutionalPML::afterElectricUpdate(EBFieldGrid& fields)
  {
    const EBGridGeometry& geometry = fields.geometry();
    const Double curlScale = geometry.dt /
      (SI::mu0 * SI::epsilon0);
    const Double scale[3] = {
      curlScale / geometry.dx,
      curlScale / geometry.dy,
      curlScale / geometry.dz
    };

    const auto correctTerm = [&](unsigned int targetComponent,
                                 unsigned int sourceComponent,
                                 unsigned int derivativeAxis,
                                 Double sign)
    {
      CPMLProfile& profile = impl_->electricProfile[derivativeAxis];
      if (profile.activeIndices.empty()) return;
      CPMLMemoryTerm& term =
        *impl_->electricTerm[targetComponent][derivativeAxis];
      YeeComponent& target = electricComponent(fields, targetComponent);
      const YeeComponent& source = magneticComponent(fields, sourceComponent);
      std::size_t begin[3] = {0, 0, 0};
      std::size_t end[3] = {target.nx(), target.ny(), target.nz()};
      if (targetComponent == 0)
        {
          begin[1] = 1;
          end[1] = geometry.ny;
          begin[2] = impl_->lowerZPhysical ? 1 : 0;
          end[2] = impl_->upperZPhysical ? geometry.nz : geometry.nz + 1;
        }
      else if (targetComponent == 1)
        {
          begin[0] = 1;
          end[0] = geometry.nx;
          begin[2] = impl_->lowerZPhysical ? 1 : 0;
          end[2] = impl_->upperZPhysical ? geometry.nz : geometry.nz + 1;
        }
      else
        {
          begin[0] = 1;
          end[0] = geometry.nx;
          begin[1] = 1;
          end[1] = geometry.ny;
        }
      forActivePoints(profile, begin, end,
        [&](std::size_t i, std::size_t j, std::size_t k)
        {
          const Double difference = derivativeAxis == 0 ?
            source(i, j, k) - source(i - 1, j, k) :
            (derivativeAxis == 1 ?
              source(i, j, k) - source(i, j - 1, k) :
              source(i, j, k) - source(i, j, k - 1));
          Double& history = term.value(i, j, k);
          const Double correction = transformedCorrection(
            profile.point[derivativeAxis == 0 ? i :
              (derivativeAxis == 1 ? j : k)], difference, history);
          target(i, j, k) += sign * scale[derivativeAxis] * correction;
        });
    };

    correctTerm(0, 2, 1,  1.0);
    correctTerm(0, 1, 2, -1.0);
    correctTerm(1, 0, 2,  1.0);
    correctTerm(1, 2, 0, -1.0);
    correctTerm(2, 1, 0,  1.0);
    correctTerm(2, 0, 1, -1.0);

    impl_->termination.afterElectricUpdate(fields);
  }

  std::size_t EBConvolutionalPML::memoryBytes() const
  {
    return impl_->memoryBytes();
  }
}
