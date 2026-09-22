#include "eb_incident.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fel
{
  EBTFSFRegion::EBTFSFRegion()
  {
    lower[0] = lower[1] = lower[2] = 0;
    upper[0] = upper[1] = upper[2] = 0;
  }

  EBTFSFRegion::EBTFSFRegion(
      std::size_t lowerX, std::size_t upperX,
      std::size_t lowerY, std::size_t upperY,
      std::size_t lowerZ, std::size_t upperZ)
  {
    lower[0] = lowerX; upper[0] = upperX;
    lower[1] = lowerY; upper[1] = upperY;
    lower[2] = lowerZ; upper[2] = upperZ;
  }

  EBMaxwellIncidentInjector::EBMaxwellIncidentInjector(
      const EBGridGeometry& geometry, const EBTFSFRegion& region)
    : EBMaxwellIncidentInjector(geometry, region, 0, geometry.nz)
  {}

  EBMaxwellIncidentInjector::EBMaxwellIncidentInjector(
      const EBGridGeometry& localGeometry,
      const EBTFSFRegion& globalRegion,
      std::size_t globalZCellOffset,
      std::size_t globalNz)
    : geometry_(localGeometry), region_(globalRegion),
      globalZCellOffset_(globalZCellOffset), globalNz_(globalNz),
      localNodeZ0_(0), localNodeZ1_(0),
      localCellZ0_(0), localCellZ1_(0),
      hasNodeIntersection_(false), hasCellIntersection_(false),
      correctLowerZMagnetic_(false), correctUpperZMagnetic_(false),
      correctLowerZElectric_(false), correctUpperZElectric_(false)
  {
    validate();

    const std::size_t slabBegin = globalZCellOffset_;
    const std::size_t slabEnd = slabBegin + geometry_.nz;
    const std::size_t nodeBegin = std::max(region_.lower[2], slabBegin);
    const std::size_t nodeEnd = std::min(region_.upper[2], slabEnd);
    if (nodeBegin <= nodeEnd)
      {
        hasNodeIntersection_ = true;
        localNodeZ0_ = nodeBegin - slabBegin;
        localNodeZ1_ = nodeEnd - slabBegin;
      }

    const std::size_t cellBegin = std::max(region_.lower[2], slabBegin);
    const std::size_t cellEnd = std::min(region_.upper[2], slabEnd);
    if (cellBegin < cellEnd)
      {
        hasCellIntersection_ = true;
        localCellZ0_ = cellBegin - slabBegin;
        localCellZ1_ = cellEnd - slabBegin;
      }

    const std::size_t lowerOutsideCell = region_.lower[2] - 1;
    const std::size_t upperOutsideCell = region_.upper[2];
    correctLowerZMagnetic_ =
      slabBegin <= lowerOutsideCell && lowerOutsideCell < slabEnd;
    correctUpperZMagnetic_ =
      slabBegin <= upperOutsideCell && upperOutsideCell < slabEnd;
    correctLowerZElectric_ =
      slabBegin <= region_.lower[2] && region_.lower[2] <= slabEnd;
    correctUpperZElectric_ =
      slabBegin <= region_.upper[2] && region_.upper[2] <= slabEnd;
  }

  const EBTFSFRegion& EBMaxwellIncidentInjector::region() const
  {
    return region_;
  }

  std::size_t EBMaxwellIncidentInjector::memoryBytes() const
  {
    return sizeof(*this);
  }

  void EBMaxwellIncidentInjector::validate() const
  {
    const std::size_t cells[3] = {geometry_.nx, geometry_.ny, globalNz_};
    for (unsigned int axis = 0; axis < 3; ++axis)
      if (region_.lower[axis] < 1 ||
          region_.upper[axis] <= region_.lower[axis] ||
          region_.upper[axis] >= cells[axis])
        throw std::invalid_argument(
          "TF/SF region must leave one scattered-field cell on every side");
    if (globalZCellOffset_ > globalNz_ ||
        geometry_.nz > globalNz_ - globalZCellOffset_)
      throw std::invalid_argument(
        "TF/SF local z slab lies outside the global E/B grid");
  }

  void EBMaxwellIncidentInjector::validateGrid(
      const EBFieldGrid& grid) const
  {
    const EBGridGeometry& actual = grid.geometry();
    if (actual.solver != EBMaxwellSolver::Yee ||
        geometry_.solver != EBMaxwellSolver::Yee)
      throw std::invalid_argument(
        "TF/SF incident injection currently supports only the Yee solver; Cowan-z requires a generalized boundary correction");
    if (actual.nx != geometry_.nx || actual.ny != geometry_.ny ||
        actual.nz != geometry_.nz || actual.dx != geometry_.dx ||
        actual.dy != geometry_.dy || actual.dz != geometry_.dz ||
        actual.dt != geometry_.dt || actual.solver != geometry_.solver)
      throw std::invalid_argument(
        "TF/SF injector geometry does not match the E/B grid");
  }

  FieldVector<Double> EBMaxwellIncidentInjector::point(
      const FieldVector<Double>& origin,
      std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k,
      Double offsetX, Double offsetY, Double offsetZ) const
  {
    FieldVector<Double> result(origin);
    result[0] += (static_cast<Double>(i) + offsetX) * geometry_.dx;
    result[1] += (static_cast<Double>(j) + offsetY) * geometry_.dy;
    result[2] += (static_cast<Double>(k) + offsetZ) * geometry_.dz;
    return result;
  }

  Double EBMaxwellIncidentInjector::incidentElectric(
      const SIFieldSourceSet& sources,
      const FieldVector<Double>& origin,
      std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k,
      Double offsetX, Double offsetY, Double offsetZ,
      unsigned int component, Double time,
      const BoostFrameTransform& frame) const
  {
    SIFieldValue fields;
    sources.sampleMaxwellIncidentBox(
      point(origin, i, j, k, offsetX, offsetY, offsetZ),
      time, frame, fields);
    return fields.electric[component];
  }

  Double EBMaxwellIncidentInjector::incidentMagnetic(
      const SIFieldSourceSet& sources,
      const FieldVector<Double>& origin,
      std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k,
      Double offsetX, Double offsetY, Double offsetZ,
      unsigned int component, Double time,
      const BoostFrameTransform& frame) const
  {
    SIFieldValue fields;
    sources.sampleMaxwellIncidentBox(
      point(origin, i, j, k, offsetX, offsetY, offsetZ),
      time, frame, fields);
    return fields.magnetic[component];
  }

  void EBMaxwellIncidentInjector::initialize(
      EBFieldGrid& grid, const SIFieldSourceSet& sources,
      const FieldVector<Double>& origin, Double timeE,
      const BoostFrameTransform& frame, bool clearFirst) const
  {
    validateGrid(grid);
    if (!std::isfinite(timeE))
      throw std::invalid_argument("TF/SF initialization time must be finite");
    if (clearFirst)
      grid.clearFields();

    SIFieldValue fields;
    FieldVector<Double> location(0.0);
    const std::size_t x0 = region_.lower[0];
    const std::size_t x1 = region_.upper[0];
    const std::size_t y0 = region_.lower[1];
    const std::size_t y1 = region_.upper[1];
    const std::size_t nodeZ0 = localNodeZ0_;
    const std::size_t nodeZ1 = localNodeZ1_;
    const std::size_t cellZ0 = localCellZ0_;
    const std::size_t cellZ1 = localCellZ1_;

    if (hasNodeIntersection_)
      {
        for (std::size_t k = nodeZ0; k <= nodeZ1; ++k)
          for (std::size_t j = y0; j <= y1; ++j)
            for (std::size_t i = x0; i < x1; ++i)
              {
                location = point(origin, i, j, k, 0.5, 0.0, 0.0);
                sources.sampleMaxwellIncidentBox(location, timeE, frame, fields);
                grid.ex()(i, j, k) += fields.electric[0];
              }
        for (std::size_t k = nodeZ0; k <= nodeZ1; ++k)
          for (std::size_t j = y0; j < y1; ++j)
            for (std::size_t i = x0; i <= x1; ++i)
              {
                location = point(origin, i, j, k, 0.0, 0.5, 0.0);
                sources.sampleMaxwellIncidentBox(location, timeE, frame, fields);
                grid.ey()(i, j, k) += fields.electric[1];
              }
      }
    if (hasCellIntersection_)
      for (std::size_t k = cellZ0; k < cellZ1; ++k)
        for (std::size_t j = y0; j <= y1; ++j)
          for (std::size_t i = x0; i <= x1; ++i)
            {
              location = point(origin, i, j, k, 0.0, 0.0, 0.5);
              sources.sampleMaxwellIncidentBox(location, timeE, frame, fields);
              grid.ez()(i, j, k) += fields.electric[2];
            }

    const Double timeB = timeE - 0.5 * geometry_.dt;
    if (hasCellIntersection_)
      {
        for (std::size_t k = cellZ0; k < cellZ1; ++k)
          for (std::size_t j = y0; j < y1; ++j)
            for (std::size_t i = x0; i <= x1; ++i)
              {
                location = point(origin, i, j, k, 0.0, 0.5, 0.5);
                sources.sampleMaxwellIncidentBox(location, timeB, frame, fields);
                grid.bx()(i, j, k) += fields.magnetic[0];
              }
        for (std::size_t k = cellZ0; k < cellZ1; ++k)
          for (std::size_t j = y0; j <= y1; ++j)
            for (std::size_t i = x0; i < x1; ++i)
              {
                location = point(origin, i, j, k, 0.5, 0.0, 0.5);
                sources.sampleMaxwellIncidentBox(location, timeB, frame, fields);
                grid.by()(i, j, k) += fields.magnetic[1];
              }
      }
    if (hasNodeIntersection_)
      for (std::size_t k = nodeZ0; k <= nodeZ1; ++k)
        for (std::size_t j = y0; j < y1; ++j)
          for (std::size_t i = x0; i < x1; ++i)
            {
              location = point(origin, i, j, k, 0.5, 0.5, 0.0);
              sources.sampleMaxwellIncidentBox(location, timeB, frame, fields);
              grid.bz()(i, j, k) += fields.magnetic[2];
            }
  }

  void EBMaxwellIncidentInjector::correctAfterMagneticUpdate(
      EBFieldGrid& grid, const SIFieldSourceSet& sources,
      const FieldVector<Double>& origin, Double timeE,
      const BoostFrameTransform& frame) const
  {
    validateGrid(grid);
    const std::size_t x0 = region_.lower[0];
    const std::size_t x1 = region_.upper[0];
    const std::size_t y0 = region_.lower[1];
    const std::size_t y1 = region_.upper[1];
    const Double dtdx = geometry_.dt / geometry_.dx;
    const Double dtdy = geometry_.dt / geometry_.dy;
    const Double dtdz = geometry_.dt / geometry_.dz;

    /* Bx: +dEz/dy - dEy/dz correction. */
    if (hasCellIntersection_)
      for (std::size_t k = localCellZ0_; k < localCellZ1_; ++k)
        for (std::size_t i = x0; i <= x1; ++i)
          {
            grid.bx()(i, y0 - 1, k) += dtdy * incidentElectric(
              sources, origin, i, y0, k, 0.0, 0.0, 0.5,
              2, timeE, frame);
            grid.bx()(i, y1, k) -= dtdy * incidentElectric(
              sources, origin, i, y1, k, 0.0, 0.0, 0.5,
              2, timeE, frame);
          }
    if (correctLowerZMagnetic_)
      {
        const std::size_t outside =
          region_.lower[2] - 1 - globalZCellOffset_;
        const std::ptrdiff_t surface =
          static_cast<std::ptrdiff_t>(region_.lower[2]) -
          static_cast<std::ptrdiff_t>(globalZCellOffset_);
        for (std::size_t j = y0; j < y1; ++j)
          for (std::size_t i = x0; i <= x1; ++i)
            grid.bx()(i, j, outside) -= dtdz * incidentElectric(
              sources, origin, i, j, surface, 0.0, 0.5, 0.0,
              1, timeE, frame);
      }
    if (correctUpperZMagnetic_)
      {
        const std::size_t outside =
          region_.upper[2] - globalZCellOffset_;
        const std::ptrdiff_t surface =
          static_cast<std::ptrdiff_t>(outside);
        for (std::size_t j = y0; j < y1; ++j)
          for (std::size_t i = x0; i <= x1; ++i)
            grid.bx()(i, j, outside) += dtdz * incidentElectric(
              sources, origin, i, j, surface, 0.0, 0.5, 0.0,
              1, timeE, frame);
      }

    /* By: +dEx/dz - dEz/dx correction. */
    if (correctLowerZMagnetic_)
      {
        const std::size_t outside =
          region_.lower[2] - 1 - globalZCellOffset_;
        const std::ptrdiff_t surface =
          static_cast<std::ptrdiff_t>(region_.lower[2]) -
          static_cast<std::ptrdiff_t>(globalZCellOffset_);
        for (std::size_t j = y0; j <= y1; ++j)
          for (std::size_t i = x0; i < x1; ++i)
            grid.by()(i, j, outside) += dtdz * incidentElectric(
              sources, origin, i, j, surface, 0.5, 0.0, 0.0,
              0, timeE, frame);
      }
    if (correctUpperZMagnetic_)
      {
        const std::size_t outside =
          region_.upper[2] - globalZCellOffset_;
        const std::ptrdiff_t surface =
          static_cast<std::ptrdiff_t>(outside);
        for (std::size_t j = y0; j <= y1; ++j)
          for (std::size_t i = x0; i < x1; ++i)
            grid.by()(i, j, outside) -= dtdz * incidentElectric(
              sources, origin, i, j, surface, 0.5, 0.0, 0.0,
              0, timeE, frame);
      }
    if (hasCellIntersection_)
      for (std::size_t k = localCellZ0_; k < localCellZ1_; ++k)
        for (std::size_t j = y0; j <= y1; ++j)
          {
            grid.by()(x0 - 1, j, k) -= dtdx * incidentElectric(
              sources, origin, x0, j, k, 0.0, 0.0, 0.5,
              2, timeE, frame);
            grid.by()(x1, j, k) += dtdx * incidentElectric(
              sources, origin, x1, j, k, 0.0, 0.0, 0.5,
              2, timeE, frame);
          }

    /* Bz: +dEy/dx - dEx/dy correction. */
    if (hasNodeIntersection_)
      {
        for (std::size_t k = localNodeZ0_; k <= localNodeZ1_; ++k)
          for (std::size_t j = y0; j < y1; ++j)
            {
              grid.bz()(x0 - 1, j, k) += dtdx * incidentElectric(
                sources, origin, x0, j, k, 0.0, 0.5, 0.0,
                1, timeE, frame);
              grid.bz()(x1, j, k) -= dtdx * incidentElectric(
                sources, origin, x1, j, k, 0.0, 0.5, 0.0,
                1, timeE, frame);
            }
        for (std::size_t k = localNodeZ0_; k <= localNodeZ1_; ++k)
          for (std::size_t i = x0; i < x1; ++i)
            {
              grid.bz()(i, y0 - 1, k) -= dtdy * incidentElectric(
                sources, origin, i, y0, k, 0.5, 0.0, 0.0,
                0, timeE, frame);
              grid.bz()(i, y1, k) += dtdy * incidentElectric(
                sources, origin, i, y1, k, 0.5, 0.0, 0.0,
                0, timeE, frame);
            }
      }
  }

  void EBMaxwellIncidentInjector::correctAfterElectricUpdate(
      EBFieldGrid& grid, const SIFieldSourceSet& sources,
      const FieldVector<Double>& origin, Double timeE,
      const BoostFrameTransform& frame) const
  {
    validateGrid(grid);
    const std::size_t x0 = region_.lower[0];
    const std::size_t x1 = region_.upper[0];
    const std::size_t y0 = region_.lower[1];
    const std::size_t y1 = region_.upper[1];
    const Double c2dt = SI::c * SI::c * geometry_.dt;
    const Double c2dtdx = c2dt / geometry_.dx;
    const Double c2dtdy = c2dt / geometry_.dy;
    const Double c2dtdz = c2dt / geometry_.dz;
    const Double timeB = timeE + 0.5 * geometry_.dt;

    /* Ex: dBz/dy - dBy/dz. */
    if (hasNodeIntersection_)
      for (std::size_t k = localNodeZ0_; k <= localNodeZ1_; ++k)
        for (std::size_t i = x0; i < x1; ++i)
          {
            grid.ex()(i, y0, k) -= c2dtdy * incidentMagnetic(
              sources, origin, i, y0 - 1, k, 0.5, 0.5, 0.0,
              2, timeB, frame);
            grid.ex()(i, y1, k) += c2dtdy * incidentMagnetic(
              sources, origin, i, y1, k, 0.5, 0.5, 0.0,
              2, timeB, frame);
          }
    if (correctLowerZElectric_)
      {
        const std::size_t surface =
          region_.lower[2] - globalZCellOffset_;
        const std::ptrdiff_t outside =
          static_cast<std::ptrdiff_t>(surface) - 1;
        for (std::size_t j = y0; j <= y1; ++j)
          for (std::size_t i = x0; i < x1; ++i)
            grid.ex()(i, j, surface) += c2dtdz * incidentMagnetic(
              sources, origin, i, j, outside, 0.5, 0.0, 0.5,
              1, timeB, frame);
      }
    if (correctUpperZElectric_)
      {
        const std::size_t surface =
          region_.upper[2] - globalZCellOffset_;
        for (std::size_t j = y0; j <= y1; ++j)
          for (std::size_t i = x0; i < x1; ++i)
            grid.ex()(i, j, surface) -= c2dtdz * incidentMagnetic(
              sources, origin, i, j,
              static_cast<std::ptrdiff_t>(surface),
              0.5, 0.0, 0.5, 1, timeB, frame);
      }

    /* Ey: dBx/dz - dBz/dx. */
    if (correctLowerZElectric_)
      {
        const std::size_t surface =
          region_.lower[2] - globalZCellOffset_;
        const std::ptrdiff_t outside =
          static_cast<std::ptrdiff_t>(surface) - 1;
        for (std::size_t j = y0; j < y1; ++j)
          for (std::size_t i = x0; i <= x1; ++i)
            grid.ey()(i, j, surface) -= c2dtdz * incidentMagnetic(
              sources, origin, i, j, outside, 0.0, 0.5, 0.5,
              0, timeB, frame);
      }
    if (correctUpperZElectric_)
      {
        const std::size_t surface =
          region_.upper[2] - globalZCellOffset_;
        for (std::size_t j = y0; j < y1; ++j)
          for (std::size_t i = x0; i <= x1; ++i)
            grid.ey()(i, j, surface) += c2dtdz * incidentMagnetic(
              sources, origin, i, j,
              static_cast<std::ptrdiff_t>(surface),
              0.0, 0.5, 0.5, 0, timeB, frame);
      }
    if (hasNodeIntersection_)
      for (std::size_t k = localNodeZ0_; k <= localNodeZ1_; ++k)
        for (std::size_t j = y0; j < y1; ++j)
          {
            grid.ey()(x0, j, k) += c2dtdx * incidentMagnetic(
              sources, origin, x0 - 1, j, k, 0.5, 0.5, 0.0,
              2, timeB, frame);
            grid.ey()(x1, j, k) -= c2dtdx * incidentMagnetic(
              sources, origin, x1, j, k, 0.5, 0.5, 0.0,
              2, timeB, frame);
          }

    /* Ez: dBy/dx - dBx/dy. */
    if (hasCellIntersection_)
      {
        for (std::size_t k = localCellZ0_; k < localCellZ1_; ++k)
          for (std::size_t j = y0; j <= y1; ++j)
            {
              grid.ez()(x0, j, k) -= c2dtdx * incidentMagnetic(
                sources, origin, x0 - 1, j, k, 0.5, 0.0, 0.5,
                1, timeB, frame);
              grid.ez()(x1, j, k) += c2dtdx * incidentMagnetic(
                sources, origin, x1, j, k, 0.5, 0.0, 0.5,
                1, timeB, frame);
            }
        for (std::size_t k = localCellZ0_; k < localCellZ1_; ++k)
          for (std::size_t i = x0; i <= x1; ++i)
            {
              grid.ez()(i, y0, k) += c2dtdy * incidentMagnetic(
                sources, origin, i, y0 - 1, k, 0.0, 0.5, 0.5,
                0, timeB, frame);
              grid.ez()(i, y1, k) -= c2dtdy * incidentMagnetic(
                sources, origin, i, y1, k, 0.0, 0.5, 0.5,
                0, timeB, frame);
            }
      }
  }

  void EBMaxwellIncidentInjector::advance(
      EBFieldGrid& grid, const SIFieldSourceSet& sources,
      const FieldVector<Double>& origin, Double timeE,
      const BoostFrameTransform& frame) const
  {
    grid.advanceMagnetic();
    correctAfterMagneticUpdate(grid, sources, origin, timeE, frame);
    grid.advanceElectric();
    correctAfterElectricUpdate(grid, sources, origin, timeE, frame);
  }
}
