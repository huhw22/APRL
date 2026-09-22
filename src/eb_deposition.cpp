#include "eb_deposition.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fel
{
  namespace
  {
    Double shape(unsigned int upper, Double coordinate)
    {
      return upper ? coordinate : 1.0 - coordinate;
    }

    Double integratedProduct(Double a0, Double a1,
                             Double b0, Double b1)
    {
      const Double averageA = 0.5 * (a0 + a1);
      const Double averageB = 0.5 * (b0 + b1);
      return averageA * averageB + (a1 - a0) * (b1 - b0) / 12.0;
    }

    Double clampUnit(Double value)
    {
      return std::max(0.0, std::min(1.0, value));
    }
  }

  CurrentDepositResult::CurrentDepositResult()
    : subsegments(0), depositedCharge(0.0)
  {}

  ChargeConservingCurrentDepositor::ChargeConservingCurrentDepositor(
      EBFieldGrid& fields, const FieldVector<Double>& originSI)
    : fields_(fields), originSI_(originSI)
  {}

  CurrentDepositResult ChargeConservingCurrentDepositor::depositSegment(
      const FieldVector<Double>& startSI,
      const FieldVector<Double>& endSI,
      Double chargeCoulomb)
  {
    validatePosition(startSI);
    validatePosition(endSI);
    if (!std::isfinite(chargeCoulomb))
      throw std::invalid_argument("Deposited particle charge must be finite");

    const NormalizedPosition start = normalize(startSI);
    const NormalizedPosition end = normalize(endSI);
    const Double dx = end.x - start.x;
    const Double dy = end.y - start.y;
    const Double dz = end.z - start.z;

    CurrentDepositResult result;
    result.depositedCharge = chargeCoulomb;
    if (dx == 0.0 && dy == 0.0 && dz == 0.0)
      return result;

    std::size_t i = containingCell(start.x, dx, fields_.geometry().nx);
    std::size_t j = containingCell(start.y, dy, fields_.geometry().ny);
    std::size_t k = containingCell(start.z, dz, fields_.geometry().nz);

    const int stepX = (dx > 0.0) ? 1 : ((dx < 0.0) ? -1 : 0);
    const int stepY = (dy > 0.0) ? 1 : ((dy < 0.0) ? -1 : 0);
    const int stepZ = (dz > 0.0) ? 1 : ((dz < 0.0) ? -1 : 0);
    const Double infinity = std::numeric_limits<Double>::infinity();

    Double tMaxX = infinity;
    Double tMaxY = infinity;
    Double tMaxZ = infinity;
    Double tDeltaX = infinity;
    Double tDeltaY = infinity;
    Double tDeltaZ = infinity;

    if (stepX > 0)
      {
        tMaxX = (static_cast<Double>(i + 1) - start.x) / dx;
        tDeltaX = 1.0 / dx;
      }
    else if (stepX < 0)
      {
        tMaxX = (static_cast<Double>(i) - start.x) / dx;
        tDeltaX = -1.0 / dx;
      }

    if (stepY > 0)
      {
        tMaxY = (static_cast<Double>(j + 1) - start.y) / dy;
        tDeltaY = 1.0 / dy;
      }
    else if (stepY < 0)
      {
        tMaxY = (static_cast<Double>(j) - start.y) / dy;
        tDeltaY = -1.0 / dy;
      }

    if (stepZ > 0)
      {
        tMaxZ = (static_cast<Double>(k + 1) - start.z) / dz;
        tDeltaZ = 1.0 / dz;
      }
    else if (stepZ < 0)
      {
        tMaxZ = (static_cast<Double>(k) - start.z) / dz;
        tDeltaZ = -1.0 / dz;
      }

    Double t = 0.0;
    const std::size_t maximumSegments = fields_.geometry().nx +
                                        fields_.geometry().ny +
                                        fields_.geometry().nz + 3;
    const Double crossingTolerance = 64.0 *
      std::numeric_limits<Double>::epsilon();

    while (t < 1.0 && result.subsegments < maximumSegments)
      {
        const Double next = std::min(1.0,
          std::min(tMaxX, std::min(tMaxY, tMaxZ)));
        if (next < t - crossingTolerance)
          throw std::runtime_error("Invalid cell crossing in current deposition");

        if (next > t + crossingTolerance)
          {
            NormalizedPosition segmentStart;
            segmentStart.x = start.x + dx * t;
            segmentStart.y = start.y + dy * t;
            segmentStart.z = start.z + dz * t;
            NormalizedPosition segmentEnd;
            segmentEnd.x = start.x + dx * next;
            segmentEnd.y = start.y + dy * next;
            segmentEnd.z = start.z + dz * next;
            depositCellSegment(i, j, k, segmentStart, segmentEnd,
                               chargeCoulomb);
            ++result.subsegments;
          }

        if (next >= 1.0 - crossingTolerance)
          break;

        if (std::abs(tMaxX - next) <= crossingTolerance)
          {
            i = static_cast<std::size_t>(static_cast<long long>(i) + stepX);
            tMaxX += tDeltaX;
          }
        if (std::abs(tMaxY - next) <= crossingTolerance)
          {
            j = static_cast<std::size_t>(static_cast<long long>(j) + stepY);
            tMaxY += tDeltaY;
          }
        if (std::abs(tMaxZ - next) <= crossingTolerance)
          {
            k = static_cast<std::size_t>(static_cast<long long>(k) + stepZ);
            tMaxZ += tDeltaZ;
          }
        t = next;
      }

    if (t < 1.0 - crossingTolerance && result.subsegments >= maximumSegments)
      throw std::runtime_error("Current deposition exceeded crossing bound");
    return result;
  }

  void ChargeConservingCurrentDepositor::depositCharge(
      const FieldVector<Double>& positionSI, Double chargeCoulomb,
      YeeComponent& rho) const
  {
    validatePosition(positionSI);
    const EBGridGeometry& g = fields_.geometry();
    if (rho.nx() != g.nx + 1 || rho.ny() != g.ny + 1 ||
        rho.nz() != g.nz + 1)
      throw std::invalid_argument(
          "Charge lattice must have (nx+1, ny+1, nz+1) vertices");

    const NormalizedPosition p = normalize(positionSI);
    const std::size_t i = containingCell(p.x, 0.0, g.nx);
    const std::size_t j = containingCell(p.y, 0.0, g.ny);
    const std::size_t k = containingCell(p.z, 0.0, g.nz);
    const Double x = clampUnit(p.x - static_cast<Double>(i));
    const Double y = clampUnit(p.y - static_cast<Double>(j));
    const Double z = clampUnit(p.z - static_cast<Double>(k));
    const Double inverseVolume = 1.0 / (g.dx * g.dy * g.dz);

    for (unsigned int c = 0; c < 2; ++c)
      for (unsigned int b = 0; b < 2; ++b)
        for (unsigned int a = 0; a < 2; ++a)
          rho(i + a, j + b, k + c) += chargeCoulomb * inverseVolume *
            shape(a, x) * shape(b, y) * shape(c, z);
  }

  Double ChargeConservingCurrentDepositor::maxContinuityResidual(
      const EBFieldGrid& fields,
      const YeeComponent& rhoBefore,
      const YeeComponent& rhoAfter)
  {
    const EBGridGeometry& g = fields.geometry();
    if (rhoBefore.nx() != g.nx + 1 || rhoBefore.ny() != g.ny + 1 ||
        rhoBefore.nz() != g.nz + 1 ||
        rhoAfter.nx() != g.nx + 1 || rhoAfter.ny() != g.ny + 1 ||
        rhoAfter.nz() != g.nz + 1)
      throw std::invalid_argument("Continuity diagnostic got invalid rho grid");

    Double maximum = 0.0;
    for (std::size_t k = 0; k <= g.nz; ++k)
      for (std::size_t j = 0; j <= g.ny; ++j)
        for (std::size_t i = 0; i <= g.nx; ++i)
          {
            const Double rightX = (i < g.nx) ? fields.jx()(i, j, k) : 0.0;
            const Double leftX = (i > 0) ? fields.jx()(i - 1, j, k) : 0.0;
            const Double rightY = (j < g.ny) ? fields.jy()(i, j, k) : 0.0;
            const Double leftY = (j > 0) ? fields.jy()(i, j - 1, k) : 0.0;
            const Double rightZ = (k < g.nz) ? fields.jz()(i, j, k) : 0.0;
            const Double leftZ = (k > 0) ? fields.jz()(i, j, k - 1) : 0.0;
            const Double residual =
              (rhoAfter(i, j, k) - rhoBefore(i, j, k)) / g.dt +
              (rightX - leftX) / g.dx +
              (rightY - leftY) / g.dy +
              (rightZ - leftZ) / g.dz;
            maximum = std::max(maximum, std::abs(residual));
          }
    return maximum;
  }

  void ChargeConservingCurrentDepositor::validatePosition(
      const FieldVector<Double>& positionSI) const
  {
    const EBGridGeometry& g = fields_.geometry();
    const Double extent[3] = {
      static_cast<Double>(g.nx) * g.dx,
      static_cast<Double>(g.ny) * g.dy,
      static_cast<Double>(g.nz) * g.dz
    };
    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        const Double relative = positionSI[axis] - originSI_[axis];
        const Double tolerance = 64.0 * std::numeric_limits<Double>::epsilon() *
          std::max(1.0, std::max(std::abs(originSI_[axis]), extent[axis]));
        if (!std::isfinite(positionSI[axis]) || relative < -tolerance ||
            relative > extent[axis] + tolerance)
          throw std::out_of_range(
              "Particle trajectory endpoint is outside the E/B domain");
      }
  }

  ChargeConservingCurrentDepositor::NormalizedPosition
  ChargeConservingCurrentDepositor::normalize(
      const FieldVector<Double>& positionSI) const
  {
    const EBGridGeometry& g = fields_.geometry();
    NormalizedPosition position;
    position.x = std::max(0.0, std::min(static_cast<Double>(g.nx),
      (positionSI[0] - originSI_[0]) / g.dx));
    position.y = std::max(0.0, std::min(static_cast<Double>(g.ny),
      (positionSI[1] - originSI_[1]) / g.dy));
    position.z = std::max(0.0, std::min(static_cast<Double>(g.nz),
      (positionSI[2] - originSI_[2]) / g.dz));
    return position;
  }

  std::size_t ChargeConservingCurrentDepositor::containingCell(
      Double coordinate, Double direction, std::size_t cells) const
  {
    if (coordinate <= 0.0)
      return 0;
    if (coordinate >= static_cast<Double>(cells))
      return cells - 1;

    const Double nearest = std::floor(coordinate);
    const Double tolerance = 64.0 *
      std::numeric_limits<Double>::epsilon() *
      std::max(1.0, std::abs(coordinate));
    if (direction < 0.0 && std::abs(coordinate - nearest) <= tolerance)
      return static_cast<std::size_t>(nearest) - 1;
    return static_cast<std::size_t>(nearest);
  }

  void ChargeConservingCurrentDepositor::depositCellSegment(
      std::size_t i, std::size_t j, std::size_t k,
      const NormalizedPosition& start,
      const NormalizedPosition& end,
      Double chargeCoulomb)
  {
    const EBGridGeometry& g = fields_.geometry();
    if (i >= g.nx || j >= g.ny || k >= g.nz)
      throw std::out_of_range("Current subsegment cell is outside the grid");

    const Double x0 = clampUnit(start.x - static_cast<Double>(i));
    const Double y0 = clampUnit(start.y - static_cast<Double>(j));
    const Double z0 = clampUnit(start.z - static_cast<Double>(k));
    const Double x1 = clampUnit(end.x - static_cast<Double>(i));
    const Double y1 = clampUnit(end.y - static_cast<Double>(j));
    const Double z1 = clampUnit(end.z - static_cast<Double>(k));
    const Double inverseVolumeTime = 1.0 /
      (g.dx * g.dy * g.dz * g.dt);

    const Double currentX = chargeCoulomb * (x1 - x0) * g.dx *
                            inverseVolumeTime;
    const Double currentY = chargeCoulomb * (y1 - y0) * g.dy *
                            inverseVolumeTime;
    const Double currentZ = chargeCoulomb * (z1 - z0) * g.dz *
                            inverseVolumeTime;

    for (unsigned int c = 0; c < 2; ++c)
      for (unsigned int b = 0; b < 2; ++b)
        fields_.jx()(i, j + b, k + c) += currentX * integratedProduct(
          shape(b, y0), shape(b, y1), shape(c, z0), shape(c, z1));

    for (unsigned int c = 0; c < 2; ++c)
      for (unsigned int a = 0; a < 2; ++a)
        fields_.jy()(i + a, j, k + c) += currentY * integratedProduct(
          shape(a, x0), shape(a, x1), shape(c, z0), shape(c, z1));

    for (unsigned int b = 0; b < 2; ++b)
      for (unsigned int a = 0; a < 2; ++a)
        fields_.jz()(i + a, j + b, k) += currentZ * integratedProduct(
          shape(a, x0), shape(a, x1), shape(b, y0), shape(b, y1));
  }
}
