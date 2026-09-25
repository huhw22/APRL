#include "eb_field.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace aprl
{
  namespace SI
  {
    const Double c = 299792458.0;
    const Double epsilon0 = 8.8541878128e-12;
    const Double mu0 = 1.25663706212e-6;
    const Double elementaryCharge = 1.602176634e-19;
    const Double electronMass = 9.1093837139e-31;
  }

  EBCowanCoefficients::EBCowanCoefficients()
    : deltaXY(0.0), deltaYZ(0.0), deltaZX(0.0)
  {
    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        alpha[axis] = 1.0;
        beta[axis] = 0.0;
        ratio[axis] = 1.0;
      }
  }

  EBCowanCoefficients EBCowanCoefficients::forZDispersion(
      Double dx, Double dy, Double dz)
  {
    if (!(dx > 0.0) || !(dy > 0.0) || !(dz > 0.0) ||
        !std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz))
      throw std::invalid_argument(
        "Cowan-z grid spacing must be positive and finite");
    const Double tolerance = 64.0 * std::numeric_limits<Double>::epsilon() *
      std::max(dx, std::max(dy, dz));
    if (dx + tolerance < dz || dy + tolerance < dz)
      throw std::invalid_argument(
        "Cowan-z requires dz to be the smallest grid spacing");

    EBCowanCoefficients result;
    result.ratio[0] = (dz / dx) * (dz / dx);
    result.ratio[1] = (dz / dy) * (dz / dy);
    result.ratio[2] = 1.0;
    const Double rx = result.ratio[0];
    const Double ry = result.ratio[1];
    const Double rz = result.ratio[2];
    const Double denominator = rx * ry + ry * rz + rz * rx;
    if (!(denominator > 0.0) || !std::isfinite(denominator))
      throw std::invalid_argument("Invalid Cowan-z aspect ratios");
    const Double productFraction = rx * ry * rz / denominator;
    result.beta[0] = rx * (1.0 - productFraction) / 8.0;
    result.beta[1] = ry * (1.0 - productFraction) / 8.0;
    result.beta[2] = rz * (1.0 - productFraction) / 8.0;
    result.deltaXY = rx * ry *
      (1.0 / 16.0 - rx * ry / (8.0 * denominator));
    result.deltaYZ = ry * rz *
      (1.0 / 16.0 - ry * rz / (8.0 * denominator));
    result.deltaZX = rz * rx *
      (1.0 / 16.0 - rz * rx / (8.0 * denominator));
    result.alpha[0] = 1.0 - 2.0 * result.beta[1] -
      2.0 * result.beta[2] - 4.0 * result.deltaYZ;
    result.alpha[1] = 1.0 - 2.0 * result.beta[2] -
      2.0 * result.beta[0] - 4.0 * result.deltaZX;
    result.alpha[2] = 1.0 - 2.0 * result.beta[0] -
      2.0 * result.beta[1] - 4.0 * result.deltaXY;
    return result;
  }

  EBGridGeometry::EBGridGeometry()
    : nx(0), ny(0), nz(0), dx(0.0), dy(0.0), dz(0.0), dt(0.0),
      solver(EBMaxwellSolver::Yee)
  {}

  EBGridGeometry::EBGridGeometry(std::size_t nxValue,
                                 std::size_t nyValue,
                                 std::size_t nzValue,
                                 Double dxValue, Double dyValue,
                                 Double dzValue, Double dtValue,
                                 EBMaxwellSolver solverValue)
    : nx(nxValue), ny(nyValue), nz(nzValue),
      dx(dxValue), dy(dyValue), dz(dzValue), dt(dtValue),
      solver(solverValue)
  {}

  EBElectricHaloView::EBElectricHaloView()
    : lowerEx(0), lowerEy(0), lowerEz(0),
      upperEx(0), upperEy(0), upperEz(0)
  {}

  YeeComponent::YeeComponent() : nx_(0), ny_(0), nz_(0), values_()
  {}

  YeeComponent::YeeComponent(std::size_t nx, std::size_t ny, std::size_t nz)
    : nx_(0), ny_(0), nz_(0), values_()
  {
    resize(nx, ny, nz);
  }

  void YeeComponent::resize(std::size_t nx, std::size_t ny, std::size_t nz)
  {
    if (nx != 0 && ny > std::numeric_limits<std::size_t>::max() / nx)
      throw std::overflow_error("Yee component size overflows size_t");
    const std::size_t xy = nx * ny;
    if (xy != 0 && nz > std::numeric_limits<std::size_t>::max() / xy)
      throw std::overflow_error("Yee component size overflows size_t");

    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
    values_.assign(xy * nz, 0.0);
  }

  void YeeComponent::fill(Double value)
  {
    std::fill(values_.begin(), values_.end(), value);
  }

  std::size_t YeeComponent::index(std::size_t i, std::size_t j,
                                  std::size_t k) const
  {
#ifdef DIRECT_EB_BOUNDS_CHECK
    if (i >= nx_ || j >= ny_ || k >= nz_)
      throw std::out_of_range("Yee component index is outside the lattice");
#endif
    return (k * ny_ + j) * nx_ + i;
  }

  Double& YeeComponent::operator()(std::size_t i, std::size_t j,
                                   std::size_t k)
  {
    return values_[index(i, j, k)];
  }

  const Double& YeeComponent::operator()(std::size_t i, std::size_t j,
                                         std::size_t k) const
  {
    return values_[index(i, j, k)];
  }

  std::size_t YeeComponent::nx() const { return nx_; }
  std::size_t YeeComponent::ny() const { return ny_; }
  std::size_t YeeComponent::nz() const { return nz_; }
  std::size_t YeeComponent::size() const { return values_.size(); }
  std::size_t YeeComponent::bytes() const
  {
    return values_.capacity() * sizeof(Double);
  }
  Double* YeeComponent::data() { return values_.empty() ? 0 : &values_[0]; }
  const Double* YeeComponent::data() const
  {
    return values_.empty() ? 0 : &values_[0];
  }

  RadiationFieldSample::RadiationFieldSample()
    : electric(0.0), magnetic(0.0), poynting(0.0), energyDensity(0.0)
  {}

  EBMemoryFootprint::EBMemoryFootprint()
    : electricBytes(0), magneticBytes(0), currentBytes(0), boundaryBytes(0)
  {}

  std::size_t EBMemoryFootprint::totalBytes() const
  {
    return electricBytes + magneticBytes + currentBytes + boundaryBytes;
  }

  EBBoundaryOperator::~EBBoundaryOperator() {}

  PerfectElectricConductorBoundary::PerfectElectricConductorBoundary(
      bool lowerZPhysical, bool upperZPhysical)
    : lowerZPhysical_(lowerZPhysical), upperZPhysical_(upperZPhysical)
  {}

  void PerfectElectricConductorBoundary::afterMagneticUpdate(
      EBFieldGrid&, const EBElectricHaloView&)
  {}

  void PerfectElectricConductorBoundary::afterElectricUpdate(
      EBFieldGrid& fields)
  {
    const EBGridGeometry& g = fields.geometry();

    /* Tangential electric field is zero on each domain face. */
    for (std::size_t k = 0; k <= g.nz; ++k)
      for (std::size_t i = 0; i < g.nx; ++i)
        {
          fields.ex()(i, 0, k) = 0.0;
          fields.ex()(i, g.ny, k) = 0.0;
        }
    if (lowerZPhysical_ || upperZPhysical_)
      for (std::size_t j = 0; j <= g.ny; ++j)
        for (std::size_t i = 0; i < g.nx; ++i)
          {
            if (lowerZPhysical_) fields.ex()(i, j, 0) = 0.0;
            if (upperZPhysical_) fields.ex()(i, j, g.nz) = 0.0;
          }

    for (std::size_t k = 0; k <= g.nz; ++k)
      for (std::size_t j = 0; j < g.ny; ++j)
        {
          fields.ey()(0, j, k) = 0.0;
          fields.ey()(g.nx, j, k) = 0.0;
        }
    if (lowerZPhysical_ || upperZPhysical_)
      for (std::size_t j = 0; j < g.ny; ++j)
        for (std::size_t i = 0; i <= g.nx; ++i)
          {
            if (lowerZPhysical_) fields.ey()(i, j, 0) = 0.0;
            if (upperZPhysical_) fields.ey()(i, j, g.nz) = 0.0;
          }

    for (std::size_t k = 0; k < g.nz; ++k)
      for (std::size_t j = 0; j <= g.ny; ++j)
        {
          fields.ez()(0, j, k) = 0.0;
          fields.ez()(g.nx, j, k) = 0.0;
        }
    for (std::size_t k = 0; k < g.nz; ++k)
      for (std::size_t i = 0; i <= g.nx; ++i)
        {
          fields.ez()(i, 0, k) = 0.0;
          fields.ez()(i, g.ny, k) = 0.0;
        }
  }

  std::size_t PerfectElectricConductorBoundary::memoryBytes() const
  {
    return 0;
  }

  EBFieldGrid::EBFieldGrid(const EBGridGeometry& geometry)
    : geometry_(geometry),
      cowan_(geometry.solver == EBMaxwellSolver::CowanZ ?
        EBCowanCoefficients::forZDispersion(
          geometry.dx, geometry.dy, geometry.dz) :
        EBCowanCoefficients()),
      ex_(geometry.nx, geometry.ny + 1, geometry.nz + 1),
      ey_(geometry.nx + 1, geometry.ny, geometry.nz + 1),
      ez_(geometry.nx + 1, geometry.ny + 1, geometry.nz),
      bx_(geometry.nx + 1, geometry.ny, geometry.nz),
      by_(geometry.nx, geometry.ny + 1, geometry.nz),
      bz_(geometry.nx, geometry.ny, geometry.nz + 1),
      jx_(geometry.nx, geometry.ny + 1, geometry.nz + 1),
      jy_(geometry.nx + 1, geometry.ny, geometry.nz + 1),
      jz_(geometry.nx + 1, geometry.ny + 1, geometry.nz),
      boundary_(new PerfectElectricConductorBoundary())
  {
    validateGeometry();
  }

  const EBGridGeometry& EBFieldGrid::geometry() const { return geometry_; }

  YeeComponent& EBFieldGrid::ex() { return ex_; }
  YeeComponent& EBFieldGrid::ey() { return ey_; }
  YeeComponent& EBFieldGrid::ez() { return ez_; }
  YeeComponent& EBFieldGrid::bx() { return bx_; }
  YeeComponent& EBFieldGrid::by() { return by_; }
  YeeComponent& EBFieldGrid::bz() { return bz_; }
  YeeComponent& EBFieldGrid::jx() { return jx_; }
  YeeComponent& EBFieldGrid::jy() { return jy_; }
  YeeComponent& EBFieldGrid::jz() { return jz_; }

  const YeeComponent& EBFieldGrid::ex() const { return ex_; }
  const YeeComponent& EBFieldGrid::ey() const { return ey_; }
  const YeeComponent& EBFieldGrid::ez() const { return ez_; }
  const YeeComponent& EBFieldGrid::bx() const { return bx_; }
  const YeeComponent& EBFieldGrid::by() const { return by_; }
  const YeeComponent& EBFieldGrid::bz() const { return bz_; }
  const YeeComponent& EBFieldGrid::jx() const { return jx_; }
  const YeeComponent& EBFieldGrid::jy() const { return jy_; }
  const YeeComponent& EBFieldGrid::jz() const { return jz_; }

  void EBFieldGrid::setBoundary(
      std::unique_ptr<EBBoundaryOperator> boundary)
  {
    if (!boundary)
      throw std::invalid_argument("E/B boundary operator cannot be null");
    boundary_ = std::move(boundary);
  }

  void EBFieldGrid::clearFields()
  {
    ex_.fill(0.0); ey_.fill(0.0); ez_.fill(0.0);
    bx_.fill(0.0); by_.fill(0.0); bz_.fill(0.0);
  }

  void EBFieldGrid::clearCurrent()
  {
    jx_.fill(0.0); jy_.fill(0.0); jz_.fill(0.0);
  }

  void EBFieldGrid::advanceMagnetic()
  {
    advanceMagnetic(EBElectricHaloView());
  }

  void EBFieldGrid::advanceMagnetic(const EBElectricHaloView& halo)
  {
    if (geometry_.solver == EBMaxwellSolver::CowanZ)
      advanceMagneticCowan(halo);
    else
      advanceMagneticYee();
    boundary_->afterMagneticUpdate(*this, halo);
  }

  void EBFieldGrid::advanceMagneticYee()
  {
    const Double dtdx = geometry_.dt / geometry_.dx;
    const Double dtdy = geometry_.dt / geometry_.dy;
    const Double dtdz = geometry_.dt / geometry_.dz;

    for (std::size_t k = 0; k < geometry_.nz; ++k)
      for (std::size_t j = 0; j < geometry_.ny; ++j)
        for (std::size_t i = 0; i <= geometry_.nx; ++i)
          bx_(i, j, k) -= dtdy * (ez_(i, j + 1, k) - ez_(i, j, k))
                        - dtdz * (ey_(i, j, k + 1) - ey_(i, j, k));

    for (std::size_t k = 0; k < geometry_.nz; ++k)
      for (std::size_t j = 0; j <= geometry_.ny; ++j)
        for (std::size_t i = 0; i < geometry_.nx; ++i)
          by_(i, j, k) -= dtdz * (ex_(i, j, k + 1) - ex_(i, j, k))
                        - dtdx * (ez_(i + 1, j, k) - ez_(i, j, k));

    for (std::size_t k = 0; k <= geometry_.nz; ++k)
      for (std::size_t j = 0; j < geometry_.ny; ++j)
        for (std::size_t i = 0; i < geometry_.nx; ++i)
          bz_(i, j, k) -= dtdx * (ey_(i + 1, j, k) - ey_(i, j, k))
                        - dtdy * (ex_(i, j + 1, k) - ex_(i, j, k));

  }

  Double EBFieldGrid::electricValue(
      const YeeComponent& component, unsigned int componentAxis,
      std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k,
      const EBElectricHaloView& halo) const
  {
    std::ptrdiff_t coordinate[3] = {i, j, k};
    const std::ptrdiff_t size[3] = {
      static_cast<std::ptrdiff_t>(component.nx()),
      static_cast<std::ptrdiff_t>(component.ny()),
      static_cast<std::ptrdiff_t>(component.nz())
    };
    Double sign = 1.0;

    for (unsigned int axis = 0; axis < 2; ++axis)
      {
        if (coordinate[axis] < 0)
          {
            if (componentAxis == axis)
              coordinate[axis] = 0;
            else
              {
                coordinate[axis] = -coordinate[axis];
                sign = -sign;
              }
          }
        else if (coordinate[axis] >= size[axis])
          {
            if (componentAxis == axis)
              coordinate[axis] = size[axis] - 1;
            else
              {
                coordinate[axis] = 2 * (size[axis] - 1) -
                                   coordinate[axis];
                sign = -sign;
              }
          }
      }

    if (coordinate[2] < 0 || coordinate[2] >= size[2])
      {
        const bool lower = coordinate[2] < 0;
        const Double* plane = 0;
        if (&component == &ex_)
          plane = lower ? halo.lowerEx : halo.upperEx;
        else if (&component == &ey_)
          plane = lower ? halo.lowerEy : halo.upperEy;
        else if (&component == &ez_)
          plane = lower ? halo.lowerEz : halo.upperEz;
        if (plane)
          return sign * plane[
            static_cast<std::size_t>(coordinate[1]) * component.nx() +
            static_cast<std::size_t>(coordinate[0])];

        if (componentAxis == 2)
          coordinate[2] = lower ? 0 : size[2] - 1;
        else
          {
            coordinate[2] = lower ? -coordinate[2] :
              2 * (size[2] - 1) - coordinate[2];
            sign = -sign;
          }
      }

    return sign * component(
      static_cast<std::size_t>(coordinate[0]),
      static_cast<std::size_t>(coordinate[1]),
      static_cast<std::size_t>(coordinate[2]));
  }

  Double EBFieldGrid::smoothedElectric(
      const YeeComponent& component, unsigned int componentAxis,
      unsigned int derivativeAxis, std::ptrdiff_t i,
      std::ptrdiff_t j, std::ptrdiff_t k,
      const EBElectricHaloView& halo) const
  {
    const unsigned int transverseA = (derivativeAxis + 1) % 3;
    const unsigned int transverseB = (derivativeAxis + 2) % 3;
    const Double delta = derivativeAxis == 0 ? cowan_.deltaYZ :
      (derivativeAxis == 1 ? cowan_.deltaZX : cowan_.deltaXY);
    std::ptrdiff_t base[3] = {i, j, k};
    Double result = cowan_.alpha[derivativeAxis] *
      electricValue(component, componentAxis, i, j, k, halo);
    for (int direction = -1; direction <= 1; direction += 2)
      {
        std::ptrdiff_t pointA[3] = {base[0], base[1], base[2]};
        std::ptrdiff_t pointB[3] = {base[0], base[1], base[2]};
        pointA[transverseA] += direction;
        pointB[transverseB] += direction;
        result += cowan_.beta[transverseA] * electricValue(
          component, componentAxis, pointA[0], pointA[1], pointA[2], halo);
        result += cowan_.beta[transverseB] * electricValue(
          component, componentAxis, pointB[0], pointB[1], pointB[2], halo);
      }
    for (int directionA = -1; directionA <= 1; directionA += 2)
      for (int directionB = -1; directionB <= 1; directionB += 2)
        {
          std::ptrdiff_t point[3] = {base[0], base[1], base[2]};
          point[transverseA] += directionA;
          point[transverseB] += directionB;
          result += delta * electricValue(
            component, componentAxis, point[0], point[1], point[2], halo);
        }
    return result;
  }

  void EBFieldGrid::advanceMagneticCowan(const EBElectricHaloView& halo)
  {
    const Double dtdx = geometry_.dt / geometry_.dx;
    const Double dtdy = geometry_.dt / geometry_.dy;
    const Double dtdz = geometry_.dt / geometry_.dz;

    for (std::size_t k = 0; k < geometry_.nz; ++k)
      for (std::size_t j = 0; j < geometry_.ny; ++j)
        for (std::size_t i = 0; i <= geometry_.nx; ++i)
          {
            const Double dEzDy = smoothedElectric(
              ez_, 2, 1, i, j + 1, k, halo) -
              smoothedElectric(ez_, 2, 1, i, j, k, halo);
            const Double dEyDz = smoothedElectric(
              ey_, 1, 2, i, j, k + 1, halo) -
              smoothedElectric(ey_, 1, 2, i, j, k, halo);
            bx_(i, j, k) -= dtdy * dEzDy - dtdz * dEyDz;
          }

    for (std::size_t k = 0; k < geometry_.nz; ++k)
      for (std::size_t j = 0; j <= geometry_.ny; ++j)
        for (std::size_t i = 0; i < geometry_.nx; ++i)
          {
            const Double dExDz = smoothedElectric(
              ex_, 0, 2, i, j, k + 1, halo) -
              smoothedElectric(ex_, 0, 2, i, j, k, halo);
            const Double dEzDx = smoothedElectric(
              ez_, 2, 0, i + 1, j, k, halo) -
              smoothedElectric(ez_, 2, 0, i, j, k, halo);
            by_(i, j, k) -= dtdz * dExDz - dtdx * dEzDx;
          }

    for (std::size_t k = 0; k <= geometry_.nz; ++k)
      for (std::size_t j = 0; j < geometry_.ny; ++j)
        for (std::size_t i = 0; i < geometry_.nx; ++i)
          {
            const Double dEyDx = smoothedElectric(
              ey_, 1, 0, i + 1, j, k, halo) -
              smoothedElectric(ey_, 1, 0, i, j, k, halo);
            const Double dExDy = smoothedElectric(
              ex_, 0, 1, i, j + 1, k, halo) -
              smoothedElectric(ex_, 0, 1, i, j, k, halo);
            bz_(i, j, k) -= dtdx * dEyDx - dtdy * dExDy;
          }
  }

  void EBFieldGrid::advanceElectric()
  {
    const Double curlScale = geometry_.dt / (SI::mu0 * SI::epsilon0);
    const Double currentScale = geometry_.dt / SI::epsilon0;

    const Double curlX = curlScale / geometry_.dx;
    const Double curlY = curlScale / geometry_.dy;
    const Double curlZ = curlScale / geometry_.dz;

    for (std::size_t k = 1; k < geometry_.nz; ++k)
      for (std::size_t j = 1; j < geometry_.ny; ++j)
        for (std::size_t i = 0; i < geometry_.nx; ++i)
          ex_(i, j, k) += curlY * (bz_(i, j, k) - bz_(i, j - 1, k))
                        - curlZ * (by_(i, j, k) - by_(i, j, k - 1))
                        - currentScale * jx_(i, j, k);

    for (std::size_t k = 1; k < geometry_.nz; ++k)
      for (std::size_t j = 0; j < geometry_.ny; ++j)
        for (std::size_t i = 1; i < geometry_.nx; ++i)
          ey_(i, j, k) += curlZ * (bx_(i, j, k) - bx_(i, j, k - 1))
                        - curlX * (bz_(i, j, k) - bz_(i - 1, j, k))
                        - currentScale * jy_(i, j, k);

    for (std::size_t k = 0; k < geometry_.nz; ++k)
      for (std::size_t j = 1; j < geometry_.ny; ++j)
        for (std::size_t i = 1; i < geometry_.nx; ++i)
          ez_(i, j, k) += curlX * (by_(i, j, k) - by_(i - 1, j, k))
                        - curlY * (bx_(i, j, k) - bx_(i, j - 1, k))
                        - currentScale * jz_(i, j, k);

    boundary_->afterElectricUpdate(*this);
  }

  void EBFieldGrid::advance()
  {
    advanceMagnetic();
    advanceElectric();
  }

  Double EBFieldGrid::courantNumber() const
  {
    return SI::c * geometry_.dt *
      std::sqrt(1.0 / (geometry_.dx * geometry_.dx) +
                1.0 / (geometry_.dy * geometry_.dy) +
                1.0 / (geometry_.dz * geometry_.dz));
  }

  Double EBFieldGrid::courantLimit() const
  {
    return 1.0 / (SI::c *
      std::sqrt(1.0 / (geometry_.dx * geometry_.dx) +
                1.0 / (geometry_.dy * geometry_.dy) +
                1.0 / (geometry_.dz * geometry_.dz)));
  }

  const EBCowanCoefficients& EBFieldGrid::cowanCoefficients() const
  {
    return cowan_;
  }

  Double EBFieldGrid::axisPhaseVelocityRatio(
      Double courantAxis, Double cellsPerWavelength)
  {
    if (!(courantAxis > 0.0) || !(courantAxis <= 1.0) ||
        !(cellsPerWavelength > 2.0))
      throw std::invalid_argument(
        "Dispersion diagnostic requires 0 < Courant <= 1 and more than two cells per wavelength");
    const Double pi = 3.1415926535897932384626433832795;
    const Double angle = pi / cellsPerWavelength;
    return cellsPerWavelength /
      (pi * courantAxis) * std::asin(courantAxis * std::sin(angle));
  }

  Double EBFieldGrid::axisGroupVelocityRatio(
      Double courantAxis, Double cellsPerWavelength)
  {
    if (!(courantAxis > 0.0) || !(courantAxis <= 1.0) ||
        !(cellsPerWavelength > 2.0))
      throw std::invalid_argument(
        "Dispersion diagnostic requires 0 < Courant <= 1 and more than two cells per wavelength");
    const Double pi = 3.1415926535897932384626433832795;
    const Double angle = pi / cellsPerWavelength;
    const Double sine = std::sin(angle);
    return std::cos(angle) /
      std::sqrt(1.0 - courantAxis * courantAxis * sine * sine);
  }

  EBMemoryFootprint EBFieldGrid::memoryFootprint() const
  {
    EBMemoryFootprint footprint;
    footprint.electricBytes = ex_.bytes() + ey_.bytes() + ez_.bytes();
    footprint.magneticBytes = bx_.bytes() + by_.bytes() + bz_.bytes();
    footprint.currentBytes = jx_.bytes() + jy_.bytes() + jz_.bytes();
    footprint.boundaryBytes = boundary_->memoryBytes();
    return footprint;
  }

  RadiationFieldSample EBFieldGrid::radiationSampleCell(
      std::size_t i, std::size_t j, std::size_t k) const
  {
    if (i >= geometry_.nx || j >= geometry_.ny || k >= geometry_.nz)
      throw std::out_of_range("Radiation sample is outside the E/B grid");

    FieldVector<Double> electric(0.0);
    electric[0] = 0.25 * (ex_(i, j, k) + ex_(i, j + 1, k) +
                          ex_(i, j, k + 1) + ex_(i, j + 1, k + 1));
    electric[1] = 0.25 * (ey_(i, j, k) + ey_(i + 1, j, k) +
                          ey_(i, j, k + 1) + ey_(i + 1, j, k + 1));
    electric[2] = 0.25 * (ez_(i, j, k) + ez_(i + 1, j, k) +
                          ez_(i, j + 1, k) + ez_(i + 1, j + 1, k));

    FieldVector<Double> magnetic(0.0);
    magnetic[0] = 0.5 * (bx_(i, j, k) + bx_(i + 1, j, k));
    magnetic[1] = 0.5 * (by_(i, j, k) + by_(i, j + 1, k));
    magnetic[2] = 0.5 * (bz_(i, j, k) + bz_(i, j, k + 1));

    return makeRadiationSample(electric, magnetic);
  }

  RadiationFieldSample EBFieldGrid::radiationSampleCellLab(
      std::size_t i, std::size_t j, std::size_t k,
      const BoostFrameTransform& frame) const
  {
    const RadiationFieldSample box = radiationSampleCell(i, j, k);
    FieldVector<Double> electricLab(0.0);
    FieldVector<Double> magneticLab(0.0);
    frame.fieldsBoxToLab(box.electric, box.magnetic,
                         electricLab, magneticLab);
    return makeRadiationSample(electricLab, magneticLab);
  }

  void EBFieldGrid::sampleFieldsPosition(
      const FieldVector<Double>& positionSI,
      const FieldVector<Double>& originSI,
      FieldVector<Double>& electric,
      FieldVector<Double>& magnetic) const
  {
    electric[0] = sampleComponent(ex_, positionSI, originSI, 0.5, 0.0, 0.0);
    electric[1] = sampleComponent(ey_, positionSI, originSI, 0.0, 0.5, 0.0);
    electric[2] = sampleComponent(ez_, positionSI, originSI, 0.0, 0.0, 0.5);
    magnetic[0] = sampleComponent(bx_, positionSI, originSI, 0.0, 0.5, 0.5);
    magnetic[1] = sampleComponent(by_, positionSI, originSI, 0.5, 0.0, 0.5);
    magnetic[2] = sampleComponent(bz_, positionSI, originSI, 0.5, 0.5, 0.0);
  }

  RadiationFieldSample EBFieldGrid::radiationSamplePositionLab(
      const FieldVector<Double>& positionBoxSI,
      const FieldVector<Double>& originBoxSI,
      const BoostFrameTransform& frame) const
  {
    FieldVector<Double> electricBox(0.0);
    FieldVector<Double> magneticBox(0.0);
    FieldVector<Double> electricLab(0.0);
    FieldVector<Double> magneticLab(0.0);
    sampleFieldsPosition(positionBoxSI, originBoxSI,
                         electricBox, magneticBox);
    frame.fieldsBoxToLab(electricBox, magneticBox,
                         electricLab, magneticLab);
    return makeRadiationSample(electricLab, magneticLab);
  }

  void EBFieldGrid::validateGeometry() const
  {
    if (geometry_.nx < 2 || geometry_.ny < 2 || geometry_.nz < 2)
      throw std::invalid_argument(
          "E/B staggered grid needs at least two cells in every direction");
    if (!(geometry_.dx > 0.0) || !(geometry_.dy > 0.0) ||
        !(geometry_.dz > 0.0) || !(geometry_.dt > 0.0) ||
        !std::isfinite(geometry_.dx) || !std::isfinite(geometry_.dy) ||
        !std::isfinite(geometry_.dz) || !std::isfinite(geometry_.dt))
      throw std::invalid_argument(
          "E/B grid spacing and time step must be positive finite SI values");
    if (geometry_.solver == EBMaxwellSolver::Yee)
      {
        if (courantNumber() >
            1.0 + 16.0 * std::numeric_limits<Double>::epsilon())
          throw std::invalid_argument(
            "E/B time step violates the three-dimensional Yee CFL limit");
      }
    else
      {
        const Double spacingTolerance = 64.0 *
          std::numeric_limits<Double>::epsilon() *
          std::max(geometry_.dx, std::max(geometry_.dy, geometry_.dz));
        if (geometry_.dx + spacingTolerance < geometry_.dz ||
            geometry_.dy + spacingTolerance < geometry_.dz)
          throw std::invalid_argument(
            "Cowan-z requires dz to be the smallest grid spacing");
        const Double axialCourant = SI::c * geometry_.dt / geometry_.dz;
        if (std::abs(axialCourant - 1.0) >
            64.0 * std::numeric_limits<Double>::epsilon())
          throw std::invalid_argument(
            "Cowan-z requires the dispersion-matched step c*dt=dz");
        const Double normalization[3] = {
          cowan_.alpha[0] + 2.0 * cowan_.beta[1] +
            2.0 * cowan_.beta[2] + 4.0 * cowan_.deltaYZ,
          cowan_.alpha[1] + 2.0 * cowan_.beta[2] +
            2.0 * cowan_.beta[0] + 4.0 * cowan_.deltaZX,
          cowan_.alpha[2] + 2.0 * cowan_.beta[0] +
            2.0 * cowan_.beta[1] + 4.0 * cowan_.deltaXY
        };
        for (unsigned int axis = 0; axis < 3; ++axis)
          if (!std::isfinite(normalization[axis]) ||
              std::abs(normalization[axis] - 1.0) > 1.0e-12)
            throw std::invalid_argument(
              "Cowan-z smoothing coefficients are not normalized");
        const Double rx = cowan_.ratio[0];
        const Double ry = cowan_.ratio[1];
        const Double rz = cowan_.ratio[2];
        const Double factorization[4] = {
          4.0 * ((cowan_.beta[2] + 2.0 * cowan_.deltaZX) / rz +
                 (cowan_.beta[1] + 2.0 * cowan_.deltaXY) / ry),
          4.0 * ((cowan_.beta[0] + 2.0 * cowan_.deltaXY) / rx +
                 (cowan_.beta[2] + 2.0 * cowan_.deltaYZ) / rz),
          4.0 * ((cowan_.beta[1] + 2.0 * cowan_.deltaYZ) / ry +
                 (cowan_.beta[0] + 2.0 * cowan_.deltaZX) / rx),
          16.0 * (cowan_.deltaYZ / (ry * rz) +
                  cowan_.deltaZX / (rz * rx) +
                  cowan_.deltaXY / (rx * ry))
        };
        for (unsigned int condition = 0; condition < 4; ++condition)
          if (!std::isfinite(factorization[condition]) ||
              std::abs(factorization[condition] - 1.0) > 1.0e-12)
            throw std::invalid_argument(
              "Cowan-z coefficients do not satisfy the stable dispersion factorization");
      }
  }

  RadiationFieldSample EBFieldGrid::makeRadiationSample(
      const FieldVector<Double>& electric,
      const FieldVector<Double>& magnetic) const
  {
    RadiationFieldSample sample;
    sample.electric = electric;
    sample.magnetic = magnetic;
    sample.poynting = cross(electric, magnetic);
    sample.poynting /= SI::mu0;
    sample.energyDensity = 0.5 *
      (SI::epsilon0 * electric.norm2() + magnetic.norm2() / SI::mu0);
    return sample;
  }

  Double EBFieldGrid::sampleComponent(
      const YeeComponent& component,
      const FieldVector<Double>& positionSI,
      const FieldVector<Double>& originSI,
      Double offsetX, Double offsetY, Double offsetZ) const
  {
    const Double lattice[3] = {
      (positionSI[0] - originSI[0]) / geometry_.dx - offsetX,
      (positionSI[1] - originSI[1]) / geometry_.dy - offsetY,
      (positionSI[2] - originSI[2]) / geometry_.dz - offsetZ
    };
    const std::size_t sizes[3] = {component.nx(), component.ny(),
                                  component.nz()};
    std::size_t lower[3];
    Double fraction[3];
    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        if (!std::isfinite(lattice[axis]) ||
            lattice[axis] < -0.5 - 64.0 * std::numeric_limits<Double>::epsilon() ||
            lattice[axis] > static_cast<Double>(sizes[axis]) - 0.5 +
                            64.0 * std::numeric_limits<Double>::epsilon())
          throw std::out_of_range("E/B interpolation position is outside the grid");
        if (lattice[axis] <= 0.0)
          {
            lower[axis] = 0;
            fraction[axis] = 0.0;
          }
        else if (lattice[axis] >= static_cast<Double>(sizes[axis] - 1))
          {
            lower[axis] = sizes[axis] - 2;
            fraction[axis] = 1.0;
          }
        else
          {
            lower[axis] = static_cast<std::size_t>(std::floor(lattice[axis]));
            fraction[axis] = lattice[axis] - static_cast<Double>(lower[axis]);
          }
      }

    Double value = 0.0;
    for (unsigned int c = 0; c < 2; ++c)
      for (unsigned int b = 0; b < 2; ++b)
        for (unsigned int a = 0; a < 2; ++a)
          {
            const Double weightX = a ? fraction[0] : 1.0 - fraction[0];
            const Double weightY = b ? fraction[1] : 1.0 - fraction[1];
            const Double weightZ = c ? fraction[2] : 1.0 - fraction[2];
            value += weightX * weightY * weightZ *
              component(lower[0] + a, lower[1] + b, lower[2] + c);
          }
    return value;
  }
}
