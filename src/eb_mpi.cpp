#include "eb_mpi.h"

#include <algorithm>
#include <climits>
#include <stdexcept>

namespace aprl
{
  namespace
  {
    const int CURRENT_UPPER_X = 310;
    const int CURRENT_LOWER_X = 311;
    const int CURRENT_UPPER_Y = 312;
    const int CURRENT_LOWER_Y = 313;
    const int MAGNETIC_UPPER = 320;
    const int MAGNETIC_LOWER = 321;
    const int ELECTRIC_UPPER = 322;
    const int ELECTRIC_LOWER = 323;
  }

  EBZSlabHaloExchange::EBZSlabHaloExchange(MPI_Comm communicator)
    : communicator_(communicator), rank_(0), size_(1),
      lowerRank_(MPI_PROC_NULL), upperRank_(MPI_PROC_NULL),
      lowerMagnetic_(), upperMagnetic_(), sendLowerMagnetic_(),
      sendUpperMagnetic_(), receiveLowerPlane_(), receiveUpperPlane_(),
      lowerElectric_(), upperElectric_(), sendLowerElectric_(),
      sendUpperElectric_()
  {
    if (communicator_ == MPI_COMM_NULL)
      throw std::invalid_argument("E/B halo communicator cannot be null");
    MPI_Comm_rank(communicator_, &rank_);
    MPI_Comm_size(communicator_, &size_);
    lowerRank_ = (rank_ > 0) ? rank_ - 1 : MPI_PROC_NULL;
    upperRank_ = (rank_ + 1 < size_) ? rank_ + 1 : MPI_PROC_NULL;
  }

  int EBZSlabHaloExchange::rank() const { return rank_; }
  int EBZSlabHaloExchange::size() const { return size_; }
  int EBZSlabHaloExchange::lowerRank() const { return lowerRank_; }
  int EBZSlabHaloExchange::upperRank() const { return upperRank_; }

  void EBZSlabHaloExchange::installPhysicalBoundaryMask(
      EBFieldGrid& fields) const
  {
    fields.setBoundary(std::unique_ptr<EBBoundaryOperator>(
      new PerfectElectricConductorBoundary(lowerRank_ == MPI_PROC_NULL,
                                           upperRank_ == MPI_PROC_NULL)));
  }

  void EBZSlabHaloExchange::sumSharedTangentialCurrent(EBFieldGrid& fields)
  {
    verifyGeometry(fields);
    if (size_ == 1) return;
    exchangeCurrentPlane(fields.jx(), CURRENT_UPPER_X, CURRENT_LOWER_X);
    exchangeCurrentPlane(fields.jy(), CURRENT_UPPER_Y, CURRENT_LOWER_Y);
  }

  void EBZSlabHaloExchange::advanceMagnetic(EBFieldGrid& fields)
  {
    verifyGeometry(fields);
    if (size_ == 1 ||
        fields.geometry().solver == EBMaxwellSolver::Yee)
      {
        fields.advanceMagnetic();
        return;
      }

    exchangeElectricPlanes(fields);
    const std::size_t exValues = fields.ex().nx() * fields.ex().ny();
    const std::size_t eyValues = fields.ey().nx() * fields.ey().ny();
    EBElectricHaloView halo;
    if (lowerRank_ != MPI_PROC_NULL)
      {
        halo.lowerEx = &lowerElectric_[0];
        halo.lowerEy = &lowerElectric_[exValues];
        halo.lowerEz = &lowerElectric_[exValues + eyValues];
      }
    if (upperRank_ != MPI_PROC_NULL)
      {
        halo.upperEx = &upperElectric_[0];
        halo.upperEy = &upperElectric_[exValues];
        halo.upperEz = &upperElectric_[exValues + eyValues];
      }
    fields.advanceMagnetic(halo);
  }

  void EBZSlabHaloExchange::advanceElectric(EBFieldGrid& fields)
  {
    verifyGeometry(fields);
    if (size_ == 1)
      {
        fields.advanceElectric();
        return;
      }

    sumSharedTangentialCurrent(fields);
    exchangeMagneticPlanes(fields);
    fields.advanceElectric();
    if (lowerRank_ != MPI_PROC_NULL) updateLowerInterface(fields);
    if (upperRank_ != MPI_PROC_NULL) updateUpperInterface(fields);
  }

  void EBZSlabHaloExchange::advance(EBFieldGrid& fields)
  {
    advanceMagnetic(fields);
    advanceElectric(fields);
  }

  void EBZSlabHaloExchange::advanceWithIncident(
      EBFieldGrid& fields,
      const EBMaxwellIncidentInjector& injector,
      const SIFieldSourceSet& sources,
      const FieldVector<Double>& localGridOriginBox,
      Double timeEBox,
      const BoostFrameTransform& frame)
  {
    /* The magnetic TF/SF correction must cross MPI with B.  The electric
     * correction must instead wait until both copies of a shared E plane have
     * received the same interface curl update. */
    advanceMagnetic(fields);
    injector.correctAfterMagneticUpdate(
      fields, sources, localGridOriginBox, timeEBox, frame);
    advanceElectric(fields);
    injector.correctAfterElectricUpdate(
      fields, sources, localGridOriginBox, timeEBox, frame);
  }

  std::size_t EBZSlabHaloExchange::memoryBytes() const
  {
    return sizeof(Double) * (lowerMagnetic_.capacity() +
      upperMagnetic_.capacity() + sendLowerMagnetic_.capacity() +
      sendUpperMagnetic_.capacity() + receiveLowerPlane_.capacity() +
      receiveUpperPlane_.capacity() + lowerElectric_.capacity() +
      upperElectric_.capacity() + sendLowerElectric_.capacity() +
      sendUpperElectric_.capacity());
  }

  void EBZSlabHaloExchange::verifyGeometry(const EBFieldGrid& fields) const
  {
    const EBGridGeometry& g = fields.geometry();
    if (g.nz < 1 || g.nx < 2 || g.ny < 2)
      throw std::invalid_argument("E/B z slab is too small for halo exchange");
  }

  void EBZSlabHaloExchange::exchangeCurrentPlane(YeeComponent& component,
                                                  int tagUpper,
                                                  int tagLower)
  {
    const std::size_t planeValues = component.nx() * component.ny();
    if (planeValues > static_cast<std::size_t>(INT_MAX))
      throw std::overflow_error("E/B current halo exceeds MPI int count");
    receiveLowerPlane_.resize(planeValues);
    receiveUpperPlane_.resize(planeValues);
    Double* lowerPlane = component.data();
    Double* upperPlane = component.data() +
                         (component.nz() - 1) * planeValues;

    MPI_Sendrecv(upperPlane, static_cast<int>(planeValues), MPI_DOUBLE,
                 upperRank_, tagUpper,
                 receiveLowerPlane_.data(), static_cast<int>(planeValues), MPI_DOUBLE,
                 lowerRank_, tagUpper, communicator_, MPI_STATUS_IGNORE);

    MPI_Sendrecv(lowerPlane, static_cast<int>(planeValues), MPI_DOUBLE,
                 lowerRank_, tagLower,
                 receiveUpperPlane_.data(), static_cast<int>(planeValues), MPI_DOUBLE,
                 upperRank_, tagLower, communicator_, MPI_STATUS_IGNORE);

    if (lowerRank_ != MPI_PROC_NULL)
      for (std::size_t value = 0; value < planeValues; ++value)
        lowerPlane[value] += receiveLowerPlane_[value];
    if (upperRank_ != MPI_PROC_NULL)
      for (std::size_t value = 0; value < planeValues; ++value)
        upperPlane[value] += receiveUpperPlane_[value];
  }

  void EBZSlabHaloExchange::exchangeMagneticPlanes(
      const EBFieldGrid& fields)
  {
    const std::size_t planeValues = fields.bx().nx() * fields.bx().ny() +
                                    fields.by().nx() * fields.by().ny();
    if (planeValues > static_cast<std::size_t>(INT_MAX))
      throw std::overflow_error("E/B magnetic halo exceeds MPI int count");
    lowerMagnetic_.resize(planeValues);
    upperMagnetic_.resize(planeValues);
    sendLowerMagnetic_.resize(planeValues);
    sendUpperMagnetic_.resize(planeValues);
    packMagneticPlane(fields, 0, sendLowerMagnetic_);
    packMagneticPlane(fields, fields.geometry().nz - 1,
                      sendUpperMagnetic_);

    MPI_Sendrecv(sendUpperMagnetic_.data(), static_cast<int>(planeValues),
                 MPI_DOUBLE, upperRank_, MAGNETIC_UPPER,
                 lowerMagnetic_.data(), static_cast<int>(planeValues),
                 MPI_DOUBLE, lowerRank_, MAGNETIC_UPPER,
                 communicator_, MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendLowerMagnetic_.data(), static_cast<int>(planeValues),
                 MPI_DOUBLE, lowerRank_, MAGNETIC_LOWER,
                 upperMagnetic_.data(), static_cast<int>(planeValues),
                 MPI_DOUBLE, upperRank_, MAGNETIC_LOWER,
                 communicator_, MPI_STATUS_IGNORE);
  }

  void EBZSlabHaloExchange::exchangeElectricPlanes(
      const EBFieldGrid& fields)
  {
    const std::size_t planeValues =
      fields.ex().nx() * fields.ex().ny() +
      fields.ey().nx() * fields.ey().ny() +
      fields.ez().nx() * fields.ez().ny();
    if (planeValues > static_cast<std::size_t>(INT_MAX))
      throw std::overflow_error("E/B electric halo exceeds MPI int count");
    lowerElectric_.resize(planeValues);
    upperElectric_.resize(planeValues);
    sendLowerElectric_.resize(planeValues);
    sendUpperElectric_.resize(planeValues);
    const std::size_t cells = fields.geometry().nz;
    packElectricPlane(fields, std::min<std::size_t>(1, cells), 0,
                      sendLowerElectric_);
    packElectricPlane(fields, cells - 1, cells - 1,
                      sendUpperElectric_);

    MPI_Sendrecv(sendUpperElectric_.data(), static_cast<int>(planeValues),
                 MPI_DOUBLE, upperRank_, ELECTRIC_UPPER,
                 lowerElectric_.data(), static_cast<int>(planeValues),
                 MPI_DOUBLE, lowerRank_, ELECTRIC_UPPER,
                 communicator_, MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendLowerElectric_.data(), static_cast<int>(planeValues),
                 MPI_DOUBLE, lowerRank_, ELECTRIC_LOWER,
                 upperElectric_.data(), static_cast<int>(planeValues),
                 MPI_DOUBLE, upperRank_, ELECTRIC_LOWER,
                 communicator_, MPI_STATUS_IGNORE);
  }

  void EBZSlabHaloExchange::updateLowerInterface(EBFieldGrid& fields)
  {
    const EBGridGeometry& g = fields.geometry();
    const std::size_t bxValues = fields.bx().nx() * fields.bx().ny();
    const Double* bxBelow = &lowerMagnetic_[0];
    const Double* byBelow = &lowerMagnetic_[bxValues];
    const Double c2dt = g.dt / (SI::mu0 * SI::epsilon0);
    const Double currentScale = g.dt / SI::epsilon0;

    for (std::size_t j = 1; j < g.ny; ++j)
      for (std::size_t i = 0; i < g.nx; ++i)
        {
          const Double bzDerivative =
            (fields.bz()(i, j, 0) - fields.bz()(i, j - 1, 0)) / g.dy;
          const Double byAbove = fields.by()(i, j, 0);
          const Double byLower = byBelow[j * g.nx + i];
          fields.ex()(i, j, 0) += c2dt *
            (bzDerivative - (byAbove - byLower) / g.dz) -
            currentScale * fields.jx()(i, j, 0);
        }

    for (std::size_t j = 0; j < g.ny; ++j)
      for (std::size_t i = 1; i < g.nx; ++i)
        {
          const Double bxAbove = fields.bx()(i, j, 0);
          const Double bxLower = bxBelow[j * (g.nx + 1) + i];
          const Double bzDerivative =
            (fields.bz()(i, j, 0) - fields.bz()(i - 1, j, 0)) / g.dx;
          fields.ey()(i, j, 0) += c2dt *
            ((bxAbove - bxLower) / g.dz - bzDerivative) -
            currentScale * fields.jy()(i, j, 0);
        }
  }

  void EBZSlabHaloExchange::updateUpperInterface(EBFieldGrid& fields)
  {
    const EBGridGeometry& g = fields.geometry();
    const std::size_t bxValues = fields.bx().nx() * fields.bx().ny();
    const Double* bxAbove = &upperMagnetic_[0];
    const Double* byAbove = &upperMagnetic_[bxValues];
    const std::size_t k = g.nz;
    const std::size_t below = g.nz - 1;
    const Double c2dt = g.dt / (SI::mu0 * SI::epsilon0);
    const Double currentScale = g.dt / SI::epsilon0;

    for (std::size_t j = 1; j < g.ny; ++j)
      for (std::size_t i = 0; i < g.nx; ++i)
        {
          const Double bzDerivative =
            (fields.bz()(i, j, k) - fields.bz()(i, j - 1, k)) / g.dy;
          const Double byUpper = byAbove[j * g.nx + i];
          const Double byLower = fields.by()(i, j, below);
          fields.ex()(i, j, k) += c2dt *
            (bzDerivative - (byUpper - byLower) / g.dz) -
            currentScale * fields.jx()(i, j, k);
        }

    for (std::size_t j = 0; j < g.ny; ++j)
      for (std::size_t i = 1; i < g.nx; ++i)
        {
          const Double bxUpper = bxAbove[j * (g.nx + 1) + i];
          const Double bxLower = fields.bx()(i, j, below);
          const Double bzDerivative =
            (fields.bz()(i, j, k) - fields.bz()(i - 1, j, k)) / g.dx;
          fields.ey()(i, j, k) += c2dt *
            ((bxUpper - bxLower) / g.dz - bzDerivative) -
            currentScale * fields.jy()(i, j, k);
        }
  }

  void EBZSlabHaloExchange::packMagneticPlane(
      const EBFieldGrid& fields, std::size_t k,
      std::vector<Double>& buffer) const
  {
    const std::size_t bxValues = fields.bx().nx() * fields.bx().ny();
    const std::size_t byValues = fields.by().nx() * fields.by().ny();
    buffer.resize(bxValues + byValues);
    std::copy(fields.bx().data() + k * bxValues,
              fields.bx().data() + (k + 1) * bxValues, buffer.begin());
    std::copy(fields.by().data() + k * byValues,
              fields.by().data() + (k + 1) * byValues,
              buffer.begin() + bxValues);
  }

  void EBZSlabHaloExchange::packElectricPlane(
      const EBFieldGrid& fields, std::size_t nodeK, std::size_t cellK,
      std::vector<Double>& buffer) const
  {
    const std::size_t exValues = fields.ex().nx() * fields.ex().ny();
    const std::size_t eyValues = fields.ey().nx() * fields.ey().ny();
    const std::size_t ezValues = fields.ez().nx() * fields.ez().ny();
    buffer.resize(exValues + eyValues + ezValues);
    std::copy(fields.ex().data() + nodeK * exValues,
              fields.ex().data() + (nodeK + 1) * exValues,
              buffer.begin());
    std::copy(fields.ey().data() + nodeK * eyValues,
              fields.ey().data() + (nodeK + 1) * eyValues,
              buffer.begin() + exValues);
    std::copy(fields.ez().data() + cellK * ezValues,
              fields.ez().data() + (cellK + 1) * ezValues,
              buffer.begin() + exValues + eyValues);
  }
}
