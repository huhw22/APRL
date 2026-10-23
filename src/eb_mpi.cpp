#include "eb_mpi.h"

#include <algorithm>
#include <climits>
#include <cstdlib>
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
    const int FILTER_UPPER = 330;
    const int FILTER_LOWER = 331;

    void requireSuccessfulTransfer(int status, MPI_Comm communicator)
    {
      if (status == MPI_SUCCESS) return;
      // A failed post/completion can leave buffers owned by MPI. Do not unwind
      // and release them when an application installed MPI_ERRORS_RETURN.
      MPI_Abort(communicator, status);
      std::abort();
    }

    class CommunicationTimer
    {
    public:
      explicit CommunicationTimer(Double& elapsed)
        : elapsed_(elapsed), start_(MPI_Wtime()) {}
      ~CommunicationTimer() { elapsed_ += MPI_Wtime() - start_; }

    private:
      Double& elapsed_;
      Double start_;
    };

    /* The field update can throw. Finish every posted transfer on either path
     * before its member buffers can be reused or destroyed. MPI stays on the
     * calling thread, as required by MPI_THREAD_FUNNELED. */
    class PendingMagneticExchange
    {
    public:
      PendingMagneticExchange(Double& elapsed, MPI_Comm communicator)
        : elapsed_(elapsed), communicator_(communicator)
      {
        for (unsigned int index = 0; index < 4; ++index)
          requests[index] = MPI_REQUEST_NULL;
      }
      ~PendingMagneticExchange() { wait(); }

      void wait()
      {
        bool pending = false;
        for (unsigned int index = 0; index < 4; ++index)
          pending = pending || requests[index] != MPI_REQUEST_NULL;
        if (!pending) return;
        CommunicationTimer timer(elapsed_);
        requireSuccessfulTransfer(
          MPI_Waitall(4, requests, MPI_STATUSES_IGNORE), communicator_);
      }

      MPI_Request requests[4];

    private:
      PendingMagneticExchange(const PendingMagneticExchange&) = delete;
      PendingMagneticExchange& operator=(const PendingMagneticExchange&) = delete;
      Double& elapsed_;
      MPI_Comm communicator_;
    };

    Double compensatedBilinear(Double minusTwo, Double minusOne,
                               Double centre, Double plusOne,
                               Double plusTwo)
    {
      return (-minusTwo + 4.0 * minusOne + 10.0 * centre +
              4.0 * plusOne - plusTwo) / 16.0;
    }
  }

  EBZSlabHaloExchange::EBZSlabHaloExchange(MPI_Comm communicator)
    : communicator_(communicator), rank_(0), size_(1),
      lowerRank_(MPI_PROC_NULL), upperRank_(MPI_PROC_NULL),
      communicationSeconds_(0.0),
      lowerMagnetic_(), upperMagnetic_(), sendLowerMagnetic_(),
      sendUpperMagnetic_(), receiveLowerPlane_(), receiveUpperPlane_(),
      filterLine_(), sendLowerFilter_(), sendUpperFilter_(),
      receiveLowerFilter_(), receiveUpperFilter_(),
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
  Double EBZSlabHaloExchange::communicationSeconds() const
  { return communicationSeconds_; }

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

  void EBZSlabHaloExchange::advanceElectric(EBFieldGrid& fields,
                                             EBCurrentFilter filter)
  {
    verifyGeometry(fields);
    if (size_ == 1)
      {
        synchronizeAndFilterCurrent(fields, filter);
        fields.advanceElectric();
        return;
      }

    synchronizeAndFilterCurrent(fields, filter);
    PendingMagneticExchange exchange(communicationSeconds_, communicator_);
    beginMagneticExchange(fields, exchange.requests);
    /* The internal E update and its boundary operator only read local B.
     * Neighbour B values are needed solely by the interface updates below. */
    fields.advanceElectric();
    exchange.wait();
    if (lowerRank_ != MPI_PROC_NULL) updateLowerInterface(fields);
    if (upperRank_ != MPI_PROC_NULL) updateUpperInterface(fields);
  }

  void EBZSlabHaloExchange::synchronizeAndFilterCurrent(
      EBFieldGrid& fields, EBCurrentFilter filter)
  {
    verifyGeometry(fields);
    sumSharedTangentialCurrent(fields);
    if (filter == EBCurrentFilter::None) return;
    if (filter != EBCurrentFilter::CompensatedBilinear)
      throw std::invalid_argument("Unsupported E/B current filter");
    if (fields.geometry().solver != EBMaxwellSolver::CowanZ)
      throw std::invalid_argument(
        "Compensated current filtering is restricted to Cowan-z runs");

    /* The separable filter is a bilinear pass followed by its standard
     * compensation pass, collapsed to the five-point stencil
     * [-1, 4, 10, 4, -1]/16.  Its one-axis response is
     * 1-sin^4(k*delta/2): fourth-order-flat at k=0 and zero at Nyquist.
     * Each one-dimensional pass leaves its two physical endpoint samples
     * untouched rather than extending a synthetic source beyond the box. */
    filterX(fields.jx()); filterY(fields.jx());
    filterX(fields.jy()); filterY(fields.jy());
    filterX(fields.jz()); filterY(fields.jz());
    if (size_ > 1) exchangeFilterPlanes(fields);
    const std::size_t jxHaloValues =
      2 * fields.jx().nx() * fields.jx().ny();
    const std::size_t jyHaloValues =
      2 * fields.jy().nx() * fields.jy().ny();
    filterZ(fields.jx(), true, 0);
    filterZ(fields.jy(), true, jxHaloValues);
    filterZ(fields.jz(), false, jxHaloValues + jyHaloValues);
  }

  void EBZSlabHaloExchange::filterX(YeeComponent& component)
  {
    const std::size_t n = component.nx();
    if (n < 5) return;
    filterLine_.resize(n);
    for (std::size_t k = 0; k < component.nz(); ++k)
      for (std::size_t j = 0; j < component.ny(); ++j)
        {
          for (std::size_t i = 0; i < n; ++i)
            filterLine_[i] = component(i, j, k);
          for (std::size_t i = 2; i + 2 < n; ++i)
            component(i, j, k) = compensatedBilinear(
              filterLine_[i - 2], filterLine_[i - 1], filterLine_[i],
              filterLine_[i + 1], filterLine_[i + 2]);
        }
  }

  void EBZSlabHaloExchange::filterY(YeeComponent& component)
  {
    const std::size_t n = component.ny();
    if (n < 5) return;
    filterLine_.resize(n);
    for (std::size_t k = 0; k < component.nz(); ++k)
      for (std::size_t i = 0; i < component.nx(); ++i)
        {
          for (std::size_t j = 0; j < n; ++j)
            filterLine_[j] = component(i, j, k);
          for (std::size_t j = 2; j + 2 < n; ++j)
            component(i, j, k) = compensatedBilinear(
              filterLine_[j - 2], filterLine_[j - 1], filterLine_[j],
              filterLine_[j + 1], filterLine_[j + 2]);
        }
  }

  void EBZSlabHaloExchange::exchangeFilterPlanes(
      const EBFieldGrid& fields)
  {
    const YeeComponent* component[3] = {
      &fields.jx(), &fields.jy(), &fields.jz()
    };
    const bool zNodeCentered[3] = {true, true, false};
    std::size_t values = 0;
    for (unsigned int field = 0; field < 3; ++field)
      values += 2 * component[field]->nx() * component[field]->ny();
    if (values > static_cast<std::size_t>(INT_MAX))
      throw std::overflow_error("E/B filter halo exceeds MPI int count");
    sendLowerFilter_.resize(values);
    sendUpperFilter_.resize(values);
    receiveLowerFilter_.resize(values);
    receiveUpperFilter_.resize(values);
    std::size_t offset = 0;
    for (unsigned int field = 0; field < 3; ++field)
      {
        const std::size_t planeValues =
          component[field]->nx() * component[field]->ny();
        const std::size_t n = component[field]->nz();
        const std::size_t lowerStart = zNodeCentered[field] ? 1 : 0;
        const std::size_t upperStart =
          zNodeCentered[field] ? n - 3 : n - 2;
        for (std::size_t layer = 0; layer < 2; ++layer)
          {
            const Double* lower = component[field]->data() +
              (lowerStart + layer) * planeValues;
            const Double* upper = component[field]->data() +
              (upperStart + layer) * planeValues;
            std::copy(lower, lower + planeValues,
                      sendLowerFilter_.begin() + offset +
                        layer * planeValues);
            std::copy(upper, upper + planeValues,
                      sendUpperFilter_.begin() + offset +
                        layer * planeValues);
          }
        offset += 2 * planeValues;
      }

    CommunicationTimer timer(communicationSeconds_);
    MPI_Sendrecv(sendUpperFilter_.data(), static_cast<int>(values), MPI_DOUBLE,
                 upperRank_, FILTER_UPPER,
                 receiveLowerFilter_.data(), static_cast<int>(values),
                 MPI_DOUBLE, lowerRank_, FILTER_UPPER,
                 communicator_, MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendLowerFilter_.data(), static_cast<int>(values), MPI_DOUBLE,
                 lowerRank_, FILTER_LOWER,
                 receiveUpperFilter_.data(), static_cast<int>(values),
                 MPI_DOUBLE, upperRank_, FILTER_LOWER,
                 communicator_, MPI_STATUS_IGNORE);
  }

  void EBZSlabHaloExchange::filterZ(YeeComponent& component,
                                     bool zNodeCentered,
                                     std::size_t haloOffset)
  {
    const std::size_t n = component.nz();
    if (n < (zNodeCentered ? 3u : 2u))
      throw std::invalid_argument(
        "Current filter needs at least two z cells per MPI rank");
    filterLine_.resize(n);
    const std::size_t planeValues = component.nx() * component.ny();

    for (std::size_t j = 0; j < component.ny(); ++j)
      for (std::size_t i = 0; i < component.nx(); ++i)
        {
          for (std::size_t k = 0; k < n; ++k)
            filterLine_[k] = component(i, j, k);
          for (std::size_t k = 0; k < n; ++k)
            {
              if ((lowerRank_ == MPI_PROC_NULL && k < 2) ||
                  (upperRank_ == MPI_PROC_NULL && k + 2 >= n))
                continue;
              const std::size_t planeIndex = j * component.nx() + i;
              const Double minusTwo = k >= 2 ? filterLine_[k - 2] :
                receiveLowerFilter_[haloOffset +
                                    (k == 0 ? 0 : 1) * planeValues +
                                    planeIndex];
              const Double minusOne = k >= 1 ? filterLine_[k - 1] :
                receiveLowerFilter_[haloOffset + planeValues + planeIndex];
              const Double plusOne = k + 1 < n ? filterLine_[k + 1] :
                receiveUpperFilter_[haloOffset + planeIndex];
              const Double plusTwo = k + 2 < n ? filterLine_[k + 2] :
                receiveUpperFilter_[haloOffset +
                                    (k + 1 == n ? 1 : 0) * planeValues +
                                    planeIndex];
              component(i, j, k) = compensatedBilinear(
                minusTwo, minusOne, filterLine_[k], plusOne, plusTwo);
            }
        }
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
      receiveUpperPlane_.capacity() + filterLine_.capacity() +
      sendLowerFilter_.capacity() + sendUpperFilter_.capacity() +
      receiveLowerFilter_.capacity() + receiveUpperFilter_.capacity() +
      lowerElectric_.capacity() +
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

    {
      CommunicationTimer timer(communicationSeconds_);
      MPI_Sendrecv(upperPlane, static_cast<int>(planeValues), MPI_DOUBLE,
                   upperRank_, tagUpper,
                   receiveLowerPlane_.data(), static_cast<int>(planeValues),
                   MPI_DOUBLE, lowerRank_, tagUpper,
                   communicator_, MPI_STATUS_IGNORE);

      MPI_Sendrecv(lowerPlane, static_cast<int>(planeValues), MPI_DOUBLE,
                   lowerRank_, tagLower,
                   receiveUpperPlane_.data(), static_cast<int>(planeValues),
                   MPI_DOUBLE, upperRank_, tagLower,
                   communicator_, MPI_STATUS_IGNORE);
    }

    if (lowerRank_ != MPI_PROC_NULL)
      for (std::size_t value = 0; value < planeValues; ++value)
        lowerPlane[value] += receiveLowerPlane_[value];
    if (upperRank_ != MPI_PROC_NULL)
      for (std::size_t value = 0; value < planeValues; ++value)
        upperPlane[value] += receiveUpperPlane_[value];
  }

  void EBZSlabHaloExchange::beginMagneticExchange(
      const EBFieldGrid& fields, MPI_Request requests[4])
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

    CommunicationTimer timer(communicationSeconds_);
    requireSuccessfulTransfer(MPI_Irecv(
              lowerMagnetic_.data(), static_cast<int>(planeValues),
              MPI_DOUBLE, lowerRank_, MAGNETIC_UPPER,
              communicator_, &requests[0]), communicator_);
    requireSuccessfulTransfer(MPI_Irecv(
              upperMagnetic_.data(), static_cast<int>(planeValues),
              MPI_DOUBLE, upperRank_, MAGNETIC_LOWER,
              communicator_, &requests[1]), communicator_);
    requireSuccessfulTransfer(MPI_Isend(
              sendUpperMagnetic_.data(), static_cast<int>(planeValues),
              MPI_DOUBLE, upperRank_, MAGNETIC_UPPER,
              communicator_, &requests[2]), communicator_);
    requireSuccessfulTransfer(MPI_Isend(
              sendLowerMagnetic_.data(), static_cast<int>(planeValues),
              MPI_DOUBLE, lowerRank_, MAGNETIC_LOWER,
              communicator_, &requests[3]), communicator_);
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

    CommunicationTimer timer(communicationSeconds_);
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
