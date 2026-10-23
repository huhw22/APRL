#ifndef DIRECT_EB_MPI_H
#define DIRECT_EB_MPI_H

#include <cstddef>
#include <vector>

#include <mpi.h>

#include "eb_field.h"
#include "eb_incident.h"

namespace aprl
{
  /* z-slab communication for the staggered E/B layout.
   *
   * Adjacent ranks share the tangential E/J vertex plane. J contributions are
   * summed on both copies before Ampere's update. The neighbouring cell-centre
   * Bx/By planes are then exchanged so both ranks independently compute the
   * same interface E value. There is no per-cell MPI call. */
  class EBZSlabHaloExchange
  {
  public:
    explicit EBZSlabHaloExchange(MPI_Comm communicator);

    int rank() const;
    int size() const;
    int lowerRank() const;
    int upperRank() const;

    void installPhysicalBoundaryMask(EBFieldGrid& fields) const;
    void sumSharedTangentialCurrent(EBFieldGrid& fields);
    void synchronizeAndFilterCurrent(EBFieldGrid& fields,
                                     EBCurrentFilter filter);
    void advanceMagnetic(EBFieldGrid& fields);
    void advanceElectric(EBFieldGrid& fields,
                         EBCurrentFilter filter = EBCurrentFilter::None);
    void advance(EBFieldGrid& fields);
    void advanceWithIncident(
        EBFieldGrid& fields,
        const EBMaxwellIncidentInjector& injector,
        const SIFieldSourceSet& sources,
        const FieldVector<Double>& localGridOriginBox,
        Double timeEBox,
        const BoostFrameTransform& frame);

    std::size_t memoryBytes() const;
    /* Time spent inside MPI communication calls, including completion waits.
     * Packing and overlapped field computation are excluded. Read differences
     * between samples to obtain an interval's communication time. */
    Double communicationSeconds() const;

  private:
    void verifyGeometry(const EBFieldGrid& fields) const;
    void exchangeCurrentPlane(YeeComponent& component,
                              int tagUpper, int tagLower);
    void filterX(YeeComponent& component);
    void filterY(YeeComponent& component);
    void filterZ(YeeComponent& component, bool zNodeCentered,
                 std::size_t haloOffset);
    void exchangeFilterPlanes(const EBFieldGrid& fields);
    void exchangeElectricPlanes(const EBFieldGrid& fields);
    void beginMagneticExchange(const EBFieldGrid& fields,
                               MPI_Request requests[4]);
    void updateLowerInterface(EBFieldGrid& fields);
    void updateUpperInterface(EBFieldGrid& fields);
    void packMagneticPlane(const EBFieldGrid& fields, std::size_t k,
                           std::vector<Double>& buffer) const;
    void packElectricPlane(const EBFieldGrid& fields,
                           std::size_t nodeK, std::size_t cellK,
                           std::vector<Double>& buffer) const;

    MPI_Comm communicator_;
    int rank_;
    int size_;
    int lowerRank_;
    int upperRank_;
    Double communicationSeconds_;

    std::vector<Double> lowerMagnetic_;
    std::vector<Double> upperMagnetic_;
    std::vector<Double> sendLowerMagnetic_;
    std::vector<Double> sendUpperMagnetic_;
    std::vector<Double> receiveLowerPlane_;
    std::vector<Double> receiveUpperPlane_;
    std::vector<Double> filterLine_;
    std::vector<Double> sendLowerFilter_;
    std::vector<Double> sendUpperFilter_;
    std::vector<Double> receiveLowerFilter_;
    std::vector<Double> receiveUpperFilter_;
    std::vector<Double> lowerElectric_;
    std::vector<Double> upperElectric_;
    std::vector<Double> sendLowerElectric_;
    std::vector<Double> sendUpperElectric_;
  };
}

#endif
