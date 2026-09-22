#ifndef DIRECT_EB_DEPOSITION_H
#define DIRECT_EB_DEPOSITION_H

#include <cstddef>

#include "eb_field.h"
#include "fieldvector.h"

namespace fel
{
  struct CurrentDepositResult
  {
    std::size_t subsegments;
    Double depositedCharge;

    CurrentDepositResult();
  };

  /* First-order, charge-conserving trajectory deposition for the Yee lattice.
   *
   * Charge lives on Yee vertices and current lives on the matching electric
   * edges. A trajectory is split exactly at every crossed cell face. Within
   * each cell the transverse shape factors are integrated analytically along
   * the straight segment. Consequently,
   *
   *   (rho[n+1] - rho[n]) / dt + discrete_divergence(J) = 0
   *
   * to round-off, including diagonal motion and multi-cell crossings.
   */
  class ChargeConservingCurrentDepositor
  {
  public:
    ChargeConservingCurrentDepositor(EBFieldGrid& fields,
                                     const FieldVector<Double>& originSI);

    CurrentDepositResult depositSegment(
        const FieldVector<Double>& startSI,
        const FieldVector<Double>& endSI,
        Double chargeCoulomb);

    /* CIC charge deposition is provided for continuity/Gauss diagnostics. The
     * rho array must have dimensions (nx+1, ny+1, nz+1) and stores C/m^3. */
    void depositCharge(const FieldVector<Double>& positionSI,
                       Double chargeCoulomb,
                       YeeComponent& rho) const;

    static Double maxContinuityResidual(
        const EBFieldGrid& fields,
        const YeeComponent& rhoBefore,
        const YeeComponent& rhoAfter);

  private:
    struct NormalizedPosition
    {
      Double x;
      Double y;
      Double z;
    };

    void validatePosition(const FieldVector<Double>& positionSI) const;
    NormalizedPosition normalize(
        const FieldVector<Double>& positionSI) const;
    std::size_t containingCell(Double coordinate, Double direction,
                               std::size_t cells) const;
    void depositCellSegment(std::size_t i, std::size_t j, std::size_t k,
                            const NormalizedPosition& start,
                            const NormalizedPosition& end,
                            Double chargeCoulomb);

    EBFieldGrid& fields_;
    FieldVector<Double> originSI_;
  };
}

#endif
