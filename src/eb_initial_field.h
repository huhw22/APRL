#ifndef DIRECT_EB_INITIAL_FIELD_H
#define DIRECT_EB_INITIAL_FIELD_H

#include <cstddef>
#include <vector>

#include <mpi.h>

#include "eb_field.h"
#include "eb_particles.h"
#include "fieldvector.h"

namespace fel
{
  struct EBGaussInitializationReport
  {
    std::size_t iterations;
    Double relativeResidual;
    Double maximumGaussResidual;
    Double totalCharge;
    std::size_t temporaryBytes;

    EBGaussInitializationReport();
  };

  /* Construct the longitudinal (curl-free) part of the initial particle
   * field without retaining a scalar potential in the evolution state.
   * Charge is deposited with the same CIC shape as the charge-conserving
   * trajectory depositor. A distributed matrix-free CG solve then enforces
   *
   *     div(E) = rho / epsilon0
   *
   * on the interior Yee vertices. The outer potential is zero, matching the
   * finite-domain field boundary; box-size/padding convergence therefore
   * remains a physical requirement. */
  class EBGaussFieldInitializer
  {
  public:
    static EBGaussInitializationReport initialize(
        EBFieldGrid& fields,
        const EBGridGeometry& globalGeometry,
        std::size_t localZOffset,
        const FieldVector<Double>& localOriginSI,
        const std::vector<RelativisticParticleSI>& particles,
        Double relativeTolerance,
        std::size_t maximumIterations,
        MPI_Comm communicator);
  };
}

#endif
