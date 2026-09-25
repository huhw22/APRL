#ifndef DIRECT_EB_INITIAL_FIELD_H
#define DIRECT_EB_INITIAL_FIELD_H

#include <cstddef>
#include <vector>

#include <mpi.h>

#include "eb_field.h"
#include "eb_particles.h"
#include "config.h"
#include "fieldvector.h"

namespace fel
{
  struct EBGaussInitializationReport
  {
    std::size_t iterations;
    Double relativeResidual;
    Double maximumGaussResidual;
    Double totalCharge;
    Double meanBetaZ;
    Double meanBetaTransverse;
    Double initialElectricEnergy;
    Double initialMagneticEnergy;
    Double rmsPositionCells[3];
    Double centrePaddingRms[3];
    std::size_t staticFieldGuardCells[3];
    std::size_t temporaryBytes;

    EBGaussInitializationReport();
  };

  /* Construct the rigid common-velocity particle E/B field without retaining
   * a scalar potential in the evolution state. Charge is deposited with the
   * production CIC shape. A distributed matrix-free relativistic-Poisson
   * solve enforces div(E)=rho/epsilon0 and B=beta x E/c on the Yee lattice.
   * For CPML, the static zero-potential surface is its inner entrance so the
   * absorbing memory and fields start consistently at zero. Physical padding
   * convergence remains mandatory. */
  class EBGaussFieldInitializer
  {
  public:
    static EBGaussInitializationReport initialize(
        EBFieldGrid& fields,
        const EBGridGeometry& globalGeometry,
        std::size_t localZOffset,
        const FieldVector<Double>& localOriginSI,
        const std::vector<RelativisticParticleSI>& particles,
        InitialSelfFieldModel model,
        const std::size_t staticFieldGuardCells[3],
        Double relativeTolerance,
        std::size_t maximumIterations,
        MPI_Comm communicator);
  };
}

#endif
