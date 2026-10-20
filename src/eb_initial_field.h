#ifndef DIRECT_EB_INITIAL_FIELD_H
#define DIRECT_EB_INITIAL_FIELD_H

#include <cstddef>
#include <vector>

#include <mpi.h>

#include "eb_field.h"
#include "eb_particles.h"
#include "config.h"
#include "fieldvector.h"

namespace aprl
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
    std::size_t transversePoissonPaddingCells;
    std::size_t temporaryBytes;

    EBGaussInitializationReport();
  };

  /* Construct the rigid common-velocity particle E/B field without retaining
   * a scalar potential in the evolution state. Charge is deposited with the
   * production CIC shape. A distributed matrix-free relativistic-Poisson
   * solve enforces div(E)=rho/epsilon0 and B=beta x E/c on the Yee lattice.
   * Normally the CPML entrance is the static zero-potential surface. The
   * opt-in transverse CPML-source solve instead pads a temporary Poisson
   * lattice and retains only physical-region E/B; evolving CPML fields and
   * memories still begin at zero. Physical padding convergence remains
   * mandatory. */
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
        Double maximumMeanTransverseBeta,
        const std::size_t staticFieldGuardCells[3],
        Double relativeTolerance,
        std::size_t maximumIterations,
        MPI_Comm communicator,
        bool includeTransverseCpmlParticles = false);
  };
}

#endif
