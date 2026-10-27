#ifndef DIRECT_EB_INITIAL_FIELD_H
#define DIRECT_EB_INITIAL_FIELD_H

#include <cstddef>
#include <vector>

#include <mpi.h>

#include "eb_field.h"
#include "eb_particles.h"
#include "boostframe.h"
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
   * The opt-in Cowan model uses -div[(D* - a zhat D_z)phi]=rho/epsilon0,
   * a=v_z^2 epsilon0 mu0, and B=v_z epsilon0 mu0 zhat x E with the same
   * rigid-translation half-step interpolation. Its 27-point distributed
   * CG solve retains one halo; field construction reuses a CG vector and
   * exchanges filtered potentials rather than allocating another volume.
   * Positive S_x, S_y and S_z-a over all transverse cosine arguments are
   * required as a conservative sufficient positive-definiteness bound.
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
        bool includeTransverseCpmlParticles = false,
        bool holdTransverseUntilActivation = false,
        const BoostFrameTransform* frame = 0);
  };
}

#endif
