#ifndef APRL_VERIFICATION_COWAN_REFERENCE_FIELD_H
#define APRL_VERIFICATION_COWAN_REFERENCE_FIELD_H

#include <cstddef>
#include <vector>

#include <mpi.h>

#include "eb_field.h"
#include "eb_particles.h"

namespace aprl {
namespace verification {

struct CowanReferenceReport {
  std::size_t transforms = 0;
  double relativeGaussL2 = 0;
  double maximumGaussResidual = 0;
  double minimumDenominator = 0;
  double minimumSmoothingSymbol = 0;
  double minimumLongitudinalSymbol = 0;
  double initializationSeconds = 0;
  double sourceCharge = 0;
};

struct CowanMagneticDivergenceReport {
  std::size_t samples = 0;
  double ordinaryMaximum = 0;          // T/m
  double modifiedMaximum = 0;          // T/m
  double ordinaryRelativeL2 = 0;
  double modifiedRelativeL2 = 0;
  double ordinaryRelativeMaximum = 0;
  double modifiedRelativeMaximum = 0;
  double normalizationRms = 0;         // T/m
  double normalizationMaximum = 0;     // T/m
};

// Cell-centre divergences at cells at least two cells from every PEC face.
// The modified divergence uses the same transverse smoothing S_i as Cowan's
// Faraday derivative. Both L2 norms divide by the L2 norm of the local scale
// sqrt(sum_i (B_i,lower^2+B_i,upper^2)/(2 h_i^2)). Maximum-relative values
// divide by the maximum of that scale. A globally zero B gives all zeros.
// Neither divergence is classified as PASS/FAIL by this diagnostic.
CowanMagneticDivergenceReport cowanMagneticDivergence(
    const EBFieldGrid& fields);

// Research-only, one rank in communicator, complete PEC box and CowanZ geometry.
// Deposits the unmodified particles with the production CIC shape. A DST-I
// solves -div(D*phi-a zhat Dz phi)=rho/epsilon0 with phi=0 on every face.
// Here D* is the actual Cowan derivative and a=v^2 epsilon0 mu0.
// At v=0 this simultaneously satisfies ordinary Gauss and curl*(E)=0,
// including the PEC parity extension used by the production Faraday update.
// At v!=0 the field is a rigid-translation, semi-discrete approximation:
// E=-D*phi+a zhat Dz phi, B=v epsilon0 mu0 zhat cross E. The B half-step
// translation uses the same linear interpolation as the existing initializer.
// It is not an exact fully discrete moving-source or discrete-CIC solution.
// The moving construction preserves the modified magnetic divergence in the
// translation-invariant interior, not necessarily the ordinary div(B).
// This is the stationary construction of Cowan et al. (2013), section V,
// equations 27--30; section II.C distinguishes the two magnetic divergences:
// https://inspirehep.net/files/dbbec6351aba54ae1a1a61f6d0b9870c
// beta is the caller's chosen common velocity/c approximation and must be
// recorded alongside the particle distribution; it is not inferred here.
// This routine overwrites E/B and clears J; it never changes particle x/p.
CowanReferenceReport initializeCowanReference(
    EBFieldGrid& fields,
    const FieldVector<double>& origin,
    const std::vector<RelativisticParticleSI>& particles,
    double beta,
    MPI_Comm communicator = MPI_COMM_WORLD);

}  // namespace verification
}  // namespace aprl
#endif
