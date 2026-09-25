#ifndef DIRECT_EB_PARTICLE_BOUNDARY_H
#define DIRECT_EB_PARTICLE_BOUNDARY_H

#include <cstddef>
#include <vector>

#include "eb_cpml.h"
#include "eb_field.h"
#include "fieldvector.h"

namespace aprl
{
  enum class ParticleBoundaryFace
  {
    LowerX = 0,
    UpperX = 1,
    LowerY = 2,
    UpperY = 3,
    LowerZ = 4,
    UpperZ = 5
  };

  const char* particleBoundaryFaceName(ParticleBoundaryFace face);

  struct ParticleBoundaryHit
  {
    FieldVector<Double> position;
    Double fraction;
    ParticleBoundaryFace face;

    ParticleBoundaryHit();
  };

  /* Geometry and matched current damping for particles entering CPML.
   * The inner CPML surface is the end of the physical particle domain. */
  class ParticleCPMLRegion
  {
  public:
    ParticleCPMLRegion(const EBGridGeometry& globalGeometry,
                       const FieldVector<Double>& globalOriginSI,
                       const EBCPMLParameters& parameters);

    bool enabled() const;
    bool containsPhysical(const FieldVector<Double>& positionSI) const;
    bool firstEntry(const FieldVector<Double>& startSI,
                    const FieldVector<Double>& endSI,
                    ParticleBoundaryHit& hit) const;

    /* Per-particle current-only damping over a straight ballistic segment.
     * This integrates sigma/epsilon0 along the actual simulation-frame path,
     * avoiding a fixed assumed particle velocity. */
    Double currentDampingFactor(
        const FieldVector<Double>& startSI,
        const FieldVector<Double>& endSI,
        Double durationSI) const;

    Double physicalLower(unsigned int axis) const;
    Double physicalUpper(unsigned int axis) const;

  private:
    Double conductivityRate(const FieldVector<Double>& positionSI) const;

    EBGridGeometry geometry_;
    FieldVector<Double> originSI_;
    EBCPMLParameters parameters_;
    Double outerLower_[3];
    Double outerUpper_[3];
    Double physicalLower_[3];
    Double physicalUpper_[3];
    Double tolerance_[3];
    Double maximumRate_[3];
  };

  /* Minimal non-physical carrier retained only while matched current is
   * attenuated inside CPML. It owns no diagnostic identity or pusher state. */
  struct ParticleCPMLCarrier
  {
    FieldVector<Double> position;
    FieldVector<Double> velocity;
    Double charge;
    Double currentWeight;

    ParticleCPMLCarrier();
  };

  /* Open, absorbing particle boundary for a z-slab decomposition.
   *
   * The in-domain part of an escaping trajectory is still deposited by the
   * charge-conserving current depositor.  The terminal CIC charge is then
   * exported through a virtual normal link on the selected outer face.  That
   * virtual current is outside the Maxwell lattice, so it is accounting for
   * the open-domain continuity equation rather than an additional field
   * source.
   *
   * Production runs retain only six counters and six charge sums.  Dense
   * face-current arrays are optional and intended solely for local continuity
   * diagnostics; disabling them has no per-face allocation or clearing cost.
   */
  class ParticleOpenBoundary
  {
  public:
    ParticleOpenBoundary(const EBGridGeometry& localGeometry,
                         const FieldVector<Double>& localOriginSI,
                         const EBGridGeometry& globalGeometry,
                         const FieldVector<Double>& globalOriginSI,
                         std::size_t localZOffset,
                         bool retainSpatialFlux = false);

    void beginStep();

    bool firstExit(const FieldVector<Double>& startSI,
                   const FieldVector<Double>& endSI,
                   ParticleBoundaryHit& hit) const;

    void depositOutgoingFlux(const ParticleBoundaryHit& hit,
                             Double chargeCoulomb);

    unsigned long long stepCount(ParticleBoundaryFace face) const;
    unsigned long long cumulativeCount(ParticleBoundaryFace face) const;
    Double stepCharge(ParticleBoundaryFace face) const;
    Double cumulativeCharge(ParticleBoundaryFace face) const;
    std::size_t memoryBytes() const;
    bool retainsSpatialFlux() const;

    /* Signed outward current density on a diagnostic face node.  The two
     * indices follow the face tangential axes: (y,z), (x,z), or (x,y).
     * Charge sign is retained; face orientation is already represented by
     * the fact that the value is outward. */
    Double outwardCurrentDensity(ParticleBoundaryFace face,
                                 std::size_t first,
                                 std::size_t second) const;

    /* Open-boundary contribution to div(J) at one local charge vertex. */
    Double outwardDivergence(std::size_t i, std::size_t j,
                             std::size_t k) const;

    Double maxContinuityResidual(const EBFieldGrid& fields,
                                 const YeeComponent& rhoBefore,
                                 const YeeComponent& rhoAfter) const;

  private:
    struct FaceFlux
    {
      std::size_t firstNodes;
      std::size_t secondNodes;
      std::vector<Double> values;

      FaceFlux();
    };

    static std::size_t faceIndex(ParticleBoundaryFace face);
    static unsigned int faceAxis(ParticleBoundaryFace face);
    static bool upperFace(ParticleBoundaryFace face);
    void validateGeometry() const;
    void validateHit(const ParticleBoundaryHit& hit) const;
    void ensureFaceStorage(ParticleBoundaryFace face);
    void depositSpatialFlux(ParticleBoundaryFace face,
                            const FieldVector<Double>& positionSI,
                            Double chargeCoulomb);
    Double faceValue(ParticleBoundaryFace face,
                     std::size_t first, std::size_t second) const;

    EBGridGeometry localGeometry_;
    FieldVector<Double> localOriginSI_;
    EBGridGeometry globalGeometry_;
    FieldVector<Double> globalOriginSI_;
    Double globalLowerSI_[3];
    Double globalUpperSI_[3];
    Double globalToleranceSI_[3];
    Double localLowerZSI_;
    Double localUpperZSI_;
    bool retainSpatialFlux_;
    FaceFlux faceFlux_[6];
    unsigned long long stepCount_[6];
    unsigned long long cumulativeCount_[6];
    Double stepCharge_[6];
    Double cumulativeCharge_[6];
  };
}

#endif
