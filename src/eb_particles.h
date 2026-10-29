#ifndef DIRECT_EB_PARTICLES_H
#define DIRECT_EB_PARTICLES_H

#include <cstdint>

#include "eb_field.h"
#include "eb_sources.h"
#include "fieldvector.h"

namespace aprl
{
  struct RelativisticParticleSI
  {
    FieldVector<Double> position;       /* m */
    FieldVector<Double> properVelocity; /* gamma * v/c */
    Double charge;                      /* C */
    Double mass;                        /* kg */
    Double weight;                      /* dimensionless diagnostic weight */
    std::uint64_t id;                   /* globally unique runtime id */
    std::uint64_t sourceId;             /* optional input-record id */
    bool frozen;                        /* upstream lab-plane source state */
    Double frozenBetaZ;                 /* cached held speed; valid while frozen */

    RelativisticParticleSI();
  };

  /* Read-only work accounting for the actual two electric kicks and magnetic
   * rotation of a Boris update. All fields are represented-particle energies
   * in joules; charge/mass already include the macro-particle weight.
   * roundoffResidual is a local arithmetic remainder, not a field/PIC energy
   * conservation residual. Source work is evaluated independently of it. */
  struct BorisWork
  {
    long double gridWork;
    long double prescribedWork;
    long double magneticRotationDeltaK;
    long double totalDeltaK;
    long double roundoffResidual;
    // Sum of absolute leaf work channels, propagated through all accumulation
    // levels (and true terminal fractions), never abs of the net history.
    long double absoluteWork;
    // Kahan state in the order above: grid, prescribed, magnetic, total,
    // roundoff, absolute. Persist together with the visible sums on restart.
    long double compensation[6];

    BorisWork();
    void clear();
    void add(const BorisWork& other, long double fraction = 1.0L);
  };

  /* Called only when explicitly supplied. The states and optional work are
   * valid only for this call; implementations that defer output must copy
   * their needed data. Parallel particle preparation requires an observer
   * private to the prepared particle, with no MPI/HDF5 from this callback. */
  class ParticleSubstepObserver
  {
  public:
    virtual ~ParticleSubstepObserver() {}
    virtual void observe(const RelativisticParticleSI& before,
                         const RelativisticParticleSI& after,
                         Double timeBefore, Double timeAfter,
                         const BorisWork* work) = 0;
  };

  class RelativisticBorisPusher
  {
  public:
    static void pushMomentum(
        FieldVector<Double>& properVelocity,
        const FieldVector<Double>& electric,
        const FieldVector<Double>& magnetic,
        Double charge, Double mass, Double timeStep);

    /* Same state arithmetic as pushMomentum. totalElectric must be the same
     * already-combined value used by the ordinary pusher; the two components
     * are used only to attribute the electric work. Overwrites work. */
    static void pushMomentumWithWork(
        FieldVector<Double>& properVelocity,
        const FieldVector<Double>& totalElectric,
        const FieldVector<Double>& totalMagnetic,
        const FieldVector<Double>& gridElectric,
        const FieldVector<Double>& prescribedElectric,
        Double charge, Double mass, Double timeStep,
        BorisWork& work);

    static void pushPosition(FieldVector<Double>& position,
                             const FieldVector<Double>& properVelocity,
                             Double timeStep);

    static void push(RelativisticParticleSI& particle,
                     const FieldVector<Double>& electric,
                     const FieldVector<Double>& magnetic,
                     Double timeStep);

    static void pushFromGrid(RelativisticParticleSI& particle,
                             const EBFieldGrid& fields,
                             const FieldVector<Double>& gridOriginSI,
                             Double timeStep);

    /* Add only PrescribedLab device fields to the self-consistent Maxwell
     * field at the particle event.  MaxwellIncident waves are deliberately
     * excluded because they already live on the grid. */
    static void pushFromGridAndPrescribedLab(
        RelativisticParticleSI& particle,
        const EBFieldGrid& fields,
        const FieldVector<Double>& gridOriginSI,
        const SIFieldSourceSet& sources,
        const BoostFrameTransform& frame,
        Double timeBoxSI,
        Double timeStep);

    /* Hold the staggered Maxwell sample fixed over one field step while
     * resolving PrescribedLab devices with smaller Boris steps. This keeps
     * MPI ownership independent of the substep count; Maxwell current and
     * detector cadence remain tied to the enclosing field step. */
    static void pushFromGridAndPrescribedLabSubcycled(
        RelativisticParticleSI& particle,
        const EBFieldGrid& fields,
        const FieldVector<Double>& gridOriginSI,
        const SIFieldSourceSet& sources,
        const BoostFrameTransform& frame,
        Double timeBoxSI,
        Double fieldTimeStep,
        unsigned int substeps,
        Double* prescribedWorkJ = 0,
        bool includeGridFields = true,
        BorisWork* borisWork = 0,
        ParticleSubstepObserver* observer = 0);
  };
}

#endif
