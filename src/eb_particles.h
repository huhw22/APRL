#ifndef DIRECT_EB_PARTICLES_H
#define DIRECT_EB_PARTICLES_H

#include <cstdint>

#include "eb_field.h"
#include "eb_sources.h"
#include "fieldvector.h"

namespace fel
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

    RelativisticParticleSI();
  };

  class RelativisticBorisPusher
  {
  public:
    static void pushMomentum(
        FieldVector<Double>& properVelocity,
        const FieldVector<Double>& electric,
        const FieldVector<Double>& magnetic,
        Double charge, Double mass, Double timeStep);

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
  };
}

#endif
