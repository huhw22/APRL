#include "eb_particles.h"

#include <cmath>
#include <stdexcept>

namespace aprl
{
  RelativisticParticleSI::RelativisticParticleSI()
    : position(0.0), properVelocity(0.0), charge(0.0), mass(0.0),
      weight(1.0), id(0), sourceId(0)
  {}

  void RelativisticBorisPusher::pushMomentum(
      FieldVector<Double>& properVelocity,
      const FieldVector<Double>& electric,
      const FieldVector<Double>& magnetic,
      Double charge, Double mass, Double timeStep)
  {
    if (!(mass > 0.0) || !std::isfinite(mass) ||
        !std::isfinite(charge) || !(timeStep > 0.0) ||
        !std::isfinite(timeStep))
      throw std::invalid_argument("Invalid SI particle pusher parameters");

    const Double electricKick = charge * timeStep /
                                (2.0 * mass * SI::c);
    FieldVector<Double> uMinus(properVelocity);
    uMinus.pmv(electricKick, electric);

    const Double gammaMinus =
      BoostFrameTransform::gammaFromProperVelocity(uMinus);
    FieldVector<Double> rotation(0.0);
    rotation.mv(charge * timeStep / (2.0 * mass * gammaMinus), magnetic);
    FieldVector<Double> uPrime = cross(uMinus, rotation);
    uPrime += uMinus;

    FieldVector<Double> s(rotation);
    s *= 2.0 / (1.0 + rotation.norm2());
    FieldVector<Double> uPlus = cross(uPrime, s);
    uPlus += uMinus;

    properVelocity = uPlus;
    properVelocity.pmv(electricKick, electric);
  }

  void RelativisticBorisPusher::pushPosition(
      FieldVector<Double>& position,
      const FieldVector<Double>& properVelocity,
      Double timeStep)
  {
    const Double gamma =
      BoostFrameTransform::gammaFromProperVelocity(properVelocity);
    position.pmv(SI::c * timeStep / gamma, properVelocity);
  }

  void RelativisticBorisPusher::push(
      RelativisticParticleSI& particle,
      const FieldVector<Double>& electric,
      const FieldVector<Double>& magnetic,
      Double timeStep)
  {
    pushMomentum(particle.properVelocity, electric, magnetic,
                 particle.charge, particle.mass, timeStep);
    pushPosition(particle.position, particle.properVelocity, timeStep);
  }

  void RelativisticBorisPusher::pushFromGrid(
      RelativisticParticleSI& particle,
      const EBFieldGrid& fields,
      const FieldVector<Double>& gridOriginSI,
      Double timeStep)
  {
    FieldVector<Double> electric(0.0);
    FieldVector<Double> magnetic(0.0);
    fields.sampleFieldsPosition(particle.position, gridOriginSI,
                                electric, magnetic);
    push(particle, electric, magnetic, timeStep);
  }

  void RelativisticBorisPusher::pushFromGridAndPrescribedLab(
      RelativisticParticleSI& particle,
      const EBFieldGrid& fields,
      const FieldVector<Double>& gridOriginSI,
      const SIFieldSourceSet& sources,
      const BoostFrameTransform& frame,
      Double timeBoxSI,
      Double timeStep)
  {
    FieldVector<Double> electric(0.0);
    FieldVector<Double> magnetic(0.0);
    fields.sampleFieldsPosition(particle.position, gridOriginSI,
                                electric, magnetic);
    sources.addPrescribedBox(particle.position, timeBoxSI, frame,
                             electric, magnetic);
    push(particle, electric, magnetic, timeStep);
  }

  void RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
      RelativisticParticleSI& particle,
      const EBFieldGrid& fields,
      const FieldVector<Double>& gridOriginSI,
      const SIFieldSourceSet& sources,
      const BoostFrameTransform& frame,
      Double timeBoxSI,
      Double fieldTimeStep,
      unsigned int substeps,
      Double* prescribedWorkJ)
  {
    if (substeps == 0 || !(fieldTimeStep > 0.0) ||
        !std::isfinite(fieldTimeStep))
      throw std::invalid_argument(
        "Particle subcycling needs positive finite dt and substeps");
    FieldVector<Double> gridElectric(0.0);
    FieldVector<Double> gridMagnetic(0.0);
    fields.sampleFieldsPosition(particle.position, gridOriginSI,
                                gridElectric, gridMagnetic);
    const Double particleDt = fieldTimeStep /
      static_cast<Double>(substeps);
    Double work = 0.0;
    for (unsigned int substep = 0; substep < substeps; ++substep)
      {
        const FieldVector<Double> startPosition(particle.position);
        FieldVector<Double> electric(gridElectric);
        FieldVector<Double> magnetic(gridMagnetic);
        FieldVector<Double> prescribedElectric(0.0);
        FieldVector<Double> prescribedMagnetic(0.0);
        const Double substepTime = std::fma(
          static_cast<Double>(substep), particleDt, timeBoxSI);
        sources.addPrescribedBox(startPosition, substepTime, frame,
                                 prescribedElectric, prescribedMagnetic);
        electric += prescribedElectric;
        magnetic += prescribedMagnetic;
        push(particle, electric, magnetic, particleDt);
        if (prescribedWorkJ)
          {
            FieldVector<Double> endElectric(0.0);
            FieldVector<Double> endMagnetic(0.0);
            sources.addPrescribedBox(particle.position,
              substepTime + particleDt, frame,
              endElectric, endMagnetic);
            FieldVector<Double> displacement(particle.position);
            displacement -= startPosition;
            prescribedElectric += endElectric;
            prescribedElectric *= 0.5;
            work = std::fma(particle.charge,
              prescribedElectric * displacement, work);
          }
      }
    if (prescribedWorkJ) *prescribedWorkJ = work;
  }
}
