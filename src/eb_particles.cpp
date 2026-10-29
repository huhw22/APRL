#include "eb_particles.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace aprl
{
  RelativisticParticleSI::RelativisticParticleSI()
    : position(0.0), properVelocity(0.0), charge(0.0), mass(0.0),
      weight(1.0), id(0), sourceId(0), frozen(false), frozenBetaZ(0.0)
  {}

  BorisWork::BorisWork() { clear(); }

  void BorisWork::clear()
  {
    gridWork = prescribedWork = magneticRotationDeltaK = 0.0L;
    totalDeltaK = roundoffResidual = 0.0L;
    absoluteWork = 0.0L;
    for (unsigned int channel = 0; channel < 6; ++channel)
      compensation[channel] = 0.0L;
  }

  void BorisWork::add(const BorisWork& other, long double fraction)
  {
    long double* destination[6] = {&gridWork, &prescribedWork,
      &magneticRotationDeltaK, &totalDeltaK, &roundoffResidual, &absoluteWork};
    const long double increment[6] = {fraction * other.gridWork,
      fraction * other.prescribedWork, fraction * other.magneticRotationDeltaK,
      fraction * other.totalDeltaK, fraction * other.roundoffResidual,
      std::abs(fraction) * other.absoluteWork};
    for (unsigned int channel = 0; channel < 6; ++channel)
      {
        const long double corrected = increment[channel] - compensation[channel];
        const long double next = *destination[channel] + corrected;
        compensation[channel] = (next - *destination[channel]) - corrected;
        *destination[channel] = next;
      }
  }

  namespace
  {
    long double auditGamma(const FieldVector<Double>& u)
    {
      long double norm = 1.0L;
      for (unsigned int axis = 0; axis < 3; ++axis)
        norm += static_cast<long double>(u[axis]) * u[axis];
      return std::sqrt(norm);
    }

    long double auditDeltaK(const FieldVector<Double>& before,
                            const FieldVector<Double>& after,
                            long double restOverGammaSum)
    {
      long double difference = 0.0L;
      for (unsigned int axis = 0; axis < 3; ++axis)
        difference +=
          (static_cast<long double>(after[axis]) - before[axis]) *
          (static_cast<long double>(after[axis]) + before[axis]);
      return restOverGammaSum * difference;
    }

    long double auditElectricKick(const FieldVector<Double>& before,
                                  const FieldVector<Double>& after,
                                  const FieldVector<Double>& electric,
                                  Double electricKick,
                                  long double restOverGammaSum)
    {
      long double dot = 0.0L;
      for (unsigned int axis = 0; axis < 3; ++axis)
        dot += (static_cast<long double>(after[axis]) + before[axis]) *
          electric[axis];
      return restOverGammaSum * static_cast<long double>(electricKick) * dot;
    }

    void pushSubstepWithWork(RelativisticParticleSI& particle,
                             const FieldVector<Double>& electric,
                             const FieldVector<Double>& magnetic,
                             const FieldVector<Double>& gridElectric,
                             const FieldVector<Double>& prescribedElectric,
                             Double dt, BorisWork* work)
    {
      if (work)
        {
          RelativisticBorisPusher::pushMomentumWithWork(
            particle.properVelocity, electric, magnetic,
            gridElectric, prescribedElectric,
            particle.charge, particle.mass, dt, *work);
          RelativisticBorisPusher::pushPosition(
            particle.position, particle.properVelocity, dt);
        }
      else
        RelativisticBorisPusher::push(particle, electric, magnetic, dt);
    }
  }

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

  void RelativisticBorisPusher::pushMomentumWithWork(
      FieldVector<Double>& properVelocity,
      const FieldVector<Double>& totalElectric,
      const FieldVector<Double>& totalMagnetic,
      const FieldVector<Double>& gridElectric,
      const FieldVector<Double>& prescribedElectric,
      Double charge, Double mass, Double timeStep,
      BorisWork& work)
  {
    if (!(mass > 0.0) || !std::isfinite(mass) ||
        !std::isfinite(charge) || !(timeStep > 0.0) ||
        !std::isfinite(timeStep))
      throw std::invalid_argument("Invalid SI particle pusher parameters");

    // Keep every state-producing expression in the same order as the
    // original pushMomentum. All extended-precision work below is a reader
    // of these states and never feeds the push or field sampling.
    const FieldVector<Double> before(properVelocity);
    const Double electricKick = charge * timeStep /
                                (2.0 * mass * SI::c);
    FieldVector<Double> uMinus(properVelocity);
    uMinus.pmv(electricKick, totalElectric);

    const Double gammaMinus =
      BoostFrameTransform::gammaFromProperVelocity(uMinus);
    FieldVector<Double> rotation(0.0);
    rotation.mv(charge * timeStep / (2.0 * mass * gammaMinus), totalMagnetic);
    FieldVector<Double> uPrime = cross(uMinus, rotation);
    uPrime += uMinus;

    FieldVector<Double> s(rotation);
    s *= 2.0 / (1.0 + rotation.norm2());
    FieldVector<Double> uPlus = cross(uPrime, s);
    uPlus += uMinus;

    properVelocity = uPlus;
    properVelocity.pmv(electricKick, totalElectric);

    work.clear();
    const long double c = SI::c;
    const long double rest = static_cast<long double>(mass) * c * c;
    const long double gammaBefore = auditGamma(before);
    const long double gammaMinusAudit = auditGamma(uMinus);
    const long double gammaPlus = auditGamma(uPlus);
    const long double gammaAfter = auditGamma(properVelocity);
    const long double firstKickScale = rest / (gammaBefore + gammaMinusAudit);
    const long double secondKickScale = rest / (gammaPlus + gammaAfter);
    work.gridWork = auditElectricKick(before, uMinus,
        gridElectric, electricKick, firstKickScale) +
      auditElectricKick(uPlus, properVelocity,
        gridElectric, electricKick, secondKickScale);
    work.prescribedWork = auditElectricKick(before, uMinus,
        prescribedElectric, electricKick, firstKickScale) +
      auditElectricKick(uPlus, properVelocity,
        prescribedElectric, electricKick, secondKickScale);
    work.magneticRotationDeltaK = auditDeltaK(uMinus, uPlus,
      rest / (gammaMinusAudit + gammaPlus));
    work.totalDeltaK = auditDeltaK(before, properVelocity,
      rest / (gammaBefore + gammaAfter));
    work.roundoffResidual = work.totalDeltaK - work.gridWork -
      work.prescribedWork - work.magneticRotationDeltaK;
    work.absoluteWork = std::abs(work.gridWork) + std::abs(work.prescribedWork) +
      std::abs(work.magneticRotationDeltaK) + std::abs(work.totalDeltaK) +
      std::abs(work.roundoffResidual);
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
      Double* prescribedWorkJ,
      bool includeGridFields,
      BorisWork* borisWork,
      ParticleSubstepObserver* observer)
  {
    if (substeps == 0 || !(fieldTimeStep > 0.0) ||
        !std::isfinite(fieldTimeStep))
      throw std::invalid_argument(
        "Particle subcycling needs positive finite dt and substeps");
    FieldVector<Double> gridElectric(0.0);
    FieldVector<Double> gridMagnetic(0.0);
    if (includeGridFields)
      {
        try
          {
            fields.sampleFieldsPosition(particle.position, gridOriginSI,
                                        gridElectric, gridMagnetic);
          }
        catch (const std::out_of_range& error)
          {
            throw std::out_of_range(std::string("particle grid-field push: ") +
                                    error.what());
          }
      }
    const Double particleDt = fieldTimeStep /
      static_cast<Double>(substeps);
    if (borisWork) borisWork->clear();
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
        if (borisWork || observer)
          {
            BorisWork substepWork;
            BorisWork* collected = borisWork ? &substepWork : 0;
            if (observer)
              {
                const RelativisticParticleSI before(particle);
                pushSubstepWithWork(particle, electric, magnetic,
                  gridElectric, prescribedElectric, particleDt, collected);
                observer->observe(before, particle, substepTime,
                  substepTime + particleDt, collected);
              }
            else
              pushSubstepWithWork(particle, electric, magnetic,
                gridElectric, prescribedElectric, particleDt, collected);
            if (borisWork) borisWork->add(substepWork);
          }
        else
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
