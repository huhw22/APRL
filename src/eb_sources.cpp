#include "eb_sources.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace aprl
{
  namespace
  {
    const Double kPi = 3.141592653589793238462643383279502884;
    const Double kTwoPi = 2.0 * kPi;
    const Double kTwoLogTwo = 1.3862943611198906188344642429163531;

    bool finiteVector(const FieldVector<Double>& value)
    {
      return std::isfinite(value[0]) && std::isfinite(value[1]) &&
             std::isfinite(value[2]);
    }

    FieldVector<Double> reversed(const FieldVector<Double>& value)
    {
      FieldVector<Double> result(value);
      result *= -1.0;
      return result;
    }

  }

  SIFieldValue::SIFieldValue() : electric(0.0), magnetic(0.0) {}

  void SIFieldValue::clear()
  {
    electric = FieldVector<Double>(0.0);
    magnetic = FieldVector<Double>(0.0);
  }

  void SIFieldValue::add(const SIFieldValue& other)
  {
    electric += other.electric;
    magnetic += other.magnetic;
  }

  SIWaveEnvelope::SIWaveEnvelope()
    : type(SIEnvelopeType::Gaussian), centerTime(0.0), duration(1.0),
      frequency(1.0), carrierPhase(0.0), risingCycles(1)
  {
    inverseGaussianSigma[0] = 1.0;
    inverseGaussianSigma[1] = 1.0;
  }

  void SIWaveEnvelope::validate() const
  {
    if (!std::isfinite(centerTime) || !(duration > 0.0) ||
        !std::isfinite(duration) || !(frequency > 0.0) ||
        !std::isfinite(frequency) || !std::isfinite(carrierPhase))
      throw std::invalid_argument("Invalid SI wave envelope");
    if ((type == SIEnvelopeType::FlatTop ||
         type == SIEnvelopeType::InverseGaussian) && risingCycles == 0)
      throw std::invalid_argument("Wave-envelope rising cycles must be nonzero");
    if (type == SIEnvelopeType::InverseGaussian &&
        (!(inverseGaussianSigma[0] > 0.0) ||
         !(inverseGaussianSigma[1] > 0.0) ||
         !std::isfinite(inverseGaussianSigma[0]) ||
         !std::isfinite(inverseGaussianSigma[1])))
      throw std::invalid_argument(
          "Inverse-Gaussian sigma values must be positive");
  }

  Double SIWaveEnvelope::value(Double time, Double additionalPhase) const
  {
    const Double delta = time - centerTime;
    if (std::abs(delta) > 10.0 * duration)
      return 0.0;

    const Double carrier =
      std::cos(kTwoPi * frequency * delta + carrierPhase + additionalPhase);
    const Double normalized = delta / duration;
    switch (type)
      {
      case SIEnvelopeType::Neumann:
        /* Dimensionless derivative-Gaussian shape. */
        return -2.0 * kTwoLogTwo * normalized * carrier *
               std::exp(-kTwoLogTwo * normalized * normalized);
      case SIEnvelopeType::Gaussian:
        return carrier * std::exp(-kTwoLogTwo * normalized * normalized);
      case SIEnvelopeType::Secant:
        return carrier / std::cosh(normalized);
      case SIEnvelopeType::FlatTop:
        if (delta < -0.5 * duration)
          return carrier * std::exp(-std::pow(
              (delta + 0.5 * duration) * frequency / risingCycles, 2.0));
        if (delta > 0.5 * duration)
          return carrier * std::exp(-std::pow(
              (delta - 0.5 * duration) * frequency / risingCycles, 2.0));
        return carrier;
      case SIEnvelopeType::InverseGaussian:
        {
          const Double compensation = std::pow(
            (1.0 + std::pow(delta / inverseGaussianSigma[0], 2.0)) *
            (1.0 + std::pow(delta / inverseGaussianSigma[1], 2.0)), 0.25);
          if (delta < -0.5 * duration)
            return carrier * compensation * std::exp(-std::pow(
              (delta + 0.5 * duration) * frequency / risingCycles, 2.0));
          if (delta > 0.5 * duration)
            return carrier * compensation * std::exp(-std::pow(
              (delta - 0.5 * duration) * frequency / risingCycles, 2.0));
          return carrier * compensation;
        }
      }
    return 0.0;
  }

  SIWaveSource::SIWaveSource()
    : profile(SIWaveProfile::Plane), position(0.0), direction(0.0),
      polarization(0.0), peakElectricField(0.0), wavelength(0.0),
      envelope(), prepared_(false)
  {
    direction[2] = 1.0;
    polarization[0] = 1.0;
    radius[0] = radius[1] = 0.0;
    order[0] = order[1] = 0;
    rayleigh_[0] = rayleigh_[1] = 0.0;
  }

  bool SIWaveSource::isGaussianProfile() const
  {
    return profile == SIWaveProfile::Gaussian ||
           profile == SIWaveProfile::StandingGaussian ||
           isSuperGaussianProfile();
  }

  bool SIWaveSource::isSuperGaussianProfile() const
  {
    return profile == SIWaveProfile::SuperGaussian ||
           profile == SIWaveProfile::StandingSuperGaussian;
  }

  bool SIWaveSource::isTruncatedProfile() const
  {
    return profile == SIWaveProfile::TruncatedPlane ||
           profile == SIWaveProfile::StandingTruncatedPlane;
  }

  bool SIWaveSource::isStandingProfile() const
  {
    return profile == SIWaveProfile::StandingPlane ||
           profile == SIWaveProfile::StandingTruncatedPlane ||
           profile == SIWaveProfile::StandingGaussian ||
           profile == SIWaveProfile::StandingSuperGaussian;
  }

  void SIWaveSource::prepare()
  {
    envelope.validate();
    if (!finiteVector(position) || !finiteVector(direction) ||
        !finiteVector(polarization) || !std::isfinite(peakElectricField) ||
        !(wavelength > 0.0) || !std::isfinite(wavelength))
      throw std::invalid_argument("Invalid SI wave source");

    const Double directionNorm = direction.norm();
    const Double polarizationNorm = polarization.norm();
    if (!(directionNorm > 0.0) || !(polarizationNorm > 0.0))
      throw std::invalid_argument(
          "Wave direction and polarization must be nonzero");
    direction /= directionNorm;
    polarization /= polarizationNorm;
    if (std::abs(direction * polarization) > 1.0e-12)
      throw std::invalid_argument(
          "Wave polarization must be perpendicular to its direction");

    if (isGaussianProfile() || isTruncatedProfile())
      {
        if (!(radius[0] > 0.0) || !(radius[1] > 0.0) ||
            !std::isfinite(radius[0]) || !std::isfinite(radius[1]))
          throw std::invalid_argument("Wave radii must be positive");
      }
    if (isSuperGaussianProfile() && (order[0] < 0 || order[1] < 0))
      throw std::invalid_argument("Super-Gaussian orders cannot be negative");

    rayleigh_[0] = kPi * radius[0] * radius[0] / wavelength;
    rayleigh_[1] = kPi * radius[1] * radius[1] / wavelength;
    prepared_ = true;
  }

  void SIWaveSource::fieldsTraveling(
      const FieldVector<Double>& positionLab, Double timeLab,
      const FieldVector<Double>& propagation,
      const FieldVector<Double>& electricDirection,
      SIFieldValue& fields) const
  {
    fields.clear();
    FieldVector<Double> relative(positionLab);
    relative -= position;
    const FieldVector<Double> secondary =
      cross(propagation, electricDirection);
    const Double longitudinal = relative * propagation;
    const Double x = relative * electricDirection;
    const Double y = relative * secondary;
    const Double retardedTime = timeLab - longitudinal / SI::c;

    if (!isGaussianProfile())
      {
        if (isTruncatedProfile() &&
            std::pow(x / radius[0], 2.0) +
            std::pow(y / radius[1], 2.0) > 1.0)
          return;
        const Double signal = envelope.value(retardedTime);
        fields.electric.mv(peakElectricField * signal, electricDirection);
        fields.magnetic.mv(peakElectricField * signal / SI::c, secondary);
        return;
      }

    const Double widthP = std::sqrt(
      1.0 + std::pow(longitudinal / rayleigh_[0], 2.0));
    const Double widthS = std::sqrt(
      1.0 + std::pow(longitudinal / rayleigh_[1], 2.0));
    const Double atanP = std::atan(longitudinal / rayleigh_[0]);
    const Double atanS = std::atan(longitudinal / rayleigh_[1]);
    const Double commonAmplitude = peakElectricField /
      std::sqrt(widthP * widthS);

    const int firstP = isSuperGaussianProfile() ? -order[0] : 0;
    const int lastP  = isSuperGaussianProfile() ?  order[0] : 0;
    const int firstS = isSuperGaussianProfile() ? -order[1] : 0;
    const int lastS  = isSuperGaussianProfile() ?  order[1] : 0;

    for (int i = firstP; i <= lastP; ++i)
      for (int j = firstS; j <= lastS; ++j)
        {
          const Double x0 = (x - i * radius[0]) / widthP;
          const Double y0 = (y - j * radius[1]) / widthS;
          if (std::abs(x0) > 4.0 * radius[0] ||
              std::abs(y0) > 4.0 * radius[1])
            continue;

          const Double phase = 0.5 * (atanP + atanS) -
            kPi * longitudinal / wavelength *
            (std::pow(x0 / rayleigh_[0], 2.0) +
             std::pow(y0 / rayleigh_[1], 2.0));
          const Double amplitude = commonAmplitude *
            std::exp(-std::pow(x0 / radius[0], 2.0) -
                     std::pow(y0 / radius[1], 2.0));

          const Double transverse =
            amplitude * envelope.value(retardedTime, phase - 0.5 * kPi);
          fields.electric.pmv(transverse, electricDirection);
          fields.magnetic.pmv(transverse / SI::c, secondary);

          const Double longitudinalE = amplitude * (-x0 / rayleigh_[0]) *
            envelope.value(retardedTime, phase + atanP);
          const Double longitudinalB = amplitude * (-y0 / rayleigh_[1]) /
            SI::c * envelope.value(retardedTime, phase + atanS);
          fields.electric.pmv(longitudinalE, propagation);
          fields.magnetic.pmv(longitudinalB, propagation);
        }
  }

  void SIWaveSource::fieldsLab(const FieldVector<Double>& positionLab,
                               Double timeLab, SIFieldValue& fields) const
  {
    if (!prepared_)
      throw std::logic_error("SI wave source was not prepared");
    if (!finiteVector(positionLab) || !std::isfinite(timeLab))
      throw std::invalid_argument("Invalid SI wave sample event");

    fieldsTraveling(positionLab, timeLab, direction, polarization, fields);
    if (isStandingProfile())
      {
        SIFieldValue backward;
        fieldsTraveling(positionLab, timeLab, reversed(direction),
                        reversed(polarization), backward);
        fields.add(backward);
      }
  }

  SIMagneticElement::SIMagneticElement()
    : type(SIMagnetType::PlanarUndulator), center(0.0), length(0.0),
      period(0.0), peakMagneticField(0.0), polarizationAngle(0.0),
      gaussianFringe(true), fringeRelativeCutoff(1.0e-9),
      waveNumber_(0.0), fringeExtent_(0.0), cosineAngle_(1.0),
      sineAngle_(0.0), prepared_(false)
  {}

  SIMagneticElement SIMagneticElement::planarUndulatorFromK(
      Double strengthParameter, Double periodSI, Double entranceZSI,
      unsigned int numberOfPeriods, Double polarizationAngleRadians)
  {
    if (!std::isfinite(strengthParameter))
      throw std::invalid_argument("Undulator K must be finite");
    SIMagneticElement result;
    result.type = SIMagnetType::PlanarUndulator;
    result.center[2] = entranceZSI;
    result.period = periodSI;
    result.length = static_cast<Double>(numberOfPeriods) * periodSI;
    result.peakMagneticField = strengthParameter * kTwoPi *
      SI::electronMass * SI::c / (SI::elementaryCharge * periodSI);
    result.polarizationAngle = polarizationAngleRadians;
    result.prepare();
    return result;
  }

  void SIMagneticElement::prepare()
  {
    if (!finiteVector(center) || !(length >= 0.0) ||
        !std::isfinite(length) || !std::isfinite(peakMagneticField) ||
        !std::isfinite(polarizationAngle))
      throw std::invalid_argument("Invalid SI magnetic element");
    if (type == SIMagnetType::PlanarUndulator &&
        (!(period > 0.0) || !std::isfinite(period)))
      throw std::invalid_argument("Undulator period must be positive");
    if (type == SIMagnetType::PlanarUndulator && gaussianFringe &&
        (!(fringeRelativeCutoff > 0.0) ||
         !(fringeRelativeCutoff < 1.0) ||
         !std::isfinite(fringeRelativeCutoff)))
      throw std::invalid_argument(
        "Undulator fringe_relative_cutoff must be in (0,1)");
    if (type == SIMagnetType::UniformDipole && !(length > 0.0))
      throw std::invalid_argument("Dipole length must be positive");

    waveNumber_ = type == SIMagnetType::PlanarUndulator
      ? kTwoPi / period : 0.0;
    fringeExtent_ = type == SIMagnetType::PlanarUndulator && gaussianFringe
      ? std::sqrt(-2.0 * std::log(fringeRelativeCutoff)) / waveNumber_
      : 0.0;
    cosineAngle_ = std::cos(polarizationAngle);
    sineAngle_ = std::sin(polarizationAngle);
    prepared_ = true;
  }

  void SIMagneticElement::fieldsLab(
      const FieldVector<Double>& positionLab, Double,
      SIFieldValue& fields) const
  {
    if (!prepared_)
      throw std::logic_error("SI magnetic element was not prepared");
    fields.clear();

    FieldVector<Double> relative(positionLab);
    relative -= center;
    const Double z = relative[2];
    if (type == SIMagnetType::UniformDipole)
      {
        if (z >= 0.0 && z <= length)
          {
            fields.magnetic[0] = peakMagneticField * cosineAngle_;
            fields.magnetic[1] = peakMagneticField * sineAngle_;
          }
        return;
      }

    Double transverseShape = 0.0;
    Double longitudinalShape = 0.0;
    if (z >= 0.0 && z <= length)
      {
        transverseShape = std::sin(waveNumber_ * z);
        longitudinalShape = std::cos(waveNumber_ * z);
      }
    else if (gaussianFringe)
      {
        const Double fringeCoordinate = z < 0.0 ? z : z - length;
        const Double distance = std::abs(fringeCoordinate);
        if (distance >= fringeExtent_) return;

        /* The legacy Gaussian end field had infinite support (and therefore
         * no usable interaction boundary).  A quintic compact-support taper
         * makes both g and dg/dz vanish at the declared edge.  Defining the
         * transverse shape from -g'/k preserves div(B)=0 through the taper. */
        const Double q = distance / fringeExtent_;
        const Double q2 = q * q;
        const Double q3 = q2 * q;
        const Double q4 = q3 * q;
        const Double q5 = q4 * q;
        const Double taper = 1.0 - 10.0 * q3 + 15.0 * q4 - 6.0 * q5;
        const Double taperDerivative =
          (-30.0 * q2 + 60.0 * q3 - 30.0 * q4) / fringeExtent_;
        const Double gaussian = std::exp(-0.5 *
          std::pow(waveNumber_ * distance, 2.0));
        longitudinalShape = gaussian * taper;
        const Double derivativeByDistance = gaussian *
          (taperDerivative - waveNumber_ * waveNumber_ * distance * taper);
        transverseShape = z < 0.0
          ? derivativeByDistance / waveNumber_
          : -derivativeByDistance / waveNumber_;
      }
    else
      return;

    const Double transverseCoordinate =
      relative[0] * cosineAngle_ + relative[1] * sineAngle_;
    const Double transverseField = peakMagneticField *
      std::cosh(waveNumber_ * transverseCoordinate) * transverseShape;
    fields.magnetic[0] = transverseField * cosineAngle_;
    fields.magnetic[1] = transverseField * sineAngle_;
    fields.magnetic[2] = peakMagneticField *
      std::sinh(waveNumber_ * transverseCoordinate) * longitudinalShape;
  }

  Double SIMagneticElement::physicalEntranceLab() const
  {
    if (!prepared_) throw std::logic_error("Magnetic element was not prepared");
    return center[2];
  }

  Double SIMagneticElement::physicalExitLab() const
  {
    if (!prepared_) throw std::logic_error("Magnetic element was not prepared");
    return center[2] + length;
  }

  Double SIMagneticElement::interactionEntranceLab() const
  {
    if (!prepared_) throw std::logic_error("Magnetic element was not prepared");
    return center[2] - fringeExtent_;
  }

  Double SIMagneticElement::interactionExitLab() const
  {
    if (!prepared_) throw std::logic_error("Magnetic element was not prepared");
    return center[2] + length + fringeExtent_;
  }

  void SIFieldSourceSet::addWave(const SIWaveSource& source,
                                 SIFieldEvolution evolution)
  {
    SIWaveSource prepared(source);
    prepared.prepare();
    if (evolution == SIFieldEvolution::MaxwellIncident)
      maxwellIncidentWaves_.push_back(prepared);
    else
      prescribedLabWaves_.push_back(prepared);
  }

  void SIFieldSourceSet::addPrescribedLabMagnet(
      const SIMagneticElement& source)
  {
    SIMagneticElement prepared(source);
    prepared.prepare();
    prescribedLabMagnets_.push_back(prepared);
  }

  std::size_t SIFieldSourceSet::maxwellIncidentWaveCount() const
  {
    return maxwellIncidentWaves_.size();
  }

  std::size_t SIFieldSourceSet::prescribedLabWaveCount() const
  {
    return prescribedLabWaves_.size();
  }

  std::size_t SIFieldSourceSet::prescribedLabMagnetCount() const
  {
    return prescribedLabMagnets_.size();
  }

  void SIFieldSourceSet::sampleMaxwellIncidentLab(
      const FieldVector<Double>& positionLab, Double timeLab,
      SIFieldValue& fields) const
  {
    fields.clear();
    SIFieldValue contribution;
    for (std::size_t i = 0; i < maxwellIncidentWaves_.size(); ++i)
      {
        maxwellIncidentWaves_[i].fieldsLab(
            positionLab, timeLab, contribution);
        fields.add(contribution);
      }
  }

  void SIFieldSourceSet::samplePrescribedLab(
      const FieldVector<Double>& positionLab, Double timeLab,
      SIFieldValue& fields) const
  {
    fields.clear();
    SIFieldValue contribution;
    for (std::size_t i = 0; i < prescribedLabWaves_.size(); ++i)
      {
        prescribedLabWaves_[i].fieldsLab(positionLab, timeLab, contribution);
        fields.add(contribution);
      }
    for (std::size_t i = 0; i < prescribedLabMagnets_.size(); ++i)
      {
        prescribedLabMagnets_[i].fieldsLab(
            positionLab, timeLab, contribution);
        fields.add(contribution);
      }
  }

  void SIFieldSourceSet::transformSampleToBox(
      const FieldVector<Double>& positionBox, Double timeBox,
      const BoostFrameTransform& frame, bool incident,
      SIFieldValue& fields) const
  {
    FieldVector<Double> positionLab(positionBox);
    Double timeLab = 0.0;
    frame.boxToLab(timeBox, positionBox[2], timeLab, positionLab[2]);

    SIFieldValue lab;
    if (incident)
      sampleMaxwellIncidentLab(positionLab, timeLab, lab);
    else
      samplePrescribedLab(positionLab, timeLab, lab);
    frame.fieldsLabToBox(lab.electric, lab.magnetic,
                         fields.electric, fields.magnetic);
  }

  void SIFieldSourceSet::sampleMaxwellIncidentBox(
      const FieldVector<Double>& positionBox, Double timeBox,
      const BoostFrameTransform& frame, SIFieldValue& fields) const
  {
    transformSampleToBox(positionBox, timeBox, frame, true, fields);
  }

  void SIFieldSourceSet::samplePrescribedBox(
      const FieldVector<Double>& positionBox, Double timeBox,
      const BoostFrameTransform& frame, SIFieldValue& fields) const
  {
    transformSampleToBox(positionBox, timeBox, frame, false, fields);
  }

  void SIFieldSourceSet::addPrescribedBox(
      const FieldVector<Double>& positionBox, Double timeBox,
      const BoostFrameTransform& frame,
      FieldVector<Double>& electricBox,
      FieldVector<Double>& magneticBox) const
  {
    SIFieldValue external;
    samplePrescribedBox(positionBox, timeBox, frame, external);
    electricBox += external.electric;
    magneticBox += external.magnetic;
  }

}
