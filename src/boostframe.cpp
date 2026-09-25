#include "boostframe.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace aprl
{
  BoostFrameTransform::BoostFrameTransform()
    : gamma_(1.0), beta_(0.0), gammaBeta_(0.0), c0_(1.0),
      labTimeOrigin_(0.0), labZOrigin_(0.0), boxZOrigin_(0.0)
  {}

  BoostFrameTransform::BoostFrameTransform(Double gamma, Double beta,
                                                 Double c0, Double zRef)
    : gamma_(1.0), beta_(0.0), gammaBeta_(0.0), c0_(1.0),
      labTimeOrigin_(0.0), labZOrigin_(0.0), boxZOrigin_(0.0)
  {
    set(gamma, beta, c0, zRef);
  }

  void BoostFrameTransform::set(Double gamma, Double beta,
                                   Double c0, Double zRef)
  {
    if (!(gamma >= 1.0) || !std::isfinite(gamma))
      throw std::invalid_argument("Lorentz gamma must be finite and >= 1");
    if (!std::isfinite(beta) || std::abs(beta) >= 1.0)
      throw std::invalid_argument("Lorentz beta must be finite and |beta| < 1");
    if (!(c0 > 0.0) || !std::isfinite(c0))
      throw std::invalid_argument("Transform light speed must be positive");

    gamma_ = gamma;
    beta_  = beta;
    gammaBeta_ = gamma * beta;
    c0_    = c0;
    labTimeOrigin_ = 0.0;
    labZOrigin_ = 0.0;
    boxZOrigin_ = zRef;
  }

  void BoostFrameTransform::setFromGamma(Double gamma, Double c0,
                                         Double zRef)
  {
    set(gamma, betaFromGamma(gamma), c0, zRef);
    gammaBeta_ = gammaBetaFromGamma(gamma);
  }

  void BoostFrameTransform::setOriginsFromGamma(
      Double gamma, Double c0, Double labTimeOrigin,
      Double labZOrigin, Double boxZOrigin)
  {
    setFromGamma(gamma, c0, boxZOrigin);
    if (!std::isfinite(labTimeOrigin) || !std::isfinite(labZOrigin) ||
        !std::isfinite(boxZOrigin))
      throw std::invalid_argument("Lorentz origins must be finite");
    labTimeOrigin_ = labTimeOrigin;
    labZOrigin_ = labZOrigin;
    boxZOrigin_ = boxZOrigin;
  }

  Double BoostFrameTransform::betaFromGamma(Double gamma)
  {
    return gammaBetaFromGamma(gamma) / gamma;
  }

  Double BoostFrameTransform::gammaBetaFromGamma(Double gamma)
  {
    if (!(gamma >= 1.0) || !std::isfinite(gamma))
      throw std::invalid_argument("Lorentz gamma must be finite and >= 1");
    return std::sqrt(gamma - 1.0) * std::sqrt(gamma + 1.0);
  }

  Double BoostFrameTransform::gammaFromProperVelocity(
      const FieldVector<Double>& properVelocity)
  {
    const Double scale = std::max(
        1.0, std::max(std::abs(properVelocity[0]),
                      std::max(std::abs(properVelocity[1]),
                               std::abs(properVelocity[2]))));
    const Double inverseScale = 1.0 / scale;
    const Double ux = properVelocity[0] * inverseScale;
    const Double uy = properVelocity[1] * inverseScale;
    const Double uz = properVelocity[2] * inverseScale;
    return scale * std::sqrt(inverseScale * inverseScale +
                             ux * ux + uy * uy + uz * uz);
  }

  namespace
  {
    Double positiveLightFrontComponent(
        const FieldVector<Double>& properVelocity)
    {
      const Double gamma =
        BoostFrameTransform::gammaFromProperVelocity(properVelocity);
      const Double transverseMassSquared = std::fma(
        properVelocity[0], properVelocity[0],
        std::fma(properVelocity[1], properVelocity[1], 1.0));
      if (!(transverseMassSquared > 0.0) ||
          !std::isfinite(transverseMassSquared))
        throw std::overflow_error(
          "Cannot form a finite transverse mass for Lorentz transformation");
      if (properVelocity[2] >= 0.0)
        return gamma + properVelocity[2];
      return transverseMassSquared / (gamma - properVelocity[2]);
    }

    void properVelocityFromPositiveLightFront(
        const FieldVector<Double>& input, Double positive,
        FieldVector<Double>& output)
    {
      const Double transverseMassSquared = std::fma(
        input[0], input[0],
        std::fma(input[1], input[1], 1.0));
      if (!(positive > 0.0) || !std::isfinite(positive))
        throw std::overflow_error(
          "Lorentz transformation produced an invalid light-front momentum");
      const Double negative = transverseMassSquared / positive;
      output[0] = input[0];
      output[1] = input[1];
      output[2] = 0.5 * (positive - negative);
      if (!std::isfinite(output[2]))
        throw std::overflow_error(
          "Lorentz transformation produced a non-finite longitudinal momentum");
    }
  }

  void BoostFrameTransform::properVelocityLabToBox(
      const FieldVector<Double>& properVelocityLab,
      FieldVector<Double>& properVelocityBox) const
  {
    const Double rapidityFactor = gamma_ + gammaBeta_;
    const Double positiveBox =
      positiveLightFrontComponent(properVelocityLab) / rapidityFactor;
    properVelocityFromPositiveLightFront(
      properVelocityLab, positiveBox, properVelocityBox);
  }

  void BoostFrameTransform::properVelocityBoxToLab(
      const FieldVector<Double>& properVelocityBox,
      FieldVector<Double>& properVelocityLab) const
  {
    const Double rapidityFactor = gamma_ + gammaBeta_;
    const Double positiveLab =
      positiveLightFrontComponent(properVelocityBox) * rapidityFactor;
    properVelocityFromPositiveLightFront(
      properVelocityBox, positiveLab, properVelocityLab);
  }

  Double BoostFrameTransform::gamma() const { return gamma_; }
  Double BoostFrameTransform::beta()  const { return beta_;  }
  Double BoostFrameTransform::gammaBeta() const { return gammaBeta_; }
  Double BoostFrameTransform::c0()    const { return c0_;    }
  Double BoostFrameTransform::zRef()  const { return boxZOrigin_; }
  Double BoostFrameTransform::labTimeOrigin() const
  { return labTimeOrigin_; }
  Double BoostFrameTransform::labZOrigin() const { return labZOrigin_; }
  Double BoostFrameTransform::boxZOrigin() const { return boxZOrigin_; }

  void BoostFrameTransform::labToBox(Double tLab, Double zLab,
                                        Double& tBox, Double& zBox) const
  {
    const Double time = tLab - labTimeOrigin_;
    const Double position = zLab - labZOrigin_;
    tBox = std::fma(-gammaBeta_ / c0_, position, gamma_ * time);
    zBox = boxZOrigin_ +
      std::fma(-gammaBeta_ * c0_, time, gamma_ * position);
  }

  void BoostFrameTransform::boxToLab(Double tBox, Double zBox,
                                        Double& tLab, Double& zLab) const
  {
    const Double zShift = zBox - boxZOrigin_;
    zLab = labZOrigin_ +
      std::fma(gammaBeta_ * c0_, tBox, gamma_ * zShift);
    tLab = labTimeOrigin_ +
      std::fma(gammaBeta_ / c0_, zShift, gamma_ * tBox);
  }

  Double BoostFrameTransform::boxZFromLabZAndBoxT(Double zLab, Double tBox) const
  {
    return boxZOrigin_ + (zLab - labZOrigin_) / gamma_ -
      beta_ * c0_ * tBox;
  }

  Double BoostFrameTransform::boxTFromLabZT(Double zLab, Double tLab) const
  {
    return std::fma(-gammaBeta_ / c0_, zLab - labZOrigin_,
                    gamma_ * (tLab - labTimeOrigin_));
  }

  Double BoostFrameTransform::labZFromBoxZT(Double zBox, Double tBox) const
  {
    return labZOrigin_ + std::fma(
      gammaBeta_ * c0_, tBox, gamma_ * (zBox - boxZOrigin_));
  }

  Double BoostFrameTransform::labTFromBoxZT(Double zBox, Double tBox) const
  {
    return labTimeOrigin_ + std::fma(
      gammaBeta_ / c0_, zBox - boxZOrigin_, gamma_ * tBox);
  }

  Double BoostFrameTransform::boxZFromLabZT(Double zLab, Double tLab) const
  {
    return boxZOrigin_ + std::fma(
      -gammaBeta_ * c0_, tLab - labTimeOrigin_,
      gamma_ * (zLab - labZOrigin_));
  }

  Double BoostFrameTransform::boxTFromBoxZLabT(Double zBox, Double tLab) const
  {
    return (tLab - labTimeOrigin_) / gamma_ -
      beta_ * (zBox - boxZOrigin_) / c0_;
  }

  Double BoostFrameTransform::labZFromBoxZLabT(Double zBox, Double tLab) const
  {
    return labZOrigin_ + (zBox - boxZOrigin_) / gamma_ +
      beta_ * c0_ * (tLab - labTimeOrigin_);
  }

  Double BoostFrameTransform::labTFromLabZBoxT(Double zLab, Double tBox) const
  {
    return labTimeOrigin_ + tBox / gamma_ +
      beta_ * (zLab - labZOrigin_) / c0_;
  }

  void BoostFrameTransform::fieldsLabToBox(
      const FieldVector<Double>& electricLab,
      const FieldVector<Double>& magneticLab,
      FieldVector<Double>& electricBox,
      FieldVector<Double>& magneticBox) const
  {
    electricBox[0] = std::fma(-gammaBeta_ * c0_, magneticLab[1],
                              gamma_ * electricLab[0]);
    electricBox[1] = std::fma( gammaBeta_ * c0_, magneticLab[0],
                              gamma_ * electricLab[1]);
    electricBox[2] = electricLab[2];

    magneticBox[0] = std::fma( gammaBeta_ / c0_, electricLab[1],
                              gamma_ * magneticLab[0]);
    magneticBox[1] = std::fma(-gammaBeta_ / c0_, electricLab[0],
                              gamma_ * magneticLab[1]);
    magneticBox[2] = magneticLab[2];
  }

  void BoostFrameTransform::fieldsBoxToLab(
      const FieldVector<Double>& electricBox,
      const FieldVector<Double>& magneticBox,
      FieldVector<Double>& electricLab,
      FieldVector<Double>& magneticLab) const
  {
    electricLab[0] = std::fma( gammaBeta_ * c0_, magneticBox[1],
                              gamma_ * electricBox[0]);
    electricLab[1] = std::fma(-gammaBeta_ * c0_, magneticBox[0],
                              gamma_ * electricBox[1]);
    electricLab[2] = electricBox[2];

    magneticLab[0] = std::fma(-gammaBeta_ / c0_, electricBox[1],
                              gamma_ * magneticBox[0]);
    magneticLab[1] = std::fma( gammaBeta_ / c0_, electricBox[0],
                              gamma_ * magneticBox[1]);
    magneticLab[2] = magneticBox[2];
  }
}
