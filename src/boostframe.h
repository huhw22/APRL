#ifndef BOOSTFRAME_H
#define BOOSTFRAME_H

#include "fieldvector.h"

namespace aprl
{
  class BoostFrameTransform
  {
    private:
      Double gamma_;
      Double beta_;
      Double gammaBeta_;
      Double c0_;
      Double labTimeOrigin_;
      Double labZOrigin_;
      Double boxZOrigin_;

    public:
      BoostFrameTransform();
      BoostFrameTransform(Double gamma, Double beta, Double c0, Double zRef);

      void set(Double gamma, Double beta, Double c0, Double zRef);
      void setFromGamma(Double gamma, Double c0, Double zRef);
      void setOriginsFromGamma(Double gamma, Double c0,
                               Double labTimeOrigin,
                               Double labZOrigin,
                               Double boxZOrigin);

      /* Stable at gamma ~= 1 and avoids the cancellation in
       * sqrt(1 - 1/gamma^2). */
      static Double betaFromGamma(Double gamma);
      static Double gammaBetaFromGamma(Double gamma);

      /* sqrt(1 + |u|^2) evaluated without overflowing when a normalized
       * proper velocity component is very large. */
      static Double gammaFromProperVelocity(
          const FieldVector<Double>& properVelocity);

      /* Longitudinal four-velocity transforms evaluated through light-front
       * components. This avoids subtracting two
       * O(gamma_boost*gamma_beam) terms when the bunch and boost are nearly
       * comoving. */
      void properVelocityLabToBox(
          const FieldVector<Double>& properVelocityLab,
          FieldVector<Double>& properVelocityBox) const;
      void properVelocityBoxToLab(
          const FieldVector<Double>& properVelocityBox,
          FieldVector<Double>& properVelocityLab) const;

      Double gamma() const;
      Double beta() const;
      Double gammaBeta() const;
      Double c0() const;
      Double zRef() const;
      Double labTimeOrigin() const;
      Double labZOrigin() const;
      Double boxZOrigin() const;

      /* Full coordinate transform: laboratory -> computational box. */
      void labToBox(Double tLab, Double zLab, Double& tBox, Double& zBox) const;

      /* Full coordinate transform: computational box -> laboratory. */
      void boxToLab(Double tBox, Double zBox, Double& tLab, Double& zLab) const;

      /* Convenience transforms used by source and output adapters. */
      Double boxZFromLabZAndBoxT(Double zLab, Double tBox) const;
      Double boxTFromLabZT(Double zLab, Double tLab) const;
      Double labZFromBoxZT(Double zBox, Double tBox) const;
      Double labTFromBoxZT(Double zBox, Double tBox) const;

      Double boxZFromLabZT(Double zLab, Double tLab) const;
      Double boxTFromBoxZLabT(Double zBox, Double tLab) const;
      Double labZFromBoxZLabT(Double zBox, Double tLab) const;
      Double labTFromLabZBoxT(Double zLab, Double tBox) const;

      /* Electromagnetic field transforms for a boost along +z.  Inputs and
       * outputs use a single unit system; the E/B core calls these in SI. */
      void fieldsLabToBox(const FieldVector<Double>& electricLab,
                          const FieldVector<Double>& magneticLab,
                          FieldVector<Double>& electricBox,
                          FieldVector<Double>& magneticBox) const;
      void fieldsBoxToLab(const FieldVector<Double>& electricBox,
                          const FieldVector<Double>& magneticBox,
                          FieldVector<Double>& electricLab,
                          FieldVector<Double>& magneticLab) const;
  };
}

#endif
