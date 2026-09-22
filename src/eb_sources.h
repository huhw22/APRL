#ifndef DIRECT_EB_SOURCES_H
#define DIRECT_EB_SOURCES_H

#include <cstddef>
#include <vector>

#include "eb_field.h"
#include "fieldvector.h"

namespace fel
{
  struct SIFieldValue
  {
    FieldVector<Double> electric; /* V/m */
    FieldVector<Double> magnetic; /* T */

    SIFieldValue();
    void clear();
    void add(const SIFieldValue& other);
  };

  enum class SIEnvelopeType
  {
    Neumann,
    Gaussian,
    Secant,
    FlatTop,
    InverseGaussian
  };

  /* Time-domain carrier and envelope in seconds and hertz.  value() is
   * dimensionless, so every wave source has one unambiguous SI E amplitude. */
  struct SIWaveEnvelope
  {
    SIEnvelopeType type;
    Double centerTime;          /* s */
    Double duration;            /* s, intensity FWHM for Gaussian */
    Double frequency;           /* Hz */
    Double carrierPhase;        /* rad */
    unsigned int risingCycles;
    Double inverseGaussianSigma[2]; /* s */

    SIWaveEnvelope();
    void validate() const;
    Double value(Double time, Double additionalPhase = 0.0) const;
  };

  enum class SIWaveProfile
  {
    Plane,
    TruncatedPlane,
    Gaussian,
    SuperGaussian,
    StandingPlane,
    StandingTruncatedPlane,
    StandingGaussian,
    StandingSuperGaussian
  };

  /* Prescribed laboratory-frame electromagnetic wave.  prepare() validates
   * and precomputes all geometry; evaluation then performs no allocation and
   * no unit conversion. */
  class SIWaveSource
  {
  public:
    SIWaveSource();

    SIWaveProfile profile;
    FieldVector<Double> position;     /* m, focus/reference point */
    FieldVector<Double> direction;    /* unit propagation vector */
    FieldVector<Double> polarization; /* unit E-polarization vector */
    Double peakElectricField;         /* V/m */
    Double wavelength;                /* m */
    Double radius[2];                  /* m */
    int order[2];
    SIWaveEnvelope envelope;

    void prepare();
    void fieldsLab(const FieldVector<Double>& positionLab,
                   Double timeLab, SIFieldValue& fields) const;

  private:
    void fieldsTraveling(const FieldVector<Double>& positionLab,
                         Double timeLab,
                         const FieldVector<Double>& propagation,
                         const FieldVector<Double>& electricDirection,
                         SIFieldValue& fields) const;
    bool isGaussianProfile() const;
    bool isSuperGaussianProfile() const;
    bool isTruncatedProfile() const;
    bool isStandingProfile() const;

    Double rayleigh_[2];
    bool prepared_;
  };

  enum class SIMagnetType
  {
    PlanarUndulator,
    UniformDipole
  };

  /* Static laboratory-frame magnetic elements.  A planar undulator includes
   * divergence-compatible Gaussian end fields; a dipole is a compact hard-
   * edge element. */
  class SIMagneticElement
  {
  public:
    SIMagneticElement();

    SIMagnetType type;
    FieldVector<Double> center; /* m; center[2] is the entrance */
    Double length;              /* m */
    Double period;              /* m; planar undulator only */
    Double peakMagneticField;   /* T */
    Double polarizationAngle;   /* rad from +x */
    bool gaussianFringe;
    Double fringeRelativeCutoff; /* raw Gaussian at compact-support edge */

    void prepare();
    void fieldsLab(const FieldVector<Double>& positionLab,
                   Double timeLab, SIFieldValue& fields) const;

    Double physicalEntranceLab() const;
    Double physicalExitLab() const;
    Double interactionEntranceLab() const;
    Double interactionExitLab() const;

    static SIMagneticElement planarUndulatorFromK(
        Double strengthParameter, Double periodSI,
        Double entranceZSI, unsigned int numberOfPeriods,
        Double polarizationAngleRadians);

  private:
    Double waveNumber_;
    Double fringeExtent_;
    Double cosineAngle_;
    Double sineAngle_;
    bool prepared_;
  };

  /* This classification describes evolution, not how a field happens to be
   * evaluated.  MaxwellIncident is injected into and advanced by the Yee
   * solver.  PrescribedLab is an externally maintained laboratory-frame field
   * which is transformed and added only when pushing particles. */
  enum class SIFieldEvolution
  {
    MaxwellIncident,
    PrescribedLab
  };

  /* Source roles remain explicit to prevent double counting. Incident waves
   * live on the Maxwell grid. Laboratory magnets use the one PrescribedLab
   * path; evaluation method does not create another coupling category. */
  class SIFieldSourceSet
  {
  public:
    void addWave(const SIWaveSource& source, SIFieldEvolution evolution);
    void addPrescribedLabMagnet(const SIMagneticElement& source);

    std::size_t maxwellIncidentWaveCount() const;
    std::size_t prescribedLabWaveCount() const;
    std::size_t prescribedLabMagnetCount() const;

    void sampleMaxwellIncidentLab(
        const FieldVector<Double>& positionLab,
        Double timeLab, SIFieldValue& fields) const;
    void samplePrescribedLab(const FieldVector<Double>& positionLab,
                             Double timeLab, SIFieldValue& fields) const;

    void sampleMaxwellIncidentBox(
        const FieldVector<Double>& positionBox,
        Double timeBox, const BoostFrameTransform& frame,
        SIFieldValue& fields) const;
    void samplePrescribedBox(const FieldVector<Double>& positionBox,
                             Double timeBox,
                             const BoostFrameTransform& frame,
                             SIFieldValue& fields) const;
    void addPrescribedBox(const FieldVector<Double>& positionBox,
                          Double timeBox,
                          const BoostFrameTransform& frame,
                          FieldVector<Double>& electricBox,
                          FieldVector<Double>& magneticBox) const;

  private:
    void transformSampleToBox(const FieldVector<Double>& positionBox,
                              Double timeBox,
                              const BoostFrameTransform& frame,
                              bool incident,
                              SIFieldValue& fields) const;

    std::vector<SIWaveSource> maxwellIncidentWaves_;
    std::vector<SIWaveSource> prescribedLabWaves_;
    std::vector<SIMagneticElement> prescribedLabMagnets_;
  };

}

#endif
