#include "eb_bunch.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fel
{
  SIBunchPlacement::SIBunchPlacement()
    : firstInteractionEntranceLab(std::numeric_limits<Double>::quiet_NaN()),
      referencePositionLab(std::numeric_limits<Double>::quiet_NaN()),
      recommendationMarginLab(0.0)
  {}

  SIBunchPlacementReport::SIBunchPlacementReport()
    : relativeHeadLab(0.0), referencePositionLab(0.0), headAfterLab(0.0),
      actualHeadDistance(0.0), recommendedMaximumReferencePosition(0.0),
      particles(0)
  {}

  SIBunchBoostReport::SIBunchBoostReport()
    : earliestLabEventTime(0.0), latestLabEventTime(0.0),
      maximumAbsoluteDriftTime(0.0)
  {}

  SIBunchPlacementReport SIBunchPreprocessor::placeRelativeLabSnapshot(
      std::vector<RelativisticParticleSI>& particles,
      const SIBunchPlacement& placement, MPI_Comm communicator)
  {
    if (communicator == MPI_COMM_NULL)
      throw std::invalid_argument("Bunch placement communicator is null");
    if (!std::isfinite(placement.firstInteractionEntranceLab) ||
        !std::isfinite(placement.referencePositionLab) ||
        !std::isfinite(placement.recommendationMarginLab) ||
        placement.recommendationMarginLab < 0.0)
      throw std::invalid_argument(
        "Relative bunch placement requires finite lab coordinates and a nonnegative margin");

    const unsigned long long localCount =
      static_cast<unsigned long long>(particles.size());
    unsigned long long globalCount = 0;
    MPI_Allreduce(&localCount, &globalCount, 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, communicator);
    if (globalCount == 0)
      throw std::invalid_argument("Cannot place an empty particle bunch");

    Double localHead = -std::numeric_limits<Double>::infinity();
    for (std::size_t i = 0; i < particles.size(); ++i)
      {
        if (!std::isfinite(particles[i].position[2]))
          throw std::invalid_argument("Bunch position must be finite");
        localHead = std::max(localHead, particles[i].position[2]);
      }
    Double globalRelativeHead = 0.0;
    MPI_Allreduce(&localHead, &globalRelativeHead, 1, MPI_DOUBLE, MPI_MAX,
                  communicator);

    SIBunchPlacementReport report;
    report.particles = globalCount;
    report.relativeHeadLab = globalRelativeHead;
    report.referencePositionLab = placement.referencePositionLab;
    report.headAfterLab = report.referencePositionLab + globalRelativeHead;
    report.actualHeadDistance = placement.firstInteractionEntranceLab -
                                report.headAfterLab;
    report.recommendedMaximumReferencePosition =
      placement.firstInteractionEntranceLab - globalRelativeHead -
      placement.recommendationMarginLab;

    for (std::size_t i = 0; i < particles.size(); ++i)
      particles[i].position[2] += report.referencePositionLab;
    return report;
  }

  SIBunchBoostReport SIBunchPreprocessor::boostLabSnapshotToBoxTimeZero(
      std::vector<RelativisticParticleSI>& particles,
      const BoostFrameTransform& frame, Double snapshotTimeLab)
  {
    if (!std::isfinite(snapshotTimeLab))
      throw std::invalid_argument("Bunch snapshot time must be finite");

    SIBunchBoostReport report;
    report.earliestLabEventTime = std::numeric_limits<Double>::infinity();
    report.latestLabEventTime = -std::numeric_limits<Double>::infinity();

    for (std::size_t index = 0; index < particles.size(); ++index)
      {
        RelativisticParticleSI& particle = particles[index];
        const Double gammaParticle =
          BoostFrameTransform::gammaFromProperVelocity(
            particle.properVelocity);
        FieldVector<Double> betaParticle(0.0);
        betaParticle.mv(1.0 / gammaParticle, particle.properVelocity);

        const Double denominator =
          1.0 - frame.beta() * betaParticle[2];
        if (!(denominator > 0.0) || !std::isfinite(denominator))
          throw std::runtime_error(
            "Particle does not intersect the requested boosted simultaneity plane");

        const Double eventTimeLab =
          (frame.labTimeOrigin() +
           frame.beta() *
             (particle.position[2] - frame.labZOrigin()) / SI::c -
           frame.beta() * betaParticle[2] * snapshotTimeLab) /
          denominator;
        const Double driftTime = eventTimeLab - snapshotTimeLab;
        FieldVector<Double> eventPosition(particle.position);
        eventPosition.pmv(SI::c * driftTime, betaParticle);

        Double eventTimeBox = 0.0;
        Double eventZBox = 0.0;
        frame.labToBox(eventTimeLab, eventPosition[2],
                       eventTimeBox, eventZBox);
        const Double timeTolerance = 128.0 *
          std::numeric_limits<Double>::epsilon() *
          std::max(1.0, std::abs(eventTimeLab));
        if (std::abs(eventTimeBox) > timeTolerance)
          throw std::runtime_error(
            "Free-drift bunch transform did not reach box time zero");

        particle.position[0] = eventPosition[0];
        particle.position[1] = eventPosition[1];
        particle.position[2] = eventZBox;
        particle.properVelocity[2] = frame.gamma() *
          (particle.properVelocity[2] -
           frame.beta() * gammaParticle);

        report.earliestLabEventTime = std::min(
          report.earliestLabEventTime, eventTimeLab);
        report.latestLabEventTime = std::max(
          report.latestLabEventTime, eventTimeLab);
        report.maximumAbsoluteDriftTime = std::max(
          report.maximumAbsoluteDriftTime, std::abs(driftTime));
      }

    if (particles.empty())
      report.earliestLabEventTime = report.latestLabEventTime = 0.0;
    return report;
  }
}
