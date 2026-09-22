#include "eb_bunch.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fel
{
  SIBunchPlacement::SIBunchPlacement()
    : mode(SIBunchPlacementMode::AbsoluteLab),
      firstElementEntranceLab(std::numeric_limits<Double>::quiet_NaN()),
      headDistanceLab(0.0)
  {}

  SIBunchPlacementReport::SIBunchPlacementReport()
    : headBeforeLab(0.0), headAfterLab(0.0),
      longitudinalTranslation(0.0), actualHeadDistance(0.0),
      particles(0)
  {}

  SIBunchBoostReport::SIBunchBoostReport()
    : earliestLabEventTime(0.0), latestLabEventTime(0.0),
      maximumAbsoluteDriftTime(0.0)
  {}

  SIBunchPlacementReport SIBunchPreprocessor::placeLabSnapshot(
      std::vector<RelativisticParticleSI>& particles,
      const SIBunchPlacement& placement, MPI_Comm communicator)
  {
    if (communicator == MPI_COMM_NULL)
      throw std::invalid_argument("Bunch placement communicator is null");
    if (placement.mode == SIBunchPlacementMode::HeadToFirstElement &&
        (!std::isfinite(placement.firstElementEntranceLab) ||
         !std::isfinite(placement.headDistanceLab) ||
         placement.headDistanceLab < 0.0))
      throw std::invalid_argument(
        "Head-to-element placement requires a finite entrance and nonnegative distance");

    const unsigned long long localCount =
      static_cast<unsigned long long>(particles.size());
    unsigned long long globalCount = 0;
    MPI_Allreduce(&localCount, &globalCount, 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, communicator);
    if (globalCount == 0)
      throw std::invalid_argument("Cannot place an empty particle bunch");

    Double localHead = -std::numeric_limits<Double>::infinity();
    for (std::size_t particle = 0; particle < particles.size(); ++particle)
      {
        const Double z = particles[particle].position[2];
        if (!std::isfinite(z))
          throw std::invalid_argument("Bunch position must be finite");
        localHead = std::max(localHead, z);
      }
    Double globalHead = 0.0;
    MPI_Allreduce(&localHead, &globalHead, 1, MPI_DOUBLE, MPI_MAX,
                  communicator);

    SIBunchPlacementReport report;
    report.headBeforeLab = globalHead;
    report.particles = globalCount;
    if (placement.mode == SIBunchPlacementMode::HeadToFirstElement)
      report.longitudinalTranslation =
        placement.firstElementEntranceLab - placement.headDistanceLab -
        globalHead;

    for (std::size_t particle = 0; particle < particles.size(); ++particle)
      particles[particle].position[2] += report.longitudinalTranslation;

    report.headAfterLab = globalHead + report.longitudinalTranslation;
    report.actualHeadDistance =
      std::isfinite(placement.firstElementEntranceLab) ?
      placement.firstElementEntranceLab - report.headAfterLab :
      std::numeric_limits<Double>::quiet_NaN();
    return report;
  }

  SIBunchBoostReport SIBunchPreprocessor::boostLabSnapshotToBoxTimeZero(
      std::vector<RelativisticParticleSI>& particles,
      const BoostFrameTransform& frame, Double snapshotTimeLab)
  {
    if (!std::isfinite(snapshotTimeLab))
      throw std::invalid_argument("Bunch snapshot time must be finite");

    SIBunchBoostReport report;
    report.earliestLabEventTime =
      std::numeric_limits<Double>::infinity();
    report.latestLabEventTime =
      -std::numeric_limits<Double>::infinity();

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
