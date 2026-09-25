#include "eb_bunch.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace aprl
{
  SIBunchPlacement::SIBunchPlacement()
    : hasFirstInteractionEntrance(false),
      firstInteractionEntranceLab(std::numeric_limits<Double>::quiet_NaN()),
      referencePositionLab(std::numeric_limits<Double>::quiet_NaN()),
      recommendationMarginLab(0.0)
  {}

  SIBunchPlacementReport::SIBunchPlacementReport()
    : relativeTailLab(0.0), relativeHeadLab(0.0),
      referencePositionLab(0.0), headAfterLab(0.0),
      actualHeadDistance(0.0), recommendedMaximumReferencePosition(0.0),
      particles(0)
  {}

  SILabPlaneProjectionReport::SILabPlaneProjectionReport()
    : inputPlaneLab(0.0), referencePositionLab(0.0),
      minimumForwardDistance(0.0), maximumForwardDistance(0.0),
      recommendedMinimumReferencePosition(0.0),
      referenceTimeOffsetLab(0.0), meanLongitudinalBeta(0.0), particles(0)
  {}

  SIBunchBoostReport::SIBunchBoostReport()
    : earliestLabEventTime(0.0), latestLabEventTime(0.0),
      maximumAbsoluteDriftTime(0.0),
      minimumLabGamma(std::numeric_limits<Double>::infinity()),
      maximumLabGamma(1.0),
      maximumRelativeMomentumRoundTripError(0.0),
      maximumRelativeGammaRoundTripError(0.0)
  {}

  SILabPlaneProjectionReport
  SIBunchPreprocessor::projectLabPlaneToSnapshot(
      std::vector<RelativisticParticleSI>& particles,
      Double inputPlaneLab, Double referencePositionLab,
      MPI_Comm communicator)
  {
    if (communicator == MPI_COMM_NULL)
      throw std::invalid_argument("Lab-plane projection communicator is null");
    if (!std::isfinite(inputPlaneLab) ||
        !std::isfinite(referencePositionLab))
      throw std::invalid_argument(
        "Lab-plane projection coordinates must be finite");

    const unsigned long long localCount =
      static_cast<unsigned long long>(particles.size());
    unsigned long long globalCount = 0;
    MPI_Allreduce(&localCount, &globalCount, 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, communicator);
    if (globalCount == 0)
      throw std::invalid_argument("Cannot project an empty particle bunch");

    Double localMinimum = std::numeric_limits<Double>::infinity();
    Double localMaximum = -std::numeric_limits<Double>::infinity();
    int localInvalidDirection = 0;
    for (std::size_t index = 0; index < particles.size(); ++index)
      {
        const RelativisticParticleSI& particle = particles[index];
        if (!(particle.properVelocity[2] > 0.0) ||
            !std::isfinite(particle.properVelocity[2]))
          localInvalidDirection = 1;
        const Double distance = referencePositionLab + particle.position[2] -
                                inputPlaneLab;
        localMinimum = std::min(localMinimum, distance);
        localMaximum = std::max(localMaximum, distance);
      }

    int globalInvalidDirection = 0;
    MPI_Allreduce(&localInvalidDirection, &globalInvalidDirection, 1,
                  MPI_INT, MPI_MAX, communicator);
    Double globalMinimum = 0.0;
    Double globalMaximum = 0.0;
    MPI_Allreduce(&localMinimum, &globalMinimum, 1, MPI_DOUBLE, MPI_MIN,
                  communicator);
    MPI_Allreduce(&localMaximum, &globalMaximum, 1, MPI_DOUBLE, MPI_MAX,
                  communicator);
    if (globalInvalidDirection)
      throw std::invalid_argument(
        "Lab-plane particle input requires every particle to have positive longitudinal proper velocity");

    const Double scale = std::max(1.0,
      std::max(std::abs(inputPlaneLab),
               std::max(std::abs(referencePositionLab),
                        std::abs(globalMinimum))));
    const Double tolerance = 128.0 *
      std::numeric_limits<Double>::epsilon() * scale;
    if (globalMinimum < -tolerance)
      {
        const Double recommendedCenter =
          referencePositionLab - globalMinimum;
        std::ostringstream message;
        message << std::setprecision(16)
          << "The reconstructed laboratory snapshot would back-propagate at "
             "least one Elegant plane record. Set "
             "beam.reference.initial_center_z to at least "
          << recommendedCenter
          << " m, move input_plane_z upstream, or correct the sign of the "
             "longitudinal offsets.";
        throw std::invalid_argument(message.str());
      }

    for (std::size_t index = 0; index < particles.size(); ++index)
      {
        RelativisticParticleSI& particle = particles[index];
        const Double distance = referencePositionLab + particle.position[2] -
                                inputPlaneLab;
        particle.position[0] +=
          particle.properVelocity[0] / particle.properVelocity[2] * distance;
        particle.position[1] +=
          particle.properVelocity[1] / particle.properVelocity[2] * distance;
      }

    SILabPlaneProjectionReport report;
    report.inputPlaneLab = inputPlaneLab;
    report.referencePositionLab = referencePositionLab;
    report.minimumForwardDistance = std::max(0.0, globalMinimum);
    report.maximumForwardDistance = globalMaximum;
    report.recommendedMinimumReferencePosition =
      referencePositionLab - globalMinimum;
    report.particles = globalCount;
    return report;
  }

  SILabPlaneProjectionReport
  SIBunchPreprocessor::projectLabPlaneEventsToSnapshot(
      std::vector<RelativisticParticleSI>& particles,
      Double inputPlaneLab, Double referencePositionLab,
      Double longitudinalOffsetLab, MPI_Comm communicator)
  {
    if (communicator == MPI_COMM_NULL)
      throw std::invalid_argument(
        "Lab-plane event projection communicator is null");
    if (!std::isfinite(inputPlaneLab) ||
        !std::isfinite(referencePositionLab) ||
        !std::isfinite(longitudinalOffsetLab))
      throw std::invalid_argument(
        "Lab-plane event projection coordinates must be finite");

    const unsigned long long localCount =
      static_cast<unsigned long long>(particles.size());
    unsigned long long globalCount = 0;
    MPI_Allreduce(&localCount, &globalCount, 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, communicator);
    if (globalCount == 0)
      throw std::invalid_argument(
        "Cannot project an empty particle event set");

    long double localWeight = 0.0L;
    long double localWeightedBetaZ = 0.0L;
    long double localWeightedBetaZTime = 0.0L;
    Double localLatestTime = -std::numeric_limits<Double>::infinity();
    int localInvalid = 0;
    for (std::size_t index = 0; index < particles.size(); ++index)
      {
        const RelativisticParticleSI& particle = particles[index];
        const Double timeOffset = particle.position[2];
        const Double gamma = std::sqrt(1.0 +
          particle.properVelocity.norm2());
        const Double betaZ = particle.properVelocity[2] / gamma;
        if (!(particle.weight > 0.0) || !std::isfinite(particle.weight) ||
            !(betaZ > 0.0) || !std::isfinite(betaZ) ||
            !std::isfinite(timeOffset))
          localInvalid = 1;
        localWeight += static_cast<long double>(particle.weight);
        localWeightedBetaZ += static_cast<long double>(particle.weight) *
          static_cast<long double>(betaZ);
        localWeightedBetaZTime +=
          static_cast<long double>(particle.weight) *
          static_cast<long double>(betaZ) *
          static_cast<long double>(timeOffset);
        localLatestTime = std::max(localLatestTime, timeOffset);
      }

    int globalInvalid = 0;
    MPI_Allreduce(&localInvalid, &globalInvalid, 1, MPI_INT, MPI_MAX,
                  communicator);
    if (globalInvalid)
      throw std::invalid_argument(
        "Lab-plane event input requires finite positive weights, finite arrival times, and positive longitudinal velocity");

    long double globalWeight = 0.0L;
    long double globalWeightedBetaZ = 0.0L;
    long double globalWeightedBetaZTime = 0.0L;
    MPI_Allreduce(&localWeight, &globalWeight, 1, MPI_LONG_DOUBLE,
                  MPI_SUM, communicator);
    MPI_Allreduce(&localWeightedBetaZ, &globalWeightedBetaZ, 1,
                  MPI_LONG_DOUBLE, MPI_SUM, communicator);
    MPI_Allreduce(&localWeightedBetaZTime, &globalWeightedBetaZTime, 1,
                  MPI_LONG_DOUBLE, MPI_SUM, communicator);
    Double globalLatestTime = 0.0;
    MPI_Allreduce(&localLatestTime, &globalLatestTime, 1, MPI_DOUBLE,
                  MPI_MAX, communicator);
    if (!(globalWeight > 0.0L) || !(globalWeightedBetaZ > 0.0L))
      throw std::invalid_argument(
        "Lab-plane event input has an invalid global weight or longitudinal velocity sum");

    const long double meanBetaZ =
      globalWeightedBetaZ / globalWeight;
    /* This choice makes the macro-weighted longitudinal centroid land
     * exactly at referencePositionLab before the optional z offset, including
     * correlations between arrival time and longitudinal velocity. */
    const long double referenceTime =
      (static_cast<long double>(referencePositionLab - inputPlaneLab) /
       static_cast<long double>(SI::c) +
       globalWeightedBetaZTime / globalWeight) / meanBetaZ;

    const Double scale = std::max(1.0,
      std::max(std::abs(inputPlaneLab), std::abs(referencePositionLab)));
    const Double distanceTolerance = 128.0 *
      std::numeric_limits<Double>::epsilon() * scale;
    if (referenceTime < static_cast<long double>(globalLatestTime) -
        static_cast<long double>(distanceTolerance / SI::c))
      {
        const long double minimumCenter =
          static_cast<long double>(inputPlaneLab) +
          static_cast<long double>(SI::c) *
          (meanBetaZ * static_cast<long double>(globalLatestTime) -
           globalWeightedBetaZTime / globalWeight);
        std::ostringstream message;
        message << std::setprecision(16)
          << "The reconstructed laboratory snapshot would require backward "
             "propagation of at least one Elegant plane event. Set "
             "beam.reference.initial_center_z to at least "
          << static_cast<Double>(minimumCenter)
          << " m, move input_plane_z upstream, or verify the Elegant t "
             "column and selected page.";
        throw std::invalid_argument(message.str());
      }

    Double localMinimum = std::numeric_limits<Double>::infinity();
    Double localMaximum = -std::numeric_limits<Double>::infinity();
    for (std::size_t index = 0; index < particles.size(); ++index)
      {
        RelativisticParticleSI& particle = particles[index];
        const Double timeOffset = particle.position[2];
        const Double gamma = std::sqrt(1.0 +
          particle.properVelocity.norm2());
        const Double betaZ = particle.properVelocity[2] / gamma;
        const Double deltaTime = static_cast<Double>(referenceTime) -
                                 timeOffset;
        Double distance = betaZ * SI::c * deltaTime;
        if (distance < 0.0 && distance >= -distanceTolerance)
          distance = 0.0;
        localMinimum = std::min(localMinimum, distance);
        localMaximum = std::max(localMaximum, distance);
        particle.position[0] +=
          particle.properVelocity[0] / particle.properVelocity[2] * distance;
        particle.position[1] +=
          particle.properVelocity[1] / particle.properVelocity[2] * distance;
        particle.position[2] = inputPlaneLab + distance -
          referencePositionLab + longitudinalOffsetLab;
      }

    Double globalMinimum = 0.0;
    Double globalMaximum = 0.0;
    MPI_Allreduce(&localMinimum, &globalMinimum, 1, MPI_DOUBLE, MPI_MIN,
                  communicator);
    MPI_Allreduce(&localMaximum, &globalMaximum, 1, MPI_DOUBLE, MPI_MAX,
                  communicator);

    SILabPlaneProjectionReport report;
    report.inputPlaneLab = inputPlaneLab;
    report.referencePositionLab = referencePositionLab;
    report.minimumForwardDistance = std::max(0.0, globalMinimum);
    report.maximumForwardDistance = globalMaximum;
    report.recommendedMinimumReferencePosition = static_cast<Double>(
      static_cast<long double>(inputPlaneLab) +
      static_cast<long double>(SI::c) *
      (meanBetaZ * static_cast<long double>(globalLatestTime) -
       globalWeightedBetaZTime / globalWeight));
    report.referenceTimeOffsetLab = static_cast<Double>(referenceTime);
    report.meanLongitudinalBeta = static_cast<Double>(meanBetaZ);
    report.particles = globalCount;
    return report;
  }

  SIBunchPlacementReport SIBunchPreprocessor::placeRelativeLabSnapshot(
      std::vector<RelativisticParticleSI>& particles,
      const SIBunchPlacement& placement, MPI_Comm communicator)
  {
    if (communicator == MPI_COMM_NULL)
      throw std::invalid_argument("Bunch placement communicator is null");
    if ((placement.hasFirstInteractionEntrance &&
         !std::isfinite(placement.firstInteractionEntranceLab)) ||
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

    Double localTail = std::numeric_limits<Double>::infinity();
    Double localHead = -std::numeric_limits<Double>::infinity();
    for (std::size_t i = 0; i < particles.size(); ++i)
      {
        if (!std::isfinite(particles[i].position[2]))
          throw std::invalid_argument("Bunch position must be finite");
        localTail = std::min(localTail, particles[i].position[2]);
        localHead = std::max(localHead, particles[i].position[2]);
      }
    Double globalRelativeTail = 0.0;
    Double globalRelativeHead = 0.0;
    MPI_Allreduce(&localTail, &globalRelativeTail, 1, MPI_DOUBLE, MPI_MIN,
                  communicator);
    MPI_Allreduce(&localHead, &globalRelativeHead, 1, MPI_DOUBLE, MPI_MAX,
                  communicator);

    SIBunchPlacementReport report;
    report.particles = globalCount;
    report.relativeTailLab = globalRelativeTail;
    report.relativeHeadLab = globalRelativeHead;
    report.referencePositionLab = placement.referencePositionLab;
    report.headAfterLab = report.referencePositionLab + globalRelativeHead;
    if (placement.hasFirstInteractionEntrance)
      {
        report.actualHeadDistance = placement.firstInteractionEntranceLab -
                                    report.headAfterLab;
        report.recommendedMaximumReferencePosition =
          placement.firstInteractionEntranceLab - globalRelativeHead -
          placement.recommendationMarginLab;
      }
    else
      {
        report.actualHeadDistance =
          std::numeric_limits<Double>::infinity();
        report.recommendedMaximumReferencePosition =
          std::numeric_limits<Double>::infinity();
      }

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
        const FieldVector<Double> properVelocityLab(
          particle.properVelocity);
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
        frame.properVelocityLabToBox(
          properVelocityLab, particle.properVelocity);

        FieldVector<Double> reconstructed(0.0);
        frame.properVelocityBoxToLab(
          particle.properVelocity, reconstructed);
        FieldVector<Double> momentumDelta(reconstructed);
        momentumDelta -= properVelocityLab;
        const Double momentumError =
          momentumDelta.norm() /
          std::max(1.0, properVelocityLab.norm());
        const Double reconstructedGamma =
          BoostFrameTransform::gammaFromProperVelocity(reconstructed);

        report.earliestLabEventTime = std::min(
          report.earliestLabEventTime, eventTimeLab);
        report.latestLabEventTime = std::max(
          report.latestLabEventTime, eventTimeLab);
        report.maximumAbsoluteDriftTime = std::max(
          report.maximumAbsoluteDriftTime, std::abs(driftTime));
        report.minimumLabGamma = std::min(
          report.minimumLabGamma, gammaParticle);
        report.maximumLabGamma = std::max(
          report.maximumLabGamma, gammaParticle);
        report.maximumRelativeMomentumRoundTripError = std::max(
          report.maximumRelativeMomentumRoundTripError, momentumError);
        report.maximumRelativeGammaRoundTripError = std::max(
          report.maximumRelativeGammaRoundTripError,
          std::abs(reconstructedGamma - gammaParticle) / gammaParticle);
      }

    if (particles.empty())
      {
        report.earliestLabEventTime = report.latestLabEventTime = 0.0;
        report.minimumLabGamma = std::numeric_limits<Double>::infinity();
      }
    return report;
  }
}
