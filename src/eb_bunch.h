#ifndef DIRECT_EB_BUNCH_H
#define DIRECT_EB_BUNCH_H

#include <vector>

#include <mpi.h>

#include "boostframe.h"
#include "eb_particles.h"

namespace fel
{
  struct SIBunchPlacement
  {
    bool hasFirstInteractionEntrance;
    Double firstInteractionEntranceLab; /* m */
    Double referencePositionLab;        /* m */
    Double recommendationMarginLab;    /* m */

    SIBunchPlacement();
  };

  struct SIBunchPlacementReport
  {
    Double relativeTailLab;
    Double relativeHeadLab;
    Double referencePositionLab;
    Double headAfterLab;
    Double actualHeadDistance;
    Double recommendedMaximumReferencePosition;
    unsigned long long particles;

    SIBunchPlacementReport();
  };

  struct SILabPlaneProjectionReport
  {
    Double inputPlaneLab;
    Double referencePositionLab;
    Double minimumForwardDistance;
    Double maximumForwardDistance;
    Double recommendedMinimumReferencePosition;
    Double referenceTimeOffsetLab;
    Double meanLongitudinalBeta;
    unsigned long long particles;

    SILabPlaneProjectionReport();
  };

  struct SIBunchBoostReport
  {
    Double earliestLabEventTime;
    Double latestLabEventTime;
    Double maximumAbsoluteDriftTime;
    Double minimumLabGamma;
    Double maximumLabGamma;
    Double maximumRelativeMomentumRoundTripError;
    Double maximumRelativeGammaRoundTripError;

    SIBunchBoostReport();
  };

  class SIBunchPreprocessor
  {
  public:
    /* Version-3 HDF5 records give x/y at a fixed laboratory observation
     * plane and a signed longitudinal offset from the bunch reference.  Form
     * the common-lab-time snapshot by advancing every record to
     * referencePositionLab + relative_z along its supplied straight-line
     * velocity.  Backward projection is rejected. */
    static SILabPlaneProjectionReport projectLabPlaneToSnapshot(
        std::vector<RelativisticParticleSI>& particles,
        Double inputPlaneLab,
        Double referencePositionLab,
        MPI_Comm communicator);

    /* Version-4 records preserve Elegant's fixed-plane arrival time instead
     * of first approximating it by a longitudinal offset.  Synchronize the
     * crossing events to one common lab time with their individual supplied
     * velocities.  The weighted longitudinal centroid is placed at
     * referencePositionLab. */
    static SILabPlaneProjectionReport projectLabPlaneEventsToSnapshot(
        std::vector<RelativisticParticleSI>& particles,
        Double inputPlaneLab,
        Double referencePositionLab,
        Double longitudinalOffsetLab,
        MPI_Comm communicator);

    /* Input z coordinates are relative to a user-defined laboratory bunch
     * centre. The optional first-interaction coordinate is a magnetic-device
     * placement constraint; detector-only and element-free propagation do not
     * invent a magnetic entrance. */
    static SIBunchPlacementReport placeRelativeLabSnapshot(
        std::vector<RelativisticParticleSI>& particles,
        const SIBunchPlacement& placement,
        MPI_Comm communicator);

    /* Map every particle from the common lab snapshot to its event on the
     * common t_box=0 simultaneity plane, then transform coordinates and
     * proper velocity.  The caller selects the anchor event through frame's
     * lab origins; production initialization anchors the downstream bunch
     * front so the large relativity-of-simultaneity span remains a diagnostic
     * rather than a fictitious required drift length. */
    static SIBunchBoostReport boostLabSnapshotToBoxTimeZero(
        std::vector<RelativisticParticleSI>& particles,
        const BoostFrameTransform& frame,
        Double snapshotTimeLab = 0.0);
  };
}

#endif
