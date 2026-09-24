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
    Double firstInteractionEntranceLab; /* m */
    Double referencePositionLab;        /* m */
    Double recommendationMarginLab;    /* m */

    SIBunchPlacement();
  };

  struct SIBunchPlacementReport
  {
    Double relativeHeadLab;
    Double referencePositionLab;
    Double headAfterLab;
    Double actualHeadDistance;
    Double recommendedMaximumReferencePosition;
    unsigned long long particles;

    SIBunchPlacementReport();
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
    /* Input z coordinates are relative to a user-defined laboratory bunch
     * centre. The physical beamline origin remains independent of the finite
     * interaction entrance created by an element fringe. */
    static SIBunchPlacementReport placeRelativeLabSnapshot(
        std::vector<RelativisticParticleSI>& particles,
        const SIBunchPlacement& placement,
        MPI_Comm communicator);

    /* Free-drift every particle from the common lab snapshot to its event on
     * the common t_box=0 simultaneity plane, then transform coordinates and
     * proper velocity. */
    static SIBunchBoostReport boostLabSnapshotToBoxTimeZero(
        std::vector<RelativisticParticleSI>& particles,
        const BoostFrameTransform& frame,
        Double snapshotTimeLab = 0.0);
  };
}

#endif
