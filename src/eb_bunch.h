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
    Double firstElementEntranceLab;    /* m */
    Double referenceDistanceLab;       /* m */
    Double recommendationMarginLab;    /* m */

    SIBunchPlacement();
  };

  struct SIBunchPlacementReport
  {
    Double relativeHeadLab;
    Double referencePositionLab;
    Double headAfterLab;
    Double actualHeadDistance;
    Double recommendedReferenceDistance;
    unsigned long long particles;

    SIBunchPlacementReport();
  };

  struct SIBunchBoostReport
  {
    Double earliestLabEventTime;
    Double latestLabEventTime;
    Double maximumAbsoluteDriftTime;

    SIBunchBoostReport();
  };

  class SIBunchPreprocessor
  {
  public:
    /* Input z coordinates are relative to a user-defined laboratory reference
     * point. referenceDistanceLab places that point upstream of the first
     * magnetic-element entrance. */
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
