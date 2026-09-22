#ifndef DIRECT_EB_BUNCH_H
#define DIRECT_EB_BUNCH_H

#include <vector>

#include <mpi.h>

#include "boostframe.h"
#include "eb_particles.h"

namespace fel
{
  enum class SIBunchPlacementMode
  {
    AbsoluteLab,
    HeadToFirstElement
  };

  struct SIBunchPlacement
  {
    SIBunchPlacementMode mode;
    Double firstElementEntranceLab; /* m */
    Double headDistanceLab;         /* m */

    SIBunchPlacement();
  };

  struct SIBunchPlacementReport
  {
    Double headBeforeLab;
    Double headAfterLab;
    Double longitudinalTranslation;
    Double actualHeadDistance;
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

  /* Converts an explicitly defined laboratory snapshot into the common
   * t_box=0 hypersurface. Particles are freely drifted from the lab snapshot
   * to their individual Lorentz-simultaneous events before coordinates and
   * proper velocities are transformed. */
  class SIBunchPreprocessor
  {
  public:
    static SIBunchPlacementReport placeLabSnapshot(
        std::vector<RelativisticParticleSI>& particles,
        const SIBunchPlacement& placement,
        MPI_Comm communicator);

    static SIBunchBoostReport boostLabSnapshotToBoxTimeZero(
        std::vector<RelativisticParticleSI>& particles,
        const BoostFrameTransform& frame,
        Double snapshotTimeLab = 0.0);
  };
}

#endif
