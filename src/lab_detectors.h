#ifndef DIRECT_EB_LAB_DETECTORS_H
#define DIRECT_EB_LAB_DETECTORS_H

#include <memory>

#include <mpi.h>

#include "boostframe.h"
#include "config.h"
#include "eb_field.h"
#include "eb_particles.h"

namespace fel
{
  /* The manager is constructed only when at least one laboratory detector is
   * configured.  All HDF5 handles live on rank zero; other ranks contribute
   * transient samples through MPI and never open detector files. */
  class LabDetectorManager
  {
  public:
    LabDetectorManager(const DetectorConfig& config,
                       const EBGridGeometry& globalGeometry,
                       const FieldVector<Double>& globalOriginBox,
                       const FieldVector<Double>& localOriginBox,
                       const BoostFrameTransform& frame,
                       MPI_Comm communicator);
    ~LabDetectorManager();

    LabDetectorManager(const LabDetectorManager&) = delete;
    LabDetectorManager& operator=(const LabDetectorManager&) = delete;

    /* Record crossings of any configured fixed lab-z particle plane during
     * one complete particle push.  No field or trajectory state is changed. */
    void captureParticleStep(const RelativisticParticleSI& before,
                             const RelativisticParticleSI& after,
                             Double timeBoxBefore,
                             Double timeBoxAfter);

    /* One detector-only collective is entered per particle step when particle
     * planes exist.  It returns immediately after a zero-count Allreduce on
     * the overwhelmingly common no-crossing steps. */
    void collectParticleCrossings();

    /* Sample due fixed-z field planes.  Only the rank owning the moving
     * boosted-frame z slice interpolates the x-y plane; rank zero is the sole
     * HDF5 writer. */
    void sampleFieldPlanes(const EBFieldGrid& localFields,
                           Double timeBox);

    void close(bool completed);

  private:
    class Impl;
    std::unique_ptr<Impl> impl_;
  };
}

#endif
