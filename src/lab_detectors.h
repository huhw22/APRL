#ifndef DIRECT_EB_LAB_DETECTORS_H
#define DIRECT_EB_LAB_DETECTORS_H

#include <memory>
#include <cstdint>
#include <string>
#include <vector>

#include <mpi.h>

#include "boostframe.h"
#include "config.h"
#include "eb_field.h"
#include "eb_particles.h"
#include "run_metadata.h"

namespace aprl
{
  /* An immutable prefix of an existing detector file at a checkpoint.  A
   * restarted run writes a new segment; reference entries in these prefixes
   * must remain available to radiation/near-field postprocessing. */
  struct DetectorCheckpointPrefix
  {
    std::string planeName;
    std::string filename;
    std::string group;
    std::uint64_t committedRecords;
    std::uint64_t validationCommittedRecords;
  };

  /* Detector-only subtraction hook for a known laboratory background such as
   * a future analytical injected laser.  It is deliberately separate from
   * EBFieldGrid: implementations are sampled only while writing a detector
   * plane and can never modify the propagated E/B state. */
  class LabFieldDetectorBackground
  {
  public:
    virtual ~LabFieldDetectorBackground() {}
    virtual const char* name() const = 0;
    virtual void sampleLab(const FieldVector<Double>& positionLab,
                           Double timeLab,
                           FieldVector<Double>& electricLab,
                           FieldVector<Double>& magneticLab) const = 0;
  };

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
                       MPI_Comm communicator,
                       const LabFieldDetectorBackground* externalBackground = NULL,
                       bool overwrite = false,
                       const RunMetadata* metadata = NULL);
    ~LabDetectorManager();

    LabDetectorManager(const LabDetectorManager&) = delete;
    LabDetectorManager& operator=(const LabDetectorManager&) = delete;

    /* Record crossings of any configured fixed lab-z particle plane during
     * one complete particle push.  No field or trajectory state is changed. */
    void captureParticleStep(const RelativisticParticleSI& before,
                             const RelativisticParticleSI& after,
                             Double timeBoxBefore,
                             Double timeBoxAfter);

    /* One detector-only collective is entered per particle step when a
     * particle plane or ballistic-reference plane exists.  It returns
     * immediately after a zero-count Allreduce on the overwhelmingly common
     * no-crossing steps. */
    void collectParticleCrossings();

    /* Sample due fixed-z field planes.  Only the rank owning the moving
     * boosted-frame z slice interpolates the x-y plane; rank zero is the sole
     * HDF5 writer. */
    void sampleFieldPlanes(const EBFieldGrid& localFields,
                           Double timeBox);

    /* Checkpoints are taken after collectParticleCrossings and the complete
     * field-step diagnostics. These methods perform no MPI operations and
     * are never called on the ordinary, checkpoint-disabled stepping path.
     * Flush publishes the old output prefix but leaves all writers open.
     * State is rank-local: the root record additionally owns pending
     * ballistic-validation entries and cumulative validation statistics. */
    void flushForCheckpoint();
    std::vector<unsigned char> checkpointState() const;
    void restoreCheckpointState(const std::vector<unsigned char>& state);
    /* Root-only prefixes, available after flushForCheckpoint; empty on
     * other ranks. Paths are the configured detector output paths. */
    std::vector<DetectorCheckpointPrefix> checkpointPrefixes() const;

    void close(bool completed);

  private:
    class Impl;
    std::unique_ptr<Impl> impl_;
  };
}

#endif
