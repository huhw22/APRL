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
  /* Affine endpoint four-momentum quadrature at a crossing. This is not an
   * on-shell particle snapshot; mass/charge already include macro weight. */
  struct EnergyPlaneCrossing
  {
    std::size_t planeIndex;
    std::int32_t direction;
    Double timeBox, timeLab, x, y;
    std::uint64_t particleId;
    long double charge;
    long double kineticBox, totalBox, cpzBox;
    long double kineticLab, totalLab, cpzLab;
  };
  /* Diagnostic-only affine four-momentum of an actual converged particle
   * substep. Interpolating its endpoint worldline is not an analytic crossing
   * solve, and the affine (gamma,u) pair need not be exactly on shell. */
  struct SubstepParticleCrossing
  {
    std::size_t planeIndex;
    std::uint64_t particleId, sourceId;
    Double timeBox, timeLab, position[3];
    Double properVelocityLab[3], properVelocityBox[3];
    Double gammaLab, gammaBox, kineticLab, kineticBox;
    Double charge, mass, weight;
    std::int32_t direction;
  };
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

    bool hasEnergyPlanes() const;
    bool hasSubstepParticlePlanes() const;
    void previewSubstepParticleStep(const RelativisticParticleSI& before,
                                   const RelativisticParticleSI& after,
                                   Double timeBoxBefore, Double timeBoxAfter,
                                   std::vector<SubstepParticleCrossing>& events) const;
    void commitSubstepParticleCrossings(
        const std::vector<SubstepParticleCrossing>& events, Double maximumTimeBox);
    void captureSubstepParticleStep(const RelativisticParticleSI& before,
                                   const RelativisticParticleSI& after,
                                   Double timeBoxBefore, Double timeBoxAfter);
    void previewEnergyParticleStep(const RelativisticParticleSI& before,
                                   const RelativisticParticleSI& after,
                                   Double timeBoxBefore, Double timeBoxAfter,
                                   std::vector<EnergyPlaneCrossing>& events) const;
    void commitEnergyParticleCrossings(
        const std::vector<EnergyPlaneCrossing>& events, Double maximumTimeBox);
    void captureEnergyParticleStep(const RelativisticParticleSI& before,
                                   const RelativisticParticleSI& after,
                                   Double timeBoxBefore, Double timeBoxAfter);

    /* One detector-only collective is entered per particle step when a
     * particle plane or ballistic-reference plane exists.  It returns
     * immediately after a zero-count Allreduce on the overwhelmingly common
     * no-crossing steps. */
    void collectParticleCrossings();

    /* Update a longitudinally translated grid after a complete-step remap,
     * before sampleFieldPlanes. Geometry and the MPI z-slab layout must stay
     * fixed; x/y origins cannot change because the laboratory pixel aperture
     * is fixed when the writers open. This performs no communication and
     * does not rebase the Lorentz frame, detector cadence or stored events. */
    void updateGridOrigins(const FieldVector<Double>& globalOriginBox,
                           const FieldVector<Double>& localOriginBox);

    /* Sample due fixed-z field planes with the global Yee interpolation
     * stencil. The owning rank obtains one adjacent Ez/Bx/By layer only when
     * the stencil crosses an internal z-slab boundary; rank zero is the sole
     * HDF5 writer. Stored E/B leapfrog times are unchanged. */
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
