#ifndef DIRECT_EB_CONFIG_H
#define DIRECT_EB_CONFIG_H

#include <cstddef>
#include <string>
#include <vector>

#include "eb_field.h"
#include "eb_cpml.h"
#include "eb_sources.h"
#include "fieldvector.h"

namespace fel
{
  struct UnitSystem
  {
    Double length;
    Double time;

    UnitSystem();
  };

  struct MeshConfig
  {
    std::size_t cells[3];
    FieldVector<Double> cellSize;
    FieldVector<Double> extent; /* derived once as cells * cellSize */
    FieldVector<Double> center;
    Double duration;
    Double boostGamma;
    unsigned int particleStepsPerUndulatorPeriod;
    unsigned int maximumParticleSubsteps;
    EBMaxwellSolver fieldSolver;

    MeshConfig();
  };

  enum class EBBoundaryType
  {
    Pec,
    Cpml
  };

  struct BoundaryConfig
  {
    EBBoundaryType type;
    EBCPMLParameters cpml;

    BoundaryConfig();
  };

  struct InitialSelfFieldConfig
  {
    bool enabled;
    Double relativeTolerance;
    std::size_t maximumIterations;

    InitialSelfFieldConfig();
  };

  enum class BeamInputType
  {
    Hdf5,
    GeneratedGaussian
  };

  struct BeamInputConfig
  {
    BeamInputType type;
    std::string file;
    Double electrons;
    FieldVector<Double> positionOffset;

    /* generated-gaussian is intentionally a small, deterministic test path. */
    std::size_t macroparticles;
    Double gamma;
    FieldVector<Double> direction;
    FieldVector<Double> center;
    FieldVector<Double> sigmaPosition;
    FieldVector<Double> sigmaProperVelocity;
    unsigned int randomSeed;

    BeamInputConfig();
  };

  struct BeamReferenceConfig
  {
    /* Lab coordinate of the particle file's relative z=0 reference at the
     * input snapshot. The first physical beamline entrance defines z=0. */
    Double initialCenterZ;

    BeamReferenceConfig();
  };

  enum class BeamlineElementRole
  {
    MagneticDevice,
    FieldDetectorPlane,
    ParticleDetectorPlane
  };

  /* Stopping uses this generic laboratory-frame extent instead of the magnet
   * container. Future field and particle detector planes can therefore join
   * the beamline as zero-length elements without changing stop semantics. */
  struct BeamlineElementExtent
  {
    BeamlineElementRole role;
    Double physicalEntrance;
    Double physicalExit;
    Double interactionEntrance;
    Double interactionExit;

    BeamlineElementExtent();
  };

  enum class StopMode
  {
    AfterLastElement,
    ReferenceCenterZ
  };

  struct StopConfig
  {
    StopMode mode;
    Double referenceZ;

    StopConfig();
  };

  enum class RuntimeMode
  {
    Interactive,
    Throughput
  };

  /* Runtime policy is global rather than belonging to one output.  This lets
   * a local field-only or detector-only test stop cleanly while the HPC path
   * avoids signal polling and durability flushes altogether. */
  struct RuntimeConfig
  {
    RuntimeMode mode;
    std::size_t stopCheckIntervalSteps;
    bool resourceReport;
    std::size_t resourceProgressIntervalSteps;
    unsigned int resourceCalibrationSteps;
    Double memorySafetyFactor;
    Double timeSafetyFactor;

    RuntimeConfig();
    bool interactive() const;
  };

  struct TrajectoryConfig
  {
    bool enabled;
    std::string directory;
    std::string basename;
    Double rhythm;
    std::size_t bufferRecords;
    unsigned int flushEverySamples;
    unsigned int compression;

    TrajectoryConfig();
  };

  struct EnergyLedgerConfig
  {
    bool enabled;
    std::string directory;
    std::string filename;
    std::size_t sampleIntervalSteps;
    std::size_t bufferRecords;
    unsigned int compression;
    Double warningRelativeTolerance;

    EnergyLedgerConfig();
  };

  /* Experimental diagnostic route for removing charged-particle current
   * upstream of a radiation field plane. It is deliberately opt-in because
   * changing a net charge inside the Maxwell domain is not continuity exact;
   * every production study must pair it with a zero-radiation baseline. */
  struct ParticleRetirementConfig
  {
    bool enabled;
    Double entranceZ;
    Double length;

    ParticleRetirementConfig();
    Double exitZ() const;
  };

  struct FieldDetectorPlaneConfig
  {
    std::string name;
    Double z;
    Double rhythm;
    std::size_t bufferSamples;
    unsigned int compression;

    /* The field plane owns a diagnostic-only, laboratory-frame reference
     * region immediately to its left.  These values are derived after the
     * input particles have been read: rho uses the full transverse mesh
     * diagonal and gamma is the largest laboratory particle gamma.  Crossing
     * the entrance records one compact straight-line reference state; it
     * never changes the physical particle or the Maxwell current. */
    bool particleBackgroundReference;
    Double referenceRho;
    Double referenceGamma;
    Double referenceDistance;
    Double referenceEntranceZ;
    std::size_t referenceBufferRecords;
    unsigned int referenceCompression;

    /* Optional small-particle validation.  When enabled, the detector also
     * captures the real particle at its sampling plane and compares that
     * crossing with the straight line launched at referenceEntranceZ. */
    bool referenceValidation;
    std::size_t referenceValidationMaximumParticles;

    /* Audit metadata copied into the field-plane file.  These values do not
     * control the retirement algorithm; the top-level configuration does. */
    bool particleRetirementCurrent;
    Double particleRetirementEntranceZ;
    Double particleRetirementExitZ;

    /* Optional detector-bound retirement design guard.  The user still sets
     * the global retirement geometry explicitly; these values only validate
     * that the C2 transition spans the requested number of cycles at the
     * lowest protected laboratory photon energy. */
    bool retirementFrequencyProtection;
    Double retirementMinimumPhotonEnergyEV;
    unsigned int retirementMinimumCycles;

    /* Derived after the laboratory particle snapshot is read.  The causal
     * guard uses the same conservative gamma*rho construction as the
     * ballistic-reference detector, but it performs no particle recording. */
    Double retirementCausalDistance;
    Double retirementCausalEntranceZ;
    Double retirementRequiredLength;
    Double retirementObservedCycles;

    FieldDetectorPlaneConfig();
  };

  struct ParticleDetectorPlaneConfig
  {
    std::string name;
    Double z;
    std::size_t bufferRecords;
    unsigned int compression;

    ParticleDetectorPlaneConfig();
  };

  /* Detector lists are empty by default.  The simulation uses that fact to
   * avoid constructing a manager, allocating buffers, or entering detector
   * MPI collectives when laboratory detector planes are not requested. */
  struct DetectorConfig
  {
    std::string directory;
    std::vector<FieldDetectorPlaneConfig> fieldPlanes;
    std::vector<ParticleDetectorPlaneConfig> particlePlanes;

    DetectorConfig();
    bool enabled() const;
  };

  struct WaveConfig
  {
    SIWaveSource source;
    SIFieldEvolution evolution;

    WaveConfig();
  };

  struct SimulationConfig
  {
    UnitSystem inputUnits;
    MeshConfig mesh;
    BoundaryConfig boundary;
    InitialSelfFieldConfig initialSelfField;
    BeamReferenceConfig reference;
    BeamInputConfig beam;
    std::vector<WaveConfig> waves;
    std::vector<SIMagneticElement> magnets;
    std::vector<BeamlineElementExtent> beamlineElements;
    RuntimeConfig runtime;
    StopConfig stop;
    TrajectoryConfig trajectory;
    EnergyLedgerConfig energyLedger;
    ParticleRetirementConfig particleRetirement;
    DetectorConfig detectors;
  };

  class YamlConfigLoader
  {
  public:
    static SimulationConfig loadFile(const std::string& filename);
  };
}

#endif
