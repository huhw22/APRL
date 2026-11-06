#ifndef DIRECT_EB_CONFIG_H
#define DIRECT_EB_CONFIG_H

#include <cstddef>
#include <string>
#include <vector>

#include "eb_field.h"
#include "eb_cpml.h"
#include "eb_sources.h"
#include "fieldvector.h"

namespace aprl
{
  enum class LongitudinalMeshMode
  {
    Explicit,
    ObservationBounded
  };

  struct LongitudinalMeshConfig
  {
    LongitudinalMeshMode mode;
    std::size_t lowerGuardCells;
    std::size_t upperGuardCells;
    Double selfFieldPaddingRms;
    Double ballisticTimeSafetyFactor;
    std::size_t alignmentCells;

    LongitudinalMeshConfig();
    bool automatic() const;
  };

  struct UnitSystem
  {
    Double length;
    Double time;

    UnitSystem();
  };

  /* Optional longitudinal window. The reference source domain is described
   * by the ordinary mesh; cells=0 selects runtime sizing of the active window.
   * Reference files follow the same working-directory path policy as output. */
  struct MovingWindowConfig
  {
    bool enabled;
    std::size_t cells;
    std::size_t lowerGuardCells;
    std::size_t upperGuardCells;
    std::size_t blockLayers;
    // Zero retains the native source-support collar. A common explicit
    // collar permits controlled source-shape comparisons at one exit plane.
    std::size_t particleCollarCells;
    // Extend the finite reference solve/fields through transverse CPML to
    // the existing outer PEC termination, instead of cutting at its entrance.
    bool referenceThroughCpml;
    bool referenceSurfaceCurrent;
    bool inletChargeMatch;
    std::string referenceDirectory;

    MovingWindowConfig();
  };

  struct MpiParticleBalanceConfig
  {
    bool enabled;
    std::size_t intervalSteps, minimumCells;
    Double minimumImprovement, frozenWeight, fieldCellWeight;
    MpiParticleBalanceConfig()
      :enabled(false),intervalSteps(5),minimumCells(3),
       minimumImprovement(.1),frozenWeight(.2),fieldCellWeight(0.) {}
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
    EBCurrentFilter currentFilter;
    /* Diagnostic control: keep evolving/depositing fields, but omit their
     * force on particles. Prescribed device fields remain active. */
    bool selfFieldFeedback;
    bool conservativeSourceCoupling;
    bool longitudinalQuintic;
    bool longitudinalNotch;
    bool longitudinalBinomial3;
    bool wideLongitudinalSource() const {return longitudinalQuintic || longitudinalNotch || longitudinalBinomial3;}
    /* Experimental finite-window boundary injection: particles outside the
     * transverse physical region remain ballistic until they enter it. */
    bool ballisticOutsideTransverse;
    LongitudinalMeshConfig longitudinal;
    MovingWindowConfig movingWindow;
    MpiParticleBalanceConfig particleBalance;

    MeshConfig();
  };

  /* Optional target-band preflight.  The target is an on-axis, forward
   * laboratory photon energy.  Only strict Nyquist failures are fatal; the
   * configurable quality levels below emit warnings and never change the
   * Maxwell or detector algorithms. */
  struct RadiationResolutionConfig
  {
    bool enabled;
    Double maximumPhotonEnergyEV;
    Double warningGridPointsPerWavelength;
    Double warningMaxwellSamplesPerCycle;
    Double warningDetectorSamplesPerCycle;

    RadiationResolutionConfig();
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

  enum class InitialSelfFieldModel
  {
    RelativisticPoisson,
    ElectrostaticPoisson,
    CowanRelativisticPoisson
  };

  struct InitialSelfFieldConfig
  {
    bool enabled;
    /* Include particles inside the transverse CPML in the one-time field
     * solve while keeping the evolving CPML field and memory at zero. */
    bool includeTransverseCpmlParticles;
    InitialSelfFieldModel model;
    Double relativeTolerance;
    Double maximumMeanTransverseBeta;
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
    /* Lab coordinate of the reconstructed snapshot's relative z=0 reference.
     * All laboratory elements and an optional Elegant input plane use the
     * same user-defined coordinate system; no beamline element is forced to
     * define z=0. */
    Double initialCenterZ;
    Double inputPlaneZ;
    bool inputPlaneZSet;
    bool frozenLabPlaneInjection;
    bool holdTransverseUntilActivation;
    Double activationPlaneZ;
    bool activationPlaneZSet;

    BeamReferenceConfig();
  };

  enum class BeamlineElementRole
  {
    MagneticDevice,
    FieldDetectorPlane,
    ParticleDetectorPlane,
    EnergyDetectorPlane
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
    bool shutdownWatchdogEnabled;
    unsigned int shutdownWatchdogTimeoutSeconds;
    unsigned int shutdownWatchdogFinalizeTimeoutSeconds;

    RuntimeConfig();
    bool interactive() const;
  };

  /* One explicit policy covers every file created by a simulation run.
   * Refusing replacement by default prevents an expensive batch rerun from
   * silently destroying an earlier result set. */
  struct OutputPolicyConfig
  {
    bool overwrite;
    std::string manifest;

    OutputPolicyConfig();
  };

  /* Optional synchronous full-state snapshots. Times are absolute box-frame
   * times in SI after YAML unit conversion; steps count completed field
   * updates. Restart provenance is filled by the loader of a saved state,
   * never accepted from the input card. */
  struct CheckpointConfig
  {
    std::string directory;
    std::vector<std::size_t> steps;
    std::vector<Double> times;
    std::string restartFrom;
    bool allowRetirementChange;
    bool resetDiagnostics;
    unsigned int ioTimeoutSeconds;
    std::size_t restartStep;
    Double restartTimeBox;

    CheckpointConfig();
    bool enabled() const;
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
    bool discreteOperatorAudit;
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

  /* Scalar transport diagnostics, independent of raw E/B and particle files. */
  struct EnergyDetectorPlaneConfig
  {
    std::string name;
    Double z;
    Double rhythm;
    std::size_t bufferSamples;
    unsigned int compression;
    EnergyDetectorPlaneConfig();
  };

  /* Detector lists are empty by default.  The simulation uses that fact to
   * avoid constructing a manager, allocating buffers, or entering detector
   * MPI collectives when laboratory detector planes are not requested. */
  struct DetectorConfig
  {
    std::string directory;
    std::vector<FieldDetectorPlaneConfig> fieldPlanes;
    std::vector<ParticleDetectorPlaneConfig> particlePlanes;
    std::vector<EnergyDetectorPlaneConfig> energyPlanes;
    /* Independent converged-substep records alongside legacy particle files. */
    bool substepParticlePlanes;
    /* Derived from the transverse CPML, never an independent input knob. */
    std::size_t physicalGuardCells[2];

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
    std::string configurationPath;
    std::string configurationText;
    std::string configurationDigest;
    UnitSystem inputUnits;
    MeshConfig mesh;
    RadiationResolutionConfig radiationResolution;
    BoundaryConfig boundary;
    InitialSelfFieldConfig initialSelfField;
    BeamReferenceConfig reference;
    BeamInputConfig beam;
    std::vector<WaveConfig> waves;
    std::vector<SIMagneticElement> magnets;
    std::vector<BeamlineElementExtent> beamlineElements;
    RuntimeConfig runtime;
    OutputPolicyConfig output;
    CheckpointConfig checkpoint;
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
