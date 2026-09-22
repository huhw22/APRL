#ifndef DIRECT_EB_CONFIG_H
#define DIRECT_EB_CONFIG_H

#include <cstddef>
#include <string>
#include <vector>

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
    FieldVector<Double> lengths;
    FieldVector<Double> resolution;
    FieldVector<Double> center;
    Double duration;
    Double boostGamma;
    unsigned int particleStepsPerUndulatorPeriod;

    MeshConfig();
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

  struct TrajectoryConfig
  {
    bool enabled;
    std::string directory;
    std::string basename;
    Double rhythm;
    bool interactive;
    std::size_t bufferRecords;
    unsigned int flushEverySamples;
    unsigned int compression;

    TrajectoryConfig();
  };

  struct FieldDetectorPlaneConfig
  {
    std::string name;
    Double z;
    Double rhythm;
    std::size_t bufferSamples;
    unsigned int compression;

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
    BeamReferenceConfig reference;
    BeamInputConfig beam;
    std::vector<WaveConfig> waves;
    std::vector<SIMagneticElement> magnets;
    std::vector<BeamlineElementExtent> beamlineElements;
    StopConfig stop;
    TrajectoryConfig trajectory;
    DetectorConfig detectors;
  };

  class YamlConfigLoader
  {
  public:
    static SimulationConfig loadFile(const std::string& filename);
  };
}

#endif
