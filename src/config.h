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
    /* Distance from the particle file's relative z=0 reference to the first
     * magnetic-element entrance, measured in the laboratory frame. */
    Double distanceToFirstMagnet;

    BeamReferenceConfig();
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
    TrajectoryConfig trajectory;
  };

  class YamlConfigLoader
  {
  public:
    static SimulationConfig loadFile(const std::string& filename);
  };
}

#endif
