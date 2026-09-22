#ifndef DIRECT_EB_CONFIG_H
#define DIRECT_EB_CONFIG_H

#include <cstddef>
#include <string>
#include <vector>

#include "eb_bunch.h"
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

  enum class BeamDistributionType
  {
    Ellipsoid,
    File,
    Manual,
    Crystal3D
  };

  struct BeamDistributionConfig
  {
    BeamDistributionType type;
    std::string profile;
    std::string generator;
    unsigned int randomSeed;
    std::size_t macroparticles;
    Double electrons;
    Double gamma;
    FieldVector<Double> direction;
    std::vector<FieldVector<Double> > positions;
    FieldVector<Double> sigmaPosition;
    FieldVector<Double> sigmaProperVelocity;
    Double transverseCutoff;
    Double longitudinalCutoff;
    std::string file;
    FieldVector<unsigned int> latticeCounts;
    FieldVector<Double> latticeConstants;
    Double bunchingFactor;
    Double bunchingPhase;
    bool shotNoise;

    BeamDistributionConfig();
  };

  struct BeamPlacementConfig
  {
    SIBunchPlacementMode mode;
    Double headDistance;

    BeamPlacementConfig();
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
    BeamPlacementConfig placement;
    std::vector<BeamDistributionConfig> beam;
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
