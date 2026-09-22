#include "config.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

namespace fel
{
  namespace
  {
    std::string lower(std::string value)
    {
      std::transform(value.begin(), value.end(), value.begin(),
        [](unsigned char character) {
          return static_cast<char>(std::tolower(character));
        });
      return value;
    }

    std::runtime_error configError(const YAML::Node& node,
                                   const std::string& message)
    {
      std::ostringstream stream;
      if (node.Mark().line >= 0)
        stream << "YAML line " << (node.Mark().line + 1) << ": ";
      stream << message;
      return std::runtime_error(stream.str());
    }

    YAML::Node required(const YAML::Node& parent, const char* key)
    {
      const YAML::Node value = parent[key];
      if (!value) throw configError(parent,
        std::string("missing required key '") + key + "'");
      return value;
    }

    Double finiteDouble(const YAML::Node& node, const std::string& name)
    {
      Double value = 0.0;
      try { value = node.as<Double>(); }
      catch (const YAML::Exception&) {
        throw configError(node, name + " must be a number");
      }
      if (!std::isfinite(value))
        throw configError(node, name + " must be finite");
      return value;
    }

    template<typename T>
    FieldVector<T> vector3(const YAML::Node& node, const std::string& name)
    {
      if (!node.IsSequence() || node.size() != 3)
        throw configError(node, name + " must contain exactly three values");
      FieldVector<T> result(static_cast<T>(0));
      for (std::size_t i = 0; i < 3; ++i)
        {
          try { result[static_cast<unsigned int>(i)] = node[i].as<T>(); }
          catch (const YAML::Exception&) {
            throw configError(node[i], name + " contains an invalid value");
          }
        }
      return result;
    }

    FieldVector<Double> finiteVector3(const YAML::Node& node,
                                      const std::string& name,
                                      Double scale = 1.0)
    {
      FieldVector<Double> result = vector3<Double>(node, name);
      for (unsigned int i = 0; i < 3; ++i)
        {
          if (!std::isfinite(result[i]))
            throw configError(node, name + " must contain finite values");
          result[i] *= scale;
        }
      return result;
    }

    Double lengthScale(const YAML::Node& node)
    {
      const std::string unit = lower(node.as<std::string>());
      if (unit == "m" || unit == "meter" || unit == "metre") return 1.0;
      if (unit == "mm" || unit == "millimeter" || unit == "millimetre") return 1.0e-3;
      if (unit == "um" || unit == "micrometer" || unit == "micrometre") return 1.0e-6;
      if (unit == "nm" || unit == "nanometer" || unit == "nanometre") return 1.0e-9;
      throw configError(node, "unsupported length unit '" + unit + "'");
    }

    Double timeScale(const YAML::Node& node)
    {
      const std::string unit = lower(node.as<std::string>());
      if (unit == "s" || unit == "second") return 1.0;
      if (unit == "ms" || unit == "millisecond") return 1.0e-3;
      if (unit == "us" || unit == "microsecond") return 1.0e-6;
      if (unit == "ns" || unit == "nanosecond") return 1.0e-9;
      if (unit == "ps" || unit == "picosecond") return 1.0e-12;
      if (unit == "fs" || unit == "femtosecond") return 1.0e-15;
      if (unit == "as" || unit == "attosecond") return 1.0e-18;
      throw configError(node, "unsupported time unit '" + unit + "'");
    }

    SIWaveProfile waveProfile(const YAML::Node& node)
    {
      const std::string value = lower(node.as<std::string>());
      if (value == "plane") return SIWaveProfile::Plane;
      if (value == "truncated-plane") return SIWaveProfile::TruncatedPlane;
      if (value == "gaussian") return SIWaveProfile::Gaussian;
      if (value == "super-gaussian") return SIWaveProfile::SuperGaussian;
      if (value == "standing-plane") return SIWaveProfile::StandingPlane;
      if (value == "standing-truncated-plane") return SIWaveProfile::StandingTruncatedPlane;
      if (value == "standing-gaussian") return SIWaveProfile::StandingGaussian;
      if (value == "standing-super-gaussian") return SIWaveProfile::StandingSuperGaussian;
      throw configError(node, "unsupported wave profile '" + value + "'");
    }

    SIEnvelopeType envelopeType(const YAML::Node& node)
    {
      const std::string value = lower(node.as<std::string>());
      if (value == "neumann") return SIEnvelopeType::Neumann;
      if (value == "gaussian") return SIEnvelopeType::Gaussian;
      if (value == "secant") return SIEnvelopeType::Secant;
      if (value == "flat-top") return SIEnvelopeType::FlatTop;
      if (value == "inverse-gaussian") return SIEnvelopeType::InverseGaussian;
      throw configError(node, "unsupported envelope type '" + value + "'");
    }

    Double electricAmplitude(const YAML::Node& node, Double frequency)
    {
      if (node["peak_electric_field_V_per_m"])
        return finiteDouble(node["peak_electric_field_V_per_m"],
                            "peak electric field");
      const Double a0 = finiteDouble(required(node, "normalized_amplitude"),
                                     "normalized amplitude");
      const Double twoPi = 6.283185307179586476925286766559;
      return a0 * SI::electronMass * SI::c * twoPi * frequency /
             SI::elementaryCharge;
    }

    BeamDistributionType beamType(const YAML::Node& node)
    {
      const std::string value = lower(node.as<std::string>());
      if (value == "ellipsoid") return BeamDistributionType::Ellipsoid;
      if (value == "file") return BeamDistributionType::File;
      if (value == "manual") return BeamDistributionType::Manual;
      if (value == "3d-crystal" || value == "crystal-3d")
        return BeamDistributionType::Crystal3D;
      throw configError(node, "unsupported beam distribution type '" + value + "'");
    }

    std::vector<FieldVector<Double> > positions(
        const YAML::Node& node, Double scale)
    {
      std::vector<FieldVector<Double> > result;
      if (!node)
        {
          result.push_back(FieldVector<Double>(0.0));
          return result;
        }
      if (!node.IsSequence())
        throw configError(node, "positions must be a sequence");
      if (node.size() == 3 && node[0].IsScalar())
        {
          result.push_back(finiteVector3(node, "position", scale));
          return result;
        }
      for (std::size_t i = 0; i < node.size(); ++i)
        result.push_back(finiteVector3(node[i], "position", scale));
      if (result.empty())
        throw configError(node, "positions cannot be empty");
      return result;
    }
  }

  UnitSystem::UnitSystem() : length(1.0), time(1.0) {}

  MeshConfig::MeshConfig()
    : lengths(0.0), resolution(0.0), center(0.0), duration(0.0),
      boostGamma(1.0), particleStepsPerUndulatorPeriod(1024)
  {}

  BeamDistributionConfig::BeamDistributionConfig()
    : type(BeamDistributionType::Ellipsoid), profile("gaussian"),
      generator("halton"), randomSeed(1), macroparticles(0), electrons(0.0),
      gamma(1.0), direction(0.0), positions(), sigmaPosition(0.0),
      sigmaProperVelocity(0.0),
      transverseCutoff(std::numeric_limits<Double>::max()),
      longitudinalCutoff(std::numeric_limits<Double>::max()),
      file(), latticeCounts(1),
      latticeConstants(0.0), bunchingFactor(0.0), bunchingPhase(0.0),
      shotNoise(false)
  {
    direction[2] = 1.0;
    positions.push_back(FieldVector<Double>(0.0));
  }

  BeamPlacementConfig::BeamPlacementConfig()
    : mode(SIBunchPlacementMode::AbsoluteLab), headDistance(0.0)
  {}

  TrajectoryConfig::TrajectoryConfig()
    : enabled(false), directory("./"), basename("trajectory"), rhythm(0.0),
      interactive(true), bufferRecords(16384), flushEverySamples(8),
      compression(0)
  {}

  WaveConfig::WaveConfig()
    : source(), evolution(SIFieldEvolution::MaxwellIncident)
  {}

  SimulationConfig YamlConfigLoader::loadFile(const std::string& filename)
  {
    YAML::Node root;
    try { root = YAML::LoadFile(filename); }
    catch (const YAML::Exception& error) {
      throw std::runtime_error(std::string("Cannot read YAML configuration: ") +
                               error.what());
    }
    if (!root.IsMap()) throw configError(root, "top level must be a map");

    SimulationConfig result;
    const YAML::Node units = root["units"];
    if (units)
      {
        result.inputUnits.length = lengthScale(required(units, "length"));
        result.inputUnits.time = timeScale(required(units, "time"));
      }

    const YAML::Node mesh = required(root, "mesh");
    result.mesh.lengths = finiteVector3(required(mesh, "lengths"),
                                        "mesh lengths", result.inputUnits.length);
    result.mesh.resolution = finiteVector3(required(mesh, "resolution"),
                                           "mesh resolution", result.inputUnits.length);
    result.mesh.center = mesh["center"] ?
      finiteVector3(mesh["center"], "mesh center", result.inputUnits.length) :
      FieldVector<Double>(0.0);
    result.mesh.duration = finiteDouble(required(mesh, "duration"),
                                        "mesh duration") * result.inputUnits.time;
    result.mesh.boostGamma = finiteDouble(required(mesh, "boost_gamma"),
                                          "boost gamma");
    if (mesh["particle_steps_per_undulator_period"])
      result.mesh.particleStepsPerUndulatorPeriod =
        mesh["particle_steps_per_undulator_period"].as<unsigned int>();

    const YAML::Node beam = required(root, "beam");
    const YAML::Node placement = required(beam, "placement");
    const std::string placementMode = lower(
      required(placement, "mode").as<std::string>());
    if (placementMode == "absolute-lab")
      result.placement.mode = SIBunchPlacementMode::AbsoluteLab;
    else if (placementMode == "head-to-first-element")
      {
        result.placement.mode = SIBunchPlacementMode::HeadToFirstElement;
        result.placement.headDistance = finiteDouble(
          required(placement, "distance"), "beam head distance") *
          result.inputUnits.length;
      }
    else
      throw configError(placement["mode"], "unsupported placement mode");

    const YAML::Node distributions = required(beam, "distributions");
    if (!distributions.IsSequence() || distributions.size() == 0)
      throw configError(distributions,
                        "beam distributions must be a nonempty sequence");
    for (std::size_t i = 0; i < distributions.size(); ++i)
      {
        const YAML::Node node = distributions[i];
        BeamDistributionConfig distribution;
        distribution.type = beamType(required(node, "type"));
        distribution.electrons = finiteDouble(required(node, "electrons"),
                                              "electron count");
        distribution.macroparticles = required(node, "macroparticles").as<std::size_t>();
        distribution.gamma = finiteDouble(required(node, "gamma"), "beam gamma");
        distribution.direction = finiteVector3(required(node, "direction"),
                                               "beam direction");
        distribution.positions = positions(node["positions"],
                                           result.inputUnits.length);
        if (node["profile"]) distribution.profile = lower(node["profile"].as<std::string>());
        if (node["generator"]) distribution.generator = lower(node["generator"].as<std::string>());
        if (node["random_seed"]) distribution.randomSeed = node["random_seed"].as<unsigned int>();
        if (node["sigma_position"])
          distribution.sigmaPosition = finiteVector3(node["sigma_position"],
            "beam position spread", result.inputUnits.length);
        if (node["sigma_proper_velocity"])
          distribution.sigmaProperVelocity = finiteVector3(
            node["sigma_proper_velocity"], "beam proper-velocity spread");
        if (node["transverse_cutoff"])
          distribution.transverseCutoff = finiteDouble(
            node["transverse_cutoff"], "transverse cutoff") *
            result.inputUnits.length;
        if (node["longitudinal_cutoff"])
          distribution.longitudinalCutoff = finiteDouble(
            node["longitudinal_cutoff"], "longitudinal cutoff") *
            result.inputUnits.length;
        if (node["file"]) distribution.file = node["file"].as<std::string>();
        if (node["lattice_counts"])
          distribution.latticeCounts = vector3<unsigned int>(
            node["lattice_counts"], "lattice counts");
        if (node["lattice_constants"])
          distribution.latticeConstants = finiteVector3(
            node["lattice_constants"], "lattice constants",
            result.inputUnits.length);
        if (node["bunching_factor"])
          distribution.bunchingFactor = finiteDouble(
            node["bunching_factor"], "bunching factor");
        if (node["bunching_phase_rad"])
          distribution.bunchingPhase = finiteDouble(
            node["bunching_phase_rad"], "bunching phase");
        if (node["shot_noise"])
          distribution.shotNoise = node["shot_noise"].as<bool>();
        result.beam.push_back(distribution);
      }

    const YAML::Node sources = root["sources"];
    if (sources && sources["incident_waves"])
      {
        const YAML::Node waves = sources["incident_waves"];
        if (!waves.IsSequence())
          throw configError(waves, "incident_waves must be a sequence");
        for (std::size_t i = 0; i < waves.size(); ++i)
          {
            const YAML::Node node = waves[i];
            WaveConfig wave;
            wave.source.profile = waveProfile(required(node, "profile"));
            wave.source.position = finiteVector3(required(node, "position"),
              "wave position", result.inputUnits.length);
            wave.source.direction = finiteVector3(required(node, "direction"),
                                                  "wave direction");
            wave.source.polarization = finiteVector3(
              required(node, "polarization"), "wave polarization");
            wave.source.wavelength = finiteDouble(
              required(node, "wavelength"), "wave wavelength") *
              result.inputUnits.length;
            wave.source.envelope.frequency = SI::c / wave.source.wavelength;
            wave.source.peakElectricField = electricAmplitude(
              node, wave.source.envelope.frequency);
            if (node["radius"])
              {
                const YAML::Node radius = node["radius"];
                if (!radius.IsSequence() || radius.size() != 2)
                  throw configError(radius, "wave radius must contain two values");
                wave.source.radius[0] = finiteDouble(radius[0], "wave radius") * result.inputUnits.length;
                wave.source.radius[1] = finiteDouble(radius[1], "wave radius") * result.inputUnits.length;
              }
            if (node["order"])
              {
                const YAML::Node order = node["order"];
                if (!order.IsSequence() || order.size() != 2)
                  throw configError(order, "wave order must contain two values");
                wave.source.order[0] = order[0].as<int>();
                wave.source.order[1] = order[1].as<int>();
              }
            const YAML::Node envelope = required(node, "envelope");
            wave.source.envelope.type = envelopeType(required(envelope, "type"));
            wave.source.envelope.centerTime = finiteDouble(
              required(envelope, "center_time"), "envelope center time") *
              result.inputUnits.time;
            wave.source.envelope.duration = finiteDouble(
              required(envelope, "duration"), "envelope duration") *
              result.inputUnits.time;
            if (envelope["carrier_phase_rad"])
              wave.source.envelope.carrierPhase = finiteDouble(
                envelope["carrier_phase_rad"], "carrier phase");
            if (envelope["rising_cycles"])
              wave.source.envelope.risingCycles = envelope["rising_cycles"].as<unsigned int>();
            if (envelope["inverse_gaussian_sigma"])
              {
                const YAML::Node sigma = envelope["inverse_gaussian_sigma"];
                if (!sigma.IsSequence() || sigma.size() != 2)
                  throw configError(sigma, "inverse Gaussian sigma must contain two values");
                wave.source.envelope.inverseGaussianSigma[0] = finiteDouble(sigma[0], "inverse Gaussian sigma") * result.inputUnits.time;
                wave.source.envelope.inverseGaussianSigma[1] = finiteDouble(sigma[1], "inverse Gaussian sigma") * result.inputUnits.time;
              }
            wave.source.prepare();
            result.waves.push_back(wave);
          }
      }

    if (sources && sources["magnetic_elements"])
      {
        const YAML::Node magnets = sources["magnetic_elements"];
        if (!magnets.IsSequence())
          throw configError(magnets, "magnetic_elements must be a sequence");
        for (std::size_t i = 0; i < magnets.size(); ++i)
          {
            const YAML::Node node = magnets[i];
            const std::string type = lower(required(node, "type").as<std::string>());
            SIMagneticElement magnet;
            const Double entrance = finiteDouble(required(node, "entrance_z"),
                                                 "magnet entrance") * result.inputUnits.length;
            const Double angle = node["polarization_angle_rad"] ?
              finiteDouble(node["polarization_angle_rad"], "magnet polarization angle") : 0.0;
            if (type == "planar-undulator")
              {
                magnet = SIMagneticElement::planarUndulatorFromK(
                  finiteDouble(required(node, "strength_parameter"), "undulator strength"),
                  finiteDouble(required(node, "period"), "undulator period") * result.inputUnits.length,
                  entrance, required(node, "periods").as<unsigned int>(), angle);
                if (node["gaussian_fringe"])
                  {
                    magnet.gaussianFringe = node["gaussian_fringe"].as<bool>();
                    magnet.prepare();
                  }
              }
            else if (type == "uniform-dipole")
              {
                magnet.type = SIMagnetType::UniformDipole;
                magnet.center[2] = entrance;
                magnet.length = finiteDouble(required(node, "length"), "dipole length") * result.inputUnits.length;
                magnet.peakMagneticField = finiteDouble(required(node, "field_T"), "dipole field");
                magnet.polarizationAngle = angle;
                magnet.gaussianFringe = false;
                magnet.prepare();
              }
            else
              throw configError(node["type"], "unsupported magnetic element type");
            result.magnets.push_back(magnet);
          }
      }

    const YAML::Node trajectory = root["trajectory"];
    if (trajectory)
      {
        result.trajectory.enabled = trajectory["enabled"] ?
          trajectory["enabled"].as<bool>() : true;
        if (trajectory["directory"]) result.trajectory.directory = trajectory["directory"].as<std::string>();
        if (trajectory["basename"]) result.trajectory.basename = trajectory["basename"].as<std::string>();
        result.trajectory.rhythm = finiteDouble(required(trajectory, "rhythm"),
                                                "trajectory rhythm") * result.inputUnits.time;
        if (trajectory["mode"])
          {
            const std::string mode = lower(trajectory["mode"].as<std::string>());
            if (mode == "interactive") result.trajectory.interactive = true;
            else if (mode == "throughput") result.trajectory.interactive = false;
            else throw configError(trajectory["mode"], "trajectory mode must be interactive or throughput");
          }
        if (trajectory["buffer_records"])
          result.trajectory.bufferRecords = trajectory["buffer_records"].as<std::size_t>();
        if (trajectory["flush_every_samples"])
          result.trajectory.flushEverySamples = trajectory["flush_every_samples"].as<unsigned int>();
        if (trajectory["compression"])
          result.trajectory.compression = trajectory["compression"].as<unsigned int>();
      }

    if (!(result.mesh.boostGamma >= 1.0))
      throw configError(mesh["boost_gamma"], "boost_gamma must be at least one");
    if (!(result.mesh.duration > 0.0))
      throw configError(mesh["duration"], "duration must be positive");
    if (result.placement.mode == SIBunchPlacementMode::HeadToFirstElement &&
        result.magnets.empty())
      throw configError(placement,
        "head-to-first-element placement requires a magnetic element");
    if (result.trajectory.enabled && !(result.trajectory.rhythm > 0.0))
      throw configError(trajectory, "enabled trajectory output requires positive rhythm");
    return result;
  }
}
