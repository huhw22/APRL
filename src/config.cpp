#include "config.h"

#include <algorithm>
#include <cctype>
#include <cmath>
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

    std::size_t positiveSize(const YAML::Node& node,
                             const std::string& name)
    {
      std::size_t value = 0;
      try { value = node.as<std::size_t>(); }
      catch (const YAML::Exception&) {
        throw configError(node, name + " must be a positive integer");
      }
      if (value == 0)
        throw configError(node, name + " must be a positive integer");
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

    bool absolutePath(const std::string& path)
    {
      return !path.empty() && path[0] == '/';
    }

    std::string resolveRelativeToConfig(const std::string& configFilename,
                                        const std::string& value)
    {
      if (value.empty() || absolutePath(value)) return value;
      const std::string::size_type separator =
        configFilename.find_last_of("/\\");
      if (separator == std::string::npos) return value;
      return configFilename.substr(0, separator + 1) + value;
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

    BeamInputType beamInputType(const YAML::Node& node)
    {
      const std::string value = lower(node.as<std::string>());
      if (value == "hdf5") return BeamInputType::Hdf5;
      if (value == "generated-gaussian")
        return BeamInputType::GeneratedGaussian;
      throw configError(node,
        "beam input type must be hdf5 or generated-gaussian");
    }
  }

  UnitSystem::UnitSystem() : length(1.0), time(1.0) {}

  MeshConfig::MeshConfig()
    : lengths(0.0), resolution(0.0), center(0.0), duration(0.0),
      boostGamma(1.0), particleStepsPerUndulatorPeriod(1024)
  {}

  BeamInputConfig::BeamInputConfig()
    : type(BeamInputType::Hdf5), file(), electrons(0.0), positionOffset(0.0),
      macroparticles(0), gamma(1.0), direction(0.0), center(0.0),
      sigmaPosition(0.0), sigmaProperVelocity(0.0), randomSeed(1)
  {
    direction[2] = 1.0;
  }

  BeamReferenceConfig::BeamReferenceConfig()
    : distanceToFirstMagnet(0.0)
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
    const YAML::Node reference = required(beam, "reference");
    result.reference.distanceToFirstMagnet = finiteDouble(
      required(reference, "distance_to_first_magnet"),
      "beam reference distance") * result.inputUnits.length;
    if (result.reference.distanceToFirstMagnet < 0.0)
      throw configError(reference["distance_to_first_magnet"],
        "distance_to_first_magnet cannot be negative");

    const YAML::Node input = required(beam, "input");
    result.beam.type = beamInputType(required(input, "type"));
    result.beam.electrons = finiteDouble(required(input, "electrons"),
                                         "electron count");
    if (!(result.beam.electrons > 0.0))
      throw configError(input["electrons"], "electron count must be positive");
    if (input["position_offset"])
      result.beam.positionOffset = finiteVector3(input["position_offset"],
        "beam position offset", result.inputUnits.length);

    if (result.beam.type == BeamInputType::Hdf5)
      {
        result.beam.file = resolveRelativeToConfig(filename,
          required(input, "file").as<std::string>());
        if (result.beam.file.empty())
          throw configError(input["file"], "particle HDF5 file cannot be empty");
      }
    else
      {
        result.beam.macroparticles = positiveSize(
          required(input, "macroparticles"), "macroparticles");
        result.beam.gamma = finiteDouble(required(input, "gamma"),
                                         "beam gamma");
        result.beam.direction = finiteVector3(required(input, "direction"),
                                               "beam direction");
        result.beam.center = input["center"] ?
          finiteVector3(input["center"], "beam center",
                        result.inputUnits.length) :
          FieldVector<Double>(0.0);
        result.beam.sigmaPosition = finiteVector3(
          required(input, "sigma_position"), "beam position spread",
          result.inputUnits.length);
        result.beam.sigmaProperVelocity = finiteVector3(
          required(input, "sigma_proper_velocity"),
          "beam proper-velocity spread");
        if (input["random_seed"])
          result.beam.randomSeed = input["random_seed"].as<unsigned int>();
        if (!(result.beam.gamma >= 1.0))
          throw configError(input["gamma"], "beam gamma must be at least one");
        if (!(result.beam.direction.norm() > 0.0))
          throw configError(input["direction"], "beam direction must be nonzero");
        for (unsigned int axis = 0; axis < 3; ++axis)
          if (result.beam.sigmaPosition[axis] < 0.0 ||
              result.beam.sigmaProperVelocity[axis] < 0.0)
            throw configError(input,
              "Gaussian standard deviations cannot be negative");
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
                  throw configError(radius,
                    "wave radius must contain two values");
                wave.source.radius[0] = finiteDouble(radius[0], "wave radius") * result.inputUnits.length;
                wave.source.radius[1] = finiteDouble(radius[1], "wave radius") * result.inputUnits.length;
              }
            if (node["order"])
              {
                const YAML::Node order = node["order"];
                if (!order.IsSequence() || order.size() != 2)
                  throw configError(order,
                    "wave order must contain two values");
                wave.source.order[0] = order[0].as<int>();
                wave.source.order[1] = order[1].as<int>();
              }
            const YAML::Node envelope = required(node, "envelope");
            wave.source.envelope.type = envelopeType(
              required(envelope, "type"));
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
              wave.source.envelope.risingCycles =
                envelope["rising_cycles"].as<unsigned int>();
            if (envelope["inverse_gaussian_sigma"])
              {
                const YAML::Node sigma = envelope["inverse_gaussian_sigma"];
                if (!sigma.IsSequence() || sigma.size() != 2)
                  throw configError(sigma,
                    "inverse Gaussian sigma must contain two values");
                wave.source.envelope.inverseGaussianSigma[0] =
                  finiteDouble(sigma[0], "inverse Gaussian sigma") *
                  result.inputUnits.time;
                wave.source.envelope.inverseGaussianSigma[1] =
                  finiteDouble(sigma[1], "inverse Gaussian sigma") *
                  result.inputUnits.time;
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
            const std::string type = lower(
              required(node, "type").as<std::string>());
            SIMagneticElement magnet;
            const Double entrance = finiteDouble(
              required(node, "entrance_z"), "magnet entrance") *
              result.inputUnits.length;
            const Double angle = node["polarization_angle_rad"] ?
              finiteDouble(node["polarization_angle_rad"],
                           "magnet polarization angle") : 0.0;
            if (type == "planar-undulator")
              {
                magnet = SIMagneticElement::planarUndulatorFromK(
                  finiteDouble(required(node, "strength_parameter"),
                               "undulator strength"),
                  finiteDouble(required(node, "period"),
                               "undulator period") * result.inputUnits.length,
                  entrance, required(node, "periods").as<unsigned int>(),
                  angle);
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
                magnet.length = finiteDouble(required(node, "length"),
                  "dipole length") * result.inputUnits.length;
                magnet.peakMagneticField = finiteDouble(
                  required(node, "field_T"), "dipole field");
                magnet.polarizationAngle = angle;
                magnet.gaussianFringe = false;
                magnet.prepare();
              }
            else
              throw configError(node["type"],
                "unsupported magnetic element type");
            result.magnets.push_back(magnet);
          }
      }

    const YAML::Node trajectory = root["trajectory"];
    if (trajectory)
      {
        result.trajectory.enabled = trajectory["enabled"] ?
          trajectory["enabled"].as<bool>() : true;
        if (trajectory["directory"])
          result.trajectory.directory =
            trajectory["directory"].as<std::string>();
        if (trajectory["basename"])
          result.trajectory.basename =
            trajectory["basename"].as<std::string>();
        if (trajectory["rhythm"])
          result.trajectory.rhythm = finiteDouble(trajectory["rhythm"],
            "trajectory rhythm") * result.inputUnits.time;
        if (trajectory["mode"])
          {
            const std::string mode = lower(
              trajectory["mode"].as<std::string>());
            if (mode == "interactive") result.trajectory.interactive = true;
            else if (mode == "throughput") result.trajectory.interactive = false;
            else throw configError(trajectory["mode"],
              "trajectory mode must be interactive or throughput");
          }
        if (trajectory["buffer_records"])
          result.trajectory.bufferRecords = positiveSize(
            trajectory["buffer_records"], "trajectory buffer_records");
        if (trajectory["flush_every_samples"])
          result.trajectory.flushEverySamples =
            trajectory["flush_every_samples"].as<unsigned int>();
        if (trajectory["compression"])
          result.trajectory.compression =
            trajectory["compression"].as<unsigned int>();
      }

    if (!(result.mesh.boostGamma >= 1.0))
      throw configError(mesh["boost_gamma"],
        "boost_gamma must be at least one");
    if (!(result.mesh.duration > 0.0))
      throw configError(mesh["duration"], "duration must be positive");
    if (result.magnets.empty())
      throw configError(sources,
        "beam reference placement requires at least one magnetic element");
    if (result.trajectory.enabled && !(result.trajectory.rhythm > 0.0))
      throw configError(trajectory,
        "enabled trajectory output requires positive rhythm");
    if (result.trajectory.compression > 9)
      throw configError(trajectory, "trajectory compression must be in [0,9]");
    return result;
  }
}
