#include "config.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <set>
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

    std::size_t nonnegativeSize(const YAML::Node& node,
                                const std::string& name)
    {
      long long value = 0;
      try { value = node.as<long long>(); }
      catch (const YAML::Exception&) {
        throw configError(node, name + " must be a nonnegative integer");
      }
      if (value < 0)
        throw configError(node, name + " must be a nonnegative integer");
      return static_cast<std::size_t>(value);
    }

    void cpmlCells(const YAML::Node& node, std::size_t result[3])
    {
      if (!node.IsSequence() || node.size() != 3)
        throw configError(node,
          "boundary cells must contain exactly three integers");
      for (unsigned int axis = 0; axis < 3; ++axis)
        {
          result[axis] = nonnegativeSize(
            node[axis], "boundary cells");
          if (result[axis] == 1)
            throw configError(node[axis],
              "an enabled CPML direction needs at least two cells");
        }
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

    RuntimeMode runtimeMode(const YAML::Node& node)
    {
      const std::string value = lower(node.as<std::string>());
      if (value == "interactive" || value == "local-test")
        return RuntimeMode::Interactive;
      if (value == "throughput" || value == "hpc")
        return RuntimeMode::Throughput;
      throw configError(node,
        "runtime mode must be interactive or throughput");
    }

    EBMaxwellSolver fieldSolver(const YAML::Node& node)
    {
      const std::string value = lower(node.as<std::string>());
      if (value == "yee") return EBMaxwellSolver::Yee;
      if (value == "cowan-z" || value == "ckc-z" || value == "ckc")
        return EBMaxwellSolver::CowanZ;
      throw configError(node,
        "field_solver must be yee or cowan-z");
    }

    EBBoundaryType boundaryType(const YAML::Node& node)
    {
      const std::string value = lower(node.as<std::string>());
      if (value == "pec") return EBBoundaryType::Pec;
      if (value == "cpml") return EBBoundaryType::Cpml;
      throw configError(node, "boundary type must be pec or cpml");
    }

    bool validDetectorName(const std::string& name)
    {
      if (name.empty()) return false;
      for (std::size_t i = 0; i < name.size(); ++i)
        {
          const unsigned char character =
            static_cast<unsigned char>(name[i]);
          if (!std::isalnum(character) && character != '-' &&
              character != '_' && character != '.')
            return false;
        }
      return name != "." && name != "..";
    }
  }

  UnitSystem::UnitSystem() : length(1.0), time(1.0) {}

  MeshConfig::MeshConfig()
    : lengths(0.0), resolution(0.0), center(0.0), duration(0.0),
      boostGamma(1.0), particleStepsPerUndulatorPeriod(1024),
      fieldSolver(EBMaxwellSolver::CowanZ)
  {}

  BoundaryConfig::BoundaryConfig()
    : type(EBBoundaryType::Pec), cpml()
  {}

  BeamInputConfig::BeamInputConfig()
    : type(BeamInputType::Hdf5), file(), electrons(0.0), positionOffset(0.0),
      macroparticles(0), gamma(1.0), direction(0.0), center(0.0),
      sigmaPosition(0.0), sigmaProperVelocity(0.0), randomSeed(1)
  {
    direction[2] = 1.0;
  }

  BeamReferenceConfig::BeamReferenceConfig()
    : initialCenterZ(0.0)
  {}

  BeamlineElementExtent::BeamlineElementExtent()
    : role(BeamlineElementRole::MagneticDevice), physicalEntrance(0.0),
      physicalExit(0.0), interactionEntrance(0.0), interactionExit(0.0)
  {}

  StopConfig::StopConfig()
    : mode(StopMode::AfterLastElement), referenceZ(0.0)
  {}

  RuntimeConfig::RuntimeConfig()
    : mode(RuntimeMode::Throughput), stopCheckIntervalSteps(16)
  {}

  bool RuntimeConfig::interactive() const
  {
    return mode == RuntimeMode::Interactive;
  }

  TrajectoryConfig::TrajectoryConfig()
    : enabled(false), directory("./"), basename("trajectory"), rhythm(0.0),
      bufferRecords(16384), flushEverySamples(8), compression(0)
  {}

  FieldDetectorPlaneConfig::FieldDetectorPlaneConfig()
    : name(), z(0.0), rhythm(0.0), bufferSamples(1), compression(0)
  {}

  ParticleDetectorPlaneConfig::ParticleDetectorPlaneConfig()
    : name(), z(0.0), bufferRecords(16384), compression(0)
  {}

  DetectorConfig::DetectorConfig()
    : directory("output/detectors"), fieldPlanes(), particlePlanes()
  {}

  bool DetectorConfig::enabled() const
  {
    return !fieldPlanes.empty() || !particlePlanes.empty();
  }

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

    const YAML::Node runtime = root["runtime"];
    bool runtimeModeExplicit = false;
    if (runtime)
      {
        if (!runtime.IsMap())
          throw configError(runtime, "runtime must be a map");
        result.runtime.mode = runtimeMode(required(runtime, "mode"));
        runtimeModeExplicit = true;
        if (runtime["stop_check_interval_steps"])
          result.runtime.stopCheckIntervalSteps = positiveSize(
            runtime["stop_check_interval_steps"],
            "runtime stop_check_interval_steps");
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
    result.mesh.fieldSolver = fieldSolver(required(mesh, "field_solver"));
    if (mesh["particle_steps_per_undulator_period"])
      result.mesh.particleStepsPerUndulatorPeriod =
        mesh["particle_steps_per_undulator_period"].as<unsigned int>();

    const YAML::Node boundary = required(root, "boundary");
    if (!boundary.IsMap())
      throw configError(boundary, "boundary must be a map");
    result.boundary.type = boundaryType(required(boundary, "type"));
    if (result.boundary.type == EBBoundaryType::Cpml)
      {
        cpmlCells(required(boundary, "cells"),
                  result.boundary.cpml.cells);
        result.boundary.cpml.polynomialOrder = finiteDouble(
          required(boundary, "polynomial_order"),
          "CPML polynomial_order");
        result.boundary.cpml.targetReflection = finiteDouble(
          required(boundary, "target_reflection"),
          "CPML target_reflection");
        result.boundary.cpml.kappaMax = finiteDouble(
          required(boundary, "kappa_max"), "CPML kappa_max");
        result.boundary.cpml.alphaFraction = finiteDouble(
          required(boundary, "alpha_fraction"),
          "CPML alpha_fraction");
        if (!result.boundary.cpml.enabled())
          throw configError(boundary["cells"],
            "CPML needs at least one nonzero direction");
        if (!(result.boundary.cpml.polynomialOrder > 0.0))
          throw configError(boundary["polynomial_order"],
            "CPML polynomial_order must be positive");
        if (!(result.boundary.cpml.targetReflection > 0.0 &&
              result.boundary.cpml.targetReflection < 1.0))
          throw configError(boundary["target_reflection"],
            "CPML target_reflection must lie strictly between zero and one");
        if (!(result.boundary.cpml.kappaMax >= 1.0))
          throw configError(boundary["kappa_max"],
            "CPML kappa_max must be at least one");
        if (!(result.boundary.cpml.alphaFraction >= 0.0))
          throw configError(boundary["alpha_fraction"],
            "CPML alpha_fraction cannot be negative");
      }

    const YAML::Node beam = required(root, "beam");
    const YAML::Node reference = required(beam, "reference");
    result.reference.initialCenterZ = finiteDouble(
      required(reference, "initial_center_z"),
      "initial beam-center z") * result.inputUnits.length;

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
                  magnet.gaussianFringe = node["gaussian_fringe"].as<bool>();
                if (node["fringe_relative_cutoff"])
                  magnet.fringeRelativeCutoff = finiteDouble(
                    node["fringe_relative_cutoff"],
                    "undulator fringe relative cutoff");
                magnet.prepare();
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
            BeamlineElementExtent extent;
            extent.role = BeamlineElementRole::MagneticDevice;
            extent.physicalEntrance = magnet.physicalEntranceLab();
            extent.physicalExit = magnet.physicalExitLab();
            extent.interactionEntrance = magnet.interactionEntranceLab();
            extent.interactionExit = magnet.interactionExitLab();
            result.beamlineElements.push_back(extent);
          }
      }

    const YAML::Node detectors = root["detectors"];
    std::set<std::string> detectorNames;
    if (detectors)
      {
        if (!detectors.IsMap())
          throw configError(detectors, "detectors must be a map");
        const bool detectorsEnabled = detectors["enabled"] ?
          detectors["enabled"].as<bool>() : true;
        if (detectorsEnabled && detectors["directory"])
          result.detectors.directory =
            detectors["directory"].as<std::string>();
        if (detectorsEnabled && result.detectors.directory.empty())
          throw configError(detectors, "detector directory cannot be empty");

        const YAML::Node fieldPlanes = detectors["field_planes"];
        if (detectorsEnabled && fieldPlanes)
          {
            if (!fieldPlanes.IsSequence())
              throw configError(fieldPlanes,
                "field_planes must be a sequence");
            for (std::size_t i = 0; i < fieldPlanes.size(); ++i)
              {
                const YAML::Node node = fieldPlanes[i];
                FieldDetectorPlaneConfig plane;
                plane.name = required(node, "name").as<std::string>();
                if (!validDetectorName(plane.name))
                  throw configError(node["name"],
                    "detector name may contain only letters, digits, '.', '-' and '_'");
                if (!detectorNames.insert(plane.name).second)
                  throw configError(node["name"],
                    "detector names must be unique");
                plane.z = finiteDouble(required(node, "z"),
                  "field detector z") * result.inputUnits.length;
                plane.rhythm = finiteDouble(required(node, "rhythm"),
                  "field detector rhythm") * result.inputUnits.time;
                if (!(plane.rhythm > 0.0))
                  throw configError(node["rhythm"],
                    "field detector rhythm must be positive");
                if (node["buffer_samples"])
                  plane.bufferSamples = positiveSize(
                    node["buffer_samples"], "field detector buffer_samples");
                if (node["compression"])
                  plane.compression = node["compression"].as<unsigned int>();
                if (plane.compression > 9)
                  throw configError(node,
                    "field detector compression must be in [0,9]");
                result.detectors.fieldPlanes.push_back(plane);

                BeamlineElementExtent extent;
                extent.role = BeamlineElementRole::FieldDetectorPlane;
                extent.physicalEntrance = extent.physicalExit = plane.z;
                extent.interactionEntrance = extent.interactionExit = plane.z;
                result.beamlineElements.push_back(extent);
              }
          }

        const YAML::Node particlePlanes = detectors["particle_planes"];
        if (detectorsEnabled && particlePlanes)
          {
            if (!particlePlanes.IsSequence())
              throw configError(particlePlanes,
                "particle_planes must be a sequence");
            for (std::size_t i = 0; i < particlePlanes.size(); ++i)
              {
                const YAML::Node node = particlePlanes[i];
                ParticleDetectorPlaneConfig plane;
                plane.name = required(node, "name").as<std::string>();
                if (!validDetectorName(plane.name))
                  throw configError(node["name"],
                    "detector name may contain only letters, digits, '.', '-' and '_'");
                if (!detectorNames.insert(plane.name).second)
                  throw configError(node["name"],
                    "detector names must be unique");
                plane.z = finiteDouble(required(node, "z"),
                  "particle detector z") * result.inputUnits.length;
                if (node["buffer_records"])
                  plane.bufferRecords = positiveSize(
                    node["buffer_records"],
                    "particle detector buffer_records");
                if (node["compression"])
                  plane.compression = node["compression"].as<unsigned int>();
                if (plane.compression > 9)
                  throw configError(node,
                    "particle detector compression must be in [0,9]");
                result.detectors.particlePlanes.push_back(plane);

                BeamlineElementExtent extent;
                extent.role = BeamlineElementRole::ParticleDetectorPlane;
                extent.physicalEntrance = extent.physicalExit = plane.z;
                extent.interactionEntrance = extent.interactionExit = plane.z;
                result.beamlineElements.push_back(extent);
              }
          }
      }

    const YAML::Node stop = required(root, "stop");
    const std::string stopMode = lower(required(stop, "mode").as<std::string>());
    if (stopMode == "after-last-element")
      result.stop.mode = StopMode::AfterLastElement;
    else if (stopMode == "reference-center-z")
      {
        result.stop.mode = StopMode::ReferenceCenterZ;
        result.stop.referenceZ = finiteDouble(required(stop, "z"),
          "stop reference-center z") * result.inputUnits.length;
      }
    else
      throw configError(stop["mode"],
        "stop mode must be after-last-element or reference-center-z");

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
            const RuntimeMode legacyMode = runtimeMode(trajectory["mode"]);
            if (runtimeModeExplicit && legacyMode != result.runtime.mode)
              throw configError(trajectory["mode"],
                "trajectory.mode conflicts with global runtime.mode");
            if (!runtimeModeExplicit) result.runtime.mode = legacyMode;
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
    if (result.beamlineElements.empty())
      throw configError(sources,
        "beam placement and stopping require at least one beamline element");
    Double firstPhysical = result.beamlineElements[0].physicalEntrance;
    Double lastInteraction = result.beamlineElements[0].interactionExit;
    for (std::size_t i = 1; i < result.beamlineElements.size(); ++i)
      {
        firstPhysical = std::min(firstPhysical,
          result.beamlineElements[i].physicalEntrance);
        lastInteraction = std::max(lastInteraction,
          result.beamlineElements[i].interactionExit);
      }
    const Double originTolerance = 64.0 *
      std::numeric_limits<Double>::epsilon() *
      std::max(1.0, std::abs(firstPhysical));
    if (std::abs(firstPhysical) > originTolerance)
      throw configError(sources,
        "the first physical beamline entrance must define lab z=0");
    if (result.stop.mode == StopMode::ReferenceCenterZ &&
        !(result.stop.referenceZ > lastInteraction))
      throw configError(stop["z"],
        "reference-center stop z must lie beyond every element interaction region");
    if (result.trajectory.enabled && !(result.trajectory.rhythm > 0.0))
      throw configError(trajectory,
        "enabled trajectory output requires positive rhythm");
    if (result.trajectory.compression > 9)
      throw configError(trajectory, "trajectory compression must be in [0,9]");
    if (result.runtime.interactive() && result.trajectory.enabled &&
        result.trajectory.flushEverySamples == 0)
      throw configError(trajectory,
        "interactive trajectory flush_every_samples must be positive");
    return result;
  }
}
