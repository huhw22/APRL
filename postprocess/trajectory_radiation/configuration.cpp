#include "configuration.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

#include "yaml-cpp/yaml.h"

namespace radiation
{
  namespace
  {
    std::string directoryName(const std::string& path)
    {
      const std::string::size_type separator = path.find_last_of("/\\");
      if (separator == std::string::npos) return ".";
      if (separator == 0) return path.substr(0, 1);
      return path.substr(0, separator);
    }

    bool isAbsolute(const std::string& path)
    {
      return !path.empty() && (path[0] == '/' || path[0] == '\\');
    }

    std::string resolvePath(const std::string& base,
                            const std::string& path)
    {
      if (isAbsolute(path)) return path;
      return base + "/" + path;
    }

    YAML::Node required(const YAML::Node& node, const char* key)
    {
      if (!node || !node[key])
        throw std::invalid_argument(std::string("Missing YAML key: ") + key);
      return node[key];
    }

    Vec3 readVector(const YAML::Node& node, const char* description)
    {
      if (!node || !node.IsSequence() || node.size() != 3)
        throw std::invalid_argument(
          std::string(description) + " must contain three values");
      return Vec3(node[0].as<double>(), node[1].as<double>(),
                  node[2].as<double>());
    }

    AxisDefinition readAxis(const YAML::Node& node,
                            const char* description,
                            bool allowLogarithmic)
    {
      if (!node || !node.IsMap())
        throw std::invalid_argument(
          std::string(description) + " must be a mapping");
      AxisDefinition axis;
      axis.minimum = required(node, "min").as<double>();
      axis.maximum = required(node, "max").as<double>();
      axis.count = required(node, "count").as<std::size_t>();
      const std::string spacing = node["spacing"] ?
        node["spacing"].as<std::string>() : "linear";
      if (spacing != "linear" && spacing != "log")
        throw std::invalid_argument(
          std::string(description) + " spacing must be linear or log");
      axis.logarithmic = spacing == "log";
      if (axis.logarithmic && !allowLogarithmic)
        throw std::invalid_argument(
          std::string(description) + " supports only linear spacing");
      if (axis.count == 0 || !std::isfinite(axis.minimum) ||
          !std::isfinite(axis.maximum) || axis.maximum < axis.minimum ||
          (axis.count > 1 && axis.maximum == axis.minimum) ||
          (axis.logarithmic && !(axis.minimum > 0.0)))
        throw std::invalid_argument(
          std::string("Invalid axis definition: ") + description);
      return axis;
    }

    std::vector<double> readDoubleList(const YAML::Node& node,
                                       const char* description)
    {
      std::vector<double> result;
      if (!node) return result;
      if (!node.IsSequence())
        throw std::invalid_argument(
          std::string(description) + " must be a sequence");
      result.reserve(node.size());
      for (std::size_t index = 0; index < node.size(); ++index)
        {
          const double value = node[index].as<double>();
          if (!(value > 0.0) || !std::isfinite(value))
            throw std::invalid_argument(
              std::string(description) + " values must be positive");
          result.push_back(value);
        }
      return result;
    }

    std::vector<Vec3> readAngleList(const YAML::Node& node,
                                    const char* description)
    {
      std::vector<Vec3> result;
      if (!node) return result;
      if (!node.IsSequence())
        throw std::invalid_argument(
          std::string(description) + " must be a sequence");
      result.reserve(node.size());
      for (std::size_t index = 0; index < node.size(); ++index)
        {
          if (!node[index].IsSequence() || node[index].size() != 2)
            throw std::invalid_argument(
              std::string(description) +
              " entries must be [theta_x_rad, theta_y_rad]");
          const double x = node[index][0].as<double>();
          const double y = node[index][1].as<double>();
          if (!std::isfinite(x) || !std::isfinite(y))
            throw std::invalid_argument(
              std::string(description) + " contains a non-finite angle");
          result.push_back(Vec3(x, y, 0.0));
        }
      return result;
    }
  }

  Vec3::Vec3() : value{0.0L, 0.0L, 0.0L} {}

  Vec3::Vec3(long double x, long double y, long double z)
    : value{x, y, z} {}

  long double& Vec3::operator[](std::size_t index)
  {
    return value[index];
  }

  long double Vec3::operator[](std::size_t index) const
  {
    return value[index];
  }

  Vec3 operator+(const Vec3& left, const Vec3& right)
  {
    return Vec3(left[0] + right[0], left[1] + right[1],
                left[2] + right[2]);
  }

  Vec3 operator-(const Vec3& left, const Vec3& right)
  {
    return Vec3(left[0] - right[0], left[1] - right[1],
                left[2] - right[2]);
  }

  Vec3 operator*(long double scalar, const Vec3& vector)
  {
    return Vec3(scalar * vector[0], scalar * vector[1],
                scalar * vector[2]);
  }

  Vec3 operator/(const Vec3& vector, long double scalar)
  {
    return (1.0L / scalar) * vector;
  }

  long double dot(const Vec3& left, const Vec3& right)
  {
    return left[0] * right[0] + left[1] * right[1] +
      left[2] * right[2];
  }

  Vec3 cross(const Vec3& left, const Vec3& right)
  {
    return Vec3(left[1] * right[2] - left[2] * right[1],
                left[2] * right[0] - left[0] * right[2],
                left[0] * right[1] - left[1] * right[0]);
  }

  long double norm(const Vec3& vector)
  {
    return std::sqrt(dot(vector, vector));
  }

  Vec3 normalized(const Vec3& vector, const char* description)
  {
    const long double length = norm(vector);
    if (!(length > 0.0L) || !std::isfinite(length))
      throw std::invalid_argument(
        std::string(description) + " must be finite and nonzero");
    return vector / length;
  }

  AxisDefinition::AxisDefinition()
    : minimum(0.0), maximum(0.0), count(0), logarithmic(false)
  {}

  std::vector<double> AxisDefinition::values() const
  {
    std::vector<double> result(count, minimum);
    if (count <= 1) return result;
    if (logarithmic)
      {
        const long double lower = std::log(
          static_cast<long double>(minimum));
        const long double upper = std::log(
          static_cast<long double>(maximum));
        for (std::size_t index = 0; index < count; ++index)
          result[index] = static_cast<double>(std::exp(
            lower + (upper - lower) * static_cast<long double>(index) /
            static_cast<long double>(count - 1)));
      }
    else
      for (std::size_t index = 0; index < count; ++index)
        result[index] = minimum + (maximum - minimum) *
          static_cast<double>(index) / static_cast<double>(count - 1);
    return result;
  }

  RadiationConfig::RadiationConfig()
    : cardPath(), shots(), requireComplete(true),
      readChunkRecords(65536), observationAxis(0.0, 0.0, 1.0),
      horizontalAxis(1.0, 0.0, 0.0),
      verticalAxis(0.0, 1.0, 0.0), distanceM(1.0), thetaX(), thetaY(),
      photonEnergyEV(), frequencyBlock(16), thetaYBlock(4),
      minimumRecordsPerParticle(3), outputFile(), compression(0),
      coherence()
  {}

  TrajectoryLoadStats::TrajectoryLoadStats()
    : inputRecords(0), uniqueRecords(0), particles(0),
      terminalEvents(0), duplicateRecords(0)
  {}

  std::size_t RadiationBlock::scalarIndex(
      std::size_t frequency, std::size_t thetaYValue,
      std::size_t thetaXValue, std::size_t polarization,
      std::size_t realOrImag) const
  {
    return (((frequency * thetaYCount + thetaYValue) * thetaXCount +
      thetaXValue) * 2 + polarization) * 2 + realOrImag;
  }

  RadiationConfig loadConfiguration(const std::string& filename)
  {
    const YAML::Node root = YAML::LoadFile(filename);
    RadiationConfig config;
    config.cardPath = filename;
    const std::string base = directoryName(filename);

    const YAML::Node input = required(root, "input");
    config.requireComplete = input["require_complete"] ?
      input["require_complete"].as<bool>() : true;
    config.readChunkRecords = input["read_chunk_records"] ?
      input["read_chunk_records"].as<std::size_t>() : 65536;
    if (config.readChunkRecords == 0)
      throw std::invalid_argument("read_chunk_records must be positive");
    const YAML::Node shots = required(input, "shots");
    if (!shots.IsSequence() || shots.size() == 0)
      throw std::invalid_argument("input.shots must be a nonempty sequence");
    std::set<std::string> shotNames;
    for (std::size_t shotIndex = 0; shotIndex < shots.size(); ++shotIndex)
      {
        ShotConfig shot;
        shot.name = required(shots[shotIndex], "name").as<std::string>();
        const YAML::Node files = required(shots[shotIndex], "files");
        if (shot.name.empty() || !files.IsSequence() || files.size() == 0)
          throw std::invalid_argument(
            "Every shot needs a name and at least one trajectory file");
        if (!shotNames.insert(shot.name).second)
          throw std::invalid_argument(
            "Every shot name must be unique: " + shot.name);
        for (std::size_t file = 0; file < files.size(); ++file)
          shot.files.push_back(resolvePath(base, files[file].as<std::string>()));
        config.shots.push_back(shot);
      }

    const YAML::Node observation = required(root, "observation");
    config.observationAxis = normalized(readVector(
      required(observation, "axis"), "observation.axis"),
      "observation.axis");
    Vec3 requestedHorizontal = readVector(
      required(observation, "horizontal"), "observation.horizontal");
    requestedHorizontal = requestedHorizontal -
      dot(requestedHorizontal, config.observationAxis) *
      config.observationAxis;
    config.horizontalAxis = normalized(
      requestedHorizontal, "projected observation.horizontal");
    config.verticalAxis = normalized(cross(
      config.observationAxis, config.horizontalAxis),
      "observation vertical basis");
    config.distanceM = required(observation, "distance_m").as<double>();
    if (!(config.distanceM > 0.0) || !std::isfinite(config.distanceM))
      throw std::invalid_argument("observation.distance_m must be positive");
    config.thetaX = readAxis(required(observation, "theta_x_rad"),
                             "observation.theta_x_rad", false);
    config.thetaY = readAxis(required(observation, "theta_y_rad"),
                             "observation.theta_y_rad", false);
    if (std::max(std::abs(config.thetaX.minimum),
                 std::abs(config.thetaX.maximum)) >= 1.4 ||
        std::max(std::abs(config.thetaY.minimum),
                 std::abs(config.thetaY.maximum)) >= 1.4)
      throw std::invalid_argument(
        "Transverse angular coordinates must stay below 1.4 rad");

    const YAML::Node spectrum = required(root, "spectrum");
    config.photonEnergyEV = readAxis(
      required(spectrum, "photon_energy_eV"),
      "spectrum.photon_energy_eV", true);
    if (!(config.photonEnergyEV.minimum > 0.0))
      throw std::invalid_argument("Photon energy must be positive");

    const YAML::Node calculation = root["calculation"];
    if (calculation)
      {
        if (calculation["frequency_block"])
          config.frequencyBlock =
            calculation["frequency_block"].as<std::size_t>();
        if (calculation["theta_y_block"])
          config.thetaYBlock =
            calculation["theta_y_block"].as<std::size_t>();
        if (calculation["minimum_records_per_particle"])
          config.minimumRecordsPerParticle = calculation[
            "minimum_records_per_particle"].as<std::size_t>();
      }
    if (config.frequencyBlock == 0 || config.thetaYBlock == 0 ||
        config.minimumRecordsPerParticle < 3)
      throw std::invalid_argument(
        "Calculation blocks must be positive and minimum records >= 3");

    const YAML::Node output = required(root, "output");
    config.outputFile = resolvePath(base,
      required(output, "file").as<std::string>());
    config.compression = output["compression"] ?
      output["compression"].as<unsigned int>() : 0;
    if (config.compression > 9)
      throw std::invalid_argument("HDF5 compression must be between 0 and 9");

    const YAML::Node coherence = root["coherence"];
    if (coherence)
      {
        config.coherence.spatialPhotonEnergyEV = readDoubleList(
          coherence["spatial_photon_energy_eV"],
          "coherence.spatial_photon_energy_eV");
        config.coherence.spatialReferenceAngles = readAngleList(
          coherence["spatial_reference_angles_rad"],
          "coherence.spatial_reference_angles_rad");
        config.coherence.temporalPhotonEnergyEV = readDoubleList(
          coherence["temporal_photon_energy_eV"],
          "coherence.temporal_photon_energy_eV");
        config.coherence.temporalReferenceAngles = readAngleList(
          coherence["temporal_reference_angles_rad"],
          "coherence.temporal_reference_angles_rad");
      }
    if (config.coherence.spatialPhotonEnergyEV.empty() !=
        config.coherence.spatialReferenceAngles.empty())
      throw std::invalid_argument(
        "Spatial coherence needs both energies and reference angles");
    if (config.coherence.temporalPhotonEnergyEV.empty() !=
        config.coherence.temporalReferenceAngles.empty())
      throw std::invalid_argument(
        "Temporal coherence needs both energies and reference angles");

    return config;
  }
}
