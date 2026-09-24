#include "configuration.h"

#include <cmath>
#include <limits>
#include <stdexcept>

#include "yaml-cpp/yaml.h"

namespace field_analysis
{
  namespace
  {
    std::string parentDirectory(const std::string& path)
    {
      const std::string::size_type separator = path.find_last_of('/');
      if (separator == std::string::npos) return std::string();
      if (separator == 0) return "/";
      return path.substr(0, separator);
    }

    std::string resolvePath(const std::string& card,
                            const std::string& path)
    {
      if (!path.empty() && path[0] == '/') return path;
      const std::string parent = parentDirectory(card);
      return parent.empty() ? path : parent + "/" + path;
    }

    YAML::Node required(const YAML::Node& parent, const char* key)
    {
      const YAML::Node value = parent[key];
      if (!value)
        throw std::runtime_error(std::string("Missing required YAML key '") +
          key + "'");
      return value;
    }

    double finiteDouble(const YAML::Node& node, const char* description)
    {
      const double value = node.as<double>();
      if (!std::isfinite(value))
        throw std::runtime_error(std::string(description) +
          " must be finite");
      return value;
    }

    std::size_t positiveSize(const YAML::Node& node,
                             const char* description)
    {
      const unsigned long long value = node.as<unsigned long long>();
      if (value == 0 || value >
          static_cast<unsigned long long>(
            std::numeric_limits<std::size_t>::max()))
        throw std::runtime_error(std::string(description) +
          " must be a positive representable integer");
      return static_cast<std::size_t>(value);
    }

    std::vector<double> positiveEnergyList(const YAML::Node& node,
                                            const char* description)
    {
      std::vector<double> result;
      if (!node) return result;
      if (!node.IsSequence())
        throw std::runtime_error(std::string(description) +
          " must be a sequence");
      for (std::size_t index = 0; index < node.size(); ++index)
        {
          const double value = finiteDouble(node[index], description);
          if (!(value > 0.0))
            throw std::runtime_error(std::string(description) +
              " values must be positive");
          result.push_back(value);
        }
      return result;
    }

    std::vector<ReferenceAngle> referenceAngles(
        const YAML::Node& node, const char* description)
    {
      std::vector<ReferenceAngle> result;
      if (!node) return result;
      if (!node.IsSequence())
        throw std::runtime_error(std::string(description) +
          " must be a sequence of [theta_x, theta_y]");
      for (std::size_t index = 0; index < node.size(); ++index)
        {
          if (!node[index].IsSequence() || node[index].size() != 2)
            throw std::runtime_error(std::string(description) +
              " entries must be [theta_x, theta_y]");
          ReferenceAngle angle;
          angle.thetaX = finiteDouble(node[index][0], description);
          angle.thetaY = finiteDouble(node[index][1], description);
          result.push_back(angle);
        }
      return result;
    }
  }

  Configuration::Configuration()
    : requireComplete(true), minimumPhotonEnergyEV(0.0),
      maximumPhotonEnergyEV(0.0), hannTimeWindows(false),
      intervalStart(0.0), intervalEnd(0.0), windowDuration(0.0),
      windowStep(0.0), transverseWindow("none"), angularPaddingY(1),
      angularPaddingX(1), frequencyBlock(8), spatialBatchPoints(32),
      fftThreads(1), maximumWorkingMiB(4096), maximumOutputGiB(64.0),
      compression(0)
  {}

  Configuration loadConfiguration(const std::string& filename)
  {
    const YAML::Node root = YAML::LoadFile(filename);
    if (!root.IsMap())
      throw std::runtime_error("Field-plane analysis card must be a map");
    const YAML::Node input = required(root, "input");
    const YAML::Node analysis = required(root, "analysis");
    const YAML::Node calculation = required(root, "calculation");
    const YAML::Node output = required(root, "output");
    if (!input.IsMap() || !analysis.IsMap() || !calculation.IsMap() ||
        !output.IsMap())
      throw std::runtime_error(
        "input, analysis, calculation, and output must be maps");

    Configuration config;
    config.cardPath = filename;
    config.fieldFile = resolvePath(filename,
      required(input, "field_file").as<std::string>());
    if (input["zero_radiation_baseline"])
      config.baselineFile = resolvePath(filename,
        input["zero_radiation_baseline"].as<std::string>());
    config.requireComplete = input["require_complete"] ?
      input["require_complete"].as<bool>() : true;
    if (config.fieldFile == config.baselineFile &&
        !config.baselineFile.empty())
      throw std::runtime_error(
        "field_file and zero_radiation_baseline must differ");

    const YAML::Node band = required(analysis, "photon_energy_band_eV");
    if (!band.IsSequence() || band.size() != 2)
      throw std::runtime_error(
        "analysis.photon_energy_band_eV must be [minimum, maximum]");
    config.minimumPhotonEnergyEV = finiteDouble(band[0],
      "minimum photon energy");
    config.maximumPhotonEnergyEV = finiteDouble(band[1],
      "maximum photon energy");
    if (!(config.minimumPhotonEnergyEV >= 0.0) ||
        !(config.maximumPhotonEnergyEV > config.minimumPhotonEnergyEV))
      throw std::runtime_error("Invalid photon-energy band");
    config.transverseWindow = analysis["transverse_window"] ?
      analysis["transverse_window"].as<std::string>() : "none";
    if (config.transverseWindow != "none" &&
        config.transverseWindow != "hann")
      throw std::runtime_error(
        "analysis.transverse_window must be none or hann");
    if (analysis["angular_zero_padding"])
      {
        const YAML::Node padding = analysis["angular_zero_padding"];
        if (!padding.IsSequence() || padding.size() != 2)
          throw std::runtime_error(
            "analysis.angular_zero_padding must be [y, x]");
        config.angularPaddingY = positiveSize(padding[0],
          "angular y padding");
        config.angularPaddingX = positiveSize(padding[1],
          "angular x padding");
      }

    if (analysis["time_windows"])
      {
        const YAML::Node windows = analysis["time_windows"];
        if (!windows.IsMap())
          throw std::runtime_error("analysis.time_windows must be a map");
        config.hannTimeWindows = windows["enabled"] ?
          windows["enabled"].as<bool>() : false;
        if (config.hannTimeWindows)
          {
            const YAML::Node interval = required(windows, "interval_s");
            if (!interval.IsSequence() || interval.size() != 2)
              throw std::runtime_error(
                "time_windows.interval_s must be [start, end]");
            config.intervalStart = finiteDouble(interval[0],
              "time-window interval start");
            config.intervalEnd = finiteDouble(interval[1],
              "time-window interval end");
            config.windowDuration = finiteDouble(
              required(windows, "duration_s"), "time-window duration");
            config.windowStep = finiteDouble(required(windows, "step_s"),
              "time-window step");
            if (!(config.intervalEnd > config.intervalStart) ||
                !(config.windowDuration > 0.0) ||
                !(config.windowStep > 0.0) ||
                config.windowDuration >
                  config.intervalEnd - config.intervalStart)
              throw std::runtime_error("Invalid Hann time-window geometry");
          }
      }

    config.frequencyBlock = calculation["frequency_block"] ?
      positiveSize(calculation["frequency_block"], "frequency_block") : 8;
    config.spatialBatchPoints = calculation["spatial_batch_points"] ?
      positiveSize(calculation["spatial_batch_points"],
        "spatial_batch_points") : 32;
    config.fftThreads = calculation["fft_threads"] ?
      positiveSize(calculation["fft_threads"], "fft_threads") : 1;
    config.maximumWorkingMiB = calculation["maximum_working_mib"] ?
      positiveSize(calculation["maximum_working_mib"],
        "maximum_working_mib") : 4096;
    config.maximumOutputGiB = calculation["maximum_output_gib"] ?
      finiteDouble(calculation["maximum_output_gib"],
        "maximum_output_gib") : 64.0;
    if (!(config.maximumOutputGiB > 0.0))
      throw std::runtime_error("maximum_output_gib must be positive");

    if (root["coherence"])
      {
        const YAML::Node coherence = root["coherence"];
        if (!coherence.IsMap())
          throw std::runtime_error("coherence must be a map");
        config.spatialPhotonEnergyEV = positiveEnergyList(
          coherence["spatial_photon_energy_eV"],
          "coherence.spatial_photon_energy_eV");
        config.spatialReferenceAngles = referenceAngles(
          coherence["spatial_reference_angles_rad"],
          "coherence.spatial_reference_angles_rad");
        config.temporalPhotonEnergyEV = positiveEnergyList(
          coherence["temporal_photon_energy_eV"],
          "coherence.temporal_photon_energy_eV");
        config.temporalReferenceAngles = referenceAngles(
          coherence["temporal_reference_angles_rad"],
          "coherence.temporal_reference_angles_rad");
        if (config.spatialPhotonEnergyEV.empty() !=
            config.spatialReferenceAngles.empty())
          throw std::runtime_error(
            "Spatial coherence needs both energies and reference angles");
        if (config.temporalPhotonEnergyEV.empty() !=
            config.temporalReferenceAngles.empty())
          throw std::runtime_error(
            "Temporal coherence needs both energies and reference angles");
      }

    config.outputFile = resolvePath(filename,
      required(output, "file").as<std::string>());
    config.compression = output["compression"] ?
      output["compression"].as<unsigned int>() : 0;
    if (config.compression > 9)
      throw std::runtime_error("output.compression must be in [0,9]");
    if (config.outputFile == config.fieldFile ||
        (!config.baselineFile.empty() &&
         config.outputFile == config.baselineFile))
      throw std::runtime_error("Output must not overwrite an input field");
    return config;
  }
}
