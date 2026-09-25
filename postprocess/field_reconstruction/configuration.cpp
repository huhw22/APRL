#include "configuration.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "yaml-cpp/yaml.h"
#include "../yaml_validation.h"

namespace reconstruction
{
  namespace
  {
    std::runtime_error cardError(const YAML::Node& node,
                                 const std::string& message)
    {
      const YAML::Mark mark = node.Mark();
      return std::runtime_error(message + " at YAML line " +
        std::to_string(mark.line + 1));
    }

    YAML::Node required(const YAML::Node& parent, const char* key)
    {
      const YAML::Node value = parent[key];
      if (!value) throw cardError(parent,
        std::string("Missing required key '") + key + "'");
      return value;
    }

    std::size_t positiveSize(const YAML::Node& node, const char* name)
    {
      const unsigned long long value = node.as<unsigned long long>();
      if (value == 0 || value >
          static_cast<unsigned long long>(
            std::numeric_limits<std::size_t>::max()))
        throw cardError(node, std::string(name) + " must be positive");
      return static_cast<std::size_t>(value);
    }

    double finiteNonnegative(const YAML::Node& node, const char* name)
    {
      const double value = node.as<double>();
      if (!(value >= 0.0) || !std::isfinite(value))
        throw cardError(node,
          std::string(name) + " must be finite and nonnegative");
      return value;
    }

    std::string parentDirectory(const std::string& path)
    {
      const std::string::size_type separator = path.find_last_of('/');
      if (separator == std::string::npos) return std::string();
      if (separator == 0) return "/";
      return path.substr(0, separator);
    }

    bool absolutePath(const std::string& path)
    {
      return !path.empty() && path[0] == '/';
    }

    std::string resolvePath(const std::string& card,
                            const std::string& path)
    {
      if (absolutePath(path)) return path;
      const std::string parent = parentDirectory(card);
      return parent.empty() ? path : parent + "/" + path;
    }
  }

  ReconstructionConfig::ReconstructionConfig()
    : cardPath(), fieldFile(), particleFile(), requireComplete(true),
      particleReadChunk(65536), fftThreads(1), paddingTime(1),
      paddingY(1), paddingX(1), transverseSmoothing(0.0),
      maximumRelativeGammaSpread(0.1), maximumRmsTransverseBeta(0.1),
      maximumOutsideChargeFraction(1.0e-3), outputFile(), compression(0),
      overwrite(false)
  {}

  ReconstructionConfig loadConfiguration(const std::string& filename)
  {
    const YAML::Node root = YAML::LoadFile(filename);
    if (!root.IsMap())
      throw cardError(root, "Field-reconstruction card must be a map");
    postprocess_common::validateMapKeys(root, "top level",
      {"input", "model", "output"});
    ReconstructionConfig config;
    config.cardPath = filename;

    const YAML::Node input = required(root, "input");
    if (!input.IsMap()) throw cardError(input, "input must be a map");
    postprocess_common::validateMapKeys(input, "input", {
      "field_file", "particle_file", "require_complete",
      "particle_read_chunk"
    });
    config.fieldFile = resolvePath(filename,
      required(input, "field_file").as<std::string>());
    config.particleFile = resolvePath(filename,
      required(input, "particle_file").as<std::string>());
    if (input["require_complete"])
      config.requireComplete = input["require_complete"].as<bool>();
    if (input["particle_read_chunk"])
      config.particleReadChunk = positiveSize(
        input["particle_read_chunk"], "particle_read_chunk");

    const YAML::Node model = required(root, "model");
    if (!model.IsMap()) throw cardError(model, "model must be a map");
    postprocess_common::validateMapKeys(model, "model", {
      "fft_threads", "padding_factor", "transverse_smoothing_m",
      "maximum_relative_gamma_spread", "maximum_rms_transverse_beta",
      "maximum_outside_charge_fraction"
    });
    if (model["fft_threads"])
      config.fftThreads = static_cast<unsigned int>(positiveSize(
        model["fft_threads"], "fft_threads"));
    const YAML::Node padding = model["padding_factor"];
    if (padding)
      {
        if (!padding.IsSequence() || padding.size() != 3)
          throw cardError(padding,
            "padding_factor must be [time, y, x]");
        config.paddingTime = positiveSize(padding[0], "time padding factor");
        config.paddingY = positiveSize(padding[1], "y padding factor");
        config.paddingX = positiveSize(padding[2], "x padding factor");
      }
    if (model["transverse_smoothing_m"])
      config.transverseSmoothing = finiteNonnegative(
        model["transverse_smoothing_m"], "transverse_smoothing_m");
    if (model["maximum_relative_gamma_spread"])
      config.maximumRelativeGammaSpread = finiteNonnegative(
        model["maximum_relative_gamma_spread"],
        "maximum_relative_gamma_spread");
    if (model["maximum_rms_transverse_beta"])
      config.maximumRmsTransverseBeta = finiteNonnegative(
        model["maximum_rms_transverse_beta"],
        "maximum_rms_transverse_beta");
    if (model["maximum_outside_charge_fraction"])
      {
        config.maximumOutsideChargeFraction = finiteNonnegative(
          model["maximum_outside_charge_fraction"],
          "maximum_outside_charge_fraction");
        if (config.maximumOutsideChargeFraction > 1.0)
          throw cardError(model["maximum_outside_charge_fraction"],
            "maximum_outside_charge_fraction must not exceed one");
      }

    const YAML::Node output = required(root, "output");
    if (!output.IsMap()) throw cardError(output, "output must be a map");
    postprocess_common::validateMapKeys(output, "output",
      {"file", "compression", "overwrite"});
    config.outputFile = resolvePath(filename,
      required(output, "file").as<std::string>());
    if (output["compression"])
      config.compression = output["compression"].as<unsigned int>();
    if (output["overwrite"])
      config.overwrite = output["overwrite"].as<bool>();
    if (config.compression > 9)
      throw cardError(output, "output compression must be in [0,9]");
    if (config.outputFile == config.fieldFile ||
        config.outputFile == config.particleFile)
      throw cardError(output,
        "output file must differ from both input files");
    postprocess_common::requireOutputAvailable(
      config.outputFile, config.overwrite);
    return config;
  }
}
