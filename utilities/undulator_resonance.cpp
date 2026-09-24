#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

namespace
{
  const double c = 299792458.0;
  const double planckEVSecond = 4.135667696e-15;

  std::string lower(std::string value)
  {
    std::transform(value.begin(), value.end(), value.begin(),
      [](unsigned char character) {
        return static_cast<char>(std::tolower(character));
      });
    return value;
  }

  YAML::Node required(const YAML::Node& parent, const char* key)
  {
    const YAML::Node value = parent[key];
    if (!value)
      throw std::runtime_error(std::string("Missing required YAML key '") +
        key + "'");
    return value;
  }

  double finite(const YAML::Node& node, const char* description)
  {
    const double value = node.as<double>();
    if (!std::isfinite(value))
      throw std::runtime_error(std::string(description) + " must be finite");
    return value;
  }

  double positive(const YAML::Node& node, const char* description)
  {
    const double value = finite(node, description);
    if (!(value > 0.0))
      throw std::runtime_error(std::string(description) + " must be positive");
    return value;
  }

  double lengthScale(const YAML::Node& units)
  {
    if (!units) return 1.0;
    const std::string unit = lower(
      required(units, "length").as<std::string>());
    if (unit == "m" || unit == "meter" || unit == "metre") return 1.0;
    if (unit == "mm" || unit == "millimeter" || unit == "millimetre")
      return 1.0e-3;
    if (unit == "um" || unit == "micrometer" || unit == "micrometre")
      return 1.0e-6;
    if (unit == "nm" || unit == "nanometer" || unit == "nanometre")
      return 1.0e-9;
    throw std::runtime_error("Unsupported input length unit '" + unit + "'");
  }

  double commandLineGamma(const char* text)
  {
    errno = 0;
    char* end = NULL;
    const double value = std::strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' ||
        !(value > 1.0) || !std::isfinite(value))
      throw std::runtime_error(
        "LAB_GAMMA must be one finite number greater than one");
    return value;
  }

  double cardGamma(const YAML::Node& root)
  {
    const YAML::Node input = root["beam"]["input"];
    if (!input || !input["gamma"])
      throw std::runtime_error(
        "The card has no beam.input.gamma; pass LAB_GAMMA as the second argument for an HDF5 beam");
    const double gamma = positive(input["gamma"], "beam gamma");
    if (!(gamma > 1.0))
      throw std::runtime_error("beam gamma must be greater than one");
    return gamma;
  }

  double delayFactor(double gamma)
  {
    const double inverseGamma2 = 1.0 / (gamma * gamma);
    const double beta = std::sqrt(1.0 - inverseGamma2);
    return inverseGamma2 / (beta * (1.0 + beta));
  }

  YAML::Node characteristicUndulator(const YAML::Node& root)
  {
    const YAML::Node magnets = root["sources"]["magnetic_elements"];
    if (!magnets || !magnets.IsSequence())
      throw std::runtime_error(
        "sources.magnetic_elements must be a nonempty sequence");
    std::vector<YAML::Node> undulators;
    std::vector<YAML::Node> marked;
    for (std::size_t index = 0; index < magnets.size(); ++index)
      {
        const YAML::Node magnet = magnets[index];
        if (lower(required(magnet, "type").as<std::string>()) !=
            "planar-undulator")
          continue;
        undulators.push_back(magnet);
        if (magnet["characteristic"] &&
            magnet["characteristic"].as<bool>())
          marked.push_back(magnet);
      }
    if (marked.size() == 1) return marked.front();
    if (marked.size() > 1)
      throw std::runtime_error(
        "More than one planar undulator is marked characteristic: true");
    if (undulators.size() == 1) return undulators.front();
    if (undulators.empty())
      throw std::runtime_error("The card contains no planar undulator");
    throw std::runtime_error(
      "The card contains multiple planar undulators; mark exactly one characteristic: true");
  }

  void printProtection(const YAML::Node& root, double gamma,
                       double inputLengthScale)
  {
    const YAML::Node planes = root["detectors"]["field_planes"];
    if (!planes || !planes.IsSequence()) return;
    const YAML::Node retirement = root["particle_retirement"];
    for (std::size_t index = 0; index < planes.size(); ++index)
      {
        const YAML::Node plane = planes[index];
        const YAML::Node protection =
          plane["retirement_frequency_protection"];
        if (!protection ||
            (protection["enabled"] &&
             !protection["enabled"].as<bool>()))
          continue;
        const double minimumEnergy = positive(
          required(protection, "minimum_photon_energy_eV"),
          "minimum protected photon energy");
        const unsigned int cycles =
          required(protection, "cycles").as<unsigned int>();
        if (cycles == 0)
          throw std::runtime_error("Protection cycles must be positive");
        const double requiredLength = static_cast<double>(cycles) *
          planckEVSecond * c / (minimumEnergy * delayFactor(gamma));
        std::cout << "\nRetirement-frequency guard for field detector '"
                  << required(plane, "name").as<std::string>() << "'\n"
                  << "  protected band begins [eV]: " << minimumEnergy << "\n"
                  << "  requested cycles N:          " << cycles << "\n"
                  << "  minimum retirement length:  " << requiredLength
                  << " m (" << requiredLength / inputLengthScale
                  << " input length units)\n";
        if (retirement && retirement["enabled"] &&
            retirement["enabled"].as<bool>() && retirement["length"])
          {
            const double actualLength = positive(
              retirement["length"], "particle retirement length") *
              inputLengthScale;
            const double actualCycles = actualLength * delayFactor(gamma) *
              minimumEnergy / (c * planckEVSecond);
            std::cout << "  configured retirement length: " << actualLength
                      << " m\n"
                      << "  configured observed cycles:   " << actualCycles
                      << "\n"
                      << "  preflight length result:      "
                      << (actualLength >= requiredLength ? "PASS" : "FAIL")
                      << "\n";
          }
      }
  }
}

int main(int argc, char** argv)
{
  if (argc < 2 || argc > 3)
    {
      std::cerr << "Usage: " << argv[0]
                << " INPUT.yaml [LAB_GAMMA]\n";
      return 2;
    }
  try
    {
      const YAML::Node root = YAML::LoadFile(argv[1]);
      if (!root.IsMap())
        throw std::runtime_error("Simulation card root must be a map");
      const double inputLengthScale = lengthScale(root["units"]);
      const double gamma = argc == 3 ?
        commandLineGamma(argv[2]) : cardGamma(root);
      const YAML::Node undulator = characteristicUndulator(root);
      const double strength = finite(
        required(undulator, "strength_parameter"), "undulator K");
      const double period = positive(
        required(undulator, "period"), "undulator period") *
        inputLengthScale;
      const double denominator = 1.0 + 0.5 * strength * strength;
      const double wavelength = period * denominator /
        (2.0 * gamma * gamma);
      const double frequency = c / wavelength;
      const double photonEnergy = planckEVSecond * frequency;
      std::cout << std::setprecision(12)
                << "Characteristic planar-undulator resonance\n"
                << "  lab gamma:                    " << gamma << "\n"
                << "  K:                            " << strength << "\n"
                << "  period [m]:                   " << period << "\n"
                << "  on-axis fundamental lambda:   " << wavelength << " m\n"
                << "  on-axis fundamental frequency:" << " " << frequency
                << " Hz\n"
                << "  photon energy:                " << photonEnergy
                << " eV\n"
                << "  interpretation: 1D cold-beam gain-centre estimate; "
                   "energy spread, emittance, space charge and 3D gain "
                   "shifts are not included.\n";
      printProtection(root, gamma, inputLengthScale);
      return 0;
    }
  catch (const std::exception& error)
    {
      std::cerr << "Undulator resonance error: " << error.what() << "\n";
      return 1;
    }
}
