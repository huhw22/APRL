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

  double betaFromGamma(double gamma)
  {
    if (!(gamma >= 1.0) || !std::isfinite(gamma))
      throw std::runtime_error("gamma must be finite and at least one");
    return std::sqrt(gamma - 1.0) * std::sqrt(gamma + 1.0) / gamma;
  }

  double configuredBoostGamma(const YAML::Node& root)
  {
    const YAML::Node mesh = required(root, "mesh");
    const double gamma = finite(required(mesh, "boost_gamma"),
      "configured boost gamma");
    if (!(gamma >= 1.0))
      throw std::runtime_error(
        "mesh.boost_gamma must be at least one");
    return gamma;
  }

  double meshCellSizeZ(const YAML::Node& root, double inputLengthScale)
  {
    const YAML::Node cells = required(required(root, "mesh"), "cell_size");
    if (!cells.IsSequence() || cells.size() != 3)
      throw std::runtime_error(
        "mesh.cell_size must contain exactly three values");
    return positive(cells[2], "mesh z cell size") * inputLengthScale;
  }

  unsigned long long meshCellsZ(const YAML::Node& root)
  {
    const YAML::Node cells = required(required(root, "mesh"), "cells");
    if (!cells.IsSequence() || cells.size() != 3)
      throw std::runtime_error(
        "mesh.cells must contain exactly three values");
    const unsigned long long value = cells[2].as<unsigned long long>();
    if (value == 0)
      throw std::runtime_error("mesh.cells z must be positive");
    return value;
  }

  unsigned int undulatorPeriods(const YAML::Node& undulator)
  {
    const unsigned int periods =
      required(undulator, "periods").as<unsigned int>();
    if (periods == 0)
      throw std::runtime_error("undulator periods must be positive");
    return periods;
  }

  void printBoostRecommendation(const YAML::Node& root,
                                double gammaLab,
                                double strength,
                                double period,
                                const YAML::Node& undulator,
                                double inputLengthScale)
  {
    const double transverseFactor = 1.0 + 0.5 * strength * strength;
    const double gammaLongitudinal =
      gammaLab / std::sqrt(transverseFactor);
    if (!(gammaLongitudinal > 1.0) ||
        !std::isfinite(gammaLongitudinal))
      throw std::runtime_error(
        "The planar-undulator longitudinal-gamma estimate is not "
        "relativistic; the paraxial recommendation is inapplicable");

    const double betaLab = betaFromGamma(gammaLab);
    const double betaLongitudinal = betaFromGamma(gammaLongitudinal);
    const double gammaBoost = configuredBoostGamma(root);
    const double betaBoost = betaFromGamma(gammaBoost);
    const unsigned int periods = undulatorPeriods(undulator);
    const double length = period * static_cast<double>(periods);
    const double dz = meshCellSizeZ(root, inputLengthScale);
    const double extentZ = dz * static_cast<double>(meshCellsZ(root));

    /* For a particle crossing a laboratory length L at mean beta_z,
     * Lorentz transformation gives Delta z' = gamma_b L
     * (1-beta_b/beta_z).  The recommendation makes this mean drift zero in
     * the constant-K undulator core. */
    const double configuredDrift = gammaBoost * length *
      (1.0 - betaBoost / betaLongitudinal);
    const double configuredDriftPerPeriod = configuredDrift /
      static_cast<double>(periods);
    const double recommendedFreeDriftPerLabMetre = gammaLongitudinal *
      (1.0 - betaLongitudinal / betaLab);

    std::cout << "\nConstant-boost recommendation for the characteristic "
                 "undulator\n"
              << "  interpretation: static B changes direction, not total "
                 "lab gamma; <u_perp^2>=K^2/2 lowers mean beta_z\n"
              << "  planar longitudinal factor:   "
              << std::sqrt(transverseFactor) << "\n"
              << "  estimated mean gamma_z:       "
              << gammaLongitudinal << "\n"
              << "  estimated mean beta_z:        "
              << betaLongitudinal << "\n"
              << "  recommended mesh boost_gamma: "
              << gammaLongitudinal << "\n"
              << "  configured mesh boost_gamma:  "
              << gammaBoost << "\n"
              << "  configured/recommended ratio: "
              << gammaBoost / gammaLongitudinal << "\n"
              << "  physical core length:         " << length << " m ("
              << periods << " periods)\n"
              << "  configured box-z drift/core:  "
              << configuredDrift << " m\n"
              << "  configured drift per period:  "
              << configuredDriftPerPeriod << " m\n"
              << "  configured drift in z cells:  "
              << configuredDrift / dz << "\n"
              << "  |drift|/box z extent:         "
              << std::abs(configuredDrift) / extentZ << "\n"
              << "  recommended-frame free-drift box shift per lab metre: "
              << recommendedFreeDriftPerLabMetre << " m/m\n"
              << "  suggested YAML: boost_gamma:  "
              << gammaLongitudinal << "\n"
              << "  scope: fixed inertial boost optimized for mean motion "
                 "inside the constant-K core; fringes, free drifts, energy "
                 "spread, emittance and collective energy change still "
                 "require longitudinal-box margin.\n";
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
    const YAML::Node detectors = root["detectors"];
    if (!detectors) return;
    const YAML::Node planes = detectors["field_planes"];
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
      printBoostRecommendation(root, gamma, strength, period, undulator,
                               inputLengthScale);
      printProtection(root, gamma, inputLengthScale);
      return 0;
    }
  catch (const std::exception& error)
    {
      std::cerr << "Undulator resonance error: " << error.what() << "\n";
      return 1;
    }
}
