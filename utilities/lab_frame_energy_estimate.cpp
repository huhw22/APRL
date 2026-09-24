#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "yaml-cpp/yaml.h"

namespace
{
  const long double kC = 299792458.0L;
  const long double kElectronCharge = 1.602176634e-19L;
  const long double kElectronMass = 9.1093837139e-31L;
  const long double kEpsilon0 = 8.8541878128e-12L;
  const long double kPlanck = 6.62607015e-34L;
  const long double kPi = 3.141592653589793238462643383279502884L;
  const long double kEulerGamma = 0.577215664901532860606512090082402431L;

  YAML::Node required(const YAML::Node& parent, const char* key)
  {
    const YAML::Node value = parent[key];
    if (!value)
      throw std::runtime_error(std::string("Missing required YAML key '") +
        key + "'");
    return value;
  }

  long double positive(const YAML::Node& parent, const char* key)
  {
    const long double value = required(parent, key).as<long double>();
    if (!(value > 0.0L) || !std::isfinite(value))
      throw std::runtime_error(std::string(key) + " must be finite and > 0");
    return value;
  }

  long double nonnegative(const YAML::Node& parent, const char* key,
                          long double fallback)
  {
    if (!parent[key]) return fallback;
    const long double value = parent[key].as<long double>();
    if (!(value >= 0.0L) || !std::isfinite(value))
      throw std::runtime_error(std::string(key) + " must be finite and >= 0");
    return value;
  }

  std::string parentDirectory(const std::string& path)
  {
    const std::string::size_type separator = path.find_last_of('/');
    if (separator == std::string::npos) return std::string();
    if (separator == 0) return "/";
    return path.substr(0, separator);
  }

  std::string resolvePath(const std::string& card, const std::string& path)
  {
    if (path.empty() || path[0] == '/') return path;
    const std::string parent = parentDirectory(card);
    return parent.empty() ? path : parent + "/" + path;
  }

  bool directoryExists(const std::string& path)
  {
    struct stat status;
    return ::stat(path.c_str(), &status) == 0 && S_ISDIR(status.st_mode);
  }

  void createDirectories(const std::string& directory)
  {
    if (directory.empty() || directory == "." || directory == "/") return;
    std::string current;
    std::size_t offset = 0;
    if (directory[0] == '/')
      {
        current = "/";
        offset = 1;
      }
    while (offset <= directory.size())
      {
        const std::string::size_type separator = directory.find('/', offset);
        const std::string part = directory.substr(offset,
          separator == std::string::npos ? std::string::npos :
          separator - offset);
        if (!part.empty())
          {
            if (!current.empty() && current[current.size() - 1] != '/')
              current += "/";
            current += part;
            if (!directoryExists(current) &&
                ::mkdir(current.c_str(), 0775) != 0 && errno != EEXIST)
              throw std::runtime_error("Cannot create directory: " + current +
                ": " + std::strerror(errno));
          }
        if (separator == std::string::npos) break;
        offset = separator + 1;
      }
  }

  struct Inputs
  {
    long double current;
    long double gamma;
    long double sigmaX;
    long double sigmaY;
    long double period;
    long double undulatorK;
    long double periods;
    long double outerRadius;
    long double coherentEnhancement;
    long double duration;
    std::vector<long double> relativeSigmaChanges;
    std::string report;
  };

  Inputs load(const std::string& filename)
  {
    const YAML::Node root = YAML::LoadFile(filename);
    const YAML::Node beam = required(root, "beam");
    const YAML::Node undulator = required(root, "undulator");
    const YAML::Node estimate = required(root, "estimate");
    if (!root.IsMap() || !beam.IsMap() || !undulator.IsMap() ||
        !estimate.IsMap())
      throw std::runtime_error(
        "beam, undulator, and estimate must be YAML maps");

    Inputs inputs;
    inputs.current = positive(beam, "average_current_A");
    inputs.gamma = positive(beam, "mean_gamma");
    if (inputs.gamma < 1.0L)
      throw std::runtime_error("mean_gamma must be >= 1");
    inputs.sigmaX = positive(beam, "sigma_x_m");
    inputs.sigmaY = positive(beam, "sigma_y_m");
    inputs.duration = nonnegative(beam, "effective_flat_top_duration_s", 0.0L);
    inputs.period = positive(undulator, "period_m");
    inputs.undulatorK = nonnegative(undulator, "strength_parameter", 0.0L);
    inputs.periods = positive(undulator, "periods");
    inputs.outerRadius = positive(estimate, "outer_radius_m");
    inputs.coherentEnhancement = positive(
      estimate, "manual_coherent_enhancement");

    const YAML::Node changes = required(estimate, "relative_sigma_changes");
    if (!changes.IsSequence() || changes.size() == 0)
      throw std::runtime_error(
        "relative_sigma_changes must be a non-empty YAML sequence");
    for (std::size_t i = 0; i < changes.size(); ++i)
      {
        const long double change = changes[i].as<long double>();
        if (!std::isfinite(change) || !(change > -1.0L))
          throw std::runtime_error(
            "every relative_sigma_change must be finite and > -1");
        inputs.relativeSigmaChanges.push_back(change);
      }
    inputs.report = root["output"] && root["output"]["report"] ?
      resolvePath(filename, root["output"]["report"].as<std::string>()) :
      std::string();
    return inputs;
  }

  void writeReport(std::ostream& out, const Inputs& in)
  {
    const long double gammaBeta = std::sqrt(in.gamma - 1.0L) *
      std::sqrt(in.gamma + 1.0L);
    const long double beta = gammaBeta / in.gamma;
    const long double sigma = std::sqrt(in.sigmaX * in.sigmaY);
    const long double aspect = std::max(in.sigmaX, in.sigmaY) /
      std::min(in.sigmaX, in.sigmaY);
    if (in.outerRadius < 5.0L * sigma)
      throw std::runtime_error(
        "outer_radius_m must be at least five equivalent rms radii for the "
        "long-round-Gaussian estimate");

    const long double lineCharge = in.current / (beta * kC);
    const long double electronsPerMetre = lineCharge / kElectronCharge;
    const long double electronRestEnergy = kElectronMass * kC * kC;
    const long double classicalRadius = kElectronCharge * kElectronCharge /
      (4.0L * kPi * kEpsilon0 * electronRestEnergy);
    const long double alfvenCurrent = 4.0L * kPi * kEpsilon0 *
      kElectronMass * kC * kC * kC / kElectronCharge;
    const long double perveance = 2.0L * in.current /
      (alfvenCurrent * beta * beta * beta * in.gamma * in.gamma * in.gamma);

    const long double kineticPerBunchMetre = electronsPerMetre *
      (in.gamma - 1.0L) * electronRestEnergy;
    const long double electrostaticScale = lineCharge * lineCharge /
      (4.0L * kPi * kEpsilon0);
    const long double shape = std::log(in.outerRadius / (2.0L * sigma)) +
      0.5L * kEulerGamma;
    const long double boundPerBunchMetre =
      (1.0L + beta * beta) * electrostaticScale * shape;
    const long double boundScale =
      (1.0L + beta * beta) * electrostaticScale;

    const long double ku = 2.0L * kPi / in.period;
    const long double spontaneousLossPerElectronPerMetre =
      classicalRadius * electronRestEnergy * in.gamma * in.gamma *
      in.undulatorK * in.undulatorK * ku * ku / 3.0L;
    const long double spontaneousLossPerElectronPerPeriod =
      spontaneousLossPerElectronPerMetre * in.period;
    const long double radiationPerBunchMetrePerPeriod = electronsPerMetre *
      spontaneousLossPerElectronPerPeriod * in.coherentEnhancement;
    const long double radiationPerBunchMetre =
      radiationPerBunchMetrePerPeriod * in.periods;
    const long double spontaneousPowerPerUndulatorMetre =
      (in.current / kElectronCharge) *
      spontaneousLossPerElectronPerMetre;
    const long double configuredPower = spontaneousPowerPerUndulatorMetre *
      in.period * in.periods * in.coherentEnhancement;

    const long double resonantWavelength = in.period *
      (1.0L + 0.5L * in.undulatorK * in.undulatorK) /
      (2.0L * in.gamma * in.gamma);
    const long double resonantPhotonEnergy =
      kPlanck * kC / (resonantWavelength * kElectronCharge);

    out << std::setprecision(18) << std::scientific;
    out << "format_version: 1\n"
        << "observer_frame: laboratory\n"
        << "external_static_undulator_work_J: 0\n"
        << "inputs:\n"
        << "  average_current_A: " << in.current << "\n"
        << "  mean_gamma: " << in.gamma << "\n"
        << "  sigma_x_m: " << in.sigmaX << "\n"
        << "  sigma_y_m: " << in.sigmaY << "\n"
        << "  undulator_period_m: " << in.period << "\n"
        << "  undulator_K: " << in.undulatorK << "\n"
        << "  configured_periods: " << in.periods << "\n"
        << "  bound_field_outer_radius_m: " << in.outerRadius << "\n"
        << "  manual_coherent_enhancement: " << in.coherentEnhancement << "\n"
        << "derived_beam:\n"
        << "  beta: " << beta << "\n"
        << "  line_charge_C_per_m: " << lineCharge << "\n"
        << "  represented_electrons_per_bunch_m: " << electronsPerMetre << "\n"
        << "  kinetic_energy_J_per_bunch_m: " << kineticPerBunchMetre << "\n"
        << "  generalized_perveance: " << perveance << "\n"
        << "  equivalent_round_sigma_m: " << sigma << "\n"
        << "  transverse_aspect_ratio: " << aspect << "\n"
        << "undulator_reference:\n"
        << "  fundamental_wavelength_m: " << resonantWavelength << "\n"
        << "  fundamental_photon_energy_eV: " << resonantPhotonEnergy << "\n"
        << "  spontaneous_loss_eV_per_electron_per_period: "
        << spontaneousLossPerElectronPerPeriod / kElectronCharge << "\n"
        << "  spontaneous_loss_eV_per_electron_per_undulator_m: "
        << spontaneousLossPerElectronPerMetre / kElectronCharge << "\n"
        << "  spontaneous_power_W_per_undulator_m: "
        << spontaneousPowerPerUndulatorMetre << "\n"
        << "  configured_radiation_power_W: " << configuredPower << "\n"
        << "lab_bound_field_reference:\n"
        << "  stored_E_plus_B_energy_J_per_bunch_m: "
        << boundPerBunchMetre << "\n"
        << "  stored_bound_to_particle_kinetic_ratio: "
        << boundPerBunchMetre / kineticPerBunchMetre << "\n"
        << "radiation_reference:\n"
        << "  energy_J_per_bunch_m_per_period: "
        << radiationPerBunchMetrePerPeriod << "\n"
        << "  configured_energy_J_per_bunch_m: "
        << radiationPerBunchMetre << "\n"
        << "  bound_store_to_configured_radiation_ratio: "
        << boundPerBunchMetre / radiationPerBunchMetre << "\n"
        << "transverse_size_change_sweep:\n";
    for (std::size_t i = 0; i < in.relativeSigmaChanges.size(); ++i)
      {
        const long double change = in.relativeSigmaChanges[i];
        const long double deltaBound = -boundScale * std::log1p(change);
        out << "  - relative_sigma_change: " << change << "\n"
            << "    delta_bound_E_plus_B_J_per_bunch_m: " << deltaBound << "\n"
            << "    abs_delta_to_radiation_per_period: "
            << std::abs(deltaBound) / radiationPerBunchMetrePerPeriod << "\n"
            << "    abs_delta_to_configured_radiation: "
            << std::abs(deltaBound) / radiationPerBunchMetre << "\n";
      }

    if (in.duration > 0.0L)
      {
        const long double bunchLength = beta * kC * in.duration;
        out << "optional_finite_bunch:\n"
            << "  effective_flat_top_duration_s: " << in.duration << "\n"
            << "  bunch_length_m: " << bunchLength << "\n"
            << "  total_charge_C: " << in.current * in.duration << "\n"
            << "  particle_kinetic_energy_J: "
            << kineticPerBunchMetre * bunchLength << "\n"
            << "  approximate_bound_E_plus_B_energy_J: "
            << boundPerBunchMetre * bunchLength << "\n"
            << "  configured_radiation_energy_J: "
            << radiationPerBunchMetre * bunchLength << "\n";
      }

    out << "interpretation:\n"
        << "  - stored bound-field energy is not its change; only the change enters the lab energy balance\n"
        << "  - positive relative_sigma_change is expansion and gives negative delta bound-field energy\n"
        << "  - spontaneous radiation is an incoherent reference unless manual_coherent_enhancement is supplied\n"
        << "  - a pre-bunched coherent prediction needs a measured bunching/form factor and cannot be inferred from current and rms size alone\n"
        << "  - the long round-Gaussian bound-field estimate uses an equivalent rms radius and a specified outer cutoff\n"
        << "  - the transverse-size sweep holds current fixed; longitudinal compression or current evolution must be evaluated separately because the bound-field scale is proportional to current squared\n";
  }
}

int main(int argc, char** argv)
{
  if (argc != 2)
    {
      std::cerr << "Usage: lab_frame_energy_estimate <estimate.yaml>\n";
      return 2;
    }
  try
    {
      const Inputs inputs = load(argv[1]);
      writeReport(std::cout, inputs);
      if (!inputs.report.empty())
        {
          createDirectories(parentDirectory(inputs.report));
          std::ofstream report(inputs.report.c_str());
          if (!report)
            throw std::runtime_error("Cannot create report: " + inputs.report);
          writeReport(report, inputs);
          report.close();
          if (!report)
            throw std::runtime_error("Cannot finish report: " + inputs.report);
        }
      return 0;
    }
  catch (const std::exception& error)
    {
      std::cerr << "lab_frame_energy_estimate: " << error.what() << "\n";
      return 1;
    }
}
