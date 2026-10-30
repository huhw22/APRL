#include "config.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

namespace aprl
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

    void validateMapKeys(const YAML::Node& node,
                         const std::string& context,
                         std::initializer_list<const char*> allowedKeys)
    {
      if (!node.IsMap())
        throw configError(node, context + " must be a map");

      std::set<std::string> allowed;
      std::ostringstream choices;
      bool first = true;
      for (std::initializer_list<const char*>::const_iterator key =
             allowedKeys.begin(); key != allowedKeys.end(); ++key)
        {
          allowed.insert(*key);
          if (!first) choices << ", ";
          choices << *key;
          first = false;
        }

      for (YAML::const_iterator entry = node.begin(); entry != node.end();
           ++entry)
        {
          std::string key;
          try { key = entry->first.as<std::string>(); }
          catch (const YAML::Exception&) {
            throw configError(entry->first,
              context + " keys must be scalar strings");
          }
          if (allowed.find(key) == allowed.end())
            throw configError(entry->first,
              context + " contains unknown key '" + key +
              "'; recognized keys are: " + choices.str());
        }
    }

    void validateConfigurationKeys(const YAML::Node& root)
    {
      validateMapKeys(root, "top level", {
        "units", "runtime", "mesh", "radiation_resolution", "boundary",
        "initial_self_field", "beam", "sources", "detectors",
        "particle_retirement", "stop", "trajectory", "energy_ledger",
        "output", "checkpoint"
      });

      const YAML::Node output = root["output"];
      if (output)
        validateMapKeys(output, "output", {"overwrite", "manifest"});

      const YAML::Node checkpoint = root["checkpoint"];
      if (checkpoint)
        validateMapKeys(checkpoint, "checkpoint", {
          "directory", "steps", "times", "restart_from",
          "allow_retirement_change", "reset_diagnostics", "io_timeout_seconds"
        });

      const YAML::Node units = root["units"];
      if (units)
        validateMapKeys(units, "units", {"length", "time"});

      const YAML::Node runtime = root["runtime"];
      if (runtime)
        {
          validateMapKeys(runtime, "runtime", {
            "mode", "stop_check_interval_steps", "resource_monitor",
            "shutdown_watchdog"
          });
          const YAML::Node resources = runtime["resource_monitor"];
          if (resources)
            validateMapKeys(resources, "runtime.resource_monitor", {
              "enabled", "progress_interval_steps", "calibration_steps",
              "memory_safety_factor", "time_safety_factor"
            });
          const YAML::Node watchdog = runtime["shutdown_watchdog"];
          if (watchdog)
            validateMapKeys(watchdog, "runtime.shutdown_watchdog", {
              "enabled", "timeout_seconds", "finalize_timeout_seconds"
            });
        }

      const YAML::Node mesh = root["mesh"];
      if (mesh)
        {
          validateMapKeys(mesh, "mesh", {
            "field_solver", "current_filter", "cells", "cell_size", "center",
            "duration", "boost_gamma", "particle_steps_per_undulator_period",
            "maximum_particle_substeps", "longitudinal", "lengths",
            "resolution", "ballistic_outside_transverse",
            "self_field_feedback", "moving_window"
          });
          const YAML::Node longitudinal = mesh["longitudinal"];
          if (longitudinal)
            validateMapKeys(longitudinal, "mesh.longitudinal", {
              "mode", "lower_guard_cells", "upper_guard_cells",
              "self_field_padding_rms", "ballistic_time_safety_factor",
              "alignment_cells"
            });
          const YAML::Node moving = mesh["moving_window"];
          if (moving)
            validateMapKeys(moving, "mesh.moving_window", {
              "enabled", "cells", "lower_guard_cells", "upper_guard_cells",
              "block_layers", "reference_directory", "reference_through_cpml", "reference_surface_current", "inlet_charge_match"
            });
        }

      const YAML::Node radiation = root["radiation_resolution"];
      if (radiation)
        validateMapKeys(radiation, "radiation_resolution", {
          "enabled", "maximum_photon_energy_eV",
          "warning_grid_points_per_wavelength",
          "warning_maxwell_samples_per_cycle",
          "warning_detector_samples_per_cycle"
        });

      const YAML::Node boundary = root["boundary"];
      if (boundary)
        validateMapKeys(boundary, "boundary", {
          "type", "cells", "polynomial_order", "target_reflection",
          "kappa_max", "alpha_fraction"
        });

      const YAML::Node initial = root["initial_self_field"];
      if (initial)
        validateMapKeys(initial, "initial_self_field", {
          "enabled", "model", "relative_tolerance", "maximum_iterations",
          "maximum_mean_transverse_beta", "include_transverse_cpml_particles"
        });

      const YAML::Node beam = root["beam"];
      if (beam)
        {
          validateMapKeys(beam, "beam", {"reference", "input"});
          const YAML::Node reference = beam["reference"];
          if (reference)
            validateMapKeys(reference, "beam.reference", {
              "initial_center_z", "input_plane_z",
              "frozen_lab_plane_injection", "activation_plane_z",
              "hold_transverse_until_activation"
            });
          const YAML::Node input = beam["input"];
          if (input)
            validateMapKeys(input, "beam.input", {
              "type", "electrons", "position_offset", "file",
              "macroparticles", "gamma", "direction", "center",
              "sigma_position", "sigma_proper_velocity", "random_seed"
            });
        }

      const YAML::Node sources = root["sources"];
      if (sources)
        {
          validateMapKeys(sources, "sources", {
            "incident_waves", "magnetic_elements"
          });
          const YAML::Node waves = sources["incident_waves"];
          if (waves)
            {
              if (!waves.IsSequence())
                throw configError(waves,
                  "sources.incident_waves must be a sequence");
              for (std::size_t i = 0; i < waves.size(); ++i)
                {
                  std::ostringstream context;
                  context << "sources.incident_waves[" << i << "]";
                  validateMapKeys(waves[i], context.str(), {
                    "profile", "position", "direction", "polarization",
                    "wavelength", "peak_electric_field_V_per_m",
                    "normalized_amplitude", "radius", "order", "envelope"
                  });
                  const YAML::Node envelope = waves[i]["envelope"];
                  if (envelope)
                    validateMapKeys(envelope, context.str() + ".envelope", {
                      "type", "center_time", "duration",
                      "carrier_phase_rad", "rising_cycles",
                      "inverse_gaussian_sigma"
                    });
                }
            }
          const YAML::Node magnets = sources["magnetic_elements"];
          if (magnets)
            {
              if (!magnets.IsSequence())
                throw configError(magnets,
                  "sources.magnetic_elements must be a sequence");
              for (std::size_t i = 0; i < magnets.size(); ++i)
                {
                  std::ostringstream context;
                  context << "sources.magnetic_elements[" << i << "]";
                  validateMapKeys(magnets[i], context.str(), {
                    "type", "entrance_z", "polarization_angle_rad",
                    "strength_parameter", "period", "periods",
                    "gaussian_fringe", "fringe_relative_cutoff", "length",
                    "field_T", "characteristic"
                  });
                }
            }
        }

      const YAML::Node detectors = root["detectors"];
      if (detectors)
        {
          validateMapKeys(detectors, "detectors", {
            "enabled", "directory", "field_planes", "particle_planes", "energy_planes"
          });
          const YAML::Node fieldPlanes = detectors["field_planes"];
          if (fieldPlanes)
            {
              if (!fieldPlanes.IsSequence())
                throw configError(fieldPlanes,
                  "detectors.field_planes must be a sequence");
              for (std::size_t i = 0; i < fieldPlanes.size(); ++i)
                {
                  std::ostringstream context;
                  context << "detectors.field_planes[" << i << "]";
                  validateMapKeys(fieldPlanes[i], context.str(), {
                    "name", "z", "rhythm", "buffer_samples", "compression",
                    "particle_background", "retirement_frequency_protection"
                  });
                  const YAML::Node background =
                    fieldPlanes[i]["particle_background"];
                  if (background)
                    {
                      validateMapKeys(background,
                        context.str() + ".particle_background", {
                          "enabled", "buffer_records", "compression",
                          "validation"
                        });
                      const YAML::Node validation = background["validation"];
                      if (validation)
                        validateMapKeys(validation,
                          context.str() + ".particle_background.validation", {
                            "enabled", "maximum_particles"
                          });
                    }
                  const YAML::Node protection =
                    fieldPlanes[i]["retirement_frequency_protection"];
                  if (protection)
                    validateMapKeys(protection,
                      context.str() + ".retirement_frequency_protection", {
                        "enabled", "minimum_photon_energy_eV", "cycles"
                      });
                }
            }
          const YAML::Node energyPlanes = detectors["energy_planes"];
          if (energyPlanes)
            {
              if (!energyPlanes.IsSequence())
                throw configError(energyPlanes, "detectors.energy_planes must be a sequence");
              for (std::size_t i = 0; i < energyPlanes.size(); ++i)
                validateMapKeys(energyPlanes[i], "detectors.energy_planes", {
                  "name", "z", "rhythm", "buffer_samples", "compression"
                });
            }
          const YAML::Node particlePlanes = detectors["particle_planes"];
          if (particlePlanes)
            {
              if (!particlePlanes.IsSequence())
                throw configError(particlePlanes,
                  "detectors.particle_planes must be a sequence");
              for (std::size_t i = 0; i < particlePlanes.size(); ++i)
                {
                  std::ostringstream context;
                  context << "detectors.particle_planes[" << i << "]";
                  validateMapKeys(particlePlanes[i], context.str(), {
                    "name", "z", "buffer_records", "compression"
                  });
                }
            }
        }

      const YAML::Node retirement = root["particle_retirement"];
      if (retirement)
        validateMapKeys(retirement, "particle_retirement", {
          "enabled", "entrance_z", "length"
        });

      const YAML::Node stop = root["stop"];
      if (stop)
        validateMapKeys(stop, "stop", {"mode", "z"});

      const YAML::Node trajectory = root["trajectory"];
      if (trajectory)
        validateMapKeys(trajectory, "trajectory", {
          "enabled", "directory", "basename", "rhythm", "mode",
          "buffer_records", "flush_every_samples", "compression"
        });

      const YAML::Node ledger = root["energy_ledger"];
      if (ledger)
        validateMapKeys(ledger, "energy_ledger", {
          "enabled", "directory", "filename", "sample_interval_steps",
          "buffer_records", "compression", "warning_relative_tolerance",
          "discrete_operator_audit"
        });
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

    std::size_t movingWindowSize(const YAML::Node& node,
                                const std::string& name)
    {
      // Do not let an unsigned YAML conversion turn a negative value into a
      // huge positive allocation. Leave the existing config paths unchanged.
      if (!node.IsScalar())
        throw configError(node, name + " must be a positive integer");
      const std::string text = node.Scalar();
      const std::string::size_type first = text.find_first_not_of(" \t\r\n");
      if (first == std::string::npos || text[first] == '-')
        throw configError(node, name + " must be a positive integer");
      return positiveSize(node, name);
    }

    std::size_t checkpointSize(const YAML::Node& node,
                               const std::string& name)
    {
      std::size_t value = 0;
      try
        {
          if (!node.IsScalar() || node.Scalar().empty() ||
              node.Scalar()[0] == '-')
            throw configError(node, name + " must be a nonnegative integer");
          value = node.as<std::size_t>();
        }
      catch (const YAML::Exception&)
        {
          throw configError(node, name +
            " must be a nonnegative integer in the supported size range");
        }
      return value;
    }

    std::string checkpointPath(const YAML::Node& node,
                               const std::string& name)
    {
      if (!node.IsScalar())
        throw configError(node, name + " must be a nonempty path");
      const std::string value = node.as<std::string>();
      if (value.find_first_not_of(" \t\r\n") == std::string::npos ||
          value.find('\0') != std::string::npos)
        throw configError(node, name + " must be a nonempty path without NUL bytes");
      return value;
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

    bool meshCells(const YAML::Node& node, std::size_t result[3])
    {
      if (!node.IsSequence() || node.size() != 3)
        throw configError(node,
          "mesh cells must contain exactly three entries");
      for (unsigned int axis = 0; axis < 3; ++axis)
        {
          if (axis == 2)
            {
              std::string token;
              try { token = lower(node[axis].as<std::string>()); }
              catch (const YAML::Exception&) {}
              if (token == "auto")
                {
                  result[axis] = 0;
                  return true;
                }
            }
          result[axis] = positiveSize(node[axis], "mesh cells");
          if (result[axis] < 3)
            throw configError(node[axis],
              "each mesh direction needs at least three cells");
          if (static_cast<long double>(result[axis]) >
              std::ldexp(1.0L, std::numeric_limits<Double>::digits))
            throw configError(node[axis],
              "mesh cell count exceeds the exactly representable double "
              "index range");
        }
      return false;
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

    EBCurrentFilter currentFilter(const YAML::Node& node)
    {
      const std::string value = lower(node.as<std::string>());
      if (value == "none") return EBCurrentFilter::None;
      if (value == "compensated-bilinear")
        return EBCurrentFilter::CompensatedBilinear;
      throw configError(node,
        "current_filter must be none or compensated-bilinear");
    }

    EBBoundaryType boundaryType(const YAML::Node& node)
    {
      const std::string value = lower(node.as<std::string>());
      if (value == "pec") return EBBoundaryType::Pec;
      if (value == "cpml") return EBBoundaryType::Cpml;
      throw configError(node, "boundary type must be pec or cpml");
    }

    InitialSelfFieldModel initialSelfFieldModel(const YAML::Node& node)
    {
      const std::string value = lower(node.as<std::string>());
      if (value == "relativistic-poisson" || value == "vay")
        return InitialSelfFieldModel::RelativisticPoisson;
      if (value == "electrostatic-poisson" || value == "legacy")
        return InitialSelfFieldModel::ElectrostaticPoisson;
      if (value == "cowan-relativistic-poisson")
        return InitialSelfFieldModel::CowanRelativisticPoisson;
      throw configError(node,
        "initial_self_field model must be relativistic-poisson, "
        "electrostatic-poisson or cowan-relativistic-poisson");
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

    std::string readTextFile(const std::string& filename)
    {
      std::ifstream stream(filename.c_str(), std::ios::in | std::ios::binary);
      if (!stream)
        throw std::runtime_error("Cannot reopen YAML configuration: " +
          filename);
      std::ostringstream contents;
      contents << stream.rdbuf();
      if (!stream.good() && !stream.eof())
        throw std::runtime_error("Cannot read YAML configuration bytes: " +
          filename);
      return contents.str();
    }

    std::string fnv1a64(const std::string& value)
    {
      std::uint64_t hash = UINT64_C(14695981039346656037);
      for (std::size_t i = 0; i < value.size(); ++i)
        {
          hash ^= static_cast<unsigned char>(value[i]);
          hash *= UINT64_C(1099511628211);
        }
      std::ostringstream result;
      result << std::hex << std::setfill('0') << std::setw(16) << hash;
      return result.str();
    }
  }

  UnitSystem::UnitSystem() : length(1.0), time(1.0) {}

  LongitudinalMeshConfig::LongitudinalMeshConfig()
    : mode(LongitudinalMeshMode::Explicit), lowerGuardCells(4),
      upperGuardCells(8), selfFieldPaddingRms(4.0),
      ballisticTimeSafetyFactor(1.1), alignmentCells(1)
  {}

  bool LongitudinalMeshConfig::automatic() const
  {
    return mode != LongitudinalMeshMode::Explicit;
  }

  MovingWindowConfig::MovingWindowConfig()
    : enabled(false), cells(0), lowerGuardCells(64), upperGuardCells(128),
      blockLayers(32), referenceThroughCpml(false), referenceSurfaceCurrent(false), inletChargeMatch(false), referenceDirectory()
  {}

  MeshConfig::MeshConfig()
    : cellSize(0.0), extent(0.0), center(0.0), duration(0.0),
      boostGamma(1.0), particleStepsPerUndulatorPeriod(1024),
      maximumParticleSubsteps(4096),
      fieldSolver(EBMaxwellSolver::CowanZ),
      currentFilter(EBCurrentFilter::None),
      selfFieldFeedback(true),
      ballisticOutsideTransverse(false)
  {
    cells[0] = cells[1] = cells[2] = 0;
  }

  RadiationResolutionConfig::RadiationResolutionConfig()
    : enabled(false), maximumPhotonEnergyEV(0.0),
      warningGridPointsPerWavelength(8.0),
      warningMaxwellSamplesPerCycle(8.0),
      warningDetectorSamplesPerCycle(8.0)
  {}

  BoundaryConfig::BoundaryConfig()
    : type(EBBoundaryType::Pec), cpml()
  {}

  InitialSelfFieldConfig::InitialSelfFieldConfig()
    : enabled(true), includeTransverseCpmlParticles(false),
      model(InitialSelfFieldModel::RelativisticPoisson),
      relativeTolerance(1.0e-10), maximumMeanTransverseBeta(1.0e-8),
      maximumIterations(10000)
  {}

  BeamInputConfig::BeamInputConfig()
    : type(BeamInputType::Hdf5), file(), electrons(0.0), positionOffset(0.0),
      macroparticles(0), gamma(1.0), direction(0.0), center(0.0),
      sigmaPosition(0.0), sigmaProperVelocity(0.0), randomSeed(1)
  {
    direction[2] = 1.0;
  }

  BeamReferenceConfig::BeamReferenceConfig()
    : initialCenterZ(0.0), inputPlaneZ(0.0), inputPlaneZSet(false),
      frozenLabPlaneInjection(false),
      holdTransverseUntilActivation(false), activationPlaneZ(0.0),
      activationPlaneZSet(false)
  {}

  BeamlineElementExtent::BeamlineElementExtent()
    : role(BeamlineElementRole::MagneticDevice), physicalEntrance(0.0),
      physicalExit(0.0), interactionEntrance(0.0), interactionExit(0.0)
  {}

  StopConfig::StopConfig()
    : mode(StopMode::AfterLastElement), referenceZ(0.0)
  {}

  RuntimeConfig::RuntimeConfig()
    : mode(RuntimeMode::Throughput), stopCheckIntervalSteps(16),
      resourceReport(true), resourceProgressIntervalSteps(1000),
      resourceCalibrationSteps(1), memorySafetyFactor(1.25),
      timeSafetyFactor(1.25), shutdownWatchdogEnabled(true),
      shutdownWatchdogTimeoutSeconds(600),
      shutdownWatchdogFinalizeTimeoutSeconds(60)
  {}

  OutputPolicyConfig::OutputPolicyConfig()
    : overwrite(false), manifest("run-manifest.yaml")
  {}

  CheckpointConfig::CheckpointConfig()
    : directory("checkpoints"), steps(), times(), restartFrom(),
      allowRetirementChange(false), resetDiagnostics(false),
      ioTimeoutSeconds(3600), restartStep(0), restartTimeBox(0.0)
  {}

  bool CheckpointConfig::enabled() const
  {
    return !steps.empty() || !times.empty();
  }

  bool RuntimeConfig::interactive() const
  {
    return mode == RuntimeMode::Interactive;
  }

  TrajectoryConfig::TrajectoryConfig()
    : enabled(false), directory("./"), basename("trajectory"), rhythm(0.0),
      bufferRecords(16384), flushEverySamples(8), compression(0)
  {}

  EnergyLedgerConfig::EnergyLedgerConfig()
    : enabled(false), discreteOperatorAudit(false), directory("output"),
      filename("energy-ledger.h5"), sampleIntervalSteps(1),
      bufferRecords(64), compression(0), warningRelativeTolerance(0.01)
  {}

  ParticleRetirementConfig::ParticleRetirementConfig()
    : enabled(false), entranceZ(0.0), length(0.0)
  {}

  Double ParticleRetirementConfig::exitZ() const
  {
    return entranceZ + length;
  }

  FieldDetectorPlaneConfig::FieldDetectorPlaneConfig()
    : name(), z(0.0), rhythm(0.0), bufferSamples(1), compression(0),
      particleBackgroundReference(false), referenceRho(0.0),
      referenceGamma(1.0), referenceDistance(0.0), referenceEntranceZ(0.0),
      referenceBufferRecords(16384), referenceCompression(0),
      referenceValidation(false),
      referenceValidationMaximumParticles(100000),
      particleRetirementCurrent(false), particleRetirementEntranceZ(0.0),
      particleRetirementExitZ(0.0),
      retirementFrequencyProtection(false),
      retirementMinimumPhotonEnergyEV(0.0), retirementMinimumCycles(0),
      retirementCausalDistance(0.0), retirementCausalEntranceZ(0.0),
      retirementRequiredLength(0.0), retirementObservedCycles(0.0)
  {}

  ParticleDetectorPlaneConfig::ParticleDetectorPlaneConfig()
    : name(), z(0.0), bufferRecords(16384), compression(0)
  {}

  EnergyDetectorPlaneConfig::EnergyDetectorPlaneConfig()
    : name(), z(0.0), rhythm(0.0), bufferSamples(128), compression(0)
  {}

  DetectorConfig::DetectorConfig()
    : directory("output/detectors"), fieldPlanes(), particlePlanes(), energyPlanes()
  { physicalGuardCells[0] = physicalGuardCells[1] = 0; }

  bool DetectorConfig::enabled() const
  {
    return !fieldPlanes.empty() || !particlePlanes.empty() || !energyPlanes.empty();
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
    validateConfigurationKeys(root);

    SimulationConfig result;
    result.configurationPath = filename;
    result.configurationText = readTextFile(filename);
    result.configurationDigest = fnv1a64(result.configurationText);

    const YAML::Node outputPolicy = root["output"];
    if (outputPolicy)
      {
        if (outputPolicy["overwrite"])
          result.output.overwrite = outputPolicy["overwrite"].as<bool>();
        if (outputPolicy["manifest"])
          result.output.manifest =
            outputPolicy["manifest"].as<std::string>();
        if (result.output.manifest.empty())
          throw configError(outputPolicy["manifest"],
            "output manifest path cannot be empty");
      }
    const YAML::Node units = root["units"];
    if (units)
      {
        result.inputUnits.length = lengthScale(required(units, "length"));
        result.inputUnits.time = timeScale(required(units, "time"));
      }

    const YAML::Node checkpoint = root["checkpoint"];
    if (checkpoint)
      {
        if (checkpoint["directory"])
          result.checkpoint.directory = checkpointPath(
            checkpoint["directory"], "checkpoint directory");
        if (checkpoint["restart_from"])
          result.checkpoint.restartFrom = checkpointPath(
            checkpoint["restart_from"], "checkpoint restart_from");
        if (checkpoint["allow_retirement_change"])
          result.checkpoint.allowRetirementChange =
            checkpoint["allow_retirement_change"].as<bool>();
        if (checkpoint["reset_diagnostics"])
          result.checkpoint.resetDiagnostics =
            checkpoint["reset_diagnostics"].as<bool>();
        if (checkpoint["io_timeout_seconds"])
          {
            const std::size_t timeout = checkpointSize(
              checkpoint["io_timeout_seconds"],
              "checkpoint io_timeout_seconds");
            if (timeout == 0 ||
                timeout > std::numeric_limits<unsigned int>::max())
              throw configError(checkpoint["io_timeout_seconds"],
                "checkpoint io_timeout_seconds must be a positive unsigned integer");
            result.checkpoint.ioTimeoutSeconds =
              static_cast<unsigned int>(timeout);
          }
        const YAML::Node steps = checkpoint["steps"];
        if (steps)
          {
            if (!steps.IsSequence())
              throw configError(steps, "checkpoint steps must be a sequence");
            for (std::size_t i = 0; i < steps.size(); ++i)
              result.checkpoint.steps.push_back(checkpointSize(
                steps[i], "checkpoint step"));
            std::sort(result.checkpoint.steps.begin(),
                      result.checkpoint.steps.end());
            result.checkpoint.steps.erase(std::unique(
              result.checkpoint.steps.begin(), result.checkpoint.steps.end()),
              result.checkpoint.steps.end());
          }
        const YAML::Node times = checkpoint["times"];
        if (times)
          {
            if (!times.IsSequence())
              throw configError(times, "checkpoint times must be a sequence");
            for (std::size_t i = 0; i < times.size(); ++i)
              {
                const Double time = finiteDouble(times[i],
                  "checkpoint time") * result.inputUnits.time;
                if (!(time >= 0.0) || !std::isfinite(time))
                  throw configError(times[i],
                    "checkpoint time must be finite and nonnegative after SI conversion");
                result.checkpoint.times.push_back(time);
              }
            std::sort(result.checkpoint.times.begin(),
                      result.checkpoint.times.end());
            result.checkpoint.times.erase(std::unique(
              result.checkpoint.times.begin(), result.checkpoint.times.end()),
              result.checkpoint.times.end());
          }
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
        const YAML::Node resources = runtime["resource_monitor"];
        if (resources)
          {
            if (!resources.IsMap())
              throw configError(resources,
                "runtime.resource_monitor must be a map");
            if (resources["enabled"])
              result.runtime.resourceReport =
                resources["enabled"].as<bool>();
            if (resources["progress_interval_steps"])
              result.runtime.resourceProgressIntervalSteps =
                resources["progress_interval_steps"].as<std::size_t>();
            if (resources["calibration_steps"])
              result.runtime.resourceCalibrationSteps =
                resources["calibration_steps"].as<unsigned int>();
            if (resources["memory_safety_factor"])
              result.runtime.memorySafetyFactor = finiteDouble(
                resources["memory_safety_factor"],
                "resource monitor memory safety factor");
            if (resources["time_safety_factor"])
              result.runtime.timeSafetyFactor = finiteDouble(
                resources["time_safety_factor"],
                "resource monitor time safety factor");
            if (result.runtime.resourceReport &&
                result.runtime.resourceCalibrationSteps == 0)
              throw configError(resources["calibration_steps"],
                "enabled resource monitor needs at least one calibration step");
            if (!(result.runtime.memorySafetyFactor >= 1.0) ||
                !(result.runtime.timeSafetyFactor >= 1.0))
              throw configError(resources,
                "resource-monitor safety factors must be at least one");
          }
        const YAML::Node watchdog = runtime["shutdown_watchdog"];
        if (watchdog)
          {
            if (!watchdog.IsMap())
              throw configError(watchdog,
                "runtime.shutdown_watchdog must be a map");
            if (watchdog["enabled"])
              result.runtime.shutdownWatchdogEnabled =
                watchdog["enabled"].as<bool>();
            if (watchdog["timeout_seconds"])
              {
                const std::size_t timeout = positiveSize(
                  watchdog["timeout_seconds"],
                  "runtime shutdown_watchdog timeout_seconds");
                if (timeout > std::numeric_limits<unsigned int>::max())
                  throw configError(watchdog["timeout_seconds"],
                    "runtime shutdown_watchdog timeout_seconds exceeds "
                    "the supported unsigned-integer range");
                result.runtime.shutdownWatchdogTimeoutSeconds =
                  static_cast<unsigned int>(timeout);
              }
            if (watchdog["finalize_timeout_seconds"])
              {
                const std::size_t timeout = positiveSize(
                  watchdog["finalize_timeout_seconds"],
                  "runtime shutdown_watchdog finalize_timeout_seconds");
                if (timeout > std::numeric_limits<unsigned int>::max())
                  throw configError(watchdog["finalize_timeout_seconds"],
                    "runtime shutdown_watchdog finalize_timeout_seconds "
                    "exceeds the supported unsigned-integer range");
                result.runtime.shutdownWatchdogFinalizeTimeoutSeconds =
                  static_cast<unsigned int>(timeout);
              }
          }
      }

    const YAML::Node mesh = required(root, "mesh");
    if (mesh["lengths"] || mesh["resolution"])
      throw configError(mesh,
        "mesh.lengths and mesh.resolution are obsolete; specify integer "
        "mesh.cells and physical mesh.cell_size instead");
    const bool automaticLongitudinal =
      meshCells(required(mesh, "cells"), result.mesh.cells);
    if (mesh["ballistic_outside_transverse"])
      result.mesh.ballisticOutsideTransverse =
        mesh["ballistic_outside_transverse"].as<bool>();
    if (mesh["self_field_feedback"])
      result.mesh.selfFieldFeedback =
        mesh["self_field_feedback"].as<bool>();
    const YAML::Node movingWindow = mesh["moving_window"];
    if (movingWindow)
      {
        if (movingWindow["enabled"])
          result.mesh.movingWindow.enabled = movingWindow["enabled"].as<bool>();
        if (movingWindow["cells"])
          result.mesh.movingWindow.cells = movingWindowSize(
            movingWindow["cells"], "mesh.moving_window cells");
        if (movingWindow["lower_guard_cells"])
          result.mesh.movingWindow.lowerGuardCells = movingWindowSize(
            movingWindow["lower_guard_cells"],
            "mesh.moving_window lower_guard_cells");
        if (movingWindow["upper_guard_cells"])
          result.mesh.movingWindow.upperGuardCells = movingWindowSize(
            movingWindow["upper_guard_cells"],
            "mesh.moving_window upper_guard_cells");
        if (movingWindow["block_layers"])
          result.mesh.movingWindow.blockLayers = movingWindowSize(
            movingWindow["block_layers"], "mesh.moving_window block_layers");
        if (movingWindow["reference_through_cpml"])
          result.mesh.movingWindow.referenceThroughCpml =
            movingWindow["reference_through_cpml"].as<bool>();
        if (movingWindow["reference_surface_current"])
          result.mesh.movingWindow.referenceSurfaceCurrent =
            movingWindow["reference_surface_current"].as<bool>();
        if(movingWindow["inlet_charge_match"])
          result.mesh.movingWindow.inletChargeMatch=movingWindow["inlet_charge_match"].as<bool>();
        if (movingWindow["reference_directory"])
          {
            if (!movingWindow["reference_directory"].IsScalar())
              throw configError(movingWindow["reference_directory"],
                "mesh.moving_window reference_directory must be a path");
            result.mesh.movingWindow.referenceDirectory =
              movingWindow["reference_directory"].as<std::string>();
            if (result.mesh.movingWindow.referenceDirectory.find('\0') !=
                std::string::npos)
              throw configError(movingWindow["reference_directory"],
                "mesh.moving_window reference_directory cannot contain NUL bytes");
          }
        if (result.mesh.movingWindow.enabled &&
            result.mesh.movingWindow.referenceDirectory.find_first_not_of(
              " \t\r\n") == std::string::npos)
          throw configError(movingWindow,
            "enabled mesh.moving_window requires a nonempty reference_directory");
        if (result.mesh.movingWindow.referenceThroughCpml && !result.mesh.movingWindow.enabled)
          throw configError(movingWindow,"reference_through_cpml requires enabled moving_window");
        if (result.mesh.movingWindow.referenceSurfaceCurrent && (!result.mesh.movingWindow.enabled || result.mesh.movingWindow.referenceThroughCpml))
          throw configError(movingWindow,"reference_surface_current requires moving_window and the original inner reference boundary");
        if(result.mesh.movingWindow.inletChargeMatch && !result.mesh.movingWindow.enabled)
          throw configError(movingWindow,"inlet_charge_match requires moving_window");
      }
    const YAML::Node longitudinal = mesh["longitudinal"];
    if (automaticLongitudinal)
      {
        if (!longitudinal)
          throw configError(mesh["cells"][2],
            "mesh.cells z is auto, so mesh.longitudinal is required");
        const std::string mode = lower(
          required(longitudinal, "mode").as<std::string>());
        if (mode != "observation-bounded")
          throw configError(longitudinal["mode"],
            "automatic longitudinal mode must be observation-bounded");
        result.mesh.longitudinal.mode =
          LongitudinalMeshMode::ObservationBounded;
        if (longitudinal["lower_guard_cells"])
          result.mesh.longitudinal.lowerGuardCells = positiveSize(
            longitudinal["lower_guard_cells"],
            "mesh.longitudinal lower_guard_cells");
        if (longitudinal["upper_guard_cells"])
          result.mesh.longitudinal.upperGuardCells = positiveSize(
            longitudinal["upper_guard_cells"],
            "mesh.longitudinal upper_guard_cells");
        if (longitudinal["self_field_padding_rms"])
          result.mesh.longitudinal.selfFieldPaddingRms = finiteDouble(
            longitudinal["self_field_padding_rms"],
            "mesh.longitudinal self_field_padding_rms");
        if (longitudinal["ballistic_time_safety_factor"])
          result.mesh.longitudinal.ballisticTimeSafetyFactor = finiteDouble(
            longitudinal["ballistic_time_safety_factor"],
            "mesh.longitudinal ballistic_time_safety_factor");
        if (longitudinal["alignment_cells"])
          result.mesh.longitudinal.alignmentCells = positiveSize(
            longitudinal["alignment_cells"],
            "mesh.longitudinal alignment_cells");
        if (result.mesh.longitudinal.lowerGuardCells < 2 ||
            result.mesh.longitudinal.upperGuardCells < 2)
          throw configError(longitudinal,
            "automatic longitudinal guards must each be at least two cells");
        if (!(result.mesh.longitudinal.selfFieldPaddingRms >= 0.0))
          throw configError(longitudinal["self_field_padding_rms"],
            "self_field_padding_rms must be nonnegative");
        if (!(result.mesh.longitudinal.ballisticTimeSafetyFactor >= 1.0))
          throw configError(longitudinal["ballistic_time_safety_factor"],
            "ballistic_time_safety_factor must be at least one");
      }
    else if (longitudinal)
      throw configError(longitudinal,
        "mesh.longitudinal is only valid when mesh.cells z is auto");
    result.mesh.cellSize = finiteVector3(required(mesh, "cell_size"),
      "mesh cell_size", result.inputUnits.length);
    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        if (!(result.mesh.cellSize[axis] > 0.0))
          throw configError(mesh["cell_size"][axis],
            "mesh cell_size values must be positive");
        if (axis == 2 && automaticLongitudinal)
          {
            result.mesh.extent[axis] = 0.0;
            continue;
          }
        const long double length =
          static_cast<long double>(result.mesh.cells[axis]) *
          static_cast<long double>(result.mesh.cellSize[axis]);
        if (!(length > 0.0L) ||
            length > static_cast<long double>(
              std::numeric_limits<Double>::max()))
          throw configError(mesh["cells"][axis],
            "mesh cells * cell_size produces an invalid physical extent");
        result.mesh.extent[axis] = static_cast<Double>(length);
      }
    result.mesh.center = mesh["center"] ?
      finiteVector3(mesh["center"], "mesh center", result.inputUnits.length) :
      FieldVector<Double>(0.0);
    result.mesh.duration = finiteDouble(required(mesh, "duration"),
      "mesh duration") * result.inputUnits.time;
    result.mesh.boostGamma = finiteDouble(required(mesh, "boost_gamma"),
      "boost gamma");
    result.mesh.fieldSolver = fieldSolver(required(mesh, "field_solver"));
    result.mesh.currentFilter = mesh["current_filter"] ?
      currentFilter(mesh["current_filter"]) : EBCurrentFilter::None;
    if(result.mesh.movingWindow.inletChargeMatch&&result.mesh.currentFilter!=EBCurrentFilter::None)
      throw configError(mesh["current_filter"],"inlet_charge_match requires unfiltered current");
    if (result.mesh.currentFilter != EBCurrentFilter::None &&
        result.mesh.fieldSolver != EBMaxwellSolver::CowanZ)
      throw configError(mesh["current_filter"],
        "current_filter is an experimental Cowan-z source treatment; use "
        "field_solver: cowan-z or set current_filter: none");
    const std::size_t particleSteps = positiveSize(
      required(mesh, "particle_steps_per_undulator_period"),
      "particle_steps_per_undulator_period");
    if (particleSteps > std::numeric_limits<unsigned int>::max())
      throw configError(mesh["particle_steps_per_undulator_period"],
        "particle_steps_per_undulator_period exceeds the unsigned-int range");
    result.mesh.particleStepsPerUndulatorPeriod =
      static_cast<unsigned int>(particleSteps);
    if (mesh["maximum_particle_substeps"])
      {
        const std::size_t maximumSubsteps = positiveSize(
          mesh["maximum_particle_substeps"],
          "maximum_particle_substeps");
        if (maximumSubsteps > std::numeric_limits<unsigned int>::max())
          throw configError(mesh["maximum_particle_substeps"],
            "maximum_particle_substeps exceeds the unsigned-int range");
        result.mesh.maximumParticleSubsteps =
          static_cast<unsigned int>(maximumSubsteps);
      }

    const YAML::Node radiationResolution = root["radiation_resolution"];
    if (radiationResolution)
      {
        if (!radiationResolution.IsMap())
          throw configError(radiationResolution,
            "radiation_resolution must be a map");
        result.radiationResolution.enabled =
          radiationResolution["enabled"] ?
          radiationResolution["enabled"].as<bool>() : true;
        if (result.radiationResolution.enabled)
          {
            result.radiationResolution.maximumPhotonEnergyEV = finiteDouble(
              required(radiationResolution, "maximum_photon_energy_eV"),
              "radiation-resolution maximum photon energy");
            if (!(result.radiationResolution.maximumPhotonEnergyEV > 0.0))
              throw configError(
                radiationResolution["maximum_photon_energy_eV"],
                "maximum_photon_energy_eV must be positive");
          }
        if (radiationResolution["warning_grid_points_per_wavelength"])
          result.radiationResolution.warningGridPointsPerWavelength =
            finiteDouble(
              radiationResolution["warning_grid_points_per_wavelength"],
              "radiation-resolution grid warning level");
        if (radiationResolution["warning_detector_samples_per_cycle"])
          result.radiationResolution.warningDetectorSamplesPerCycle =
            finiteDouble(
              radiationResolution["warning_detector_samples_per_cycle"],
              "radiation-resolution detector warning level");
        if (radiationResolution["warning_maxwell_samples_per_cycle"])
          result.radiationResolution.warningMaxwellSamplesPerCycle =
            finiteDouble(
              radiationResolution["warning_maxwell_samples_per_cycle"],
              "radiation-resolution Maxwell warning level");
        if (!(result.radiationResolution.warningGridPointsPerWavelength >=
              2.0) ||
            !(result.radiationResolution.warningMaxwellSamplesPerCycle >=
              2.0) ||
            !(result.radiationResolution.warningDetectorSamplesPerCycle >=
              2.0))
          throw configError(radiationResolution,
            "radiation-resolution warning levels must be at least the "
            "two-sample Nyquist limit");
      }

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

    const YAML::Node initialSelfField = root["initial_self_field"];
    if (initialSelfField)
      {
        if (!initialSelfField.IsMap())
          throw configError(initialSelfField,
            "initial_self_field must be a map");
        if (initialSelfField["enabled"])
          result.initialSelfField.enabled =
            initialSelfField["enabled"].as<bool>();
        if (initialSelfField["include_transverse_cpml_particles"])
          result.initialSelfField.includeTransverseCpmlParticles =
            initialSelfField["include_transverse_cpml_particles"].as<bool>();
        if (initialSelfField["model"])
          result.initialSelfField.model = initialSelfFieldModel(
            initialSelfField["model"]);
        if (initialSelfField["relative_tolerance"])
          result.initialSelfField.relativeTolerance = finiteDouble(
            initialSelfField["relative_tolerance"],
            "initial self-field relative_tolerance");
        if (initialSelfField["maximum_mean_transverse_beta"])
          result.initialSelfField.maximumMeanTransverseBeta = finiteDouble(
            initialSelfField["maximum_mean_transverse_beta"],
            "initial self-field maximum_mean_transverse_beta");
        if (initialSelfField["maximum_iterations"])
          result.initialSelfField.maximumIterations = positiveSize(
            initialSelfField["maximum_iterations"],
            "initial self-field maximum_iterations");
        if (!(result.initialSelfField.relativeTolerance > 0.0) ||
            !(result.initialSelfField.relativeTolerance < 1.0))
          throw configError(initialSelfField["relative_tolerance"],
            "initial self-field relative_tolerance must lie between zero and one");
        if (!(result.initialSelfField.maximumMeanTransverseBeta >= 0.0) ||
            !(result.initialSelfField.maximumMeanTransverseBeta < 1.0))
          throw configError(initialSelfField["maximum_mean_transverse_beta"],
            "initial self-field maximum_mean_transverse_beta must lie in [0, 1)");
      }

    const YAML::Node beam = required(root, "beam");
    const YAML::Node reference = required(beam, "reference");
    result.reference.initialCenterZ = finiteDouble(
      required(reference, "initial_center_z"),
      "initial beam-center z") * result.inputUnits.length;
    if (reference["input_plane_z"])
      {
        result.reference.inputPlaneZ = finiteDouble(
          reference["input_plane_z"], "particle input-plane z") *
          result.inputUnits.length;
        result.reference.inputPlaneZSet = true;
      }
    if (reference["frozen_lab_plane_injection"])
      result.reference.frozenLabPlaneInjection =
        reference["frozen_lab_plane_injection"].as<bool>();
    if (reference["hold_transverse_until_activation"])
      result.reference.holdTransverseUntilActivation =
        reference["hold_transverse_until_activation"].as<bool>();
    if (result.reference.holdTransverseUntilActivation &&
        !result.reference.frozenLabPlaneInjection)
      throw configError(reference["hold_transverse_until_activation"],
        "hold_transverse_until_activation requires frozen_lab_plane_injection: true");
    if (reference["activation_plane_z"])
      {
        result.reference.activationPlaneZ = finiteDouble(
          reference["activation_plane_z"],
          "particle activation-plane z") * result.inputUnits.length;
        result.reference.activationPlaneZSet = true;
      }
    if (result.reference.frozenLabPlaneInjection &&
        !result.reference.inputPlaneZSet)
      throw configError(reference["frozen_lab_plane_injection"],
        "frozen_lab_plane_injection requires beam.reference.input_plane_z");
    if (result.reference.activationPlaneZSet &&
        !result.reference.frozenLabPlaneInjection)
      throw configError(reference["activation_plane_z"],
        "activation_plane_z requires frozen_lab_plane_injection: true");
    if (result.reference.frozenLabPlaneInjection &&
        !result.reference.activationPlaneZSet)
      result.reference.activationPlaneZ = result.reference.inputPlaneZ;
    if (result.reference.activationPlaneZSet &&
        result.reference.activationPlaneZ > result.reference.inputPlaneZ)
      throw configError(reference["activation_plane_z"],
        "activation_plane_z must be at or upstream of input_plane_z");
    const YAML::Node input = required(beam, "input");
    result.beam.type = beamInputType(required(input, "type"));
    result.beam.electrons = finiteDouble(required(input, "electrons"),
                                         "electron count");
    if (!(result.beam.electrons > 0.0))
      throw configError(input["electrons"], "electron count must be positive");
    if (input["position_offset"])
      result.beam.positionOffset = finiteVector3(input["position_offset"],
        "beam position offset", result.inputUnits.length);
    if (result.reference.frozenLabPlaneInjection &&
        result.beam.positionOffset[2] != 0.0)
      throw configError(input["position_offset"],
        "frozen_lab_plane_injection requires beam.input.position_offset z = 0; shift input_plane_z instead");

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
                const YAML::Node reference = node["particle_background"];
                if (reference)
                  {
                    if (!reference.IsMap())
                      throw configError(reference,
                        "field detector particle_background must be a map");
                    if (reference["enabled"])
                      plane.particleBackgroundReference =
                        reference["enabled"].as<bool>();
                    if (reference["buffer_records"])
                      plane.referenceBufferRecords = positiveSize(
                        reference["buffer_records"],
                        "field detector particle_background buffer_records");
                    if (reference["compression"])
                      plane.referenceCompression =
                        reference["compression"].as<unsigned int>();
                    if (plane.referenceCompression > 9)
                      throw configError(reference,
                        "field detector particle_background compression must be in [0,9]");
                    const YAML::Node validation = reference["validation"];
                    if (validation)
                      {
                        if (!validation.IsMap())
                          throw configError(validation,
                            "field detector particle_background validation must be a map");
                        if (validation["enabled"])
                          plane.referenceValidation =
                            validation["enabled"].as<bool>();
                        if (validation["maximum_particles"])
                          plane.referenceValidationMaximumParticles =
                            positiveSize(validation["maximum_particles"],
                              "field detector particle_background validation maximum_particles");
                      }
                    if (plane.referenceValidation &&
                        !plane.particleBackgroundReference)
                      throw configError(reference,
                        "field detector particle_background validation requires enabled: true");
                  }
                const YAML::Node protection =
                  node["retirement_frequency_protection"];
                if (protection)
                  {
                    if (!protection.IsMap())
                      throw configError(protection,
                        "field detector retirement_frequency_protection must be a map");
                    plane.retirementFrequencyProtection =
                      protection["enabled"] ?
                      protection["enabled"].as<bool>() : true;
                    if (plane.retirementFrequencyProtection)
                      {
                        plane.retirementMinimumPhotonEnergyEV = finiteDouble(
                          required(protection,
                            "minimum_photon_energy_eV"),
                          "retirement protection minimum photon energy");
                        if (!(plane.retirementMinimumPhotonEnergyEV > 0.0))
                          throw configError(
                            protection["minimum_photon_energy_eV"],
                            "retirement protection minimum photon energy must be positive");
                        const std::size_t cycles = positiveSize(
                          required(protection, "cycles"),
                          "retirement protection cycles");
                        if (cycles >
                            std::numeric_limits<unsigned int>::max())
                          throw configError(protection["cycles"],
                            "retirement protection cycles are too large");
                        plane.retirementMinimumCycles =
                          static_cast<unsigned int>(cycles);
                        if (plane.particleBackgroundReference)
                          throw configError(protection,
                            "retirement_frequency_protection and particle_background are alternative field-detector settings");
                      }
                  }
                result.detectors.fieldPlanes.push_back(plane);

                BeamlineElementExtent extent;
                extent.role = BeamlineElementRole::FieldDetectorPlane;
                extent.physicalEntrance = extent.physicalExit = plane.z;
                extent.interactionEntrance = extent.interactionExit = plane.z;
                result.beamlineElements.push_back(extent);
              }
          }

        const YAML::Node energyPlanes = detectors["energy_planes"];
        if (detectorsEnabled && energyPlanes)
          {
            for (std::size_t i = 0; i < energyPlanes.size(); ++i)
              {
                const YAML::Node node = energyPlanes[i];
                EnergyDetectorPlaneConfig plane;
                plane.name = required(node, "name").as<std::string>();
                if (!validDetectorName(plane.name))
                  throw configError(node["name"], "invalid energy detector name");
                if (!detectorNames.insert(plane.name).second)
                  throw configError(node["name"], "detector names must be unique");
                plane.z = finiteDouble(required(node, "z"), "energy detector z") * result.inputUnits.length;
                plane.rhythm = finiteDouble(required(node, "rhythm"), "energy detector rhythm") * result.inputUnits.time;
                if (!(plane.rhythm > 0.0))
                  throw configError(node, "energy detector rhythm must be positive");
                if (node["buffer_samples"])
                  plane.bufferSamples = positiveSize(node["buffer_samples"], "energy detector buffer_samples");
                if (node["compression"])
                  plane.compression = node["compression"].as<unsigned int>();
                if (plane.compression > 9)
                  throw configError(node, "energy detector compression must be in [0,9]");
                result.detectors.energyPlanes.push_back(plane);
                BeamlineElementExtent extent;
                extent.role = BeamlineElementRole::EnergyDetectorPlane;
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

    const YAML::Node retirement = root["particle_retirement"];
    if (retirement)
      {
        if (!retirement.IsMap())
          throw configError(retirement, "particle_retirement must be a map");
        result.particleRetirement.enabled = retirement["enabled"] ?
          retirement["enabled"].as<bool>() : true;
        if (result.particleRetirement.enabled)
          {
            result.particleRetirement.entranceZ = finiteDouble(
              required(retirement, "entrance_z"),
              "particle retirement entrance_z") *
              result.inputUnits.length;
            result.particleRetirement.length = finiteDouble(
              required(retirement, "length"),
              "particle retirement length") * result.inputUnits.length;
            if (!(result.particleRetirement.length > 0.0))
              throw configError(retirement["length"],
                "particle retirement length must be positive");
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

    const YAML::Node energyLedger = root["energy_ledger"];
    if (energyLedger)
      {
        if (!energyLedger.IsMap())
          throw configError(energyLedger, "energy_ledger must be a map");
        result.energyLedger.enabled = energyLedger["enabled"] ?
          energyLedger["enabled"].as<bool>() : true;
        if (energyLedger["discrete_operator_audit"])
          result.energyLedger.discreteOperatorAudit = energyLedger["discrete_operator_audit"].as<bool>();
        if (result.energyLedger.discreteOperatorAudit && !result.energyLedger.enabled)
          throw configError(energyLedger,"discrete_operator_audit requires enabled energy_ledger");
        if (energyLedger["directory"])
          result.energyLedger.directory =
            energyLedger["directory"].as<std::string>();
        if (energyLedger["filename"])
          result.energyLedger.filename =
            energyLedger["filename"].as<std::string>();
        if (energyLedger["sample_interval_steps"])
          result.energyLedger.sampleIntervalSteps = positiveSize(
            energyLedger["sample_interval_steps"],
            "energy_ledger sample_interval_steps");
        if (energyLedger["buffer_records"])
          result.energyLedger.bufferRecords = positiveSize(
            energyLedger["buffer_records"],
            "energy_ledger buffer_records");
        if (energyLedger["compression"])
          result.energyLedger.compression =
            energyLedger["compression"].as<unsigned int>();
        if (energyLedger["warning_relative_tolerance"])
          result.energyLedger.warningRelativeTolerance = finiteDouble(
            energyLedger["warning_relative_tolerance"],
            "energy_ledger warning_relative_tolerance");
        if (result.energyLedger.enabled &&
            (result.energyLedger.directory.empty() ||
             result.energyLedger.filename.empty()))
          throw configError(energyLedger,
            "enabled energy_ledger needs nonempty directory and filename");
        if (result.energyLedger.filename.find('/') != std::string::npos ||
            result.energyLedger.filename.find('\\') != std::string::npos)
          throw configError(energyLedger["filename"],
            "energy_ledger filename must be a basename; use directory for its path");
        if (result.energyLedger.compression > 9)
          throw configError(energyLedger["compression"],
            "energy_ledger compression must be in [0,9]");
        if (!(result.energyLedger.warningRelativeTolerance > 0.0))
          throw configError(energyLedger["warning_relative_tolerance"],
            "energy_ledger warning_relative_tolerance must be positive");
      }

    if (!(result.mesh.boostGamma >= 1.0))
      throw configError(mesh["boost_gamma"],
        "boost_gamma must be at least one");
    if (!(result.mesh.duration > 0.0))
      throw configError(mesh["duration"], "duration must be positive");
    if (result.mesh.movingWindow.enabled)
      {
        if (result.mesh.fieldSolver != EBMaxwellSolver::CowanZ ||
            !result.initialSelfField.enabled ||
            result.initialSelfField.model !=
              InitialSelfFieldModel::CowanRelativisticPoisson)
          throw configError(movingWindow,
            "moving_window requires Cowan-z and enabled "
            "initial_self_field.model: cowan-relativistic-poisson");
        if (result.mesh.currentFilter != EBCurrentFilter::None)
          throw configError(movingWindow,
            "moving_window requires current_filter: none; its causal exit "
            "does not support a longitudinally extended current filter");
        if (result.boundary.type != EBBoundaryType::Cpml ||
            result.boundary.cpml.cells[0] < 2 ||
            result.boundary.cpml.cells[1] < 2 ||
            result.boundary.cpml.cells[2] != 0)
          throw configError(movingWindow,
            "moving_window requires at least two transverse x/y CPML cells "
            "and zero z CPML");
        if (!result.reference.frozenLabPlaneInjection ||
            !result.reference.holdTransverseUntilActivation)
          throw configError(movingWindow,
            "moving_window requires frozen_lab_plane_injection and "
            "hold_transverse_until_activation");
        if (result.checkpoint.enabled() || !result.checkpoint.restartFrom.empty())
          throw configError(movingWindow,
            "moving_window checkpoint writing/restart is not supported");
        if (result.particleRetirement.enabled || !result.waves.empty() ||
            result.mesh.ballisticOutsideTransverse ||
            result.initialSelfField.includeTransverseCpmlParticles)
          throw configError(movingWindow,
            "moving_window does not support particle retirement, incident "
            "waves, ballistic_outside_transverse or transverse CPML source padding");
      }
    if (result.mesh.ballisticOutsideTransverse &&
        result.particleRetirement.enabled)
      throw configError(mesh["ballistic_outside_transverse"],
        "ballistic transverse boundary injection is incompatible with "
        "particle retirement");
    if (result.mesh.ballisticOutsideTransverse &&
        (result.boundary.type != EBBoundaryType::Cpml ||
         result.boundary.cpml.cells[0] < 2 ||
         result.boundary.cpml.cells[1] < 2))
      throw configError(mesh["ballistic_outside_transverse"],
        "ballistic transverse boundary injection requires at least two "
        "CPML cells on both x and y");
    if (result.initialSelfField.model ==
          InitialSelfFieldModel::CowanRelativisticPoisson &&
        result.mesh.fieldSolver != EBMaxwellSolver::CowanZ)
      throw configError(initialSelfField["model"],
        "initial_self_field.model: cowan-relativistic-poisson requires "
        "mesh.field_solver: cowan-z");
    if (result.initialSelfField.includeTransverseCpmlParticles &&
        (!result.initialSelfField.enabled ||
         !result.mesh.ballisticOutsideTransverse ||
         result.boundary.type != EBBoundaryType::Cpml))
      throw configError(initialSelfField["include_transverse_cpml_particles"],
        "initial_self_field.include_transverse_cpml_particles requires "
        "initial_self_field.enabled, mesh.ballistic_outside_transverse "
        "and CPML");
    if (result.mesh.longitudinal.automatic())
      {
        if (!result.mesh.movingWindow.enabled &&
            (result.boundary.type != EBBoundaryType::Cpml ||
             result.boundary.cpml.cells[2] < 2))
          throw configError(mesh["longitudinal"],
            "observation-bounded longitudinal sizing requires z CPML");
        if (!result.waves.empty())
          throw configError(mesh["longitudinal"],
            "observation-bounded longitudinal sizing does not yet support incident waves");
      }
    if (result.stop.mode == StopMode::AfterLastElement &&
        result.beamlineElements.empty())
      throw configError(stop["mode"],
        "after-last-element requires at least one enabled magnetic element "
        "or detector plane; for element-free propagation use "
        "reference-center-z with a finite downstream z");
    if (result.stop.mode == StopMode::ReferenceCenterZ)
      {
        Double requiredBeyond = result.reference.initialCenterZ;
        for (std::size_t i = 0; i < result.beamlineElements.size(); ++i)
          requiredBeyond = std::max(requiredBeyond,
            result.beamlineElements[i].interactionExit);
        if (!(result.stop.referenceZ > requiredBeyond))
          throw configError(stop["z"],
            "reference-center stop z must lie downstream of the initial "
            "beam centre and every element interaction region");
      }
    if (result.particleRetirement.enabled)
      {
        Double lastMagneticExit =
          -std::numeric_limits<Double>::infinity();
        for (std::size_t i = 0; i < result.magnets.size(); ++i)
          lastMagneticExit = std::max(lastMagneticExit,
            result.magnets[i].interactionExitLab());
        const Double margin = result.mesh.boostGamma *
          result.mesh.cellSize[2];
        if (!result.magnets.empty() &&
            result.particleRetirement.entranceZ <
              lastMagneticExit + margin)
          {
            std::ostringstream message;
            message << "particle_retirement.entrance_z must be at least "
                    << (lastMagneticExit + margin) /
                         result.inputUnits.length
                    << " in the configured length unit: the retirement "
                       "layer must start one lab-equivalent z cell beyond "
                       "the final magnetic interaction region";
            throw configError(retirement["entrance_z"], message.str());
          }
        if (result.detectors.fieldPlanes.empty())
          throw configError(retirement,
            "enabled particle_retirement requires a downstream field plane");
        Double lastFieldPlane =
          -std::numeric_limits<Double>::infinity();
        for (std::size_t i = 0;
             i < result.detectors.fieldPlanes.size(); ++i)
          {
            FieldDetectorPlaneConfig& plane =
              result.detectors.fieldPlanes[i];
            if (plane.particleBackgroundReference)
              throw configError(retirement,
                "particle_retirement and field-plane particle_background "
                "are alternative routes and cannot be enabled together");
            if (!(plane.z > result.particleRetirement.exitZ() + margin))
              {
                std::ostringstream message;
                message << "field detector '" << plane.name
                        << "' must be at least one lab-equivalent z cell "
                           "downstream of the particle-retirement exit at "
                        << result.particleRetirement.exitZ() /
                             result.inputUnits.length;
                throw configError(retirement, message.str());
              }
            plane.particleRetirementCurrent = true;
            plane.particleRetirementEntranceZ =
              result.particleRetirement.entranceZ;
            plane.particleRetirementExitZ =
              result.particleRetirement.exitZ();
            lastFieldPlane = std::max(lastFieldPlane, plane.z);
          }
        if (result.stop.mode != StopMode::ReferenceCenterZ)
          throw configError(stop,
            "particle_retirement requires stop.mode: reference-center-z so "
            "the run does not stop as soon as active particles retire");
        if (!(result.stop.referenceZ > lastFieldPlane + margin))
          {
            std::ostringstream message;
            message << "particle_retirement requires stop.z beyond the "
                       "last field plane; use at least "
                    << (lastFieldPlane + margin) /
                         result.inputUnits.length
                    << " in the configured length unit";
            throw configError(stop["z"], message.str());
          }
      }

    std::size_t protectedRetirementDetectors = 0;
    for (std::size_t detector = 0;
         detector < result.detectors.fieldPlanes.size(); ++detector)
      {
        const FieldDetectorPlaneConfig& plane =
          result.detectors.fieldPlanes[detector];
        if (!plane.retirementFrequencyProtection) continue;
        ++protectedRetirementDetectors;
        if (!result.particleRetirement.enabled)
          throw configError(detectors,
            "retirement_frequency_protection requires enabled top-level particle_retirement");

        const Double scale = std::max(1.0, std::abs(plane.z));
        const Double tolerance = 128.0 *
          std::numeric_limits<Double>::epsilon() * scale;
        for (std::size_t element = 0;
             element < result.beamlineElements.size(); ++element)
          if (result.beamlineElements[element].physicalExit >
              plane.z + tolerance)
            {
              std::ostringstream message;
              message << "frequency-protected field detector '" << plane.name
                      << "' at z=" << plane.z / result.inputUnits.length
                      << " must be the final beamline element; element "
                      << element << " extends downstream to "
                      << result.beamlineElements[element].physicalExit /
                           result.inputUnits.length
                      << " in the configured length unit. Ballistic analytic "
                         "particle_background detectors do not have this "
                         "downstream-element restriction.";
              throw configError(detectors, message.str());
            }
      }
    if (protectedRetirementDetectors > 1)
      throw configError(detectors,
        "only one field detector may bind the global particle_retirement through retirement_frequency_protection");

    if (result.trajectory.enabled && !(result.trajectory.rhythm > 0.0))
      throw configError(trajectory,
        "enabled trajectory output requires positive rhythm");
    if (result.trajectory.compression > 9)
      throw configError(trajectory, "trajectory compression must be in [0,9]");
    if (result.runtime.interactive() && result.trajectory.enabled &&
        result.trajectory.flushEverySamples == 0)
      throw configError(trajectory,
        "interactive trajectory flush_every_samples must be positive");
    if (!result.detectors.energyPlanes.empty())
      {
        std::set<std::string> outputStems;
        for (std::size_t i = 0; i < result.detectors.fieldPlanes.size(); ++i)
          {
            const FieldDetectorPlaneConfig& plane = result.detectors.fieldPlanes[i];
            outputStems.insert(plane.name);
            if (plane.particleBackgroundReference)
              outputStems.insert(plane.name + "-ballistic-reference");
          }
        for (std::size_t i = 0; i < result.detectors.particlePlanes.size(); ++i)
          outputStems.insert(result.detectors.particlePlanes[i].name);
        for (std::size_t i = 0; i < result.detectors.energyPlanes.size(); ++i)
          if (!outputStems.insert(result.detectors.energyPlanes[i].name + "-energy").second)
            throw configError(detectors, "energy plane output filename collides with another detector file");
      }
    if (result.boundary.type == EBBoundaryType::Cpml)
      for (unsigned int axis = 0; axis < 2; ++axis)
        result.detectors.physicalGuardCells[axis] = result.boundary.cpml.cells[axis];
    if (!result.detectors.energyPlanes.empty() &&
        (result.checkpoint.enabled() || !result.checkpoint.restartFrom.empty()))
      throw configError(detectors,
        "energy_planes do not support checkpoint/restart; use an independent diagnostics run");
    return result;
  }
}
