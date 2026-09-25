#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "hdf5.h"
#include "yaml-cpp/yaml.h"
#include "../yaml_validation.h"

namespace
{
  const long double kC = 299792458.0L;

  void requireHandle(hid_t handle, const std::string& message)
  {
    if (handle < 0) throw std::runtime_error(message);
  }

  void requireStatus(herr_t status, const std::string& message)
  {
    if (status < 0) throw std::runtime_error(message);
  }

  YAML::Node required(const YAML::Node& parent, const char* key)
  {
    const YAML::Node value = parent[key];
    if (!value)
      throw std::runtime_error(std::string("Missing required YAML key '") +
        key + "'");
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
    if (!path.empty() && path[0] == '/') return path;
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

  struct Config
  {
    std::string signalEntry;
    std::string signalExit;
    std::string baselineEntry;
    std::string baselineExit;
    std::string signalEntryField;
    std::string signalExitField;
    std::string baselineEntryField;
    std::string baselineExitField;
    std::string fieldAnalysis;
    std::string report;
    bool requireComplete;
    bool requireAllParticles;
    bool overwrite;
  };

  Config loadConfig(const std::string& filename)
  {
    const YAML::Node root = YAML::LoadFile(filename);
    if (!root.IsMap())
      throw std::runtime_error("Energy-closure card must be a YAML map");
    postprocess_common::validateMapKeys(root, "top level",
      {"input", "output"});
    const YAML::Node input = required(root, "input");
    const YAML::Node output = required(root, "output");
    if (!input.IsMap() || !output.IsMap())
      throw std::runtime_error("input and output must be YAML maps");
    postprocess_common::validateMapKeys(input, "input", {
      "signal_entry_particle_plane", "signal_exit_particle_plane",
      "baseline_entry_particle_plane", "baseline_exit_particle_plane",
      "signal_entry_field_reconstruction",
      "signal_exit_field_reconstruction",
      "baseline_entry_field_reconstruction",
      "baseline_exit_field_reconstruction", "field_analysis",
      "require_complete", "require_all_particles"
    });
    postprocess_common::validateMapKeys(output, "output",
      {"report", "overwrite"});

    Config config;
    config.signalEntry = resolvePath(filename,
      required(input, "signal_entry_particle_plane").as<std::string>());
    config.signalExit = resolvePath(filename,
      required(input, "signal_exit_particle_plane").as<std::string>());
    config.fieldAnalysis = resolvePath(filename,
      required(input, "field_analysis").as<std::string>());
    config.requireComplete = input["require_complete"] ?
      input["require_complete"].as<bool>() : true;
    config.requireAllParticles = input["require_all_particles"] ?
      input["require_all_particles"].as<bool>() : true;
    const bool hasBaselineEntry =
      static_cast<bool>(input["baseline_entry_particle_plane"]);
    const bool hasBaselineExit =
      static_cast<bool>(input["baseline_exit_particle_plane"]);
    if (hasBaselineEntry != hasBaselineExit)
      throw std::runtime_error(
        "baseline entry and exit particle planes must be supplied together");
    if (hasBaselineEntry)
      {
        config.baselineEntry = resolvePath(filename,
          input["baseline_entry_particle_plane"].as<std::string>());
        config.baselineExit = resolvePath(filename,
          input["baseline_exit_particle_plane"].as<std::string>());
      }
    const bool hasSignalEntryField =
      static_cast<bool>(input["signal_entry_field_reconstruction"]);
    const bool hasSignalExitField =
      static_cast<bool>(input["signal_exit_field_reconstruction"]);
    if (hasSignalEntryField != hasSignalExitField)
      throw std::runtime_error(
        "signal entry and exit field reconstructions must be supplied together");
    const bool hasBaselineEntryField =
      static_cast<bool>(input["baseline_entry_field_reconstruction"]);
    const bool hasBaselineExitField =
      static_cast<bool>(input["baseline_exit_field_reconstruction"]);
    if (hasBaselineEntryField != hasBaselineExitField)
      throw std::runtime_error(
        "baseline entry and exit field reconstructions must be supplied together");
    if (hasSignalEntryField)
      {
        config.signalEntryField = resolvePath(filename,
          input["signal_entry_field_reconstruction"].as<std::string>());
        config.signalExitField = resolvePath(filename,
          input["signal_exit_field_reconstruction"].as<std::string>());
      }
    if (hasBaselineEntryField)
      {
        if (!hasBaselineEntry || !hasSignalEntryField)
          throw std::runtime_error(
            "baseline field reconstructions require both particle baseline "
            "planes and signal field reconstructions");
        config.baselineEntryField = resolvePath(filename,
          input["baseline_entry_field_reconstruction"].as<std::string>());
        config.baselineExitField = resolvePath(filename,
          input["baseline_exit_field_reconstruction"].as<std::string>());
      }
    if (hasSignalEntryField && hasBaselineEntry && !hasBaselineEntryField)
      throw std::runtime_error(
        "a matched particle baseline requires matched baseline field "
        "reconstructions for the lab control-volume report");
    config.report = resolvePath(filename,
      required(output, "report").as<std::string>());
    config.overwrite = output["overwrite"] ?
      output["overwrite"].as<bool>() : false;
    postprocess_common::requireOutputAvailable(
      config.report, config.overwrite);
    return config;
  }

  std::uint64_t readUnsignedDataset(hid_t group, const char* name)
  {
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
    requireHandle(dataset, std::string("Missing dataset: ") + name);
    std::uint64_t value = 0;
    const herr_t status = H5Dread(dataset, H5T_NATIVE_UINT64, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &value);
    H5Dclose(dataset);
    requireStatus(status, std::string("Cannot read dataset: ") + name);
    return value;
  }

  unsigned char readByteDataset(hid_t group, const char* name)
  {
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
    requireHandle(dataset, std::string("Missing dataset: ") + name);
    unsigned char value = 0;
    const herr_t status = H5Dread(dataset, H5T_NATIVE_UCHAR, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &value);
    H5Dclose(dataset);
    requireStatus(status, std::string("Cannot read dataset: ") + name);
    return value;
  }

  double readDoubleAttribute(hid_t object, const char* name)
  {
    hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
    requireHandle(attribute, std::string("Missing attribute: ") + name);
    double value = 0.0;
    const herr_t status = H5Aread(attribute, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attribute);
    requireStatus(status, std::string("Cannot read attribute: ") + name);
    return value;
  }

  struct ParticleRecord
  {
    std::uint64_t particleId;
    std::uint64_t sourceId;
    double time;
    double position[3];
    double properVelocity[3];
    double charge;
    double mass;
    double weight;
    std::int32_t direction;
  };

  hid_t createParticleMemoryType()
  {
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(ParticleRecord));
    requireHandle(type, "Cannot create particle memory type");
    hsize_t three[1] = {3};
    hid_t vector = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, three);
    requireHandle(vector, "Cannot create particle vector type");
    requireStatus(H5Tinsert(type, "particle_id",
      HOFFSET(ParticleRecord, particleId), H5T_NATIVE_UINT64),
      "Cannot map particle_id");
    requireStatus(H5Tinsert(type, "source_id",
      HOFFSET(ParticleRecord, sourceId), H5T_NATIVE_UINT64),
      "Cannot map source_id");
    requireStatus(H5Tinsert(type, "time_s", HOFFSET(ParticleRecord, time),
      H5T_NATIVE_DOUBLE), "Cannot map time_s");
    requireStatus(H5Tinsert(type, "position_m",
      HOFFSET(ParticleRecord, position), vector), "Cannot map position_m");
    requireStatus(H5Tinsert(type, "proper_velocity",
      HOFFSET(ParticleRecord, properVelocity), vector),
      "Cannot map proper_velocity");
    requireStatus(H5Tinsert(type, "charge_C",
      HOFFSET(ParticleRecord, charge), H5T_NATIVE_DOUBLE),
      "Cannot map charge_C");
    requireStatus(H5Tinsert(type, "mass_kg",
      HOFFSET(ParticleRecord, mass), H5T_NATIVE_DOUBLE),
      "Cannot map mass_kg");
    requireStatus(H5Tinsert(type, "weight",
      HOFFSET(ParticleRecord, weight), H5T_NATIVE_DOUBLE),
      "Cannot map weight");
    requireStatus(H5Tinsert(type, "crossing_direction",
      HOFFSET(ParticleRecord, direction), H5T_NATIVE_INT32),
      "Cannot map crossing_direction");
    H5Tclose(vector);
    return type;
  }

  struct Plane
  {
    std::map<std::uint64_t, ParticleRecord> downstream;
    std::size_t upstream;
  };

  Plane readPlane(const std::string& filename, bool requireComplete)
  {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open particle plane: " + filename);
    hid_t group = -1;
    hid_t dataset = -1;
    hid_t memoryType = -1;
    try
      {
        group = H5Gopen2(file, "/particle_plane", H5P_DEFAULT);
        requireHandle(group, "Missing /particle_plane in " + filename);
        if (requireComplete && readByteDataset(group, "complete") == 0)
          throw std::runtime_error("Incomplete particle plane: " + filename);
        const std::uint64_t committed =
          readUnsignedDataset(group, "committed_records");
        dataset = H5Dopen2(group, "records", H5P_DEFAULT);
        requireHandle(dataset, "Missing particle-plane records: " + filename);
        hid_t space = H5Dget_space(dataset);
        requireHandle(space, "Cannot inspect particle-plane records");
        hsize_t dimensions[1] = {0};
        const int rank = H5Sget_simple_extent_dims(space, dimensions, 0);
        H5Sclose(space);
        if (rank != 1 || committed > dimensions[0])
          throw std::runtime_error("Invalid committed particle-plane prefix: " +
            filename);
        std::vector<ParticleRecord> records(
          static_cast<std::size_t>(committed));
        memoryType = createParticleMemoryType();
        if (!records.empty())
          {
            hid_t fileSpace = H5Dget_space(dataset);
            const hsize_t count[1] = {committed};
            const hsize_t start[1] = {0};
            requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
              start, 0, count, 0), "Cannot select committed records");
            hid_t memorySpace = H5Screate_simple(1, count, 0);
            requireHandle(memorySpace, "Cannot create record memory space");
            const herr_t status = H5Dread(dataset, memoryType, memorySpace,
              fileSpace, H5P_DEFAULT, &records[0]);
            H5Sclose(memorySpace);
            H5Sclose(fileSpace);
            requireStatus(status, "Cannot read particle-plane records: " +
              filename);
          }

        Plane plane;
        plane.upstream = 0;
        for (std::size_t i = 0; i < records.size(); ++i)
          {
            const ParticleRecord& record = records[i];
            if (record.direction <= 0)
              {
                ++plane.upstream;
                continue;
              }
            if (!plane.downstream.insert(std::make_pair(
                  record.particleId, record)).second)
              throw std::runtime_error(
                "Duplicate downstream crossing for particle " +
                std::to_string(record.particleId) + " in " + filename);
          }
        H5Tclose(memoryType);
        H5Dclose(dataset);
        H5Gclose(group);
        H5Fclose(file);
        return plane;
      }
    catch (...)
      {
        if (memoryType >= 0) H5Tclose(memoryType);
        if (dataset >= 0) H5Dclose(dataset);
        if (group >= 0) H5Gclose(group);
        H5Fclose(file);
        throw;
      }
  }

  class CompensatedSum
  {
  public:
    CompensatedSum() : sum_(0.0L), correction_(0.0L) {}

    void add(long double value)
    {
      const long double updated = sum_ + value;
      if (std::abs(sum_) >= std::abs(value))
        correction_ += (sum_ - updated) + value;
      else
        correction_ += (value - updated) + sum_;
      sum_ = updated;
    }

    long double value() const { return sum_ + correction_; }

  private:
    long double sum_;
    long double correction_;
  };

  struct Losses
  {
    std::map<std::uint64_t, long double> byParticle;
    long double total;
    long double positive;
    long double negative;
    long double meanGammaEntry;
    long double meanGammaExit;
    long double sigmaGammaEntry;
    long double sigmaGammaExit;
    std::size_t missingEntry;
    std::size_t missingExit;
  };

  long double relativeDifference(long double first, long double second)
  {
    return std::abs(first - second) /
      std::max(std::max(std::abs(first), std::abs(second)), 1.0e-300L);
  }

  Losses calculateLosses(const Plane& entry, const Plane& exit,
                         bool requireAll, const std::string& label)
  {
    Losses result;
    result.total = result.positive = result.negative = 0.0L;
    result.meanGammaEntry = result.meanGammaExit = 0.0L;
    result.sigmaGammaEntry = result.sigmaGammaExit = 0.0L;
    result.missingEntry = result.missingExit = 0;
    CompensatedSum total;
    CompensatedSum positive;
    CompensatedSum negative;
    CompensatedSum representedMass;
    CompensatedSum entryGammaMoment;
    CompensatedSum exitGammaMoment;
    for (std::map<std::uint64_t, ParticleRecord>::const_iterator iterator =
         entry.downstream.begin(); iterator != entry.downstream.end();
         ++iterator)
      {
        const std::map<std::uint64_t, ParticleRecord>::const_iterator found =
          exit.downstream.find(iterator->first);
        if (found == exit.downstream.end())
          {
            ++result.missingExit;
            continue;
          }
        const ParticleRecord& before = iterator->second;
        const ParticleRecord& after = found->second;
        if (!(before.mass > 0.0) || !(after.mass > 0.0) ||
            !std::isfinite(before.mass) || !std::isfinite(after.mass))
          throw std::runtime_error("Invalid macro-particle mass in " + label);
        if (relativeDifference(before.mass, after.mass) > 1.0e-12L ||
            relativeDifference(before.charge, after.charge) > 1.0e-12L)
          throw std::runtime_error(
            "Macro-particle mass or charge changed between planes in " + label);

        long double beforeSquared = 0.0L;
        long double afterSquared = 0.0L;
        for (std::size_t axis = 0; axis < 3; ++axis)
          {
            const long double uBefore = before.properVelocity[axis];
            const long double uAfter = after.properVelocity[axis];
            beforeSquared += uBefore * uBefore;
            afterSquared += uAfter * uAfter;
          }
        const long double gammaBefore = std::sqrt(1.0L + beforeSquared);
        const long double gammaAfter = std::sqrt(1.0L + afterSquared);
        const long double deltaGamma = (beforeSquared - afterSquared) /
          (gammaBefore + gammaAfter);
        const long double energy = static_cast<long double>(before.mass) *
          kC * kC * deltaGamma;
        result.byParticle[iterator->first] = energy;
        total.add(energy);
        if (energy >= 0.0L) positive.add(energy);
        else negative.add(energy);
        const long double weight = before.mass;
        representedMass.add(weight);
        entryGammaMoment.add(weight * gammaBefore);
        exitGammaMoment.add(weight * gammaAfter);
      }
    for (std::map<std::uint64_t, ParticleRecord>::const_iterator iterator =
         exit.downstream.begin(); iterator != exit.downstream.end(); ++iterator)
      if (entry.downstream.find(iterator->first) == entry.downstream.end())
        ++result.missingEntry;
    if (requireAll && (result.missingEntry != 0 || result.missingExit != 0))
      {
        std::ostringstream message;
        message << label << " particle planes do not form a complete pair: "
                << result.missingEntry << " missing entries, "
                << result.missingExit << " missing exits";
        throw std::runtime_error(message.str());
      }
    result.total = total.value();
    result.positive = positive.value();
    result.negative = negative.value();
    if (representedMass.value() > 0.0L)
      {
        result.meanGammaEntry =
          entryGammaMoment.value() / representedMass.value();
        result.meanGammaExit =
          exitGammaMoment.value() / representedMass.value();
        CompensatedSum entryVarianceMoment;
        CompensatedSum exitVarianceMoment;
        for (std::map<std::uint64_t, ParticleRecord>::const_iterator iterator =
             entry.downstream.begin(); iterator != entry.downstream.end();
             ++iterator)
          {
            const std::map<std::uint64_t, ParticleRecord>::const_iterator found =
              exit.downstream.find(iterator->first);
            if (found == exit.downstream.end()) continue;
            long double entrySquared = 0.0L;
            long double exitSquared = 0.0L;
            for (std::size_t axis = 0; axis < 3; ++axis)
              {
                const long double entryU = iterator->second.properVelocity[axis];
                const long double exitU = found->second.properVelocity[axis];
                entrySquared += entryU * entryU;
                exitSquared += exitU * exitU;
              }
            const long double entryGamma = std::sqrt(1.0L + entrySquared);
            const long double exitGamma = std::sqrt(1.0L + exitSquared);
            const long double entryDelta = entryGamma - result.meanGammaEntry;
            const long double exitDelta = exitGamma - result.meanGammaExit;
            const long double weight = iterator->second.mass;
            entryVarianceMoment.add(weight * entryDelta * entryDelta);
            exitVarianceMoment.add(weight * exitDelta * exitDelta);
          }
        result.sigmaGammaEntry = std::sqrt(std::max(0.0L,
          entryVarianceMoment.value() / representedMass.value()));
        result.sigmaGammaExit = std::sqrt(std::max(0.0L,
          exitVarianceMoment.value() / representedMass.value()));
      }
    return result;
  }

  struct FieldEnergy
  {
    long double energy;
    double minimumPhotonEnergy;
    double maximumPhotonEnergy;
    double nyquistPhotonEnergy;
    std::string source;
  };

  struct PlaneFieldEnergy
  {
    long double rawSigned;
    long double cleanedSigned;
  };

  PlaneFieldEnergy readPlaneFieldEnergy(const std::string& filename,
                                        bool requireComplete)
  {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open field reconstruction: " + filename);
    hid_t group = H5Gopen2(file, "/reconstructed_field", H5P_DEFAULT);
    if (group < 0)
      {
        H5Fclose(file);
        throw std::runtime_error(
          "Missing /reconstructed_field in " + filename);
      }
    try
      {
        if (requireComplete && readByteDataset(group, "complete") == 0)
          throw std::runtime_error(
            "Incomplete field reconstruction: " + filename);
        PlaneFieldEnergy result;
        result.rawSigned = readDoubleAttribute(group, "raw_signed_energy_J");
        result.cleanedSigned = readDoubleAttribute(
          group, "cleaned_signed_energy_J");
        if (!std::isfinite(static_cast<double>(result.rawSigned)) ||
            !std::isfinite(static_cast<double>(result.cleanedSigned)))
          throw std::runtime_error(
            "Invalid signed energy in field reconstruction: " + filename);
        H5Gclose(group);
        H5Fclose(file);
        return result;
      }
    catch (...)
      {
        H5Gclose(group);
        H5Fclose(file);
        throw;
      }
  }

  struct LabControlVolume
  {
    PlaneFieldEnergy signalEntry;
    PlaneFieldEnergy signalExit;
    PlaneFieldEnergy baselineEntry;
    PlaneFieldEnergy baselineExit;
    bool hasBaseline;
  };

  FieldEnergy readFieldEnergy(const std::string& filename,
                              bool requireComplete)
  {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open field analysis: " + filename);
    const bool fieldPlane = H5Lexists(file, "/field_plane_analysis",
      H5P_DEFAULT) > 0;
    const char* groupPath = fieldPlane ? "/field_plane_analysis" :
      "/far_field";
    hid_t group = H5Gopen2(file, groupPath, H5P_DEFAULT);
    if (group < 0)
      {
        H5Fclose(file);
        throw std::runtime_error(
          "Missing /field_plane_analysis or /far_field in " + filename);
      }
    try
      {
        if (requireComplete && readByteDataset(group, "complete") == 0)
          throw std::runtime_error("Incomplete field analysis: " + filename);
        FieldEnergy result;
        result.energy = readDoubleAttribute(group, "mean_band_energy_J");
        result.minimumPhotonEnergy = readDoubleAttribute(group,
          "band_min_photon_energy_eV");
        result.maximumPhotonEnergy = readDoubleAttribute(group,
          "band_max_photon_energy_eV");
        result.nyquistPhotonEnergy = fieldPlane ? readDoubleAttribute(group,
          "nyquist_photon_energy_eV") :
          std::numeric_limits<double>::quiet_NaN();
        result.source = fieldPlane ? "field_plane_analysis" :
          "trajectory_far_field";
        if (!(result.energy >= 0.0L) || !std::isfinite(
              static_cast<double>(result.energy)))
          throw std::runtime_error("Invalid field energy in " + filename);
        H5Gclose(group);
        H5Fclose(file);
        return result;
      }
    catch (...)
      {
        H5Gclose(group);
        H5Fclose(file);
        throw;
      }
  }

  long double matchedDifference(const Losses& signal,
                                const Losses& baseline,
                                bool requireAll)
  {
    CompensatedSum difference;
    std::size_t missing = 0;
    for (std::map<std::uint64_t, long double>::const_iterator iterator =
         signal.byParticle.begin(); iterator != signal.byParticle.end();
         ++iterator)
      {
        const std::map<std::uint64_t, long double>::const_iterator found =
          baseline.byParticle.find(iterator->first);
        if (found == baseline.byParticle.end())
          {
            ++missing;
            continue;
          }
        difference.add(iterator->second - found->second);
      }
    for (std::map<std::uint64_t, long double>::const_iterator iterator =
         baseline.byParticle.begin(); iterator != baseline.byParticle.end();
         ++iterator)
      if (signal.byParticle.find(iterator->first) == signal.byParticle.end())
        ++missing;
    if (requireAll && missing != 0)
      throw std::runtime_error(
        "Signal and baseline particle pairs have different particle IDs");
    return difference.value();
  }

  void writeReport(std::ostream& output, const Config& config,
                   const Plane& signalEntry, const Plane& signalExit,
                   const Losses& signal, const Plane* baselineEntry,
                   const Plane* baselineExit, const Losses* baseline,
                   long double correctedLoss, const FieldEnergy& field,
                   const LabControlVolume* controlVolume)
  {
    const long double residual = correctedLoss - field.energy;
    const long double ratio = correctedLoss > 0.0L ?
      field.energy / correctedLoss :
      std::numeric_limits<long double>::quiet_NaN();
    const long double baselineMagnitudeRatio = baseline &&
      std::abs(correctedLoss) > 0.0L ?
      std::abs(baseline->total) / std::abs(correctedLoss) :
      std::numeric_limits<long double>::quiet_NaN();
    const long double unresolvedFraction = correctedLoss > 0.0L ?
      residual / correctedLoss :
      std::numeric_limits<long double>::quiet_NaN();
    output << std::setprecision(18) << std::scientific;
    output << "format_version: 3\n"
           << "observer_frame: laboratory\n"
           << "particle_hypersurface: fixed laboratory-z detector planes\n"
           << "boosted_runtime_ledger_used: false\n"
           << "static_magnetic_device_lab_work_J: 0\n"
           << "method: per-particle stable delta-gamma with compensated sums\n"
           << "signal:\n"
           << "  entry_file: " << config.signalEntry << "\n"
           << "  exit_file: " << config.signalExit << "\n"
           << "  paired_downstream_particles: "
           << signal.byParticle.size() << "\n"
           << "  ignored_upstream_entry_crossings: "
           << signalEntry.upstream << "\n"
           << "  ignored_upstream_exit_crossings: "
           << signalExit.upstream << "\n"
           << "  particle_energy_loss_J: " << signal.total << "\n"
           << "  positive_particle_loss_J: " << signal.positive << "\n"
           << "  negative_particle_loss_J: " << signal.negative << "\n"
           << "  entry_mean_gamma_lab: " << signal.meanGammaEntry << "\n"
           << "  exit_mean_gamma_lab: " << signal.meanGammaExit << "\n"
           << "  entry_sigma_gamma_lab: " << signal.sigmaGammaEntry << "\n"
           << "  exit_sigma_gamma_lab: " << signal.sigmaGammaExit << "\n"
           << "  delta_sigma_gamma_lab: "
           << signal.sigmaGammaExit - signal.sigmaGammaEntry << "\n";
    if (baseline && baselineEntry && baselineExit)
      output << "baseline:\n"
             << "  entry_file: " << config.baselineEntry << "\n"
             << "  exit_file: " << config.baselineExit << "\n"
             << "  paired_downstream_particles: "
             << baseline->byParticle.size() << "\n"
             << "  ignored_upstream_entry_crossings: "
             << baselineEntry->upstream << "\n"
             << "  ignored_upstream_exit_crossings: "
             << baselineExit->upstream << "\n"
             << "  particle_energy_loss_J: " << baseline->total << "\n"
             << "  positive_particle_loss_J: " << baseline->positive << "\n"
             << "  negative_particle_loss_J: " << baseline->negative << "\n"
             << "  entry_mean_gamma_lab: " << baseline->meanGammaEntry << "\n"
             << "  exit_mean_gamma_lab: " << baseline->meanGammaExit << "\n"
             << "  entry_sigma_gamma_lab: " << baseline->sigmaGammaEntry << "\n"
             << "  exit_sigma_gamma_lab: " << baseline->sigmaGammaExit << "\n"
             << "  delta_sigma_gamma_lab: "
             << baseline->sigmaGammaExit - baseline->sigmaGammaEntry << "\n";
    output << "closure:\n"
           << "  baseline_corrected_particle_loss_J: " << correctedLoss << "\n"
           << "  radiation_source: " << field.source << "\n"
           << "  forward_radiation_band_energy_J: " << field.energy << "\n"
           << "  radiation_to_particle_loss_ratio: " << ratio << "\n"
           << "  particle_minus_field_J: " << residual << "\n"
           << "  radiation_band_min_eV: " << field.minimumPhotonEnergy << "\n"
           << "  radiation_band_max_eV: " << field.maximumPhotonEnergy << "\n"
           << "  field_plane_nyquist_eV: " << field.nyquistPhotonEnergy << "\n"
           << "lab_energy_budget:\n"
           << "  observed_signal_particle_loss_J: " << signal.total << "\n";
    if (baseline)
      output << "  zero_magnet_control_particle_loss_J: "
             << baseline->total << "\n";
    output << "  device_associated_particle_loss_J: " << correctedLoss << "\n"
           << "  collected_forward_radiation_band_J: " << field.energy << "\n"
           << "  unresolved_device_associated_energy_J: " << residual << "\n"
           << "  forward_fraction_of_device_associated_loss: " << ratio << "\n"
           << "  unresolved_fraction_of_device_associated_loss: "
           << unresolvedFraction << "\n"
           << "  abs_control_exchange_to_abs_device_associated_loss: "
           << baselineMagnitudeRatio << "\n"
           << "  has_zero_magnet_control: "
           << (baseline ? "true" : "false") << "\n"
           << "  near_field_change_directly_measured: false\n"
           << "  energy_spread_is_a_reservoir: false\n"
           << "interpretation:\n"
           << "  - the zero-magnet particle change is a space-charge plus numerical control, not a direct stored-near-field measurement\n"
           << "  - the unresolved term contains differential bound-field change, side or backward radiation, missed aperture or frequency, and numerical residual\n"
           << "  - sigma_gamma describes redistribution inside particle kinetic energy and is not added as a separate energy term\n"
           << "  - only a lab closed-surface flux or equal-lab-time 3D field diagnostic can separate total radiation from retained bound field without this residual\n";
    if (controlVolume)
      {
        const long double signalRaw =
          controlVolume->signalExit.rawSigned -
          controlVolume->signalEntry.rawSigned;
        const long double signalCleaned =
          controlVolume->signalExit.cleanedSigned -
          controlVolume->signalEntry.cleanedSigned;
        const long double baselineRaw = controlVolume->hasBaseline ?
          controlVolume->baselineExit.rawSigned -
            controlVolume->baselineEntry.rawSigned : 0.0L;
        const long double baselineCleaned = controlVolume->hasBaseline ?
          controlVolume->baselineExit.cleanedSigned -
            controlVolume->baselineEntry.cleanedSigned : 0.0L;
        const long double matchedRaw = signalRaw - baselineRaw;
        const long double matchedCleaned = signalCleaned - baselineCleaned;
        const long double rawResidual = correctedLoss - matchedRaw;
        const long double rawFraction = correctedLoss > 0.0L ?
          matchedRaw / correctedLoss :
          std::numeric_limits<long double>::quiet_NaN();
        output << "lab_longitudinal_control_volume:\n"
               << "  sign_convention: positive signed Poynting flux is plus-z\n"
               << "  signal_entry_file: " << config.signalEntryField << "\n"
               << "  signal_exit_file: " << config.signalExitField << "\n"
               << "  signal_entry_raw_signed_J: "
               << controlVolume->signalEntry.rawSigned << "\n"
               << "  signal_exit_raw_signed_J: "
               << controlVolume->signalExit.rawSigned << "\n"
               << "  signal_exit_minus_entry_raw_signed_J: "
               << signalRaw << "\n"
               << "  signal_exit_minus_entry_cleaned_signed_J: "
               << signalCleaned << "\n";
        if (controlVolume->hasBaseline)
          output << "  baseline_entry_file: "
                 << config.baselineEntryField << "\n"
                 << "  baseline_exit_file: "
                 << config.baselineExitField << "\n"
                 << "  baseline_entry_raw_signed_J: "
                 << controlVolume->baselineEntry.rawSigned << "\n"
                 << "  baseline_exit_raw_signed_J: "
                 << controlVolume->baselineExit.rawSigned << "\n"
                 << "  baseline_exit_minus_entry_raw_signed_J: "
                 << baselineRaw << "\n"
                 << "  baseline_exit_minus_entry_cleaned_signed_J: "
                 << baselineCleaned << "\n";
        output << "  matched_exit_minus_entry_raw_signed_J: "
               << matchedRaw << "\n"
               << "  matched_exit_minus_entry_cleaned_signed_J: "
               << matchedCleaned << "\n"
               << "  raw_signed_to_particle_loss_ratio: "
               << rawFraction << "\n"
               << "  particle_minus_raw_signed_J: " << rawResidual << "\n"
               << "  raw_signed_minus_forward_band_J: "
               << matchedRaw - field.energy << "\n"
               << "  transverse_flux_and_stored_field_included: false\n"
               << "  interpretation: raw closure tests the resolved longitudinal lab control volume; raw-minus-forward-band contains bound or non-propagating field, excluded frequencies or angles, and decomposition cross terms\n";
      }
  }
}

int main(int argc, char** argv)
{
  if (argc != 2)
    {
      std::cerr << "Usage: energy_closure <closure.yaml>\n";
      return 2;
    }
  try
    {
      const Config config = loadConfig(argv[1]);
      const Plane signalEntry = readPlane(config.signalEntry,
        config.requireComplete);
      const Plane signalExit = readPlane(config.signalExit,
        config.requireComplete);
      const Losses signal = calculateLosses(signalEntry, signalExit,
        config.requireAllParticles, "signal");

      Plane baselineEntry;
      Plane baselineExit;
      Losses baseline;
      Losses* baselinePointer = 0;
      Plane* baselineEntryPointer = 0;
      Plane* baselineExitPointer = 0;
      long double correctedLoss = signal.total;
      if (!config.baselineEntry.empty())
        {
          baselineEntry = readPlane(config.baselineEntry,
            config.requireComplete);
          baselineExit = readPlane(config.baselineExit,
            config.requireComplete);
          baseline = calculateLosses(baselineEntry, baselineExit,
            config.requireAllParticles, "baseline");
          correctedLoss = matchedDifference(signal, baseline,
            config.requireAllParticles);
          baselinePointer = &baseline;
          baselineEntryPointer = &baselineEntry;
          baselineExitPointer = &baselineExit;
        }
      const FieldEnergy field = readFieldEnergy(config.fieldAnalysis,
        config.requireComplete);
      LabControlVolume controlVolume;
      LabControlVolume* controlVolumePointer = 0;
      if (!config.signalEntryField.empty())
        {
          controlVolume.signalEntry = readPlaneFieldEnergy(
            config.signalEntryField, config.requireComplete);
          controlVolume.signalExit = readPlaneFieldEnergy(
            config.signalExitField, config.requireComplete);
          controlVolume.hasBaseline = !config.baselineEntryField.empty();
          if (controlVolume.hasBaseline)
            {
              controlVolume.baselineEntry = readPlaneFieldEnergy(
                config.baselineEntryField, config.requireComplete);
              controlVolume.baselineExit = readPlaneFieldEnergy(
                config.baselineExitField, config.requireComplete);
            }
          controlVolumePointer = &controlVolume;
        }

      createDirectories(parentDirectory(config.report));
      std::ofstream report(config.report.c_str());
      if (!report)
        throw std::runtime_error("Cannot create report: " + config.report);
      writeReport(report, config, signalEntry, signalExit, signal,
        baselineEntryPointer, baselineExitPointer, baselinePointer,
        correctedLoss, field, controlVolumePointer);
      report.close();
      if (!report)
        throw std::runtime_error("Cannot finish report: " + config.report);

      writeReport(std::cout, config, signalEntry, signalExit, signal,
        baselineEntryPointer, baselineExitPointer, baselinePointer,
        correctedLoss, field, controlVolumePointer);
      return 0;
    }
  catch (const std::exception& error)
    {
      std::cerr << "Energy closure failed: " << error.what() << "\n";
      return 1;
    }
}
