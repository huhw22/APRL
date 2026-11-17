#include "energy_ledger.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "hdf5.h"

namespace
{
  struct ParticlePlaneRecord
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

  void require(bool condition, const std::string& message)
  {
    if (!condition) throw std::runtime_error(message);
  }

  void requireHandle(hid_t handle, const std::string& message)
  {
    if (handle < 0) throw std::runtime_error(message);
  }

  void requireStatus(herr_t status, const std::string& message)
  {
    if (status < 0) throw std::runtime_error(message);
  }

  std::uint64_t readUnsignedScalar(hid_t group, const char* name)
  {
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
    requireHandle(dataset, std::string("Cannot open scalar dataset ") + name);
    std::uint64_t value = 0;
    const herr_t status = H5Dread(dataset, H5T_NATIVE_UINT64,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    H5Dclose(dataset);
    requireStatus(status, std::string("Cannot read scalar dataset ") + name);
    return value;
  }

  unsigned char readByteScalar(hid_t group, const char* name)
  {
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
    requireHandle(dataset, std::string("Cannot open scalar dataset ") + name);
    unsigned char value = 0;
    const herr_t status = H5Dread(dataset, H5T_NATIVE_UCHAR,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    H5Dclose(dataset);
    requireStatus(status, std::string("Cannot read scalar dataset ") + name);
    return value;
  }

  double readDoubleAttribute(hid_t object, const char* name)
  {
    hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
    requireHandle(attribute, std::string("Cannot open attribute ") + name);
    double value = 0.0;
    const herr_t status = H5Aread(attribute, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attribute);
    requireStatus(status, std::string("Cannot read attribute ") + name);
    return value;
  }

  std::string readStringAttribute(hid_t object, const char* name)
  {
    hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
    requireHandle(attribute, std::string("Cannot open attribute ") + name);
    hid_t type = H5Aget_type(attribute);
    requireHandle(type, std::string("Cannot inspect attribute ") + name);
    const std::size_t size = H5Tget_size(type);
    std::vector<char> buffer(size + 1, '\0');
    const herr_t status = H5Aread(attribute, type, buffer.data());
    H5Tclose(type);
    H5Aclose(attribute);
    requireStatus(status, std::string("Cannot read attribute ") + name);
    return std::string(buffer.data());
  }

  std::vector<hsize_t> dimensions(hid_t dataset)
  {
    hid_t space = H5Dget_space(dataset);
    requireHandle(space, "Cannot get dataset dataspace");
    const int rank = H5Sget_simple_extent_ndims(space);
    require(rank >= 0, "Cannot get dataset rank");
    std::vector<hsize_t> result(static_cast<std::size_t>(rank), 0);
    if (rank > 0)
      requireStatus(H5Sget_simple_extent_dims(space, result.data(), nullptr),
        "Cannot get dataset dimensions");
    H5Sclose(space);
    return result;
  }

  hid_t createParticleMemoryType()
  {
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(ParticlePlaneRecord));
    requireHandle(type, "Cannot create particle-plane memory type");
    hsize_t three[1] = {3};
    hid_t vectorType = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, three);
    requireHandle(vectorType, "Cannot create particle-plane vector type");
    try
      {
        requireStatus(H5Tinsert(type, "particle_id",
          HOFFSET(ParticlePlaneRecord, particleId), H5T_NATIVE_UINT64),
          "Cannot map particle_id");
        requireStatus(H5Tinsert(type, "source_id",
          HOFFSET(ParticlePlaneRecord, sourceId), H5T_NATIVE_UINT64),
          "Cannot map source_id");
        requireStatus(H5Tinsert(type, "time_s",
          HOFFSET(ParticlePlaneRecord, time), H5T_NATIVE_DOUBLE),
          "Cannot map time_s");
        requireStatus(H5Tinsert(type, "position_m",
          HOFFSET(ParticlePlaneRecord, position), vectorType),
          "Cannot map position_m");
        requireStatus(H5Tinsert(type, "proper_velocity",
          HOFFSET(ParticlePlaneRecord, properVelocity), vectorType),
          "Cannot map proper velocity");
        requireStatus(H5Tinsert(type, "charge_C",
          HOFFSET(ParticlePlaneRecord, charge), H5T_NATIVE_DOUBLE),
          "Cannot map charge_C");
        requireStatus(H5Tinsert(type, "mass_kg",
          HOFFSET(ParticlePlaneRecord, mass), H5T_NATIVE_DOUBLE),
          "Cannot map mass_kg");
        requireStatus(H5Tinsert(type, "weight",
          HOFFSET(ParticlePlaneRecord, weight), H5T_NATIVE_DOUBLE),
          "Cannot map weight");
        requireStatus(H5Tinsert(type, "crossing_direction",
          HOFFSET(ParticlePlaneRecord, direction), H5T_NATIVE_INT32),
          "Cannot map direction");
      }
    catch (...)
      {
        H5Tclose(vectorType);
        H5Tclose(type);
        throw;
      }
    H5Tclose(vectorType);
    return type;
  }

  std::vector<ParticlePlaneRecord> readParticlePlane(
      const std::string& filename, bool requireComplete)
  {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open particle-plane file: " + filename);
    hid_t group = H5Gopen2(file, "/particle_plane", H5P_DEFAULT);
    if (group < 0)
      {
        H5Fclose(file);
        throw std::runtime_error("Cannot open /particle_plane in " + filename);
      }
    try
      {
        if (requireComplete)
          require(readByteScalar(group, "complete") == 1,
            "Particle-plane output is not marked complete: " + filename);
        const std::uint64_t committed =
          readUnsignedScalar(group, "committed_records");
        hid_t dataset = H5Dopen2(group, "records", H5P_DEFAULT);
        requireHandle(dataset, "Cannot open particle-plane records");
        const std::vector<hsize_t> shape = dimensions(dataset);
        require(shape.size() == 1 && committed <= shape[0],
          "Particle-plane commit marker exceeds dataset extent");
        std::vector<ParticlePlaneRecord> records(
          static_cast<std::size_t>(committed));
        if (committed > 0)
          {
            hid_t memoryType = createParticleMemoryType();
            const herr_t status = H5Dread(dataset, memoryType, H5S_ALL,
              H5S_ALL, H5P_DEFAULT, records.data());
            H5Tclose(memoryType);
            requireStatus(status, "Cannot read particle-plane records");
          }
        H5Dclose(dataset);
        H5Gclose(group);
        H5Fclose(file);
        return records;
      }
    catch (...)
      {
        H5Gclose(group);
        H5Fclose(file);
        throw;
      }
  }

  bool finiteParticle(const ParticlePlaneRecord& record)
  {
    if (!std::isfinite(record.time) || !std::isfinite(record.charge) ||
        !std::isfinite(record.mass) || !std::isfinite(record.weight))
      return false;
    for (int component = 0; component < 3; ++component)
      if (!std::isfinite(record.position[component]) ||
          !std::isfinite(record.properVelocity[component])) return false;
    return true;
  }

  bool nearlyEqual(double left, double right,
                   double absoluteTolerance, double relativeTolerance)
  {
    const double scale = std::max(std::fabs(left), std::fabs(right));
    return std::fabs(left - right) <=
      absoluteTolerance + relativeTolerance * scale;
  }

  int checkDetector(const std::string& fieldFilename,
                    const std::string& particleFilename,
                    std::uint64_t expectedParticles,
                    double minimumInterval,
                    const std::string& manifestFilename,
                    const std::string& transverseLimit = "1e-08")
  {
    hid_t file = H5Fopen(fieldFilename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open field-plane file: " + fieldFilename);
    hid_t group = H5Gopen2(file, "/field_plane", H5P_DEFAULT);
    if (group < 0)
      {
        H5Fclose(file);
        throw std::runtime_error("Cannot open /field_plane");
      }
    require(readByteScalar(group, "complete") == 1,
      "Field-plane output is not marked complete");
    const std::string runId = readStringAttribute(group, "run_id");
    require(!runId.empty(), "Field-plane run_id is empty");
    require(!readStringAttribute(group,
      "configuration_digest_fnv1a64").empty(),
      "Field-plane configuration digest is empty");
    require(!readStringAttribute(group, "source_revision").empty(),
      "Field-plane source revision is empty");
    require(readStringAttribute(group, "particle_background_status") ==
      "disabled; raw total field retained",
      "Omitted particle_background did not remain disabled");
    require(readDoubleAttribute(group, "reference_distance_m") == 0.0,
      "Omitted particle_background created a reference region");
    const std::uint64_t committed =
      readUnsignedScalar(group, "committed_samples");
    require(committed >= 2, "Field-plane cadence test needs at least 2 samples");
    hid_t timeDataset = H5Dopen2(group, "time_s", H5P_DEFAULT);
    requireHandle(timeDataset, "Cannot open field-plane times");
    const std::vector<hsize_t> timeShape = dimensions(timeDataset);
    require(timeShape.size() == 1 && timeShape[0] == committed,
      "Field-plane time extent does not equal committed_samples");
    hid_t fieldDataset = H5Dopen2(group, "fields", H5P_DEFAULT);
    requireHandle(fieldDataset, "Cannot open field-plane fields");
    const std::vector<hsize_t> fieldShape = dimensions(fieldDataset);
    require(fieldShape.size() == 3 && fieldShape[0] == committed,
      "Field-plane field extent does not equal committed_samples");
    std::vector<double> times(static_cast<std::size_t>(committed), 0.0);
    requireStatus(H5Dread(timeDataset, H5T_NATIVE_DOUBLE, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, times.data()), "Cannot read field-plane times");
    double minimumObservedInterval = std::numeric_limits<double>::infinity();
    for (std::size_t sample = 0; sample < times.size(); ++sample)
      {
        require(std::isfinite(times[sample]), "Non-finite detector time");
        if (sample > 0)
          {
            minimumObservedInterval = std::min(minimumObservedInterval,
              times[sample] - times[sample - 1]);
            const double scheduled = times[0] +
              static_cast<double>(sample) * minimumInterval;
            require(times[sample] + minimumInterval * 1.0e-10 >= scheduled,
              "Field-plane samples run ahead of the requested cadence");
          }
      }
    std::cout << "field cadence: samples=" << committed
              << ", minimum interval_s=" << minimumObservedInterval
              << ", requested_s=" << minimumInterval << ", times_s=";
    for (std::size_t sample = 0; sample < times.size(); ++sample)
      std::cout << (sample == 0 ? "[" : ",") << times[sample];
    std::cout << "]\n";
    require(minimumObservedInterval > 0.0,
      "Field-plane times are not strictly increasing");
    H5Dclose(fieldDataset);
    H5Dclose(timeDataset);
    H5Gclose(group);
    H5Fclose(file);

    std::string referenceFilename = fieldFilename;
    if (referenceFilename.size() >= 3 &&
        referenceFilename.substr(referenceFilename.size() - 3) == ".h5")
      referenceFilename.erase(referenceFilename.size() - 3);
    referenceFilename += "-ballistic-reference.h5";
    std::ifstream referenceFile(referenceFilename.c_str(),
      std::ios::in | std::ios::binary);
    require(!referenceFile.good(),
      "Omitted particle_background created a companion reference file");

    std::vector<ParticlePlaneRecord> particles =
      readParticlePlane(particleFilename, true);
    require(particles.size() == expectedParticles,
      "Particle detector did not record the expected crossing count");
    hid_t particleFile = H5Fopen(particleFilename.c_str(), H5F_ACC_RDONLY,
      H5P_DEFAULT);
    requireHandle(particleFile, "Cannot reopen particle-plane file");
    hid_t particleGroup = H5Gopen2(particleFile, "/particle_plane",
      H5P_DEFAULT);
    requireHandle(particleGroup, "Cannot reopen /particle_plane");
    const double planeZ = readDoubleAttribute(particleGroup, "plane_z_m");
    require(readStringAttribute(particleGroup, "run_id") == runId,
      "Field and particle detector files have different run_id values");
    H5Gclose(particleGroup);
    H5Fclose(particleFile);
    for (std::size_t particle = 0; particle < particles.size(); ++particle)
      {
        require(finiteParticle(particles[particle]),
          "Particle detector contains non-finite data");
        require(particles[particle].direction == 1,
          "Expected a downstream detector crossing");
        require(nearlyEqual(particles[particle].position[2], planeZ,
          1.0e-15, 1.0e-12), "Particle crossing is not on the detector plane");
      }
    std::ifstream manifest(manifestFilename.c_str());
    require(static_cast<bool>(manifest), "Cannot open run manifest");
    const std::string manifestText((std::istreambuf_iterator<char>(manifest)),
      std::istreambuf_iterator<char>());
    require(manifestText.find(runId) != std::string::npos,
      "Run manifest does not contain the detector run_id");
    require(manifestText.find("configuration_yaml: |") != std::string::npos,
      "Run manifest does not preserve the input YAML");
    require(manifestText.find("resolved_initial_self_field:\n") !=
              std::string::npos &&
            manifestText.find("  maximum_mean_transverse_beta: " + transverseLimit + "\n") !=
              std::string::npos,
      "Run manifest does not preserve the resolved transverse-beta limit");
    std::cout << "detector: samples=" << committed
              << ", particle crossings=" << particles.size() << '\n';
    return EXIT_SUCCESS;
  }

  void requireCompleteGroup(const std::string& filename,
                            const char* groupName,
                            const char* finiteDataset)
  {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open post-process output: " + filename);
    hid_t group = H5Gopen2(file, groupName, H5P_DEFAULT);
    if (group < 0)
      {
        H5Fclose(file);
        throw std::runtime_error("Missing expected group " +
          std::string(groupName) + " in " + filename);
      }
    require(readByteScalar(group, "complete") == 1,
      "Post-process output is not marked complete: " + filename);
    hid_t dataset = H5Dopen2(group, finiteDataset, H5P_DEFAULT);
    requireHandle(dataset, "Missing numerical result " +
      std::string(finiteDataset) + " in " + filename);
    const std::vector<hsize_t> shape = dimensions(dataset);
    hsize_t count = 1;
    for (std::size_t axis = 0; axis < shape.size(); ++axis)
      count *= shape[axis];
    require(count > 0, "Post-process numerical result is empty: " + filename);
    std::vector<double> values(static_cast<std::size_t>(count), 0.0);
    requireStatus(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, values.data()), "Cannot read numerical result: " + filename);
    for (std::size_t i = 0; i < values.size(); ++i)
      require(std::isfinite(values[i]),
        "Post-process numerical result contains NaN/Inf: " + filename);
    H5Dclose(dataset);
    H5Gclose(group);
    H5Fclose(file);
  }

  int checkPostprocess(const std::string& reconstruction,
                       const std::string& fieldAnalysis,
                       const std::string& trajectory,
                       const std::string& power,
                       const std::string& closure)
  {
    requireCompleteGroup(reconstruction, "/reconstructed_field", "time_s");
    requireCompleteGroup(fieldAnalysis, "/field_plane_analysis",
      "mean_energy_spectrum_J_per_eV");
    requireCompleteGroup(trajectory, "/far_field", "band_energy_J");
    requireCompleteGroup(power, "/power_comparison",
      "signal_forward_power_W");
    std::ifstream report(closure.c_str());
    require(static_cast<bool>(report), "Cannot open energy-closure report");
    const std::string text((std::istreambuf_iterator<char>(report)),
      std::istreambuf_iterator<char>());
    require(text.find("particle_energy_loss_J") != std::string::npos &&
            text.find("forward_radiation_band_energy_J") != std::string::npos,
      "Energy-closure report is missing required budget fields");
    std::cout << "post-process chain: all outputs complete and finite\n";
    return EXIT_SUCCESS;
  }

  std::string readGroupStringAttribute(const std::string& filename,
                                       const char* groupName,
                                       const char* attributeName)
  {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open HDF5 output: " + filename);
    hid_t group = H5Gopen2(file, groupName, H5P_DEFAULT);
    if (group < 0)
      {
        H5Fclose(file);
        throw std::runtime_error("Missing expected group " +
          std::string(groupName) + " in " + filename);
      }
    try
      {
        const std::string value = readStringAttribute(group, attributeName);
        H5Gclose(group);
        H5Fclose(file);
        return value;
      }
    catch (...)
      {
        H5Gclose(group);
        H5Fclose(file);
        throw;
      }
  }

  int checkFieldRoutes(const std::string& rawAnalysis,
                       const std::string& reconstruction,
                       const std::string& reconstructedAnalysis)
  {
    requireCompleteGroup(rawAnalysis, "/field_plane_analysis",
      "mean_energy_spectrum_J_per_eV");
    requireCompleteGroup(reconstruction, "/reconstructed_field", "time_s");
    requireCompleteGroup(reconstructedAnalysis, "/field_plane_analysis",
      "mean_energy_spectrum_J_per_eV");
    require(readGroupStringAttribute(rawAnalysis, "/field_plane_analysis",
      "input_field_kind") == "field_plane",
      "Direct field-analysis route did not consume a raw field plane");
    require(readGroupStringAttribute(reconstruction, "/reconstructed_field",
      "particle_input_role") == "ballistic_reference",
      "Free-drift reconstruction did not consume a ballistic reference");
    require(readGroupStringAttribute(reconstructedAnalysis,
      "/field_plane_analysis", "input_field_kind") ==
      "reconstructed_field",
      "Reconstructed field-analysis route did not consume reconstructed data");
    std::cout << "field routes: raw field and ballistic reconstruction are "
                 "both analyzable\n";
    return EXIT_SUCCESS;
  }

  int compareParticles(const std::string& leftFilename,
                       const std::string& rightFilename,
                       std::size_t expectedRecords)
  {
    std::vector<ParticlePlaneRecord> left =
      readParticlePlane(leftFilename, true);
    std::vector<ParticlePlaneRecord> right =
      readParticlePlane(rightFilename, true);
    const auto order = [](const ParticlePlaneRecord& first,
                          const ParticlePlaneRecord& second) {
      if (first.sourceId != second.sourceId)
        return first.sourceId < second.sourceId;
      if (first.particleId != second.particleId)
        return first.particleId < second.particleId;
      return first.direction < second.direction;
    };
    std::sort(left.begin(), left.end(), order);
    std::sort(right.begin(), right.end(), order);
    require(left.size() == right.size(),
      "MPI runs produced different detector record counts");
    require(left.size() == expectedRecords,
      "MPI consistency runs did not produce the expected crossings");
    double maximumTimeDifference = 0.0;
    double maximumPositionDifference = 0.0;
    double maximumVelocityDifference = 0.0;
    for (std::size_t record = 0; record < left.size(); ++record)
      {
        const ParticlePlaneRecord& a = left[record];
        const ParticlePlaneRecord& b = right[record];
        require(a.sourceId == b.sourceId && a.particleId == b.particleId &&
          a.direction == b.direction, "MPI runs changed particle identity/order");
        maximumTimeDifference = std::max(maximumTimeDifference,
          std::fabs(a.time - b.time));
        for (int component = 0; component < 3; ++component)
          {
            maximumPositionDifference = std::max(maximumPositionDifference,
              std::fabs(a.position[component] - b.position[component]));
            maximumVelocityDifference = std::max(maximumVelocityDifference,
              std::fabs(a.properVelocity[component] -
                        b.properVelocity[component]));
          }
        require(nearlyEqual(a.charge, b.charge, 1.0e-30, 2.0e-14) &&
          nearlyEqual(a.mass, b.mass, 1.0e-40, 2.0e-14) &&
          nearlyEqual(a.weight, b.weight, 1.0e-14, 2.0e-14),
          "MPI runs disagree on particle metadata");
      }
    std::cout << "mpi differences: time_s=" << maximumTimeDifference
              << ", position_m=" << maximumPositionDifference
              << ", proper_velocity=" << maximumVelocityDifference << '\n';
    require(maximumTimeDifference <= 2.0e-22,
      "MPI runs disagree on crossing time");
    require(maximumPositionDifference <= 2.0e-14,
      "MPI runs disagree on crossing position");
    require(maximumVelocityDifference <= 5.0e-8,
      "MPI runs disagree on proper velocity");
    std::cout << "mpi consistency: " << left.size()
              << " crossing records agree\n";
    return EXIT_SUCCESS;
  }

  bool finiteLedger(const aprl::EnergyLedgerRecord& record)
  {
    const double values[] = {
      record.timeBox, record.activeRepresentedElectrons,
      record.particleKineticActive, record.particleTotalActive,
      record.particleKineticRemoved, record.particleTotalRemoved,
      record.fieldEnergyInterior, record.outwardFieldEnergy,
      record.prescribedWork, record.boundarySourceWork,
      record.boundarySourceAbsWork,
      record.balanceResidualRaw, record.balanceResidual,
      record.relativeBalanceExchange, record.relativeBalanceInitial,
      record.fieldFractionKinetic, record.fieldFractionIncludingRest,
      record.meanGammaLab, record.sigmaGammaLab,
      record.relativeEnergySpreadLab, record.minimumGammaLab,
      record.maximumGammaLab, record.linearChirpGammaPerM,
      record.uncorrelatedSigmaGammaLab, record.meanGammaBox,
      record.sigmaGammaBox
    };
    for (double value : values)
      if (!std::isfinite(value)) return false;
    for (double value : record.outwardFieldEnergyFace)
      if (!std::isfinite(value)) return false;
    return true;
  }

  void requireDurationTransaction(const std::string& filename,
      const char* groupName, const char* countName, const char* recordsName,
      const std::vector<hsize_t>& expectedShape, const std::string& runId)
  {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open bounded output: " + filename);
    hid_t group = H5Gopen2(file, groupName, H5P_DEFAULT);
    requireHandle(group, "Cannot open bounded output group: " + filename);
    require(readByteScalar(group, "complete") == 0,
      "Duration-limited output must remain marked incomplete: " + filename);
    require(readUnsignedScalar(group, countName) == expectedShape[0],
      "Bounded output commit count differs: " + filename);
    require(readStringAttribute(group, "run_id") == runId,
      "Bounded outputs have inconsistent run identities: " + filename);
    hid_t records = H5Dopen2(group, recordsName, H5P_DEFAULT);
    requireHandle(records, "Cannot open bounded records: " + filename);
    require(dimensions(records) == expectedShape,
      "Bounded output extent differs from committed records: " + filename);
    H5Dclose(records);
    H5Gclose(group);
    H5Fclose(file);
  }

  int checkDurationClose(const std::string& directory)
  {
    // The fixture advances eight complete field steps. All buffers exceed
    // their final sample counts so explicit close must flush pending data.
    const std::size_t steps = 8;
    const std::size_t samples = steps + 1;
    const std::string fieldName = directory + "/detectors/duration-field.h5";
    hid_t file = H5Fopen(fieldName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open delayed field output");
    hid_t group = H5Gopen2(file, "/field_plane", H5P_DEFAULT);
    requireHandle(group, "Cannot open delayed field group");
    const std::string runId = readStringAttribute(group, "run_id");
    require(!runId.empty(), "Bounded run identity is empty");
    hid_t timeDataset = H5Dopen2(group, "time_s", H5P_DEFAULT);
    requireHandle(timeDataset, "Cannot open bounded field time axis");
    require(dimensions(timeDataset) == std::vector<hsize_t>{samples},
      "Bounded field time extent differs");
    std::vector<double> times(samples);
    requireStatus(H5Dread(timeDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, times.data()), "Cannot read bounded field time axis");
    for (std::size_t i = 0; i < samples; ++i)
      require(std::isfinite(times[i]) && (i == 0 || times[i] > times[i - 1]),
        "Bounded field times are not finite and strictly increasing");
    H5Dclose(timeDataset);
    struct FieldRecord { double electric[3]; double magnetic[3]; };
    hid_t memory = H5Tcreate(H5T_COMPOUND, sizeof(FieldRecord));
    requireHandle(memory, "Cannot create bounded field memory type");
    const hsize_t three[1] = {3};
    hid_t vector = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, three);
    requireHandle(vector, "Cannot create bounded field vector type");
    requireStatus(H5Tinsert(memory, "electric_V_per_m",
      HOFFSET(FieldRecord, electric), vector), "Cannot map bounded electric field");
    requireStatus(H5Tinsert(memory, "magnetic_T",
      HOFFSET(FieldRecord, magnetic), vector), "Cannot map bounded magnetic field");
    H5Tclose(vector);
    hid_t fields = H5Dopen2(group, "fields", H5P_DEFAULT);
    requireHandle(fields, "Cannot open bounded fields");
    require(dimensions(fields) == std::vector<hsize_t>{samples, 8, 8},
      "Unexpected bounded field dimensions");
    std::vector<FieldRecord> values(samples * 8 * 8);
    requireStatus(H5Dread(fields, memory, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      values.data()), "Cannot read bounded field values");
    for (const FieldRecord& value : values)
      for (int axis = 0; axis < 3; ++axis)
        require(std::isfinite(value.electric[axis]) &&
                std::isfinite(value.magnetic[axis]), "Nonfinite bounded field");
    H5Dclose(fields);
    H5Tclose(memory);
    H5Gclose(group);
    H5Fclose(file);
    requireDurationTransaction(fieldName, "/field_plane", "committed_samples",
      "fields", {samples, 8, 8}, runId);

    const char* particleNames[] = {
      "duration-particles", "duration-empty-near", "duration-empty-far"};
    for (std::size_t index = 0; index < 3; ++index)
      {
        const std::string name = directory + "/detectors/" +
          particleNames[index] + ".h5";
        const std::size_t expected = index == 0 ? 4 : 0;
        requireDurationTransaction(name, "/particle_plane", "committed_records",
          "records", {expected}, runId);
        const std::vector<ParticlePlaneRecord> records = readParticlePlane(name, false);
        require(records.size() == expected, "Unexpected bounded particle cohort");
        for (const ParticlePlaneRecord& record : records)
          require(finiteParticle(record), "Nonfinite bounded particle record");
      }

    const std::string ledgerName = directory + "/ledger/energy-ledger.h5";
    requireDurationTransaction(ledgerName, "/energy_ledger", "committed_records",
      "records", {samples}, runId);
    file = H5Fopen(ledgerName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot reopen bounded energy ledger");
    hid_t records = H5Dopen2(file, "/energy_ledger/records", H5P_DEFAULT);
    requireHandle(records, "Cannot open bounded ledger records");
    memory = aprl::EnergyLedgerWriter::createMemoryRecordType();
    std::vector<aprl::EnergyLedgerRecord> ledger(samples);
    requireStatus(H5Dread(records, memory, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      ledger.data()), "Cannot read bounded ledger");
    H5Tclose(memory);
    H5Dclose(records);
    H5Fclose(file);
    const double dt = 0.1e-6 / 299792458.0;
    for (std::size_t i = 0; i < samples; ++i)
      require(ledger[i].step == i && finiteLedger(ledger[i]) &&
        nearlyEqual(ledger[i].timeBox, static_cast<double>(i) * dt, 1e-29, 1e-12),
        "Bounded ledger is nonfinite or has a missing complete step");
    std::cout << "duration-close: all five files readable, complete=0, field="
              << samples << ", particle cohort=4, empty planes=2, ledger="
              << samples << ", final step=" << steps << '\n';
    return EXIT_SUCCESS;
  }

  int checkLedger(const std::string& filename)
  {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open energy ledger: " + filename);
    hid_t group = H5Gopen2(file, "/energy_ledger", H5P_DEFAULT);
    if (group < 0)
      {
        H5Fclose(file);
        throw std::runtime_error("Cannot open /energy_ledger");
      }
    require(readByteScalar(group, "complete") == 1,
      "Energy ledger is not marked complete");
    const std::uint64_t committed =
      readUnsignedScalar(group, "committed_records");
    require(committed >= 2, "Energy ledger needs at least two samples");
    hid_t dataset = H5Dopen2(group, "records", H5P_DEFAULT);
    requireHandle(dataset, "Cannot open energy-ledger records");
    const std::vector<hsize_t> shape = dimensions(dataset);
    require(shape.size() == 1 && shape[0] == committed,
      "Energy-ledger extent does not equal committed_records");
    std::vector<aprl::EnergyLedgerRecord> records(
      static_cast<std::size_t>(committed));
    hid_t memoryType = aprl::EnergyLedgerWriter::createMemoryRecordType();
    requireStatus(H5Dread(dataset, memoryType, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, records.data()), "Cannot read energy-ledger records");
    H5Tclose(memoryType);
    H5Dclose(dataset);
    H5Gclose(group);
    H5Fclose(file);

    double maximumRelativeInitial = 0.0;
    double maximumPrescribedWork = 0.0;
    for (std::size_t sample = 0; sample < records.size(); ++sample)
      {
        require(finiteLedger(records[sample]),
          "Energy ledger contains non-finite data");
        require(records[sample].closureValid == 1,
          "K=0 ledger unexpectedly invalidated closure");
        require(records[sample].boundaryInjectionEntries == 0 &&
                records[sample].boundarySourceWork == 0.0 &&
                records[sample].boundarySourceAbsWork == 0.0 &&
                records[sample].balanceResidualRaw ==
                  records[sample].balanceResidual,
          "Default energy ledger unexpectedly contains boundary source work");
        if (sample > 0)
          require(records[sample].step > records[sample - 1].step &&
            records[sample].timeBox > records[sample - 1].timeBox,
            "Energy-ledger samples are not strictly ordered");
        maximumRelativeInitial = std::max(maximumRelativeInitial,
          std::fabs(records[sample].relativeBalanceInitial));
        maximumPrescribedWork = std::max(maximumPrescribedWork,
          std::fabs(records[sample].prescribedWork));
      }
    require(maximumPrescribedWork <= 1.0e-24,
      "K=0 baseline accumulated unexpected prescribed work");
    require(maximumRelativeInitial <= 1.0e-5,
      "K=0 baseline violates the lightweight energy-ledger sanity bound");
    std::cout << "energy ledger: samples=" << committed
              << ", max |relative balance/initial|="
              << maximumRelativeInitial << '\n';
    return EXIT_SUCCESS;
  }

  aprl::EnergyLedgerRecord readFinalLedger(const std::string& filename)
  {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open energy ledger: " + filename);
    hid_t group = H5Gopen2(file, "/energy_ledger", H5P_DEFAULT);
    requireHandle(group, "Cannot open energy-ledger group: " + filename);
    require(readByteScalar(group, "complete") == 1,
      "Boundary-injection ledger is incomplete");
    const std::uint64_t committed =
      readUnsignedScalar(group, "committed_records");
    require(committed > 1, "Boundary-injection ledger lacks step records");
    hid_t dataset = H5Dopen2(group, "records", H5P_DEFAULT);
    requireHandle(dataset, "Cannot open energy-ledger records");
    std::vector<aprl::EnergyLedgerRecord> records(
      static_cast<std::size_t>(committed));
    hid_t memoryType = aprl::EnergyLedgerWriter::createMemoryRecordType();
    requireStatus(H5Dread(dataset, memoryType, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, records.data()),
      "Cannot read boundary-injection energy ledger");
    H5Tclose(memoryType);
    H5Dclose(dataset);
    H5Gclose(group);
    H5Fclose(file);
    require(records.front().boundaryInjectionEntries == 0,
      "Boundary injection was counted before the first field step");
    require(finiteLedger(records.back()),
      "Boundary-injection ledger contains non-finite data");
    return records.back();
  }

  int checkBoundaryLedger(const std::string& oneRank,
                          const std::string& twoRank)
  {
    const aprl::EnergyLedgerRecord one = readFinalLedger(oneRank);
    const aprl::EnergyLedgerRecord two = readFinalLedger(twoRank);
    const double values[] = {
      one.boundarySourceWork, one.balanceResidualRaw,
      one.balanceResidual
    };
    double scale = 1.0e-30;
    for (double value : values)
      scale = std::max(scale, std::abs(value));
    require(one.closureValid == 0 && two.closureValid == 0,
      "Approximate boundary injection must invalidate exact closure");
    require(one.boundaryInjectionEntries == 1 &&
            two.boundaryInjectionEntries == 1,
      "Expected exactly one injected entry in both MPI layouts");
    require(std::abs(one.boundarySourceWork) > 1.0e-30,
      "Boundary source work was not measured");
    require(one.boundarySourceAbsWork > 1.0e-30 &&
            one.boundarySourceAbsWork + 1.0e-40 >=
              std::abs(one.boundarySourceWork),
      "Absolute boundary-source work did not capture entry activity");
    require(std::abs(one.balanceResidual) > 1.0e-30,
      "Boundary correction unexpectedly forced the residual to zero");
    require(std::abs(one.balanceResidualRaw -
                     one.boundarySourceWork -
                     one.balanceResidual) <= 1.0e-10 * scale,
      "Corrected boundary ledger does not match its separate source term");
    require(std::abs(one.boundarySourceWork - two.boundarySourceWork) <=
              1.0e-8 * scale &&
            std::abs(one.boundarySourceAbsWork -
                     two.boundarySourceAbsWork) <= 1.0e-8 * scale &&
            std::abs(one.balanceResidual - two.balanceResidual) <=
              1.0e-8 * scale,
      "Boundary-source energy depends on the MPI z-slab count");
    std::cout << "boundary ledger: source_J="
              << one.boundarySourceWork
              << " raw_residual_J=" << one.balanceResidualRaw
              << " corrected_residual_J=" << one.balanceResidual << '\n';
    return EXIT_SUCCESS;
  }
}

int main(int argc, char** argv)
{
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  try
    {
      require(argc >= 2,
        "usage: aprl_output_checks <detector|compare-particles|ledger|postprocess|field-routes> ...");
      const std::string command(argv[1]);
      if (command == "detector")
        {
          require(argc == 7,
            "detector needs FIELD_H5 PARTICLE_H5 COUNT MIN_INTERVAL_S MANIFEST");
          return checkDetector(argv[2], argv[3],
            static_cast<std::uint64_t>(std::stoull(argv[4])),
            std::stod(argv[5]), argv[6]);
        }
      if (command == "release-detector")
        {
          require(argc == 3, "release-detector needs RUN_DIRECTORY");
          const std::string path(argv[2]);
          return checkDetector(path + "/detectors/field.h5",
            path + "/detectors/particles.h5", 21, 1e-16,
            path + "/run-manifest.yaml", "0.01");
        }
      if (command == "compare-particles")
        {
          require(argc == 5,
            "compare-particles needs LEFT_H5 RIGHT_H5 EXPECTED_RECORDS");
          return compareParticles(argv[2], argv[3],
            static_cast<std::size_t>(std::stoull(argv[4])));
        }
      if (command == "ledger")
        {
          require(argc == 3, "ledger needs ENERGY_LEDGER_H5");
          return checkLedger(argv[2]);
        }
      if (command == "duration-close")
        {
          require(argc == 3, "duration-close needs RUN_DIRECTORY");
          return checkDurationClose(argv[2]);
        }
      if (command == "boundary-ledger")
        {
          require(argc == 4,
            "boundary-ledger needs ONE_RANK_H5 TWO_RANK_H5");
          return checkBoundaryLedger(argv[2], argv[3]);
        }
      if (command == "postprocess")
        {
          require(argc == 7,
            "postprocess needs RECONSTRUCTION ANALYSIS TRAJECTORY POWER CLOSURE");
          return checkPostprocess(argv[2], argv[3], argv[4], argv[5], argv[6]);
        }
      if (command == "field-routes")
        {
          require(argc == 5,
            "field-routes needs RAW_ANALYSIS RECONSTRUCTION RECONSTRUCTED_ANALYSIS");
          return checkFieldRoutes(argv[2], argv[3], argv[4]);
        }
      throw std::runtime_error("unknown output-check command: " + command);
    }
  catch (const std::exception& error)
    {
      std::cerr << "verification error: " << error.what() << '\n';
      return EXIT_FAILURE;
    }
}
