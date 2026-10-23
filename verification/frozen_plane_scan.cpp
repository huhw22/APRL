#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>
#include <sys/stat.h>
#include <unistd.h>

#include "config.h"
#include "eb_deposition.h"
#include "eb_bunch.h"
#include "energy_ledger.h"
#include "simulation.h"

namespace
{
  struct InputRecord
  {
    double position[3];
    double properVelocity[3];
    std::uint64_t sourceId;
    double macroWeight;
  };

  struct DetectorRecord
  {
    std::uint64_t sourceId;
    double time;
    double position[3];
  };

  struct FieldRecord
  {
    double electric[3];
    double magnetic[3];
  };

  void check(bool valid, const std::string& message)
  {
    if (!valid) throw std::runtime_error(message);
  }

  void makeDirectory(const std::string& path)
  {
    check(::mkdir(path.c_str(), 0755) == 0 ||
      access(path.c_str(), F_OK) == 0,
      "Cannot create directory " + path);
  }

  double maximumAbsolute(const aprl::YeeComponent& component)
  {
    double result = 0.0;
    for (std::size_t index = 0; index < component.size(); ++index)
      result = std::max(result, std::abs(component.data()[index]));
    return result;
  }

  void checkHeldReleaseCurrent()
  {
    const aprl::EBGridGeometry geometry(16, 16, 16,
      1.0e-6, 1.0e-6, 1.0e-6, 1.0e-15,
      aprl::EBMaxwellSolver::Yee);
    aprl::EBFieldGrid fields(geometry);
    aprl::ChargeConservingCurrentDepositor depositor(
      fields, aprl::FieldVector<double>(0.0));
    aprl::YeeComponent before(geometry.nx + 1, geometry.ny + 1,
                              geometry.nz + 1);
    aprl::YeeComponent after(geometry.nx + 1, geometry.ny + 1,
                             geometry.nz + 1);
    before.fill(0.0);
    after.fill(0.0);
    aprl::FieldVector<double> start(0.0);
    start[0] = 6.2e-6;
    start[1] = 7.1e-6;
    start[2] = 4.3e-6;
    aprl::FieldVector<double> release(start);
    release[2] = 6.7e-6;
    aprl::FieldVector<double> end(release);
    end[0] += 0.8e-6;
    end[1] -= 0.55e-6;
    end[2] += 1.6e-6;
    const double releaseFraction =
      (release[2] - start[2]) / (end[2] - start[2]);
    check(releaseFraction > 0.0 && releaseFraction < 1.0 &&
          std::abs(releaseFraction - 0.6) < 1.0e-14,
      "Held release must occur strictly inside the field step");
    const double charge = -2.0e-12;
    depositor.depositCharge(start, charge, before);
    fields.clearCurrent();
    depositor.depositSegment(start, release, charge);
    check(maximumAbsolute(fields.jx()) == 0.0 &&
          maximumAbsolute(fields.jy()) == 0.0 &&
          maximumAbsolute(fields.jz()) > 0.0,
      "Frozen leg generated transverse current before release");
    depositor.depositSegment(release, end, charge);
    check(maximumAbsolute(fields.jx()) > 0.0 &&
          maximumAbsolute(fields.jy()) > 0.0,
      "Active leg lost released transverse current");
    depositor.depositCharge(end, charge, after);
    const double residual =
      aprl::ChargeConservingCurrentDepositor::maxContinuityResidual(
        fields, before, after);
    const double scale = std::max(
      maximumAbsolute(before), maximumAbsolute(after)) / geometry.dt;
    check(scale > 0.0 && residual / scale < 3.0e-13,
      "Held-release segmented source violates discrete continuity");
    std::cout << "held_release_current_continuity_relative="
              << residual / scale << " frozen_jx_jy=0"
              << " released_jx_max=" << maximumAbsolute(fields.jx())
              << " released_jy_max=" << maximumAbsolute(fields.jy())
              << '\n';
  }

  void checkHighBoostPlaneGeometry()
  {
    std::vector<aprl::RelativisticParticleSI> frozen(2);
    for (std::size_t i = 0; i < frozen.size(); ++i)
      {
        frozen[i].position[0] = (i == 0 ? -1.0 : 1.0) * 1.0e-4;
        frozen[i].position[1] = 0.0;
        frozen[i].position[2] = (i == 0 ? -1.0 : 1.0) * 6.7e-14;
        frozen[i].properVelocity[0] =
          i == 0 ? 0.117417 : 0.5;
        frozen[i].properVelocity[1] = i == 0 ? 0.0 : 0.2;
        frozen[i].properVelocity[2] = 1174.17;
        frozen[i].weight = 1.0;
        frozen[i].charge = -aprl::SI::elementaryCharge;
        frozen[i].mass = aprl::SI::electronMass;
      }
    std::vector<aprl::RelativisticParticleSI> held(frozen);
    const double inputPlane = -0.05;
    const double center = -0.04995;
    aprl::SIBunchPreprocessor::projectLabPlaneEventsToSnapshot(
      frozen, inputPlane, center, 0.0, MPI_COMM_WORLD);
    aprl::SIBunchPreprocessor::projectLabPlaneEventsToSnapshot(
      held, inputPlane, center, 0.0, MPI_COMM_WORLD, true);
    aprl::SIBunchPlacement placement;
    placement.referencePositionLab = center;
    const aprl::SIBunchPlacementReport placed =
      aprl::SIBunchPreprocessor::placeRelativeLabSnapshot(
        frozen, placement, MPI_COMM_WORLD);
    aprl::SIBunchPreprocessor::placeRelativeLabSnapshot(
      held, placement, MPI_COMM_WORLD);
    aprl::BoostFrameTransform frame;
    frame.setOriginsFromGamma(929.130549, aprl::SI::c, 0.0,
      placed.headAfterLab, 0.0);
    aprl::SIBunchPreprocessor::boostLabSnapshotToBoxTimeZero(
      frozen, frame, 0.0);
    double heldInputTimes[2] = {};
    double heldBetas[2] = {};
    double heldFrameOriginTime = std::numeric_limits<double>::infinity();
    const double heldDzBox = 1.25e-6;
    for (std::size_t i = 0; i < held.size(); ++i)
      {
        const double labGamma =
          aprl::BoostFrameTransform::gammaFromProperVelocity(
            held[i].properVelocity);
        const double labBetaZ = held[i].properVelocity[2] / labGamma;
        heldInputTimes[i] =
          (inputPlane - held[i].position[2]) /
          (aprl::SI::c * labBetaZ);
        aprl::FieldVector<double> boxProperVelocity(0.0);
        frame.properVelocityLabToBox(
          held[i].properVelocity, boxProperVelocity);
        const double betaHeld = frame.heldLongitudinalBetaZ(
          boxProperVelocity);
        heldBetas[i] = betaHeld;
        heldFrameOriginTime = std::min(heldFrameOriginTime,
          heldInputTimes[i] - 2.0 * heldDzBox /
            (frame.gamma() * aprl::SI::c *
              (frame.beta() + betaHeld)));
      }
    aprl::BoostFrameTransform heldFrame;
    check(std::abs(heldBetas[0] - heldBetas[1]) < 1.0e-13,
      "Held longitudinal speed depends on transverse momentum");
    heldFrame.setOriginsFromGamma(929.130549, aprl::SI::c,
      heldFrameOriginTime, inputPlane, 0.0);
    aprl::SIBunchPreprocessor::boostLabSnapshotToBoxTimeZero(
      held, heldFrame, 0.0, true, inputPlane);
    double maximumVirtualShift = 0.0;
    double maximumInputPlaneError = 0.0;
    double maximumHeldInputTimeError = 0.0;
    for (std::size_t i = 0; i < frozen.size(); ++i)
      {
        const double inputX = (i == 0 ? -1.0 : 1.0) * 1.0e-4;
        maximumVirtualShift = std::max(maximumVirtualShift,
          std::abs(frozen[i].position[0] - inputX));
        const double gamma =
          aprl::BoostFrameTransform::gammaFromProperVelocity(
            frozen[i].properVelocity);
        const double betaBoxZ = frozen[i].properVelocity[2] / gamma;
        const double labZ = frame.labZFromBoxZT(
          frozen[i].position[2], 0.0);
        const double labZPerBoxTime = frame.gamma() * aprl::SI::c *
          (frame.beta() + betaBoxZ);
        const double crossingTime =
          (inputPlane - labZ) / labZPerBoxTime;
        const double crossingX = frozen[i].position[0] +
          aprl::SI::c * crossingTime *
          frozen[i].properVelocity[0] / gamma;
        maximumInputPlaneError = std::max(maximumInputPlaneError,
          std::abs(crossingX - inputX));
        check(held[i].position[0] == inputX,
          "Held-transverse boost changed its input-plane x");
        const double heldBetaBox = heldFrame.heldLongitudinalBetaZ(
          held[i].properVelocity);
        const double heldLabZ = heldFrame.labZFromBoxZT(
          held[i].position[2], 0.0);
        check(heldLabZ < inputPlane &&
              held[i].position[2] < -1.9 * heldDzBox,
          "Held high-boost particle starts past its activation plane");
        const double heldCrossingBoxTime =
          (inputPlane - heldLabZ) /
          (heldFrame.gamma() * aprl::SI::c *
            (heldFrame.beta() + heldBetaBox));
        const double heldCrossingZ = held[i].position[2] +
          aprl::SI::c * heldBetaBox * heldCrossingBoxTime;
        maximumHeldInputTimeError = std::max(
          maximumHeldInputTimeError,
          std::abs(heldFrame.labTFromBoxZT(
            heldCrossingZ, heldCrossingBoxTime) - heldInputTimes[i]));
      }
    check(maximumVirtualShift > 1.0e-3,
      "High-boost transform did not exercise millimetre-scale transverse drift");
    check(maximumInputPlaneError < 1.0e-10,
      "Ballistic high-boost trajectory misses its input-plane x");
    check(maximumHeldInputTimeError < 1.0e-18,
      "Held high-boost pseudo-trajectory misses its input-plane arrival time");
    std::cout << "high_boost_virtual_x_shift_mm="
              << maximumVirtualShift * 1.0e3
              << " high_boost_input_plane_max_error_nm="
              << maximumInputPlaneError * 1.0e9
              << " held_high_boost_input_time_max_error_fs="
              << maximumHeldInputTimeError * 1.0e15 << '\n';
  }

  hid_t inputType(bool fileType)
  {
    hid_t type = H5Tcreate(H5T_COMPOUND,
      fileType ? 8 * 8 : sizeof(InputRecord));
    check(type >= 0, "Cannot create input datatype");
    hsize_t two[1] = {2};
    hsize_t three[1] = {3};
    hid_t xy = H5Tarray_create2(
      fileType ? H5T_IEEE_F64LE : H5T_NATIVE_DOUBLE, 1, two);
    hid_t velocity = H5Tarray_create2(
      fileType ? H5T_IEEE_F64LE : H5T_NATIVE_DOUBLE, 1, three);
    check(xy >= 0 && velocity >= 0, "Cannot create vector datatypes");
    const std::size_t position = fileType ? 0 : HOFFSET(InputRecord, position);
    const std::size_t proper = fileType ? 24 : HOFFSET(InputRecord, properVelocity);
    const std::size_t id = fileType ? 48 : HOFFSET(InputRecord, sourceId);
    const std::size_t weight = fileType ? 56 : HOFFSET(InputRecord, macroWeight);
    check(H5Tinsert(type, "plane_position_m", position, xy) >= 0 &&
      H5Tinsert(type, "arrival_time_offset_s", position + 16,
        fileType ? H5T_IEEE_F64LE : H5T_NATIVE_DOUBLE) >= 0 &&
      H5Tinsert(type, "proper_velocity", proper, velocity) >= 0 &&
      H5Tinsert(type, "source_id", id,
        fileType ? H5T_STD_U64LE : H5T_NATIVE_UINT64) >= 0 &&
      H5Tinsert(type, "macro_weight", weight,
        fileType ? H5T_IEEE_F64LE : H5T_NATIVE_DOUBLE) >= 0,
      "Cannot define input columns");
    H5Tclose(xy);
    H5Tclose(velocity);
    return type;
  }

  std::vector<InputRecord> makeInput(const std::string& path,
                                     double transverseScale = 0.004)
  {
    std::vector<InputRecord> records(4096);
    for (std::size_t i = 0; i < records.size(); ++i)
      {
        InputRecord& record = records[i];
        const int ix = static_cast<int>(i % 16);
        const int iy = static_cast<int>((i / 16) % 16);
        const int it = static_cast<int>((i / 256) % 16);
        record.position[0] = (ix - 7.5) * 1.0e-6;
        record.position[1] = (iy - 7.5) * 1.0e-6;
        record.position[2] = (it - 7.5) * 0.8e-15;
        record.properVelocity[0] = (ix - 7.5) * transverseScale;
        record.properVelocity[1] = (iy - 7.5) * transverseScale;
        record.properVelocity[2] = 10.0;
        record.sourceId = i + 1;
        record.macroWeight = 1.0;
      }
    hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC,
      H5P_DEFAULT, H5P_DEFAULT);
    check(file >= 0, "Cannot create HDF5 fixture");
    hid_t group = H5Gcreate2(file, "/particles", H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);
    check(group >= 0, "Cannot create particle group");
    hid_t scalar = H5Screate(H5S_SCALAR);
    hid_t attribute = H5Acreate2(group, "format_version", H5T_STD_I32LE,
      scalar, H5P_DEFAULT, H5P_DEFAULT);
    int version = 4;
    check(attribute >= 0 &&
      H5Awrite(attribute, H5T_NATIVE_INT, &version) >= 0,
      "Cannot write HDF5 version");
    H5Aclose(attribute);
    H5Sclose(scalar);
    hsize_t count[1] = {records.size()};
    hid_t space = H5Screate_simple(1, count, NULL);
    hid_t fileType = inputType(true);
    hid_t memoryType = inputType(false);
    hid_t dataset = H5Dcreate2(group, "records", fileType, space,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    check(dataset >= 0 && H5Dwrite(dataset, memoryType, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &records[0]) >= 0,
      "Cannot write HDF5 particle records");
    H5Dclose(dataset);
    H5Tclose(fileType);
    H5Tclose(memoryType);
    H5Sclose(space);
    H5Gclose(group);
    H5Fclose(file);
    return records;
  }

  void writeCard(const std::string& path, const std::string& input,
                 bool frozen, bool feedback, bool initialField,
                 bool automaticZ, bool holdTransverse,
                 bool activationAtInput)
  {
    std::ofstream out(path.c_str());
    check(static_cast<bool>(out), "Cannot write scan card");
    out << "units: {length: micrometer, time: femtosecond}\n"
        << "runtime: {mode: throughput, resource_monitor: {enabled: false}}\n"
        << "output: {overwrite: true, manifest: run-manifest.yaml}\n"
        << "mesh:\n  field_solver: cowan-z\n  current_filter: none\n"
        << "  cells: [32, 32, "
        << (automaticZ ? "auto" : "160")
        << "]\n  cell_size: [2, 2, 1]\n"
        << "  center: [0, 0, 0]\n  duration: 250\n"
        << "  boost_gamma: 2\n  particle_steps_per_undulator_period: 1\n"
        << "  self_field_feedback: " << (feedback ? "true" : "false") << "\n"
        << (automaticZ ?
          "  longitudinal: {mode: observation-bounded, lower_guard_cells: 4, upper_guard_cells: 8, self_field_padding_rms: 4, ballistic_time_safety_factor: 1.1, alignment_cells: 1}\n" : "")
        << "boundary:\n  type: cpml\n  cells: [4, 4, 4]\n"
        << "  polynomial_order: 3\n  target_reflection: 1.0e-8\n"
        << "  kappa_max: 8\n  alpha_fraction: 0\n"
        << "initial_self_field:\n  enabled: "
        << (initialField ? "true" : "false") << "\n"
        << "  model: relativistic-poisson\n"
        << "  maximum_mean_transverse_beta: 0.1\n"
        << "beam:\n  reference:\n    input_plane_z: 0\n"
        << "    initial_center_z: 8\n"
        << "    frozen_lab_plane_injection: "
        << (frozen ? "true" : "false") << "\n"
        << (frozen ? (activationAtInput ?
              "    activation_plane_z: 0\n" :
              "    activation_plane_z: -1\n") : "")
        << (holdTransverse ?
              "    hold_transverse_until_activation: true\n" : "")
        << "  input:\n    type: hdf5\n    file: " << input << "\n"
        << "    electrons: 1000000\n"
        << "    position_offset: [0, 0, 0]\n"
        << "sources: {incident_waves: [], magnetic_elements: []}\n"
        << "detectors:\n  enabled: true\n  directory: detectors\n"
        << "  field_planes:\n"
        << "    - {name: self-field, z: 12, rhythm: 5, buffer_samples: 32, compression: 0}\n"
        << "  particle_planes:\n"
        << "    - {name: upstream, z: -1, buffer_records: 4096, compression: 0}\n"
        << "    - {name: input, z: 0, buffer_records: 4096, compression: 0}\n"
        << "    - {name: near, z: 12, buffer_records: 4096, compression: 0}\n"
        << "    - {name: far, z: 20, buffer_records: 4096, compression: 0}\n"
        << "stop: {mode: after-last-element}\n"
        << "trajectory: {enabled: false}\n"
        << "energy_ledger:\n  enabled: true\n  directory: ledger\n"
        << "  filename: energy-ledger.h5\n"
        << "  sample_interval_steps: 1\n"
        << "  buffer_records: 128\n  compression: 0\n";
    check(static_cast<bool>(out), "Cannot finish scan card");
  }

  std::vector<DetectorRecord> readDetector(const std::string& path)
  {
    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    check(file >= 0, "Cannot open detector " + path);
    hid_t dataset = H5Dopen2(file, "/particle_plane/records", H5P_DEFAULT);
    check(dataset >= 0, "Cannot open detector records");
    hid_t space = H5Dget_space(dataset);
    hsize_t count[1] = {0};
    check(H5Sget_simple_extent_dims(space, count, NULL) == 1,
      "Cannot inspect detector records");
    std::vector<DetectorRecord> records(count[0]);
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(DetectorRecord));
    hsize_t three[1] = {3};
    hid_t vectorType = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, three);
    check(H5Tinsert(type, "source_id", HOFFSET(DetectorRecord, sourceId),
      H5T_NATIVE_UINT64) >= 0 &&
      H5Tinsert(type, "time_s", HOFFSET(DetectorRecord, time),
      H5T_NATIVE_DOUBLE) >= 0 &&
      H5Tinsert(type, "position_m", HOFFSET(DetectorRecord, position),
      vectorType) >= 0, "Cannot define detector read type");
    if (!records.empty())
      check(H5Dread(dataset, type, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &records[0]) >= 0, "Cannot read detector records");
    H5Tclose(vectorType);
    H5Tclose(type);
    H5Sclose(space);
    H5Dclose(dataset);
    H5Fclose(file);
    return records;
  }

  struct FieldSummary
  {
    double initial;
    double maximum;
    double final;
    double maximumStepIncrease;
  };

  void readFieldAmplitude(const std::string& path,
                          double& maximumElectric,
                          double& maximumLightScaledMagnetic)
  {
    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    check(file >= 0, "Cannot open field detector");
    hid_t dataset = H5Dopen2(file, "/field_plane/fields", H5P_DEFAULT);
    check(dataset >= 0, "Cannot open field detector samples");
    hid_t space = H5Dget_space(dataset);
    hsize_t dimensions[3] = {0, 0, 0};
    check(H5Sget_simple_extent_dims(space, dimensions, NULL) == 3 &&
      dimensions[0] > 0, "Field detector has no samples");
    const std::size_t count = static_cast<std::size_t>(dimensions[0] *
      dimensions[1] * dimensions[2]);
    std::vector<FieldRecord> records(count);
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(FieldRecord));
    hsize_t three[1] = {3};
    hid_t vectorType = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, three);
    check(H5Tinsert(type, "electric_V_per_m",
      HOFFSET(FieldRecord, electric), vectorType) >= 0 &&
      H5Tinsert(type, "magnetic_T", HOFFSET(FieldRecord, magnetic),
      vectorType) >= 0, "Cannot define field detector read type");
    check(H5Dread(dataset, type, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &records[0]) >= 0, "Cannot read field detector");
    maximumElectric = 0.0;
    maximumLightScaledMagnetic = 0.0;
    for (std::size_t i = 0; i < records.size(); ++i)
      {
        double electric2 = 0.0;
        double magnetic2 = 0.0;
        for (unsigned int axis = 0; axis < 3; ++axis)
          {
            electric2 += records[i].electric[axis] *
              records[i].electric[axis];
            magnetic2 += records[i].magnetic[axis] *
              records[i].magnetic[axis];
          }
        maximumElectric = std::max(maximumElectric,
          std::sqrt(electric2));
        maximumLightScaledMagnetic = std::max(
          maximumLightScaledMagnetic,
          aprl::SI::c * std::sqrt(magnetic2));
      }
    H5Tclose(vectorType);
    H5Tclose(type);
    H5Sclose(space);
    H5Dclose(dataset);
    H5Fclose(file);
    check(std::isfinite(maximumElectric) &&
      std::isfinite(maximumLightScaledMagnetic) &&
      maximumElectric < 1.0e11 &&
      maximumLightScaledMagnetic < 1.0e11,
      "Free-drift self-field amplitude is nonfinite or pathological");
  }

  FieldSummary readFieldEnergy(const std::string& path)
  {
    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    check(file >= 0, "Cannot open energy ledger");
    hid_t dataset = H5Dopen2(file, "/energy_ledger/records", H5P_DEFAULT);
    check(dataset >= 0, "Cannot open energy ledger records");
    hid_t space = H5Dget_space(dataset);
    hsize_t count[1] = {0};
    check(H5Sget_simple_extent_dims(space, count, NULL) == 1 && count[0] > 2,
      "Energy ledger has too few samples");
    std::vector<aprl::EnergyLedgerRecord> records(count[0]);
    hid_t type = aprl::EnergyLedgerWriter::createMemoryRecordType();
    check(H5Dread(dataset, type, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &records[0]) >= 0, "Cannot read energy ledger");
    FieldSummary result = {records.front().fieldEnergyInterior, 0.0,
                           records.back().fieldEnergyInterior, 0.0};
    for (std::size_t i = 0; i < records.size(); ++i)
      {
        check(std::isfinite(records[i].fieldEnergyInterior) &&
          records[i].fieldEnergyInterior >= 0.0,
          "Nonfinite or negative field energy");
        result.maximum = std::max(result.maximum,
          records[i].fieldEnergyInterior);
        if (i != 0)
          result.maximumStepIncrease = std::max(
            result.maximumStepIncrease,
            records[i].fieldEnergyInterior -
              records[i - 1].fieldEnergyInterior);
      }
    H5Tclose(type);
    H5Sclose(space);
    H5Dclose(dataset);
    H5Fclose(file);
    return result;
  }

  void checkIntegratedHeldFrameOrigin(
      const std::string& output,
      const std::vector<InputRecord>& input)
  {
    const std::string label = "origin_lab_t=";
    const std::size_t start = output.find(label);
    check(start != std::string::npos,
      "Held simulation did not report its laboratory frame origin");
    std::istringstream value(output.substr(start + label.size()));
    double actual = 0.0;
    value >> actual;
    check(static_cast<bool>(value) && std::isfinite(actual),
      "Held simulation reported an invalid laboratory frame origin");

    const double referenceZ = 8.0e-6;
    const double dz = 1.0e-6;
    const double boostGamma = 2.0;
    long double weightedBeta = 0.0L;
    long double weightedBetaTime = 0.0L;
    for (std::size_t i = 0; i < input.size(); ++i)
      {
        const double ux = input[i].properVelocity[0];
        const double uy = input[i].properVelocity[1];
        const double uz = input[i].properVelocity[2];
        const double gamma = std::sqrt(1.0 + ux * ux + uy * uy + uz * uz);
        const double betaZ = uz / gamma;
        weightedBeta += betaZ;
        weightedBetaTime +=
          static_cast<long double>(betaZ) * input[i].position[2];
      }
    const long double meanBeta = weightedBeta / input.size();
    const double referenceTime = static_cast<double>(
      (static_cast<long double>(referenceZ / aprl::SI::c) +
       weightedBetaTime / input.size()) / meanBeta);
    aprl::BoostFrameTransform probe;
    probe.setOriginsFromGamma(boostGamma, aprl::SI::c, 0.0, 0.0, 0.0);
    double expected = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < input.size(); ++i)
      {
        aprl::FieldVector<double> lab(0.0);
        for (unsigned int axis = 0; axis < 3; ++axis)
          lab[axis] = input[i].properVelocity[axis];
        aprl::FieldVector<double> box(0.0);
        probe.properVelocityLabToBox(lab, box);
        const double betaHeld = probe.heldLongitudinalBetaZ(box);
        const double margin = 2.0 * dz /
          (boostGamma * aprl::SI::c * (probe.beta() + betaHeld));
        expected = std::min(expected,
          input[i].position[2] - referenceTime - margin);
      }
    check(std::abs(actual - expected) < 1.0e-19,
      "Held frame origin double-counted the placed laboratory centre");
  }

  void runCase(const std::string& directory, const std::string& input,
               const std::vector<InputRecord>& expected,
               bool frozen, bool feedback, bool initialField,
               bool automaticZ = false,
               bool holdTransverse = false,
               bool activationAtInput = false,
               bool strongTransverse = false)
  {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) makeDirectory(directory);
    MPI_Barrier(MPI_COMM_WORLD);
    const std::string card = directory + "/run.yaml";
    if (rank == 0)
      writeCard(card, input, frozen, feedback, initialField,
                automaticZ, holdTransverse, activationAtInput);
    MPI_Barrier(MPI_COMM_WORLD);
    const std::string oldDirectory = []() {
      char buffer[4096];
      check(getcwd(buffer, sizeof(buffer)) != NULL, "Cannot read cwd");
      return std::string(buffer);
    }();
    check(chdir(directory.c_str()) == 0, "Cannot enter run directory");
    std::ostringstream startupLog;
    std::streambuf* originalOutput = NULL;
    if (rank == 0 && holdTransverse && activationAtInput)
      originalOutput = std::cout.rdbuf(startupLog.rdbuf());
    try
      {
        aprl::SimulationConfig config = aprl::YamlConfigLoader::loadFile(card);
        aprl::Simulation simulation(config, MPI_COMM_WORLD);
        simulation.solve();
      }
    catch (...)
      {
        if (originalOutput)
          std::cout.rdbuf(originalOutput);
        if (chdir(oldDirectory.c_str()) != 0)
          std::cerr << "Cannot restore cwd after scan failure\n";
        throw;
      }
    if (originalOutput)
      {
        std::cout.rdbuf(originalOutput);
        std::cout << startupLog.str();
        checkIntegratedHeldFrameOrigin(startupLog.str(), expected);
      }
    check(chdir(oldDirectory.c_str()) == 0, "Cannot restore cwd");
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank != 0)
      {
        MPI_Barrier(MPI_COMM_WORLD);
        return;
      }

    std::ifstream manifest((directory + "/run-manifest.yaml").c_str());
    check(static_cast<bool>(manifest), "Frozen-plane run manifest is missing");
    const std::string manifestText(
      (std::istreambuf_iterator<char>(manifest)),
      std::istreambuf_iterator<char>());
    const std::string startupMode = !frozen ? "front-anchored" :
      holdTransverse ? "held" : "ballistic";
    check(manifestText.find("resolved_beam_reference:\n  startup_mode: " +
      startupMode + "\n") != std::string::npos,
      "Frozen-plane run manifest has the wrong startup mode");

    if (holdTransverse)
      {
        const std::vector<DetectorRecord> upstream = readDetector(
          directory + "/detectors/upstream.h5");
        check(upstream.size() == expected.size(),
          "Held upstream plane lost particles");
        double maximumFrozenPositionError = 0.0;
        for (std::size_t i = 0; i < upstream.size(); ++i)
          {
            check(upstream[i].sourceId >= 1 &&
              upstream[i].sourceId <= expected.size(),
              "Upstream detector source ID is invalid");
            const InputRecord& inputRecord =
              expected[upstream[i].sourceId - 1];
            for (unsigned int axis = 0; axis < 2; ++axis)
              maximumFrozenPositionError = std::max(
                maximumFrozenPositionError,
                std::abs(upstream[i].position[axis] -
                  inputRecord.position[axis]));
          }
        check(maximumFrozenPositionError < 1.0e-12,
          "Held particle moved transversely before activation");
        std::cout << directory << " upstream_held_x_y_max_error_nm="
                  << maximumFrozenPositionError * 1.0e9 << '\n';
      }

    const char* names[3] = {"input", "near", "far"};
    const double planes[3] = {0.0, 12.0e-6, 20.0e-6};
    for (unsigned int plane = frozen ? 0 : 1; plane < 3; ++plane)
      {
        const std::vector<DetectorRecord> found = readDetector(
          directory + "/detectors/" + names[plane] + ".h5");
        check(found.size() == expected.size(),
          std::string(names[plane]) + " plane lost particles");
        double maxError = 0.0;
        double sumSquares = 0.0;
        for (std::size_t i = 0; i < found.size(); ++i)
          {
            check(found[i].sourceId >= 1 &&
              found[i].sourceId <= expected.size(),
              "Detector source ID is invalid");
            const InputRecord& inputRecord =
              expected[found[i].sourceId - 1];
            for (unsigned int axis = 0; axis < 2; ++axis)
              {
                const double ballistic = inputRecord.position[axis] +
                  inputRecord.properVelocity[axis] /
                  inputRecord.properVelocity[2] * planes[plane];
                const double error = found[i].position[axis] - ballistic;
                maxError = std::max(maxError, std::abs(error));
                sumSquares += error * error;
              }
          }
        const double rms = std::sqrt(sumSquares / (2 * found.size()));
        std::cout << directory << " " << names[plane]
                  << " ballistic_rms_um=" << rms * 1.0e6
                  << " ballistic_max_um=" << maxError * 1.0e6 << '\n';
        if (frozen && !feedback)
          check(maxError < 1.0e-10,
            "Frozen no-feedback orbit differs excessively from input-plane ballistic motion");
        else
          check(maxError < 2.0e-6,
            "Self-field case has a pathological transverse orbit deviation");
      }
    const FieldSummary field = readFieldEnergy(
      directory + "/ledger/energy-ledger.h5");
    double peakElectric = 0.0;
    double peakLightScaledMagnetic = 0.0;
    readFieldAmplitude(directory + "/detectors/self-field.h5",
      peakElectric, peakLightScaledMagnetic);
    if (initialField)
      check(field.initial > 0.0, "Initial self-field is absent");
    else
      check(field.initial == 0.0,
        "Initial self-field was not disabled");
    const double growth = field.initial > 0.0 ?
      field.maximum / field.initial : 0.0;
    std::cout << directory << " field_initial_J=" << field.initial
              << " field_max_J=" << field.maximum
              << " field_final_J=" << field.final;
    if (initialField)
      std::cout << " peak_to_initial=" << growth
                << " max_step_increase_to_initial="
                << field.maximumStepIncrease / field.initial;
    else
      std::cout << " peak_to_initial=undefined(initial_EB_zero)";
    std::cout << " peak_E_V_per_m=" << peakElectric
              << " peak_cB_V_per_m=" << peakLightScaledMagnetic << '\n';
    if (initialField)
      check(std::isfinite(growth) && growth < 2.0,
        "Free-drift self-field energy grew pathologically");
    else
      check(std::isfinite(field.maximum) && field.maximum < 1.0e-9,
        "No-initial-field startup produced pathological field energy");
    if (holdTransverse && initialField)
      check(growth < (strongTransverse ? 1.25 : 1.05) &&
            peakElectric < (strongTransverse ? 5.0e8 : 1.0e8) &&
            peakLightScaledMagnetic <
              (strongTransverse ? 5.0e8 : 1.0e8) &&
            field.maximumStepIncrease <
              (strongTransverse ? 0.15 : 0.05) * field.initial,
        "Held-transverse release produced a self-field spike");
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  try
    {
      check(argc == 2, "Usage: frozen_plane_scan <work-directory>");
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      checkHeldReleaseCurrent();
      checkHighBoostPlaneGeometry();
      const std::string work = argv[1];
      if (rank == 0) makeDirectory(work);
      MPI_Barrier(MPI_COMM_WORLD);
      const std::string input = work + "/particles.h5";
      const std::vector<InputRecord> records = rank == 0 ?
        makeInput(input) : std::vector<InputRecord>();
      MPI_Barrier(MPI_COMM_WORLD);
      runCase(work + "/frozen-initial-field", input, records,
        true, true, true);
      runCase(work + "/frozen-no-initial-field", input, records,
        true, true, false);
      runCase(work + "/unfrozen-initial-field", input, records,
        false, true, true);
      runCase(work + "/frozen-no-feedback-control", input, records,
        true, false, true);
      runCase(work + "/frozen-coupled-auto-z", input, records,
        true, true, true, true);
      runCase(work + "/ballistic-at-input", input, records,
        true, true, true, false, false, true);
      runCase(work + "/held-initial-field", input, records,
        true, true, true, false, true, true);
      runCase(work + "/held-no-initial-field", input, records,
        true, true, false, false, true, true);
      runCase(work + "/held-no-feedback-control", input, records,
        true, false, true, false, true, true);
      const std::string strongInput = work + "/particles-strong-transverse.h5";
      const std::vector<InputRecord> strongRecords = rank == 0 ?
        makeInput(strongInput, 0.4) : std::vector<InputRecord>();
      MPI_Barrier(MPI_COMM_WORLD);
      runCase(work + "/held-strong-transverse", strongInput,
        strongRecords, true, true, true, false, true, true, true);
      runCase(work + "/held-strong-transverse-no-feedback", strongInput,
        strongRecords, true, false, true, false, true, true, true);
      MPI_Finalize();
      return 0;
    }
  catch (const std::exception& error)
    {
      std::cerr << "Frozen-plane scan failed: " << error.what() << '\n';
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }
}
