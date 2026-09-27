#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>

namespace
{
  const double pi = 3.141592653589793238462643383279502884;
  const double speedOfLight = 299792458.0;

  struct ParticleRecord
  {
    double position[3];
    double properVelocity[3];
    std::uint64_t sourceId;
    double macroWeight;
  };

  struct Options
  {
    std::string output;
    double electrons;
    double requestedLength;
    double resonantWavelength;
    double gamma;
    double sigmaGamma;
    double sigmaX;
    double sigmaY;
    double emittanceNX;
    double emittanceNY;
    double twissAlphaX;
    double twissAlphaY;
    unsigned int harmonics;
    std::uint64_t beamletsPerWavelength;
    std::uint64_t seed;
    std::uint64_t chunkRecords;
    unsigned int compression;
    bool overwrite;
    bool selfTest;

    Options()
      : electrons(std::numeric_limits<double>::quiet_NaN()),
        requestedLength(std::numeric_limits<double>::quiet_NaN()),
        resonantWavelength(std::numeric_limits<double>::quiet_NaN()),
        gamma(std::numeric_limits<double>::quiet_NaN()), sigmaGamma(0.0),
        sigmaX(std::numeric_limits<double>::quiet_NaN()),
        sigmaY(std::numeric_limits<double>::quiet_NaN()),
        emittanceNX(0.0), emittanceNY(0.0), twissAlphaX(0.0),
        twissAlphaY(0.0), harmonics(3), beamletsPerWavelength(64),
        seed(1), chunkRecords(65536), compression(0), overwrite(false),
        selfTest(false)
    {}
  };

  class DeterministicRandom
  {
  public:
    explicit DeterministicRandom(std::uint64_t seed) : state_(seed) {}

    std::uint64_t integer()
    {
      state_ += UINT64_C(0x9e3779b97f4a7c15);
      std::uint64_t z = state_;
      z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
      z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
      return z ^ (z >> 31);
    }

    double uniformOpen()
    {
      const std::uint64_t mantissa = integer() >> 11;
      return (static_cast<double>(mantissa) + 0.5) /
        9007199254740992.0;
    }

    double normal()
    {
      const double radius = std::sqrt(-2.0 * std::log(uniformOpen()));
      const double angle = 2.0 * pi * uniformOpen();
      return radius * std::cos(angle);
    }

  private:
    std::uint64_t state_;
  };

  void requireHandle(hid_t handle, const std::string& message)
  {
    if (handle < 0) throw std::runtime_error(message);
  }

  void requireStatus(herr_t status, const std::string& message)
  {
    if (status < 0) throw std::runtime_error(message);
  }

  double parseDouble(const char* text, const std::string& option)
  {
    errno = 0;
    char* end = NULL;
    const double value = std::strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !std::isfinite(value))
      throw std::invalid_argument("Invalid value for " + option + ": " + text);
    return value;
  }

  std::uint64_t parseUnsigned(const char* text, const std::string& option)
  {
    if (text[0] == '-')
      throw std::invalid_argument("Invalid value for " + option + ": " + text);
    errno = 0;
    char* end = NULL;
    const unsigned long long value = std::strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0')
      throw std::invalid_argument("Invalid value for " + option + ": " + text);
    return static_cast<std::uint64_t>(value);
  }

  void usage(std::ostream& output)
  {
    output <<
      "Usage:\n"
      "  sase_bunch_generator --output FILE --electrons N --length-m L\n"
      "    --resonant-wavelength-m LAMBDA --gamma GAMMA\n"
      "    --sigma-x-m SX --sigma-y-m SY [options]\n\n"
      "Required options:\n"
      "  --output FILE                    APRL HDF5-v2 snapshot\n"
      "  --electrons N                    represented physical electrons\n"
      "  --length-m L                     requested top-hat bunch length\n"
      "  --resonant-wavelength-m LAMBDA   wavelength used for FEL phases\n"
      "  --gamma GAMMA                    mean laboratory Lorentz factor\n"
      "  --sigma-x-m SX                   rms Gaussian x size\n"
      "  --sigma-y-m SY                   rms Gaussian y size\n\n"
      "Optional phase-space controls:\n"
      "  --sigma-gamma VALUE              rms uncorrelated gamma spread [0]\n"
      "  --emittance-n-x-m VALUE          normalized rms x emittance [0]\n"
      "  --emittance-n-y-m VALUE          normalized rms y emittance [0]\n"
      "  --twiss-alpha-x VALUE            x Twiss alpha at snapshot [0]\n"
      "  --twiss-alpha-y VALUE            y Twiss alpha at snapshot [0]\n\n"
      "Shot-noise controls:\n"
      "  --harmonics M                    load correct noise through M [3]\n"
      "  --beamlets-per-wavelength B      independent beamlets per cell [64]\n"
      "  --seed SEED                      deterministic 64-bit seed [1]\n\n"
      "I/O controls:\n"
      "  --chunk-records N                HDF5 streaming chunk [65536]\n"
      "  --compression LEVEL              gzip level 0-9 [0]\n"
      "  --overwrite                      replace an existing output\n"
      "  --self-test                      run algorithm checks, write no file\n"
      "  --help                           show this message\n";
  }

  Options parseOptions(int argc, char** argv)
  {
    Options options;
    for (int i = 1; i < argc; ++i)
      {
        const std::string name(argv[i]);
        if (name == "--help")
          {
            usage(std::cout);
            std::exit(EXIT_SUCCESS);
          }
        if (name == "--overwrite")
          {
            options.overwrite = true;
            continue;
          }
        if (name == "--self-test")
          {
            options.selfTest = true;
            continue;
          }
        if (i + 1 >= argc)
          throw std::invalid_argument("Missing value after " + name);
        const char* value = argv[++i];
        if (name == "--output") options.output = value;
        else if (name == "--electrons")
          options.electrons = parseDouble(value, name);
        else if (name == "--length-m")
          options.requestedLength = parseDouble(value, name);
        else if (name == "--resonant-wavelength-m")
          options.resonantWavelength = parseDouble(value, name);
        else if (name == "--gamma")
          options.gamma = parseDouble(value, name);
        else if (name == "--sigma-gamma")
          options.sigmaGamma = parseDouble(value, name);
        else if (name == "--sigma-x-m")
          options.sigmaX = parseDouble(value, name);
        else if (name == "--sigma-y-m")
          options.sigmaY = parseDouble(value, name);
        else if (name == "--emittance-n-x-m")
          options.emittanceNX = parseDouble(value, name);
        else if (name == "--emittance-n-y-m")
          options.emittanceNY = parseDouble(value, name);
        else if (name == "--twiss-alpha-x")
          options.twissAlphaX = parseDouble(value, name);
        else if (name == "--twiss-alpha-y")
          options.twissAlphaY = parseDouble(value, name);
        else if (name == "--harmonics")
          options.harmonics = static_cast<unsigned int>(
            parseUnsigned(value, name));
        else if (name == "--beamlets-per-wavelength")
          options.beamletsPerWavelength = parseUnsigned(value, name);
        else if (name == "--seed")
          options.seed = parseUnsigned(value, name);
        else if (name == "--chunk-records")
          options.chunkRecords = parseUnsigned(value, name);
        else if (name == "--compression")
          options.compression = static_cast<unsigned int>(
            parseUnsigned(value, name));
        else
          throw std::invalid_argument("Unknown option: " + name);
      }
    return options;
  }

  void validateOptions(const Options& options)
  {
    if (options.selfTest) return;
    if (options.output.empty())
      throw std::invalid_argument("--output is required");
    if (!(options.electrons > 0.0))
      throw std::invalid_argument("--electrons must be positive");
    if (!(options.requestedLength > 0.0))
      throw std::invalid_argument("--length-m must be positive");
    if (!(options.resonantWavelength > 0.0))
      throw std::invalid_argument(
        "--resonant-wavelength-m must be positive");
    if (!(options.gamma > 1.0))
      throw std::invalid_argument("--gamma must be greater than one");
    if (!(options.sigmaGamma >= 0.0) || !(options.sigmaX > 0.0) ||
        !(options.sigmaY > 0.0) || !(options.emittanceNX >= 0.0) ||
        !(options.emittanceNY >= 0.0))
      throw std::invalid_argument(
        "Beam sizes must be positive and spreads/emittances nonnegative");
    if (options.harmonics == 0 || options.harmonics > 64)
      throw std::invalid_argument("--harmonics must be in [1,64]");
    if (options.beamletsPerWavelength == 0)
      throw std::invalid_argument(
        "--beamlets-per-wavelength must be positive");
    if (options.chunkRecords == 0)
      throw std::invalid_argument("--chunk-records must be positive");
    if (options.compression > 9)
      throw std::invalid_argument("--compression must be in [0,9]");
  }

  double wrapPhase(double phase)
  {
    phase = std::fmod(phase, 2.0 * pi);
    if (phase < 0.0) phase += 2.0 * pi;
    return phase;
  }

  std::vector<double> noisyPhases(unsigned int harmonics,
      std::uint64_t beamletIndex, std::uint64_t beamletCount,
      double electronsPerBeamlet, DeterministicRandom& random)
  {
    std::vector<double> a(harmonics + 1, 0.0);
    std::vector<double> b(harmonics + 1, 0.0);
    for (unsigned int m = 1; m <= harmonics; ++m)
      {
        const double sigma = std::sqrt(2.0 /
          (static_cast<double>(m) * static_cast<double>(m) *
           electronsPerBeamlet));
        a[m] = sigma * random.normal();
        b[m] = sigma * random.normal();
      }

    // Keep the highest loaded noise harmonic strictly below the quiet-start
    // cancellation order.  This also suppresses the adjacent unphysical
    // harmonic that can couple back into the requested highest harmonic.
    const unsigned int quietOrder = harmonics + 1U;
    const std::uint64_t members = 2U * quietOrder;
    const double base = (static_cast<double>(beamletIndex) + 0.5) * pi /
      (static_cast<double>(quietOrder) *
       static_cast<double>(beamletCount));
    std::vector<double> phases(static_cast<std::size_t>(members));
    for (std::uint64_t j = 0; j < members; ++j)
      {
        const double quiet = base + static_cast<double>(j) * pi /
          static_cast<double>(quietOrder);
        double displacement = 0.0;
        for (unsigned int m = 1; m <= harmonics; ++m)
          displacement += a[m] * std::cos(m * quiet) +
                          b[m] * std::sin(m * quiet);
        phases[static_cast<std::size_t>(j)] =
          wrapPhase(quiet + displacement);
      }
    return phases;
  }

  hid_t memoryRecordType()
  {
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(ParticleRecord));
    requireHandle(type, "Cannot create HDF5 memory record type");
    hsize_t dimensions[1] = {3};
    hid_t vector = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dimensions);
    requireHandle(vector, "Cannot create HDF5 memory vector type");
    requireStatus(H5Tinsert(type, "position_m",
      HOFFSET(ParticleRecord, position), vector),
      "Cannot add position_m to HDF5 memory record");
    requireStatus(H5Tinsert(type, "proper_velocity",
      HOFFSET(ParticleRecord, properVelocity), vector),
      "Cannot add proper_velocity to HDF5 memory record");
    requireStatus(H5Tinsert(type, "source_id",
      HOFFSET(ParticleRecord, sourceId), H5T_NATIVE_UINT64),
      "Cannot add source_id to HDF5 memory record");
    requireStatus(H5Tinsert(type, "macro_weight",
      HOFFSET(ParticleRecord, macroWeight), H5T_NATIVE_DOUBLE),
      "Cannot add macro_weight to HDF5 memory record");
    H5Tclose(vector);
    return type;
  }

  hid_t fileRecordType()
  {
    const std::size_t f64 = 8;
    hid_t type = H5Tcreate(H5T_COMPOUND, 8 * f64);
    requireHandle(type, "Cannot create HDF5 file record type");
    hsize_t dimensions[1] = {3};
    hid_t vector = H5Tarray_create2(H5T_IEEE_F64LE, 1, dimensions);
    requireHandle(vector, "Cannot create HDF5 file vector type");
    requireStatus(H5Tinsert(type, "position_m", 0, vector),
      "Cannot add position_m to HDF5 file record");
    requireStatus(H5Tinsert(type, "proper_velocity", 3 * f64, vector),
      "Cannot add proper_velocity to HDF5 file record");
    requireStatus(H5Tinsert(type, "source_id", 6 * f64, H5T_STD_U64LE),
      "Cannot add source_id to HDF5 file record");
    requireStatus(H5Tinsert(type, "macro_weight", 7 * f64,
      H5T_IEEE_F64LE), "Cannot add macro_weight to HDF5 file record");
    H5Tclose(vector);
    return type;
  }

  void writeIntAttribute(hid_t object, const char* name, int value)
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create integer attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_STD_I32LE, space,
      H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute ") + name);
    requireStatus(H5Awrite(attribute, H5T_NATIVE_INT, &value),
      std::string("Cannot write attribute ") + name);
    H5Aclose(attribute);
  }

  void writeUInt64Attribute(hid_t object, const char* name,
                            std::uint64_t value)
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create uint64 attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_STD_U64LE, space,
      H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute ") + name);
    requireStatus(H5Awrite(attribute, H5T_NATIVE_UINT64, &value),
      std::string("Cannot write attribute ") + name);
    H5Aclose(attribute);
  }

  void writeDoubleAttribute(hid_t object, const char* name, double value)
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create floating attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_IEEE_F64LE, space,
      H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute ") + name);
    requireStatus(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &value),
      std::string("Cannot write attribute ") + name);
    H5Aclose(attribute);
  }

  void writeDoubleVectorAttribute(hid_t object, const char* name,
                                  const std::vector<double>& values)
  {
    const hsize_t dimensions[1] = {
      static_cast<hsize_t>(values.size())};
    hid_t space = H5Screate_simple(1, dimensions, NULL);
    requireHandle(space, "Cannot create vector attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_IEEE_F64LE, space,
      H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute ") + name);
    requireStatus(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &values[0]),
      std::string("Cannot write attribute ") + name);
    H5Aclose(attribute);
  }

  void writeStringAttribute(hid_t object, const char* name,
                            const std::string& value)
  {
    hid_t type = H5Tcopy(H5T_C_S1);
    requireHandle(type, "Cannot create HDF5 string type");
    requireStatus(H5Tset_size(type, value.size() + 1),
      "Cannot size HDF5 string type");
    requireStatus(H5Tset_strpad(type, H5T_STR_NULLTERM),
      "Cannot set HDF5 string padding");
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create string attribute space");
    hid_t attribute = H5Acreate2(object, name, type, space,
      H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute ") + name);
    requireStatus(H5Awrite(attribute, type, value.c_str()),
      std::string("Cannot write attribute ") + name);
    H5Aclose(attribute);
    H5Tclose(type);
  }

  void writeBuffer(hid_t dataset, hid_t memoryType,
      const std::vector<ParticleRecord>& buffer, std::uint64_t offset)
  {
    if (buffer.empty()) return;
    hid_t fileSpace = H5Dget_space(dataset);
    requireHandle(fileSpace, "Cannot inspect particle output dataspace");
    const hsize_t start[1] = {static_cast<hsize_t>(offset)};
    const hsize_t count[1] = {static_cast<hsize_t>(buffer.size())};
    requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
      start, NULL, count, NULL), "Cannot select particle output hyperslab");
    hid_t memorySpace = H5Screate_simple(1, count, NULL);
    requireHandle(memorySpace, "Cannot create particle memory dataspace");
    const herr_t status = H5Dwrite(dataset, memoryType, memorySpace,
      fileSpace, H5P_DEFAULT, &buffer[0]);
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
    requireStatus(status, "Cannot write particle output records");
  }

  void sampleTransverse(const Options& options,
      DeterministicRandom& random, double& x, double& y,
      double& ux, double& uy)
  {
    const double x1 = random.normal();
    const double x2 = random.normal();
    const double y1 = random.normal();
    const double y2 = random.normal();
    x = options.sigmaX * x1;
    y = options.sigmaY * y1;
    ux = options.emittanceNX == 0.0 ? 0.0 :
      options.emittanceNX / options.sigmaX *
      (-options.twissAlphaX * x1 + x2);
    uy = options.emittanceNY == 0.0 ? 0.0 :
      options.emittanceNY / options.sigmaY *
      (-options.twissAlphaY * y1 + y2);
  }

  double sampleGamma(const Options& options, double ux, double uy,
                     DeterministicRandom& random)
  {
    const double minimum = std::sqrt(1.0 + ux * ux + uy * uy);
    for (unsigned int attempt = 0; attempt < 1000; ++attempt)
      {
        const double value = options.gamma +
          options.sigmaGamma * random.normal();
        if (value > minimum) return value;
      }
    throw std::runtime_error(
      "Cannot sample a physical gamma; reduce spreads or increase --gamma");
  }

  void validateWrittenFile(const std::string& filename,
                           std::uint64_t expectedRecords)
  {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot reopen generated output");
    hid_t group = H5Gopen2(file, "/particles", H5P_DEFAULT);
    requireHandle(group, "Generated output is missing /particles");
    hid_t attribute = H5Aopen(group, "format_version", H5P_DEFAULT);
    requireHandle(attribute, "Generated output is missing format_version");
    int version = 0;
    requireStatus(H5Aread(attribute, H5T_NATIVE_INT, &version),
      "Cannot read generated format_version");
    H5Aclose(attribute);
    if (version != 2)
      throw std::runtime_error("Generated output format_version is not 2");
    attribute = H5Aopen(group, "complete", H5P_DEFAULT);
    requireHandle(attribute, "Generated output is missing complete");
    int complete = 0;
    requireStatus(H5Aread(attribute, H5T_NATIVE_INT, &complete),
      "Cannot read generated complete marker");
    H5Aclose(attribute);
    if (complete != 1)
      throw std::runtime_error("Generated output is incomplete");
    hid_t dataset = H5Dopen2(group, "records", H5P_DEFAULT);
    requireHandle(dataset, "Generated output is missing records");
    hid_t space = H5Dget_space(dataset);
    requireHandle(space, "Cannot inspect generated records");
    hsize_t dimensions[1] = {0};
    if (H5Sget_simple_extent_ndims(space) != 1)
      throw std::runtime_error("Generated records are not one-dimensional");
    H5Sget_simple_extent_dims(space, dimensions, NULL);
    H5Sclose(space);
    H5Dclose(dataset);
    H5Gclose(group);
    H5Fclose(file);
    if (dimensions[0] != static_cast<hsize_t>(expectedRecords))
      throw std::runtime_error("Generated record count changed on reopen");
  }

  void generate(const Options& options)
  {
    const double rawCells = options.requestedLength /
      options.resonantWavelength;
    const std::uint64_t wavelengthCells = static_cast<std::uint64_t>(
      std::floor(rawCells + 0.5));
    if (wavelengthCells == 0)
      throw std::invalid_argument(
        "Bunch length must contain at least one resonant wavelength");
    const double effectiveLength = wavelengthCells *
      options.resonantWavelength;
    const std::uint64_t members = 2U * (options.harmonics + 1U);
    if (options.beamletsPerWavelength >
        std::numeric_limits<std::uint64_t>::max() / members ||
        wavelengthCells > std::numeric_limits<std::uint64_t>::max() /
          (options.beamletsPerWavelength * members))
      throw std::overflow_error("Requested macroparticle count overflows uint64");
    const std::uint64_t totalRecords = wavelengthCells *
      options.beamletsPerWavelength * members;
    if (totalRecords > static_cast<std::uint64_t>(
        std::numeric_limits<hsize_t>::max()))
      throw std::overflow_error("Requested macroparticle count exceeds HDF5");

    const double electronsPerCell = options.electrons /
      static_cast<double>(wavelengthCells);
    const double electronsPerBeamlet = electronsPerCell /
      static_cast<double>(options.beamletsPerWavelength);
    if (!(electronsPerBeamlet > 1.0))
      throw std::invalid_argument(
        "Each beamlet must represent more than one physical electron");
    if (electronsPerBeamlet < 100.0)
      std::cerr << "warning: only " << electronsPerBeamlet
                << " physical electrons per beamlet; Fawley's small-phase "
                   "approximation is marginal\n";

    const unsigned int createFlag = options.overwrite ?
      H5F_ACC_TRUNC : H5F_ACC_EXCL;
    hid_t file = H5Fcreate(options.output.c_str(), createFlag,
      H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(file, "Cannot create output file (use --overwrite only "
      "when replacement is intentional): " + options.output);
    hid_t group = H5Gcreate2(file, "/particles", H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(group, "Cannot create /particles group");

    writeIntAttribute(group, "format_version", 2);
    writeIntAttribute(group, "complete", 0);
    writeStringAttribute(group, "generator", "aprl_sase_bunch_generator");
    writeStringAttribute(group, "coordinate_frame",
      "common-laboratory-time-relative-snapshot");
    writeStringAttribute(group, "longitudinal_distribution",
      "uniform-top-hat-in-integer-resonant-wavelength-cells");
    writeStringAttribute(group, "transverse_distribution",
      "gaussian-beamlets-with-twiss-correlated-proper-velocity");
    writeStringAttribute(group, "shot_noise_model",
      "Fawley-2002-multibeamlet-multiharmonic-quiet-start");
    writeStringAttribute(group, "macro_weight_definition",
      "unit-relative-weight-normalized-by-APRL-to-beam.input.electrons");
    writeStringAttribute(group, "position_unit", "m");
    writeStringAttribute(group, "proper_velocity_unit", "gamma*v/c");
    writeStringAttribute(group, "source_id_definition",
      "one-based-cell-major-beamlet-major-member-major-index");
    writeDoubleAttribute(group, "represented_electrons", options.electrons);
    writeDoubleAttribute(group, "requested_length_m",
      options.requestedLength);
    writeDoubleAttribute(group, "effective_length_m", effectiveLength);
    writeDoubleAttribute(group, "resonant_wavelength_m",
      options.resonantWavelength);
    writeDoubleAttribute(group, "electrons_per_wavelength_cell",
      electronsPerCell);
    writeDoubleAttribute(group, "electrons_per_beamlet",
      electronsPerBeamlet);
    writeDoubleAttribute(group, "gamma_mean", options.gamma);
    writeDoubleAttribute(group, "sigma_gamma", options.sigmaGamma);
    writeDoubleAttribute(group, "sigma_x_m", options.sigmaX);
    writeDoubleAttribute(group, "sigma_y_m", options.sigmaY);
    writeDoubleAttribute(group, "emittance_n_x_m", options.emittanceNX);
    writeDoubleAttribute(group, "emittance_n_y_m", options.emittanceNY);
    writeDoubleAttribute(group, "twiss_alpha_x", options.twissAlphaX);
    writeDoubleAttribute(group, "twiss_alpha_y", options.twissAlphaY);
    writeUInt64Attribute(group, "wavelength_cells", wavelengthCells);
    writeUInt64Attribute(group, "beamlets_per_wavelength",
      options.beamletsPerWavelength);
    writeUInt64Attribute(group, "macroparticles_per_beamlet", members);
    writeIntAttribute(group, "quiet_start_harmonic_order",
      static_cast<int>(options.harmonics + 1U));
    writeUInt64Attribute(group, "macroparticles", totalRecords);
    writeUInt64Attribute(group, "random_seed", options.seed);
    writeIntAttribute(group, "shot_noise_harmonics",
      static_cast<int>(options.harmonics));

    const hsize_t dimensions[1] = {static_cast<hsize_t>(totalRecords)};
    hid_t dataSpace = H5Screate_simple(1, dimensions, NULL);
    requireHandle(dataSpace, "Cannot create particle dataset space");
    hid_t fileType = fileRecordType();
    hid_t creation = H5Pcreate(H5P_DATASET_CREATE);
    requireHandle(creation, "Cannot create particle dataset properties");
    const hsize_t chunk[1] = {static_cast<hsize_t>(
      std::min(options.chunkRecords, totalRecords))};
    requireStatus(H5Pset_chunk(creation, 1, chunk),
      "Cannot set particle dataset chunk");
    if (options.compression > 0)
      requireStatus(H5Pset_deflate(creation, options.compression),
        "Cannot enable particle dataset compression");
    hid_t dataset = H5Dcreate2(group, "records", fileType, dataSpace,
      H5P_DEFAULT, creation, H5P_DEFAULT);
    H5Pclose(creation);
    H5Sclose(dataSpace);
    H5Tclose(fileType);
    requireHandle(dataset, "Cannot create /particles/records");
    hid_t memoryType = memoryRecordType();

    std::vector<ParticleRecord> buffer;
    buffer.reserve(static_cast<std::size_t>(
      std::min(options.chunkRecords, totalRecords)));
    std::vector<double> globalReal(options.harmonics, 0.0);
    std::vector<double> globalImag(options.harmonics, 0.0);
    std::vector<double> meanScaledLocalPower(options.harmonics, 0.0);
    DeterministicRandom random(options.seed);
    std::uint64_t written = 0;
    std::uint64_t sourceId = 1;
    double betaZSum = 0.0;

    for (std::uint64_t cell = 0; cell < wavelengthCells; ++cell)
      {
        std::vector<double> cellReal(options.harmonics, 0.0);
        std::vector<double> cellImag(options.harmonics, 0.0);
        for (std::uint64_t beamlet = 0;
             beamlet < options.beamletsPerWavelength; ++beamlet)
          {
            double x = 0.0, y = 0.0, ux = 0.0, uy = 0.0;
            sampleTransverse(options, random, x, y, ux, uy);
            const double gamma = sampleGamma(options, ux, uy, random);
            const double uz = std::sqrt(
              gamma * gamma - 1.0 - ux * ux - uy * uy);
            betaZSum += static_cast<double>(members) * uz / gamma;
            const std::vector<double> phases = noisyPhases(
              options.harmonics, beamlet,
              options.beamletsPerWavelength,
              electronsPerBeamlet, random);
            for (std::uint64_t member = 0; member < members; ++member)
              {
                ParticleRecord record = {};
                const double phase = phases[static_cast<std::size_t>(member)];
                record.position[0] = x;
                record.position[1] = y;
                record.position[2] = -0.5 * effectiveLength +
                  (static_cast<double>(cell) + phase / (2.0 * pi)) *
                  options.resonantWavelength;
                record.properVelocity[0] = ux;
                record.properVelocity[1] = uy;
                record.properVelocity[2] = uz;
                record.sourceId = sourceId++;
                record.macroWeight = 1.0;
                buffer.push_back(record);
                for (unsigned int m = 1; m <= options.harmonics; ++m)
                  {
                    const double real = std::cos(m * phase);
                    const double imag = std::sin(m * phase);
                    cellReal[m - 1] += real;
                    cellImag[m - 1] += imag;
                    globalReal[m - 1] += real;
                    globalImag[m - 1] += imag;
                  }
                if (buffer.size() == buffer.capacity())
                  {
                    writeBuffer(dataset, memoryType, buffer, written);
                    written += static_cast<std::uint64_t>(buffer.size());
                    buffer.clear();
                  }
              }
          }
        const double recordsPerCell = static_cast<double>(
          options.beamletsPerWavelength * members);
        for (unsigned int m = 0; m < options.harmonics; ++m)
          {
            cellReal[m] /= recordsPerCell;
            cellImag[m] /= recordsPerCell;
            meanScaledLocalPower[m] += electronsPerCell *
              (cellReal[m] * cellReal[m] +
               cellImag[m] * cellImag[m]);
          }
      }
    writeBuffer(dataset, memoryType, buffer, written);
    written += static_cast<std::uint64_t>(buffer.size());
    H5Tclose(memoryType);
    H5Dclose(dataset);
    if (written != totalRecords || sourceId != totalRecords + 1)
      throw std::runtime_error("Internal generated record count mismatch");

    std::vector<double> observedGlobal(options.harmonics, 0.0);
    std::vector<double> expectedGlobal(options.harmonics,
      1.0 / std::sqrt(options.electrons));
    for (unsigned int m = 0; m < options.harmonics; ++m)
      {
        globalReal[m] /= static_cast<double>(totalRecords);
        globalImag[m] /= static_cast<double>(totalRecords);
        observedGlobal[m] = std::sqrt(globalReal[m] * globalReal[m] +
                                      globalImag[m] * globalImag[m]);
        meanScaledLocalPower[m] /= static_cast<double>(wavelengthCells);
      }
    const double meanBetaZ = betaZSum / static_cast<double>(totalRecords);
    const double meanCurrent = options.electrons * 1.602176634e-19 *
      speedOfLight * meanBetaZ / effectiveLength;
    writeDoubleAttribute(group, "mean_beta_z", meanBetaZ);
    writeDoubleAttribute(group, "mean_current_A", meanCurrent);
    writeDoubleVectorAttribute(group, "observed_global_bunching_abs",
      observedGlobal);
    writeDoubleVectorAttribute(group, "expected_global_bunching_rms",
      expectedGlobal);
    writeDoubleVectorAttribute(group,
      "mean_electrons_per_cell_times_local_bunching_squared",
      meanScaledLocalPower);
    H5Adelete(group, "complete");
    writeIntAttribute(group, "complete", 1);
    requireStatus(H5Fflush(file, H5F_SCOPE_GLOBAL),
      "Cannot flush generated output");
    H5Gclose(group);
    H5Fclose(file);

    validateWrittenFile(options.output, totalRecords);

    std::cout << std::setprecision(12)
      << "Generated APRL HDF5-v2 SASE input\n"
      << "  output: " << options.output << "\n"
      << "  requested/effective length [m]: "
      << options.requestedLength << " / " << effectiveLength << "\n"
      << "  wavelength cells: " << wavelengthCells << "\n"
      << "  macroparticles: " << totalRecords << "\n"
      << "  represented electrons: " << options.electrons << "\n"
      << "  mean beta_z/current [A]: " << meanBetaZ << " / "
      << meanCurrent << "\n"
      << "  electrons/beamlet: " << electronsPerBeamlet << "\n";
    for (unsigned int m = 0; m < options.harmonics; ++m)
      std::cout << "  harmonic " << (m + 1)
                << ": |b_global|=" << observedGlobal[m]
                << ", expected_rms=" << expectedGlobal[m]
                << ", mean(N_cell*|b_cell|^2)="
                << meanScaledLocalPower[m] << "\n";
    std::cout << "Use beam.input.electrons=" << options.electrons
              << " in the APRL run card.\n";
  }

  void selfTest()
  {
    for (unsigned int harmonics = 1; harmonics <= 8; ++harmonics)
      {
        const unsigned int members = 2U * harmonics;
        for (unsigned int m = 1; m <= harmonics; ++m)
          {
            double real = 0.0;
            double imag = 0.0;
            for (unsigned int j = 0; j < members; ++j)
              {
                const double phase = 0.173 + j * pi / harmonics;
                real += std::cos(m * phase);
                imag += std::sin(m * phase);
              }
            if (std::sqrt(real * real + imag * imag) / members > 1.0e-13)
              throw std::runtime_error(
                "Quiet-start harmonic cancellation self-test failed");
          }
      }

    const unsigned int harmonics = 4;
    const std::uint64_t beamlets = 32;
    const std::uint64_t cells = 4096;
    const double electronsPerCell = 2.0e7;
    const double electronsPerBeamlet = electronsPerCell / beamlets;
    const std::uint64_t members = 2U * (harmonics + 1U);
    DeterministicRandom random(UINT64_C(0x4d595df4d0f33173));
    std::vector<double> scaledPower(harmonics, 0.0);
    for (std::uint64_t cell = 0; cell < cells; ++cell)
      {
        std::vector<double> real(harmonics, 0.0);
        std::vector<double> imag(harmonics, 0.0);
        for (std::uint64_t beamlet = 0; beamlet < beamlets; ++beamlet)
          {
            const std::vector<double> phases = noisyPhases(harmonics,
              beamlet, beamlets, electronsPerBeamlet, random);
            for (std::size_t j = 0; j < phases.size(); ++j)
              for (unsigned int m = 1; m <= harmonics; ++m)
                {
                  real[m - 1] += std::cos(m * phases[j]);
                  imag[m - 1] += std::sin(m * phases[j]);
                }
          }
        const double count = static_cast<double>(beamlets * members);
        for (unsigned int m = 0; m < harmonics; ++m)
          scaledPower[m] += electronsPerCell *
            (real[m] * real[m] + imag[m] * imag[m]) / (count * count);
      }
    for (unsigned int m = 0; m < harmonics; ++m)
      {
        scaledPower[m] /= static_cast<double>(cells);
        if (std::abs(scaledPower[m] - 1.0) > 0.10)
          {
            std::ostringstream message;
            message << "Shot-noise variance self-test failed at harmonic "
                    << (m + 1) << ": " << scaledPower[m];
            throw std::runtime_error(message.str());
          }
      }
    std::cout << "SASE bunch generator self-test passed";
    for (unsigned int m = 0; m < harmonics; ++m)
      std::cout << " h" << (m + 1) << "=" << scaledPower[m];
    std::cout << "\n";
  }
}

int main(int argc, char** argv)
{
  try
    {
      const Options options = parseOptions(argc, argv);
      validateOptions(options);
      if (options.selfTest) selfTest();
      else generate(options);
      return EXIT_SUCCESS;
    }
  catch (const std::exception& error)
    {
      std::cerr << "sase_bunch_generator: " << error.what() << "\n";
      return EXIT_FAILURE;
    }
}
