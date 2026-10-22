#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>

namespace
{
  const double pi = 3.141592653589793238462643383279502884;
  const double speedOfLight = 299792458.0;
  const double electronCharge = 1.602176634e-19;
  const double electronRestEnergyMeV = 0.51099895000;

  struct ParticleRecord
  {
    double planePosition[2];
    double arrivalTimeOffset;
    double properVelocity[3];
    std::uint64_t sourceId;
    double macroWeight;
  };
  static_assert(sizeof(ParticleRecord) == 64,
    "The synthetic v4 record must occupy 64 bytes");

  struct Options
  {
    std::string output;
    double length;
    double current;
    double macroWeightElectrons;
    double centralMomentumMeVc;
    double emitX;
    double emitY;
    double betaX;
    double betaY;
    double alphaX;
    double alphaY;
    double sigmaDp;
    double cutoffSigma;
    std::uint64_t seed;
    std::uint64_t chunkRecords;
    unsigned int compression;
    bool overwrite;
    bool plan;
    std::string verify;

    Options()
      : length(10.0e-6), current(50.0), macroWeightElectrons(1.0),
        centralMomentumMeVc(600.0), emitX(1.0e-9), emitY(1.0e-9),
        betaX(4.527431766), betaY(4.771248310), alphaX(0.0),
        alphaY(0.0), sigmaDp(1.3e-4), cutoffSigma(5.0),
        seed(UINT64_C(20261019)), chunkRecords(65536), compression(0),
        overwrite(false), plan(false)
    {}
  };

  class DeterministicRandom
  {
  public:
    explicit DeterministicRandom(std::uint64_t seed) : state_(seed),
      spare_(0.0), hasSpare_(false) {}

    std::uint64_t integer()
    {
      state_ += UINT64_C(0x9e3779b97f4a7c15);
      std::uint64_t value = state_;
      value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
      value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
      return value ^ (value >> 31);
    }

    double uniformOpen()
    {
      return (static_cast<double>(integer() >> 11) + 0.5) /
        9007199254740992.0;
    }

    double normal()
    {
      if (hasSpare_)
        {
          hasSpare_ = false;
          return spare_;
        }
      const double radius = std::sqrt(-2.0 * std::log(uniformOpen()));
      const double angle = 2.0 * pi * uniformOpen();
      spare_ = radius * std::sin(angle);
      hasSpare_ = true;
      return radius * std::cos(angle);
    }

    double truncatedNormal(double cutoff)
    {
      double sample = 0.0;
      do sample = normal(); while (std::abs(sample) > cutoff);
      return sample;
    }

  private:
    std::uint64_t state_;
    double spare_;
    bool hasSpare_;
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
      throw std::invalid_argument("Negative value for " + option);
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
      "Usage: uniform_lab_plane_bunch_generator --output FILE [options]\n"
      "       uniform_lab_plane_bunch_generator --plan [options]\n\n"
      "Independent APRL HDF5-v4 fixed-plane input; no Elegant required.\n"
      "Defaults reproduce a 10 um, 50 A, 600 MeV/c synthetic bunch:\n"
      "  --length-m L                    Full uniform z width [1e-5]\n"
      "  --current-A I                   Requested current [50]\n"
      "  --macro-weight-electrons W      Electrons per macro-particle [1]\n"
      "  --central-momentum-MeV-c P      Central momentum [600]\n"
      "  --emit-x-m EX                   Geometric horizontal emittance [1e-9]\n"
      "  --emit-y-m EY                   Geometric vertical emittance [1e-9]\n"
      "  --beta-x-m BX                   Horizontal Twiss beta [4.527431766]\n"
      "  --beta-y-m BY                   Vertical Twiss beta [4.771248310]\n"
      "  --alpha-x AX                    Horizontal Twiss alpha [0]\n"
      "  --alpha-y AY                    Vertical Twiss alpha [0]\n"
      "  --sigma-dp VALUE                RMS relative momentum spread [1.3e-4]\n"
      "  --cutoff-sigma VALUE            Independent Gaussian cutoff [5]\n"
      "  --seed N                        Deterministic RNG seed [20261019]\n"
      "  --chunk-records N               Streaming HDF5 chunk [65536]\n"
      "  --compression LEVEL            gzip level 0-9 [0]\n"
      "  --overwrite                    Replace output explicitly\n"
      "  --plan                         Print size/current plan only\n"
      "  --verify FILE                  Read back and validate an existing file\n"
      "  --help                         Print this message\n";
  }

  Options parseOptions(int argc, char** argv)
  {
    Options options;
    for (int index = 1; index < argc; ++index)
      {
        const std::string option(argv[index]);
        if (option == "--help")
          {
            usage(std::cout);
            std::exit(EXIT_SUCCESS);
          }
        if (option == "--overwrite")
          {
            options.overwrite = true;
            continue;
          }
        if (option == "--plan")
          {
            options.plan = true;
            continue;
          }
        if (index + 1 >= argc)
          throw std::invalid_argument("Missing value after " + option);
        const char* value = argv[++index];
        if (option == "--output") options.output = value;
        else if (option == "--verify") options.verify = value;
        else if (option == "--length-m")
          options.length = parseDouble(value, option);
        else if (option == "--current-A")
          options.current = parseDouble(value, option);
        else if (option == "--macro-weight-electrons")
          options.macroWeightElectrons = parseDouble(value, option);
        else if (option == "--central-momentum-MeV-c")
          options.centralMomentumMeVc = parseDouble(value, option);
        else if (option == "--emit-x-m")
          options.emitX = parseDouble(value, option);
        else if (option == "--emit-y-m")
          options.emitY = parseDouble(value, option);
        else if (option == "--beta-x-m")
          options.betaX = parseDouble(value, option);
        else if (option == "--beta-y-m")
          options.betaY = parseDouble(value, option);
        else if (option == "--alpha-x")
          options.alphaX = parseDouble(value, option);
        else if (option == "--alpha-y")
          options.alphaY = parseDouble(value, option);
        else if (option == "--sigma-dp")
          options.sigmaDp = parseDouble(value, option);
        else if (option == "--cutoff-sigma")
          options.cutoffSigma = parseDouble(value, option);
        else if (option == "--seed")
          options.seed = parseUnsigned(value, option);
        else if (option == "--chunk-records")
          options.chunkRecords = parseUnsigned(value, option);
        else if (option == "--compression")
          {
            const std::uint64_t level = parseUnsigned(value, option);
            if (level > 9)
              throw std::invalid_argument("--compression must be in [0,9]");
            options.compression = static_cast<unsigned int>(level);
          }
        else
          throw std::invalid_argument("Unknown option: " + option);
      }
    return options;
  }

  std::uint64_t particleCount(const Options& options, double& betaReference)
  {
    const double p = options.centralMomentumMeVc / electronRestEnergyMeV;
    betaReference = p / std::hypot(1.0, p);
    const long double target =
      static_cast<long double>(options.current) * options.length /
      (static_cast<long double>(electronCharge) * speedOfLight *
       betaReference * options.macroWeightElectrons);
    if (!(target >= 0.5L) ||
        target > static_cast<long double>(
          std::numeric_limits<long long>::max()))
      throw std::invalid_argument(
        "Current, length and macro weight imply an invalid particle count");
    return static_cast<std::uint64_t>(std::llround(target));
  }

  void validateOptions(const Options& options)
  {
    if (!options.verify.empty())
      {
        if (options.plan || !options.output.empty())
          throw std::invalid_argument(
            "--verify cannot be combined with --plan or --output");
        return;
      }
    if (!options.plan && options.output.empty())
      throw std::invalid_argument("--output is required unless --plan is set");
    if (!(options.length > 0.0) || !(options.current > 0.0) ||
        !(options.macroWeightElectrons > 0.0) ||
        !(options.centralMomentumMeVc > 0.0) ||
        !(options.emitX > 0.0) || !(options.emitY > 0.0) ||
        !(options.betaX > 0.0) || !(options.betaY > 0.0) ||
        !(options.sigmaDp >= 0.0) || !(options.cutoffSigma > 0.0) ||
        options.chunkRecords == 0 || options.chunkRecords > 1048576 ||
        options.compression > 9)
      throw std::invalid_argument("Invalid positive beam or output control");
    const double p = options.centralMomentumMeVc / electronRestEnergyMeV;
    const double transverse = p * options.sigmaDp * options.cutoffSigma;
    if (!std::isfinite(p) || !std::isfinite(transverse) ||
        options.sigmaDp * options.cutoffSigma >= 1.0)
      throw std::invalid_argument(
        "Momentum spread/cutoff must leave positive forward momentum");
    double betaReference = 0.0;
    (void)particleCount(options, betaReference);
  }

  struct Sample
  {
    double x;
    double y;
    double z;
    double xp;
    double yp;
    double dp;
    double uz;
    double ux;
    double uy;
    double betaZ;
    double crossingTime;
  };

  Sample sampleParticle(const Options& options, DeterministicRandom& random)
  {
    const double rx = random.truncatedNormal(options.cutoffSigma);
    const double rpx = random.truncatedNormal(options.cutoffSigma);
    const double ry = random.truncatedNormal(options.cutoffSigma);
    const double rpy = random.truncatedNormal(options.cutoffSigma);
    const double rdp = random.truncatedNormal(options.cutoffSigma);
    Sample sample = {};
    sample.z = options.length * (random.uniformOpen() - 0.5);
    sample.x = std::sqrt(options.emitX * options.betaX) * rx;
    sample.y = std::sqrt(options.emitY * options.betaY) * ry;
    sample.xp = std::sqrt(options.emitX / options.betaX) *
      (-options.alphaX * rx + rpx);
    sample.yp = std::sqrt(options.emitY / options.betaY) *
      (-options.alphaY * ry + rpy);
    sample.dp = options.sigmaDp * rdp;
    const double p = options.centralMomentumMeVc /
      electronRestEnergyMeV * (1.0 + sample.dp);
    const double slopes = std::sqrt(1.0 + sample.xp * sample.xp +
      sample.yp * sample.yp);
    sample.uz = p / slopes;
    sample.ux = sample.xp * sample.uz;
    sample.uy = sample.yp * sample.uz;
    sample.betaZ = sample.uz / std::hypot(1.0, p);
    sample.crossingTime = -sample.z / (speedOfLight * sample.betaZ);
    return sample;
  }

  struct Statistics
  {
    long double time;
    long double betaZ;
    long double x;
    long double y;
    long double xp;
    long double yp;
    long double x2;
    long double y2;
    long double xp2;
    long double yp2;

    Statistics() : time(0), betaZ(0), x(0), y(0), xp(0), yp(0),
      x2(0), y2(0), xp2(0), yp2(0) {}

    void add(const Sample& sample)
    {
      time += sample.crossingTime;
      betaZ += sample.betaZ;
      x += sample.x;
      y += sample.y;
      xp += sample.xp;
      yp += sample.yp;
      x2 += static_cast<long double>(sample.x) * sample.x;
      y2 += static_cast<long double>(sample.y) * sample.y;
      xp2 += static_cast<long double>(sample.xp) * sample.xp;
      yp2 += static_cast<long double>(sample.yp) * sample.yp;
    }
  };

  double rms(long double sum, long double squares, std::uint64_t count)
  {
    const long double mean = sum / count;
    return std::sqrt(static_cast<double>(std::max(
      0.0L, squares / count - mean * mean)));
  }

  hid_t memoryRecordType()
  {
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(ParticleRecord));
    requireHandle(type, "Cannot create HDF5 memory record type");
    hsize_t two[1] = {2};
    hsize_t three[1] = {3};
    hid_t xy = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, two);
    hid_t velocity = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, three);
    requireHandle(xy, "Cannot create plane coordinate type");
    requireHandle(velocity, "Cannot create velocity type");
    requireStatus(H5Tinsert(type, "plane_position_m",
      HOFFSET(ParticleRecord, planePosition), xy),
      "Cannot add plane_position_m");
    requireStatus(H5Tinsert(type, "arrival_time_offset_s",
      HOFFSET(ParticleRecord, arrivalTimeOffset), H5T_NATIVE_DOUBLE),
      "Cannot add arrival_time_offset_s");
    requireStatus(H5Tinsert(type, "proper_velocity",
      HOFFSET(ParticleRecord, properVelocity), velocity),
      "Cannot add proper_velocity");
    requireStatus(H5Tinsert(type, "source_id",
      HOFFSET(ParticleRecord, sourceId), H5T_NATIVE_UINT64),
      "Cannot add source_id");
    requireStatus(H5Tinsert(type, "macro_weight",
      HOFFSET(ParticleRecord, macroWeight), H5T_NATIVE_DOUBLE),
      "Cannot add macro_weight");
    H5Tclose(xy);
    H5Tclose(velocity);
    return type;
  }

  hid_t fileRecordType()
  {
    hid_t type = H5Tcreate(H5T_COMPOUND, 64);
    requireHandle(type, "Cannot create HDF5 file record type");
    hsize_t two[1] = {2};
    hsize_t three[1] = {3};
    hid_t xy = H5Tarray_create2(H5T_IEEE_F64LE, 1, two);
    hid_t velocity = H5Tarray_create2(H5T_IEEE_F64LE, 1, three);
    requireHandle(xy, "Cannot create file plane coordinate type");
    requireHandle(velocity, "Cannot create file velocity type");
    requireStatus(H5Tinsert(type, "plane_position_m", 0, xy),
      "Cannot add file plane_position_m");
    requireStatus(H5Tinsert(type, "arrival_time_offset_s", 16,
      H5T_IEEE_F64LE), "Cannot add file arrival_time_offset_s");
    requireStatus(H5Tinsert(type, "proper_velocity", 24, velocity),
      "Cannot add file proper_velocity");
    requireStatus(H5Tinsert(type, "source_id", 48, H5T_STD_U64LE),
      "Cannot add file source_id");
    requireStatus(H5Tinsert(type, "macro_weight", 56,
      H5T_IEEE_F64LE), "Cannot add file macro_weight");
    H5Tclose(xy);
    H5Tclose(velocity);
    return type;
  }

  void writeIntAttribute(hid_t object, const char* name, int value)
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create integer attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_STD_I32LE,
      space, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create ") + name);
    requireStatus(H5Awrite(attribute, H5T_NATIVE_INT, &value),
      std::string("Cannot write ") + name);
    H5Aclose(attribute);
  }

  void writeDoubleAttribute(hid_t object, const char* name, double value)
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create double attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_IEEE_F64LE,
      space, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create ") + name);
    requireStatus(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &value),
      std::string("Cannot write ") + name);
    H5Aclose(attribute);
  }

  void writeUInt64Attribute(hid_t object, const char* name,
                            std::uint64_t value)
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create uint64 attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_STD_U64LE,
      space, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create ") + name);
    requireStatus(H5Awrite(attribute, H5T_NATIVE_UINT64, &value),
      std::string("Cannot write ") + name);
    H5Aclose(attribute);
  }

  void writeStringAttribute(hid_t object, const char* name,
                            const std::string& value)
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create string attribute space");
    hid_t type = H5Tcopy(H5T_C_S1);
    requireHandle(type, "Cannot create string attribute type");
    requireStatus(H5Tset_size(type, value.size() + 1),
      "Cannot size string attribute type");
    requireStatus(H5Tset_strpad(type, H5T_STR_NULLTERM),
      "Cannot set string attribute padding");
    hid_t attribute = H5Acreate2(object, name, type,
      space, H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(attribute, std::string("Cannot create ") + name);
    requireStatus(H5Awrite(attribute, type, value.c_str()),
      std::string("Cannot write ") + name);
    H5Aclose(attribute);
    H5Tclose(type);
    H5Sclose(space);
  }

  void writeChunk(hid_t dataset, hid_t memoryType,
                  const std::vector<ParticleRecord>& records,
                  std::uint64_t offset)
  {
    const hsize_t start[1] = {static_cast<hsize_t>(offset)};
    const hsize_t count[1] = {static_cast<hsize_t>(records.size())};
    hid_t fileSpace = H5Dget_space(dataset);
    requireHandle(fileSpace, "Cannot obtain HDF5 file space");
    requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
      start, NULL, count, NULL), "Cannot select HDF5 output block");
    hid_t memorySpace = H5Screate_simple(1, count, NULL);
    requireHandle(memorySpace, "Cannot create HDF5 memory block");
    requireStatus(H5Dwrite(dataset, memoryType, memorySpace, fileSpace,
      H5P_DEFAULT, &records[0]), "Cannot write HDF5 output block");
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
  }

  int readIntAttribute(hid_t object, const char* name)
  {
    hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
    requireHandle(attribute, std::string("Missing attribute ") + name);
    int value = 0;
    requireStatus(H5Aread(attribute, H5T_NATIVE_INT, &value),
      std::string("Cannot read attribute ") + name);
    H5Aclose(attribute);
    return value;
  }

  double readDoubleAttribute(hid_t object, const char* name)
  {
    hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
    requireHandle(attribute, std::string("Missing attribute ") + name);
    double value = 0.0;
    requireStatus(H5Aread(attribute, H5T_NATIVE_DOUBLE, &value),
      std::string("Cannot read attribute ") + name);
    H5Aclose(attribute);
    return value;
  }

  std::uint64_t readUInt64Attribute(hid_t object, const char* name)
  {
    hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
    requireHandle(attribute, std::string("Missing attribute ") + name);
    std::uint64_t value = 0;
    requireStatus(H5Aread(attribute, H5T_NATIVE_UINT64, &value),
      std::string("Cannot read attribute ") + name);
    H5Aclose(attribute);
    return value;
  }

  void verifyFile(const std::string& filename)
  {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open synthetic bunch for verification: " +
      filename);
    hid_t group = H5Gopen2(file, "/particles", H5P_DEFAULT);
    requireHandle(group, "Missing /particles group");
    if (readIntAttribute(group, "format_version") != 4 ||
        readIntAttribute(group, "complete") != 1)
      throw std::runtime_error("Synthetic bunch is incomplete or not v4");
    const std::uint64_t expected =
      readUInt64Attribute(group, "macro_particles");
    const double electrons =
      readDoubleAttribute(group, "represented_electrons");
    const double weight =
      readDoubleAttribute(group, "macro_weight_electrons");
    const double length =
      readDoubleAttribute(group, "requested_full_width_m");
    const double current =
      readDoubleAttribute(group, "realized_mean_current_A");
    if (expected == 0 || !(electrons > 0.0) || !(weight > 0.0) ||
        !(length > 0.0) || !(current > 0.0) ||
        std::abs(electrons / (static_cast<double>(expected) * weight) - 1.0) >
          1.0e-12)
      throw std::runtime_error("Synthetic bunch metadata is inconsistent");
    hid_t dataset = H5Dopen2(group, "records", H5P_DEFAULT);
    requireHandle(dataset, "Missing /particles/records dataset");
    hid_t fileSpace = H5Dget_space(dataset);
    requireHandle(fileSpace, "Cannot inspect synthetic bunch records");
    hsize_t dimensions[1] = {0};
    if (H5Sget_simple_extent_ndims(fileSpace) != 1 ||
        H5Sget_simple_extent_dims(fileSpace, dimensions, NULL) < 0 ||
        dimensions[0] != expected)
      throw std::runtime_error("Synthetic bunch record count mismatch");
    H5Sclose(fileSpace);
    hid_t memoryType = memoryRecordType();
    std::vector<ParticleRecord> records(65536);
    long double sumOffset = 0.0L;
    long double sumBetaZ = 0.0L;
    double minimumZ = std::numeric_limits<double>::infinity();
    double maximumZ = -std::numeric_limits<double>::infinity();
    for (std::uint64_t offset = 0; offset < expected;)
      {
        const std::uint64_t rows = std::min<std::uint64_t>(
          records.size(), expected - offset);
        const hsize_t start[1] = {static_cast<hsize_t>(offset)};
        const hsize_t count[1] = {static_cast<hsize_t>(rows)};
        hid_t space = H5Dget_space(dataset);
        requireHandle(space, "Cannot obtain verification file space");
        requireStatus(H5Sselect_hyperslab(space, H5S_SELECT_SET,
          start, NULL, count, NULL), "Cannot select verification block");
        hid_t memorySpace = H5Screate_simple(1, count, NULL);
        requireHandle(memorySpace, "Cannot create verification memory space");
        requireStatus(H5Dread(dataset, memoryType, memorySpace, space,
          H5P_DEFAULT, &records[0]), "Cannot read synthetic record block");
        H5Sclose(memorySpace);
        H5Sclose(space);
        for (std::uint64_t index = 0; index < rows; ++index)
          {
            const ParticleRecord& record = records[index];
            const double ux = record.properVelocity[0];
            const double uy = record.properVelocity[1];
            const double uz = record.properVelocity[2];
            const double gamma = std::hypot(std::hypot(1.0, ux),
              std::hypot(uy, uz));
            const double betaZ = uz / gamma;
            const double z = -speedOfLight * betaZ *
              record.arrivalTimeOffset;
            if (record.sourceId != offset + index + 1 ||
                record.macroWeight != 1.0 ||
                !(uz > 0.0) || !std::isfinite(betaZ) ||
                !std::isfinite(record.planePosition[0]) ||
                !std::isfinite(record.planePosition[1]) ||
                !std::isfinite(record.arrivalTimeOffset))
              throw std::runtime_error(
                "Synthetic bunch contains an invalid particle record");
            sumOffset += record.arrivalTimeOffset;
            sumBetaZ += betaZ;
            minimumZ = std::min(minimumZ, z);
            maximumZ = std::max(maximumZ, z);
          }
        offset += rows;
      }
    H5Tclose(memoryType);
    H5Dclose(dataset);
    H5Gclose(group);
    H5Fclose(file);
    const double meanOffset = static_cast<double>(sumOffset / expected);
    const double measuredCurrent = electrons * electronCharge *
      static_cast<double>(sumBetaZ / expected) * speedOfLight / length;
    if (std::abs(meanOffset) > 1.0e-22 ||
        std::abs(measuredCurrent / current - 1.0) > 1.0e-10 ||
        maximumZ - minimumZ > 1.01 * length)
      throw std::runtime_error(
        "Synthetic bunch centered time, current or width check failed");
    std::cout << std::setprecision(12)
      << "Verified independent HDF5-v4 bunch\n"
      << "  file: " << filename << "\n"
      << "  rows: " << expected << "\n"
      << "  represented electrons: " << electrons << "\n"
      << "  mean arrival offset [s]: " << meanOffset << "\n"
      << "  reconstructed z span [m]: "
      << minimumZ << " to " << maximumZ << "\n"
      << "  measured current [A]: " << measuredCurrent << '\n';
  }

  void generate(const Options& options, std::uint64_t count,
                double betaReference)
  {
    std::cout << "Generating " << count << " independent particles in two "
                 "bounded-memory passes..." << std::endl;
    DeterministicRandom firstPass(options.seed);
    Statistics statistics;
    for (std::uint64_t index = 0; index < count; ++index)
      {
        statistics.add(sampleParticle(options, firstPass));
        if ((index + 1) % UINT64_C(2000000) == 0)
          std::cout << "  statistics: " << index + 1 << " / " << count
                    << std::endl;
      }
    const double referenceTime =
      static_cast<double>(statistics.time / count);
    const double meanBetaZ =
      static_cast<double>(statistics.betaZ / count);
    const double representedElectrons =
      static_cast<double>(count) * options.macroWeightElectrons;
    const double charge = representedElectrons * electronCharge;
    const double realizedCurrentReference = charge * betaReference *
      speedOfLight / options.length;
    const double realizedCurrentMean = charge * meanBetaZ *
      speedOfLight / options.length;

    hid_t file = H5Fcreate(options.output.c_str(),
      options.overwrite ? H5F_ACC_TRUNC : H5F_ACC_EXCL,
      H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(file, "Cannot create output file (use --overwrite only "
      "for an intentional replacement): " + options.output);
    hid_t group = H5Gcreate2(file, "/particles", H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(group, "Cannot create /particles group");
    writeIntAttribute(group, "format_version", 4);
    writeIntAttribute(group, "complete", 0);
    writeStringAttribute(group, "generator",
      "aprl_random_uniform_lab_plane_bunch_generator");
    writeStringAttribute(group, "coordinate_frame",
      "fixed-lab-plane-crossing-events");
    writeStringAttribute(group, "arrival_time_reference",
      "mean-crossing-time-of-synthetic-common-lab-time-snapshot");
    writeStringAttribute(group, "longitudinal_distribution",
      "random-uniform-in-common-lab-time-z");
    writeStringAttribute(group, "transverse_distribution",
      "independent-truncated-gaussian-twiss");
    writeStringAttribute(group, "plane_position_unit", "m");
    writeStringAttribute(group, "arrival_time_offset_unit", "s");
    writeStringAttribute(group, "proper_velocity_unit", "gamma*v/c");
    writeStringAttribute(group, "source_id_definition",
      "one-based-generation-order");
    writeStringAttribute(group, "macro_weight_definition",
      "unit-relative-weight; each row represents --macro-weight-electrons");
    writeDoubleAttribute(group, "requested_current_A", options.current);
    writeDoubleAttribute(group, "realized_reference_current_A",
      realizedCurrentReference);
    writeDoubleAttribute(group, "realized_mean_current_A",
      realizedCurrentMean);
    writeDoubleAttribute(group, "represented_electrons",
      representedElectrons);
    writeDoubleAttribute(group, "represented_charge_C", charge);
    writeDoubleAttribute(group, "macro_weight_electrons",
      options.macroWeightElectrons);
    writeDoubleAttribute(group, "requested_full_width_m", options.length);
    writeDoubleAttribute(group, "central_momentum_MeV_c",
      options.centralMomentumMeVc);
    writeDoubleAttribute(group, "source_arrival_time_reference_s",
      referenceTime);
    writeDoubleAttribute(group, "emit_x_m", options.emitX);
    writeDoubleAttribute(group, "emit_y_m", options.emitY);
    writeDoubleAttribute(group, "beta_x_m", options.betaX);
    writeDoubleAttribute(group, "beta_y_m", options.betaY);
    writeDoubleAttribute(group, "alpha_x", options.alphaX);
    writeDoubleAttribute(group, "alpha_y", options.alphaY);
    writeDoubleAttribute(group, "sigma_dp", options.sigmaDp);
    writeDoubleAttribute(group, "cutoff_sigma", options.cutoffSigma);
    writeDoubleAttribute(group, "mean_beta_z", meanBetaZ);
    writeDoubleAttribute(group, "rms_x_m",
      rms(statistics.x, statistics.x2, count));
    writeDoubleAttribute(group, "rms_y_m",
      rms(statistics.y, statistics.y2, count));
    writeDoubleAttribute(group, "rms_xp",
      rms(statistics.xp, statistics.xp2, count));
    writeDoubleAttribute(group, "rms_yp",
      rms(statistics.yp, statistics.yp2, count));
    writeUInt64Attribute(group, "random_seed", options.seed);
    writeUInt64Attribute(group, "macro_particles", count);

    const hsize_t dimensions[1] = {static_cast<hsize_t>(count)};
    hid_t space = H5Screate_simple(1, dimensions, NULL);
    requireHandle(space, "Cannot create HDF5 output dataspace");
    hid_t fileType = fileRecordType();
    hid_t properties = H5Pcreate(H5P_DATASET_CREATE);
    requireHandle(properties, "Cannot create HDF5 dataset properties");
    const hsize_t chunk[1] = {static_cast<hsize_t>(
      std::min(options.chunkRecords, count))};
    requireStatus(H5Pset_chunk(properties, 1, chunk),
      "Cannot set HDF5 chunk shape");
    if (options.compression > 0)
      requireStatus(H5Pset_deflate(properties, options.compression),
        "Cannot enable HDF5 gzip compression");
    hid_t dataset = H5Dcreate2(group, "records", fileType, space,
      H5P_DEFAULT, properties, H5P_DEFAULT);
    requireHandle(dataset, "Cannot create /particles/records");
    H5Pclose(properties);
    H5Tclose(fileType);
    H5Sclose(space);
    hid_t memoryType = memoryRecordType();

    DeterministicRandom secondPass(options.seed);
    std::vector<ParticleRecord> records;
    records.reserve(static_cast<std::size_t>(options.chunkRecords));
    std::uint64_t written = 0;
    for (std::uint64_t index = 0; index < count; ++index)
      {
        const Sample sample = sampleParticle(options, secondPass);
        ParticleRecord record = {};
        record.planePosition[0] = sample.x - sample.xp * sample.z;
        record.planePosition[1] = sample.y - sample.yp * sample.z;
        record.arrivalTimeOffset = sample.crossingTime - referenceTime;
        record.properVelocity[0] = sample.ux;
        record.properVelocity[1] = sample.uy;
        record.properVelocity[2] = sample.uz;
        record.sourceId = index + 1;
        record.macroWeight = 1.0;
        records.push_back(record);
        if (records.size() == options.chunkRecords || index + 1 == count)
          {
            writeChunk(dataset, memoryType, records, written);
            written += records.size();
            records.clear();
            if (written % UINT64_C(2000000) < options.chunkRecords ||
                written == count)
              std::cout << "  HDF5 records: " << written << " / " << count
                        << std::endl;
          }
      }
    H5Tclose(memoryType);
    H5Dclose(dataset);
    requireStatus(H5Adelete(group, "complete"),
      "Cannot replace completion attribute");
    writeIntAttribute(group, "complete", 1);
    requireStatus(H5Fflush(file, H5F_SCOPE_GLOBAL),
      "Cannot flush synthetic bunch");
    H5Gclose(group);
    H5Fclose(file);

    std::cout << std::setprecision(12)
      << "Generated independent random-uniform APRL HDF5-v4 bunch\n"
      << "  output: " << options.output << "\n"
      << "  full uniform z width [m]: " << options.length << "\n"
      << "  macro-particles: " << count << "\n"
      << "  electrons per macro-particle: " << options.macroWeightElectrons
      << "\n"
      << "  represented electrons: " << representedElectrons << "\n"
      << "  represented charge [C]: " << charge << "\n"
      << "  requested/reference/mean-beta current [A]: "
      << options.current << " / " << realizedCurrentReference
      << " / " << realizedCurrentMean << "\n"
      << "  rms x/y [um]: "
      << 1.0e6 * rms(statistics.x, statistics.x2, count) << " / "
      << 1.0e6 * rms(statistics.y, statistics.y2, count) << "\n"
      << "  rms x'/y' [urad]: "
      << 1.0e6 * rms(statistics.xp, statistics.xp2, count) << " / "
      << 1.0e6 * rms(statistics.yp, statistics.yp2, count) << "\n"
      << "Use beam.input.electrons=" << representedElectrons
      << " and input_plane_z=activation_plane_z at the synthetic plane.\n";
  }
}

int main(int argc, char** argv)
{
  try
    {
      const Options options = parseOptions(argc, argv);
      validateOptions(options);
      if (!options.verify.empty())
        {
          verifyFile(options.verify);
          return EXIT_SUCCESS;
        }
      double betaReference = 0.0;
      const std::uint64_t count = particleCount(options, betaReference);
      const double representedElectrons =
        static_cast<double>(count) * options.macroWeightElectrons;
      const double realizedCurrent = representedElectrons * electronCharge *
        betaReference * speedOfLight / options.length;
      if (options.plan)
        {
          std::cout << std::setprecision(12)
            << "Independent random-uniform HDF5-v4 bunch plan\n"
            << "  full uniform z width [m]: " << options.length << "\n"
            << "  requested current [A]: " << options.current << "\n"
            << "  macro weight [electrons]: "
            << options.macroWeightElectrons << "\n"
            << "  macro-particles: " << count << "\n"
            << "  represented electrons: " << representedElectrons << "\n"
            << "  realized reference current [A]: "
            << realizedCurrent << "\n"
            << "  raw record bytes: " << count * UINT64_C(64) << "\n";
        }
      else
        generate(options, count, betaReference);
      return EXIT_SUCCESS;
    }
  catch (const std::exception& error)
    {
      std::cerr << "uniform_lab_plane_bunch_generator: "
                << error.what() << '\n';
      return EXIT_FAILURE;
    }
}
