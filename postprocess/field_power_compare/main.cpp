#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "fftw3.h"
#include "hdf5.h"
#include "yaml-cpp/yaml.h"

namespace
{
  const double kMu0 = 1.25663706212e-6;
  const double kEpsilon0 = 8.8541878128e-12;
  const double kC = 299792458.0;
  const double kPi = 3.141592653589793238462643383279502884;
  const double kPlanckEVSecond = 4.135667696e-15;

  void requireHandle(hid_t handle, const std::string& message)
  {
    if (handle < 0) throw std::runtime_error(message);
  }

  void requireStatus(herr_t status, const std::string& message)
  {
    if (status < 0) throw std::runtime_error(message);
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

  YAML::Node required(const YAML::Node& parent, const char* key)
  {
    const YAML::Node value = parent[key];
    if (!value) throw std::runtime_error(
      std::string("Missing required YAML key '") + key + "'");
    return value;
  }

  struct Config
  {
    std::string signalFile;
    std::string baselineFile;
    std::string trajectoryFile;
    bool requireComplete;
    double minimumPhotonEnergyEV;
    double maximumPhotonEnergyEV;
    std::size_t batchPoints;
    std::string outputFile;
    unsigned int compression;
  };

  Config loadConfig(const std::string& filename)
  {
    const YAML::Node root = YAML::LoadFile(filename);
    if (!root.IsMap())
      throw std::runtime_error("Field-power comparison card must be a map");
    const YAML::Node input = required(root, "input");
    const YAML::Node analysis = required(root, "analysis");
    const YAML::Node output = required(root, "output");
    if (!input.IsMap() || !analysis.IsMap() || !output.IsMap())
      throw std::runtime_error("input, analysis, and output must be maps");
    Config config;
    config.signalFile = resolvePath(filename,
      required(input, "signal_field").as<std::string>());
    config.baselineFile = resolvePath(filename,
      required(input, "zero_radiation_baseline").as<std::string>());
    config.trajectoryFile = input["trajectory_far_field"] ?
      resolvePath(filename,
        input["trajectory_far_field"].as<std::string>()) : std::string();
    config.requireComplete = input["require_complete"] ?
      input["require_complete"].as<bool>() : true;
    const YAML::Node band = required(analysis, "photon_energy_band_eV");
    if (!band.IsSequence() || band.size() != 2)
      throw std::runtime_error(
        "analysis.photon_energy_band_eV must be [minimum, maximum]");
    config.minimumPhotonEnergyEV = band[0].as<double>();
    config.maximumPhotonEnergyEV = band[1].as<double>();
    if (!(config.minimumPhotonEnergyEV >= 0.0) ||
        !(config.maximumPhotonEnergyEV > config.minimumPhotonEnergyEV) ||
        !std::isfinite(config.minimumPhotonEnergyEV) ||
        !std::isfinite(config.maximumPhotonEnergyEV))
      throw std::runtime_error("Invalid photon-energy band");
    config.batchPoints = analysis["spatial_batch_points"] ?
      analysis["spatial_batch_points"].as<std::size_t>() : 16;
    if (config.batchPoints == 0)
      throw std::runtime_error("spatial_batch_points must be positive");
    config.outputFile = resolvePath(filename,
      required(output, "file").as<std::string>());
    config.compression = output["compression"] ?
      output["compression"].as<unsigned int>() : 0;
    if (config.compression > 9)
      throw std::runtime_error("output.compression must be in [0,9]");
    if (config.signalFile == config.baselineFile ||
        config.outputFile == config.signalFile ||
        config.outputFile == config.baselineFile)
      throw std::runtime_error("Signal, baseline, and output paths must differ");
    return config;
  }

  double readDoubleAttribute(hid_t object, const char* name)
  {
    hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
    requireHandle(attribute, std::string("Missing HDF5 attribute: ") + name);
    double value = 0.0;
    const herr_t status = H5Aread(attribute, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attribute);
    requireStatus(status, std::string("Cannot read HDF5 attribute: ") + name);
    return value;
  }

  unsigned long long readUnsignedAttribute(hid_t object, const char* name)
  {
    hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
    requireHandle(attribute, std::string("Missing HDF5 attribute: ") + name);
    unsigned long long value = 0;
    const herr_t status = H5Aread(attribute,
      H5T_NATIVE_ULLONG, &value);
    H5Aclose(attribute);
    requireStatus(status, std::string("Cannot read HDF5 attribute: ") + name);
    return value;
  }

  unsigned long long readUnsignedDataset(hid_t group, const char* name)
  {
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
    requireHandle(dataset, std::string("Missing HDF5 dataset: ") + name);
    unsigned long long value = 0;
    const herr_t status = H5Dread(dataset, H5T_NATIVE_ULLONG,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    H5Dclose(dataset);
    requireStatus(status, std::string("Cannot read HDF5 dataset: ") + name);
    return value;
  }

  unsigned char readByteDataset(hid_t group, const char* name)
  {
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
    requireHandle(dataset, std::string("Missing HDF5 dataset: ") + name);
    unsigned char value = 0;
    const herr_t status = H5Dread(dataset, H5T_NATIVE_UCHAR,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    H5Dclose(dataset);
    requireStatus(status, std::string("Cannot read HDF5 dataset: ") + name);
    return value;
  }

  struct FieldPoint
  {
    double electric[3];
    double magnetic[3];
  };

  hid_t createFieldType()
  {
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(FieldPoint));
    requireHandle(type, "Cannot create field memory type");
    hsize_t three[1] = {3};
    hid_t vector = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, three);
    requireHandle(vector, "Cannot create vector memory type");
    requireStatus(H5Tinsert(type, "electric_V_per_m",
      HOFFSET(FieldPoint, electric), vector), "Cannot map electric field");
    requireStatus(H5Tinsert(type, "magnetic_T",
      HOFFSET(FieldPoint, magnetic), vector), "Cannot map magnetic field");
    H5Tclose(vector);
    return type;
  }

  class FieldFile
  {
  public:
    FieldFile(const std::string& filename, bool requireComplete)
      : filename_(filename), file_(-1), group_(-1), fields_(-1),
        memoryType_(-1), samples_(0), nx_(0), ny_(0), planeZ_(0.0),
        xFirst_(0.0), yFirst_(0.0), dx_(0.0), dy_(0.0), times_()
    {
      file_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      requireHandle(file_, "Cannot open field file: " + filename);
      try
        {
          group_ = H5Gopen2(file_, "/field_plane", H5P_DEFAULT);
          requireHandle(group_, "Missing /field_plane in " + filename);
          if (requireComplete && readByteDataset(group_, "complete") == 0)
            throw std::runtime_error("Incomplete field file: " + filename);
          samples_ = static_cast<std::size_t>(
            readUnsignedDataset(group_, "committed_samples"));
          nx_ = static_cast<std::size_t>(readUnsignedAttribute(group_, "nx"));
          ny_ = static_cast<std::size_t>(readUnsignedAttribute(group_, "ny"));
          planeZ_ = readDoubleAttribute(group_, "plane_z_m");
          xFirst_ = readDoubleAttribute(group_, "x_first_m");
          yFirst_ = readDoubleAttribute(group_, "y_first_m");
          dx_ = readDoubleAttribute(group_, "dx_m");
          dy_ = readDoubleAttribute(group_, "dy_m");
          if (samples_ < 3 || nx_ == 0 || ny_ == 0 ||
              !(dx_ > 0.0) || !(dy_ > 0.0))
            throw std::runtime_error("Invalid field-plane geometry: " + filename);
          hid_t time = H5Dopen2(group_, "time_s", H5P_DEFAULT);
          requireHandle(time, "Missing field time axis: " + filename);
          times_.resize(samples_);
          hsize_t start[1] = {0};
          hsize_t count[1] = {static_cast<hsize_t>(samples_)};
          hid_t fileSpace = H5Dget_space(time);
          requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
            start, NULL, count, NULL), "Cannot select field time axis");
          hid_t memorySpace = H5Screate_simple(1, count, NULL);
          const herr_t status = H5Dread(time, H5T_NATIVE_DOUBLE,
            memorySpace, fileSpace, H5P_DEFAULT, &times_[0]);
          H5Sclose(memorySpace);
          H5Sclose(fileSpace);
          H5Dclose(time);
          requireStatus(status, "Cannot read field time axis");
          fields_ = H5Dopen2(group_, "fields", H5P_DEFAULT);
          requireHandle(fields_, "Missing field samples: " + filename);
          memoryType_ = createFieldType();
        }
      catch (...)
        {
          close();
          throw;
        }
    }

    ~FieldFile() { close(); }

    void readBlock(std::size_t y, std::size_t x, std::size_t count,
                   std::vector<FieldPoint>& values) const
    {
      if (y >= ny_ || x + count > nx_ || count == 0)
        throw std::out_of_range("Field block is outside detector plane");
      values.resize(samples_ * count);
      hsize_t start[3] = {0, static_cast<hsize_t>(y),
                          static_cast<hsize_t>(x)};
      hsize_t size[3] = {static_cast<hsize_t>(samples_), 1,
                         static_cast<hsize_t>(count)};
      hid_t fileSpace = H5Dget_space(fields_);
      requireHandle(fileSpace, "Cannot inspect field samples");
      requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
        start, NULL, size, NULL), "Cannot select field block");
      hid_t memorySpace = H5Screate_simple(3, size, NULL);
      const herr_t status = H5Dread(fields_, memoryType_, memorySpace,
        fileSpace, H5P_DEFAULT, &values[0]);
      H5Sclose(memorySpace);
      H5Sclose(fileSpace);
      requireStatus(status, "Cannot read field block");
    }

    void close()
    {
      if (memoryType_ >= 0) H5Tclose(memoryType_);
      if (fields_ >= 0) H5Dclose(fields_);
      if (group_ >= 0) H5Gclose(group_);
      if (file_ >= 0) H5Fclose(file_);
      memoryType_ = fields_ = group_ = file_ = -1;
    }

    const std::string& filename() const { return filename_; }
    std::size_t samples() const { return samples_; }
    std::size_t nx() const { return nx_; }
    std::size_t ny() const { return ny_; }
    double planeZ() const { return planeZ_; }
    double xFirst() const { return xFirst_; }
    double yFirst() const { return yFirst_; }
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    const std::vector<double>& times() const { return times_; }

  private:
    FieldFile(const FieldFile&);
    FieldFile& operator=(const FieldFile&);
    std::string filename_;
    hid_t file_;
    hid_t group_;
    hid_t fields_;
    hid_t memoryType_;
    std::size_t samples_;
    std::size_t nx_;
    std::size_t ny_;
    double planeZ_;
    double xFirst_;
    double yFirst_;
    double dx_;
    double dy_;
    std::vector<double> times_;
  };

  bool nearlyEqual(double left, double right)
  {
    const double scale = std::max(std::numeric_limits<double>::min(),
      std::max(std::abs(left), std::abs(right)));
    return std::abs(left - right) <=
      256.0 * std::numeric_limits<double>::epsilon() * scale;
  }

  void validateMatched(const FieldFile& signal, const FieldFile& baseline)
  {
    if (signal.samples() != baseline.samples() ||
        signal.nx() != baseline.nx() || signal.ny() != baseline.ny() ||
        !nearlyEqual(signal.planeZ(), baseline.planeZ()) ||
        !nearlyEqual(signal.xFirst(), baseline.xFirst()) ||
        !nearlyEqual(signal.yFirst(), baseline.yFirst()) ||
        !nearlyEqual(signal.dx(), baseline.dx()) ||
        !nearlyEqual(signal.dy(), baseline.dy()))
      throw std::runtime_error(
        "Signal and zero-radiation baseline have different detector geometry");
    for (std::size_t sample = 0; sample < signal.samples(); ++sample)
      if (!nearlyEqual(signal.times()[sample], baseline.times()[sample]))
        throw std::runtime_error(
          "Signal and zero-radiation baseline have different sample times");
  }

  class BandFilter
  {
  public:
    BandFilter(std::size_t samples, double timeStep,
               double minimumEV, double maximumEV)
      : samples_(samples), input_(fftw_alloc_real(samples)),
        output_(fftw_alloc_real(samples)),
        spectrum_(fftw_alloc_complex(samples / 2 + 1)),
        forward_(NULL), inverse_(NULL), keptBins_(0)
    {
      if (!input_ || !output_ || !spectrum_)
        throw std::bad_alloc();
      forward_ = fftw_plan_dft_r2c_1d(static_cast<int>(samples_), input_,
        spectrum_, FFTW_ESTIMATE);
      inverse_ = fftw_plan_dft_c2r_1d(static_cast<int>(samples_), spectrum_,
        output_, FFTW_ESTIMATE);
      if (!forward_ || !inverse_)
        throw std::runtime_error("Cannot create FFTW band-filter plans");
      keep_.resize(samples_ / 2 + 1, false);
      for (std::size_t bin = 0; bin < keep_.size(); ++bin)
        {
          const double energyEV = kPlanckEVSecond *
            static_cast<double>(bin) /
            (static_cast<double>(samples_) * timeStep);
          keep_[bin] = energyEV >= minimumEV && energyEV <= maximumEV;
          if (keep_[bin]) ++keptBins_;
        }
      if (keptBins_ == 0)
        throw std::runtime_error(
          "The requested photon-energy band contains no discrete FFT bins");
    }

    ~BandFilter()
    {
      if (forward_) fftw_destroy_plan(forward_);
      if (inverse_) fftw_destroy_plan(inverse_);
      fftw_free(spectrum_);
      fftw_free(output_);
      fftw_free(input_);
    }

    void apply(const std::vector<double>& input, std::vector<double>& output)
    {
      if (input.size() != samples_)
        throw std::invalid_argument("Band-filter series has wrong length");
      std::copy(input.begin(), input.end(), input_);
      fftw_execute(forward_);
      for (std::size_t bin = 0; bin < keep_.size(); ++bin)
        if (!keep_[bin])
          spectrum_[bin][0] = spectrum_[bin][1] = 0.0;
      fftw_execute(inverse_);
      output.resize(samples_);
      const double normalization = 1.0 / static_cast<double>(samples_);
      for (std::size_t sample = 0; sample < samples_; ++sample)
        output[sample] = output_[sample] * normalization;
    }

    std::size_t keptBins() const { return keptBins_; }

  private:
    BandFilter(const BandFilter&);
    BandFilter& operator=(const BandFilter&);
    std::size_t samples_;
    double* input_;
    double* output_;
    fftw_complex* spectrum_;
    fftw_plan forward_;
    fftw_plan inverse_;
    std::vector<bool> keep_;
    std::size_t keptBins_;
  };

  struct PowerSeries
  {
    std::vector<double> signedPower;
    std::vector<double> forwardPower;

    explicit PowerSeries(std::size_t samples)
      : signedPower(samples, 0.0), forwardPower(samples, 0.0)
    {}
  };

  double component(const FieldPoint& value, unsigned int component)
  {
    if (component == 0) return value.electric[0];
    if (component == 1) return value.electric[1];
    if (component == 2) return value.magnetic[0];
    return value.magnetic[1];
  }

  double scenarioComponent(unsigned int scenario,
                           const FieldPoint& signal,
                           const FieldPoint& baseline,
                           unsigned int fieldComponent)
  {
    if (scenario == 0) return component(signal, fieldComponent);
    if (scenario == 1) return component(baseline, fieldComponent);
    return component(signal, fieldComponent) -
      component(baseline, fieldComponent);
  }

  void accumulatePoynting(PowerSeries& power, std::size_t sample,
                          double ex, double ey, double bx, double by,
                          double area)
  {
    const double value = (ex * by - ey * bx) * area / kMu0;
    power.signedPower[sample] += value;
    power.forwardPower[sample] += std::max(0.0, value);
  }

  double integrate(const std::vector<double>& time,
                   const std::vector<double>& values)
  {
    double result = 0.0;
    for (std::size_t sample = 1; sample < time.size(); ++sample)
      result += 0.5 * (values[sample - 1] + values[sample]) *
        (time[sample] - time[sample - 1]);
    return result;
  }

  double peak(const std::vector<double>& values)
  {
    return *std::max_element(values.begin(), values.end());
  }

  void resampleLinear(const std::vector<double>& sourceTime,
                      const std::vector<double>& source,
                      const std::vector<double>& targetTime,
                      std::vector<double>& target)
  {
    target.resize(targetTime.size());
    std::size_t upper = 1;
    for (std::size_t sample = 0; sample < targetTime.size(); ++sample)
      {
        while (upper + 1 < sourceTime.size() &&
               sourceTime[upper] < targetTime[sample])
          ++upper;
        if (targetTime[sample] <= sourceTime.front())
          {
            target[sample] = source.front();
            continue;
          }
        if (targetTime[sample] >= sourceTime.back())
          {
            target[sample] = source.back();
            continue;
          }
        const std::size_t lower = upper - 1;
        const double fraction = (targetTime[sample] - sourceTime[lower]) /
          (sourceTime[upper] - sourceTime[lower]);
        target[sample] = source[lower] + fraction *
          (source[upper] - source[lower]);
      }
  }

  void writeStringAttribute(hid_t object, const char* name,
                            const std::string& value)
  {
    hid_t type = H5Tcopy(H5T_C_S1);
    requireStatus(H5Tset_size(type, value.size() + 1),
      "Cannot size string attribute");
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attribute = H5Acreate2(object, name, type, space,
      H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute: ") + name);
    const herr_t status = H5Awrite(attribute, type, value.c_str());
    H5Aclose(attribute);
    H5Tclose(type);
    requireStatus(status, std::string("Cannot write attribute: ") + name);
  }

  void writeDoubleAttribute(hid_t object, const char* name, double value)
  {
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attribute = H5Acreate2(object, name, H5T_IEEE_F64LE, space,
      H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute: ") + name);
    const herr_t status = H5Awrite(attribute, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attribute);
    requireStatus(status, std::string("Cannot write attribute: ") + name);
  }

  void writeVector(hid_t group, const char* name,
                   const std::vector<double>& values,
                   unsigned int compression, const char* unit)
  {
    hsize_t dimensions[1] = {static_cast<hsize_t>(values.size())};
    hid_t space = H5Screate_simple(1, dimensions, NULL);
    hid_t creation = H5P_DEFAULT;
    hid_t ownedCreation = -1;
    if (compression > 0)
      {
        ownedCreation = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk[1] = {static_cast<hsize_t>(
          std::min<std::size_t>(values.size(), 4096))};
        requireStatus(H5Pset_chunk(ownedCreation, 1, chunk),
          "Cannot set output chunk");
        requireStatus(H5Pset_shuffle(ownedCreation),
          "Cannot set output shuffle");
        requireStatus(H5Pset_deflate(ownedCreation, compression),
          "Cannot set output compression");
        creation = ownedCreation;
      }
    hid_t dataset = H5Dcreate2(group, name, H5T_IEEE_F64LE, space,
      H5P_DEFAULT, creation, H5P_DEFAULT);
    if (ownedCreation >= 0) H5Pclose(ownedCreation);
    H5Sclose(space);
    requireHandle(dataset, std::string("Cannot create dataset: ") + name);
    requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &values[0]), std::string("Cannot write dataset: ") + name);
    writeStringAttribute(dataset, "unit", unit);
    H5Dclose(dataset);
  }

  struct ComplexValue
  {
    double real;
    double imag;
  };

  struct TrajectoryPower
  {
    double energy;
    double peakPower;
    double period;
    std::vector<double> relativeTime;
    std::vector<double> power;

    TrajectoryPower()
      : energy(0.0), peakPower(0.0), period(0.0),
        relativeTime(), power()
    {}
  };

  std::vector<double> readDoubleVector(hid_t group, const char* name)
  {
    hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
    requireHandle(dataset, std::string("Missing HDF5 dataset: ") + name);
    hid_t space = H5Dget_space(dataset);
    hsize_t count = 0;
    if (H5Sget_simple_extent_ndims(space) != 1 ||
        H5Sget_simple_extent_dims(space, &count, NULL) != 1 || count == 0)
      {
        H5Sclose(space);
        H5Dclose(dataset);
        throw std::runtime_error(std::string("Invalid HDF5 axis: ") + name);
      }
    H5Sclose(space);
    std::vector<double> values(static_cast<std::size_t>(count));
    const herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &values[0]);
    H5Dclose(dataset);
    requireStatus(status, std::string("Cannot read HDF5 axis: ") + name);
    return values;
  }

  double trajectorySolidAngleWeight(const std::vector<double>& thetaX,
                                    const std::vector<double>& thetaY,
                                    std::size_t y, std::size_t x)
  {
    const double dx = (thetaX.back() - thetaX.front()) /
      static_cast<double>(thetaX.size() - 1);
    const double dy = (thetaY.back() - thetaY.front()) /
      static_cast<double>(thetaY.size() - 1);
    const double tx = std::tan(thetaX[x]);
    const double ty = std::tan(thetaY[y]);
    const double jacobian = (1.0 + tx * tx) * (1.0 + ty * ty) /
      std::pow(1.0 + tx * tx + ty * ty, 1.5);
    const double endpointX =
      (x == 0 || x + 1 == thetaX.size()) ? 0.5 : 1.0;
    const double endpointY =
      (y == 0 || y + 1 == thetaY.size()) ? 0.5 : 1.0;
    return dx * dy * jacobian * endpointX * endpointY;
  }

  TrajectoryPower readTrajectoryPower(const std::string& filename,
                                      double minimumEV, double maximumEV,
                                      bool requireComplete)
  {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    requireHandle(file, "Cannot open trajectory far-field file: " + filename);
    hid_t group = H5Gopen2(file, "/far_field", H5P_DEFAULT);
    if (group < 0)
      {
        H5Fclose(file);
        throw std::runtime_error("Missing /far_field in " + filename);
      }
    try
      {
        if (requireComplete && readByteDataset(group, "complete") == 0)
          throw std::runtime_error("Incomplete trajectory far-field file");
        if (readUnsignedAttribute(group, "angular_integral_available") == 0)
          throw std::runtime_error(
            "Trajectory far-field angular grid cannot be integrated");
        const double storedMinimum = readDoubleAttribute(
          group, "band_min_photon_energy_eV");
        const double storedMaximum = readDoubleAttribute(
          group, "band_max_photon_energy_eV");
        if (!nearlyEqual(storedMinimum, minimumEV) ||
            !nearlyEqual(storedMaximum, maximumEV))
          throw std::runtime_error(
            "Trajectory and field photon-energy bands differ");
        if (readUnsignedAttribute(group, "ensemble_shots") != 1)
          throw std::runtime_error(
            "Trajectory time-power reconstruction currently requires exactly one shot");
        TrajectoryPower result;
        result.energy = readDoubleAttribute(group, "mean_band_energy_J");
        const double observerDistance = readDoubleAttribute(
          group, "observer_distance_m");
        const std::vector<double> omega = readDoubleVector(
          group, "omega_rad_per_s");
        const std::vector<double> thetaX = readDoubleVector(
          group, "theta_x_rad");
        const std::vector<double> thetaY = readDoubleVector(
          group, "theta_y_rad");
        if (omega.size() < 2 || thetaX.size() < 2 || thetaY.size() < 2)
          throw std::runtime_error(
            "Trajectory power needs at least two frequency and angle samples");
        const double deltaOmega = (omega.back() - omega.front()) /
          static_cast<double>(omega.size() - 1);
        for (std::size_t frequency = 1; frequency < omega.size(); ++frequency)
          if (std::abs((omega[frequency] - omega[frequency - 1]) -
                       deltaOmega) > 1.0e-8 * deltaOmega)
            throw std::runtime_error(
              "Trajectory time-power reconstruction requires linear frequency spacing");

        hid_t dataset = H5Dopen2(group, "electric_field_spectral",
          H5P_DEFAULT);
        requireHandle(dataset,
          "Missing trajectory electric_field_spectral dataset");
        hid_t fieldSpace = H5Dget_space(dataset);
        hsize_t dimensions[5] = {};
        if (H5Sget_simple_extent_ndims(fieldSpace) != 5 ||
            H5Sget_simple_extent_dims(fieldSpace, dimensions, NULL) != 5 ||
            dimensions[0] != 1 ||
            dimensions[1] != static_cast<hsize_t>(omega.size()) ||
            dimensions[2] != static_cast<hsize_t>(thetaY.size()) ||
            dimensions[3] != static_cast<hsize_t>(thetaX.size()) ||
            dimensions[4] != 2)
          {
            H5Sclose(fieldSpace);
            H5Dclose(dataset);
            throw std::runtime_error(
              "Unexpected trajectory spectral-field dimensions");
          }
        H5Sclose(fieldSpace);
        hid_t complexType = H5Tcreate(H5T_COMPOUND, sizeof(ComplexValue));
        requireHandle(complexType, "Cannot create trajectory complex type");
        requireStatus(H5Tinsert(complexType, "real",
          HOFFSET(ComplexValue, real), H5T_NATIVE_DOUBLE),
          "Cannot map trajectory complex real component");
        requireStatus(H5Tinsert(complexType, "imag",
          HOFFSET(ComplexValue, imag), H5T_NATIVE_DOUBLE),
          "Cannot map trajectory complex imaginary component");
        const std::size_t fieldValues = omega.size() * thetaY.size() *
          thetaX.size() * 2;
        std::vector<ComplexValue> field(fieldValues);
        const herr_t fieldStatus = H5Dread(dataset, complexType,
          H5S_ALL, H5S_ALL, H5P_DEFAULT, &field[0]);
        H5Tclose(complexType);
        H5Dclose(dataset);
        requireStatus(fieldStatus,
          "Cannot read trajectory spectral electric field");

        result.period = 2.0 * kPi / deltaOmega;
        const std::size_t samples = std::max<std::size_t>(
          256, 4 * omega.size());
        const double timeStep = result.period /
          static_cast<double>(samples);
        result.power.assign(samples, 0.0);
        const double inverseScale = deltaOmega / kPi;
        const double powerScale = 0.5 * kEpsilon0 * kC *
          observerDistance * observerDistance;
        for (std::size_t y = 0; y < thetaY.size(); ++y)
          for (std::size_t x = 0; x < thetaX.size(); ++x)
            {
              const double solidAngle = trajectorySolidAngleWeight(
                thetaX, thetaY, y, x);
              for (std::size_t sample = 0; sample < samples; ++sample)
                {
                  const double time = timeStep * static_cast<double>(sample);
                  double magnitudeSquared = 0.0;
                  for (std::size_t polarization = 0;
                       polarization < 2; ++polarization)
                    {
                      double real = 0.0;
                      double imag = 0.0;
                      for (std::size_t frequency = 0;
                           frequency < omega.size(); ++frequency)
                        {
                          const std::size_t index =
                            (((frequency * thetaY.size() + y) *
                              thetaX.size() + x) * 2 + polarization);
                          const double phase = std::remainder(
                            omega[frequency] * time, 2.0 * kPi);
                          const double cosine = std::cos(phase);
                          const double sine = std::sin(phase);
                          real += field[index].real * cosine +
                                  field[index].imag * sine;
                          imag += field[index].imag * cosine -
                                  field[index].real * sine;
                        }
                      magnitudeSquared += inverseScale * inverseScale *
                        (real * real + imag * imag);
                    }
                  result.power[sample] += powerScale * solidAngle *
                    magnitudeSquared;
                }
            }
        double periodicEnergy = 0.0;
        for (std::size_t sample = 0; sample < samples; ++sample)
          periodicEnergy += result.power[sample] * timeStep;
        if (!(periodicEnergy > 0.0) || !std::isfinite(periodicEnergy))
          throw std::runtime_error(
            "Invalid reconstructed trajectory time-domain energy");
        const double normalization = result.energy / periodicEnergy;
        for (std::size_t sample = 0; sample < samples; ++sample)
          result.power[sample] *= normalization;
        const std::size_t peakIndex = static_cast<std::size_t>(
          std::max_element(result.power.begin(), result.power.end()) -
          result.power.begin());
        std::vector<double> centered(samples);
        result.relativeTime.resize(samples);
        const std::size_t middle = samples / 2;
        for (std::size_t sample = 0; sample < samples; ++sample)
          {
            const std::size_t original =
              (peakIndex + sample + samples - middle) % samples;
            centered[sample] = result.power[original];
            result.relativeTime[sample] =
              (static_cast<double>(sample) -
               static_cast<double>(middle)) * timeStep;
          }
        result.power.swap(centered);
        result.peakPower = peak(result.power);
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

  double safeRatio(double numerator, double denominator)
  {
    return denominator != 0.0 ? numerator / denominator :
      std::numeric_limits<double>::quiet_NaN();
  }
}

int main(int argc, char** argv)
{
  if (argc != 2)
    {
      std::cerr << "Usage: field_power_compare <configuration.yaml>\n";
      return 2;
    }
  try
    {
      const Config config = loadConfig(argv[1]);
      FieldFile signal(config.signalFile, config.requireComplete);
      FieldFile baseline(config.baselineFile, config.requireComplete);
      validateMatched(signal, baseline);
      const std::vector<double>& time = signal.times();
      const double dt = (time.back() - time.front()) /
        static_cast<double>(time.size() - 1);
      if (!(dt > 0.0) || !std::isfinite(dt))
        throw std::runtime_error("Invalid detector sample interval");
      for (std::size_t sample = 1; sample < time.size(); ++sample)
        if (!(time[sample] > time[sample - 1]))
          throw std::runtime_error(
            "Detector laboratory sample times must be strictly increasing");
      std::vector<double> bandTime(time.size());
      for (std::size_t sample = 0; sample < time.size(); ++sample)
        bandTime[sample] = time.front() + dt * static_cast<double>(sample);
      const double nyquistEV = kPlanckEVSecond / (2.0 * dt);
      if (config.maximumPhotonEnergyEV > nyquistEV * (1.0 + 1.0e-12))
        throw std::runtime_error(
          "Requested photon-energy maximum exceeds detector Nyquist energy");

      PowerSeries raw[3] = {
        PowerSeries(time.size()), PowerSeries(time.size()),
        PowerSeries(time.size())
      };
      PowerSeries band[3] = {
        PowerSeries(time.size()), PowerSeries(time.size()),
        PowerSeries(time.size())
      };
      BandFilter filter(time.size(), dt, config.minimumPhotonEnergyEV,
                        config.maximumPhotonEnergyEV);
      const double area = signal.dx() * signal.dy();
      std::vector<FieldPoint> signalBlock;
      std::vector<FieldPoint> baselineBlock;
      std::vector<double> series(time.size());
      std::vector<double> uniformSeries(time.size());
      std::vector<double> filtered[4];

      for (std::size_t y = 0; y < signal.ny(); ++y)
        for (std::size_t x = 0; x < signal.nx();
             x += config.batchPoints)
          {
            const std::size_t count = std::min(config.batchPoints,
                                               signal.nx() - x);
            signal.readBlock(y, x, count, signalBlock);
            baseline.readBlock(y, x, count, baselineBlock);
            for (std::size_t point = 0; point < count; ++point)
              {
                for (std::size_t sample = 0; sample < time.size(); ++sample)
                  {
                    const std::size_t index = sample * count + point;
                    for (unsigned int scenario = 0; scenario < 3; ++scenario)
                      accumulatePoynting(raw[scenario], sample,
                        scenarioComponent(scenario, signalBlock[index],
                          baselineBlock[index], 0),
                        scenarioComponent(scenario, signalBlock[index],
                          baselineBlock[index], 1),
                        scenarioComponent(scenario, signalBlock[index],
                          baselineBlock[index], 2),
                        scenarioComponent(scenario, signalBlock[index],
                          baselineBlock[index], 3), area);
                  }
                for (unsigned int scenario = 0; scenario < 3; ++scenario)
                  {
                    for (unsigned int componentIndex = 0;
                         componentIndex < 4; ++componentIndex)
                      {
                        for (std::size_t sample = 0;
                             sample < time.size(); ++sample)
                          {
                            const std::size_t index = sample * count + point;
                            series[sample] = scenarioComponent(scenario,
                              signalBlock[index], baselineBlock[index],
                              componentIndex);
                          }
                        resampleLinear(time, series, bandTime, uniformSeries);
                        filter.apply(uniformSeries,
                                     filtered[componentIndex]);
                      }
                    for (std::size_t sample = 0;
                         sample < time.size(); ++sample)
                      accumulatePoynting(band[scenario], sample,
                        filtered[0][sample], filtered[1][sample],
                        filtered[2][sample], filtered[3][sample], area);
                  }
              }
          }

      const char* names[3] = {"signal", "baseline", "difference"};
      double rawEnergy[3] = {};
      double bandEnergy[3] = {};
      double rawPeak[3] = {};
      double bandPeak[3] = {};
      for (unsigned int scenario = 0; scenario < 3; ++scenario)
        {
          rawEnergy[scenario] = integrate(time, raw[scenario].forwardPower);
          bandEnergy[scenario] = integrate(
            bandTime, band[scenario].forwardPower);
          rawPeak[scenario] = peak(raw[scenario].forwardPower);
          bandPeak[scenario] = peak(band[scenario].forwardPower);
        }

      TrajectoryPower trajectory;
      if (!config.trajectoryFile.empty())
        trajectory = readTrajectoryPower(config.trajectoryFile,
          config.minimumPhotonEnergyEV, config.maximumPhotonEnergyEV,
          config.requireComplete);

      createDirectories(parentDirectory(config.outputFile));
      hid_t file = H5Fcreate(config.outputFile.c_str(), H5F_ACC_TRUNC,
        H5P_DEFAULT, H5P_DEFAULT);
      requireHandle(file, "Cannot create output: " + config.outputFile);
      hid_t group = H5Gcreate2(file, "/power_comparison", H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT);
      requireHandle(group, "Cannot create /power_comparison");
      try
        {
          writeVector(group, "time_s", time, config.compression, "s");
          writeVector(group, "band_time_s", bandTime,
            config.compression, "s");
          for (unsigned int scenario = 0; scenario < 3; ++scenario)
            {
              const std::string prefix(names[scenario]);
              writeVector(group, (prefix + "_signed_power_W").c_str(),
                raw[scenario].signedPower, config.compression, "W");
              writeVector(group, (prefix + "_forward_power_W").c_str(),
                raw[scenario].forwardPower, config.compression, "W");
              writeVector(group,
                (prefix + "_band_signed_power_W").c_str(),
                band[scenario].signedPower, config.compression, "W");
              writeVector(group,
                (prefix + "_band_forward_power_W").c_str(),
                band[scenario].forwardPower, config.compression, "W");
              writeDoubleAttribute(group,
                (prefix + "_forward_energy_J").c_str(),
                rawEnergy[scenario]);
              writeDoubleAttribute(group,
                (prefix + "_band_forward_energy_J").c_str(),
                bandEnergy[scenario]);
              writeDoubleAttribute(group,
                (prefix + "_peak_forward_power_W").c_str(),
                rawPeak[scenario]);
              writeDoubleAttribute(group,
                (prefix + "_band_peak_forward_power_W").c_str(),
                bandPeak[scenario]);
            }
          writeDoubleAttribute(group, "photon_energy_min_eV",
            config.minimumPhotonEnergyEV);
          writeDoubleAttribute(group, "photon_energy_max_eV",
            config.maximumPhotonEnergyEV);
          writeDoubleAttribute(group, "nyquist_photon_energy_eV", nyquistEV);
          writeDoubleAttribute(group,
            "baseline_to_difference_band_energy_ratio",
            safeRatio(bandEnergy[1], bandEnergy[2]));
          writeDoubleAttribute(group,
            "baseline_to_difference_band_peak_power_ratio",
            safeRatio(bandPeak[1], bandPeak[2]));
          writeDoubleAttribute(group, "plane_z_m", signal.planeZ());
          writeDoubleAttribute(group, "detector_area_m2",
            area * static_cast<double>(signal.nx() * signal.ny()));
          if (!config.trajectoryFile.empty())
            {
              writeVector(group, "trajectory_relative_time_s",
                trajectory.relativeTime, config.compression, "s");
              writeVector(group, "trajectory_band_power_W",
                trajectory.power, config.compression, "W");
              writeDoubleAttribute(group, "trajectory_band_energy_J",
                trajectory.energy);
              writeDoubleAttribute(group,
                "trajectory_band_peak_power_W", trajectory.peakPower);
              writeDoubleAttribute(group,
                "trajectory_reconstruction_period_s", trajectory.period);
              writeDoubleAttribute(group,
                "field_to_trajectory_band_energy_ratio",
                safeRatio(bandEnergy[2], trajectory.energy));
              writeDoubleAttribute(group,
                "field_to_trajectory_band_peak_power_ratio",
                safeRatio(bandPeak[2], trajectory.peakPower));
              writeStringAttribute(group,
                "trajectory_comparison_requirement",
                "trajectory angular grid must represent the same accepted radiation cone as the finite field plane");
              writeStringAttribute(group,
                "trajectory_power_definition",
                "cycle-averaged analytic-signal power reconstructed from the uniformly spaced positive-frequency far field; periodic time origin is shifted to center its peak and normalization is set by trajectory_band_energy_J");
            }
          writeStringAttribute(group, "difference_definition",
            "Poynting flux of E_signal-E_baseline and B_signal-B_baseline; powers are not subtracted");
          writeStringAttribute(group, "band_filter_definition",
            "linear resampling to a uniform time grid followed by rectangular finite-record DFT; bins outside requested photon-energy band are zeroed before inverse transform");
          writeStringAttribute(group, "forward_power_definition",
            "integral max((E cross B)_z/mu0,0) dx dy over stored cell centres");
          hsize_t scalar[1] = {1};
          hid_t space = H5Screate_simple(1, scalar, NULL);
          hid_t complete = H5Dcreate2(group, "complete", H5T_STD_U8LE,
            space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          H5Sclose(space);
          const unsigned char one = 1;
          requireStatus(H5Dwrite(complete, H5T_NATIVE_UCHAR, H5S_ALL,
            H5S_ALL, H5P_DEFAULT, &one), "Cannot commit comparison output");
          H5Dclose(complete);
          H5Gclose(group);
          H5Fclose(file);
        }
      catch (...)
        {
          H5Gclose(group);
          H5Fclose(file);
          throw;
        }

      std::cout << std::setprecision(10)
                << "Matched field-power comparison complete\n"
                << "  samples=" << time.size() << ", plane="
                << signal.nx() << "x" << signal.ny()
                << ", Nyquist=" << nyquistEV << " eV, retained FFT bins="
                << filter.keptBins() << "\n"
                << "  band [" << config.minimumPhotonEnergyEV << ", "
                << config.maximumPhotonEnergyEV << "] eV\n"
                << "  signal:     peak=" << bandPeak[0]
                << " W, energy=" << bandEnergy[0] << " J\n"
                << "  baseline:   peak=" << bandPeak[1]
                << " W, energy=" << bandEnergy[1] << " J\n"
                << "  difference: peak=" << bandPeak[2]
                << " W, energy=" << bandEnergy[2] << " J\n"
                << "  baseline/difference: peak="
                << safeRatio(bandPeak[1], bandPeak[2])
                << ", energy=" << safeRatio(bandEnergy[1], bandEnergy[2])
                << "\n";
      if (!config.trajectoryFile.empty())
        std::cout << "  trajectory: peak=" << trajectory.peakPower
                  << " W, band energy=" << trajectory.energy << " J\n"
                  << "  field/trajectory: peak="
                  << safeRatio(bandPeak[2], trajectory.peakPower)
                  << ", energy="
                  << safeRatio(bandEnergy[2], trajectory.energy) << "\n";
      return 0;
    }
  catch (const std::exception& error)
    {
      std::cerr << "Field-power comparison failure: " << error.what()
                << "\n";
      return 1;
    }
}
