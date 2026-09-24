#include "field_plane_analysis.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "fftw3.h"
#include "hdf5.h"

namespace field_analysis
{
  namespace
  {
    const double kPi = 3.141592653589793238462643383279502884;
    const double kC = 299792458.0;
    const double kEpsilon0 = 8.8541878128e-12;
    const double kPlanckEVSecond = 4.135667696e-15;

    void requireHandle(hid_t handle, const std::string& message)
    {
      if (handle < 0) throw std::runtime_error(message);
    }

    void requireStatus(herr_t status, const std::string& message)
    {
      if (status < 0) throw std::runtime_error(message);
    }

    std::size_t checkedProduct(std::size_t left, std::size_t right,
                               const char* description)
    {
      if (right != 0 && left >
          std::numeric_limits<std::size_t>::max() / right)
        throw std::overflow_error(std::string("Size overflow in ") +
          description);
      return left * right;
    }

    std::size_t checkedSum(std::size_t left, std::size_t right,
                           const char* description)
    {
      if (left > std::numeric_limits<std::size_t>::max() - right)
        throw std::overflow_error(std::string("Size overflow in ") +
          description);
      return left + right;
    }

    std::string parentDirectory(const std::string& path)
    {
      const std::string::size_type separator = path.find_last_of('/');
      if (separator == std::string::npos) return std::string();
      if (separator == 0) return "/";
      return path.substr(0, separator);
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
          const std::string::size_type separator =
            directory.find('/', offset);
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
                throw std::runtime_error("Cannot create directory: " +
                  current + ": " + std::strerror(errno));
            }
          if (separator == std::string::npos) break;
          offset = separator + 1;
        }
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

    std::uint64_t readUnsignedAttribute(hid_t object, const char* name)
    {
      hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
      requireHandle(attribute, std::string("Missing HDF5 attribute: ") + name);
      std::uint64_t value = 0;
      const herr_t status = H5Aread(attribute, H5T_NATIVE_UINT64, &value);
      H5Aclose(attribute);
      requireStatus(status, std::string("Cannot read HDF5 attribute: ") + name);
      return value;
    }

    std::uint64_t readUnsignedDataset(hid_t group, const char* name)
    {
      hid_t dataset = H5Dopen2(group, name, H5P_DEFAULT);
      requireHandle(dataset, std::string("Missing HDF5 dataset: ") + name);
      std::uint64_t value = 0;
      const herr_t status = H5Dread(dataset, H5T_NATIVE_UINT64,
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

    hid_t createFieldMemoryType()
    {
      hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(FieldPoint));
      requireHandle(type, "Cannot create field memory type");
      hsize_t three[1] = {3};
      hid_t vector = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, three);
      requireHandle(vector, "Cannot create field vector type");
      requireStatus(H5Tinsert(type, "electric_V_per_m",
        HOFFSET(FieldPoint, electric), vector), "Cannot map electric field");
      requireStatus(H5Tinsert(type, "magnetic_T",
        HOFFSET(FieldPoint, magnetic), vector), "Cannot map magnetic field");
      H5Tclose(vector);
      return type;
    }

    class FieldInput
    {
    public:
      FieldInput(const std::string& path, bool requireComplete)
        : path_(path), file_(-1), group_(-1), rawFields_(-1),
          reconstructedElectric_(-1), fieldType_(-1), reconstructed_(false),
          samples_(0), nx_(0), ny_(0), planeZ_(0.0), xFirst_(0.0),
          yFirst_(0.0), dx_(0.0), dy_(0.0), times_()
      {
        file_ = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        requireHandle(file_, "Cannot open field input: " + path);
        try
          {
            if (H5Lexists(file_, "/field_plane", H5P_DEFAULT) > 0)
              openRaw(requireComplete);
            else if (H5Lexists(file_, "/reconstructed_field",
                               H5P_DEFAULT) > 0)
              openReconstructed(requireComplete);
            else
              throw std::runtime_error(
                "Field input has neither /field_plane nor /reconstructed_field: " +
                path);
            validate();
          }
        catch (...)
          {
            close();
            throw;
          }
      }

      ~FieldInput() { close(); }

      void readElectricBlock(std::size_t y, std::size_t x,
                             std::size_t count,
                             std::vector<double>& values) const
      {
        if (y >= ny_ || count == 0 || x + count > nx_)
          throw std::out_of_range("Electric block is outside field plane");
        values.resize(checkedProduct(checkedProduct(samples_, count,
          "electric block"), 3, "electric block"));
        if (!reconstructed_)
          {
            std::vector<FieldPoint> points(checkedProduct(samples_, count,
              "raw field block"));
            hsize_t start[3] = {0, static_cast<hsize_t>(y),
              static_cast<hsize_t>(x)};
            hsize_t size[3] = {static_cast<hsize_t>(samples_), 1,
              static_cast<hsize_t>(count)};
            hid_t fileSpace = H5Dget_space(rawFields_);
            requireHandle(fileSpace, "Cannot inspect raw field dataset");
            requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
              start, NULL, size, NULL), "Cannot select raw field block");
            hid_t memorySpace = H5Screate_simple(3, size, NULL);
            requireHandle(memorySpace, "Cannot create raw field memory space");
            const herr_t status = H5Dread(rawFields_, fieldType_, memorySpace,
              fileSpace, H5P_DEFAULT, &points[0]);
            H5Sclose(memorySpace);
            H5Sclose(fileSpace);
            requireStatus(status, "Cannot read raw field block");
            for (std::size_t index = 0; index < points.size(); ++index)
              for (unsigned int component = 0; component < 3; ++component)
                values[3 * index + component] =
                  points[index].electric[component];
          }
        else
          {
            hsize_t start[4] = {0, static_cast<hsize_t>(y),
              static_cast<hsize_t>(x), 0};
            hsize_t size[4] = {static_cast<hsize_t>(samples_), 1,
              static_cast<hsize_t>(count), 3};
            hid_t fileSpace = H5Dget_space(reconstructedElectric_);
            requireHandle(fileSpace,
              "Cannot inspect reconstructed electric dataset");
            requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
              start, NULL, size, NULL),
              "Cannot select reconstructed electric block");
            hid_t memorySpace = H5Screate_simple(4, size, NULL);
            requireHandle(memorySpace,
              "Cannot create reconstructed electric memory space");
            const herr_t status = H5Dread(reconstructedElectric_,
              H5T_NATIVE_DOUBLE, memorySpace, fileSpace, H5P_DEFAULT,
              &values[0]);
            H5Sclose(memorySpace);
            H5Sclose(fileSpace);
            requireStatus(status,
              "Cannot read reconstructed electric block");
          }
      }

      const std::string& path() const { return path_; }
      const std::vector<double>& times() const { return times_; }
      std::size_t samples() const { return samples_; }
      std::size_t nx() const { return nx_; }
      std::size_t ny() const { return ny_; }
      double planeZ() const { return planeZ_; }
      double xFirst() const { return xFirst_; }
      double yFirst() const { return yFirst_; }
      double dx() const { return dx_; }
      double dy() const { return dy_; }
      bool reconstructed() const { return reconstructed_; }

    private:
      FieldInput(const FieldInput&);
      FieldInput& operator=(const FieldInput&);

      void openRaw(bool requireComplete)
      {
        group_ = H5Gopen2(file_, "/field_plane", H5P_DEFAULT);
        requireHandle(group_, "Cannot open /field_plane");
        if (requireComplete && readByteDataset(group_, "complete") == 0)
          throw std::runtime_error("Incomplete field input: " + path_);
        samples_ = static_cast<std::size_t>(
          readUnsignedDataset(group_, "committed_samples"));
        readGeometry();
        rawFields_ = H5Dopen2(group_, "fields", H5P_DEFAULT);
        requireHandle(rawFields_, "Missing /field_plane/fields");
        fieldType_ = createFieldMemoryType();
      }

      void openReconstructed(bool requireComplete)
      {
        reconstructed_ = true;
        group_ = H5Gopen2(file_, "/reconstructed_field", H5P_DEFAULT);
        requireHandle(group_, "Cannot open /reconstructed_field");
        if (requireComplete && readByteDataset(group_, "complete") == 0)
          throw std::runtime_error("Incomplete reconstructed input: " + path_);
        samples_ = static_cast<std::size_t>(
          readUnsignedAttribute(group_, "samples"));
        readGeometry();
        reconstructedElectric_ = H5Dopen2(group_, "electric_V_per_m",
          H5P_DEFAULT);
        requireHandle(reconstructedElectric_,
          "Missing reconstructed electric field");
      }

      void readGeometry()
      {
        nx_ = static_cast<std::size_t>(readUnsignedAttribute(group_, "nx"));
        ny_ = static_cast<std::size_t>(readUnsignedAttribute(group_, "ny"));
        planeZ_ = readDoubleAttribute(group_, "plane_z_m");
        xFirst_ = readDoubleAttribute(group_, "x_first_m");
        yFirst_ = readDoubleAttribute(group_, "y_first_m");
        dx_ = readDoubleAttribute(group_, "dx_m");
        dy_ = readDoubleAttribute(group_, "dy_m");
        hid_t time = H5Dopen2(group_, "time_s", H5P_DEFAULT);
        requireHandle(time, "Missing field time axis");
        times_.resize(samples_);
        hsize_t start[1] = {0};
        hsize_t count[1] = {static_cast<hsize_t>(samples_)};
        hid_t fileSpace = H5Dget_space(time);
        requireHandle(fileSpace, "Cannot inspect field time axis");
        requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
          start, NULL, count, NULL), "Cannot select committed time prefix");
        hid_t memorySpace = H5Screate_simple(1, count, NULL);
        requireHandle(memorySpace, "Cannot create time-axis memory space");
        const herr_t status = H5Dread(time, H5T_NATIVE_DOUBLE,
          memorySpace, fileSpace, H5P_DEFAULT,
          times_.empty() ? NULL : &times_[0]);
        H5Sclose(memorySpace);
        H5Sclose(fileSpace);
        H5Dclose(time);
        requireStatus(status, "Cannot read field time axis");
      }

      void validate() const
      {
        if (samples_ < 4 || nx_ == 0 || ny_ == 0 ||
            !(dx_ > 0.0) || !(dy_ > 0.0) || !std::isfinite(planeZ_) ||
            !std::isfinite(xFirst_) || !std::isfinite(yFirst_))
          throw std::runtime_error("Invalid field geometry: " + path_);
        for (std::size_t index = 0; index < times_.size(); ++index)
          if (!std::isfinite(times_[index]) ||
              (index > 0 && !(times_[index] > times_[index - 1])))
            throw std::runtime_error(
              "Field times must be finite and strictly increasing: " + path_);
      }

      void close()
      {
        if (fieldType_ >= 0) H5Tclose(fieldType_);
        if (reconstructedElectric_ >= 0) H5Dclose(reconstructedElectric_);
        if (rawFields_ >= 0) H5Dclose(rawFields_);
        if (group_ >= 0) H5Gclose(group_);
        if (file_ >= 0) H5Fclose(file_);
        fieldType_ = reconstructedElectric_ = rawFields_ = group_ = file_ = -1;
      }

      std::string path_;
      hid_t file_;
      hid_t group_;
      hid_t rawFields_;
      hid_t reconstructedElectric_;
      hid_t fieldType_;
      bool reconstructed_;
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
      const double scale = std::max(1.0,
        std::max(std::abs(left), std::abs(right)));
      return std::abs(left - right) <=
        256.0 * std::numeric_limits<double>::epsilon() * scale;
    }

    void validateMatched(const FieldInput& field, const FieldInput& baseline)
    {
      if (field.samples() != baseline.samples() ||
          field.nx() != baseline.nx() || field.ny() != baseline.ny() ||
          !nearlyEqual(field.planeZ(), baseline.planeZ()) ||
          !nearlyEqual(field.xFirst(), baseline.xFirst()) ||
          !nearlyEqual(field.yFirst(), baseline.yFirst()) ||
          !nearlyEqual(field.dx(), baseline.dx()) ||
          !nearlyEqual(field.dy(), baseline.dy()))
        throw std::runtime_error(
          "Signal and zero-radiation baseline geometries differ");
      for (std::size_t index = 0; index < field.samples(); ++index)
        if (!nearlyEqual(field.times()[index], baseline.times()[index]))
          throw std::runtime_error(
            "Signal and zero-radiation baseline times differ");
    }

    struct TimeWindow
    {
      std::size_t start;
      double center;
    };

    struct TimeModel
    {
      double start;
      double dt;
      double maximumIrregularity;
      std::size_t segmentSamples;
      std::vector<double> uniformTimes;
      std::vector<std::size_t> interpolationLower;
      std::vector<double> interpolationFraction;
      std::vector<double> windowWeight;
      std::vector<TimeWindow> windows;
    };

    TimeModel makeTimeModel(const FieldInput& field,
                            const Configuration& config)
    {
      TimeModel model;
      model.start = field.times().front();
      model.dt = (field.times().back() - field.times().front()) /
        static_cast<double>(field.samples() - 1);
      if (!(model.dt > 0.0) || !std::isfinite(model.dt))
        throw std::runtime_error("Cannot derive a positive uniform time step");
      model.uniformTimes.resize(field.samples());
      model.interpolationLower.resize(field.samples());
      model.interpolationFraction.resize(field.samples());
      model.maximumIrregularity = 0.0;
      std::size_t upper = 1;
      for (std::size_t sample = 0; sample < field.samples(); ++sample)
        {
          const double target = std::fma(static_cast<double>(sample),
            model.dt, model.start);
          model.uniformTimes[sample] = target;
          model.maximumIrregularity = std::max(model.maximumIrregularity,
            std::abs(field.times()[sample] - target) / model.dt);
          while (upper + 1 < field.samples() &&
                 field.times()[upper] < target)
            ++upper;
          if (target <= field.times().front())
            {
              model.interpolationLower[sample] = 0;
              model.interpolationFraction[sample] = 0.0;
            }
          else if (target >= field.times().back())
            {
              model.interpolationLower[sample] = field.samples() - 2;
              model.interpolationFraction[sample] = 1.0;
            }
          else
            {
              const std::size_t lower = upper - 1;
              model.interpolationLower[sample] = lower;
              model.interpolationFraction[sample] =
                (target - field.times()[lower]) /
                (field.times()[upper] - field.times()[lower]);
            }
        }

      if (!config.hannTimeWindows)
        {
          model.segmentSamples = field.samples();
          TimeWindow window;
          window.start = 0;
          window.center = 0.5 *
            (model.uniformTimes.front() + model.uniformTimes.back());
          model.windows.push_back(window);
          model.windowWeight.assign(model.segmentSamples, 1.0);
          return model;
        }

      model.segmentSamples = static_cast<std::size_t>(std::llround(
        config.windowDuration / model.dt));
      const std::size_t stepSamples = static_cast<std::size_t>(std::llround(
        config.windowStep / model.dt));
      if (model.segmentSamples < 4 || model.segmentSamples > field.samples() ||
          stepSamples == 0)
        throw std::runtime_error(
          "Hann duration/step quantizes to an invalid sample count; increase "
          "the duration or step, or refine detector cadence");
      const long long first = static_cast<long long>(std::ceil(
        (config.intervalStart - model.start) / model.dt - 1.0e-12));
      const long long lastExclusive = static_cast<long long>(std::floor(
        (config.intervalEnd - model.start) / model.dt + 1.0e-12)) + 1;
      for (long long begin = std::max<long long>(0, first);
           begin + static_cast<long long>(model.segmentSamples) <=
             std::min<long long>(static_cast<long long>(field.samples()),
                                 lastExclusive);
           begin += static_cast<long long>(stepSamples))
        {
          TimeWindow window;
          window.start = static_cast<std::size_t>(begin);
          window.center = model.start +
            (static_cast<double>(window.start) +
             0.5 * static_cast<double>(model.segmentSamples - 1)) * model.dt;
          model.windows.push_back(window);
        }
      if (model.windows.empty())
        throw std::runtime_error(
          "No complete Hann window fits in the configured interval");

      model.windowWeight.resize(model.segmentSamples);
      long double squared = 0.0L;
      for (std::size_t sample = 0; sample < model.segmentSamples; ++sample)
        {
          const double value = 0.5 * (1.0 - std::cos(
            2.0 * kPi * static_cast<double>(sample) /
            static_cast<double>(model.segmentSamples - 1)));
          model.windowWeight[sample] = value;
          squared += static_cast<long double>(value) * value;
        }
      const double correction = std::sqrt(
        static_cast<double>(model.segmentSamples) /
        static_cast<double>(squared));
      for (std::size_t sample = 0; sample < model.windowWeight.size(); ++sample)
        model.windowWeight[sample] *= correction;
      return model;
    }

    struct FrequencyModel
    {
      std::vector<std::size_t> bins;
      std::vector<double> frequency;
      std::vector<double> photonEnergy;
      double binWidthEV;
      double nyquistEV;
    };

    FrequencyModel makeFrequencyModel(const TimeModel& time,
                                      const Configuration& config)
    {
      FrequencyModel model;
      const double denominator =
        static_cast<double>(time.segmentSamples) * time.dt;
      model.binWidthEV = kPlanckEVSecond / denominator;
      model.nyquistEV = kPlanckEVSecond * 0.5 / time.dt;
      if (config.maximumPhotonEnergyEV > model.nyquistEV *
          (1.0 + 64.0 * std::numeric_limits<double>::epsilon()))
        throw std::runtime_error(
          "Requested photon-energy maximum exceeds the field-plane Nyquist "
          "energy; decrease detector rhythm/Maxwell step or lower the band");
      for (std::size_t bin = 1; bin <= time.segmentSamples / 2; ++bin)
        {
          const double frequency = static_cast<double>(bin) / denominator;
          const double energy = kPlanckEVSecond * frequency;
          if (energy >= config.minimumPhotonEnergyEV &&
              energy <= config.maximumPhotonEnergyEV)
            {
              model.bins.push_back(bin);
              model.frequency.push_back(frequency);
              model.photonEnergy.push_back(energy);
            }
        }
      if (model.bins.empty())
        throw std::runtime_error(
          "Requested photon-energy band contains no positive FFT bins");
      return model;
    }

    std::vector<std::size_t> mapEnergies(
        const std::vector<double>& requested,
        const FrequencyModel& frequencies,
        const char* description)
    {
      std::vector<std::size_t> result;
      for (std::size_t request = 0; request < requested.size(); ++request)
        {
          if (requested[request] < frequencies.photonEnergy.front() -
                0.5 * frequencies.binWidthEV ||
              requested[request] > frequencies.photonEnergy.back() +
                0.5 * frequencies.binWidthEV)
            throw std::runtime_error(std::string(description) +
              " lies outside the retained FFT band");
          std::size_t nearest = 0;
          double distance = std::abs(requested[request] -
            frequencies.photonEnergy[0]);
          for (std::size_t index = 1;
               index < frequencies.photonEnergy.size(); ++index)
            {
              const double candidate = std::abs(requested[request] -
                frequencies.photonEnergy[index]);
              if (candidate < distance)
                {
                  distance = candidate;
                  nearest = index;
                }
            }
          result.push_back(nearest);
        }
      return result;
    }

    struct ComplexValue
    {
      double real;
      double imaginary;
    };

    hid_t createComplexType()
    {
      hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(ComplexValue));
      requireHandle(type, "Cannot create complex datatype");
      requireStatus(H5Tinsert(type, "real", HOFFSET(ComplexValue, real),
        H5T_IEEE_F64LE), "Cannot add complex real member");
      requireStatus(H5Tinsert(type, "imag", HOFFSET(ComplexValue, imaginary),
        H5T_IEEE_F64LE), "Cannot add complex imaginary member");
      return type;
    }

    void writeStringAttribute(hid_t object, const char* name,
                              const std::string& value)
    {
      hid_t type = H5Tcopy(H5T_C_S1);
      requireHandle(type, "Cannot create string datatype");
      requireStatus(H5Tset_size(type, value.size() + 1),
        "Cannot size string datatype");
      hid_t space = H5Screate(H5S_SCALAR);
      requireHandle(space, "Cannot create string attribute space");
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
      requireHandle(space, "Cannot create scalar attribute space");
      hid_t attribute = H5Acreate2(object, name, H5T_IEEE_F64LE, space,
        H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(space);
      requireHandle(attribute, std::string("Cannot create attribute: ") + name);
      const herr_t status = H5Awrite(attribute, H5T_NATIVE_DOUBLE, &value);
      H5Aclose(attribute);
      requireStatus(status, std::string("Cannot write attribute: ") + name);
    }

    void writeUnsignedAttribute(hid_t object, const char* name,
                                std::uint64_t value)
    {
      hid_t space = H5Screate(H5S_SCALAR);
      requireHandle(space, "Cannot create scalar attribute space");
      hid_t attribute = H5Acreate2(object, name, H5T_STD_U64LE, space,
        H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(space);
      requireHandle(attribute, std::string("Cannot create attribute: ") + name);
      const herr_t status = H5Awrite(attribute, H5T_NATIVE_UINT64, &value);
      H5Aclose(attribute);
      requireStatus(status, std::string("Cannot write attribute: ") + name);
    }

    hid_t createDataset(hid_t group, const char* name, hid_t type,
                        const std::vector<hsize_t>& dimensions,
                        const std::vector<hsize_t>& chunks,
                        unsigned int compression)
    {
      hid_t space = H5Screate_simple(static_cast<int>(dimensions.size()),
        &dimensions[0], NULL);
      requireHandle(space, std::string("Cannot create dataset space: ") + name);
      hid_t creation = H5P_DEFAULT;
      if (compression > 0)
        {
          creation = H5Pcreate(H5P_DATASET_CREATE);
          requireHandle(creation,
            std::string("Cannot create dataset properties: ") + name);
          requireStatus(H5Pset_chunk(creation, static_cast<int>(chunks.size()),
            &chunks[0]), std::string("Cannot set chunks: ") + name);
          requireStatus(H5Pset_deflate(creation, compression),
            std::string("Cannot set compression: ") + name);
        }
      hid_t dataset = H5Dcreate2(group, name, type, space, H5P_DEFAULT,
        creation, H5P_DEFAULT);
      if (creation != H5P_DEFAULT) H5Pclose(creation);
      H5Sclose(space);
      requireHandle(dataset, std::string("Cannot create dataset: ") + name);
      return dataset;
    }

    void writeVector(hid_t group, const char* name,
                     const std::vector<double>& values, const char* unit)
    {
      hsize_t dimension[1] = {static_cast<hsize_t>(values.size())};
      hid_t space = H5Screate_simple(1, dimension, NULL);
      requireHandle(space, std::string("Cannot create axis space: ") + name);
      hid_t dataset = H5Dcreate2(group, name, H5T_IEEE_F64LE, space,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(space);
      requireHandle(dataset, std::string("Cannot create axis: ") + name);
      requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, values.empty() ? NULL : &values[0]),
        std::string("Cannot write axis: ") + name);
      if (unit != NULL) writeStringAttribute(dataset, "unit", unit);
      H5Dclose(dataset);
    }

    void writeHyperslab(hid_t dataset, hid_t memoryType,
                        const std::vector<hsize_t>& start,
                        const std::vector<hsize_t>& count, const void* data,
                        const char* description)
    {
      hid_t fileSpace = H5Dget_space(dataset);
      requireHandle(fileSpace, std::string("Cannot inspect ") + description);
      requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
        &start[0], NULL, &count[0], NULL),
        std::string("Cannot select ") + description);
      hid_t memorySpace = H5Screate_simple(static_cast<int>(count.size()),
        &count[0], NULL);
      requireHandle(memorySpace,
        std::string("Cannot create memory space for ") + description);
      const herr_t status = H5Dwrite(dataset, memoryType, memorySpace,
        fileSpace, H5P_DEFAULT, data);
      H5Sclose(memorySpace);
      H5Sclose(fileSpace);
      requireStatus(status, std::string("Cannot write ") + description);
    }

    std::vector<double> requestedAngles(
        const std::vector<ReferenceAngle>& angles)
    {
      std::vector<double> result(angles.size() * 2);
      for (std::size_t index = 0; index < angles.size(); ++index)
        {
          result[2 * index] = angles[index].thetaX;
          result[2 * index + 1] = angles[index].thetaY;
        }
      return result;
    }

    class OutputFile
    {
    public:
      OutputFile(const Configuration& config, const FieldInput& field,
                 const TimeModel& time, const FrequencyModel& frequencies,
                 std::size_t ny, std::size_t nx,
                 const std::vector<double>& ky,
                 const std::vector<double>& kx,
                 const std::vector<std::size_t>& spatialFrequency,
                 const std::vector<std::size_t>& temporalFrequency,
                 double workingMiB, double outputGiB)
        : file_(-1), group_(-1), spatialGroup_(-1), temporalGroup_(-1),
          complexType_(-1), complete_(-1), meanField_(-1), meanIntensity_(-1),
          coherentIntensity_(-1), fluctuationIntensity_(-1), stokes_(-1),
          solidAngle_(-1), meanSpectrum_(-1), sampleSpectrum_(-1),
          coherentFraction_(-1), ensembleGram_(-1), globalCoherence_(-1),
          spatialCsd_(-1), spatialMu2_(-1),
          spatialActualAngles_(-1), spatialActualIndices_(-1),
          temporalCsd_(-1), temporalMu2_(-1), temporalActualAngles_(-1),
          config_(config), ny_(ny), nx_(nx), samples_(time.windows.size()),
          frequencies_(frequencies.photonEnergy.size()), open_(false)
      {
        createDirectories(parentDirectory(config.outputFile));
        file_ = H5Fcreate(config.outputFile.c_str(), H5F_ACC_TRUNC,
          H5P_DEFAULT, H5P_DEFAULT);
        requireHandle(file_, "Cannot create field-plane analysis output: " +
          config.outputFile);
        try
          {
            group_ = H5Gcreate2(file_, "/field_plane_analysis", H5P_DEFAULT,
              H5P_DEFAULT, H5P_DEFAULT);
            requireHandle(group_, "Cannot create /field_plane_analysis");
            complexType_ = createComplexType();
            hid_t scalar = H5Screate(H5S_SCALAR);
            complete_ = H5Dcreate2(group_, "complete", H5T_STD_U8LE, scalar,
              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(scalar);
            requireHandle(complete_, "Cannot create completion marker");
            const unsigned char incomplete = 0;
            requireStatus(H5Dwrite(complete_, H5T_NATIVE_UCHAR,
              H5S_ALL, H5S_ALL, H5P_DEFAULT, &incomplete),
              "Cannot initialize completion marker");

            writeUnsignedAttribute(group_, "format_version", 1);
            writeStringAttribute(group_, "input_field_file", field.path());
            writeStringAttribute(group_, "input_field_kind",
              field.reconstructed() ? "reconstructed_field" : "field_plane");
            writeStringAttribute(group_, "zero_radiation_baseline",
              config.baselineFile.empty() ? "none" : config.baselineFile);
            writeStringAttribute(group_, "radiation_model",
              "forward propagating vacuum angular spectrum from transverse projection of the laboratory electric field");
            writeStringAttribute(group_, "time_resampling",
              "linear interpolation onto a uniform grid spanning the committed input interval");
            writeStringAttribute(group_, "time_ensemble",
              config.hannTimeWindows ?
              "overlapping Hann windows; sqrt(N/sum(w^2)) energy normalization" :
              "one rectangular full-record sample");
            writeStringAttribute(group_, "transverse_window",
              config.transverseWindow);
            writeStringAttribute(group_, "angular_coordinates",
              "kx,ky are exact FFT axes; for each photon energy nx=kx/k, ny=ky/k, nz=sqrt(1-nx^2-ny^2), theta_x=atan2(nx,nz), theta_y=atan2(ny,nz); zero solid-angle weight marks nonpropagating bins");
            writeStringAttribute(group_, "stokes_convention",
              "I=|Eh|^2+|Ev|^2,Q=|Eh|^2-|Ev|^2,U=2Re(Eh conj(Ev)),V=-2Im(Eh conj(Ev)); scaled to J/eV/sr");
            writeStringAttribute(group_, "coherence_contract",
              "CSD is an ensemble mean over configured time windows; one sample is necessarily rank one");
            writeDoubleAttribute(group_, "plane_z_m", field.planeZ());
            writeDoubleAttribute(group_, "x_first_m", field.xFirst());
            writeDoubleAttribute(group_, "y_first_m", field.yFirst());
            writeDoubleAttribute(group_, "dx_m", field.dx());
            writeDoubleAttribute(group_, "dy_m", field.dy());
            writeDoubleAttribute(group_, "uniform_dt_s", time.dt);
            writeDoubleAttribute(group_, "maximum_time_irregularity_steps",
              time.maximumIrregularity);
            writeDoubleAttribute(group_, "frequency_bin_width_eV",
              frequencies.binWidthEV);
            writeDoubleAttribute(group_, "band_min_photon_energy_eV",
              config.minimumPhotonEnergyEV);
            writeDoubleAttribute(group_, "band_max_photon_energy_eV",
              config.maximumPhotonEnergyEV);
            writeDoubleAttribute(group_, "nyquist_photon_energy_eV",
              frequencies.nyquistEV);
            writeDoubleAttribute(group_, "estimated_working_mib", workingMiB);
            writeDoubleAttribute(group_, "estimated_output_gib", outputGiB);
            writeUnsignedAttribute(group_, "input_samples", field.samples());
            writeUnsignedAttribute(group_, "segment_samples",
              time.segmentSamples);
            writeUnsignedAttribute(group_, "ensemble_samples", samples_);
            writeUnsignedAttribute(group_, "ny", ny_);
            writeUnsignedAttribute(group_, "nx", nx_);
            writeUnsignedAttribute(group_, "angular_padding_y",
              config.angularPaddingY);
            writeUnsignedAttribute(group_, "angular_padding_x",
              config.angularPaddingX);

            writeVector(group_, "photon_energy_eV",
              frequencies.photonEnergy, "eV");
            writeVector(group_, "frequency_Hz", frequencies.frequency, "Hz");
            std::vector<double> omega(frequencies.frequency.size());
            std::vector<double> wavelength(frequencies.frequency.size());
            for (std::size_t index = 0; index < omega.size(); ++index)
              {
                omega[index] = 2.0 * kPi * frequencies.frequency[index];
                wavelength[index] = kC / frequencies.frequency[index];
              }
            writeVector(group_, "omega_rad_per_s", omega, "rad/s");
            writeVector(group_, "wavelength_m", wavelength, "m");
            writeVector(group_, "kx_rad_per_m", kx, "rad/m");
            writeVector(group_, "ky_rad_per_m", ky, "rad/m");
            std::vector<double> centers(time.windows.size());
            for (std::size_t index = 0; index < centers.size(); ++index)
              centers[index] = time.windows[index].center;
            writeVector(group_, "sample_center_time_s", centers, "s");

            const std::vector<hsize_t> scalarDimensions{
              static_cast<hsize_t>(frequencies_),
              static_cast<hsize_t>(ny_), static_cast<hsize_t>(nx_)};
            const std::vector<hsize_t> scalarChunks{1,
              static_cast<hsize_t>(std::min<std::size_t>(16, ny_)),
              static_cast<hsize_t>(std::min<std::size_t>(64, nx_))};
            std::vector<hsize_t> fieldDimensions(scalarDimensions);
            fieldDimensions.push_back(2);
            std::vector<hsize_t> fieldChunks(scalarChunks);
            fieldChunks.push_back(2);
            std::vector<hsize_t> stokesDimensions(scalarDimensions);
            stokesDimensions.push_back(4);
            std::vector<hsize_t> stokesChunks(scalarChunks);
            stokesChunks.push_back(4);
            meanField_ = createDataset(group_,
              "mean_angular_electric_field_spectral", complexType_,
              fieldDimensions, fieldChunks, config.compression);
            writeStringAttribute(meanField_, "unit", "V s m");
            meanIntensity_ = createDataset(group_,
              "mean_angular_spectral_energy_density_J_per_eV_sr",
              H5T_IEEE_F64LE, scalarDimensions, scalarChunks,
              config.compression);
            writeStringAttribute(meanIntensity_, "unit", "J/eV/sr");
            coherentIntensity_ = createDataset(group_,
              "coherent_angular_spectral_energy_density_J_per_eV_sr",
              H5T_IEEE_F64LE, scalarDimensions, scalarChunks,
              config.compression);
            writeStringAttribute(coherentIntensity_, "unit", "J/eV/sr");
            fluctuationIntensity_ = createDataset(group_,
              "fluctuation_angular_spectral_energy_density_J_per_eV_sr",
              H5T_IEEE_F64LE, scalarDimensions, scalarChunks,
              config.compression);
            writeStringAttribute(fluctuationIntensity_, "unit", "J/eV/sr");
            stokes_ = createDataset(group_,
              "mean_stokes_spectral_energy_density_J_per_eV_sr",
              H5T_IEEE_F64LE, stokesDimensions, stokesChunks,
              config.compression);
            writeStringAttribute(stokes_, "unit", "J/eV/sr");
            solidAngle_ = createDataset(group_, "solid_angle_weight_sr",
              H5T_IEEE_F64LE, scalarDimensions, scalarChunks,
              config.compression);
            writeStringAttribute(solidAngle_, "unit", "sr");

            const std::vector<hsize_t> spectrumDimensions{
              static_cast<hsize_t>(frequencies_)};
            const std::vector<hsize_t> spectrumChunks{
              static_cast<hsize_t>(std::min<std::size_t>(256, frequencies_))};
            meanSpectrum_ = createDataset(group_,
              "mean_energy_spectrum_J_per_eV", H5T_IEEE_F64LE,
              spectrumDimensions, spectrumChunks, config.compression);
            writeStringAttribute(meanSpectrum_, "unit", "J/eV");
            coherentFraction_ = createDataset(group_,
              "coherent_fraction_spectrum", H5T_IEEE_F64LE,
              spectrumDimensions, spectrumChunks, config.compression);
            globalCoherence_ = createDataset(group_,
              "global_degree_of_transverse_coherence", H5T_IEEE_F64LE,
              spectrumDimensions, spectrumChunks, config.compression);
            const std::vector<hsize_t> sampleSpectrumDimensions{
              static_cast<hsize_t>(samples_),
              static_cast<hsize_t>(frequencies_)};
            const std::vector<hsize_t> sampleSpectrumChunks{1,
              static_cast<hsize_t>(std::min<std::size_t>(256, frequencies_))};
            sampleSpectrum_ = createDataset(group_,
              "sample_energy_spectrum_J_per_eV", H5T_IEEE_F64LE,
              sampleSpectrumDimensions, sampleSpectrumChunks,
              config.compression);
            writeStringAttribute(sampleSpectrum_, "unit", "J/eV");
            const std::vector<hsize_t> gramDimensions{
              static_cast<hsize_t>(frequencies_),
              static_cast<hsize_t>(samples_),
              static_cast<hsize_t>(samples_)};
            const std::vector<hsize_t> gramChunks{1,
              static_cast<hsize_t>(std::min<std::size_t>(32, samples_)),
              static_cast<hsize_t>(std::min<std::size_t>(32, samples_))};
            ensembleGram_ = createDataset(group_,
              "ensemble_gram_matrix_J_per_eV", complexType_, gramDimensions,
              gramChunks, config.compression);
            writeStringAttribute(ensembleGram_, "unit", "J/eV");
            writeStringAttribute(ensembleGram_, "eigenvalue_contract",
              "nonzero eigenvalues equal those of the angular-polarization CSD; matrix includes 1/ensemble_samples normalization and solid-angle quadrature");

            if (!spatialFrequency.empty())
              createSpatial(config, frequencies, spatialFrequency);
            if (!temporalFrequency.empty())
              createTemporal(config, frequencies, temporalFrequency);
            open_ = true;
          }
        catch (...)
          {
            close(false);
            throw;
          }
      }

      ~OutputFile() { close(false); }

      void writeFrequency(std::size_t frequency,
                          const std::vector<ComplexValue>& meanField,
                          const std::vector<double>& meanIntensity,
                          const std::vector<double>& coherent,
                          const std::vector<double>& fluctuation,
                          const std::vector<double>& stokes,
                          const std::vector<double>& solidAngle,
                          const std::vector<double>& sampleSpectrum,
                          const std::vector<ComplexValue>& ensembleGram,
                          double meanSpectrum, double coherentFraction,
                          double globalCoherence)
      {
        const std::vector<hsize_t> scalarStart{
          static_cast<hsize_t>(frequency), 0, 0};
        const std::vector<hsize_t> scalarCount{1,
          static_cast<hsize_t>(ny_), static_cast<hsize_t>(nx_)};
        writeHyperslab(meanIntensity_, H5T_NATIVE_DOUBLE, scalarStart,
          scalarCount, &meanIntensity[0], "mean angular intensity");
        writeHyperslab(coherentIntensity_, H5T_NATIVE_DOUBLE, scalarStart,
          scalarCount, &coherent[0], "coherent angular intensity");
        writeHyperslab(fluctuationIntensity_, H5T_NATIVE_DOUBLE, scalarStart,
          scalarCount, &fluctuation[0], "fluctuation angular intensity");
        writeHyperslab(solidAngle_, H5T_NATIVE_DOUBLE, scalarStart,
          scalarCount, &solidAngle[0], "solid-angle weights");
        std::vector<hsize_t> fieldStart(scalarStart);
        fieldStart.push_back(0);
        std::vector<hsize_t> fieldCount(scalarCount);
        fieldCount.push_back(2);
        writeHyperslab(meanField_, complexType_, fieldStart, fieldCount,
          &meanField[0], "mean angular field");
        std::vector<hsize_t> stokesStart(scalarStart);
        stokesStart.push_back(0);
        std::vector<hsize_t> stokesCount(scalarCount);
        stokesCount.push_back(4);
        writeHyperslab(stokes_, H5T_NATIVE_DOUBLE, stokesStart, stokesCount,
          &stokes[0], "mean Stokes spectrum");
        const std::vector<hsize_t> scalarAxisStart{
          static_cast<hsize_t>(frequency)};
        const std::vector<hsize_t> scalarAxisCount{1};
        writeHyperslab(meanSpectrum_, H5T_NATIVE_DOUBLE, scalarAxisStart,
          scalarAxisCount, &meanSpectrum, "mean energy spectrum");
        writeHyperslab(coherentFraction_, H5T_NATIVE_DOUBLE, scalarAxisStart,
          scalarAxisCount, &coherentFraction,
          "coherent fraction spectrum");
        writeHyperslab(globalCoherence_, H5T_NATIVE_DOUBLE, scalarAxisStart,
          scalarAxisCount, &globalCoherence,
          "global transverse coherence");
        const std::vector<hsize_t> sampleStart{0,
          static_cast<hsize_t>(frequency)};
        const std::vector<hsize_t> sampleCount{
          static_cast<hsize_t>(samples_), 1};
        writeHyperslab(sampleSpectrum_, H5T_NATIVE_DOUBLE, sampleStart,
          sampleCount, &sampleSpectrum[0], "sample energy spectrum");
        const std::vector<hsize_t> gramStart{
          static_cast<hsize_t>(frequency), 0, 0};
        const std::vector<hsize_t> gramCount{1,
          static_cast<hsize_t>(samples_), static_cast<hsize_t>(samples_)};
        writeHyperslab(ensembleGram_, complexType_, gramStart, gramCount,
          &ensembleGram[0], "ensemble Gram matrix");
      }

      void writeSpatial(std::size_t selected,
                        const std::vector<double>& actualAngles,
                        const std::vector<std::uint64_t>& actualIndices,
                        const std::vector<ComplexValue>& csd,
                        const std::vector<double>& mu2,
                        std::size_t references)
      {
        const std::vector<hsize_t> angleStart{
          static_cast<hsize_t>(selected), 0, 0};
        const std::vector<hsize_t> angleCount{1,
          static_cast<hsize_t>(references), 2};
        writeHyperslab(spatialActualAngles_, H5T_NATIVE_DOUBLE, angleStart,
          angleCount, &actualAngles[0], "spatial actual reference angles");
        writeHyperslab(spatialActualIndices_, H5T_NATIVE_UINT64, angleStart,
          angleCount, &actualIndices[0], "spatial actual reference indices");
        const std::vector<hsize_t> muStart{
          static_cast<hsize_t>(selected), 0, 0, 0};
        const std::vector<hsize_t> muCount{1,
          static_cast<hsize_t>(references), static_cast<hsize_t>(ny_),
          static_cast<hsize_t>(nx_)};
        writeHyperslab(spatialMu2_, H5T_NATIVE_DOUBLE, muStart, muCount,
          &mu2[0], "spatial degree of coherence");
        std::vector<hsize_t> csdStart(muStart);
        csdStart.push_back(0);
        csdStart.push_back(0);
        std::vector<hsize_t> csdCount(muCount);
        csdCount.push_back(2);
        csdCount.push_back(2);
        writeHyperslab(spatialCsd_, complexType_, csdStart, csdCount,
          &csd[0], "spatial cross-spectral density");
      }

      void writeTemporal(const std::vector<double>& actualAngles,
                         const std::vector<ComplexValue>& csd,
                         const std::vector<double>& mu2,
                         std::size_t references,
                         std::size_t selectedFrequencies)
      {
        const std::vector<hsize_t> angleStart{0, 0, 0};
        const std::vector<hsize_t> angleCount{
          static_cast<hsize_t>(selectedFrequencies),
          static_cast<hsize_t>(references), 2};
        writeHyperslab(temporalActualAngles_, H5T_NATIVE_DOUBLE,
          angleStart, angleCount, &actualAngles[0],
          "temporal actual reference angles");
        const std::vector<hsize_t> muStart{0, 0, 0};
        const std::vector<hsize_t> muCount{
          static_cast<hsize_t>(references),
          static_cast<hsize_t>(selectedFrequencies),
          static_cast<hsize_t>(selectedFrequencies)};
        writeHyperslab(temporalMu2_, H5T_NATIVE_DOUBLE, muStart, muCount,
          &mu2[0], "temporal degree of coherence");
        const std::vector<hsize_t> csdStart{0, 0, 0, 0, 0};
        const std::vector<hsize_t> csdCount{
          static_cast<hsize_t>(references),
          static_cast<hsize_t>(selectedFrequencies),
          static_cast<hsize_t>(selectedFrequencies), 2, 2};
        writeHyperslab(temporalCsd_, complexType_, csdStart, csdCount,
          &csd[0], "two-frequency cross-spectral density");
      }

      void finish(const std::vector<double>& sampleBandEnergy,
                  double meanBandEnergy, double bandCoherentFraction)
      {
        writeVector(group_, "sample_band_energy_J", sampleBandEnergy, "J");
        writeDoubleAttribute(group_, "mean_band_energy_J", meanBandEnergy);
        writeDoubleAttribute(group_, "band_coherent_fraction",
          bandCoherentFraction);
        close(true);
      }

    private:
      OutputFile(const OutputFile&);
      OutputFile& operator=(const OutputFile&);

      void createSpatial(const Configuration& config,
                         const FrequencyModel& frequencies,
                         const std::vector<std::size_t>& selected)
      {
        spatialGroup_ = H5Gcreate2(file_, "/spatial_coherence", H5P_DEFAULT,
          H5P_DEFAULT, H5P_DEFAULT);
        requireHandle(spatialGroup_, "Cannot create /spatial_coherence");
        std::vector<double> actualEnergy(selected.size());
        for (std::size_t index = 0; index < selected.size(); ++index)
          actualEnergy[index] = frequencies.photonEnergy[selected[index]];
        writeVector(spatialGroup_, "requested_photon_energy_eV",
          config.spatialPhotonEnergyEV, "eV");
        writeVector(spatialGroup_, "actual_photon_energy_eV", actualEnergy,
          "eV");
        writeAngleDataset(spatialGroup_, "requested_reference_angles_rad",
          requestedAngles(config.spatialReferenceAngles),
          config.spatialReferenceAngles.size());
        const std::size_t references = config.spatialReferenceAngles.size();
        const std::vector<hsize_t> angleDimensions{
          static_cast<hsize_t>(selected.size()),
          static_cast<hsize_t>(references), 2};
        const std::vector<hsize_t> angleChunks{1,
          static_cast<hsize_t>(references), 2};
        spatialActualAngles_ = createDataset(spatialGroup_,
          "actual_reference_angles_rad", H5T_IEEE_F64LE,
          angleDimensions, angleChunks, config.compression);
        spatialActualIndices_ = createDataset(spatialGroup_,
          "actual_reference_indices_y_x", H5T_STD_U64LE,
          angleDimensions, angleChunks, config.compression);
        const std::vector<hsize_t> muDimensions{
          static_cast<hsize_t>(selected.size()),
          static_cast<hsize_t>(references), static_cast<hsize_t>(ny_),
          static_cast<hsize_t>(nx_)};
        const std::vector<hsize_t> muChunks{1, 1,
          static_cast<hsize_t>(std::min<std::size_t>(8, ny_)),
          static_cast<hsize_t>(std::min<std::size_t>(32, nx_))};
        spatialMu2_ = createDataset(spatialGroup_,
          "spectral_degree_of_coherence_squared", H5T_IEEE_F64LE,
          muDimensions, muChunks, config.compression);
        std::vector<hsize_t> csdDimensions(muDimensions);
        csdDimensions.push_back(2);
        csdDimensions.push_back(2);
        std::vector<hsize_t> csdChunks(muChunks);
        csdChunks.push_back(2);
        csdChunks.push_back(2);
        spatialCsd_ = createDataset(spatialGroup_,
          "angular_cross_spectral_density_J_per_eV_sr", complexType_,
          csdDimensions, csdChunks, config.compression);
        writeStringAttribute(spatialCsd_, "unit", "J/eV/sr");
      }

      void createTemporal(const Configuration& config,
                          const FrequencyModel& frequencies,
                          const std::vector<std::size_t>& selected)
      {
        temporalGroup_ = H5Gcreate2(file_, "/temporal_coherence", H5P_DEFAULT,
          H5P_DEFAULT, H5P_DEFAULT);
        requireHandle(temporalGroup_, "Cannot create /temporal_coherence");
        std::vector<double> actualEnergy(selected.size());
        for (std::size_t index = 0; index < selected.size(); ++index)
          actualEnergy[index] = frequencies.photonEnergy[selected[index]];
        writeVector(temporalGroup_, "requested_photon_energy_eV",
          config.temporalPhotonEnergyEV, "eV");
        writeVector(temporalGroup_, "actual_photon_energy_eV", actualEnergy,
          "eV");
        writeAngleDataset(temporalGroup_, "requested_reference_angles_rad",
          requestedAngles(config.temporalReferenceAngles),
          config.temporalReferenceAngles.size());
        const std::size_t references = config.temporalReferenceAngles.size();
        const std::size_t selectedCount = selected.size();
        const std::vector<hsize_t> angleDimensions{
          static_cast<hsize_t>(selectedCount),
          static_cast<hsize_t>(references), 2};
        const std::vector<hsize_t> angleChunks{1,
          static_cast<hsize_t>(references), 2};
        temporalActualAngles_ = createDataset(temporalGroup_,
          "actual_reference_angles_rad", H5T_IEEE_F64LE,
          angleDimensions, angleChunks, config.compression);
        const std::vector<hsize_t> muDimensions{
          static_cast<hsize_t>(references),
          static_cast<hsize_t>(selectedCount),
          static_cast<hsize_t>(selectedCount)};
        const std::vector<hsize_t> muChunks{1,
          static_cast<hsize_t>(std::min<std::size_t>(32, selectedCount)),
          static_cast<hsize_t>(std::min<std::size_t>(32, selectedCount))};
        temporalMu2_ = createDataset(temporalGroup_,
          "spectral_degree_of_coherence_squared", H5T_IEEE_F64LE,
          muDimensions, muChunks, config.compression);
        std::vector<hsize_t> csdDimensions(muDimensions);
        csdDimensions.push_back(2);
        csdDimensions.push_back(2);
        std::vector<hsize_t> csdChunks(muChunks);
        csdChunks.push_back(2);
        csdChunks.push_back(2);
        temporalCsd_ = createDataset(temporalGroup_,
          "two_frequency_cross_spectral_density_J_per_eV_sr", complexType_,
          csdDimensions, csdChunks, config.compression);
        writeStringAttribute(temporalCsd_, "unit", "J/eV/sr");
      }

      void writeAngleDataset(hid_t group, const char* name,
                             const std::vector<double>& values,
                             std::size_t references)
      {
        hsize_t dimensions[2] = {static_cast<hsize_t>(references), 2};
        hid_t space = H5Screate_simple(2, dimensions, NULL);
        requireHandle(space, "Cannot create reference-angle space");
        hid_t dataset = H5Dcreate2(group, name, H5T_IEEE_F64LE, space,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(space);
        requireHandle(dataset, "Cannot create reference-angle dataset");
        requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
          H5P_DEFAULT, values.empty() ? NULL : &values[0]),
          "Cannot write reference angles");
        writeStringAttribute(dataset, "unit", "rad");
        H5Dclose(dataset);
      }

      void close(bool completed)
      {
        if (file_ < 0) return;
        if (complete_ >= 0)
          {
            const unsigned char value = completed ? 1 : 0;
            H5Dwrite(complete_, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
              H5P_DEFAULT, &value);
          }
        H5Fflush(file_, H5F_SCOPE_GLOBAL);
        if (temporalActualAngles_ >= 0) H5Dclose(temporalActualAngles_);
        if (temporalMu2_ >= 0) H5Dclose(temporalMu2_);
        if (temporalCsd_ >= 0) H5Dclose(temporalCsd_);
        if (spatialActualIndices_ >= 0) H5Dclose(spatialActualIndices_);
        if (spatialActualAngles_ >= 0) H5Dclose(spatialActualAngles_);
        if (spatialMu2_ >= 0) H5Dclose(spatialMu2_);
        if (spatialCsd_ >= 0) H5Dclose(spatialCsd_);
        if (coherentFraction_ >= 0) H5Dclose(coherentFraction_);
        if (globalCoherence_ >= 0) H5Dclose(globalCoherence_);
        if (ensembleGram_ >= 0) H5Dclose(ensembleGram_);
        if (sampleSpectrum_ >= 0) H5Dclose(sampleSpectrum_);
        if (meanSpectrum_ >= 0) H5Dclose(meanSpectrum_);
        if (solidAngle_ >= 0) H5Dclose(solidAngle_);
        if (stokes_ >= 0) H5Dclose(stokes_);
        if (fluctuationIntensity_ >= 0) H5Dclose(fluctuationIntensity_);
        if (coherentIntensity_ >= 0) H5Dclose(coherentIntensity_);
        if (meanIntensity_ >= 0) H5Dclose(meanIntensity_);
        if (meanField_ >= 0) H5Dclose(meanField_);
        if (complete_ >= 0) H5Dclose(complete_);
        if (temporalGroup_ >= 0) H5Gclose(temporalGroup_);
        if (spatialGroup_ >= 0) H5Gclose(spatialGroup_);
        if (group_ >= 0) H5Gclose(group_);
        if (complexType_ >= 0) H5Tclose(complexType_);
        H5Fclose(file_);
        temporalActualAngles_ = temporalMu2_ = temporalCsd_ = -1;
        spatialActualIndices_ = spatialActualAngles_ = spatialMu2_ =
          spatialCsd_ = -1;
        ensembleGram_ = globalCoherence_ = coherentFraction_ =
          sampleSpectrum_ = meanSpectrum_ = solidAngle_ = stokes_ =
          fluctuationIntensity_ = coherentIntensity_ = meanIntensity_ =
          meanField_ = complete_ = -1;
        temporalGroup_ = spatialGroup_ = group_ = complexType_ = file_ = -1;
        open_ = false;
      }

      hid_t file_;
      hid_t group_;
      hid_t spatialGroup_;
      hid_t temporalGroup_;
      hid_t complexType_;
      hid_t complete_;
      hid_t meanField_;
      hid_t meanIntensity_;
      hid_t coherentIntensity_;
      hid_t fluctuationIntensity_;
      hid_t stokes_;
      hid_t solidAngle_;
      hid_t meanSpectrum_;
      hid_t sampleSpectrum_;
      hid_t coherentFraction_;
      hid_t ensembleGram_;
      hid_t globalCoherence_;
      hid_t spatialCsd_;
      hid_t spatialMu2_;
      hid_t spatialActualAngles_;
      hid_t spatialActualIndices_;
      hid_t temporalCsd_;
      hid_t temporalMu2_;
      hid_t temporalActualAngles_;
      Configuration config_;
      std::size_t ny_;
      std::size_t nx_;
      std::size_t samples_;
      std::size_t frequencies_;
      bool open_;
    };

    struct ProjectionGeometry
    {
      std::vector<std::size_t> sourceIndex;
      std::vector<double> horizontal;
      std::vector<double> vertical;
      std::vector<double> factor;
      std::vector<double> solidAngle;
      std::vector<double> thetaX;
      std::vector<double> thetaY;
    };

    ProjectionGeometry makeProjectionGeometry(
        double frequency, std::size_t temporalBin,
        std::size_t segmentSamples,
        std::size_t ny, std::size_t nx,
        const std::vector<double>& ky, const std::vector<double>& kx)
    {
      ProjectionGeometry geometry;
      const std::size_t points = checkedProduct(ny, nx,
        "projection geometry");
      geometry.sourceIndex.resize(points);
      geometry.horizontal.assign(checkedProduct(points, 3,
        "horizontal projection basis"), 0.0);
      geometry.vertical.assign(checkedProduct(points, 3,
        "vertical projection basis"), 0.0);
      geometry.factor.assign(points, 0.0);
      geometry.solidAngle.assign(points, 0.0);
      geometry.thetaX.assign(points, 0.0);
      geometry.thetaY.assign(points, 0.0);
      const double waveNumber = 2.0 * kPi * frequency / kC;
      const double deltaKx = nx > 1 ? std::abs(kx[1] - kx[0]) : 0.0;
      const double deltaKy = ny > 1 ? std::abs(ky[1] - ky[0]) : 0.0;
      const bool selfConjugate = temporalBin == 0 ||
        (segmentSamples % 2 == 0 && temporalBin == segmentSamples / 2);
      const double oneSided = selfConjugate ? 1.0 : 2.0;
      for (std::size_t y = 0; y < ny; ++y)
        {
          const long long qy = static_cast<long long>(y) -
            static_cast<long long>(ny / 2);
          const std::size_t sourceY = static_cast<std::size_t>(
            (qy + static_cast<long long>(ny)) % static_cast<long long>(ny));
          for (std::size_t x = 0; x < nx; ++x)
            {
              const long long qx = static_cast<long long>(x) -
                static_cast<long long>(nx / 2);
              const std::size_t sourceX = static_cast<std::size_t>(
                (qx + static_cast<long long>(nx)) %
                static_cast<long long>(nx));
              const std::size_t point = y * nx + x;
              geometry.sourceIndex[point] = sourceY * nx + sourceX;
              const double directionX = kx[x] / waveNumber;
              const double directionY = ky[y] / waveNumber;
              const double transverse2 = directionX * directionX +
                directionY * directionY;
              if (!(transverse2 < 1.0 - 1.0e-14)) continue;
              const double directionZ = std::sqrt(1.0 - transverse2);
              const double horizontalNorm = std::sqrt(std::max(0.0,
                1.0 - directionX * directionX));
              if (!(horizontalNorm > 0.0)) continue;
              double* horizontal = &geometry.horizontal[3 * point];
              double* vertical = &geometry.vertical[3 * point];
              horizontal[0] = (1.0 - directionX * directionX) /
                horizontalNorm;
              horizontal[1] = -directionX * directionY / horizontalNorm;
              horizontal[2] = -directionX * directionZ / horizontalNorm;
              vertical[0] = directionY * horizontal[2] -
                directionZ * horizontal[1];
              vertical[1] = directionZ * horizontal[0] -
                directionX * horizontal[2];
              vertical[2] = directionX * horizontal[1] -
                directionY * horizontal[0];
              geometry.thetaX[point] = std::atan2(directionX, directionZ);
              geometry.thetaY[point] = std::atan2(directionY, directionZ);
              geometry.factor[point] = oneSided * kEpsilon0 * kC *
                waveNumber * waveNumber * directionZ * directionZ /
                (4.0 * kPi * kPi * kPlanckEVSecond);
              geometry.solidAngle[point] = deltaKx * deltaKy /
                (waveNumber * waveNumber * directionZ);
            }
        }
      return geometry;
    }

    std::size_t nearestReference(const ProjectionGeometry& geometry,
                                 const ReferenceAngle& requested)
    {
      const double tangentX = std::tan(requested.thetaX);
      const double tangentY = std::tan(requested.thetaY);
      const double norm = std::sqrt(1.0 + tangentX * tangentX +
        tangentY * tangentY);
      const double requestedX = tangentX / norm;
      const double requestedY = tangentY / norm;
      std::size_t nearest = geometry.factor.size();
      double minimum = std::numeric_limits<double>::infinity();
      for (std::size_t point = 0; point < geometry.factor.size(); ++point)
        {
          if (!(geometry.factor[point] > 0.0)) continue;
          const double actualTangentX = std::tan(geometry.thetaX[point]);
          const double actualTangentY = std::tan(geometry.thetaY[point]);
          const double actualNorm = std::sqrt(1.0 +
            actualTangentX * actualTangentX +
            actualTangentY * actualTangentY);
          const double differenceX = actualTangentX / actualNorm - requestedX;
          const double differenceY = actualTangentY / actualNorm - requestedY;
          const double distance = differenceX * differenceX +
            differenceY * differenceY;
          if (distance < minimum)
            {
              minimum = distance;
              nearest = point;
            }
        }
      if (nearest == geometry.factor.size())
        throw std::runtime_error(
          "Angular grid contains no forward propagating FFT bin");
      return nearest;
    }

    std::complex<double> fftValue(const fftw_complex& value)
    {
      return std::complex<double>(value[0], value[1]);
    }

    ComplexValue storedComplex(const std::complex<double>& value)
    {
      ComplexValue result;
      result.real = value.real();
      result.imaginary = value.imag();
      return result;
    }

    double transverseWeight(std::size_t index, std::size_t count,
                            const std::string& window)
    {
      if (window == "none" || count == 1) return 1.0;
      return 0.5 * (1.0 - std::cos(2.0 * kPi *
        static_cast<double>(index) / static_cast<double>(count - 1)));
    }

    std::size_t temporalAmplitudeIndex(std::size_t sample,
                                       std::size_t reference,
                                       std::size_t frequency,
                                       std::size_t polarization,
                                       std::size_t references,
                                       std::size_t frequencies)
    {
      return (((sample * references + reference) * frequencies +
        frequency) * 2 + polarization);
    }

    std::size_t temporalCsdIndex(std::size_t reference,
                                 std::size_t firstFrequency,
                                 std::size_t secondFrequency,
                                 std::size_t firstPolarization,
                                 std::size_t secondPolarization,
                                 std::size_t frequencies)
    {
      return (((((reference * frequencies + firstFrequency) * frequencies +
        secondFrequency) * 2 + firstPolarization) * 2) +
        secondPolarization);
    }
  }

  AnalysisSummary::AnalysisSummary()
    : inputSamples(0), ensembleSamples(0), frequencyBins(0),
      uniformTimeStep(0.0), nyquistPhotonEnergyEV(0.0),
      meanBandEnergy(0.0), bandCoherentFraction(0.0),
      estimatedWorkingMiB(0.0), estimatedOutputGiB(0.0)
  {}

  AnalysisSummary analyze(const Configuration& config)
  {
    FieldInput field(config.fieldFile, config.requireComplete);
    FieldInput* baselinePointer = NULL;
    if (!config.baselineFile.empty())
      {
        baselinePointer = new FieldInput(config.baselineFile,
          config.requireComplete);
        validateMatched(field, *baselinePointer);
      }
    struct BaselineGuard
    {
      FieldInput*& pointer;
      explicit BaselineGuard(FieldInput*& value) : pointer(value) {}
      ~BaselineGuard() { delete pointer; }
    } baselineGuard(baselinePointer);

    const TimeModel time = makeTimeModel(field, config);
    const FrequencyModel frequencies = makeFrequencyModel(time, config);
    const std::vector<std::size_t> spatialFrequency = mapEnergies(
      config.spatialPhotonEnergyEV, frequencies,
      "Spatial-coherence photon energy");
    const std::vector<std::size_t> temporalFrequency = mapEnergies(
      config.temporalPhotonEnergyEV, frequencies,
      "Temporal-coherence photon energy");
    const std::size_t ny = checkedProduct(field.ny(),
      config.angularPaddingY, "padded angular y size");
    const std::size_t nx = checkedProduct(field.nx(),
      config.angularPaddingX, "padded angular x size");
    if (ny > static_cast<std::size_t>(INT_MAX) ||
        nx > static_cast<std::size_t>(INT_MAX) ||
        time.segmentSamples > static_cast<std::size_t>(INT_MAX))
      throw std::runtime_error("FFTW dimensions exceed the supported int range");
    const std::size_t points = checkedProduct(ny, nx,
      "padded angular plane");
    const std::size_t samples = time.windows.size();
    const std::size_t blockLimit = std::min(config.frequencyBlock,
      frequencies.photonEnergy.size());
    if (config.fftThreads > static_cast<std::size_t>(INT_MAX))
      throw std::runtime_error("fft_threads exceeds the supported int range");

    std::vector<double> kx(nx);
    std::vector<double> ky(ny);
    for (std::size_t x = 0; x < nx; ++x)
      {
        const long long q = static_cast<long long>(x) -
          static_cast<long long>(nx / 2);
        kx[x] = 2.0 * kPi * static_cast<double>(q) /
          (static_cast<double>(nx) * field.dx());
      }
    for (std::size_t y = 0; y < ny; ++y)
      {
        const long long q = static_cast<long long>(y) -
          static_cast<long long>(ny / 2);
        ky[y] = 2.0 * kPi * static_cast<double>(q) /
          (static_cast<double>(ny) * field.dy());
      }

    std::size_t apertureValues = checkedProduct(samples, blockLimit,
      "aperture samples and frequencies");
    apertureValues = checkedProduct(apertureValues, 3,
      "aperture components");
    apertureValues = checkedProduct(apertureValues, points,
      "aperture field");
    std::size_t workingBytes = checkedProduct(apertureValues,
      sizeof(fftw_complex), "aperture bytes");
    std::size_t polarizedValues = checkedProduct(samples, points,
      "polarized samples");
    polarizedValues = checkedProduct(polarizedValues, 2,
      "polarized components");
    workingBytes = checkedSum(workingBytes, checkedProduct(
      polarizedValues, sizeof(std::complex<double>), "polarized bytes"),
      "working memory");
    std::size_t perPointBytes = 0;
    perPointBytes = checkedSum(perPointBytes, 2 * sizeof(ComplexValue),
      "per-point output memory");
    perPointBytes = checkedSum(perPointBytes, 8 * sizeof(double),
      "per-point output memory");
    perPointBytes = checkedSum(perPointBytes, 10 * sizeof(double) +
      sizeof(std::size_t), "projection geometry memory");
    workingBytes = checkedSum(workingBytes, checkedProduct(points,
      perPointBytes, "per-point work arrays"), "working memory");
    if (!spatialFrequency.empty())
      {
        std::size_t spatialSlice = checkedProduct(
          config.spatialReferenceAngles.size(), points,
          "spatial CSD slice");
        spatialSlice = checkedProduct(spatialSlice,
          4 * sizeof(ComplexValue) + sizeof(double),
          "spatial CSD bytes");
        workingBytes = checkedSum(workingBytes, spatialSlice,
          "working memory");
      }
    std::size_t temporalValues = checkedProduct(samples,
      config.temporalReferenceAngles.size(), "temporal reference samples");
    temporalValues = checkedProduct(temporalValues,
      temporalFrequency.size(), "temporal selected frequencies");
    temporalValues = checkedProduct(temporalValues, 2,
      "temporal polarization");
    workingBytes = checkedSum(workingBytes, checkedProduct(temporalValues,
      sizeof(std::complex<double>), "temporal amplitude bytes"),
      "working memory");
    const std::size_t gramValues = checkedProduct(samples, samples,
      "ensemble Gram matrix");
    workingBytes = checkedSum(workingBytes, checkedProduct(gramValues,
      sizeof(std::complex<double>) + sizeof(ComplexValue),
      "ensemble Gram work bytes"), "working memory");
    const std::size_t readPoints = std::min(config.spatialBatchPoints,
      field.nx());
    std::size_t inputValues = checkedProduct(field.samples(), readPoints,
      "input read block");
    inputValues = checkedProduct(inputValues, 3, "input read block");
    inputValues = checkedProduct(inputValues,
      baselinePointer == NULL ? 1 : 2, "signal/baseline read blocks");
    workingBytes = checkedSum(workingBytes, checkedProduct(inputValues,
      sizeof(double), "input read-block bytes"), "working memory");
    workingBytes = checkedSum(workingBytes, checkedProduct(
      field.samples() + time.segmentSamples +
      (time.segmentSamples / 2 + 1) * 2,
      sizeof(double), "temporal work bytes"), "working memory");
    const double workingMiB = 1.25 * static_cast<double>(workingBytes) /
      (1024.0 * 1024.0);
    if (workingMiB > static_cast<double>(config.maximumWorkingMiB))
      throw std::runtime_error(
        "Estimated working memory exceeds calculation.maximum_working_mib; "
        "reduce frequency_block, angular zero padding, or the number of time windows");

    long double outputBytes = static_cast<long double>(
      frequencies.photonEnergy.size()) * points * 96.0L;
    outputBytes += static_cast<long double>(samples) *
      frequencies.photonEnergy.size() * sizeof(double);
    outputBytes += static_cast<long double>(frequencies.photonEnergy.size()) *
      samples * samples * sizeof(ComplexValue);
    outputBytes += static_cast<long double>(spatialFrequency.size()) *
      config.spatialReferenceAngles.size() * points * 72.0L;
    outputBytes += static_cast<long double>(config.temporalReferenceAngles.size()) *
      temporalFrequency.size() * temporalFrequency.size() *
      (4.0L * sizeof(ComplexValue) + sizeof(double));
    const double outputGiB = static_cast<double>(outputBytes /
      (1024.0L * 1024.0L * 1024.0L));
    if (outputGiB > config.maximumOutputGiB)
      throw std::runtime_error(
        "Estimated HDF5 output exceeds calculation.maximum_output_gib; "
        "reduce angular padding, photon band, or spatial coherence references");

    std::cout << std::setprecision(10)
              << "Field input: samples=" << field.samples()
              << ", grid=" << field.nx() << "x" << field.ny()
              << ", uniform dt=" << time.dt << " s"
              << ", maximum timing irregularity="
              << time.maximumIrregularity << " steps\n"
              << "Ensemble: " << samples << " sample(s), segment="
              << time.segmentSamples << " points, retained FFT bins="
              << frequencies.photonEnergy.size() << ", angular FFT="
              << nx << "x" << ny << "\n"
              << "Resource guards: working estimate=" << workingMiB
              << " MiB, output estimate=" << outputGiB << " GiB\n";

    if (fftw_init_threads() == 0)
      throw std::runtime_error("Cannot initialize FFTW threads");
    fftw_plan_with_nthreads(static_cast<int>(config.fftThreads));
    double* temporalInput = fftw_alloc_real(time.segmentSamples);
    fftw_complex* temporalOutput = fftw_alloc_complex(
      time.segmentSamples / 2 + 1);
    if (!temporalInput || !temporalOutput) throw std::bad_alloc();
    fftw_plan temporalPlan = fftw_plan_dft_r2c_1d(
      static_cast<int>(time.segmentSamples), temporalInput, temporalOutput,
      FFTW_ESTIMATE);
    if (!temporalPlan)
      throw std::runtime_error("Cannot create temporal FFT plan");

    OutputFile output(config, field, time, frequencies, ny, nx, ky, kx,
      spatialFrequency, temporalFrequency, workingMiB, outputGiB);
    std::vector<double> sampleBandEnergy(samples, 0.0);
    double coherentBandEnergy = 0.0;
    std::vector<std::complex<double> > temporalAmplitude(temporalValues,
      std::complex<double>(0.0, 0.0));
    std::vector<double> temporalActualAngles(checkedProduct(
      checkedProduct(temporalFrequency.size(),
        config.temporalReferenceAngles.size(),
        "temporal actual angles"), 2, "temporal actual angles"), 0.0);

    const std::size_t offsetY = (ny - field.ny()) / 2;
    const std::size_t offsetX = (nx - field.nx()) / 2;
    const double paddedXFirst = field.xFirst() -
      static_cast<double>(offsetX) * field.dx();
    const double paddedYFirst = field.yFirst() -
      static_cast<double>(offsetY) * field.dy();
    std::vector<double> signalBlock;
    std::vector<double> baselineBlock;
    std::vector<double> uniformSeries(field.samples());

    for (std::size_t frequencyOffset = 0;
         frequencyOffset < frequencies.photonEnergy.size();
         frequencyOffset += config.frequencyBlock)
      {
        const std::size_t blockCount = std::min(config.frequencyBlock,
          frequencies.photonEnergy.size() - frequencyOffset);
        std::size_t blockValues = checkedProduct(samples, blockCount,
          "frequency-block samples");
        blockValues = checkedProduct(blockValues, 3,
          "frequency-block components");
        blockValues = checkedProduct(blockValues, points,
          "frequency-block aperture");
        fftw_complex* aperture = fftw_alloc_complex(blockValues);
        if (!aperture) throw std::bad_alloc();
        std::memset(aperture, 0, blockValues * sizeof(fftw_complex));

        for (std::size_t y = 0; y < field.ny(); ++y)
          for (std::size_t xStart = 0; xStart < field.nx();
               xStart += config.spatialBatchPoints)
            {
              const std::size_t count = std::min(config.spatialBatchPoints,
                field.nx() - xStart);
              field.readElectricBlock(y, xStart, count, signalBlock);
              if (baselinePointer != NULL)
                baselinePointer->readElectricBlock(y, xStart, count,
                  baselineBlock);
              for (std::size_t localPoint = 0; localPoint < count;
                   ++localPoint)
                {
                  const std::size_t originalX = xStart + localPoint;
                  const std::size_t paddedPoint = (y + offsetY) * nx +
                    originalX + offsetX;
                  const double spatialWeight = transverseWeight(y,
                    field.ny(), config.transverseWindow) *
                    transverseWeight(originalX, field.nx(),
                      config.transverseWindow);
                  for (unsigned int component = 0; component < 3; ++component)
                    {
                      for (std::size_t target = 0;
                           target < field.samples(); ++target)
                        {
                          const std::size_t lower =
                            time.interpolationLower[target];
                          const double fraction =
                            time.interpolationFraction[target];
                          const std::size_t lowerIndex =
                            (lower * count + localPoint) * 3 + component;
                          const std::size_t upperIndex =
                            ((lower + 1) * count + localPoint) * 3 + component;
                          const double signal = signalBlock[lowerIndex] +
                            fraction * (signalBlock[upperIndex] -
                              signalBlock[lowerIndex]);
                          double baseline = 0.0;
                          if (baselinePointer != NULL)
                            baseline = baselineBlock[lowerIndex] + fraction *
                              (baselineBlock[upperIndex] -
                               baselineBlock[lowerIndex]);
                          uniformSeries[target] = signal - baseline;
                        }
                      for (std::size_t sample = 0; sample < samples; ++sample)
                        {
                          const TimeWindow& window = time.windows[sample];
                          for (std::size_t point = 0;
                               point < time.segmentSamples; ++point)
                            temporalInput[point] =
                              uniformSeries[window.start + point] *
                              time.windowWeight[point];
                          fftw_execute(temporalPlan);
                          const double segmentStartTime = time.start +
                            static_cast<double>(window.start) * time.dt;
                          for (std::size_t localFrequency = 0;
                               localFrequency < blockCount; ++localFrequency)
                            {
                              const std::size_t globalFrequency =
                                frequencyOffset + localFrequency;
                              const std::size_t bin =
                                frequencies.bins[globalFrequency];
                              const double phaseAngle = -2.0 * kPi *
                                frequencies.frequency[globalFrequency] *
                                segmentStartTime;
                              const std::complex<double> phase =
                                std::polar(1.0, phaseAngle);
                              const std::complex<double> value =
                                fftValue(temporalOutput[bin]) * time.dt *
                                spatialWeight * phase;
                              const std::size_t index =
                                (((sample * blockCount + localFrequency) * 3 +
                                   component) * points + paddedPoint);
                              aperture[index][0] = value.real();
                              aperture[index][1] = value.imag();
                            }
                        }
                    }
                }
            }

        fftw_plan spatialPlan = fftw_plan_dft_2d(static_cast<int>(ny),
          static_cast<int>(nx), aperture, aperture, FFTW_FORWARD,
          FFTW_ESTIMATE);
        if (!spatialPlan)
          {
            fftw_free(aperture);
            throw std::runtime_error("Cannot create spatial FFT plan");
          }
        for (std::size_t sample = 0; sample < samples; ++sample)
          for (std::size_t localFrequency = 0;
               localFrequency < blockCount; ++localFrequency)
            for (unsigned int component = 0; component < 3; ++component)
              {
                fftw_complex* plane = aperture +
                  (((sample * blockCount + localFrequency) * 3 + component) *
                   points);
                fftw_execute_dft(spatialPlan, plane, plane);
              }

        for (std::size_t localFrequency = 0;
             localFrequency < blockCount; ++localFrequency)
          {
            const std::size_t globalFrequency = frequencyOffset +
              localFrequency;
            const double frequency = frequencies.frequency[globalFrequency];
            const ProjectionGeometry geometry = makeProjectionGeometry(
              frequency, frequencies.bins[globalFrequency],
              time.segmentSamples, ny, nx, ky, kx);
            std::vector<std::complex<double> > polarized(polarizedValues,
              std::complex<double>(0.0, 0.0));
            for (std::size_t point = 0; point < points; ++point)
              {
                if (!(geometry.factor[point] > 0.0)) continue;
                const std::size_t source = geometry.sourceIndex[point];
                const std::complex<double> originPhase = std::polar(
                  field.dx() * field.dy(),
                  -(kx[point % nx] * paddedXFirst +
                    ky[point / nx] * paddedYFirst));
                const double* horizontal = &geometry.horizontal[3 * point];
                const double* vertical = &geometry.vertical[3 * point];
                for (std::size_t sample = 0; sample < samples; ++sample)
                  {
                    std::complex<double> electric[3];
                    for (unsigned int component = 0; component < 3; ++component)
                      {
                        const std::size_t index =
                          (((sample * blockCount + localFrequency) * 3 +
                             component) * points + source);
                        electric[component] = fftValue(aperture[index]) *
                          originPhase;
                      }
                    std::complex<double> h(0.0, 0.0);
                    std::complex<double> v(0.0, 0.0);
                    for (unsigned int component = 0; component < 3; ++component)
                      {
                        h += horizontal[component] * electric[component];
                        v += vertical[component] * electric[component];
                      }
                    polarized[(sample * points + point) * 2] = h;
                    polarized[(sample * points + point) * 2 + 1] = v;
                  }
              }

            std::vector<ComplexValue> meanField(points * 2);
            std::vector<double> meanIntensity(points, 0.0);
            std::vector<double> coherent(points, 0.0);
            std::vector<double> fluctuation(points, 0.0);
            std::vector<double> stokes(points * 4, 0.0);
            std::vector<double> sampleSpectrum(samples, 0.0);
            double meanSpectrum = 0.0;
            double coherentSpectrum = 0.0;
            for (std::size_t point = 0; point < points; ++point)
              {
                std::complex<double> meanH(0.0, 0.0);
                std::complex<double> meanV(0.0, 0.0);
                double meanI = 0.0;
                double meanQ = 0.0;
                double meanU = 0.0;
                double meanVStokes = 0.0;
                for (std::size_t sample = 0; sample < samples; ++sample)
                  {
                    const std::complex<double>& h =
                      polarized[(sample * points + point) * 2];
                    const std::complex<double>& v =
                      polarized[(sample * points + point) * 2 + 1];
                    const double h2 = std::norm(h);
                    const double v2 = std::norm(v);
                    const std::complex<double> cross = h * std::conj(v);
                    meanH += h;
                    meanV += v;
                    meanI += h2 + v2;
                    meanQ += h2 - v2;
                    meanU += 2.0 * cross.real();
                    meanVStokes += -2.0 * cross.imag();
                    sampleSpectrum[sample] += geometry.factor[point] *
                      (h2 + v2) * geometry.solidAngle[point];
                  }
                const double inverseSamples = 1.0 /
                  static_cast<double>(samples);
                meanH *= inverseSamples;
                meanV *= inverseSamples;
                meanI *= inverseSamples;
                meanQ *= inverseSamples;
                meanU *= inverseSamples;
                meanVStokes *= inverseSamples;
                meanField[2 * point] = storedComplex(meanH);
                meanField[2 * point + 1] = storedComplex(meanV);
                meanIntensity[point] = geometry.factor[point] * meanI;
                coherent[point] = geometry.factor[point] *
                  (std::norm(meanH) + std::norm(meanV));
                fluctuation[point] = std::max(0.0,
                  meanIntensity[point] - coherent[point]);
                stokes[4 * point] = meanIntensity[point];
                stokes[4 * point + 1] = geometry.factor[point] * meanQ;
                stokes[4 * point + 2] = geometry.factor[point] * meanU;
                stokes[4 * point + 3] =
                  geometry.factor[point] * meanVStokes;
                meanSpectrum += meanIntensity[point] *
                  geometry.solidAngle[point];
                coherentSpectrum += coherent[point] *
                  geometry.solidAngle[point];
              }
            const double coherentFraction = meanSpectrum > 0.0 ?
              coherentSpectrum / meanSpectrum : 0.0;
            std::vector<ComplexValue> ensembleGram(gramValues);
            long double gramNorm2 = 0.0L;
            for (std::size_t firstSample = 0; firstSample < samples;
                 ++firstSample)
              for (std::size_t secondSample = 0; secondSample < samples;
                   ++secondSample)
                {
                  std::complex<double> overlap(0.0, 0.0);
                  for (std::size_t point = 0; point < points; ++point)
                    for (unsigned int polarization = 0; polarization < 2;
                         ++polarization)
                      overlap += std::conj(polarized[
                        (firstSample * points + point) * 2 + polarization]) *
                        polarized[(secondSample * points + point) * 2 +
                          polarization] * geometry.factor[point] *
                        geometry.solidAngle[point];
                  overlap /= static_cast<double>(samples);
                  ensembleGram[firstSample * samples + secondSample] =
                    storedComplex(overlap);
                  gramNorm2 += static_cast<long double>(std::norm(overlap));
                }
            const double globalCoherence = meanSpectrum > 0.0 ?
              std::min(1.0, static_cast<double>(gramNorm2) /
                (meanSpectrum * meanSpectrum)) : 0.0;
            output.writeFrequency(globalFrequency, meanField, meanIntensity,
              coherent, fluctuation, stokes, geometry.solidAngle,
              sampleSpectrum, ensembleGram, meanSpectrum, coherentFraction,
              globalCoherence);
            for (std::size_t sample = 0; sample < samples; ++sample)
              sampleBandEnergy[sample] += sampleSpectrum[sample] *
                frequencies.binWidthEV;
            coherentBandEnergy += coherentSpectrum * frequencies.binWidthEV;

            for (std::size_t selected = 0;
                 selected < spatialFrequency.size(); ++selected)
              if (spatialFrequency[selected] == globalFrequency)
                {
                  const std::size_t references =
                    config.spatialReferenceAngles.size();
                  std::vector<double> actualAngles(references * 2);
                  std::vector<std::uint64_t> actualIndices(references * 2);
                  std::vector<ComplexValue> csd(checkedProduct(
                    checkedProduct(references, points, "spatial CSD"), 4,
                    "spatial CSD"));
                  std::vector<double> mu2(checkedProduct(references, points,
                    "spatial coherence"), 0.0);
                  for (std::size_t reference = 0; reference < references;
                       ++reference)
                    {
                      const std::size_t referencePoint = nearestReference(
                        geometry, config.spatialReferenceAngles[reference]);
                      actualAngles[2 * reference] =
                        geometry.thetaX[referencePoint];
                      actualAngles[2 * reference + 1] =
                        geometry.thetaY[referencePoint];
                      actualIndices[2 * reference] = referencePoint / nx;
                      actualIndices[2 * reference + 1] = referencePoint % nx;
                      double referenceIntensity = 0.0;
                      for (std::size_t sample = 0; sample < samples; ++sample)
                        referenceIntensity += std::norm(
                          polarized[(sample * points + referencePoint) * 2]) +
                          std::norm(polarized[
                            (sample * points + referencePoint) * 2 + 1]);
                      referenceIntensity *= geometry.factor[referencePoint] /
                        static_cast<double>(samples);
                      for (std::size_t point = 0; point < points; ++point)
                        {
                          double numerator = 0.0;
                          for (unsigned int a = 0; a < 2; ++a)
                            for (unsigned int b = 0; b < 2; ++b)
                              {
                                std::complex<double> value(0.0, 0.0);
                                for (std::size_t sample = 0; sample < samples;
                                     ++sample)
                                  value += std::conj(polarized[
                                    (sample * points + referencePoint) * 2 + a]) *
                                    polarized[(sample * points + point) * 2 + b];
                                value *= std::sqrt(
                                  geometry.factor[referencePoint] *
                                  geometry.factor[point]) /
                                  static_cast<double>(samples);
                                csd[((reference * points + point) * 2 + a) *
                                  2 + b] = storedComplex(value);
                                numerator += std::norm(value);
                              }
                          const double denominator = referenceIntensity *
                            meanIntensity[point];
                          mu2[reference * points + point] = denominator > 0.0 ?
                            std::min(1.0, numerator / denominator) : 0.0;
                        }
                    }
                  output.writeSpatial(selected, actualAngles, actualIndices,
                    csd, mu2, references);
                }

            for (std::size_t selected = 0;
                 selected < temporalFrequency.size(); ++selected)
              if (temporalFrequency[selected] == globalFrequency)
                for (std::size_t reference = 0;
                     reference < config.temporalReferenceAngles.size();
                     ++reference)
                  {
                    const std::size_t referencePoint = nearestReference(
                      geometry, config.temporalReferenceAngles[reference]);
                    temporalActualAngles[
                      (selected * config.temporalReferenceAngles.size() +
                       reference) * 2] = geometry.thetaX[referencePoint];
                    temporalActualAngles[
                      (selected * config.temporalReferenceAngles.size() +
                       reference) * 2 + 1] = geometry.thetaY[referencePoint];
                    const double scale = std::sqrt(
                      geometry.factor[referencePoint]);
                    for (std::size_t sample = 0; sample < samples; ++sample)
                      for (unsigned int polarization = 0; polarization < 2;
                           ++polarization)
                        temporalAmplitude[temporalAmplitudeIndex(sample,
                          reference, selected, polarization,
                          config.temporalReferenceAngles.size(),
                          temporalFrequency.size())] = scale * polarized[
                            (sample * points + referencePoint) * 2 +
                            polarization];
                  }
          }
        fftw_destroy_plan(spatialPlan);
        fftw_free(aperture);
        std::cout << "Processed field FFT bins " << frequencyOffset + 1
                  << "--" << frequencyOffset + blockCount << " of "
                  << frequencies.photonEnergy.size() << "\n";
      }

    if (!temporalFrequency.empty())
      {
        const std::size_t references =
          config.temporalReferenceAngles.size();
        const std::size_t selectedCount = temporalFrequency.size();
        const std::size_t csdValues = checkedProduct(checkedProduct(
          checkedProduct(references, selectedCount, "temporal CSD"),
          selectedCount, "temporal CSD"), 4, "temporal CSD");
        std::vector<ComplexValue> csd(csdValues);
        std::vector<double> mu2(checkedProduct(checkedProduct(references,
          selectedCount, "temporal coherence"), selectedCount,
          "temporal coherence"), 0.0);
        for (std::size_t reference = 0; reference < references; ++reference)
          for (std::size_t first = 0; first < selectedCount; ++first)
            for (std::size_t second = 0; second < selectedCount; ++second)
              {
                double firstIntensity = 0.0;
                double secondIntensity = 0.0;
                for (std::size_t sample = 0; sample < samples; ++sample)
                  for (unsigned int polarization = 0; polarization < 2;
                       ++polarization)
                    {
                      firstIntensity += std::norm(temporalAmplitude[
                        temporalAmplitudeIndex(sample, reference, first,
                          polarization, references, selectedCount)]);
                      secondIntensity += std::norm(temporalAmplitude[
                        temporalAmplitudeIndex(sample, reference, second,
                          polarization, references, selectedCount)]);
                    }
                firstIntensity /= static_cast<double>(samples);
                secondIntensity /= static_cast<double>(samples);
                double numerator = 0.0;
                for (unsigned int a = 0; a < 2; ++a)
                  for (unsigned int b = 0; b < 2; ++b)
                    {
                      std::complex<double> value(0.0, 0.0);
                      for (std::size_t sample = 0; sample < samples; ++sample)
                        value += std::conj(temporalAmplitude[
                          temporalAmplitudeIndex(sample, reference, first, a,
                            references, selectedCount)]) *
                          temporalAmplitude[temporalAmplitudeIndex(sample,
                            reference, second, b, references, selectedCount)];
                      value /= static_cast<double>(samples);
                      csd[temporalCsdIndex(reference, first, second, a, b,
                        selectedCount)] = storedComplex(value);
                      numerator += std::norm(value);
                    }
                const double denominator = firstIntensity * secondIntensity;
                mu2[(reference * selectedCount + first) * selectedCount +
                  second] = denominator > 0.0 ?
                  std::min(1.0, numerator / denominator) : 0.0;
              }
        output.writeTemporal(temporalActualAngles, csd, mu2, references,
          selectedCount);
      }

    const double meanBandEnergy = std::accumulate(sampleBandEnergy.begin(),
      sampleBandEnergy.end(), 0.0) / static_cast<double>(samples);
    const double bandCoherentFraction = meanBandEnergy > 0.0 ?
      coherentBandEnergy / meanBandEnergy : 0.0;
    output.finish(sampleBandEnergy, meanBandEnergy, bandCoherentFraction);
    fftw_destroy_plan(temporalPlan);
    fftw_free(temporalOutput);
    fftw_free(temporalInput);
    fftw_cleanup_threads();

    AnalysisSummary summary;
    summary.inputSamples = field.samples();
    summary.ensembleSamples = samples;
    summary.frequencyBins = frequencies.photonEnergy.size();
    summary.uniformTimeStep = time.dt;
    summary.nyquistPhotonEnergyEV = frequencies.nyquistEV;
    summary.meanBandEnergy = meanBandEnergy;
    summary.bandCoherentFraction = bandCoherentFraction;
    summary.estimatedWorkingMiB = workingMiB;
    summary.estimatedOutputGiB = outputGiB;
    return summary;
  }
}
