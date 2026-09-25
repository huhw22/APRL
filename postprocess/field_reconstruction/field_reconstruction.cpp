#include "field_reconstruction.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "fftw3.h"
#include "hdf5.h"

namespace reconstruction
{
  namespace
  {
    const double kPi = 3.141592653589793238462643383279502884;
    const double kC = 299792458.0;
    const double kEpsilon0 = 8.8541878128e-12;
    const double kMu0 = 1.25663706212e-6;

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
                throw std::runtime_error("Cannot create output directory: " +
                  current + ": " + std::strerror(errno));
            }
          if (separator == std::string::npos) break;
          offset = separator + 1;
        }
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

    void writeStringAttribute(hid_t object, const char* name,
                              const std::string& value)
    {
      hid_t type = H5Tcopy(H5T_C_S1);
      requireHandle(type, "Cannot create string attribute type");
      requireStatus(H5Tset_size(type, value.size() + 1),
        "Cannot size string attribute");
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
      requireHandle(space, "Cannot create integer attribute space");
      hid_t attribute = H5Acreate2(object, name, H5T_STD_U64LE, space,
        H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(space);
      requireHandle(attribute, std::string("Cannot create attribute: ") + name);
      const herr_t status = H5Awrite(attribute, H5T_NATIVE_UINT64, &value);
      H5Aclose(attribute);
      requireStatus(status, std::string("Cannot write attribute: ") + name);
    }

    struct FieldPoint
    {
      double electric[3];
      double magnetic[3];
    };

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
      requireStatus(H5Tinsert(type, "time_s",
        HOFFSET(ParticleRecord, time), H5T_NATIVE_DOUBLE),
        "Cannot map particle time");
      requireStatus(H5Tinsert(type, "position_m",
        HOFFSET(ParticleRecord, position), vector),
        "Cannot map particle position");
      requireStatus(H5Tinsert(type, "proper_velocity",
        HOFFSET(ParticleRecord, properVelocity), vector),
        "Cannot map particle proper velocity");
      requireStatus(H5Tinsert(type, "charge_C",
        HOFFSET(ParticleRecord, charge), H5T_NATIVE_DOUBLE),
        "Cannot map particle charge");
      requireStatus(H5Tinsert(type, "mass_kg",
        HOFFSET(ParticleRecord, mass), H5T_NATIVE_DOUBLE),
        "Cannot map particle mass");
      requireStatus(H5Tinsert(type, "weight",
        HOFFSET(ParticleRecord, weight), H5T_NATIVE_DOUBLE),
        "Cannot map particle weight");
      requireStatus(H5Tinsert(type, "crossing_direction",
        HOFFSET(ParticleRecord, direction), H5T_NATIVE_INT32),
        "Cannot map particle direction");
      H5Tclose(vector);
      return type;
    }

    class FieldInput
    {
    public:
      FieldInput(const std::string& filename, bool requireComplete)
        : filename_(filename), file_(-1), group_(-1), time_(-1), fields_(-1),
          memoryType_(-1), samples_(0), nx_(0), ny_(0), planeZ_(0.0),
          xFirst_(0.0), yFirst_(0.0), dx_(0.0), dy_(0.0), times_()
      {
        file_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        requireHandle(file_, "Cannot open raw field file: " + filename);
        try
          {
            group_ = H5Gopen2(file_, "/field_plane", H5P_DEFAULT);
            requireHandle(group_, "Missing /field_plane in " + filename);
            if (requireComplete && readByteDataset(group_, "complete") == 0)
              throw std::runtime_error(
                "Raw field file is incomplete: " + filename);
            samples_ = static_cast<std::size_t>(
              readUnsignedDataset(group_, "committed_samples"));
            nx_ = static_cast<std::size_t>(
              readUnsignedAttribute(group_, "nx"));
            ny_ = static_cast<std::size_t>(
              readUnsignedAttribute(group_, "ny"));
            planeZ_ = readDoubleAttribute(group_, "plane_z_m");
            xFirst_ = readDoubleAttribute(group_, "x_first_m");
            yFirst_ = readDoubleAttribute(group_, "y_first_m");
            dx_ = readDoubleAttribute(group_, "dx_m");
            dy_ = readDoubleAttribute(group_, "dy_m");
            if (samples_ < 2 || nx_ == 0 || ny_ == 0 ||
                !(dx_ > 0.0) || !(dy_ > 0.0))
              throw std::runtime_error("Invalid raw field-plane geometry");
            time_ = H5Dopen2(group_, "time_s", H5P_DEFAULT);
            fields_ = H5Dopen2(group_, "fields", H5P_DEFAULT);
            requireHandle(time_, "Missing field time dataset");
            requireHandle(fields_, "Missing field dataset");
            memoryType_ = createFieldMemoryType();
            times_.resize(samples_);
            hsize_t start[1] = {0};
            hsize_t count[1] = {static_cast<hsize_t>(samples_)};
            hid_t fileSpace = H5Dget_space(time_);
            requireHandle(fileSpace, "Cannot inspect field time dataset");
            requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
              start, NULL, count, NULL), "Cannot select field times");
            hid_t memorySpace = H5Screate_simple(1, count, NULL);
            requireHandle(memorySpace, "Cannot create field time memory space");
            const herr_t status = H5Dread(time_, H5T_NATIVE_DOUBLE,
              memorySpace, fileSpace, H5P_DEFAULT, &times_[0]);
            H5Sclose(memorySpace);
            H5Sclose(fileSpace);
            requireStatus(status, "Cannot read field sample times");
          }
        catch (...)
          {
            close();
            throw;
          }
      }

      ~FieldInput() { close(); }

      void close()
      {
        if (memoryType_ >= 0) H5Tclose(memoryType_);
        if (fields_ >= 0) H5Dclose(fields_);
        if (time_ >= 0) H5Dclose(time_);
        if (group_ >= 0) H5Gclose(group_);
        if (file_ >= 0) H5Fclose(file_);
        memoryType_ = fields_ = time_ = group_ = file_ = -1;
      }

      void readSample(std::size_t sample,
                      std::vector<FieldPoint>& values) const
      {
        if (sample >= samples_) throw std::out_of_range("Field sample index");
        values.resize(nx_ * ny_);
        hsize_t start[3] = {static_cast<hsize_t>(sample), 0, 0};
        hsize_t count[3] = {1, static_cast<hsize_t>(ny_),
                            static_cast<hsize_t>(nx_)};
        hid_t fileSpace = H5Dget_space(fields_);
        requireHandle(fileSpace, "Cannot inspect raw field dataset");
        requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
          start, NULL, count, NULL), "Cannot select raw field sample");
        hsize_t memoryDimensions[2] = {count[1], count[2]};
        hid_t memorySpace = H5Screate_simple(2, memoryDimensions, NULL);
        requireHandle(memorySpace, "Cannot create raw field memory space");
        const herr_t status = H5Dread(fields_, memoryType_, memorySpace,
          fileSpace, H5P_DEFAULT, &values[0]);
        H5Sclose(memorySpace);
        H5Sclose(fileSpace);
        requireStatus(status, "Cannot read raw field sample");
      }

      void readInterpolated(double time,
                            std::vector<FieldPoint>& values) const
      {
        if (time < times_.front() || time > times_.back())
          throw std::out_of_range("Interpolated field time is outside input");
        const std::vector<double>::const_iterator upper =
          std::lower_bound(times_.begin(), times_.end(), time);
        if (upper == times_.begin())
          {
            readSample(0, values);
            return;
          }
        if (upper == times_.end())
          {
            readSample(samples_ - 1, values);
            return;
          }
        const std::size_t second = static_cast<std::size_t>(
          upper - times_.begin());
        const double scale = std::max(std::numeric_limits<double>::min(),
          std::max(std::abs(time), std::abs(times_[second])));
        if (std::abs(time - times_[second]) <=
            64.0 * std::numeric_limits<double>::epsilon() * scale)
          {
            readSample(second, values);
            return;
          }
        const std::size_t first = second - 1;
        std::vector<FieldPoint> left;
        std::vector<FieldPoint> right;
        readSample(first, left);
        readSample(second, right);
        const double fraction = (time - times_[first]) /
          (times_[second] - times_[first]);
        values.resize(left.size());
        for (std::size_t point = 0; point < values.size(); ++point)
          for (unsigned int component = 0; component < 3; ++component)
            {
              values[point].electric[component] =
                left[point].electric[component] + fraction *
                (right[point].electric[component] -
                 left[point].electric[component]);
              values[point].magnetic[component] =
                left[point].magnetic[component] + fraction *
                (right[point].magnetic[component] -
                 left[point].magnetic[component]);
            }
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
      FieldInput(const FieldInput&);
      FieldInput& operator=(const FieldInput&);
      std::string filename_;
      hid_t file_;
      hid_t group_;
      hid_t time_;
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

    class ParticleInput
    {
    public:
      ParticleInput(const std::string& filename, bool requireComplete,
                    double fieldPlaneZ)
        : filename_(filename), file_(-1), group_(-1), records_(-1),
          memoryType_(-1), committed_(0), planeZ_(0.0),
          fieldPlaneZ_(fieldPlaneZ), ballisticReference_(false)
      {
        file_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        requireHandle(file_, "Cannot open particle-plane file: " + filename);
        try
          {
            if (H5Lexists(file_, "/ballistic_reference", H5P_DEFAULT) > 0)
              {
                group_ = H5Gopen2(file_, "/ballistic_reference", H5P_DEFAULT);
                ballisticReference_ = true;
              }
            else if (H5Lexists(file_, "/particle_plane", H5P_DEFAULT) > 0)
              group_ = H5Gopen2(file_, "/particle_plane", H5P_DEFAULT);
            else
              throw std::runtime_error(
                "Particle file has neither /ballistic_reference nor /particle_plane");
            requireHandle(group_, "Cannot open particle crossing group");
            if (requireComplete && readByteDataset(group_, "complete") == 0)
              throw std::runtime_error(
                "Particle crossing file is incomplete: " + filename);
            committed_ = readUnsignedDataset(group_, "committed_records");
            planeZ_ = readDoubleAttribute(group_, "plane_z_m");
            if (ballisticReference_)
              {
                const double declaredDetector = readDoubleAttribute(
                  group_, "field_detector_z_m");
                const double tolerance = 128.0 *
                  std::numeric_limits<double>::epsilon() *
                  std::max(1.0, std::abs(fieldPlaneZ_));
                if (std::abs(declaredDetector - fieldPlaneZ_) > tolerance)
                  throw std::runtime_error(
                    "Ballistic reference belongs to a different field plane");
              }
            else
              {
                const double tolerance = 128.0 *
                  std::numeric_limits<double>::epsilon() *
                  std::max(1.0, std::abs(fieldPlaneZ_));
                if (std::abs(planeZ_ - fieldPlaneZ_) > tolerance)
                  throw std::runtime_error(
                    "Particle detector must be colocated with the field plane");
              }
            records_ = H5Dopen2(group_, "records", H5P_DEFAULT);
            requireHandle(records_, "Missing particle crossing records");
            memoryType_ = createParticleMemoryType();
          }
        catch (...)
          {
            close();
            throw;
          }
      }

      ~ParticleInput() { close(); }

      void close()
      {
        if (memoryType_ >= 0) H5Tclose(memoryType_);
        if (records_ >= 0) H5Dclose(records_);
        if (group_ >= 0) H5Gclose(group_);
        if (file_ >= 0) H5Fclose(file_);
        memoryType_ = records_ = group_ = file_ = -1;
      }

      void read(std::uint64_t offset, std::size_t count,
                std::vector<ParticleRecord>& values) const
      {
        if (offset + count > committed_)
          throw std::out_of_range("Particle record range");
        values.resize(count);
        if (count == 0) return;
        hsize_t start[1] = {static_cast<hsize_t>(offset)};
        hsize_t size[1] = {static_cast<hsize_t>(count)};
        hid_t fileSpace = H5Dget_space(records_);
        requireHandle(fileSpace, "Cannot inspect particle records");
        requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
          start, NULL, size, NULL), "Cannot select particle records");
        hid_t memorySpace = H5Screate_simple(1, size, NULL);
        requireHandle(memorySpace, "Cannot create particle memory space");
        const herr_t status = H5Dread(records_, memoryType_, memorySpace,
          fileSpace, H5P_DEFAULT, &values[0]);
        H5Sclose(memorySpace);
        H5Sclose(fileSpace);
        requireStatus(status, "Cannot read particle crossing records");
      }

      std::uint64_t committed() const { return committed_; }
      double planeZ() const { return planeZ_; }
      bool ballisticReference() const { return ballisticReference_; }
      const std::string& filename() const { return filename_; }

    private:
      ParticleInput(const ParticleInput&);
      ParticleInput& operator=(const ParticleInput&);
      std::string filename_;
      hid_t file_;
      hid_t group_;
      hid_t records_;
      hid_t memoryType_;
      std::uint64_t committed_;
      double planeZ_;
      double fieldPlaneZ_;
      bool ballisticReference_;
    };

    struct ProjectedParticle
    {
      std::uint64_t id;
      double time;
      double x;
      double y;
      double velocity[3];
      double gamma;
      double charge;
    };

    bool projectParticle(const ParticleRecord& record,
                         double sourcePlaneZ,
                         double fieldPlaneZ,
                         ProjectedParticle& result)
    {
      if (record.direction != 1 || record.particleId == 0 ||
          !std::isfinite(record.time) || !std::isfinite(record.charge))
        return false;
      double u2 = 0.0;
      for (unsigned int component = 0; component < 3; ++component)
        {
          if (!std::isfinite(record.position[component]) ||
              !std::isfinite(record.properVelocity[component]))
            return false;
          u2 += record.properVelocity[component] *
            record.properVelocity[component];
        }
      const double gamma = std::sqrt(1.0 + u2);
      if (!(gamma >= 1.0) || !std::isfinite(gamma)) return false;
      result.id = record.particleId;
      result.gamma = gamma;
      result.charge = record.charge;
      for (unsigned int component = 0; component < 3; ++component)
        result.velocity[component] =
          kC * record.properVelocity[component] / gamma;
      if (!(result.velocity[2] > 0.0) ||
          !std::isfinite(result.velocity[2])) return false;
      const double flightTime =
        (fieldPlaneZ - sourcePlaneZ) / result.velocity[2];
      if (!(flightTime >= 0.0) || !std::isfinite(flightTime)) return false;
      result.time = record.time + flightTime;
      result.x = record.position[0] + result.velocity[0] * flightTime;
      result.y = record.position[1] + result.velocity[1] * flightTime;
      return std::isfinite(result.time) && std::isfinite(result.x) &&
             std::isfinite(result.y);
    }

    struct ParticleStatistics
    {
      unsigned long long accepted;
      unsigned long long duplicates;
      unsigned long long invalid;
      long double absoluteCharge;
      long double gamma;
      long double gammaSquared;
      long double betaTransverseSquared;
      long double denominatorCoefficient;

      ParticleStatistics()
        : accepted(0), duplicates(0), invalid(0), absoluteCharge(0.0L),
          gamma(0.0L), gammaSquared(0.0L),
          betaTransverseSquared(0.0L), denominatorCoefficient(0.0L)
      {}
    };

    class DenseParticleIdSet
    {
    public:
      bool insert(std::uint64_t id)
      {
        const std::size_t maximumWords = 67108864;
        const std::uint64_t word64 = id / 64;
        if (word64 >= static_cast<std::uint64_t>(maximumWords))
          throw std::runtime_error(
            "Particle ID exceeds the dense-ID reconstruction limit");
        const std::size_t word = static_cast<std::size_t>(word64);
        if (word >= words_.size())
          {
            const std::size_t doubled = words_.empty() ? 1024 :
              checkedProduct(words_.size(), 2, "particle-ID bitmap");
            const std::size_t requested = word + 1;
            words_.resize(std::max(requested,
              std::min(doubled, maximumWords)), UINT64_C(0));
          }
        const std::uint64_t mask = UINT64_C(1) << (id % 64);
        const bool inserted = (words_[word] & mask) == 0;
        words_[word] |= mask;
        return inserted;
      }

    private:
      std::vector<std::uint64_t> words_;
    };

    void forEachUniqueParticle(const ParticleInput& input,
                               std::size_t chunkRecords,
                               double fieldPlaneZ,
                               const std::function<void(
                                 const ProjectedParticle&)>& callback,
                               unsigned long long& duplicates,
                               unsigned long long& invalid)
    {
      DenseParticleIdSet seen;
      std::vector<ParticleRecord> buffer;
      for (std::uint64_t offset = 0; offset < input.committed();)
        {
          const std::size_t count = static_cast<std::size_t>(
            std::min<std::uint64_t>(chunkRecords,
              input.committed() - offset));
          input.read(offset, count, buffer);
          for (std::size_t index = 0; index < buffer.size(); ++index)
            {
              if (buffer[index].direction != 1) continue;
              ProjectedParticle particle;
              if (!projectParticle(buffer[index], input.planeZ(),
                                   fieldPlaneZ, particle))
                {
                  ++invalid;
                  continue;
                }
              if (!seen.insert(particle.id))
                {
                  ++duplicates;
                  continue;
                }
              callback(particle);
            }
          offset += count;
        }
    }

    struct PaddedGrid
    {
      std::size_t nt;
      std::size_t ny;
      std::size_t nx;
      std::size_t timeOffset;
      std::size_t yOffset;
      std::size_t xOffset;
      std::size_t realSize;
      std::size_t complexSize;
      double timeFirst;
      double yFirst;
      double xFirst;
      double dt;
      double dy;
      double dx;
      std::vector<double> outputTimes;
    };

    PaddedGrid makePaddedGrid(const FieldInput& field,
                              const ReconstructionConfig& config)
    {
      PaddedGrid grid;
      grid.nt = checkedProduct(field.samples(), config.paddingTime,
        "padded time grid");
      grid.ny = checkedProduct(field.ny(), config.paddingY,
        "padded y grid");
      grid.nx = checkedProduct(field.nx(), config.paddingX,
        "padded x grid");
      if (grid.nt > static_cast<std::size_t>(INT_MAX) ||
          grid.ny > static_cast<std::size_t>(INT_MAX) ||
          grid.nx > static_cast<std::size_t>(INT_MAX))
        throw std::overflow_error("Padded FFT dimension exceeds INT_MAX");
      grid.timeOffset = (grid.nt - field.samples()) / 2;
      grid.yOffset = (grid.ny - field.ny()) / 2;
      grid.xOffset = (grid.nx - field.nx()) / 2;
      grid.realSize = checkedProduct(checkedProduct(grid.nt, grid.ny,
        "padded real grid"), grid.nx, "padded real grid");
      grid.complexSize = checkedProduct(checkedProduct(grid.nt, grid.ny,
        "padded spectral grid"), grid.nx / 2 + 1,
        "padded spectral grid");
      const std::vector<double>& times = field.times();
      grid.dt = (times.back() - times.front()) /
        static_cast<double>(times.size() - 1);
      if (!(grid.dt > 0.0) || !std::isfinite(grid.dt))
        throw std::runtime_error("Field sample times are not increasing");
      for (std::size_t index = 1; index < times.size(); ++index)
        if (!(times[index] > times[index - 1]))
          throw std::runtime_error(
            "Field sample times must be strictly increasing");
      grid.dx = field.dx();
      grid.dy = field.dy();
      grid.timeFirst = times.front() - grid.timeOffset * grid.dt;
      grid.yFirst = field.yFirst() - grid.yOffset * grid.dy;
      grid.xFirst = field.xFirst() - grid.xOffset * grid.dx;
      grid.outputTimes.resize(field.samples());
      for (std::size_t index = 0; index < field.samples(); ++index)
        grid.outputTimes[index] = times.front() + index * grid.dt;
      return grid;
    }

    struct AxisWeights
    {
      std::size_t index[2];
      double weight[2];
    };

    bool interpolationWeights(double coordinate, std::size_t size,
                              AxisWeights& result)
    {
      if (!std::isfinite(coordinate) || coordinate < 0.0 ||
          coordinate > static_cast<double>(size - 1)) return false;
      const std::size_t lower = static_cast<std::size_t>(
        std::floor(coordinate));
      const std::size_t upper = std::min(lower + 1, size - 1);
      const double fraction = coordinate - static_cast<double>(lower);
      result.index[0] = lower;
      result.index[1] = upper;
      result.weight[0] = 1.0 - fraction;
      result.weight[1] = fraction;
      if (upper == lower)
        {
          result.weight[0] = 1.0;
          result.weight[1] = 0.0;
        }
      return true;
    }

    bool depositParticle(const ProjectedParticle& particle,
                         const PaddedGrid& grid,
                         double* charge,
                         double* chargeOverVelocity,
                         double* chargeOverGamma2Velocity2)
    {
      AxisWeights time;
      AxisWeights y;
      AxisWeights x;
      if (!interpolationWeights((particle.time - grid.timeFirst) / grid.dt,
                                grid.nt, time) ||
          !interpolationWeights((particle.y - grid.yFirst) / grid.dy,
                                grid.ny, y) ||
          !interpolationWeights((particle.x - grid.xFirst) / grid.dx,
                                grid.nx, x)) return false;
      const double inverseCell = 1.0 / (grid.dt * grid.dy * grid.dx);
      const double valueQ = particle.charge * inverseCell;
      const double valueQV = valueQ / particle.velocity[2];
      const double valueQGammaV = valueQ /
        (particle.gamma * particle.gamma *
         particle.velocity[2] * particle.velocity[2]);
      for (unsigned int a = 0; a < 2; ++a)
        for (unsigned int b = 0; b < 2; ++b)
          for (unsigned int c = 0; c < 2; ++c)
            {
              const double shape = time.weight[a] * y.weight[b] * x.weight[c];
              if (shape == 0.0) continue;
              const std::size_t index =
                (time.index[a] * grid.ny + y.index[b]) * grid.nx +
                x.index[c];
              charge[index] += valueQ * shape;
              chargeOverVelocity[index] += valueQV * shape;
              chargeOverGamma2Velocity2[index] += valueQGammaV * shape;
            }
      return true;
    }

    double signedWaveNumber(std::size_t index, std::size_t size,
                            double spacing)
    {
      const long long signedIndex = index <= size / 2 ?
        static_cast<long long>(index) :
        static_cast<long long>(index) - static_cast<long long>(size);
      return 2.0 * kPi * static_cast<double>(signedIndex) /
        (static_cast<double>(size) * spacing);
    }

    enum FieldComponent
    {
      ElectricX = 0,
      ElectricY = 1,
      ElectricZ = 2,
      MagneticX = 3,
      MagneticY = 4,
      MagneticZ = 5
    };

    void multiplyByImaginary(fftw_complex& destination,
                             const fftw_complex& source,
                             double factor)
    {
      destination[0] = -factor * source[1];
      destination[1] = factor * source[0];
    }

    void backgroundSpectrum(FieldComponent component,
                            const PaddedGrid& grid,
                            const fftw_complex* charge,
                            const fftw_complex* chargeOverVelocity,
                            const fftw_complex* chargeOverGamma2Velocity2,
                            double denominatorCoefficient,
                            double smoothing,
                            fftw_complex* destination)
    {
      const std::size_t kxCount = grid.nx / 2 + 1;
      for (std::size_t it = 0; it < grid.nt; ++it)
        {
          const double omega = signedWaveNumber(it, grid.nt, grid.dt);
          for (std::size_t iy = 0; iy < grid.ny; ++iy)
            {
              const double ky = signedWaveNumber(iy, grid.ny, grid.dy);
              for (std::size_t ix = 0; ix < kxCount; ++ix)
                {
                  const double kx = 2.0 * kPi * static_cast<double>(ix) /
                    (static_cast<double>(grid.nx) * grid.dx);
                  const std::size_t index =
                    (it * grid.ny + iy) * kxCount + ix;
                  const double denominator = kx * kx + ky * ky +
                    omega * omega * denominatorCoefficient;
                  if (!(denominator > 0.0))
                    {
                      destination[index][0] = 0.0;
                      destination[index][1] = 0.0;
                      continue;
                    }
                  const double filter = smoothing > 0.0 ?
                    std::exp(-0.5 * smoothing * smoothing *
                      (kx * kx + ky * ky)) : 1.0;
                  double factor = filter / (kEpsilon0 * denominator);
                  const fftw_complex* source = charge;
                  switch (component)
                    {
                    case ElectricX:
                      factor *= -kx;
                      source = chargeOverVelocity;
                      break;
                    case ElectricY:
                      factor *= -ky;
                      source = chargeOverVelocity;
                      break;
                    case ElectricZ:
                      factor *= omega;
                      source = chargeOverGamma2Velocity2;
                      break;
                    case MagneticX:
                      factor *= ky / (kC * kC);
                      source = charge;
                      break;
                    case MagneticY:
                      factor *= -kx / (kC * kC);
                      source = charge;
                      break;
                    case MagneticZ:
                      factor = 0.0;
                      source = charge;
                      break;
                    }
                  multiplyByImaginary(destination[index], source[index], factor);
                }
            }
        }
    }

    class OutputFile
    {
    public:
      OutputFile(const ReconstructionConfig& config,
                 const FieldInput& field,
                 const ParticleInput& particles,
                 const PaddedGrid& grid,
                 const ReconstructionSummary& summary,
                 double denominatorCoefficient)
        : file_(-1), group_(-1), time_(-1), electric_(-1), magnetic_(-1),
          complete_(-1), nx_(field.nx()), ny_(field.ny()),
          samples_(field.samples())
      {
        createDirectories(parentDirectory(config.outputFile));
        file_ = H5Fcreate(config.outputFile.c_str(),
          config.overwrite ? H5F_ACC_TRUNC : H5F_ACC_EXCL,
          H5P_DEFAULT, H5P_DEFAULT);
        requireHandle(file_, "Cannot create reconstruction output: " +
          config.outputFile);
        try
          {
            group_ = H5Gcreate2(file_, "/reconstructed_field", H5P_DEFAULT,
              H5P_DEFAULT, H5P_DEFAULT);
            requireHandle(group_, "Cannot create reconstructed field group");
            hsize_t timeDimensions[1] = {static_cast<hsize_t>(samples_)};
            hid_t timeSpace = H5Screate_simple(1, timeDimensions, NULL);
            time_ = H5Dcreate2(group_, "time_s", H5T_IEEE_F64LE,
              timeSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(timeSpace);
            requireHandle(time_, "Cannot create reconstructed time dataset");
            requireStatus(H5Dwrite(time_, H5T_NATIVE_DOUBLE, H5S_ALL,
              H5S_ALL, H5P_DEFAULT, &grid.outputTimes[0]),
              "Cannot write reconstructed field times");

            hsize_t dimensions[4] = {
              static_cast<hsize_t>(samples_), static_cast<hsize_t>(ny_),
              static_cast<hsize_t>(nx_), 3
            };
            hid_t space = H5Screate_simple(4, dimensions, NULL);
            requireHandle(space, "Cannot create reconstructed field space");
            hid_t creation = H5Pcreate(H5P_DATASET_CREATE);
            requireHandle(creation,
              "Cannot create reconstructed field properties");
            hsize_t chunk[4] = {1,
              static_cast<hsize_t>(std::min<std::size_t>(ny_, 64)),
              static_cast<hsize_t>(std::min<std::size_t>(nx_, 64)), 3};
            requireStatus(H5Pset_chunk(creation, 4, chunk),
              "Cannot set reconstructed field chunks");
            if (config.compression > 0)
              requireStatus(H5Pset_deflate(creation, config.compression),
                "Cannot set reconstructed field compression");
            electric_ = H5Dcreate2(group_, "electric_V_per_m",
              H5T_IEEE_F64LE, space, H5P_DEFAULT, creation, H5P_DEFAULT);
            magnetic_ = H5Dcreate2(group_, "magnetic_T",
              H5T_IEEE_F64LE, space, H5P_DEFAULT, creation, H5P_DEFAULT);
            H5Pclose(creation);
            H5Sclose(space);
            requireHandle(electric_, "Cannot create reconstructed E dataset");
            requireHandle(magnetic_, "Cannot create reconstructed B dataset");

            hid_t scalar = H5Screate(H5S_SCALAR);
            complete_ = H5Dcreate2(group_, "complete", H5T_STD_U8LE,
              scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(scalar);
            requireHandle(complete_, "Cannot create reconstruction marker");
            const unsigned char incomplete = 0;
            requireStatus(H5Dwrite(complete_, H5T_NATIVE_UCHAR, H5S_ALL,
              H5S_ALL, H5P_DEFAULT, &incomplete),
              "Cannot initialize reconstruction marker");

            writeUnsignedAttribute(group_, "format_version", 1);
            writeStringAttribute(group_, "raw_field_file", field.filename());
            writeStringAttribute(group_, "particle_crossing_file",
              particles.filename());
            writeStringAttribute(group_, "particle_input_role",
              particles.ballisticReference() ? "ballistic_reference" :
              "particle_plane");
            writeStringAttribute(group_, "model",
              "mean-denominator uniform-longitudinal-velocity Maxwell field from CIC crossing density");
            writeStringAttribute(group_, "raw_magnetic_time_stagger",
              "preserved from input; reconstructed background is evaluated at labelled sample time");
            writeStringAttribute(group_, "time_resampling",
              "linear interpolation of raw E/B onto a uniform grid spanning the original committed time range");
            writeDoubleAttribute(group_, "plane_z_m", field.planeZ());
            writeDoubleAttribute(group_, "x_first_m", field.xFirst());
            writeDoubleAttribute(group_, "y_first_m", field.yFirst());
            writeDoubleAttribute(group_, "dx_m", field.dx());
            writeDoubleAttribute(group_, "dy_m", field.dy());
            writeDoubleAttribute(group_, "dt_s", grid.dt);
            writeDoubleAttribute(group_, "transverse_smoothing_m",
              config.transverseSmoothing);
            writeDoubleAttribute(group_, "effective_inverse_gamma2_vz2",
              denominatorCoefficient);
            writeDoubleAttribute(group_, "mean_gamma", summary.meanGamma);
            writeDoubleAttribute(group_, "relative_gamma_spread",
              summary.relativeGammaSpread);
            writeDoubleAttribute(group_, "rms_transverse_beta",
              summary.rmsTransverseBeta);
            writeDoubleAttribute(group_, "outside_charge_fraction",
              summary.outsideChargeFraction);
            writeUnsignedAttribute(group_, "nx", nx_);
            writeUnsignedAttribute(group_, "ny", ny_);
            writeUnsignedAttribute(group_, "samples", samples_);
            writeUnsignedAttribute(group_, "accepted_particles",
              summary.acceptedParticles);
            writeUnsignedAttribute(group_, "duplicate_particles",
              summary.duplicateParticles);
            writeUnsignedAttribute(group_, "invalid_particles",
              summary.invalidParticles);
            writeUnsignedAttribute(group_, "padding_time", config.paddingTime);
            writeUnsignedAttribute(group_, "padding_y", config.paddingY);
            writeUnsignedAttribute(group_, "padding_x", config.paddingX);
          }
        catch (...)
          {
            close(false);
            throw;
          }
      }

      ~OutputFile() { close(false); }

      void writeComponent(FieldComponent component, const double* padded,
                          const PaddedGrid& grid, double scale)
      {
        hid_t dataset = component <= ElectricZ ? electric_ : magnetic_;
        const std::size_t componentIndex = component <= ElectricZ ?
          static_cast<std::size_t>(component) :
          static_cast<std::size_t>(component - MagneticX);
        std::vector<double> plane(nx_ * ny_);
        for (std::size_t sample = 0; sample < samples_; ++sample)
          {
            for (std::size_t y = 0; y < ny_; ++y)
              for (std::size_t x = 0; x < nx_; ++x)
                {
                  const std::size_t paddedIndex =
                    ((sample + grid.timeOffset) * grid.ny +
                     (y + grid.yOffset)) * grid.nx + x + grid.xOffset;
                  plane[y * nx_ + x] = padded[paddedIndex] * scale;
                }
            hsize_t start[4] = {static_cast<hsize_t>(sample), 0, 0,
              static_cast<hsize_t>(componentIndex)};
            hsize_t count[4] = {1, static_cast<hsize_t>(ny_),
              static_cast<hsize_t>(nx_), 1};
            hid_t fileSpace = H5Dget_space(dataset);
            requireHandle(fileSpace, "Cannot inspect reconstructed dataset");
            requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
              start, NULL, count, NULL),
              "Cannot select reconstructed component");
            hsize_t memoryDimensions[2] = {
              static_cast<hsize_t>(ny_), static_cast<hsize_t>(nx_)
            };
            hid_t memorySpace = H5Screate_simple(2, memoryDimensions, NULL);
            requireHandle(memorySpace,
              "Cannot create reconstructed component memory space");
            const herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE,
              memorySpace, fileSpace, H5P_DEFAULT, &plane[0]);
            H5Sclose(memorySpace);
            H5Sclose(fileSpace);
            requireStatus(status, "Cannot write reconstructed component");
          }
      }

      void readSample(std::size_t sample,
                      std::vector<double>& electric,
                      std::vector<double>& magnetic) const
      {
        electric.resize(nx_ * ny_ * 3);
        magnetic.resize(nx_ * ny_ * 3);
        hsize_t start[4] = {static_cast<hsize_t>(sample), 0, 0, 0};
        hsize_t count[4] = {1, static_cast<hsize_t>(ny_),
          static_cast<hsize_t>(nx_), 3};
        hsize_t memoryDimensions[3] = {count[1], count[2], 3};
        hid_t memorySpace = H5Screate_simple(3, memoryDimensions, NULL);
        requireHandle(memorySpace, "Cannot create cleaned sample memory space");
        hid_t electricSpace = H5Dget_space(electric_);
        requireHandle(electricSpace, "Cannot inspect cleaned E dataset");
        requireStatus(H5Sselect_hyperslab(electricSpace, H5S_SELECT_SET,
          start, NULL, count, NULL), "Cannot select cleaned E sample");
        requireStatus(H5Dread(electric_, H5T_NATIVE_DOUBLE, memorySpace,
          electricSpace, H5P_DEFAULT, &electric[0]),
          "Cannot read cleaned E sample");
        H5Sclose(electricSpace);
        hid_t magneticSpace = H5Dget_space(magnetic_);
        requireHandle(magneticSpace, "Cannot inspect cleaned B dataset");
        requireStatus(H5Sselect_hyperslab(magneticSpace, H5S_SELECT_SET,
          start, NULL, count, NULL), "Cannot select cleaned B sample");
        requireStatus(H5Dread(magnetic_, H5T_NATIVE_DOUBLE, memorySpace,
          magneticSpace, H5P_DEFAULT, &magnetic[0]),
          "Cannot read cleaned B sample");
        H5Sclose(magneticSpace);
        H5Sclose(memorySpace);
      }

      void writePowerDiagnostics(const std::vector<double>& rawSigned,
                                 const std::vector<double>& cleanSigned,
                                 const std::vector<double>& rawForward,
                                 const std::vector<double>& cleanForward,
                                 double rawSignedEnergy,
                                 double cleanSignedEnergy,
                                 double rawEnergy,
                                 double cleanEnergy,
                                 double relativeChange)
      {
        writeVector("raw_signed_power_W", rawSigned);
        writeVector("cleaned_signed_power_W", cleanSigned);
        writeVector("raw_forward_power_W", rawForward);
        writeVector("cleaned_forward_power_W", cleanForward);
        writeDoubleAttribute(group_, "raw_signed_energy_J", rawSignedEnergy);
        writeDoubleAttribute(group_, "cleaned_signed_energy_J",
          cleanSignedEnergy);
        writeDoubleAttribute(group_, "raw_forward_energy_J", rawEnergy);
        writeDoubleAttribute(group_, "cleaned_forward_energy_J", cleanEnergy);
        writeDoubleAttribute(group_, "relative_forward_energy_change",
          relativeChange);
        writeStringAttribute(group_, "forward_power_definition",
          "integral max((E cross B)_z/mu0,0) dx dy over stored cell centres");
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
        if (complete_ >= 0) H5Dclose(complete_);
        if (magnetic_ >= 0) H5Dclose(magnetic_);
        if (electric_ >= 0) H5Dclose(electric_);
        if (time_ >= 0) H5Dclose(time_);
        if (group_ >= 0) H5Gclose(group_);
        H5Fclose(file_);
        complete_ = magnetic_ = electric_ = time_ = group_ = file_ = -1;
      }

    private:
      void writeVector(const char* name, const std::vector<double>& values)
      {
        hsize_t dimensions[1] = {static_cast<hsize_t>(values.size())};
        hid_t space = H5Screate_simple(1, dimensions, NULL);
        requireHandle(space, "Cannot create diagnostic vector space");
        hid_t dataset = H5Dcreate2(group_, name, H5T_IEEE_F64LE, space,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(space);
        requireHandle(dataset, std::string("Cannot create diagnostic: ") + name);
        const herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE,
          H5S_ALL, H5S_ALL, H5P_DEFAULT, &values[0]);
        H5Dclose(dataset);
        requireStatus(status, std::string("Cannot write diagnostic: ") + name);
      }

      OutputFile(const OutputFile&);
      OutputFile& operator=(const OutputFile&);
      hid_t file_;
      hid_t group_;
      hid_t time_;
      hid_t electric_;
      hid_t magnetic_;
      hid_t complete_;
      std::size_t nx_;
      std::size_t ny_;
      std::size_t samples_;
    };

    void loadRawComponent(const FieldInput& field,
                          const PaddedGrid& grid,
                          FieldComponent component,
                          double* destination)
    {
      std::fill(destination, destination + grid.realSize, 0.0);
      std::vector<FieldPoint> sample;
      for (std::size_t time = 0; time < field.samples(); ++time)
        {
          field.readInterpolated(grid.outputTimes[time], sample);
          for (std::size_t y = 0; y < field.ny(); ++y)
            for (std::size_t x = 0; x < field.nx(); ++x)
              {
                const FieldPoint& point = sample[y * field.nx() + x];
                const double value = component <= ElectricZ ?
                  point.electric[static_cast<unsigned int>(component)] :
                  point.magnetic[static_cast<unsigned int>(component - MagneticX)];
                const std::size_t padded =
                  ((time + grid.timeOffset) * grid.ny +
                   (y + grid.yOffset)) * grid.nx + x + grid.xOffset;
                destination[padded] = value;
              }
        }
    }

    double integrateTrapezoid(const std::vector<double>& time,
                              const std::vector<double>& value)
    {
      double result = 0.0;
      for (std::size_t index = 1; index < time.size(); ++index)
        result += 0.5 * (value[index - 1] + value[index]) *
          (time[index] - time[index - 1]);
      return result;
    }
  }

  ReconstructionSummary::ReconstructionSummary()
    : inputRecords(0), acceptedParticles(0), duplicateParticles(0),
      invalidParticles(0), outsideChargeFraction(0.0), meanGamma(0.0),
      relativeGammaSpread(0.0), rmsTransverseBeta(0.0),
      rawSignedEnergy(0.0), cleanedSignedEnergy(0.0),
      rawForwardEnergy(0.0), cleanedForwardEnergy(0.0),
      relativeForwardEnergyChange(0.0)
  {}

  ReconstructionSummary reconstructField(
      const ReconstructionConfig& config)
  {
    FieldInput field(config.fieldFile, config.requireComplete);
    ParticleInput particles(config.particleFile, config.requireComplete,
                            field.planeZ());
    const PaddedGrid grid = makePaddedGrid(field, config);

    ParticleStatistics statistics;
    unsigned long long duplicates = 0;
    unsigned long long invalid = 0;
    forEachUniqueParticle(particles, config.particleReadChunk, field.planeZ(),
      [&](const ProjectedParticle& particle)
      {
        const long double weight = std::abs(particle.charge);
        const long double betaX = particle.velocity[0] / kC;
        const long double betaY = particle.velocity[1] / kC;
        statistics.absoluteCharge += weight;
        statistics.gamma += weight * particle.gamma;
        statistics.gammaSquared += weight * particle.gamma * particle.gamma;
        statistics.betaTransverseSquared += weight *
          (betaX * betaX + betaY * betaY);
        statistics.denominatorCoefficient += weight /
          (particle.gamma * particle.gamma *
           particle.velocity[2] * particle.velocity[2]);
        ++statistics.accepted;
      }, duplicates, invalid);
    if (statistics.accepted == 0 || !(statistics.absoluteCharge > 0.0L))
      throw std::runtime_error("No valid downstream particle crossings found");

    ReconstructionSummary summary;
    summary.inputRecords = particles.committed();
    summary.acceptedParticles = statistics.accepted;
    summary.duplicateParticles = duplicates;
    summary.invalidParticles = invalid;
    summary.meanGamma = static_cast<double>(statistics.gamma /
      statistics.absoluteCharge);
    const long double gammaVariance = std::max(0.0L,
      statistics.gammaSquared / statistics.absoluteCharge -
      static_cast<long double>(summary.meanGamma) * summary.meanGamma);
    summary.relativeGammaSpread = static_cast<double>(
      std::sqrt(gammaVariance) / summary.meanGamma);
    summary.rmsTransverseBeta = static_cast<double>(std::sqrt(
      statistics.betaTransverseSquared / statistics.absoluteCharge));
    const double denominatorCoefficient = static_cast<double>(
      statistics.denominatorCoefficient / statistics.absoluteCharge);
    if (summary.relativeGammaSpread > config.maximumRelativeGammaSpread)
      throw std::runtime_error(
        "Particle gamma spread exceeds the configured common-kernel limit");
    if (summary.rmsTransverseBeta > config.maximumRmsTransverseBeta)
      throw std::runtime_error(
        "Particle transverse velocity exceeds the longitudinal model limit");

    if (fftw_init_threads() == 0)
      throw std::runtime_error("Cannot initialize FFTW threads");
    fftw_plan_with_nthreads(static_cast<int>(config.fftThreads));
    double* sourceQ = fftw_alloc_real(grid.realSize);
    double* sourceQV = fftw_alloc_real(grid.realSize);
    double* sourceQGammaV = fftw_alloc_real(grid.realSize);
    fftw_complex* spectrumQ = fftw_alloc_complex(grid.complexSize);
    fftw_complex* spectrumQV = fftw_alloc_complex(grid.complexSize);
    fftw_complex* spectrumQGammaV = fftw_alloc_complex(grid.complexSize);
    double* workReal = fftw_alloc_real(grid.realSize);
    fftw_complex* workSpectrum = fftw_alloc_complex(grid.complexSize);
    fftw_complex* background = fftw_alloc_complex(grid.complexSize);
    if (!sourceQ || !sourceQV || !sourceQGammaV || !spectrumQ ||
        !spectrumQV || !spectrumQGammaV || !workReal || !workSpectrum ||
        !background)
      throw std::bad_alloc();
    std::fill(sourceQ, sourceQ + grid.realSize, 0.0);
    std::fill(sourceQV, sourceQV + grid.realSize, 0.0);
    std::fill(sourceQGammaV, sourceQGammaV + grid.realSize, 0.0);

    long double outsideCharge = 0.0L;
    duplicates = invalid = 0;
    forEachUniqueParticle(particles, config.particleReadChunk, field.planeZ(),
      [&](const ProjectedParticle& particle)
      {
        if (!depositParticle(particle, grid, sourceQ, sourceQV,
                             sourceQGammaV))
          outsideCharge += std::abs(particle.charge);
      }, duplicates, invalid);
    summary.outsideChargeFraction = static_cast<double>(
      outsideCharge / statistics.absoluteCharge);
    if (summary.outsideChargeFraction >
        config.maximumOutsideChargeFraction)
      throw std::runtime_error(
        "Too much particle charge lies outside the padded reconstruction grid");

    fftw_plan planQ = fftw_plan_dft_r2c_3d(
      static_cast<int>(grid.nt), static_cast<int>(grid.ny),
      static_cast<int>(grid.nx), sourceQ, spectrumQ, FFTW_ESTIMATE);
    fftw_plan planQV = fftw_plan_dft_r2c_3d(
      static_cast<int>(grid.nt), static_cast<int>(grid.ny),
      static_cast<int>(grid.nx), sourceQV, spectrumQV, FFTW_ESTIMATE);
    fftw_plan planQGammaV = fftw_plan_dft_r2c_3d(
      static_cast<int>(grid.nt), static_cast<int>(grid.ny),
      static_cast<int>(grid.nx), sourceQGammaV, spectrumQGammaV,
      FFTW_ESTIMATE);
    fftw_plan rawForward = fftw_plan_dft_r2c_3d(
      static_cast<int>(grid.nt), static_cast<int>(grid.ny),
      static_cast<int>(grid.nx), workReal, workSpectrum, FFTW_ESTIMATE);
    fftw_plan cleanedBackward = fftw_plan_dft_c2r_3d(
      static_cast<int>(grid.nt), static_cast<int>(grid.ny),
      static_cast<int>(grid.nx), workSpectrum, workReal, FFTW_ESTIMATE);
    if (!planQ || !planQV || !planQGammaV || !rawForward || !cleanedBackward)
      throw std::runtime_error("Cannot create FFTW reconstruction plans");
    fftw_execute(planQ);
    fftw_execute(planQV);
    fftw_execute(planQGammaV);
    fftw_free(sourceQ);
    fftw_free(sourceQV);
    fftw_free(sourceQGammaV);
    sourceQ = sourceQV = sourceQGammaV = NULL;

    OutputFile output(config, field, particles, grid, summary,
                      denominatorCoefficient);
    const double inverseTransform = 1.0 /
      static_cast<double>(grid.realSize);
    for (int componentValue = static_cast<int>(ElectricX);
         componentValue <= static_cast<int>(MagneticY); ++componentValue)
      {
        const FieldComponent component =
          static_cast<FieldComponent>(componentValue);
        loadRawComponent(field, grid, component, workReal);
        fftw_execute(rawForward);
        backgroundSpectrum(component, grid, spectrumQ, spectrumQV,
          spectrumQGammaV, denominatorCoefficient,
          config.transverseSmoothing, background);
        for (std::size_t index = 0; index < grid.complexSize; ++index)
          {
            workSpectrum[index][0] -= background[index][0];
            workSpectrum[index][1] -= background[index][1];
          }
        fftw_execute(cleanedBackward);
        output.writeComponent(component, workReal, grid, inverseTransform);
      }
    loadRawComponent(field, grid, MagneticZ, workReal);
    output.writeComponent(MagneticZ, workReal, grid, 1.0);

    std::vector<double> rawSigned(field.samples(), 0.0);
    std::vector<double> cleanSigned(field.samples(), 0.0);
    std::vector<double> rawForwardPower(field.samples(), 0.0);
    std::vector<double> cleanForwardPower(field.samples(), 0.0);
    std::vector<FieldPoint> raw;
    std::vector<double> cleanElectric;
    std::vector<double> cleanMagnetic;
    const double area = field.dx() * field.dy();
    for (std::size_t sample = 0; sample < field.samples(); ++sample)
      {
        field.readInterpolated(grid.outputTimes[sample], raw);
        output.readSample(sample, cleanElectric, cleanMagnetic);
        for (std::size_t point = 0; point < raw.size(); ++point)
          {
            const double rawFlux = (raw[point].electric[0] *
              raw[point].magnetic[1] - raw[point].electric[1] *
              raw[point].magnetic[0]) / kMu0;
            const double cleanFlux = (cleanElectric[point * 3] *
              cleanMagnetic[point * 3 + 1] -
              cleanElectric[point * 3 + 1] *
              cleanMagnetic[point * 3]) / kMu0;
            rawSigned[sample] += rawFlux * area;
            cleanSigned[sample] += cleanFlux * area;
            rawForwardPower[sample] += std::max(0.0, rawFlux) * area;
            cleanForwardPower[sample] += std::max(0.0, cleanFlux) * area;
          }
      }
    summary.rawSignedEnergy = integrateTrapezoid(grid.outputTimes, rawSigned);
    summary.cleanedSignedEnergy = integrateTrapezoid(
      grid.outputTimes, cleanSigned);
    summary.rawForwardEnergy = integrateTrapezoid(
      grid.outputTimes, rawForwardPower);
    summary.cleanedForwardEnergy = integrateTrapezoid(
      grid.outputTimes, cleanForwardPower);
    summary.relativeForwardEnergyChange = summary.rawForwardEnergy == 0.0 ?
      0.0 : (summary.cleanedForwardEnergy - summary.rawForwardEnergy) /
        summary.rawForwardEnergy;
    output.writePowerDiagnostics(rawSigned, cleanSigned, rawForwardPower,
      cleanForwardPower, summary.rawSignedEnergy,
      summary.cleanedSignedEnergy, summary.rawForwardEnergy,
      summary.cleanedForwardEnergy, summary.relativeForwardEnergyChange);
    output.close(true);

    fftw_destroy_plan(planQ);
    fftw_destroy_plan(planQV);
    fftw_destroy_plan(planQGammaV);
    fftw_destroy_plan(rawForward);
    fftw_destroy_plan(cleanedBackward);
    fftw_free(spectrumQ);
    fftw_free(spectrumQV);
    fftw_free(spectrumQGammaV);
    fftw_free(workReal);
    fftw_free(workSpectrum);
    fftw_free(background);
    fftw_cleanup_threads();
    return summary;
  }
}
