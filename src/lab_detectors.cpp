#include "lab_detectors.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "hdf5.h"

#include "runtime_util.h"

namespace fel
{
  namespace
  {
    const int PARTICLE_EVENT_TAG = 8701;
    const int FIELD_SAMPLE_TAG = 8702;

    void requireHandle(hid_t handle, const std::string& message)
    {
      if (handle < 0) throw std::runtime_error(message);
    }

    void requireStatus(herr_t status, const std::string& message)
    {
      if (status < 0) throw std::runtime_error(message);
    }

    int checkedByteCount(std::size_t records, std::size_t recordBytes,
                         const char* description)
    {
      if (recordBytes == 0 ||
          records > static_cast<std::size_t>(INT_MAX) / recordBytes)
        throw std::overflow_error(std::string(description) +
          " exceeds one MPI message; reduce the detector batch size");
      return static_cast<int>(records * recordBytes);
    }

    void writeStringAttribute(hid_t object, const char* name,
                              const std::string& value)
    {
      hid_t type = H5Tcopy(H5T_C_S1);
      requireHandle(type, std::string("Cannot create detector string type: ") + name);
      try
        {
          requireStatus(H5Tset_size(type, value.size() + 1),
            std::string("Cannot size detector string attribute: ") + name);
          requireStatus(H5Tset_strpad(type, H5T_STR_NULLTERM),
            std::string("Cannot set detector string padding: ") + name);
          hid_t space = H5Screate(H5S_SCALAR);
          requireHandle(space, "Cannot create detector attribute space");
          hid_t attribute = H5Acreate2(object, name, type, space,
            H5P_DEFAULT, H5P_DEFAULT);
          H5Sclose(space);
          requireHandle(attribute,
            std::string("Cannot create detector attribute: ") + name);
          const herr_t status = H5Awrite(attribute, type, value.c_str());
          H5Aclose(attribute);
          requireStatus(status,
            std::string("Cannot write detector attribute: ") + name);
        }
      catch (...)
        {
          H5Tclose(type);
          throw;
        }
      H5Tclose(type);
    }

    void writeDoubleAttribute(hid_t object, const char* name, double value)
    {
      hid_t space = H5Screate(H5S_SCALAR);
      requireHandle(space, "Cannot create detector scalar attribute space");
      hid_t attribute = H5Acreate2(object, name, H5T_IEEE_F64LE, space,
        H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(space);
      requireHandle(attribute,
        std::string("Cannot create detector attribute: ") + name);
      const herr_t status = H5Awrite(attribute, H5T_NATIVE_DOUBLE, &value);
      H5Aclose(attribute);
      requireStatus(status,
        std::string("Cannot write detector attribute: ") + name);
    }

    void writeUnsignedAttribute(hid_t object, const char* name,
                                std::uint64_t value)
    {
      hid_t space = H5Screate(H5S_SCALAR);
      requireHandle(space, "Cannot create detector integer attribute space");
      hid_t attribute = H5Acreate2(object, name, H5T_STD_U64LE, space,
        H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(space);
      requireHandle(attribute,
        std::string("Cannot create detector attribute: ") + name);
      const herr_t status = H5Awrite(attribute, H5T_NATIVE_UINT64, &value);
      H5Aclose(attribute);
      requireStatus(status,
        std::string("Cannot write detector attribute: ") + name);
    }

    void writeUnsignedScalar(hid_t dataset, std::uint64_t value,
                             const char* description)
    {
      requireStatus(H5Dwrite(dataset, H5T_NATIVE_UINT64,
        H5S_ALL, H5S_ALL, H5P_DEFAULT, &value), description);
    }

    void writeByteScalar(hid_t dataset, unsigned char value,
                         const char* description)
    {
      requireStatus(H5Dwrite(dataset, H5T_NATIVE_UCHAR,
        H5S_ALL, H5S_ALL, H5P_DEFAULT, &value), description);
    }

    struct FieldPointRecord
    {
      double electric[3];
      double magnetic[3];
    };

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

    struct ParticleCrossingWire
    {
      std::uint32_t detectorIndex;
      ParticlePlaneRecord record;
    };

    hid_t createFieldMemoryType()
    {
      hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(FieldPointRecord));
      requireHandle(type, "Cannot create field-detector memory datatype");
      hsize_t vectorSize[1] = {3};
      hid_t vectorType = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, vectorSize);
      if (vectorType < 0)
        {
          H5Tclose(type);
          throw std::runtime_error(
            "Cannot create field-detector vector datatype");
        }
      try
        {
          requireStatus(H5Tinsert(type, "electric_V_per_m",
            HOFFSET(FieldPointRecord, electric), vectorType),
            "Cannot add field-detector electric field");
          requireStatus(H5Tinsert(type, "magnetic_T",
            HOFFSET(FieldPointRecord, magnetic), vectorType),
            "Cannot add field-detector magnetic field");
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

    hid_t createFieldFileType()
    {
      hid_t type = H5Tcreate(H5T_COMPOUND, 6 * 8);
      requireHandle(type, "Cannot create field-detector file datatype");
      hsize_t vectorSize[1] = {3};
      hid_t vectorType = H5Tarray_create2(H5T_IEEE_F64LE, 1, vectorSize);
      if (vectorType < 0)
        {
          H5Tclose(type);
          throw std::runtime_error(
            "Cannot create field-detector file vector datatype");
        }
      try
        {
          requireStatus(H5Tinsert(type, "electric_V_per_m", 0, vectorType),
            "Cannot add field-detector file electric field");
          requireStatus(H5Tinsert(type, "magnetic_T", 3 * 8, vectorType),
            "Cannot add field-detector file magnetic field");
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

    hid_t createParticleMemoryType()
    {
      hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(ParticlePlaneRecord));
      requireHandle(type, "Cannot create particle-detector memory datatype");
      hsize_t vectorSize[1] = {3};
      hid_t vectorType = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, vectorSize);
      if (vectorType < 0)
        {
          H5Tclose(type);
          throw std::runtime_error(
            "Cannot create particle-detector vector datatype");
        }
      try
        {
          requireStatus(H5Tinsert(type, "particle_id",
            HOFFSET(ParticlePlaneRecord, particleId), H5T_NATIVE_UINT64),
            "Cannot add particle-detector particle_id");
          requireStatus(H5Tinsert(type, "source_id",
            HOFFSET(ParticlePlaneRecord, sourceId), H5T_NATIVE_UINT64),
            "Cannot add particle-detector source_id");
          requireStatus(H5Tinsert(type, "time_s",
            HOFFSET(ParticlePlaneRecord, time), H5T_NATIVE_DOUBLE),
            "Cannot add particle-detector time");
          requireStatus(H5Tinsert(type, "position_m",
            HOFFSET(ParticlePlaneRecord, position), vectorType),
            "Cannot add particle-detector position");
          requireStatus(H5Tinsert(type, "proper_velocity",
            HOFFSET(ParticlePlaneRecord, properVelocity), vectorType),
            "Cannot add particle-detector proper velocity");
          requireStatus(H5Tinsert(type, "charge_C",
            HOFFSET(ParticlePlaneRecord, charge), H5T_NATIVE_DOUBLE),
            "Cannot add particle-detector charge");
          requireStatus(H5Tinsert(type, "mass_kg",
            HOFFSET(ParticlePlaneRecord, mass), H5T_NATIVE_DOUBLE),
            "Cannot add particle-detector mass");
          requireStatus(H5Tinsert(type, "weight",
            HOFFSET(ParticlePlaneRecord, weight), H5T_NATIVE_DOUBLE),
            "Cannot add particle-detector weight");
          requireStatus(H5Tinsert(type, "crossing_direction",
            HOFFSET(ParticlePlaneRecord, direction), H5T_NATIVE_INT32),
            "Cannot add particle-detector crossing direction");
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

    hid_t createParticleFileType()
    {
      const std::size_t u64 = 8;
      const std::size_t f64 = 8;
      const std::size_t fileSize = 2 * u64 + 10 * f64 + 4;
      hid_t type = H5Tcreate(H5T_COMPOUND, fileSize);
      requireHandle(type, "Cannot create particle-detector file datatype");
      hsize_t vectorSize[1] = {3};
      hid_t vectorType = H5Tarray_create2(H5T_IEEE_F64LE, 1, vectorSize);
      if (vectorType < 0)
        {
          H5Tclose(type);
          throw std::runtime_error(
            "Cannot create particle-detector file vector datatype");
        }
      try
        {
          std::size_t offset = 0;
          requireStatus(H5Tinsert(type, "particle_id", offset, H5T_STD_U64LE),
            "Cannot add particle-detector file particle_id");
          offset += u64;
          requireStatus(H5Tinsert(type, "source_id", offset, H5T_STD_U64LE),
            "Cannot add particle-detector file source_id");
          offset += u64;
          requireStatus(H5Tinsert(type, "time_s", offset, H5T_IEEE_F64LE),
            "Cannot add particle-detector file time");
          offset += f64;
          requireStatus(H5Tinsert(type, "position_m", offset, vectorType),
            "Cannot add particle-detector file position");
          offset += 3 * f64;
          requireStatus(H5Tinsert(type, "proper_velocity", offset, vectorType),
            "Cannot add particle-detector file proper velocity");
          offset += 3 * f64;
          requireStatus(H5Tinsert(type, "charge_C", offset, H5T_IEEE_F64LE),
            "Cannot add particle-detector file charge");
          offset += f64;
          requireStatus(H5Tinsert(type, "mass_kg", offset, H5T_IEEE_F64LE),
            "Cannot add particle-detector file mass");
          offset += f64;
          requireStatus(H5Tinsert(type, "weight", offset, H5T_IEEE_F64LE),
            "Cannot add particle-detector file weight");
          offset += f64;
          requireStatus(H5Tinsert(type, "crossing_direction", offset,
            H5T_STD_I32LE),
            "Cannot add particle-detector file crossing direction");
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

    class FieldPlaneWriter
    {
    public:
      FieldPlaneWriter()
        : open_(false), nx_(0), ny_(0), bufferLimit_(0), committed_(0),
          times_(), fields_(), file_(-1), group_(-1), timeDataset_(-1),
          fieldDataset_(-1), committedDataset_(-1), completeDataset_(-1),
          memoryType_(-1), fileType_(-1)
      {}

      ~FieldPlaneWriter()
      {
        try { close(false); }
        catch (...) { closeHandles(); }
      }

      void open(const std::string& filename,
                const FieldDetectorPlaneConfig& config,
                const EBGridGeometry& geometry,
                const FieldVector<Double>& origin)
      {
        if (open_) throw std::runtime_error("Field detector is already open");
        nx_ = geometry.nx;
        ny_ = geometry.ny;
        bufferLimit_ = config.bufferSamples;
        committed_ = 0;
        times_.reserve(bufferLimit_);
        if (nx_ != 0 && ny_ > std::numeric_limits<std::size_t>::max() / nx_)
          throw std::overflow_error("Field detector plane dimensions overflow");
        const std::size_t points = nx_ * ny_;
        if (bufferLimit_ > std::numeric_limits<std::size_t>::max() / points)
          throw std::overflow_error("Field detector buffer dimensions overflow");
        fields_.reserve(bufferLimit_ * points);

        file_ = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,
          H5P_DEFAULT, H5P_DEFAULT);
        requireHandle(file_, "Cannot create field-detector file: " + filename);
        try
          {
            group_ = H5Gcreate2(file_, "/field_plane", H5P_DEFAULT,
              H5P_DEFAULT, H5P_DEFAULT);
            requireHandle(group_, "Cannot create field-detector group");
            memoryType_ = createFieldMemoryType();
            fileType_ = createFieldFileType();

            hsize_t timeInitial[1] = {0};
            hsize_t timeMaximum[1] = {H5S_UNLIMITED};
            hid_t timeSpace = H5Screate_simple(1, timeInitial, timeMaximum);
            requireHandle(timeSpace, "Cannot create field-detector time space");
            hid_t timeCreation = H5Pcreate(H5P_DATASET_CREATE);
            requireHandle(timeCreation,
              "Cannot create field-detector time properties");
            hsize_t timeChunk[1] = {
              static_cast<hsize_t>(std::max<std::size_t>(1, bufferLimit_))
            };
            requireStatus(H5Pset_chunk(timeCreation, 1, timeChunk),
              "Cannot set field-detector time chunk");
            timeDataset_ = H5Dcreate2(group_, "time_s", H5T_IEEE_F64LE,
              timeSpace, H5P_DEFAULT, timeCreation, H5P_DEFAULT);
            H5Pclose(timeCreation);
            H5Sclose(timeSpace);
            requireHandle(timeDataset_,
              "Cannot create field-detector time dataset");

            hsize_t fieldInitial[3] = {
              0, static_cast<hsize_t>(ny_), static_cast<hsize_t>(nx_)
            };
            hsize_t fieldMaximum[3] = {
              H5S_UNLIMITED, static_cast<hsize_t>(ny_),
              static_cast<hsize_t>(nx_)
            };
            hid_t fieldSpace = H5Screate_simple(3, fieldInitial, fieldMaximum);
            requireHandle(fieldSpace,
              "Cannot create field-detector sample space");
            hid_t creation = H5Pcreate(H5P_DATASET_CREATE);
            requireHandle(creation,
              "Cannot create field-detector sample properties");
            hsize_t chunk[3] = {
              1,
              static_cast<hsize_t>(std::min<std::size_t>(ny_, 256)),
              static_cast<hsize_t>(std::min<std::size_t>(nx_, 256))
            };
            requireStatus(H5Pset_chunk(creation, 3, chunk),
              "Cannot set field-detector sample chunk");
            if (config.compression > 0)
              {
                requireStatus(H5Pset_shuffle(creation),
                  "Cannot enable field-detector byte shuffle");
                requireStatus(H5Pset_deflate(creation, config.compression),
                  "Cannot enable field-detector compression");
              }
            fieldDataset_ = H5Dcreate2(group_, "fields", fileType_,
              fieldSpace, H5P_DEFAULT, creation, H5P_DEFAULT);
            H5Pclose(creation);
            H5Sclose(fieldSpace);
            requireHandle(fieldDataset_,
              "Cannot create field-detector sample dataset");

            hid_t scalar = H5Screate(H5S_SCALAR);
            requireHandle(scalar,
              "Cannot create field-detector scalar space");
            committedDataset_ = H5Dcreate2(group_, "committed_samples",
              H5T_STD_U64LE, scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            completeDataset_ = H5Dcreate2(group_, "complete", H5T_STD_U8LE,
              scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(scalar);
            requireHandle(committedDataset_,
              "Cannot create field-detector commit marker");
            requireHandle(completeDataset_,
              "Cannot create field-detector completion marker");
            writeUnsignedScalar(committedDataset_, 0,
              "Cannot initialize field-detector commit marker");
            writeByteScalar(completeDataset_, 0,
              "Cannot initialize field-detector completion marker");

            writeUnsignedAttribute(group_, "format_version", 1);
            writeStringAttribute(group_, "name", config.name);
            writeStringAttribute(group_, "frame", "laboratory");
            writeStringAttribute(group_, "write_mode", "rank-zero-throughput");
            writeStringAttribute(group_, "x_y_sampling", "cell-centres");
            writeStringAttribute(group_, "magnetic_time_stagger",
              "Yee B is sampled at its stored leapfrog half time");
            writeDoubleAttribute(group_, "plane_z_m", config.z);
            writeDoubleAttribute(group_, "x_first_m", origin[0] + 0.5 * geometry.dx);
            writeDoubleAttribute(group_, "y_first_m", origin[1] + 0.5 * geometry.dy);
            writeDoubleAttribute(group_, "dx_m", geometry.dx);
            writeDoubleAttribute(group_, "dy_m", geometry.dy);
            writeUnsignedAttribute(group_, "nx", nx_);
            writeUnsignedAttribute(group_, "ny", ny_);
            writeUnsignedAttribute(group_, "compression_level", config.compression);
            open_ = true;
          }
        catch (...)
          {
            closeHandles();
            times_.clear();
            fields_.clear();
            throw;
          }
      }

      void append(double time, const std::vector<FieldPointRecord>& plane)
      {
        if (!open_) throw std::runtime_error("Field detector is not open");
        if (!std::isfinite(time) || plane.size() != nx_ * ny_)
          throw std::invalid_argument("Invalid field-detector plane sample");
        times_.push_back(time);
        fields_.insert(fields_.end(), plane.begin(), plane.end());
        if (times_.size() >= bufferLimit_) flush();
      }

      void flush()
      {
        if (!open_ || times_.empty()) return;
        const hsize_t samples = static_cast<hsize_t>(times_.size());
        const hsize_t timeStart[1] = {static_cast<hsize_t>(committed_)};
        const hsize_t timeCount[1] = {samples};
        const hsize_t timeExtent[1] = {timeStart[0] + samples};
        requireStatus(H5Dset_extent(timeDataset_, timeExtent),
          "Cannot extend field-detector time dataset");
        hid_t timeFileSpace = H5Dget_space(timeDataset_);
        requireHandle(timeFileSpace,
          "Cannot select field-detector time file space");
        requireStatus(H5Sselect_hyperslab(timeFileSpace, H5S_SELECT_SET,
          timeStart, NULL, timeCount, NULL),
          "Cannot select field-detector time batch");
        hid_t timeMemorySpace = H5Screate_simple(1, timeCount, NULL);
        requireHandle(timeMemorySpace,
          "Cannot create field-detector time memory space");
        const herr_t timeStatus = H5Dwrite(timeDataset_, H5T_NATIVE_DOUBLE,
          timeMemorySpace, timeFileSpace, H5P_DEFAULT, &times_[0]);
        H5Sclose(timeMemorySpace);
        H5Sclose(timeFileSpace);
        requireStatus(timeStatus, "Cannot write field-detector times");

        const hsize_t fieldStart[3] = {
          static_cast<hsize_t>(committed_), 0, 0
        };
        const hsize_t fieldCount[3] = {
          samples, static_cast<hsize_t>(ny_), static_cast<hsize_t>(nx_)
        };
        const hsize_t fieldExtent[3] = {
          fieldStart[0] + samples,
          static_cast<hsize_t>(ny_), static_cast<hsize_t>(nx_)
        };
        requireStatus(H5Dset_extent(fieldDataset_, fieldExtent),
          "Cannot extend field-detector field dataset");
        hid_t fieldFileSpace = H5Dget_space(fieldDataset_);
        requireHandle(fieldFileSpace,
          "Cannot select field-detector field file space");
        requireStatus(H5Sselect_hyperslab(fieldFileSpace, H5S_SELECT_SET,
          fieldStart, NULL, fieldCount, NULL),
          "Cannot select field-detector field batch");
        hid_t fieldMemorySpace = H5Screate_simple(3, fieldCount, NULL);
        requireHandle(fieldMemorySpace,
          "Cannot create field-detector field memory space");
        const herr_t fieldStatus = H5Dwrite(fieldDataset_, memoryType_,
          fieldMemorySpace, fieldFileSpace, H5P_DEFAULT, &fields_[0]);
        H5Sclose(fieldMemorySpace);
        H5Sclose(fieldFileSpace);
        requireStatus(fieldStatus, "Cannot write field-detector batch");

        committed_ += static_cast<std::uint64_t>(times_.size());
        writeUnsignedScalar(committedDataset_, committed_,
          "Cannot update field-detector commit marker");
        times_.clear();
        fields_.clear();
      }

      void close(bool completed)
      {
        if (!open_)
          {
            closeHandles();
            return;
          }
        flush();
        writeByteScalar(completeDataset_, completed ? 1 : 0,
          "Cannot write field-detector completion marker");
        requireStatus(H5Fflush(file_, H5F_SCOPE_GLOBAL),
          "Cannot finalize field-detector file");
        open_ = false;
        closeHandles();
        times_.clear();
        fields_.clear();
      }

    private:
      void closeHandles() noexcept
      {
        if (timeDataset_ >= 0) H5Dclose(timeDataset_);
        if (fieldDataset_ >= 0) H5Dclose(fieldDataset_);
        if (committedDataset_ >= 0) H5Dclose(committedDataset_);
        if (completeDataset_ >= 0) H5Dclose(completeDataset_);
        if (memoryType_ >= 0) H5Tclose(memoryType_);
        if (fileType_ >= 0) H5Tclose(fileType_);
        if (group_ >= 0) H5Gclose(group_);
        if (file_ >= 0) H5Fclose(file_);
        timeDataset_ = fieldDataset_ = committedDataset_ = completeDataset_ = -1;
        memoryType_ = fileType_ = group_ = file_ = -1;
      }

      bool open_;
      std::size_t nx_;
      std::size_t ny_;
      std::size_t bufferLimit_;
      std::uint64_t committed_;
      std::vector<double> times_;
      std::vector<FieldPointRecord> fields_;
      hid_t file_;
      hid_t group_;
      hid_t timeDataset_;
      hid_t fieldDataset_;
      hid_t committedDataset_;
      hid_t completeDataset_;
      hid_t memoryType_;
      hid_t fileType_;
    };

    class ParticlePlaneWriter
    {
    public:
      ParticlePlaneWriter()
        : open_(false), bufferLimit_(0), committed_(0), buffer_(), file_(-1),
          group_(-1), records_(-1), committedDataset_(-1),
          completeDataset_(-1), memoryType_(-1), fileType_(-1)
      {}

      ~ParticlePlaneWriter()
      {
        try { close(false); }
        catch (...) { closeHandles(); }
      }

      void open(const std::string& filename,
                const ParticleDetectorPlaneConfig& config,
                int mpiSize)
      {
        if (open_) throw std::runtime_error("Particle detector is already open");
        bufferLimit_ = config.bufferRecords;
        committed_ = 0;
        buffer_.reserve(bufferLimit_);
        file_ = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,
          H5P_DEFAULT, H5P_DEFAULT);
        requireHandle(file_, "Cannot create particle-detector file: " + filename);
        try
          {
            group_ = H5Gcreate2(file_, "/particle_plane", H5P_DEFAULT,
              H5P_DEFAULT, H5P_DEFAULT);
            requireHandle(group_, "Cannot create particle-detector group");
            memoryType_ = createParticleMemoryType();
            fileType_ = createParticleFileType();

            hsize_t initial[1] = {0};
            hsize_t maximum[1] = {H5S_UNLIMITED};
            hid_t space = H5Screate_simple(1, initial, maximum);
            requireHandle(space, "Cannot create particle-detector record space");
            hid_t creation = H5Pcreate(H5P_DATASET_CREATE);
            requireHandle(creation,
              "Cannot create particle-detector record properties");
            hsize_t chunk[1] = {static_cast<hsize_t>(bufferLimit_)};
            requireStatus(H5Pset_chunk(creation, 1, chunk),
              "Cannot set particle-detector record chunk");
            if (config.compression > 0)
              {
                requireStatus(H5Pset_shuffle(creation),
                  "Cannot enable particle-detector byte shuffle");
                requireStatus(H5Pset_deflate(creation, config.compression),
                  "Cannot enable particle-detector compression");
              }
            records_ = H5Dcreate2(group_, "records", fileType_, space,
              H5P_DEFAULT, creation, H5P_DEFAULT);
            H5Pclose(creation);
            H5Sclose(space);
            requireHandle(records_,
              "Cannot create particle-detector record dataset");

            hid_t scalar = H5Screate(H5S_SCALAR);
            requireHandle(scalar,
              "Cannot create particle-detector scalar space");
            committedDataset_ = H5Dcreate2(group_, "committed_records",
              H5T_STD_U64LE, scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            completeDataset_ = H5Dcreate2(group_, "complete", H5T_STD_U8LE,
              scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(scalar);
            requireHandle(committedDataset_,
              "Cannot create particle-detector commit marker");
            requireHandle(completeDataset_,
              "Cannot create particle-detector completion marker");
            writeUnsignedScalar(committedDataset_, 0,
              "Cannot initialize particle-detector commit marker");
            writeByteScalar(completeDataset_, 0,
              "Cannot initialize particle-detector completion marker");

            writeUnsignedAttribute(group_, "format_version", 1);
            writeStringAttribute(group_, "name", config.name);
            writeStringAttribute(group_, "frame", "laboratory");
            writeStringAttribute(group_, "write_mode", "rank-zero-throughput");
            writeStringAttribute(group_, "proper_velocity_unit", "gamma*v/c");
            writeStringAttribute(group_, "crossing_direction_contract",
              "+1 downstream, -1 upstream");
            writeDoubleAttribute(group_, "plane_z_m", config.z);
            writeUnsignedAttribute(group_, "mpi_size", mpiSize);
            writeUnsignedAttribute(group_, "compression_level", config.compression);
            open_ = true;
          }
        catch (...)
          {
            closeHandles();
            buffer_.clear();
            throw;
          }
      }

      void append(const ParticlePlaneRecord& record)
      {
        if (!open_) throw std::runtime_error("Particle detector is not open");
        buffer_.push_back(record);
        if (buffer_.size() >= bufferLimit_) flush();
      }

      void flush()
      {
        if (!open_ || buffer_.empty()) return;
        const hsize_t start[1] = {static_cast<hsize_t>(committed_)};
        const hsize_t count[1] = {static_cast<hsize_t>(buffer_.size())};
        const hsize_t extent[1] = {start[0] + count[0]};
        requireStatus(H5Dset_extent(records_, extent),
          "Cannot extend particle-detector records");
        hid_t fileSpace = H5Dget_space(records_);
        requireHandle(fileSpace,
          "Cannot select particle-detector file space");
        requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
          start, NULL, count, NULL),
          "Cannot select particle-detector record batch");
        hid_t memorySpace = H5Screate_simple(1, count, NULL);
        requireHandle(memorySpace,
          "Cannot create particle-detector memory space");
        const herr_t status = H5Dwrite(records_, memoryType_, memorySpace,
          fileSpace, H5P_DEFAULT, &buffer_[0]);
        H5Sclose(memorySpace);
        H5Sclose(fileSpace);
        requireStatus(status, "Cannot write particle-detector record batch");
        committed_ += static_cast<std::uint64_t>(buffer_.size());
        writeUnsignedScalar(committedDataset_, committed_,
          "Cannot update particle-detector commit marker");
        buffer_.clear();
      }

      void close(bool completed)
      {
        if (!open_)
          {
            closeHandles();
            return;
          }
        flush();
        writeByteScalar(completeDataset_, completed ? 1 : 0,
          "Cannot write particle-detector completion marker");
        requireStatus(H5Fflush(file_, H5F_SCOPE_GLOBAL),
          "Cannot finalize particle-detector file");
        open_ = false;
        closeHandles();
        buffer_.clear();
      }

    private:
      void closeHandles() noexcept
      {
        if (records_ >= 0) H5Dclose(records_);
        if (committedDataset_ >= 0) H5Dclose(committedDataset_);
        if (completeDataset_ >= 0) H5Dclose(completeDataset_);
        if (memoryType_ >= 0) H5Tclose(memoryType_);
        if (fileType_ >= 0) H5Tclose(fileType_);
        if (group_ >= 0) H5Gclose(group_);
        if (file_ >= 0) H5Fclose(file_);
        records_ = committedDataset_ = completeDataset_ = -1;
        memoryType_ = fileType_ = group_ = file_ = -1;
      }

      bool open_;
      std::size_t bufferLimit_;
      std::uint64_t committed_;
      std::vector<ParticlePlaneRecord> buffer_;
      hid_t file_;
      hid_t group_;
      hid_t records_;
      hid_t committedDataset_;
      hid_t completeDataset_;
      hid_t memoryType_;
      hid_t fileType_;
    };

    struct FieldPlaneRuntime
    {
      bool sampled;
      double nextTimeLab;

      FieldPlaneRuntime()
        : sampled(false), nextTimeLab(0.0)
      {}
    };
  }

  class LabDetectorManager::Impl
  {
  public:
    Impl(const DetectorConfig& config,
         const EBGridGeometry& globalGeometry,
         const FieldVector<Double>& globalOriginBox,
         const FieldVector<Double>& localOriginBox,
         const BoostFrameTransform& frame,
         MPI_Comm communicator)
      : config_(config), globalGeometry_(globalGeometry),
        globalOriginBox_(globalOriginBox),
        localOriginBox_(localOriginBox), frame_(frame),
        communicator_(communicator), rank_(0), size_(1), closed_(false),
        fieldRuntime_(config.fieldPlanes.size()), localParticleEvents_(),
        fieldWriters_(), particleWriters_()
    {
      MPI_Comm_rank(communicator_, &rank_);
      MPI_Comm_size(communicator_, &size_);
      if (!config_.enabled())
        throw std::invalid_argument(
          "LabDetectorManager requires at least one detector plane");

      createDirectories(config_.directory, communicator_);
      int failed = 0;
      std::string failureMessage;
      if (rank_ == 0)
        {
          try
            {
              fieldWriters_.reserve(config_.fieldPlanes.size());
              for (std::size_t i = 0; i < config_.fieldPlanes.size(); ++i)
                {
                  std::unique_ptr<FieldPlaneWriter> writer(
                    new FieldPlaneWriter());
                  writer->open(joinPath(config_.directory,
                    config_.fieldPlanes[i].name + ".h5"),
                    config_.fieldPlanes[i], globalGeometry_, globalOriginBox_);
                  fieldWriters_.push_back(std::move(writer));
                }
              particleWriters_.reserve(config_.particlePlanes.size());
              for (std::size_t i = 0; i < config_.particlePlanes.size(); ++i)
                {
                  std::unique_ptr<ParticlePlaneWriter> writer(
                    new ParticlePlaneWriter());
                  writer->open(joinPath(config_.directory,
                    config_.particlePlanes[i].name + ".h5"),
                    config_.particlePlanes[i], size_);
                  particleWriters_.push_back(std::move(writer));
                }
            }
          catch (const std::exception& error)
            {
              failed = 1;
              failureMessage = error.what();
            }
        }
      MPI_Bcast(&failed, 1, MPI_INT, 0, communicator_);
      if (failed)
        {
          if (rank_ == 0)
            throw std::runtime_error(failureMessage);
          throw std::runtime_error(
            "Rank zero could not initialize laboratory detector output");
        }
    }

    ~Impl()
    {
      try { close(false); }
      catch (...) {}
    }

    void captureParticleStep(const RelativisticParticleSI& before,
                             const RelativisticParticleSI& after,
                             Double timeBoxBefore,
                             Double timeBoxAfter)
    {
      if (config_.particlePlanes.empty()) return;
      Double beforeTimeLab = 0.0;
      Double beforeZLab = 0.0;
      Double afterTimeLab = 0.0;
      Double afterZLab = 0.0;
      frame_.boxToLab(timeBoxBefore, before.position[2],
                      beforeTimeLab, beforeZLab);
      frame_.boxToLab(timeBoxAfter, after.position[2],
                      afterTimeLab, afterZLab);

      for (std::size_t detector = 0;
           detector < config_.particlePlanes.size(); ++detector)
        {
          const double planeZ = config_.particlePlanes[detector].z;
          const double beforeDelta = beforeZLab - planeZ;
          const double afterDelta = afterZLab - planeZ;
          std::int32_t direction = 0;
          if (beforeDelta < 0.0 && afterDelta >= 0.0) direction = 1;
          else if (beforeDelta > 0.0 && afterDelta <= 0.0) direction = -1;
          else continue;

          const double denominator = afterZLab - beforeZLab;
          if (denominator == 0.0) continue;
          const double fraction = std::max(0.0, std::min(1.0,
            (planeZ - beforeZLab) / denominator));

          ParticleCrossingWire event = {};
          event.detectorIndex = static_cast<std::uint32_t>(detector);
          event.record.particleId = after.id;
          event.record.sourceId = after.sourceId;
          event.record.time = beforeTimeLab +
            fraction * (afterTimeLab - beforeTimeLab);
          event.record.position[0] = before.position[0] +
            fraction * (after.position[0] - before.position[0]);
          event.record.position[1] = before.position[1] +
            fraction * (after.position[1] - before.position[1]);
          event.record.position[2] = planeZ;

          FieldVector<Double> properVelocityBox(0.0);
          for (unsigned int component = 0; component < 3; ++component)
            properVelocityBox[component] = before.properVelocity[component] +
              fraction * (after.properVelocity[component] -
                          before.properVelocity[component]);
          const double gammaBox =
            BoostFrameTransform::gammaFromProperVelocity(properVelocityBox);
          event.record.properVelocity[0] = properVelocityBox[0];
          event.record.properVelocity[1] = properVelocityBox[1];
          event.record.properVelocity[2] = frame_.gamma() *
            (properVelocityBox[2] + frame_.beta() * gammaBox);
          event.record.charge = after.charge;
          event.record.mass = after.mass;
          event.record.weight = after.weight;
          event.record.direction = direction;
          localParticleEvents_.push_back(event);
        }
    }

    void collectParticleCrossings()
    {
      if (config_.particlePlanes.empty()) return;
      const unsigned long long localCount =
        static_cast<unsigned long long>(localParticleEvents_.size());
      unsigned long long globalCount = 0;
      MPI_Allreduce(&localCount, &globalCount, 1, MPI_UNSIGNED_LONG_LONG,
                    MPI_SUM, communicator_);
      if (globalCount == 0) return;

      std::vector<unsigned long long> counts;
      if (rank_ == 0) counts.resize(static_cast<std::size_t>(size_));
      MPI_Gather(&localCount, 1, MPI_UNSIGNED_LONG_LONG,
        rank_ == 0 ? &counts[0] : NULL, 1, MPI_UNSIGNED_LONG_LONG,
        0, communicator_);

      const std::size_t maximumChunk = std::max<std::size_t>(1,
        std::min<std::size_t>(1048576,
          static_cast<std::size_t>(INT_MAX) /
          sizeof(ParticleCrossingWire)));
      if (rank_ == 0)
        {
          dispatch(localParticleEvents_);
          std::vector<ParticleCrossingWire> receive;
          for (int source = 1; source < size_; ++source)
            {
              unsigned long long remaining = counts[source];
              while (remaining > 0)
                {
                  const std::size_t chunk = static_cast<std::size_t>(
                    std::min<unsigned long long>(remaining, maximumChunk));
                  receive.resize(chunk);
                  MPI_Recv(&receive[0], checkedByteCount(chunk,
                    sizeof(ParticleCrossingWire), "Particle detector receive"),
                    MPI_BYTE, source, PARTICLE_EVENT_TAG,
                    communicator_, MPI_STATUS_IGNORE);
                  dispatch(receive);
                  remaining -= static_cast<unsigned long long>(chunk);
                }
            }
        }
      else
        {
          std::size_t offset = 0;
          while (offset < localParticleEvents_.size())
            {
              const std::size_t chunk = std::min(maximumChunk,
                localParticleEvents_.size() - offset);
              MPI_Send(&localParticleEvents_[offset], checkedByteCount(chunk,
                sizeof(ParticleCrossingWire), "Particle detector send"),
                MPI_BYTE, 0, PARTICLE_EVENT_TAG, communicator_);
              offset += chunk;
            }
        }
      localParticleEvents_.clear();
    }

    void sampleFieldPlanes(const EBFieldGrid& localFields, Double timeBox)
    {
      for (std::size_t detector = 0;
           detector < config_.fieldPlanes.size(); ++detector)
        {
          const FieldDetectorPlaneConfig& plane =
            config_.fieldPlanes[detector];
          const double boxZ = frame_.boxZFromLabZAndBoxT(plane.z, timeBox);
          const int owner = ownerRank(boxZ);
          if (owner < 0) continue;
          const double timeLab =
            frame_.labTFromLabZBoxT(plane.z, timeBox);
          FieldPlaneRuntime& runtime = fieldRuntime_[detector];
          const double tolerance = 64.0 *
            std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::abs(timeLab));
          if (runtime.sampled &&
              timeLab < runtime.nextTimeLab - tolerance)
            continue;

          std::vector<FieldPointRecord> sample;
          if (rank_ == owner)
            samplePlane(localFields, boxZ, sample);
          const std::size_t points = globalGeometry_.nx * globalGeometry_.ny;
          if (owner != 0)
            {
              const int bytes = checkedByteCount(points,
                sizeof(FieldPointRecord), "Field detector plane");
              if (rank_ == owner)
                MPI_Send(&sample[0], bytes, MPI_BYTE, 0,
                         FIELD_SAMPLE_TAG, communicator_);
              else if (rank_ == 0)
                {
                  sample.resize(points);
                  MPI_Recv(&sample[0], bytes, MPI_BYTE, owner,
                    FIELD_SAMPLE_TAG, communicator_, MPI_STATUS_IGNORE);
                }
            }
          if (rank_ == 0)
            fieldWriters_[detector]->append(timeLab, sample);

          if (!runtime.sampled)
            {
              runtime.sampled = true;
              runtime.nextTimeLab = timeLab + plane.rhythm;
            }
          else
            {
              do { runtime.nextTimeLab += plane.rhythm; }
              while (runtime.nextTimeLab <= timeLab + tolerance);
            }
        }
    }

    void close(bool completed)
    {
      if (closed_) return;
      if (rank_ == 0)
        {
          for (std::size_t i = 0; i < fieldWriters_.size(); ++i)
            fieldWriters_[i]->close(completed);
          for (std::size_t i = 0; i < particleWriters_.size(); ++i)
            particleWriters_[i]->close(completed);
        }
      closed_ = true;
    }

  private:
    int ownerRank(double boxZ) const
    {
      const double coordinate =
        (boxZ - globalOriginBox_[2]) / globalGeometry_.dz;
      const double tolerance = 64.0 *
        std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(coordinate));
      if (coordinate < -tolerance ||
          coordinate > static_cast<double>(globalGeometry_.nz) + tolerance)
        return -1;
      const std::size_t cell =
        coordinate >= static_cast<double>(globalGeometry_.nz) ?
        globalGeometry_.nz - 1 :
        static_cast<std::size_t>(std::max(0.0, std::floor(coordinate)));
      const std::size_t base = globalGeometry_.nz /
        static_cast<std::size_t>(size_);
      const std::size_t remainder = globalGeometry_.nz %
        static_cast<std::size_t>(size_);
      const std::size_t longCells = (base + 1) * remainder;
      if (cell < longCells)
        return static_cast<int>(cell / (base + 1));
      return static_cast<int>(remainder + (cell - longCells) / base);
    }

    void samplePlane(const EBFieldGrid& localFields, double boxZ,
                     std::vector<FieldPointRecord>& sample) const
    {
      sample.resize(globalGeometry_.nx * globalGeometry_.ny);
      FieldVector<Double> position(0.0);
      position[2] = boxZ;
      for (std::size_t j = 0; j < globalGeometry_.ny; ++j)
        {
          position[1] = globalOriginBox_[1] +
            (static_cast<double>(j) + 0.5) * globalGeometry_.dy;
          for (std::size_t i = 0; i < globalGeometry_.nx; ++i)
            {
              position[0] = globalOriginBox_[0] +
                (static_cast<double>(i) + 0.5) * globalGeometry_.dx;
              const RadiationFieldSample fields =
                localFields.radiationSamplePositionLab(
                  position, localOriginBox_, frame_);
              FieldPointRecord& point =
                sample[j * globalGeometry_.nx + i];
              for (unsigned int component = 0; component < 3; ++component)
                {
                  point.electric[component] = fields.electric[component];
                  point.magnetic[component] = fields.magnetic[component];
                }
            }
        }
    }

    void dispatch(const std::vector<ParticleCrossingWire>& events)
    {
      for (std::size_t i = 0; i < events.size(); ++i)
        {
          const std::size_t detector = events[i].detectorIndex;
          if (detector >= particleWriters_.size())
            throw std::runtime_error(
              "Received an invalid particle-detector index");
          particleWriters_[detector]->append(events[i].record);
        }
    }

    DetectorConfig config_;
    EBGridGeometry globalGeometry_;
    FieldVector<Double> globalOriginBox_;
    FieldVector<Double> localOriginBox_;
    BoostFrameTransform frame_;
    MPI_Comm communicator_;
    int rank_;
    int size_;
    bool closed_;
    std::vector<FieldPlaneRuntime> fieldRuntime_;
    std::vector<ParticleCrossingWire> localParticleEvents_;
    std::vector<std::unique_ptr<FieldPlaneWriter> > fieldWriters_;
    std::vector<std::unique_ptr<ParticlePlaneWriter> > particleWriters_;
  };

  LabDetectorManager::LabDetectorManager(
      const DetectorConfig& config,
      const EBGridGeometry& globalGeometry,
      const FieldVector<Double>& globalOriginBox,
      const FieldVector<Double>& localOriginBox,
      const BoostFrameTransform& frame,
      MPI_Comm communicator)
    : impl_(new Impl(config, globalGeometry, globalOriginBox,
                     localOriginBox, frame, communicator))
  {}

  LabDetectorManager::~LabDetectorManager() {}

  void LabDetectorManager::captureParticleStep(
      const RelativisticParticleSI& before,
      const RelativisticParticleSI& after,
      Double timeBoxBefore,
      Double timeBoxAfter)
  {
    impl_->captureParticleStep(before, after, timeBoxBefore, timeBoxAfter);
  }

  void LabDetectorManager::collectParticleCrossings()
  {
    impl_->collectParticleCrossings();
  }

  void LabDetectorManager::sampleFieldPlanes(
      const EBFieldGrid& localFields, Double timeBox)
  {
    impl_->sampleFieldPlanes(localFields, timeBox);
  }

  void LabDetectorManager::close(bool completed)
  {
    impl_->close(completed);
  }
}
