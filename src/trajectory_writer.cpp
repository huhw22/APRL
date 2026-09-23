#include "trajectory_writer.h"

#include <cmath>
#include <stdexcept>

namespace fel
{
  namespace
  {
    void requireHandle(hid_t handle, const std::string& message)
    {
      if (handle < 0) throw std::runtime_error(message);
    }

    void requireStatus(herr_t status, const std::string& message)
    {
      if (status < 0) throw std::runtime_error(message);
    }
  }

  TrajectoryRecord::TrajectoryRecord()
    : particleId(0), sourceId(0), time(0.0),
      charge(0.0), weight(1.0),
      event(static_cast<std::uint8_t>(TrajectoryEvent::Sample)),
      boundaryFace(-1)
  {
    position[0] = position[1] = position[2] = 0.0;
    properVelocity[0] = properVelocity[1] = properVelocity[2] = 0.0;
  }

  TrajectoryWriter::TrajectoryWriter()
    : isOpen_(false), swmrEnabled_(false), durableProgress_(true),
      bufferLimit_(0),
      committedRecords_(0), buffer_(), file_(-1), group_(-1),
      recordDataset_(-1), committedDataset_(-1), completeDataset_(-1),
      memoryRecordType_(-1), fileRecordType_(-1)
  {}

  TrajectoryWriter::~TrajectoryWriter()
  {
    try
      {
        close(false);
      }
    catch (...)
      {
        closeHandles();
      }
  }

  hid_t TrajectoryWriter::createMemoryRecordType() const
  {
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(TrajectoryRecord));
    requireHandle(type, "Cannot create trajectory memory datatype");
    try
      {
        requireStatus(H5Tinsert(type, "particle_id",
                                HOFFSET(TrajectoryRecord, particleId),
                                H5T_NATIVE_UINT64),
                      "Cannot add trajectory particle_id field");
        requireStatus(H5Tinsert(type, "source_id",
                                HOFFSET(TrajectoryRecord, sourceId),
                                H5T_NATIVE_UINT64),
                      "Cannot add trajectory source_id field");
        requireStatus(H5Tinsert(type, "time_s",
                                HOFFSET(TrajectoryRecord, time),
                                H5T_NATIVE_DOUBLE),
                      "Cannot add trajectory time field");

        hsize_t vectorSize[1] = {3};
        hid_t vectorType = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, vectorSize);
        requireHandle(vectorType, "Cannot create trajectory vector datatype");
        requireStatus(H5Tinsert(type, "position_m",
                                HOFFSET(TrajectoryRecord, position),
                                vectorType),
                      "Cannot add trajectory position field");
        requireStatus(H5Tinsert(type, "proper_velocity",
                                HOFFSET(TrajectoryRecord, properVelocity),
                                vectorType),
                      "Cannot add trajectory momentum field");
        H5Tclose(vectorType);

        requireStatus(H5Tinsert(type, "charge_C",
                                HOFFSET(TrajectoryRecord, charge),
                                H5T_NATIVE_DOUBLE),
                      "Cannot add trajectory charge field");
        requireStatus(H5Tinsert(type, "weight",
                                HOFFSET(TrajectoryRecord, weight),
                                H5T_NATIVE_DOUBLE),
                      "Cannot add trajectory weight field");
        requireStatus(H5Tinsert(type, "event_type",
                                HOFFSET(TrajectoryRecord, event),
                                H5T_NATIVE_UCHAR),
                      "Cannot add trajectory event field");
        requireStatus(H5Tinsert(type, "boundary_face",
                                HOFFSET(TrajectoryRecord, boundaryFace),
                                H5T_NATIVE_SCHAR),
                      "Cannot add trajectory boundary face field");
      }
    catch (...)
      {
        H5Tclose(type);
        throw;
      }
    return type;
  }

  hid_t TrajectoryWriter::createFileRecordType() const
  {
    const std::size_t u64 = 8;
    const std::size_t f64 = 8;
    const std::size_t fileSize = 2 * u64 + 9 * f64 + 2;
    hid_t type = H5Tcreate(H5T_COMPOUND, fileSize);
    requireHandle(type, "Cannot create trajectory file datatype");
    try
      {
        std::size_t offset = 0;
        requireStatus(H5Tinsert(type, "particle_id", offset, H5T_STD_U64LE),
                      "Cannot add trajectory file particle_id");
        offset += u64;
        requireStatus(H5Tinsert(type, "source_id", offset, H5T_STD_U64LE),
                      "Cannot add trajectory file source_id");
        offset += u64;
        requireStatus(H5Tinsert(type, "time_s", offset, H5T_IEEE_F64LE),
                      "Cannot add trajectory file time");
        offset += f64;

        hsize_t vectorSize[1] = {3};
        hid_t vectorType = H5Tarray_create2(H5T_IEEE_F64LE, 1, vectorSize);
        requireHandle(vectorType, "Cannot create trajectory file vector type");
        requireStatus(H5Tinsert(type, "position_m", offset, vectorType),
                      "Cannot add trajectory file position");
        offset += 3 * f64;
        requireStatus(H5Tinsert(type, "proper_velocity", offset, vectorType),
                      "Cannot add trajectory file momentum");
        offset += 3 * f64;
        H5Tclose(vectorType);

        requireStatus(H5Tinsert(type, "charge_C", offset, H5T_IEEE_F64LE),
                      "Cannot add trajectory file charge");
        offset += f64;
        requireStatus(H5Tinsert(type, "weight", offset, H5T_IEEE_F64LE),
                      "Cannot add trajectory file weight");
        offset += f64;
        requireStatus(H5Tinsert(type, "event_type", offset, H5T_STD_U8LE),
                      "Cannot add trajectory file event");
        offset += 1;
        requireStatus(H5Tinsert(type, "boundary_face", offset, H5T_STD_I8LE),
                      "Cannot add trajectory file boundary face");
      }
    catch (...)
      {
        H5Tclose(type);
        throw;
      }
    return type;
  }

  void TrajectoryWriter::open(const std::string& filename,
                              int mpiRank, int mpiSize,
                              std::size_t bufferRecords,
                              unsigned int compressionLevel,
                              bool durableProgress)
  {
    if (isOpen_) throw std::runtime_error("Trajectory file is already open");
    if (filename.empty())
      throw std::invalid_argument("Trajectory filename cannot be empty");
    if (mpiRank < 0 || mpiSize < 1 || mpiRank >= mpiSize)
      throw std::invalid_argument("Invalid trajectory MPI metadata");
    if (bufferRecords == 0)
      throw std::invalid_argument("Trajectory buffer must contain records");
    if (compressionLevel > 9)
      throw std::invalid_argument("Trajectory compression must be in [0,9]");

    bufferLimit_ = bufferRecords;
    durableProgress_ = durableProgress;
    committedRecords_ = 0;
    buffer_.clear();
    buffer_.reserve(bufferLimit_);

    hid_t access = H5Pcreate(H5P_FILE_ACCESS);
    requireHandle(access, "Cannot create trajectory file-access properties");
#if H5_VERSION_GE(1, 10, 0)
    requireStatus(H5Pset_libver_bounds(access, H5F_LIBVER_LATEST,
                                       H5F_LIBVER_LATEST),
                  "Cannot select crash-readable HDF5 format");
#endif
    file_ = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,
                      H5P_DEFAULT, access);
    H5Pclose(access);
    requireHandle(file_, "Cannot create trajectory file: " + filename);

    try
      {
        group_ = H5Gcreate2(file_, "/trajectory", H5P_DEFAULT,
                            H5P_DEFAULT, H5P_DEFAULT);
        requireHandle(group_, "Cannot create trajectory group");
        memoryRecordType_ = createMemoryRecordType();
        fileRecordType_ = createFileRecordType();

        hsize_t initial[1] = {0};
        hsize_t maximum[1] = {H5S_UNLIMITED};
        hid_t space = H5Screate_simple(1, initial, maximum);
        requireHandle(space, "Cannot create trajectory record space");
        hid_t creation = H5Pcreate(H5P_DATASET_CREATE);
        if (creation < 0)
          {
            H5Sclose(space);
            throw std::runtime_error("Cannot create trajectory chunk properties");
          }
        hsize_t chunk[1] = {static_cast<hsize_t>(bufferLimit_)};
        requireStatus(H5Pset_chunk(creation, 1, chunk),
                      "Cannot set trajectory chunk size");
        requireStatus(H5Pset_shuffle(creation),
                      "Cannot enable trajectory byte shuffle");
        if (compressionLevel > 0)
          requireStatus(H5Pset_deflate(creation, compressionLevel),
                        "Cannot enable trajectory compression");
        recordDataset_ = H5Dcreate2(group_, "records", fileRecordType_,
                                    space, H5P_DEFAULT, creation,
                                    H5P_DEFAULT);
        H5Pclose(creation);
        H5Sclose(space);
        requireHandle(recordDataset_, "Cannot create trajectory record dataset");

        hid_t scalar = H5Screate(H5S_SCALAR);
        requireHandle(scalar, "Cannot create trajectory scalar space");
        committedDataset_ = H5Dcreate2(
          group_, "committed_records", H5T_STD_U64LE, scalar,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        completeDataset_ = H5Dcreate2(
          group_, "complete", H5T_STD_U8LE, scalar,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(scalar);
        requireHandle(committedDataset_,
                      "Cannot create trajectory commit marker");
        requireHandle(completeDataset_,
                      "Cannot create trajectory completion marker");
        writeUnsignedScalar(committedDataset_, 0);
        writeByteScalar(completeDataset_, 0);

        writeIntAttribute(group_, "format_version", 2);
        writeIntAttribute(group_, "mpi_rank", mpiRank);
        writeIntAttribute(group_, "mpi_size", mpiSize);
        writeIntAttribute(group_, "compression_level",
                          static_cast<int>(compressionLevel));
        writeStringAttribute(group_, "write_mode",
          durableProgress_ ? "interactive" : "throughput");
        writeStringAttribute(group_, "frame", "laboratory");
        writeStringAttribute(group_, "position_unit", "m");
        writeStringAttribute(group_, "time_unit", "s");
        writeStringAttribute(group_, "charge_unit", "C");
        writeStringAttribute(group_, "proper_velocity_unit", "gamma*v/c");
        writeStringAttribute(group_, "event_type_definition",
          "0=periodic_sample,1=cpml_entry,2=outer_domain_exit,3=particle_retirement_entry");
        writeStringAttribute(group_, "boundary_face_definition",
          "-1=none,0=x-,1=x+,2=y-,3=y+,4=z-,5=z+");
        writeStringAttribute(group_, "reader_contract",
          "read only records[0:committed_records]");

        if (durableProgress_)
          requireStatus(H5Fflush(file_, H5F_SCOPE_GLOBAL),
                        "Cannot initialize trajectory file on disk");
#if H5_VERSION_GE(1, 10, 0)
        swmrEnabled_ = durableProgress_ &&
          H5Fstart_swmr_write(file_) >= 0;
#else
        swmrEnabled_ = false;
#endif
        isOpen_ = true;
      }
    catch (...)
      {
        closeHandles();
        buffer_.clear();
        throw;
      }
  }

  void TrajectoryWriter::append(const TrajectoryRecord& record)
  {
    if (!isOpen_)
      throw std::runtime_error("Trajectory append requires an open file");
    if (!std::isfinite(record.time) || !std::isfinite(record.charge) ||
        !std::isfinite(record.weight))
      throw std::invalid_argument("Trajectory scalar is not finite");
    const std::uint8_t retirement =
      static_cast<std::uint8_t>(TrajectoryEvent::RetirementEntry);
    if (record.event > retirement ||
        record.boundaryFace < -1 || record.boundaryFace > 5 ||
        (record.event ==
           static_cast<std::uint8_t>(TrajectoryEvent::Sample) &&
         record.boundaryFace != -1) ||
        (record.event != retirement && record.event !=
           static_cast<std::uint8_t>(TrajectoryEvent::Sample) &&
         record.boundaryFace < 0))
      throw std::invalid_argument("Invalid trajectory event metadata");
    for (unsigned int component = 0; component < 3; ++component)
      if (!std::isfinite(record.position[component]) ||
          !std::isfinite(record.properVelocity[component]))
        throw std::invalid_argument("Trajectory vector is not finite");

    buffer_.push_back(record);
    if (buffer_.size() >= bufferLimit_) flush();
  }

  void TrajectoryWriter::flush()
  {
    if (!isOpen_)
      throw std::runtime_error("Trajectory flush requires an open file");
    if (buffer_.empty()) return;

    const hsize_t start[1] = {static_cast<hsize_t>(committedRecords_)};
    const hsize_t count[1] = {static_cast<hsize_t>(buffer_.size())};
    const hsize_t extent[1] = {start[0] + count[0]};
    requireStatus(H5Dset_extent(recordDataset_, extent),
                  "Cannot extend trajectory record dataset");
    hid_t fileSpace = H5Dget_space(recordDataset_);
    requireHandle(fileSpace, "Cannot select trajectory file space");
    try
      {
        requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
                                          start, NULL, count, NULL),
                      "Cannot select trajectory record batch");
        hid_t memorySpace = H5Screate_simple(1, count, NULL);
        requireHandle(memorySpace, "Cannot create trajectory memory space");
        const herr_t status = H5Dwrite(
          recordDataset_, memoryRecordType_, memorySpace, fileSpace,
          H5P_DEFAULT, &buffer_[0]);
        H5Sclose(memorySpace);
        requireStatus(status, "Cannot write trajectory record batch");
      }
    catch (...)
      {
        H5Sclose(fileSpace);
        throw;
      }
    H5Sclose(fileSpace);

    if (durableProgress_)
      {
#if H5_VERSION_GE(1, 10, 0)
        requireStatus(H5Dflush(recordDataset_),
                      "Cannot flush trajectory record data");
#endif
        requireStatus(H5Fflush(file_, H5F_SCOPE_GLOBAL),
                      "Cannot make trajectory record data durable");
      }

    const std::uint64_t newCommitted =
      committedRecords_ + static_cast<std::uint64_t>(buffer_.size());
    writeUnsignedScalar(committedDataset_, newCommitted);
    if (durableProgress_)
      {
#if H5_VERSION_GE(1, 10, 0)
        requireStatus(H5Dflush(committedDataset_),
                      "Cannot flush trajectory commit marker");
#endif
        requireStatus(H5Fflush(file_, H5F_SCOPE_GLOBAL),
                      "Cannot make trajectory commit marker durable");
      }
    committedRecords_ = newCommitted;
    buffer_.clear();
  }

  void TrajectoryWriter::close(bool completed)
  {
    if (!isOpen_)
      {
        closeHandles();
        return;
      }
    flush();
    writeByteScalar(completeDataset_, completed ? 1 : 0);
#if H5_VERSION_GE(1, 10, 0)
    if (durableProgress_)
      requireStatus(H5Dflush(completeDataset_),
                    "Cannot flush trajectory completion marker");
#endif
    requireStatus(H5Fflush(file_, H5F_SCOPE_GLOBAL),
                  "Cannot finalize trajectory file");
    isOpen_ = false;
    closeHandles();
    buffer_.clear();
  }

  void TrajectoryWriter::writeUnsignedScalar(
      hid_t dataset, std::uint64_t value) const
  {
    requireStatus(H5Dwrite(dataset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, &value),
                  "Cannot write trajectory unsigned scalar");
  }

  void TrajectoryWriter::writeByteScalar(
      hid_t dataset, unsigned char value) const
  {
    requireStatus(H5Dwrite(dataset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, &value),
                  "Cannot write trajectory byte scalar");
  }

  void TrajectoryWriter::writeIntAttribute(
      hid_t object, const char* name, int value) const
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, std::string("Cannot create attribute space: ") + name);
    hid_t attribute = H5Acreate2(object, name, H5T_STD_I32LE, space,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute: ") + name);
    const herr_t status = H5Awrite(attribute, H5T_NATIVE_INT, &value);
    H5Aclose(attribute);
    requireStatus(status, std::string("Cannot write attribute: ") + name);
  }

  void TrajectoryWriter::writeStringAttribute(
      hid_t object, const char* name, const std::string& value) const
  {
    hid_t type = H5Tcopy(H5T_C_S1);
    requireHandle(type, std::string("Cannot create string type: ") + name);
    requireStatus(H5Tset_size(type, value.size() + 1),
                  std::string("Cannot size string type: ") + name);
    requireStatus(H5Tset_strpad(type, H5T_STR_NULLTERM),
                  std::string("Cannot set string padding: ") + name);
    hid_t space = H5Screate(H5S_SCALAR);
    if (space < 0)
      {
        H5Tclose(type);
        throw std::runtime_error(std::string("Cannot create attribute space: ") + name);
      }
    hid_t attribute = H5Acreate2(object, name, type, space,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    if (attribute < 0)
      {
        H5Tclose(type);
        throw std::runtime_error(std::string("Cannot create attribute: ") + name);
      }
    const herr_t status = H5Awrite(attribute, type, value.c_str());
    H5Aclose(attribute);
    H5Tclose(type);
    requireStatus(status, std::string("Cannot write attribute: ") + name);
  }

  void TrajectoryWriter::closeHandles() noexcept
  {
    if (recordDataset_ >= 0) H5Dclose(recordDataset_);
    if (committedDataset_ >= 0) H5Dclose(committedDataset_);
    if (completeDataset_ >= 0) H5Dclose(completeDataset_);
    if (memoryRecordType_ >= 0) H5Tclose(memoryRecordType_);
    if (fileRecordType_ >= 0) H5Tclose(fileRecordType_);
    if (group_ >= 0) H5Gclose(group_);
    if (file_ >= 0) H5Fclose(file_);
    recordDataset_ = committedDataset_ = completeDataset_ = -1;
    memoryRecordType_ = fileRecordType_ = -1;
    group_ = file_ = -1;
    isOpen_ = false;
    swmrEnabled_ = false;
  }

  bool TrajectoryWriter::isOpen() const { return isOpen_; }
  bool TrajectoryWriter::swmrEnabled() const { return swmrEnabled_; }
  std::size_t TrajectoryWriter::bufferedRecords() const
  {
    return buffer_.size();
  }
  std::uint64_t TrajectoryWriter::committedRecords() const
  {
    return committedRecords_;
  }
  std::size_t TrajectoryWriter::memoryBytes() const
  {
    return buffer_.capacity() * sizeof(TrajectoryRecord);
  }
}
