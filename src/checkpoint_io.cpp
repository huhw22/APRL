#include "checkpoint_io.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

#include <hdf5.h>

namespace aprl
{
  namespace
  {
    enum BlockKind { ByteBlock = 1, DoubleBlock = 2,
                     ParticleBlock = 3, CarrierBlock = 4 };
    const std::size_t headerCount = 12;
    const unsigned long long checkpointMagic = 0x4150524c434b5031ULL;
    const unsigned long long checkpointVersion = 1;
    const std::size_t transferChunkBytes = 64U * 1024U * 1024U;
    const std::size_t particleFileBytes = 10 * 8 + 2 * 8 + 1;
    const std::size_t carrierFileBytes = 8 * 8;

    class Handle
    {
    public:
      typedef herr_t (*Closer)(hid_t);
      Handle(hid_t value, Closer closer) : value_(value), closer_(closer) {}
      ~Handle() { if (value_ >= 0) closer_(value_); }
      Handle(const Handle&) = delete;
      Handle& operator=(const Handle&) = delete;
      operator hid_t() const { return value_; }
    private:
      hid_t value_;
      Closer closer_;
    };

    bool validName(const std::string& name)
    {
      if (name.empty() || name.front() == '/' || name.back() == '/')
        return false;
      std::size_t start = 0;
      while (start < name.size())
        {
          const std::size_t end = name.find('/', start);
          const std::string part = name.substr(start,
            end == std::string::npos ? std::string::npos : end - start);
          if (part.empty() || part == "." || part == "..") return false;
          for (std::size_t i = 0; i < part.size(); ++i)
            {
              const unsigned char c = static_cast<unsigned char>(part[i]);
              if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                    (c >= '0' && c <= '9') || c == '_' || c == '-' ||
                    c == '.')) return false;
            }
          if (end == std::string::npos) break;
          start = end + 1;
        }
      return true;
    }

    template<typename Record, typename Field>
    std::size_t fieldOffset(const Record& record, const Field& field)
    {
      return static_cast<std::size_t>(
        reinterpret_cast<const unsigned char*>(&field) -
        reinterpret_cast<const unsigned char*>(&record));
    }

    std::size_t fileWidth(unsigned long long kind)
    {
      switch (kind)
        {
        case ByteBlock: return 1;
        case DoubleBlock: return 8;
        case ParticleBlock: return particleFileBytes;
        case CarrierBlock: return carrierFileBytes;
        default: return 0;
        }
    }
  }

  class CheckpointArchive::Impl
  {
  public:
    struct Block
    {
      unsigned long long kind;
      unsigned long long count;
      unsigned long long offset;
      unsigned long long total;
      unsigned long long maximum;
      Block() : kind(0), count(0), offset(0), total(0), maximum(0) {}
    };

    Impl(MPI_Comm communicator, const std::string& path, bool writing)
      : communicator_(communicator), path_(path), partial_(path + ".partial"),
        writing_(writing), closed_(false), rank_(0), size_(0), file_(-1),
        transferred_(0)
    {
      if (communicator_ == MPI_COMM_NULL)
        throw std::invalid_argument("Checkpoint communicator is null");
      MPI_Comm_rank(communicator_, &rank_);
      MPI_Comm_size(communicator_, &size_);
      sameString(path_, "archive path");
      unsigned long long mode = writing ? 1 : 0, rootMode = mode;
      MPI_Bcast(&rootMode, 1, MPI_UNSIGNED_LONG_LONG, 0, communicator_);
      require(mode == rootMode, "Checkpoint access mode differs between ranks");
      require(!path_.empty() && path_.find('\0') == std::string::npos &&
        !(path_.size() >= 8 && path_.substr(path_.size() - 8) == ".partial"),
        "Checkpoint path must be a nonempty final filename, not .partial");
#ifndef H5_HAVE_PARALLEL
      require(size_ == 1, "Multi-rank checkpoint requires parallel HDF5");
#endif
      if (writing_)
        {
          bool available = true;
          if (rank_ == 0)
            {
              struct stat status;
              errno = 0;
              available = ::lstat(path_.c_str(), &status) != 0 && errno == ENOENT;
              errno = 0;
              available = available &&
                (::lstat(partial_.c_str(), &status) != 0 && errno == ENOENT);
            }
          require(available,
            "Refusing existing or inaccessible checkpoint final/partial path: " + path_);
        }
      Handle access(H5Pcreate(H5P_FILE_ACCESS), H5Pclose);
      require(access >= 0, "Cannot create checkpoint file access properties");
#ifdef H5_HAVE_PARALLEL
      require(H5Pset_fapl_mpio(access, communicator_, MPI_INFO_NULL) >= 0,
              "Cannot enable parallel checkpoint file I/O");
#endif
      file_ = writing_ ? H5Fcreate(partial_.c_str(), H5F_ACC_EXCL,
        H5P_DEFAULT, access) : H5Fopen(path_.c_str(), H5F_ACC_RDONLY, access);
      require(file_ >= 0, "Cannot open checkpoint: " + path_);
      if (writing_)
        {
          makeHeader(false);
          const hsize_t extent[1] = {headerCount};
          Handle space(H5Screate_simple(1, extent, NULL), H5Sclose);
          require(space >= 0, "Cannot create checkpoint header space");
          Handle header(H5Dcreate2(file_, "/header", H5T_STD_U64LE, space,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
          require(header >= 0, "Cannot create checkpoint header");
          rootTransfer(header, H5T_NATIVE_ULLONG, header_, headerCount, true);
          Handle blocks(H5Gcreate2(file_, "/blocks", H5P_DEFAULT,
            H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
          require(blocks >= 0, "Cannot create checkpoint blocks group");
        }
      else
        {
          Handle header(H5Dopen2(file_, "/header", H5P_DEFAULT), H5Dclose);
          require(header >= 0, "Checkpoint header is missing");
          checkExtent(header, headerCount);
          require(H5Dread(header, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, header_) >= 0, "Cannot read checkpoint header");
          unsigned long long expected[headerCount];
          std::copy(header_, header_ + headerCount, expected);
          makeHeader(true);
          require(std::equal(header_, header_ + headerCount, expected),
            "Checkpoint schema/completion/MPI count/native ABI does not match this run");
        }
    }

    ~Impl() {}

    void require(bool success, const std::string& message) const
    {
      const int failed = success ? 0 : 1;
      int anyFailed = 0;
      MPI_Allreduce(&failed, &anyFailed, 1, MPI_INT, MPI_MAX, communicator_);
      if (anyFailed) throw std::runtime_error(message);
    }

    void sameString(const std::string& value, const char* context) const
    {
      unsigned long long length = rank_ == 0 ? value.size() : 0;
      MPI_Bcast(&length, 1, MPI_UNSIGNED_LONG_LONG, 0, communicator_);
      require(length <= static_cast<unsigned long long>(INT_MAX),
              "Checkpoint name exceeds MPI count range");
      std::string expected;
      bool allocated = true;
      try { expected.resize(static_cast<std::size_t>(length)); }
      catch (...) { allocated = false; }
      require(allocated, "Cannot allocate checkpoint name metadata");
      if (rank_ == 0) expected = value;
      if (length) MPI_Bcast(&expected[0], static_cast<int>(length), MPI_CHAR,
                            0, communicator_);
      require(value == expected,
              std::string("Checkpoint ") + context + " differs between ranks");
    }

    void begin(const std::string& name, bool writing) const
    {
      require(!closed_ && writing_ == writing,
              "Checkpoint archive is closed or has the wrong access mode");
      sameString(name, "block name");
      require(validName(name), "Invalid checkpoint block name: " + name);
    }

    void makeHeader(bool complete)
    {
      const std::uint16_t endian = 1;
      header_[0] = checkpointMagic;
      header_[1] = checkpointVersion;
      header_[2] = complete ? 1 : 0;
      header_[3] = static_cast<unsigned long long>(size_);
      header_[4] = sizeof(double);
      header_[5] = sizeof(std::size_t);
      header_[6] = sizeof(long double);
      header_[7] = std::numeric_limits<long double>::digits;
      header_[8] = std::numeric_limits<long double>::max_exponent;
      header_[9] = *reinterpret_cast<const unsigned char*>(&endian);
      header_[10] = std::numeric_limits<double>::is_iec559 ? 1 : 0;
      header_[11] = sizeof(bool);
    }

    void checkExtent(hid_t dataset, unsigned long long expected) const
    {
      Handle space(H5Dget_space(dataset), H5Sclose);
      require(space >= 0, "Cannot inspect checkpoint dataspace");
      require(H5Sget_simple_extent_ndims(space) == 1,
              "Checkpoint dataset must be one-dimensional");
      hsize_t extent[1] = {0};
      require(H5Sget_simple_extent_dims(space, extent, NULL) >= 0 &&
        extent[0] == expected, "Checkpoint dataset extent is inconsistent");
    }

    void rootTransfer(hid_t dataset, hid_t type, void* data,
                      std::size_t count, bool writing) const
    {
      Handle space(H5Dget_space(dataset), H5Sclose);
      require(space >= 0, "Cannot inspect checkpoint metadata space");
      const hsize_t extent[1] = {static_cast<hsize_t>(count)};
      Handle memory(rank_ == 0 ? H5Screate_simple(1, extent, NULL) :
        H5Screate(H5S_NULL), H5Sclose);
      require(memory >= 0, "Cannot create checkpoint metadata memory space");
      require(rank_ == 0 || H5Sselect_none(space) >= 0,
              "Cannot select empty checkpoint metadata space");
      Handle transfer(H5Pcreate(H5P_DATASET_XFER), H5Pclose);
      require(transfer >= 0, "Cannot create checkpoint metadata transfer");
#ifdef H5_HAVE_PARALLEL
      require(H5Pset_dxpl_mpio(transfer, H5FD_MPIO_COLLECTIVE) >= 0,
              "Cannot enable collective checkpoint metadata transfer");
#endif
      unsigned char empty = 0;
      void* buffer = rank_ == 0 ? data : &empty;
      const herr_t status = writing ? H5Dwrite(dataset, type, memory, space,
        transfer, buffer) : H5Dread(dataset, type, memory, space, transfer, buffer);
      require(status >= 0, "Checkpoint metadata transfer failed");
    }

    std::vector<unsigned long long> counts(std::size_t local) const
    {
      std::vector<unsigned long long> result;
      bool allocated = true;
      try { result.resize(static_cast<std::size_t>(size_)); }
      catch (...) { allocated = false; }
      require(allocated, "Cannot allocate checkpoint rank counts");
      const unsigned long long localCount = static_cast<unsigned long long>(local);
      MPI_Allgather(&localCount, 1, MPI_UNSIGNED_LONG_LONG, result.data(), 1,
                    MPI_UNSIGNED_LONG_LONG, communicator_);
      return result;
    }

    Block describe(unsigned long long kind,
                   const std::vector<unsigned long long>& rankCounts) const
    {
      Block result;
      result.kind = kind;
      result.count = rankCounts[static_cast<std::size_t>(rank_)];
      bool fits = fileWidth(kind) != 0;
      for (int rank = 0; rank < size_; ++rank)
        {
          const unsigned long long count = rankCounts[static_cast<std::size_t>(rank)];
          if (rank == rank_) result.offset = result.total;
          fits = fits && count <= std::numeric_limits<unsigned long long>::max() -
            result.total;
          if (fits) result.total += count;
          result.maximum = std::max(result.maximum, count);
        }
      fits = fits && result.total <= std::numeric_limits<hsize_t>::max() &&
        result.count <= std::numeric_limits<std::size_t>::max() / fileWidth(kind);
      require(fits, "Checkpoint block sizes overflow address or file space");
      return result;
    }

    void createParents(const std::string& path) const
    {
      for (std::size_t slash = path.find('/', 1); slash != std::string::npos;
           slash = path.find('/', slash + 1))
        {
          const std::string parent = path.substr(0, slash);
          const htri_t exists = H5Lexists(file_, parent.c_str(), H5P_DEFAULT);
          require(exists >= 0, "Cannot inspect checkpoint parent group");
          if (!exists)
            {
              Handle group(H5Gcreate2(file_, parent.c_str(), H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
              require(group >= 0, "Cannot create checkpoint parent group");
            }
        }
    }

    void write(const std::string& name, const void* data, std::size_t count,
               unsigned long long kind, hid_t memoryType, hid_t fileType,
               std::size_t memoryWidth)
    {
      begin(name, true);
      require(count == 0 || data != NULL, "Null checkpoint data pointer");
      require(count <= std::numeric_limits<std::size_t>::max() / memoryWidth,
              "Checkpoint local memory size overflow");
      std::vector<unsigned long long> rankCounts = counts(count);
      const Block block = describe(kind, rankCounts);
      const std::string path = "/blocks/" + name;
      createParents(path);
      require(H5Lexists(file_, path.c_str(), H5P_DEFAULT) == 0,
              "Checkpoint block already exists: " + name);
      Handle group(H5Gcreate2(file_, path.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                             H5P_DEFAULT), H5Gclose);
      require(group >= 0, "Cannot create checkpoint block: " + name);
      const hsize_t one[1] = {1};
      Handle kindSpace(H5Screate_simple(1, one, NULL), H5Sclose);
      require(kindSpace >= 0, "Cannot create checkpoint kind space");
      Handle kindData(H5Dcreate2(group, "kind", H5T_STD_U64LE, kindSpace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
      require(kindData >= 0, "Cannot create checkpoint kind");
      rootTransfer(kindData, H5T_NATIVE_ULLONG, &kind, 1, true);
      const hsize_t ranks[1] = {static_cast<hsize_t>(size_)};
      Handle countSpace(H5Screate_simple(1, ranks, NULL), H5Sclose);
      require(countSpace >= 0, "Cannot create checkpoint rank-count space");
      Handle countData(H5Dcreate2(group, "rank_counts", H5T_STD_U64LE,
        countSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
      require(countData >= 0, "Cannot create checkpoint rank counts");
      rootTransfer(countData, H5T_NATIVE_ULLONG, rankCounts.data(),
                   rankCounts.size(), true);
      const hsize_t extent[1] = {static_cast<hsize_t>(block.total)};
      Handle dataSpace(H5Screate_simple(1, extent, NULL), H5Sclose);
      require(dataSpace >= 0, "Cannot create checkpoint data space");
      Handle dataset(H5Dcreate2(group, "data", fileType, dataSpace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
      require(dataset >= 0, "Cannot create checkpoint dataset: " + name);
      transfer(dataset, const_cast<void*>(data), block, memoryType,
               memoryWidth, true);
    }

    Block inspect(const std::string& name)
    {
      begin(name, false);
      const std::string path = "/blocks/" + name;
      Handle group(H5Gopen2(file_, path.c_str(), H5P_DEFAULT), H5Gclose);
      require(group >= 0, "Checkpoint block is missing: " + name);
      Handle kindData(H5Dopen2(group, "kind", H5P_DEFAULT), H5Dclose);
      require(kindData >= 0, "Checkpoint block kind is missing");
      checkExtent(kindData, 1);
      unsigned long long kind = 0;
      require(H5Dread(kindData, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &kind) >= 0, "Cannot read checkpoint block kind");
      Handle countData(H5Dopen2(group, "rank_counts", H5P_DEFAULT), H5Dclose);
      require(countData >= 0, "Checkpoint rank counts are missing");
      checkExtent(countData, static_cast<unsigned long long>(size_));
      std::vector<unsigned long long> rankCounts = counts(0);
      require(H5Dread(countData, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, rankCounts.data()) >= 0, "Cannot read checkpoint rank counts");
      const Block block = describe(kind, rankCounts);
      Handle dataset(H5Dopen2(group, "data", H5P_DEFAULT), H5Dclose);
      require(dataset >= 0, "Checkpoint block data is missing");
      checkExtent(dataset, block.total);
      Handle storedType(H5Dget_type(dataset), H5Tclose);
      require(storedType >= 0 && H5Tget_size(storedType) == fileWidth(kind),
              "Checkpoint dataset type size is inconsistent");
      return block;
    }

    void read(const std::string& name, void* data, std::size_t count,
              unsigned long long kind, hid_t memoryType, hid_t expectedFileType,
              std::size_t memoryWidth)
    {
      const Block block = inspect(name);
      require(block.kind == kind && block.count == count,
              "Checkpoint block kind/count does not match: " + name);
      require(count == 0 || data != NULL, "Null checkpoint restore pointer");
      require(count <= std::numeric_limits<std::size_t>::max() / memoryWidth,
              "Checkpoint restore memory size overflow");
      const std::string path = "/blocks/" + name + "/data";
      Handle dataset(H5Dopen2(file_, path.c_str(), H5P_DEFAULT), H5Dclose);
      require(dataset >= 0, "Cannot open checkpoint block data");
      Handle storedType(H5Dget_type(dataset), H5Tclose);
      require(storedType >= 0 && H5Tequal(storedType, expectedFileType) > 0,
              "Checkpoint dataset datatype does not match its schema");
      transfer(dataset, data, block, memoryType, memoryWidth, false);
    }

    void transfer(hid_t dataset, void* data, const Block& block,
                  hid_t memoryType, std::size_t memoryWidth, bool writing)
    {
      const unsigned long long chunk = std::max<std::size_t>(1,
        transferChunkBytes / std::max(memoryWidth, fileWidth(block.kind)));
      const unsigned long long rounds = std::max<unsigned long long>(1,
        block.maximum / chunk + (block.maximum % chunk != 0));
      Handle fileSpace(H5Dget_space(dataset), H5Sclose);
      require(fileSpace >= 0, "Cannot inspect checkpoint transfer space");
      Handle transferProperties(H5Pcreate(H5P_DATASET_XFER), H5Pclose);
      require(transferProperties >= 0, "Cannot create checkpoint transfer properties");
#ifdef H5_HAVE_PARALLEL
      require(H5Pset_dxpl_mpio(transferProperties, H5FD_MPIO_COLLECTIVE) >= 0,
              "Cannot enable collective checkpoint data transfer");
#endif
      for (unsigned long long round = 0; round < rounds; ++round)
        {
          const unsigned long long localStart = std::min(block.count, round * chunk);
          const hsize_t amount[1] = {static_cast<hsize_t>(
            std::min(chunk, block.count - localStart))};
          const hsize_t start[1] = {static_cast<hsize_t>(block.offset + localStart)};
          const herr_t selected = amount[0] ?
            H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start, NULL, amount, NULL) :
            H5Sselect_none(fileSpace);
          require(selected >= 0, "Cannot select checkpoint data hyperslab");
          Handle memory(amount[0] ? H5Screate_simple(1, amount, NULL) :
            H5Screate(H5S_NULL), H5Sclose);
          require(memory >= 0, "Cannot create checkpoint transfer memory space");
          unsigned char empty = 0;
          void* buffer = amount[0] ? static_cast<void*>(
            static_cast<unsigned char*>(data) +
            static_cast<std::size_t>(localStart) * memoryWidth) : &empty;
          const herr_t status = writing ? H5Dwrite(dataset, memoryType, memory,
            fileSpace, transferProperties, buffer) : H5Dread(dataset,
            memoryType, memory, fileSpace, transferProperties, buffer);
          require(status >= 0, "Parallel checkpoint data transfer failed");
          const std::uint64_t bytes = static_cast<std::uint64_t>(amount[0]) *
            static_cast<std::uint64_t>(fileWidth(block.kind));
          require(bytes <= std::numeric_limits<std::uint64_t>::max() - transferred_,
                  "Checkpoint transfer accounting overflow");
          transferred_ += bytes;
        }
    }

    hid_t particleType(bool file)
    {
      const bool trueValue = true;
      const bool falseValue = false;
      require(sizeof(bool) == 1 &&
        *reinterpret_cast<const unsigned char*>(&trueValue) == 1 &&
        *reinterpret_cast<const unsigned char*>(&falseValue) == 0,
        "Checkpoint particles require one-byte zero/one bool storage");
      RelativisticParticleSI particle;
      const hid_t type = H5Tcreate(H5T_COMPOUND,
        file ? particleFileBytes : sizeof(RelativisticParticleSI));
      require(type >= 0, "Cannot create checkpoint particle datatype");
      std::size_t offset = 0;
      const auto insert = [&](const char* name, const void* address,
                               hid_t native, hid_t disk, std::size_t width) {
        const std::size_t memoryOffset = static_cast<std::size_t>(
          static_cast<const unsigned char*>(address) -
          reinterpret_cast<const unsigned char*>(&particle));
        require(H5Tinsert(type, name, file ? offset : memoryOffset,
          file ? disk : native) >= 0, "Cannot define checkpoint particle field");
        offset += width;
      };
      insert("x", &particle.position[0], H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, 8);
      insert("y", &particle.position[1], H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, 8);
      insert("z", &particle.position[2], H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, 8);
      insert("ux", &particle.properVelocity[0], H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, 8);
      insert("uy", &particle.properVelocity[1], H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, 8);
      insert("uz", &particle.properVelocity[2], H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, 8);
      insert("charge", &particle.charge, H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, 8);
      insert("mass", &particle.mass, H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, 8);
      insert("weight", &particle.weight, H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, 8);
      insert("id", &particle.id, H5T_NATIVE_UINT64, H5T_STD_U64LE, 8);
      insert("source_id", &particle.sourceId, H5T_NATIVE_UINT64, H5T_STD_U64LE, 8);
      insert("frozen", &particle.frozen, H5T_NATIVE_UCHAR, H5T_STD_U8LE, 1);
      insert("frozen_beta_z", &particle.frozenBetaZ, H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, 8);
      return type;
    }

    hid_t carrierType(bool file)
    {
      ParticleCPMLCarrier carrier;
      const hid_t type = H5Tcreate(H5T_COMPOUND,
        file ? carrierFileBytes : sizeof(ParticleCPMLCarrier));
      require(type >= 0, "Cannot create checkpoint carrier datatype");
      const char* names[8] = {"x", "y", "z", "vx", "vy", "vz",
                              "charge", "current_weight"};
      const double* fields[8] = {&carrier.position[0], &carrier.position[1],
        &carrier.position[2], &carrier.velocity[0], &carrier.velocity[1],
        &carrier.velocity[2], &carrier.charge, &carrier.currentWeight};
      for (std::size_t index = 0; index < 8; ++index)
        require(H5Tinsert(type, names[index], file ? index * 8 :
          fieldOffset(carrier, *fields[index]), file ? H5T_IEEE_F64LE :
          H5T_NATIVE_DOUBLE) >= 0, "Cannot define checkpoint carrier field");
      return type;
    }

    template<typename Record>
    void resize(std::vector<Record>& records, unsigned long long count)
    {
      require(count <= std::numeric_limits<std::size_t>::max() / sizeof(Record),
              "Checkpoint restored record count exceeds memory address range");
      bool allocated = true;
      try { records.resize(static_cast<std::size_t>(count)); }
      catch (...) { allocated = false; }
      require(allocated, "Cannot allocate checkpoint restored records");
    }

    void close()
    {
      require(!closed_, "Checkpoint archive is already closed");
      const herr_t result = H5Fclose(file_);
      file_ = -1;
      closed_ = true;
      require(result >= 0, "Cannot close checkpoint on all ranks");
    }

    void commit()
    {
      require(writing_ && !closed_, "Only an open checkpoint writer can commit");
      makeHeader(true);
      {
        Handle header(H5Dopen2(file_, "/header", H5P_DEFAULT), H5Dclose);
        require(header >= 0, "Cannot reopen checkpoint completion header");
        rootTransfer(header, H5T_NATIVE_ULLONG, header_, headerCount, true);
      }
      require(H5Fflush(file_, H5F_SCOPE_GLOBAL) >= 0,
              "Cannot flush checkpoint before publication");
      close();
      /* close() includes a collective success check. Publication only starts
       * once every rank has closed the shared file. link() atomically refuses
       * any existing final path, unlike rename(), including a racing writer. */
      bool published = true;
      if (rank_ == 0) published = ::link(partial_.c_str(), path_.c_str()) == 0;
      require(published, "Cannot publish checkpoint without replacing an existing file: " + path_);
      bool removed = true;
      if (rank_ == 0) removed = ::unlink(partial_.c_str()) == 0;
      require(removed,
              "Checkpoint published, but its partial hard link could not be removed");
    }

    MPI_Comm communicator_;
    std::string path_;
    std::string partial_;
    bool writing_;
    bool closed_;
    int rank_;
    int size_;
    hid_t file_;
    unsigned long long header_[headerCount];
    std::uint64_t transferred_;
  };

  CheckpointArchive::CheckpointArchive(MPI_Comm communicator,
      const std::string& path, bool writing)
    : impl_(new Impl(communicator, path, writing)) {}
  CheckpointArchive::~CheckpointArchive() = default;

  void CheckpointArchive::writeBytes(const std::string& name,
      const void* data, std::size_t bytes)
  {
    impl_->write(name, data, bytes, ByteBlock, H5T_NATIVE_UCHAR,
                  H5T_STD_U8LE, 1);
  }
  void CheckpointArchive::readBytes(const std::string& name,
      void* data, std::size_t bytes)
  {
    impl_->read(name, data, bytes, ByteBlock, H5T_NATIVE_UCHAR,
                 H5T_STD_U8LE, 1);
  }
  std::vector<unsigned char> CheckpointArchive::readLocalBytes(
      const std::string& name)
  {
    const Impl::Block block = impl_->inspect(name);
    impl_->require(block.kind == ByteBlock, "Checkpoint block is not native bytes");
    std::vector<unsigned char> result;
    impl_->resize(result, block.count);
    readBytes(name, result.empty() ? NULL : result.data(), result.size());
    return result;
  }
  void CheckpointArchive::writeDoubles(const std::string& name,
      const double* data, std::size_t count)
  {
    impl_->write(name, data, count, DoubleBlock, H5T_NATIVE_DOUBLE,
                  H5T_IEEE_F64LE, sizeof(double));
  }
  void CheckpointArchive::readDoubles(const std::string& name,
      double* data, std::size_t count)
  {
    impl_->read(name, data, count, DoubleBlock, H5T_NATIVE_DOUBLE,
                 H5T_IEEE_F64LE, sizeof(double));
  }
  void CheckpointArchive::writeParticles(const std::string& name,
      const std::vector<RelativisticParticleSI>& particles)
  {
    Handle memory(impl_->particleType(false), H5Tclose);
    Handle disk(impl_->particleType(true), H5Tclose);
    impl_->write(name, particles.empty() ? NULL : particles.data(),
      particles.size(), ParticleBlock, memory, disk, sizeof(RelativisticParticleSI));
  }
  void CheckpointArchive::readParticles(const std::string& name,
      std::vector<RelativisticParticleSI>& particles)
  {
    const Impl::Block block = impl_->inspect(name);
    impl_->require(block.kind == ParticleBlock, "Checkpoint block is not particles");
    impl_->resize(particles, block.count);
    Handle memory(impl_->particleType(false), H5Tclose);
    Handle disk(impl_->particleType(true), H5Tclose);
    impl_->read(name, particles.empty() ? NULL : particles.data(), particles.size(),
      ParticleBlock, memory, disk, sizeof(RelativisticParticleSI));
    bool validFrozen = true;
    for (std::size_t index = 0; index < particles.size(); ++index)
      {
        const unsigned char raw = *reinterpret_cast<const unsigned char*>(
          &particles[index].frozen);
        validFrozen = validFrozen && raw <= 1;
      }
    impl_->require(validFrozen, "Checkpoint particle frozen flag is invalid");
  }
  void CheckpointArchive::writeCarriers(const std::string& name,
      const std::vector<ParticleCPMLCarrier>& carriers)
  {
    Handle memory(impl_->carrierType(false), H5Tclose);
    Handle disk(impl_->carrierType(true), H5Tclose);
    impl_->write(name, carriers.empty() ? NULL : carriers.data(), carriers.size(),
      CarrierBlock, memory, disk, sizeof(ParticleCPMLCarrier));
  }
  void CheckpointArchive::readCarriers(const std::string& name,
      std::vector<ParticleCPMLCarrier>& carriers)
  {
    const Impl::Block block = impl_->inspect(name);
    impl_->require(block.kind == CarrierBlock, "Checkpoint block is not carriers");
    impl_->resize(carriers, block.count);
    Handle memory(impl_->carrierType(false), H5Tclose);
    Handle disk(impl_->carrierType(true), H5Tclose);
    impl_->read(name, carriers.empty() ? NULL : carriers.data(), carriers.size(),
      CarrierBlock, memory, disk, sizeof(ParticleCPMLCarrier));
  }
  std::size_t CheckpointArchive::localSize(const std::string& name)
  {
    const Impl::Block block = impl_->inspect(name);
    return static_cast<std::size_t>(block.count) * fileWidth(block.kind);
  }
  std::uint64_t CheckpointArchive::localTransferredBytes() const
  { return impl_->transferred_; }
  void CheckpointArchive::commit() { impl_->commit(); }
  void CheckpointArchive::close() { impl_->close(); }
}
