#include "trajectory_reader.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "hdf5.h"

namespace radiation
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

    int readIntAttribute(hid_t object, const char* name)
    {
      hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
      requireHandle(attribute, std::string("Missing HDF5 attribute: ") + name);
      int value = 0;
      const herr_t status = H5Aread(attribute, H5T_NATIVE_INT, &value);
      H5Aclose(attribute);
      requireStatus(status, std::string("Cannot read HDF5 attribute: ") + name);
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

    hid_t createMemoryRecordType(int formatVersion)
    {
      hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(TrajectoryRecord));
      requireHandle(type, "Cannot create trajectory memory datatype");
      hsize_t vectorSize[1] = {3};
      hid_t vectorType = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, vectorSize);
      requireHandle(vectorType, "Cannot create trajectory vector datatype");
      try
        {
          requireStatus(H5Tinsert(type, "particle_id",
            HOFFSET(TrajectoryRecord, particleId), H5T_NATIVE_UINT64),
            "Cannot map particle_id");
          requireStatus(H5Tinsert(type, "source_id",
            HOFFSET(TrajectoryRecord, sourceId), H5T_NATIVE_UINT64),
            "Cannot map source_id");
          requireStatus(H5Tinsert(type, "time_s",
            HOFFSET(TrajectoryRecord, time), H5T_NATIVE_DOUBLE),
            "Cannot map time_s");
          requireStatus(H5Tinsert(type, "position_m",
            HOFFSET(TrajectoryRecord, position), vectorType),
            "Cannot map position_m");
          requireStatus(H5Tinsert(type, "proper_velocity",
            HOFFSET(TrajectoryRecord, properVelocity), vectorType),
            "Cannot map proper_velocity");
          requireStatus(H5Tinsert(type, "charge_C",
            HOFFSET(TrajectoryRecord, charge), H5T_NATIVE_DOUBLE),
            "Cannot map charge_C");
          requireStatus(H5Tinsert(type, "weight",
            HOFFSET(TrajectoryRecord, weight), H5T_NATIVE_DOUBLE),
            "Cannot map weight");
          if (formatVersion >= 2)
            {
              requireStatus(H5Tinsert(type, "event_type",
                HOFFSET(TrajectoryRecord, event), H5T_NATIVE_UCHAR),
                "Cannot map event_type");
              requireStatus(H5Tinsert(type, "boundary_face",
                HOFFSET(TrajectoryRecord, boundaryFace), H5T_NATIVE_SCHAR),
                "Cannot map boundary_face");
            }
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

    void validateRecord(const TrajectoryRecord& record,
                        const std::string& filename)
    {
      if (record.particleId == 0 || !std::isfinite(record.time) ||
          !std::isfinite(record.charge) || !std::isfinite(record.weight) ||
          record.event > 2 || record.boundaryFace < -1 ||
          record.boundaryFace > 5)
        throw std::runtime_error(
          "Invalid trajectory record in " + filename);
      for (unsigned int component = 0; component < 3; ++component)
        if (!std::isfinite(record.position[component]) ||
            !std::isfinite(record.properVelocity[component]))
          throw std::runtime_error(
            "Non-finite trajectory vector in " + filename);
    }

    void readFile(const std::string& filename, bool requireComplete,
                  std::size_t chunkRecords,
                  std::vector<TrajectoryRecord>& destination)
    {
      hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      requireHandle(file, "Cannot open trajectory file: " + filename);
      hid_t group = -1;
      hid_t dataset = -1;
      hid_t memoryType = -1;
      try
        {
          group = H5Gopen2(file, "/trajectory", H5P_DEFAULT);
          requireHandle(group, "Missing /trajectory group in " + filename);
          const int formatVersion = readIntAttribute(
            group, "format_version");
          if (formatVersion != 1 && formatVersion != 2)
            throw std::runtime_error(
              "Trajectory radiation supports format_version 1 or 2: " +
              filename);
          if (requireComplete && readByteDataset(group, "complete") == 0)
            throw std::runtime_error(
              "Trajectory file is marked incomplete: " + filename);
          const std::uint64_t committed = readUnsignedDataset(
            group, "committed_records");
          dataset = H5Dopen2(group, "records", H5P_DEFAULT);
          requireHandle(dataset, "Missing trajectory records in " + filename);
          hid_t space = H5Dget_space(dataset);
          requireHandle(space, "Cannot inspect trajectory records in " + filename);
          hsize_t extent = 0;
          const int dimensions = H5Sget_simple_extent_dims(space, &extent, NULL);
          H5Sclose(space);
          if (dimensions != 1 || committed > extent)
            throw std::runtime_error(
              "Invalid committed trajectory prefix in " + filename);
          memoryType = createMemoryRecordType(formatVersion);
          std::vector<TrajectoryRecord> buffer;
          for (std::uint64_t offset = 0; offset < committed;)
            {
              const std::uint64_t remaining = committed - offset;
              const std::size_t count = static_cast<std::size_t>(
                std::min<std::uint64_t>(remaining, chunkRecords));
              buffer.resize(count);
              hid_t fileSpace = H5Dget_space(dataset);
              requireHandle(fileSpace, "Cannot select trajectory file space");
              hsize_t start[1] = {static_cast<hsize_t>(offset)};
              hsize_t size[1] = {static_cast<hsize_t>(count)};
              requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
                start, NULL, size, NULL), "Cannot select trajectory chunk");
              hid_t memorySpace = H5Screate_simple(1, size, NULL);
              requireHandle(memorySpace, "Cannot create trajectory memory space");
              const herr_t status = H5Dread(dataset, memoryType, memorySpace,
                fileSpace, H5P_DEFAULT, buffer.empty() ? NULL : &buffer[0]);
              H5Sclose(memorySpace);
              H5Sclose(fileSpace);
              requireStatus(status, "Cannot read trajectory chunk: " + filename);
              for (std::size_t index = 0; index < buffer.size(); ++index)
                {
                  if (formatVersion == 1)
                    {
                      buffer[index].event = 0;
                      buffer[index].boundaryFace = -1;
                    }
                  validateRecord(buffer[index], filename);
                }
              destination.insert(destination.end(), buffer.begin(), buffer.end());
              offset += count;
            }
          H5Tclose(memoryType);
          H5Dclose(dataset);
          H5Gclose(group);
          H5Fclose(file);
        }
      catch (...)
        {
          if (memoryType >= 0) H5Tclose(memoryType);
          if (dataset >= 0) H5Dclose(dataset);
          if (group >= 0) H5Gclose(group);
          H5Fclose(file);
          throw;
        }
    }

    int checkedByteCount(std::size_t records)
    {
      if (records > static_cast<std::size_t>(INT_MAX) /
                    sizeof(TrajectoryRecord))
        throw std::overflow_error(
          "One MPI trajectory exchange exceeds the MPI int byte limit");
      return static_cast<int>(records * sizeof(TrajectoryRecord));
    }

    bool recordLess(const TrajectoryRecord& left,
                    const TrajectoryRecord& right)
    {
      if (left.particleId != right.particleId)
        return left.particleId < right.particleId;
      if (left.time != right.time) return left.time < right.time;
      return left.event < right.event;
    }

    bool sameEventTime(double left, double right)
    {
      const double scale = std::max(std::numeric_limits<double>::min(),
        std::max(std::abs(left), std::abs(right)));
      return std::abs(left - right) <=
        64.0 * std::numeric_limits<double>::epsilon() * scale;
    }
  }

  std::vector<ParticleTrajectory> loadShotDistributed(
      const ShotConfig& shot, bool requireComplete,
      std::size_t chunkRecords, MPI_Comm communicator,
      TrajectoryLoadStats& globalStats)
  {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    std::vector<TrajectoryRecord> localInput;
    for (std::size_t file = 0; file < shot.files.size(); ++file)
      if (static_cast<int>(file % static_cast<std::size_t>(size)) == rank)
        readFile(shot.files[file], requireComplete, chunkRecords, localInput);

    const unsigned long long localInputCount =
      static_cast<unsigned long long>(localInput.size());
    MPI_Allreduce(&localInputCount, &globalStats.inputRecords, 1,
      MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

    std::vector<TrajectoryRecord> received;
    if (size == 1)
      {
        // Avoid MPI's signed-int byte-count limit on the small-server path.
        received.swap(localInput);
      }
    else
      {
        std::vector<int> sendCounts(static_cast<std::size_t>(size), 0);
        for (std::size_t index = 0; index < localInput.size(); ++index)
          {
            const int owner = static_cast<int>(localInput[index].particleId %
              static_cast<std::uint64_t>(size));
            if (sendCounts[static_cast<std::size_t>(owner)] == INT_MAX)
              throw std::overflow_error(
                "Too many trajectory records for one MPI peer");
            ++sendCounts[static_cast<std::size_t>(owner)];
          }
        std::vector<int> receiveCounts(static_cast<std::size_t>(size), 0);
        MPI_Alltoall(&sendCounts[0], 1, MPI_INT,
                     &receiveCounts[0], 1, MPI_INT, communicator);
        std::vector<int> sendDisplacements(static_cast<std::size_t>(size), 0);
        std::vector<int> receiveDisplacements(
          static_cast<std::size_t>(size), 0);
        for (int process = 1; process < size; ++process)
          {
            sendDisplacements[process] = sendDisplacements[process - 1] +
              sendCounts[process - 1];
            receiveDisplacements[process] =
              receiveDisplacements[process - 1] +
              receiveCounts[process - 1];
          }
        const std::size_t receiveTotal = static_cast<std::size_t>(
          receiveDisplacements[size - 1] + receiveCounts[size - 1]);
        std::vector<TrajectoryRecord> send(localInput.size());
        std::vector<int> cursor(sendDisplacements);
        for (std::size_t index = 0; index < localInput.size(); ++index)
          {
            const int owner = static_cast<int>(
              localInput[index].particleId %
              static_cast<std::uint64_t>(size));
            send[static_cast<std::size_t>(
              cursor[static_cast<std::size_t>(owner)]++)] = localInput[index];
          }
        localInput.clear();
        localInput.shrink_to_fit();
        received.resize(receiveTotal);
        std::vector<int> sendBytes(sendCounts);
        std::vector<int> receiveBytes(receiveCounts);
        std::vector<int> sendByteDisplacements(sendDisplacements);
        std::vector<int> receiveByteDisplacements(receiveDisplacements);
        for (int process = 0; process < size; ++process)
          {
            sendBytes[process] = checkedByteCount(
              static_cast<std::size_t>(sendCounts[process]));
            receiveBytes[process] = checkedByteCount(
              static_cast<std::size_t>(receiveCounts[process]));
            sendByteDisplacements[process] = checkedByteCount(
              static_cast<std::size_t>(sendDisplacements[process]));
            receiveByteDisplacements[process] = checkedByteCount(
              static_cast<std::size_t>(receiveDisplacements[process]));
          }
        MPI_Alltoallv(send.empty() ? NULL : &send[0],
          &sendBytes[0], &sendByteDisplacements[0], MPI_BYTE,
          received.empty() ? NULL : &received[0],
          &receiveBytes[0], &receiveByteDisplacements[0], MPI_BYTE,
          communicator);
      }
    std::sort(received.begin(), received.end(), recordLess);

    TrajectoryLoadStats localStats;
    std::vector<ParticleTrajectory> trajectories;
    for (std::size_t index = 0; index < received.size();)
      {
        const std::uint64_t particleId = received[index].particleId;
        ParticleTrajectory trajectory;
        trajectory.particleId = particleId;
        trajectory.sourceId = received[index].sourceId;
        trajectory.charge = received[index].charge;
        bool terminalSeen = false;
        while (index < received.size() &&
               received[index].particleId == particleId)
          {
            const TrajectoryRecord candidate = received[index++];
            if (candidate.sourceId != trajectory.sourceId ||
                std::abs(candidate.charge - trajectory.charge) >
                  64.0 * std::numeric_limits<double>::epsilon() *
                  std::max(std::numeric_limits<double>::min(),
                           std::abs(trajectory.charge)))
              throw std::runtime_error(
                "A particle changes source or charge across trajectory files");
            if (!trajectory.records.empty() && sameEventTime(
                  trajectory.records.back().time, candidate.time))
              {
                ++localStats.duplicateRecords;
                if (candidate.event != 0)
                  trajectory.records.back() = candidate;
                continue;
              }
            if (terminalSeen)
              throw std::runtime_error(
                "A particle has records after a terminal trajectory event");
            if (candidate.event != 0)
              {
                terminalSeen = true;
                ++localStats.terminalEvents;
              }
            trajectory.records.push_back(candidate);
          }
        localStats.uniqueRecords += trajectory.records.size();
        ++localStats.particles;
        trajectories.push_back(trajectory);
      }

    unsigned long long localValues[4] = {
      localStats.uniqueRecords, localStats.particles,
      localStats.terminalEvents, localStats.duplicateRecords
    };
    unsigned long long globalValues[4] = {};
    MPI_Allreduce(localValues, globalValues, 4,
      MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);
    globalStats.uniqueRecords = globalValues[0];
    globalStats.particles = globalValues[1];
    globalStats.terminalEvents = globalValues[2];
    globalStats.duplicateRecords = globalValues[3];
    return trajectories;
  }
}
