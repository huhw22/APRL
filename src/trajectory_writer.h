#ifndef DIRECT_EB_TRAJECTORY_WRITER_H
#define DIRECT_EB_TRAJECTORY_WRITER_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "hdf5.h"

namespace fel
{
  /* Compact radiation-postprocessing record in the laboratory frame.
   * Positions and proper velocities remain float64 because small phase errors
   * are amplified in coherent-radiation reconstruction. */
  struct TrajectoryRecord
  {
    std::uint64_t particleId;
    std::uint64_t sourceId;
    double time;
    double position[3];
    double properVelocity[3];
    double charge;
    double weight;

    TrajectoryRecord();
  };

  /* Buffered append-only HDF5 writer.
   *
   * Each MPI rank owns a separate file, so appending never performs a global
   * synchronization.  Interactive mode makes each committed batch durable in
   * two phases.  Throughput mode keeps the same schema but defers durability
   * to normal close.  Readers must ignore any physical dataset tail beyond
   * committed_records. */
  class TrajectoryWriter
  {
  public:
    TrajectoryWriter();
    ~TrajectoryWriter();

    TrajectoryWriter(const TrajectoryWriter&) = delete;
    TrajectoryWriter& operator=(const TrajectoryWriter&) = delete;

    void open(const std::string& filename,
              int mpiRank, int mpiSize,
              std::size_t bufferRecords = 16384,
              unsigned int compressionLevel = 0,
              bool durableProgress = true);
    void append(const TrajectoryRecord& record);
    void flush();

    /* completed=false marks an intentionally interrupted or partial run. */
    void close(bool completed = false);

    bool isOpen() const;
    bool swmrEnabled() const;
    std::size_t bufferedRecords() const;
    std::uint64_t committedRecords() const;
    std::size_t memoryBytes() const;

  private:
    hid_t createMemoryRecordType() const;
    hid_t createFileRecordType() const;
    void writeUnsignedScalar(hid_t dataset, std::uint64_t value) const;
    void writeByteScalar(hid_t dataset, unsigned char value) const;
    void writeIntAttribute(hid_t object, const char* name, int value) const;
    void writeStringAttribute(hid_t object, const char* name,
                              const std::string& value) const;
    void closeHandles() noexcept;

    bool isOpen_;
    bool swmrEnabled_;
    bool durableProgress_;
    std::size_t bufferLimit_;
    std::uint64_t committedRecords_;
    std::vector<TrajectoryRecord> buffer_;

    hid_t file_;
    hid_t group_;
    hid_t recordDataset_;
    hid_t committedDataset_;
    hid_t completeDataset_;
    hid_t memoryRecordType_;
    hid_t fileRecordType_;
  };
}

#endif
