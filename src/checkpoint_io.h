#ifndef APRL_CHECKPOINT_IO_H
#define APRL_CHECKPOINT_IO_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <mpi.h>

#include "eb_particles.h"
#include "particle_boundary.h"

namespace aprl
{
  /* Every operation, including localSize(), is collective and must be called
   * in the same order on every rank. Data remains distributed. A writer owns
   * path + ".partial" until commit() closes all ranks and publishes path.
   * Existing final or partial files are never overwritten.
   *
   * The destructor deliberately performs no collective HDF5 close: an error
   * may unwind only one rank, whose caller must then abort the MPI job. Call
   * commit() on a successful writer, or close() on a successful reader.
   * The caller supplies a local, MPI-free timeout guard around archive I/O.
   * Native byte blocks are for ABI-checked state, not a portable wire format.
   */
  class CheckpointArchive
  {
  public:
    CheckpointArchive(MPI_Comm communicator, const std::string& path,
                      bool writing);
    ~CheckpointArchive();
    CheckpointArchive(const CheckpointArchive&) = delete;
    CheckpointArchive& operator=(const CheckpointArchive&) = delete;

    void writeBytes(const std::string& name, const void* data,
                    std::size_t bytes);
    void readBytes(const std::string& name, void* data, std::size_t bytes);
    std::vector<unsigned char> readLocalBytes(const std::string& name);
    void writeDoubles(const std::string& name, const double* data,
                      std::size_t count);
    void readDoubles(const std::string& name, double* data, std::size_t count);
    void writeParticles(const std::string& name,
                        const std::vector<RelativisticParticleSI>& particles);
    void readParticles(const std::string& name,
                       std::vector<RelativisticParticleSI>& particles);
    void writeCarriers(const std::string& name,
                       const std::vector<ParticleCPMLCarrier>& carriers);
    void readCarriers(const std::string& name,
                      std::vector<ParticleCPMLCarrier>& carriers);

    /* Number of packed on-disk bytes belonging to this rank. For particles
     * and carriers, use their typed read methods to allocate native records. */
    std::size_t localSize(const std::string& name);
    std::uint64_t localTransferredBytes() const;
    void commit();
    /* A writer closed without commit retains its unpublished partial file. */
    void close();

  private:
    class Impl;
    std::unique_ptr<Impl> impl_;
  };
}

#endif
