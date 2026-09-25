#ifndef DIRECT_EB_PARTICLE_FILE_H
#define DIRECT_EB_PARTICLE_FILE_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <mpi.h>

#include "eb_particles.h"

namespace fel
{
  struct ParticleInputRecord
  {
    Double position[3];
    Double properVelocity[3];
    std::uint64_t sourceId;
    Double macroWeight;

    ParticleInputRecord();
  };

  class ParticleHdf5File
  {
  public:
    static const int legacyFormatVersion = 1;
    static const int snapshotFormatVersion = 2;
    static const int formatVersion = 3;

    /* Reads one contiguous hyperslab per MPI rank. Parallel HDF5 uses one
     * collective file open and one collective dataset read. */
    static std::vector<RelativisticParticleSI> readDistributed(
        const std::string& filename,
        Double totalElectrons,
        const FieldVector<Double>& positionOffsetSI,
        MPI_Comm communicator,
        unsigned long long& globalRecords,
        int& inputFormatVersion);
  };
}

#endif
