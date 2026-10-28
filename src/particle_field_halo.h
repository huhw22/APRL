#ifndef APRL_PARTICLE_FIELD_HALO_H
#define APRL_PARTICLE_FIELD_HALO_H

#include "eb_field.h"
#include <mpi.h>
#include <vector>

namespace aprl
{
  /* Supplies the missing half-cell Ez/Bx/By interpolation neighbours for a
   * fixed z-slab decomposition. No field values or push formula are changed.
   * The caller owns the field-step lifetime: install after B, clear after all
   * particle pushes (including transfers), before E/diagnostics/window shift.
   * The installed grid must be cleared before this buffer owner is destroyed. */
  class EBParticleFieldHalo
  {
  public:
    EBParticleFieldHalo(const EBGridGeometry& geometry, MPI_Comm communicator);
    EBParticleFieldHalo(const EBParticleFieldHalo&) = delete;
    EBParticleFieldHalo& operator=(const EBParticleFieldHalo&) = delete;

    void exchangeAndInstall(EBFieldGrid& fields);
    void clear(EBFieldGrid& fields);
    std::size_t memoryBytes() const;
    Double communicationSeconds() const { return communicationSeconds_; }
    std::size_t exchanges() const { return exchanges_; }
    std::size_t sentBytes() const { return sentBytes_; }
    std::size_t receivedBytes() const { return receivedBytes_; }

  private:
    void verifyGeometry(const EBFieldGrid& fields) const;
    void pack(const EBFieldGrid& fields, std::size_t k);

    EBGridGeometry geometry_;
    MPI_Comm communicator_;
    int lowerRank_, upperRank_;
    std::size_t ezCount_, bxCount_, planeCount_;
    std::vector<Double> send_, lower_, upper_;
    EBFieldGrid* installed_;
    Double communicationSeconds_;
    std::size_t exchanges_, sentBytes_, receivedBytes_;
  };
}
#endif
