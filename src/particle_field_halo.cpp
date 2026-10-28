#include "particle_field_halo.h"

#include <algorithm>
#include <climits>
#include <limits>
#include <stdexcept>

namespace aprl
{
  namespace
  {
    std::size_t planeSize(std::size_t nx, std::size_t ny)
    {
      if (nx == 0 || ny == 0 ||
          nx > static_cast<std::size_t>(INT_MAX) / ny)
        throw std::overflow_error("Particle interpolation plane exceeds MPI count");
      return nx * ny;
    }
  }

  EBParticleFieldHalo::EBParticleFieldHalo(
      const EBGridGeometry& geometry, MPI_Comm communicator)
    : geometry_(geometry), communicator_(communicator),
      lowerRank_(MPI_PROC_NULL), upperRank_(MPI_PROC_NULL),
      ezCount_(0), bxCount_(0), planeCount_(0), installed_(0),
      communicationSeconds_(0.0), exchanges_(0), sentBytes_(0),
      receivedBytes_(0)
  {
    if (geometry.nx < 2 || geometry.ny < 2 || geometry.nz < 2 ||
        geometry.nx == std::numeric_limits<std::size_t>::max() ||
        geometry.ny == std::numeric_limits<std::size_t>::max())
      throw std::invalid_argument("Particle interpolation halo needs a valid slab");
    ezCount_ = planeSize(geometry.nx + 1, geometry.ny + 1);
    bxCount_ = planeSize(geometry.nx + 1, geometry.ny);
    const std::size_t byCount = planeSize(geometry.nx, geometry.ny + 1);
    if (ezCount_ > static_cast<std::size_t>(INT_MAX) - bxCount_ ||
        ezCount_ + bxCount_ > static_cast<std::size_t>(INT_MAX) - byCount)
      throw std::overflow_error("Packed particle interpolation planes exceed MPI count");
    planeCount_ = ezCount_ + bxCount_ + byCount;
    int rank = 0, size = 0;
    MPI_Comm_rank(communicator_, &rank);
    MPI_Comm_size(communicator_, &size);
    lowerRank_ = rank == 0 ? MPI_PROC_NULL : rank - 1;
    upperRank_ = rank + 1 == size ? MPI_PROC_NULL : rank + 1;
    // nx/ny and native spacing must match. Longitudinal slab lengths may differ.
    const Double local[6] = {static_cast<Double>(geometry.nx),
      static_cast<Double>(geometry.ny), geometry.dx, geometry.dy,
      geometry.dz, geometry.dt};
    Double minimum[6], maximum[6];
    MPI_Allreduce(local, minimum, 6, MPI_DOUBLE, MPI_MIN, communicator_);
    MPI_Allreduce(local, maximum, 6, MPI_DOUBLE, MPI_MAX, communicator_);
    for (unsigned i = 0; i < 6; ++i)
      if (minimum[i] != maximum[i])
        throw std::invalid_argument("Particle interpolation halo slab geometry differs across ranks");
    if (size > 1) send_.resize(planeCount_);
    if (lowerRank_ != MPI_PROC_NULL) lower_.resize(planeCount_);
    if (upperRank_ != MPI_PROC_NULL) upper_.resize(planeCount_);
  }

  void EBParticleFieldHalo::verifyGeometry(const EBFieldGrid& fields) const
  {
    const EBGridGeometry& g = fields.geometry();
    if (g.nx != geometry_.nx || g.ny != geometry_.ny || g.nz != geometry_.nz ||
        g.dx != geometry_.dx || g.dy != geometry_.dy ||
        g.dz != geometry_.dz || g.dt != geometry_.dt ||
        g.solver != geometry_.solver)
      throw std::invalid_argument("Particle interpolation halo grid changed");
  }

  void EBParticleFieldHalo::pack(const EBFieldGrid& fields, std::size_t k)
  {
    const YeeComponent* c[3] = {&fields.ez(), &fields.bx(), &fields.by()};
    std::size_t offset = 0;
    for (unsigned a = 0; a < 3; ++a)
      {
        const std::size_t count = c[a]->nx() * c[a]->ny();
        std::copy(c[a]->data() + k * count,
                  c[a]->data() + (k + 1) * count, send_.begin() + offset);
        offset += count;
      }
  }

  void EBParticleFieldHalo::exchangeAndInstall(EBFieldGrid& fields)
  {
    verifyGeometry(fields);
    if (installed_ || fields.hasParticleInterpolationHalo())
      throw std::logic_error("Particle interpolation halo is already installed");
    EBParticleInterpolationHaloView view;
    if (!send_.empty())
      {
        // Send the uppermost cell centre to the upper rank, receive the lower
        // neighbour's uppermost centre. Then perform the opposite direction.
        pack(fields, geometry_.nz - 1);
        Double started = MPI_Wtime();
        const int first = MPI_Sendrecv(send_.data(), static_cast<int>(planeCount_),
          MPI_DOUBLE, upperRank_, 19301, lower_.empty() ? 0 : lower_.data(),
          lower_.empty() ? 0 : static_cast<int>(planeCount_),
          MPI_DOUBLE, lowerRank_, 19301,
          communicator_, MPI_STATUS_IGNORE);
        communicationSeconds_ += MPI_Wtime() - started;
        if (first != MPI_SUCCESS)
          throw std::runtime_error("Particle lower interpolation halo exchange failed");
        pack(fields, 0);
        started = MPI_Wtime();
        const int second = MPI_Sendrecv(send_.data(), static_cast<int>(planeCount_),
          MPI_DOUBLE, lowerRank_, 19302, upper_.empty() ? 0 : upper_.data(),
          upper_.empty() ? 0 : static_cast<int>(planeCount_),
          MPI_DOUBLE, upperRank_, 19302,
          communicator_, MPI_STATUS_IGNORE);
        communicationSeconds_ += MPI_Wtime() - started;
        if (second != MPI_SUCCESS)
          throw std::runtime_error("Particle upper interpolation halo exchange failed");
        const std::size_t peers = (lowerRank_ != MPI_PROC_NULL ? 1 : 0) +
          (upperRank_ != MPI_PROC_NULL ? 1 : 0);
        sentBytes_ += peers * planeCount_ * sizeof(Double);
        receivedBytes_ += peers * planeCount_ * sizeof(Double);
      }
    if (!lower_.empty())
      {
        view.lowerEz = lower_.data();
        view.lowerBx = lower_.data() + ezCount_;
        view.lowerBy = lower_.data() + ezCount_ + bxCount_;
      }
    if (!upper_.empty())
      {
        view.upperEz = upper_.data();
        view.upperBx = upper_.data() + ezCount_;
        view.upperBy = upper_.data() + ezCount_ + bxCount_;
      }
    fields.setParticleInterpolationHalo(view);
    installed_ = &fields;
    ++exchanges_;
  }

  void EBParticleFieldHalo::clear(EBFieldGrid& fields)
  {
    if (installed_ != &fields)
      throw std::logic_error("Particle interpolation halo clear has no matching install");
    fields.clearParticleInterpolationHalo();
    installed_ = 0;
  }

  std::size_t EBParticleFieldHalo::memoryBytes() const
  {
    return (send_.capacity() + lower_.capacity() + upper_.capacity()) * sizeof(Double);
  }
}
