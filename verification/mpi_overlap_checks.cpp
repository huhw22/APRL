#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

#include <mpi.h>

#include "eb_cpml.h"
#include "eb_field.h"
#include "eb_mpi.h"

namespace
{
  void require(bool condition, const char* message)
  {
    if (!condition) throw std::runtime_error(message);
  }

  void initialize(aprl::EBFieldGrid& fields, std::size_t offset)
  {
    aprl::YeeComponent* component[] = {
      &fields.ex(), &fields.ey(), &fields.ez(),
      &fields.bx(), &fields.by(), &fields.bz()
    };
    for (unsigned int axis = 0; axis < 6; ++axis)
      for (std::size_t k = 0; k < component[axis]->nz(); ++k)
        for (std::size_t j = 0; j < component[axis]->ny(); ++j)
          for (std::size_t i = 0; i < component[axis]->nx(); ++i)
            (*component[axis])(i, j, k) =
              (std::sin(0.13 * static_cast<double>(i + 2 * j + 3 * (k + offset))
                 + 0.7 * axis) + 0.25 * std::cos(0.21 * (i + j))) /
              (axis < 3 ? 1.0 : aprl::SI::c);
  }

  void setCurrent(aprl::EBFieldGrid& fields, std::size_t offset,
                  std::size_t globalNz, unsigned int step)
  {
    aprl::YeeComponent* component[] = {
      &fields.jx(), &fields.jy(), &fields.jz()
    };
    for (unsigned int axis = 0; axis < 3; ++axis)
      for (std::size_t k = 0; k < component[axis]->nz(); ++k)
        {
          // Reconstruct one global deposited current when the shared
          // tangential vertex plane is summed by neighbouring slabs.
          const bool shared = axis < 2 &&
            ((k == 0 && offset > 0) ||
             (k == fields.geometry().nz && k + offset < globalNz));
          for (std::size_t j = 0; j < component[axis]->ny(); ++j)
            for (std::size_t i = 0; i < component[axis]->nx(); ++i)
              (*component[axis])(i, j, k) = (shared ? 0.5 : 1.0) *
                (0.02 + 0.01 * std::sin(0.17 * (i + 3 * j + k + offset)
                   + 0.31 * step + axis));
        }
  }

  double difference(const aprl::EBFieldGrid& local,
                    const aprl::EBFieldGrid& reference,
                    std::size_t offset)
  {
    const aprl::YeeComponent* found[] = {
      &local.ex(), &local.ey(), &local.ez(),
      &local.bx(), &local.by(), &local.bz()
    };
    const aprl::YeeComponent* expected[] = {
      &reference.ex(), &reference.ey(), &reference.ez(),
      &reference.bx(), &reference.by(), &reference.bz()
    };
    double error = 0.0;
    for (unsigned int axis = 0; axis < 6; ++axis)
      for (std::size_t k = 0; k < found[axis]->nz(); ++k)
        for (std::size_t j = 0; j < found[axis]->ny(); ++j)
          for (std::size_t i = 0; i < found[axis]->nx(); ++i)
            {
              const double scale = axis < 3 ? 1.0 : aprl::SI::c;
              const double value = scale * (*found[axis])(i, j, k);
              const double target = scale * (*expected[axis])(i, j, k + offset);
              require(std::isfinite(value) && std::isfinite(target),
                "Nonfinite field in MPI overlap comparison");
              error = std::max(error, std::abs(value - target) /
                std::max(1.0, std::abs(target)));
            }
    return error;
  }

  void compareWithWholeGrid(int rank, int size, aprl::EBMaxwellSolver solver,
                           aprl::EBCurrentFilter filter, bool cpml)
  {
    const std::size_t globalNz = 6 * static_cast<std::size_t>(size) + 1;
    const std::size_t localNz = 6 + (rank == 0 ? 1 : 0);
    const std::size_t offset = 6 * static_cast<std::size_t>(rank) +
      (rank == 0 ? 0 : 1);
    const double dz = 1.0e-5;
    const double dt = dz / aprl::SI::c *
      (solver == aprl::EBMaxwellSolver::Yee ? 0.4 : 1.0);
    aprl::EBGridGeometry global(8, 7, globalNz, 2.0e-5, 2.5e-5, dz, dt, solver);
    aprl::EBGridGeometry slab(global);
    slab.nz = localNz;
    aprl::EBFieldGrid reference(global);
    aprl::EBFieldGrid local(slab);
    aprl::EBZSlabHaloExchange halo(MPI_COMM_WORLD);
    aprl::EBZSlabHaloExchange serialHalo(MPI_COMM_SELF);
    if (cpml)
      {
        aprl::EBCPMLParameters parameters;
        for (unsigned int axis = 0; axis < 3; ++axis)
          parameters.cells[axis] = 2;
        reference.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(
          new aprl::EBConvolutionalPML(global, globalNz, 0, parameters, true, true)));
        local.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(
          new aprl::EBConvolutionalPML(slab, globalNz, offset, parameters,
            rank == 0, rank + 1 == size)));
      }
    else
      {
        halo.installPhysicalBoundaryMask(local);
        serialHalo.installPhysicalBoundaryMask(reference);
      }
    initialize(reference, 0);
    initialize(local, offset);
    double localError = 0.0;
    double previousCommunication = halo.communicationSeconds();
    for (unsigned int step = 0; step < 12; ++step)
      {
        setCurrent(reference, 0, globalNz, step);
        setCurrent(local, offset, globalNz, step);
        serialHalo.advanceMagnetic(reference);
        serialHalo.advanceElectric(reference, filter);
        halo.advanceMagnetic(local);
        halo.advanceElectric(local, filter);
        localError = std::max(localError, difference(local, reference, offset));
        require(halo.communicationSeconds() >= previousCommunication,
          "Halo communication timer decreased");
        previousCommunication = halo.communicationSeconds();
      }
    double globalError = 0.0;
    MPI_Allreduce(&localError, &globalError, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    require(globalError < 5.0e-12,
      "Overlapped slab update disagrees with the unsplit whole-grid solver");
    require(halo.communicationSeconds() > 0.0,
      "Halo communication timer did not record MPI calls");
    require(serialHalo.communicationSeconds() == 0.0,
      "Single-rank halo recorded communication that never occurred");
    if (rank == 0)
      std::cout << "overlap solver=" << (solver == aprl::EBMaxwellSolver::Yee ?
        "yee" : "cowan-z") << " filter=" << static_cast<int>(filter)
        << " cpml=" << cpml << " max_relative_error=" << globalError << '\n';
  }

  class ThrowingBoundary : public aprl::PerfectElectricConductorBoundary
  {
  public:
    ThrowingBoundary(bool lower, bool upper, bool shouldThrow)
      : aprl::PerfectElectricConductorBoundary(lower, upper),
        shouldThrow_(shouldThrow) {}
    void afterElectricUpdate(aprl::EBFieldGrid& fields) override
    {
      if (shouldThrow_) throw std::runtime_error("Injected field-update failure");
      aprl::PerfectElectricConductorBoundary::afterElectricUpdate(fields);
    }
  private:
    bool shouldThrow_;
  };

  void checkExceptionCompletion(int rank, int size)
  {
    aprl::EBGridGeometry geometry(8, 7, 6, 2.0e-5, 2.5e-5, 1.0e-5,
      0.4e-5 / aprl::SI::c);
    aprl::EBFieldGrid fields(geometry);
    aprl::EBZSlabHaloExchange halo(MPI_COMM_WORLD);
    initialize(fields, 6 * static_cast<std::size_t>(rank));
    fields.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(
      new ThrowingBoundary(rank == 0, rank + 1 == size, rank == 0)));
    bool caught = false;
    try { halo.advanceElectric(fields); }
    catch (const std::runtime_error&) { caught = true; }
    require(caught == (rank == 0), "Injected field-update exception was lost");
    // Reusing the same buffers and tags after the failure must not retain a
    // request from the failed update or consume a stale message.
    halo.installPhysicalBoundaryMask(fields);
    fields.clearFields();
    fields.clearCurrent();
    halo.advanceMagnetic(fields);
    halo.advanceElectric(fields);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

int main(int argc, char** argv)
{
  int provided = MPI_THREAD_SINGLE;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  try
    {
      require(size >= 2, "MPI overlap checks require at least two ranks");
      require(provided >= MPI_THREAD_FUNNELED, "MPI FUNNELED support is required");
      for (unsigned int boundary = 0; boundary < 2; ++boundary)
        {
          compareWithWholeGrid(rank, size, aprl::EBMaxwellSolver::Yee,
            aprl::EBCurrentFilter::None, boundary != 0);
          compareWithWholeGrid(rank, size, aprl::EBMaxwellSolver::CowanZ,
            aprl::EBCurrentFilter::None, boundary != 0);
          compareWithWholeGrid(rank, size, aprl::EBMaxwellSolver::CowanZ,
            aprl::EBCurrentFilter::CompensatedBilinear, boundary != 0);
        }
      checkExceptionCompletion(rank, size);
      if (rank == 0) std::cout << "MPI overlap checks passed\n";
      MPI_Finalize();
      return 0;
    }
  catch (const std::exception& error)
    {
      std::cerr << "rank=" << rank << " " << error.what() << '\n';
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }
}
