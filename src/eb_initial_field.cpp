#include "eb_initial_field.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

#include "eb_deposition.h"

namespace fel
{
  namespace
  {
    const int RHO_UPPER = 910;
    const int RHO_LOWER = 911;
    const int SCALAR_UPPER = 912;
    const int SCALAR_LOWER = 913;
    const int GAUSS_EZ_UPPER = 914;

    class ScalarSlab
    {
    public:
      ScalarSlab(std::size_t nxNodes, std::size_t nyNodes,
                 std::size_t ownedZNodes)
        : nx_(nxNodes), ny_(nyNodes), ownedZ_(ownedZNodes),
          values_(checkedValues(nxNodes, nyNodes, ownedZNodes + 2), 0.0)
      {}

      Double& operator()(std::size_t i, std::size_t j, std::size_t p)
      {
        return values_[(p * ny_ + j) * nx_ + i];
      }

      const Double& operator()(std::size_t i, std::size_t j,
                               std::size_t p) const
      {
        return values_[(p * ny_ + j) * nx_ + i];
      }

      Double* plane(std::size_t p)
      {
        return &values_[p * planeValues()];
      }

      std::size_t ownedZ() const { return ownedZ_; }
      std::size_t planeValues() const { return nx_ * ny_; }
      std::size_t bytes() const { return values_.capacity() * sizeof(Double); }

      void fill(Double value)
      {
        std::fill(values_.begin(), values_.end(), value);
      }

    private:
      static std::size_t checkedValues(std::size_t nx, std::size_t ny,
                                       std::size_t nz)
      {
        if (nx != 0 && ny > std::numeric_limits<std::size_t>::max() / nx)
          throw std::overflow_error(
            "Initial-field scalar plane overflows size_t");
        const std::size_t plane = nx * ny;
        if (plane != 0 && nz >
            std::numeric_limits<std::size_t>::max() / plane)
          throw std::overflow_error(
            "Initial-field scalar slab overflows size_t");
        return plane * nz;
      }

      std::size_t nx_;
      std::size_t ny_;
      std::size_t ownedZ_;
      std::vector<Double> values_;
    };

    void communicatorRanks(MPI_Comm communicator, int& rank, int& size)
    {
      if (communicator == MPI_COMM_NULL)
        throw std::invalid_argument(
          "Initial Gauss field communicator cannot be null");
      MPI_Comm_rank(communicator, &rank);
      MPI_Comm_size(communicator, &size);
    }

    int checkedPlaneCount(std::size_t values)
    {
      if (values > static_cast<std::size_t>(INT_MAX))
        throw std::overflow_error(
          "Initial-field scalar halo exceeds MPI int count");
      return static_cast<int>(values);
    }

    void exchangeScalarGhosts(ScalarSlab& scalar, int rank, int size,
                              MPI_Comm communicator)
    {
      const int lower = rank > 0 ? rank - 1 : MPI_PROC_NULL;
      const int upper = rank + 1 < size ? rank + 1 : MPI_PROC_NULL;
      const int count = checkedPlaneCount(scalar.planeValues());
      MPI_Sendrecv(scalar.plane(scalar.ownedZ()), count, MPI_DOUBLE,
                   upper, SCALAR_UPPER,
                   scalar.plane(0), count, MPI_DOUBLE,
                   lower, SCALAR_UPPER, communicator, MPI_STATUS_IGNORE);
      MPI_Sendrecv(scalar.plane(1), count, MPI_DOUBLE,
                   lower, SCALAR_LOWER,
                   scalar.plane(scalar.ownedZ() + 1), count, MPI_DOUBLE,
                   upper, SCALAR_LOWER, communicator, MPI_STATUS_IGNORE);
      if (lower == MPI_PROC_NULL)
        std::fill(scalar.plane(0),
                  scalar.plane(0) + scalar.planeValues(), 0.0);
      if (upper == MPI_PROC_NULL)
        std::fill(scalar.plane(scalar.ownedZ() + 1),
                  scalar.plane(scalar.ownedZ() + 1) +
                    scalar.planeValues(), 0.0);
    }

    void sumSharedChargePlanes(YeeComponent& rho, int rank, int size,
                               MPI_Comm communicator)
    {
      if (size == 1) return;
      const int lower = rank > 0 ? rank - 1 : MPI_PROC_NULL;
      const int upper = rank + 1 < size ? rank + 1 : MPI_PROC_NULL;
      const std::size_t planeValues = rho.nx() * rho.ny();
      const int count = checkedPlaneCount(planeValues);
      std::vector<Double> receiveLower(planeValues, 0.0);
      std::vector<Double> receiveUpper(planeValues, 0.0);
      Double* lowerPlane = rho.data();
      Double* upperPlane = rho.data() + (rho.nz() - 1) * planeValues;

      MPI_Sendrecv(upperPlane, count, MPI_DOUBLE, upper, RHO_UPPER,
                   receiveLower.data(), count, MPI_DOUBLE, lower, RHO_UPPER,
                   communicator, MPI_STATUS_IGNORE);
      MPI_Sendrecv(lowerPlane, count, MPI_DOUBLE, lower, RHO_LOWER,
                   receiveUpper.data(), count, MPI_DOUBLE, upper, RHO_LOWER,
                   communicator, MPI_STATUS_IGNORE);
      if (lower != MPI_PROC_NULL)
        for (std::size_t value = 0; value < planeValues; ++value)
          lowerPlane[value] += receiveLower[value];
      if (upper != MPI_PROC_NULL)
        for (std::size_t value = 0; value < planeValues; ++value)
          upperPlane[value] += receiveUpper[value];
    }

    std::unique_ptr<YeeComponent> depositCharge(
        EBFieldGrid& fields,
        const FieldVector<Double>& localOriginSI,
        const std::vector<RelativisticParticleSI>& particles,
        int rank, int size, MPI_Comm communicator)
    {
      const EBGridGeometry& local = fields.geometry();
      std::unique_ptr<YeeComponent> rho(new YeeComponent(
        local.nx + 1, local.ny + 1, local.nz + 1));
      ChargeConservingCurrentDepositor depositor(fields, localOriginSI);
      for (std::size_t index = 0; index < particles.size(); ++index)
        depositor.depositCharge(particles[index].position,
                                particles[index].charge, *rho);
      sumSharedChargePlanes(*rho, rank, size, communicator);
      return rho;
    }

    bool interiorVertex(std::size_t i, std::size_t j,
                        std::size_t globalK,
                        const EBGridGeometry& global)
    {
      return i > 0 && i < global.nx &&
             j > 0 && j < global.ny &&
             globalK > 0 && globalK < global.nz;
    }

    void applyNegativeLaplacian(
        ScalarSlab& input, ScalarSlab& output,
        const EBGridGeometry& global, std::size_t localZOffset,
        int rank, int size, MPI_Comm communicator)
    {
      exchangeScalarGhosts(input, rank, size, communicator);
      output.fill(0.0);
      const Double cx = 1.0 / (global.dx * global.dx);
      const Double cy = 1.0 / (global.dy * global.dy);
      const Double cz = 1.0 / (global.dz * global.dz);
      const Double diagonal = 2.0 * (cx + cy + cz);
      for (std::size_t p = 1; p <= input.ownedZ(); ++p)
        {
          const std::size_t globalK = localZOffset + p - 1;
          if (globalK == 0 || globalK == global.nz) continue;
          for (std::size_t j = 1; j < global.ny; ++j)
            for (std::size_t i = 1; i < global.nx; ++i)
              output(i, j, p) = diagonal * input(i, j, p) -
                cx * (input(i - 1, j, p) + input(i + 1, j, p)) -
                cy * (input(i, j - 1, p) + input(i, j + 1, p)) -
                cz * (input(i, j, p - 1) + input(i, j, p + 1));
        }
    }

    long double globalDot(
        const ScalarSlab& left, const ScalarSlab& right,
        const EBGridGeometry& global, std::size_t localZOffset,
        MPI_Comm communicator)
    {
      long double local = 0.0L;
      for (std::size_t p = 1; p <= left.ownedZ(); ++p)
        {
          const std::size_t globalK = localZOffset + p - 1;
          if (globalK == 0 || globalK == global.nz) continue;
          for (std::size_t j = 1; j < global.ny; ++j)
            for (std::size_t i = 1; i < global.nx; ++i)
              local += static_cast<long double>(left(i, j, p)) *
                       static_cast<long double>(right(i, j, p));
        }
      long double result = 0.0L;
      MPI_Allreduce(&local, &result, 1, MPI_LONG_DOUBLE, MPI_SUM,
                    communicator);
      return result;
    }

    void updateInterior(ScalarSlab& target, const ScalarSlab& source,
                        Double scale, const EBGridGeometry& global,
                        std::size_t localZOffset)
    {
      for (std::size_t p = 1; p <= target.ownedZ(); ++p)
        {
          const std::size_t globalK = localZOffset + p - 1;
          if (globalK == 0 || globalK == global.nz) continue;
          for (std::size_t j = 1; j < global.ny; ++j)
            for (std::size_t i = 1; i < global.nx; ++i)
              target(i, j, p) += scale * source(i, j, p);
        }
    }

    void combineDirection(ScalarSlab& direction,
                          const ScalarSlab& residual, Double beta,
                          const EBGridGeometry& global,
                          std::size_t localZOffset)
    {
      for (std::size_t p = 1; p <= direction.ownedZ(); ++p)
        {
          const std::size_t globalK = localZOffset + p - 1;
          for (std::size_t j = 0; j <= global.ny; ++j)
            for (std::size_t i = 0; i <= global.nx; ++i)
              direction(i, j, p) =
                interiorVertex(i, j, globalK, global) ?
                residual(i, j, p) + beta * direction(i, j, p) : 0.0;
        }
    }

    void electricFromPotential(EBFieldGrid& fields, ScalarSlab& potential,
                               int rank, int size, MPI_Comm communicator)
    {
      exchangeScalarGhosts(potential, rank, size, communicator);
      const EBGridGeometry& g = fields.geometry();
      for (std::size_t k = 0; k <= g.nz; ++k)
        for (std::size_t j = 0; j <= g.ny; ++j)
          for (std::size_t i = 0; i < g.nx; ++i)
            fields.ex()(i, j, k) =
              -(potential(i + 1, j, k + 1) -
                potential(i, j, k + 1)) / g.dx;
      for (std::size_t k = 0; k <= g.nz; ++k)
        for (std::size_t j = 0; j < g.ny; ++j)
          for (std::size_t i = 0; i <= g.nx; ++i)
            fields.ey()(i, j, k) =
              -(potential(i, j + 1, k + 1) -
                potential(i, j, k + 1)) / g.dy;
      for (std::size_t k = 0; k < g.nz; ++k)
        for (std::size_t j = 0; j <= g.ny; ++j)
          for (std::size_t i = 0; i <= g.nx; ++i)
            fields.ez()(i, j, k) =
              -(potential(i, j, k + 2) -
                potential(i, j, k + 1)) / g.dz;
    }

    void validateGauss(
        const EBFieldGrid& fields, const YeeComponent& rho,
        const EBGridGeometry& global, std::size_t localZOffset,
        MPI_Comm communicator, Double& relative, Double& maximum)
    {
      const EBGridGeometry& local = fields.geometry();
      int rank = 0;
      int size = 1;
      communicatorRanks(communicator, rank, size);
      const int lower = rank > 0 ? rank - 1 : MPI_PROC_NULL;
      const int upper = rank + 1 < size ? rank + 1 : MPI_PROC_NULL;
      const std::size_t ezPlaneValues = fields.ez().nx() * fields.ez().ny();
      const int ezCount = checkedPlaneCount(ezPlaneValues);
      std::vector<Double> lowerEz(ezPlaneValues, 0.0);
      const Double* upperEz = fields.ez().data() +
        (local.nz - 1) * ezPlaneValues;
      MPI_Sendrecv(upperEz, ezCount, MPI_DOUBLE, upper, GAUSS_EZ_UPPER,
                   lowerEz.data(), ezCount, MPI_DOUBLE, lower,
                   GAUSS_EZ_UPPER, communicator, MPI_STATUS_IGNORE);
      long double localResidual2 = 0.0L;
      long double localSource2 = 0.0L;
      Double localMaximum = 0.0;
      for (std::size_t k = 0; k < local.nz; ++k)
        {
          const std::size_t globalK = localZOffset + k;
          if (globalK == 0 || globalK == global.nz) continue;
          for (std::size_t j = 1; j < global.ny; ++j)
            for (std::size_t i = 1; i < global.nx; ++i)
              {
                const Double leftZ = k > 0 ?
                  fields.ez()(i, j, k - 1) :
                  lowerEz[j * fields.ez().nx() + i];
                const Double divergence =
                  (fields.ex()(i, j, k) -
                   fields.ex()(i - 1, j, k)) / local.dx +
                  (fields.ey()(i, j, k) -
                   fields.ey()(i, j - 1, k)) / local.dy +
                  (fields.ez()(i, j, k) -
                   leftZ) / local.dz;
                const Double source = rho(i, j, k) / SI::epsilon0;
                const Double error = divergence - source;
                localResidual2 += static_cast<long double>(error) * error;
                localSource2 += static_cast<long double>(source) * source;
                localMaximum = std::max(localMaximum, std::abs(error));
              }
        }
      long double residual2 = 0.0L;
      long double source2 = 0.0L;
      MPI_Allreduce(&localResidual2, &residual2, 1, MPI_LONG_DOUBLE,
                    MPI_SUM, communicator);
      MPI_Allreduce(&localSource2, &source2, 1, MPI_LONG_DOUBLE,
                    MPI_SUM, communicator);
      MPI_Allreduce(&localMaximum, &maximum, 1, MPI_DOUBLE, MPI_MAX,
                    communicator);
      relative = source2 > 0.0L ?
        static_cast<Double>(std::sqrt(residual2 / source2)) : 0.0;
    }
  }

  EBGaussInitializationReport::EBGaussInitializationReport()
    : iterations(0), relativeResidual(0.0), maximumGaussResidual(0.0),
      totalCharge(0.0), temporaryBytes(0)
  {}

  EBGaussInitializationReport EBGaussFieldInitializer::initialize(
      EBFieldGrid& fields,
      const EBGridGeometry& globalGeometry,
      std::size_t localZOffset,
      const FieldVector<Double>& localOriginSI,
      const std::vector<RelativisticParticleSI>& particles,
      Double relativeTolerance,
      std::size_t maximumIterations,
      MPI_Comm communicator)
  {
    if (!(relativeTolerance > 0.0) || !(relativeTolerance < 1.0) ||
        !std::isfinite(relativeTolerance) || maximumIterations == 0)
      throw std::invalid_argument(
        "Invalid initial Gauss-field convergence controls");
    const EBGridGeometry& local = fields.geometry();
    if (local.nx != globalGeometry.nx || local.ny != globalGeometry.ny ||
        localZOffset + local.nz > globalGeometry.nz)
      throw std::invalid_argument(
        "Initial Gauss-field local slab does not match global geometry");

    int rank = 0;
    int size = 1;
    communicatorRanks(communicator, rank, size);
    fields.clearFields();

    ScalarSlab potential(globalGeometry.nx + 1, globalGeometry.ny + 1,
                         local.nz);
    ScalarSlab residual(globalGeometry.nx + 1, globalGeometry.ny + 1,
                        local.nz);
    ScalarSlab direction(globalGeometry.nx + 1, globalGeometry.ny + 1,
                         local.nz);
    ScalarSlab image(globalGeometry.nx + 1, globalGeometry.ny + 1,
                     local.nz);

    EBGaussInitializationReport report;
    report.temporaryBytes = potential.bytes() + residual.bytes() +
      direction.bytes() + image.bytes();

    std::unique_ptr<YeeComponent> rho = depositCharge(
      fields, localOriginSI, particles, rank, size, communicator);
    long double localCharge = 0.0L;
    for (std::size_t index = 0; index < particles.size(); ++index)
      localCharge += static_cast<long double>(particles[index].charge);
    long double totalCharge = 0.0L;
    MPI_Allreduce(&localCharge, &totalCharge, 1, MPI_LONG_DOUBLE, MPI_SUM,
                  communicator);
    report.totalCharge = static_cast<Double>(totalCharge);

    long double localInteriorCharge = 0.0L;
    const long double vertexVolume =
      static_cast<long double>(globalGeometry.dx) *
      static_cast<long double>(globalGeometry.dy) *
      static_cast<long double>(globalGeometry.dz);
    for (std::size_t p = 1; p <= local.nz; ++p)
      {
        const std::size_t globalK = localZOffset + p - 1;
        for (std::size_t j = 0; j <= globalGeometry.ny; ++j)
          for (std::size_t i = 0; i <= globalGeometry.nx; ++i)
            if (interiorVertex(i, j, globalK, globalGeometry))
              {
                residual(i, j, p) =
                  (*rho)(i, j, p - 1) / SI::epsilon0;
                localInteriorCharge +=
                  static_cast<long double>((*rho)(i, j, p - 1)) *
                  vertexVolume;
              }
      }
    long double interiorCharge = 0.0L;
    MPI_Allreduce(&localInteriorCharge, &interiorCharge, 1,
                  MPI_LONG_DOUBLE, MPI_SUM, communicator);
    if (std::abs((interiorCharge - totalCharge) / totalCharge) > 1.0e-12L)
      throw std::runtime_error(
        "Initial CIC charge cloud touches the zero-potential boundary; "
        "move every particle at least one complete cell inward or enlarge "
        "the mesh/CPML padding");
    direction = residual;
    rho.reset();

    const long double sourceNorm2 = globalDot(
      residual, residual, globalGeometry, localZOffset, communicator);
    if (!(sourceNorm2 > 0.0L) ||
        !std::isfinite(static_cast<Double>(sourceNorm2)))
      throw std::runtime_error(
        "Initial particle charge produced an invalid Poisson source norm");

    long double residualNorm2 = sourceNorm2;
    for (std::size_t iteration = 0;
         iteration < maximumIterations; ++iteration)
      {
        applyNegativeLaplacian(direction, image, globalGeometry,
          localZOffset, rank, size, communicator);
        const long double denominator = globalDot(
          direction, image, globalGeometry, localZOffset, communicator);
        if (!(denominator > 0.0L) ||
            !std::isfinite(static_cast<Double>(denominator)))
          throw std::runtime_error(
            "Initial Gauss-field CG lost positive definiteness");
        const Double alpha = static_cast<Double>(
          residualNorm2 / denominator);
        updateInterior(potential, direction, alpha, globalGeometry,
                       localZOffset);
        updateInterior(residual, image, -alpha, globalGeometry,
                       localZOffset);
        const long double nextNorm2 = globalDot(
          residual, residual, globalGeometry, localZOffset, communicator);
        report.iterations = iteration + 1;
        report.relativeResidual = static_cast<Double>(
          std::sqrt(nextNorm2 / sourceNorm2));
        if (report.relativeResidual <= relativeTolerance)
          {
            residualNorm2 = nextNorm2;
            break;
          }
        if (!(nextNorm2 >= 0.0L) ||
            !std::isfinite(static_cast<Double>(nextNorm2)))
          throw std::runtime_error(
            "Initial Gauss-field CG produced a non-finite residual");
        const Double beta = static_cast<Double>(
          nextNorm2 / residualNorm2);
        combineDirection(direction, residual, beta, globalGeometry,
                         localZOffset);
        residualNorm2 = nextNorm2;
      }

    if (report.relativeResidual > relativeTolerance)
      throw std::runtime_error(
        "Initial Gauss-field solve did not converge; increase "
        "initial_self_field.maximum_iterations, relax its relative_tolerance, "
        "or reduce the grid aspect ratio");

    electricFromPotential(fields, potential, rank, size, communicator);
    std::unique_ptr<YeeComponent> validationRho = depositCharge(
      fields, localOriginSI, particles, rank, size, communicator);
    validateGauss(fields, *validationRho, globalGeometry, localZOffset,
                  communicator, report.relativeResidual,
                  report.maximumGaussResidual);
    if (report.relativeResidual > 4.0 * relativeTolerance)
      throw std::runtime_error(
        "Initial E field failed the post-solve discrete Gauss check");
    return report;
  }
}
