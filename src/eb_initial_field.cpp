#include "eb_initial_field.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

#include "eb_deposition.h"

namespace aprl
{
  namespace
  {
    const int RHO_UPPER = 910;
    const int RHO_LOWER = 911;
    const int SCALAR_UPPER = 912;
    const int SCALAR_LOWER = 913;
    const int GAUSS_EZ_UPPER = 914;

    struct ScalarDomain
    {
      std::size_t lower[3];
      std::size_t upper[3];
    };

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
                        const ScalarDomain& domain)
    {
      return i > domain.lower[0] && i < domain.upper[0] &&
             j > domain.lower[1] && j < domain.upper[1] &&
             globalK > domain.lower[2] && globalK < domain.upper[2];
    }

    void applyNegativeLaplacian(
        ScalarSlab& input, ScalarSlab& output,
        const EBGridGeometry& global, const ScalarDomain& domain,
        std::size_t localZOffset,
        Double longitudinalCoefficient,
        int rank, int size, MPI_Comm communicator)
    {
      exchangeScalarGhosts(input, rank, size, communicator);
      output.fill(0.0);
      const Double cx = 1.0 / (global.dx * global.dx);
      const Double cy = 1.0 / (global.dy * global.dy);
      const Double cz = longitudinalCoefficient /
        (global.dz * global.dz);
      const Double diagonal = 2.0 * (cx + cy + cz);
      for (std::size_t p = 1; p <= input.ownedZ(); ++p)
        {
          const std::size_t globalK = localZOffset + p - 1;
          if (globalK <= domain.lower[2] ||
              globalK >= domain.upper[2]) continue;
          for (std::size_t j = domain.lower[1] + 1;
               j < domain.upper[1]; ++j)
            for (std::size_t i = domain.lower[0] + 1;
                 i < domain.upper[0]; ++i)
              output(i, j, p) = diagonal * input(i, j, p) -
                cx * (input(i - 1, j, p) + input(i + 1, j, p)) -
                cy * (input(i, j - 1, p) + input(i, j + 1, p)) -
                cz * (input(i, j, p - 1) + input(i, j, p + 1));
        }
    }

    long double globalDot(
        const ScalarSlab& left, const ScalarSlab& right,
        const ScalarDomain& domain, std::size_t localZOffset,
        MPI_Comm communicator)
    {
      long double local = 0.0L;
      for (std::size_t p = 1; p <= left.ownedZ(); ++p)
        {
          const std::size_t globalK = localZOffset + p - 1;
          if (globalK <= domain.lower[2] ||
              globalK >= domain.upper[2]) continue;
          for (std::size_t j = domain.lower[1] + 1;
               j < domain.upper[1]; ++j)
            for (std::size_t i = domain.lower[0] + 1;
                 i < domain.upper[0]; ++i)
              local += static_cast<long double>(left(i, j, p)) *
                       static_cast<long double>(right(i, j, p));
        }
      long double result = 0.0L;
      MPI_Allreduce(&local, &result, 1, MPI_LONG_DOUBLE, MPI_SUM,
                    communicator);
      return result;
    }

    void updateInterior(ScalarSlab& target, const ScalarSlab& source,
                        Double scale, const ScalarDomain& domain,
                        std::size_t localZOffset)
    {
      for (std::size_t p = 1; p <= target.ownedZ(); ++p)
        {
          const std::size_t globalK = localZOffset + p - 1;
          if (globalK <= domain.lower[2] ||
              globalK >= domain.upper[2]) continue;
          for (std::size_t j = domain.lower[1] + 1;
               j < domain.upper[1]; ++j)
            for (std::size_t i = domain.lower[0] + 1;
                 i < domain.upper[0]; ++i)
              target(i, j, p) += scale * source(i, j, p);
        }
    }

    void combineDirection(ScalarSlab& direction,
                          const ScalarSlab& residual, Double beta,
                          const EBGridGeometry& global,
                          const ScalarDomain& domain,
                          std::size_t localZOffset)
    {
      for (std::size_t p = 1; p <= direction.ownedZ(); ++p)
        {
          const std::size_t globalK = localZOffset + p - 1;
          for (std::size_t j = 0; j <= global.ny; ++j)
            for (std::size_t i = 0; i <= global.nx; ++i)
              direction(i, j, p) =
                interiorVertex(i, j, globalK, domain) ?
                residual(i, j, p) + beta * direction(i, j, p) : 0.0;
        }
    }

    void fieldsFromPotential(EBFieldGrid& fields, ScalarSlab& potential,
                             Double meanBetaZ,
                             Double longitudinalCoefficient,
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
              -longitudinalCoefficient *
              (potential(i, j, k + 2) -
               potential(i, j, k + 1)) / g.dz;

      /* For a charge distribution translating rigidly along z, the Vay
       * relativistic-Poisson construction gives B=beta x E/c.  B is stored
       * at the leapfrog time -dt/2 while E is at time zero. Rigid translation
       * maps the earlier B face to E(z+beta*c*dt/2, t=0); linear z
       * interpolation therefore also supplies the required half-time
       * staggering without allocating another field copy. */
      const Double magneticScale = meanBetaZ / SI::c;
      const Double upperWeight = 0.5 + 0.5 * meanBetaZ * SI::c *
        g.dt / g.dz;
      const Double lowerWeight = 1.0 - upperWeight;
      if (upperWeight < 0.0 || upperWeight > 1.0)
        throw std::runtime_error(
          "Initial B half-step shift exceeds one longitudinal cell");
      for (std::size_t k = 0; k < g.nz; ++k)
        for (std::size_t j = 0; j < g.ny; ++j)
          for (std::size_t i = 0; i <= g.nx; ++i)
            fields.bx()(i, j, k) = -magneticScale *
              (lowerWeight * fields.ey()(i, j, k) +
               upperWeight * fields.ey()(i, j, k + 1));
      for (std::size_t k = 0; k < g.nz; ++k)
        for (std::size_t j = 0; j <= g.ny; ++j)
          for (std::size_t i = 0; i < g.nx; ++i)
            fields.by()(i, j, k) = magneticScale *
              (lowerWeight * fields.ex()(i, j, k) +
               upperWeight * fields.ex()(i, j, k + 1));
      fields.bz().fill(0.0);
    }

    void meanVelocity(const std::vector<RelativisticParticleSI>& particles,
                      MPI_Comm communicator, Double& betaZ,
                      Double& betaTransverse)
    {
      long double local[4] = {};
      for (std::size_t index = 0; index < particles.size(); ++index)
        {
          const RelativisticParticleSI& particle = particles[index];
          const Double gamma = BoostFrameTransform::gammaFromProperVelocity(
            particle.properVelocity);
          const long double weight = static_cast<long double>(particle.mass);
          local[0] += weight;
          for (unsigned int axis = 0; axis < 3; ++axis)
            local[axis + 1] += weight * static_cast<long double>(
              particle.properVelocity[axis] / gamma);
        }
      long double global[4] = {};
      MPI_Allreduce(local, global, 4, MPI_LONG_DOUBLE, MPI_SUM,
                    communicator);
      if (!(global[0] > 0.0L))
        throw std::runtime_error(
          "Initial relativistic-Poisson field requires positive represented mass");
      const Double betaX = static_cast<Double>(global[1] / global[0]);
      const Double betaY = static_cast<Double>(global[2] / global[0]);
      betaZ = static_cast<Double>(global[3] / global[0]);
      betaTransverse = std::hypot(betaX, betaY);
      if (!std::isfinite(betaZ) || !std::isfinite(betaTransverse) ||
          betaTransverse >= 1.0 || std::abs(betaZ) >= 1.0)
        throw std::runtime_error(
          "Initial particle mean velocity is not a finite subluminal vector");
    }

    void initialFieldEnergy(const EBFieldGrid& fields,
                            MPI_Comm communicator,
                            Double& electric, Double& magnetic)
    {
      const EBGridGeometry& g = fields.geometry();
      long double localElectric = 0.0L;
      long double localMagnetic = 0.0L;
      const long double volume = static_cast<long double>(g.dx) *
        static_cast<long double>(g.dy) * static_cast<long double>(g.dz);
      for (std::size_t k = 0; k < g.nz; ++k)
        for (std::size_t j = 0; j < g.ny; ++j)
          for (std::size_t i = 0; i < g.nx; ++i)
            {
              const RadiationFieldSample sample =
                fields.radiationSampleCell(i, j, k);
              localElectric += 0.5L *
                static_cast<long double>(SI::epsilon0) *
                static_cast<long double>(sample.electric.norm2()) * volume;
              localMagnetic += 0.5L /
                static_cast<long double>(SI::mu0) *
                static_cast<long double>(sample.magnetic.norm2()) * volume;
            }
      long double globalElectric = 0.0L;
      long double globalMagnetic = 0.0L;
      MPI_Allreduce(&localElectric, &globalElectric, 1, MPI_LONG_DOUBLE,
                    MPI_SUM, communicator);
      MPI_Allreduce(&localMagnetic, &globalMagnetic, 1, MPI_LONG_DOUBLE,
                    MPI_SUM, communicator);
      electric = static_cast<Double>(globalElectric);
      magnetic = static_cast<Double>(globalMagnetic);
    }

    void distributionResolution(
        const std::vector<RelativisticParticleSI>& particles,
        const EBGridGeometry& global, const ScalarDomain& domain,
        std::size_t localZOffset,
        const FieldVector<Double>& localOriginSI,
        MPI_Comm communicator, EBGaussInitializationReport& report)
    {
      long double local[7] = {};
      for (std::size_t index = 0; index < particles.size(); ++index)
        {
          const long double weight = static_cast<long double>(
            particles[index].mass);
          local[0] += weight;
          for (unsigned int axis = 0; axis < 3; ++axis)
            {
              const long double position = static_cast<long double>(
                particles[index].position[axis]);
              local[axis + 1] += weight * position;
              local[axis + 4] += weight * position * position;
            }
        }
      long double globalMoment[7] = {};
      MPI_Allreduce(local, globalMoment, 7, MPI_LONG_DOUBLE, MPI_SUM,
                    communicator);
      if (!(globalMoment[0] > 0.0L))
        throw std::runtime_error(
          "Initial field cannot diagnose an empty particle distribution");
      const Double spacing[3] = {global.dx, global.dy, global.dz};
      FieldVector<Double> globalOrigin(localOriginSI);
      globalOrigin[2] -= static_cast<Double>(localZOffset) * global.dz;
      for (unsigned int axis = 0; axis < 3; ++axis)
        {
          const long double mean = globalMoment[axis + 1] /
            globalMoment[0];
          const long double variance = std::max(0.0L,
            globalMoment[axis + 4] / globalMoment[0] - mean * mean);
          const Double rms = static_cast<Double>(std::sqrt(variance));
          report.rmsPositionCells[axis] = rms / spacing[axis];
          if (rms > 0.0)
            {
              const Double lower = globalOrigin[axis] +
                static_cast<Double>(domain.lower[axis]) * spacing[axis];
              const Double upper = globalOrigin[axis] +
                static_cast<Double>(domain.upper[axis]) * spacing[axis];
              report.centrePaddingRms[axis] = std::min(
                static_cast<Double>(mean) - lower,
                upper - static_cast<Double>(mean)) / rms;
            }
          else
            report.centrePaddingRms[axis] = 0.0;
        }
    }

    void validateGauss(
        const EBFieldGrid& fields, const YeeComponent& rho,
        const ScalarDomain& domain,
        std::size_t localZOffset,
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
          if (globalK <= domain.lower[2] ||
              globalK >= domain.upper[2]) continue;
          for (std::size_t j = domain.lower[1] + 1;
               j < domain.upper[1]; ++j)
            for (std::size_t i = domain.lower[0] + 1;
                 i < domain.upper[0]; ++i)
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
      totalCharge(0.0), meanBetaZ(0.0), meanBetaTransverse(0.0),
      initialElectricEnergy(0.0), initialMagneticEnergy(0.0),
      temporaryBytes(0)
  {
    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        rmsPositionCells[axis] = 0.0;
        centrePaddingRms[axis] = 0.0;
        staticFieldGuardCells[axis] = 0;
      }
  }

  EBGaussInitializationReport EBGaussFieldInitializer::initialize(
      EBFieldGrid& fields,
      const EBGridGeometry& globalGeometry,
      std::size_t localZOffset,
      const FieldVector<Double>& localOriginSI,
      const std::vector<RelativisticParticleSI>& particles,
      InitialSelfFieldModel model,
      const std::size_t staticFieldGuardCells[3],
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

    ScalarDomain domain;
    const std::size_t globalCells[3] = {
      globalGeometry.nx, globalGeometry.ny, globalGeometry.nz
    };
    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        if (2 * staticFieldGuardCells[axis] + 2 >= globalCells[axis])
          throw std::invalid_argument(
            "Initial static-field guard leaves fewer than two physical cells");
        domain.lower[axis] = staticFieldGuardCells[axis];
        domain.upper[axis] = globalCells[axis] - staticFieldGuardCells[axis];
      }

    int rank = 0;
    int size = 1;
    communicatorRanks(communicator, rank, size);
    fields.clearFields();

    Double meanBetaZ = 0.0;
    Double meanBetaTransverse = 0.0;
    meanVelocity(particles, communicator, meanBetaZ,
                 meanBetaTransverse);
    if (model == InitialSelfFieldModel::RelativisticPoisson &&
        meanBetaTransverse > 1.0e-8)
      throw std::runtime_error(
        "relativistic-poisson initialization currently supports a common "
        "mean velocity along the boost z axis only; rotate/recenter the "
        "input beam or use a closer longitudinal boost frame");
    const Double appliedBetaZ =
      model == InitialSelfFieldModel::RelativisticPoisson ?
      meanBetaZ : 0.0;
    const Double longitudinalCoefficient =
      std::fma(-appliedBetaZ, appliedBetaZ, 1.0);
    if (!(longitudinalCoefficient >
          64.0 * std::numeric_limits<Double>::epsilon()))
      throw std::runtime_error(
        "Initial relativistic-Poisson operator is ill-conditioned because "
        "the bunch remains too relativistic in the simulation frame; choose "
        "boost_gamma closer to the bunch mean gamma");

    ScalarSlab potential(globalGeometry.nx + 1, globalGeometry.ny + 1,
                         local.nz);
    ScalarSlab residual(globalGeometry.nx + 1, globalGeometry.ny + 1,
                        local.nz);
    ScalarSlab direction(globalGeometry.nx + 1, globalGeometry.ny + 1,
                         local.nz);
    ScalarSlab image(globalGeometry.nx + 1, globalGeometry.ny + 1,
                     local.nz);

    EBGaussInitializationReport report;
    report.meanBetaZ = meanBetaZ;
    report.meanBetaTransverse = meanBetaTransverse;
    for (unsigned int axis = 0; axis < 3; ++axis)
      report.staticFieldGuardCells[axis] = staticFieldGuardCells[axis];
    distributionResolution(particles, globalGeometry, domain,
      localZOffset, localOriginSI, communicator, report);
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
            if (interiorVertex(i, j, globalK, domain))
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
        "Initial CIC charge cloud touches the static zero-potential boundary; "
        "move every particle at least one complete cell inside the CPML "
        "entrance or enlarge the physical mesh padding");
    direction = residual;
    rho.reset();

    const long double sourceNorm2 = globalDot(
      residual, residual, domain, localZOffset, communicator);
    if (!(sourceNorm2 > 0.0L) ||
        !std::isfinite(static_cast<Double>(sourceNorm2)))
      throw std::runtime_error(
        "Initial particle charge produced an invalid Poisson source norm");

    long double residualNorm2 = sourceNorm2;
    for (std::size_t iteration = 0;
         iteration < maximumIterations; ++iteration)
      {
        applyNegativeLaplacian(direction, image, globalGeometry, domain,
          localZOffset, longitudinalCoefficient,
          rank, size, communicator);
        const long double denominator = globalDot(
          direction, image, domain, localZOffset, communicator);
        if (!(denominator > 0.0L) ||
            !std::isfinite(static_cast<Double>(denominator)))
          throw std::runtime_error(
            "Initial Gauss-field CG lost positive definiteness");
        const Double alpha = static_cast<Double>(
          residualNorm2 / denominator);
        updateInterior(potential, direction, alpha, domain,
                       localZOffset);
        updateInterior(residual, image, -alpha, domain,
                       localZOffset);
        const long double nextNorm2 = globalDot(
          residual, residual, domain, localZOffset, communicator);
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
        combineDirection(direction, residual, beta, globalGeometry, domain,
                         localZOffset);
        residualNorm2 = nextNorm2;
      }

    if (report.relativeResidual > relativeTolerance)
      throw std::runtime_error(
        "Initial Gauss-field solve did not converge; increase "
        "initial_self_field.maximum_iterations, relax its relative_tolerance, "
        "or reduce the grid aspect ratio");

    fieldsFromPotential(fields, potential, appliedBetaZ,
      longitudinalCoefficient, rank, size, communicator);
    std::unique_ptr<YeeComponent> validationRho = depositCharge(
      fields, localOriginSI, particles, rank, size, communicator);
    validateGauss(fields, *validationRho, domain, localZOffset,
                  communicator, report.relativeResidual,
                  report.maximumGaussResidual);
    if (report.relativeResidual > 4.0 * relativeTolerance)
      throw std::runtime_error(
        "Initial E field failed the post-solve discrete Gauss check");
    initialFieldEnergy(fields, communicator, report.initialElectricEnergy,
                       report.initialMagneticEnergy);
    return report;
  }
}
