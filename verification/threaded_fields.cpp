#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

#include "compute_threads.h"
#include "eb_cpml.h"
#include "eb_field.h"

namespace
{
  void require(bool condition, const std::string& message)
  {
    if (!condition) throw std::runtime_error(message);
  }

  std::vector<aprl::YeeComponent*> components(aprl::EBFieldGrid& fields)
  {
    return {&fields.ex(), &fields.ey(), &fields.ez(),
            &fields.bx(), &fields.by(), &fields.bz(),
            &fields.jx(), &fields.jy(), &fields.jz()};
  }

  void seed(aprl::EBFieldGrid& fields)
  {
    const std::vector<aprl::YeeComponent*> all = components(fields);
    for (std::size_t component = 0; component < all.size(); ++component)
      for (std::size_t index = 0; index < all[component]->size(); ++index)
        {
          const double scale = component >= 3 && component < 6 ?
            1.0 / aprl::SI::c : 1.0;
          all[component]->data()[index] = scale *
            (std::sin(0.017 * static_cast<double>(index + 3 * component)) +
             0.2 * std::cos(0.031 * static_cast<double>(index + component)));
        }
  }

  void equalFields(aprl::EBFieldGrid& left, aprl::EBFieldGrid& right,
                   const std::string& phase)
  {
    const std::vector<aprl::YeeComponent*> a = components(left);
    const std::vector<aprl::YeeComponent*> b = components(right);
    for (std::size_t component = 0; component < a.size(); ++component)
      {
        require(a[component]->size() == b[component]->size(),
          "Mismatched component sizes");
        for (std::size_t index = 0; index < a[component]->size(); ++index)
          require(std::isfinite(a[component]->data()[index]),
            "Non-finite serial reference field");
        require(std::memcmp(a[component]->data(), b[component]->data(),
                  a[component]->size() * sizeof(double)) == 0,
          phase + ": threaded field differs bitwise from serial component " +
          std::to_string(component));
      }
  }

  // Preserve the general PEC/halo sampler as an independent reference for
  // the optimized interior stencil.  The reference always checks boundaries.
  double referenceElectric(const aprl::YeeComponent& component,
                           unsigned int componentAxis,
                           std::ptrdiff_t i, std::ptrdiff_t j,
                           std::ptrdiff_t k,
                           const aprl::EBElectricHaloView& halo)
  {
    std::ptrdiff_t coordinate[3] = {i, j, k};
    const std::ptrdiff_t size[3] = {
      static_cast<std::ptrdiff_t>(component.nx()),
      static_cast<std::ptrdiff_t>(component.ny()),
      static_cast<std::ptrdiff_t>(component.nz())
    };
    double sign = 1.0;
    for (unsigned int axis = 0; axis < 2; ++axis)
      {
        if (coordinate[axis] < 0)
          {
            if (componentAxis == axis)
              coordinate[axis] = 0;
            else
              {
                coordinate[axis] = -coordinate[axis];
                sign = -sign;
              }
          }
        else if (coordinate[axis] >= size[axis])
          {
            if (componentAxis == axis)
              coordinate[axis] = size[axis] - 1;
            else
              {
                coordinate[axis] = 2 * (size[axis] - 1) - coordinate[axis];
                sign = -sign;
              }
          }
      }
    if (coordinate[2] < 0 || coordinate[2] >= size[2])
      {
        const bool lower = coordinate[2] < 0;
        const double* plane = componentAxis == 0 ?
          (lower ? halo.lowerEx : halo.upperEx) :
          (componentAxis == 1 ? (lower ? halo.lowerEy : halo.upperEy) :
           (lower ? halo.lowerEz : halo.upperEz));
        if (plane)
          return sign * plane[
            static_cast<std::size_t>(coordinate[1]) * component.nx() +
            static_cast<std::size_t>(coordinate[0])];
        if (componentAxis == 2)
          coordinate[2] = lower ? 0 : size[2] - 1;
        else
          {
            coordinate[2] = lower ? -coordinate[2] :
              2 * (size[2] - 1) - coordinate[2];
            sign = -sign;
          }
      }
    return sign * component(static_cast<std::size_t>(coordinate[0]),
                            static_cast<std::size_t>(coordinate[1]),
                            static_cast<std::size_t>(coordinate[2]));
  }

  double referenceSmoothed(const aprl::EBCowanCoefficients& coefficients,
                           const aprl::YeeComponent& component,
                           unsigned int componentAxis,
                           unsigned int derivativeAxis,
                           std::ptrdiff_t i, std::ptrdiff_t j,
                           std::ptrdiff_t k,
                           const aprl::EBElectricHaloView& halo)
  {
    const unsigned int transverseA = (derivativeAxis + 1) % 3;
    const unsigned int transverseB = (derivativeAxis + 2) % 3;
    const double delta = derivativeAxis == 0 ? coefficients.deltaYZ :
      (derivativeAxis == 1 ? coefficients.deltaZX : coefficients.deltaXY);
    const std::ptrdiff_t base[3] = {i, j, k};
    double result = coefficients.alpha[derivativeAxis] *
      referenceElectric(component, componentAxis, i, j, k, halo);
    for (int direction = -1; direction <= 1; direction += 2)
      {
        std::ptrdiff_t pointA[3] = {base[0], base[1], base[2]};
        std::ptrdiff_t pointB[3] = {base[0], base[1], base[2]};
        pointA[transverseA] += direction;
        pointB[transverseB] += direction;
        result += coefficients.beta[transverseA] * referenceElectric(
          component, componentAxis, pointA[0], pointA[1], pointA[2], halo);
        result += coefficients.beta[transverseB] * referenceElectric(
          component, componentAxis, pointB[0], pointB[1], pointB[2], halo);
      }
    for (int directionA = -1; directionA <= 1; directionA += 2)
      for (int directionB = -1; directionB <= 1; directionB += 2)
        {
          std::ptrdiff_t point[3] = {base[0], base[1], base[2]};
          point[transverseA] += directionA;
          point[transverseB] += directionB;
          result += delta * referenceElectric(
            component, componentAxis, point[0], point[1], point[2], halo);
        }
    return result;
  }

  void referenceMagnetic(aprl::EBFieldGrid& fields,
                         const aprl::EBElectricHaloView& halo)
  {
    const aprl::EBGridGeometry& geometry = fields.geometry();
    const aprl::EBCowanCoefficients& coefficients = fields.cowanCoefficients();
    const double dtdx = geometry.dt / geometry.dx;
    const double dtdy = geometry.dt / geometry.dy;
    const double dtdz = geometry.dt / geometry.dz;
    for (std::size_t k = 0; k < geometry.nz; ++k)
      for (std::size_t j = 0; j < geometry.ny; ++j)
        for (std::size_t i = 0; i <= geometry.nx; ++i)
          {
            const double dEzDy = referenceSmoothed(
              coefficients, fields.ez(), 2, 1, i, j + 1, k, halo) -
              referenceSmoothed(coefficients, fields.ez(), 2, 1,
                                i, j, k, halo);
            const double dEyDz = referenceSmoothed(
              coefficients, fields.ey(), 1, 2, i, j, k + 1, halo) -
              referenceSmoothed(coefficients, fields.ey(), 1, 2,
                                i, j, k, halo);
            fields.bx()(i, j, k) -= dtdy * dEzDy - dtdz * dEyDz;
          }
    for (std::size_t k = 0; k < geometry.nz; ++k)
      for (std::size_t j = 0; j <= geometry.ny; ++j)
        for (std::size_t i = 0; i < geometry.nx; ++i)
          {
            const double dExDz = referenceSmoothed(
              coefficients, fields.ex(), 0, 2, i, j, k + 1, halo) -
              referenceSmoothed(coefficients, fields.ex(), 0, 2,
                                i, j, k, halo);
            const double dEzDx = referenceSmoothed(
              coefficients, fields.ez(), 2, 0, i + 1, j, k, halo) -
              referenceSmoothed(coefficients, fields.ez(), 2, 0,
                                i, j, k, halo);
            fields.by()(i, j, k) -= dtdz * dExDz - dtdx * dEzDx;
          }
    for (std::size_t k = 0; k <= geometry.nz; ++k)
      for (std::size_t j = 0; j < geometry.ny; ++j)
        for (std::size_t i = 0; i < geometry.nx; ++i)
          {
            const double dEyDx = referenceSmoothed(
              coefficients, fields.ey(), 1, 0, i + 1, j, k, halo) -
              referenceSmoothed(coefficients, fields.ey(), 1, 0,
                                i, j, k, halo);
            const double dExDy = referenceSmoothed(
              coefficients, fields.ex(), 0, 1, i, j + 1, k, halo) -
              referenceSmoothed(coefficients, fields.ex(), 0, 1,
                                i, j, k, halo);
            fields.bz()(i, j, k) -= dtdx * dEyDx - dtdy * dExDy;
          }
  }

  void compareCowanReference(int provided, int threads,
                             std::size_t nx, std::size_t ny, std::size_t nz,
                             unsigned int haloMask)
  {
    const aprl::EBGridGeometry geometry(
      nx, ny, nz, 0.12, 0.10, 0.08, 0.08 / aprl::SI::c,
      aprl::EBMaxwellSolver::CowanZ);
    aprl::EBFieldGrid expected(geometry);
    aprl::EBFieldGrid actual(geometry);
    seed(expected);
    seed(actual);
    const aprl::YeeComponent* electric[3] = {
      &expected.ex(), &expected.ey(), &expected.ez()
    };
    std::vector<double> lower[3];
    std::vector<double> upper[3];
    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        const std::size_t size = electric[axis]->nx() * electric[axis]->ny();
        lower[axis].resize(size);
        upper[axis].resize(size);
        for (std::size_t index = 0; index < size; ++index)
          {
            lower[axis][index] = 0.3 + std::sin(
              0.013 * static_cast<double>(index + 11 * axis));
            upper[axis][index] = -0.7 + std::cos(
              0.021 * static_cast<double>(index + 17 * axis));
          }
      }
    aprl::EBElectricHaloView halo;
    if (haloMask & 1)
      {
        halo.lowerEx = lower[0].data();
        halo.lowerEy = lower[1].data();
        halo.lowerEz = lower[2].data();
      }
    if (haloMask & 2)
      {
        halo.upperEx = upper[0].data();
        halo.upperEy = upper[1].data();
        halo.upperEz = upper[2].data();
      }
    referenceMagnetic(expected, halo);
    aprl::configureComputeThreads(provided, threads);
    actual.advanceMagnetic(halo);
    equalFields(expected, actual, "general Cowan reference " +
      std::to_string(nx) + "x" + std::to_string(ny) + "x" +
      std::to_string(nz) + " halo=" + std::to_string(haloMask) +
      " threads=" + std::to_string(threads));
  }

  void compare(int provided, int threads, aprl::EBMaxwellSolver solver,
               bool cpml, bool useHalo)
  {
    const aprl::EBGridGeometry geometry(
      17, 15, 12, 0.12, 0.10, 0.08,
      (solver == aprl::EBMaxwellSolver::CowanZ ? 0.08 : 0.02) / aprl::SI::c,
      solver);
    aprl::EBFieldGrid serial(geometry);
    aprl::EBFieldGrid parallel(geometry);
    seed(serial);
    seed(parallel);
    if (cpml)
      {
        aprl::EBCPMLParameters parameters;
        parameters.cells[0] = parameters.cells[1] = 2;
        // Synthetic interior halos cover the Cowan MPI-face read path;
        // z CPML is used only when the tested faces are physical.
        parameters.cells[2] = useHalo ? 0 : 2;
        serial.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(
          new aprl::EBConvolutionalPML(geometry, geometry.nz, 0, parameters,
                                      !useHalo, !useHalo)));
        parallel.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(
          new aprl::EBConvolutionalPML(geometry, geometry.nz, 0, parameters,
                                      !useHalo, !useHalo)));
      }
    std::vector<double> haloEx(serial.ex().nx() * serial.ex().ny());
    std::vector<double> haloEy(serial.ey().nx() * serial.ey().ny());
    std::vector<double> haloEz(serial.ez().nx() * serial.ez().ny());
    for (std::size_t index = 0; index < haloEx.size(); ++index)
      haloEx[index] = std::sin(0.02 * static_cast<double>(index));
    for (std::size_t index = 0; index < haloEy.size(); ++index)
      haloEy[index] = std::cos(0.03 * static_cast<double>(index));
    for (std::size_t index = 0; index < haloEz.size(); ++index)
      haloEz[index] = std::sin(0.04 * static_cast<double>(index));
    aprl::EBElectricHaloView halo;
    if (useHalo)
      {
        halo.lowerEx = halo.upperEx = haloEx.data();
        halo.lowerEy = halo.upperEy = haloEy.data();
        halo.lowerEz = halo.upperEz = haloEz.data();
      }
    for (unsigned int step = 0; step < 8; ++step)
      {
        aprl::configureComputeThreads(provided, 1);
        serial.advanceMagnetic(halo);
        aprl::configureComputeThreads(provided, threads);
        parallel.advanceMagnetic(halo);
        equalFields(serial, parallel, "magnetic update");

        aprl::configureComputeThreads(provided, 1);
        serial.advanceElectric();
        aprl::configureComputeThreads(provided, threads);
        parallel.advanceElectric();
        equalFields(serial, parallel, "electric update");

        aprl::configureComputeThreads(provided, 1);
        serial.clearCurrent();
        aprl::configureComputeThreads(provided, threads);
        parallel.clearCurrent();
        equalFields(serial, parallel, "current clearing");
      }
  }
}

int main(int argc, char** argv)
{
  int provided = MPI_THREAD_SINGLE;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  int status = 0;
  try
    {
      if (argc == 2 && std::string(argv[1]) == "--environment")
        {
          aprl::configureComputeThreads(provided);
          std::cout << "compute_threads=" << aprl::computeThreadCount()
                    << std::endl;
        }
      else
        {
          require(aprl::computeThreadCount() == 1,
            "Compute threads must default to one before configuration");
          bool insufficientRejected = false;
          try { aprl::configureComputeThreads(MPI_THREAD_SINGLE, 2); }
          catch (const std::exception&) { insufficientRejected = true; }
          require(insufficientRejected,
            "Insufficient MPI thread support was silently accepted");

          if (!aprl::computeThreadingAvailable())
            {
              bool unavailableRejected = false;
              try { aprl::configureComputeThreads(provided, 2); }
              catch (const std::exception&) { unavailableRejected = true; }
              require(unavailableRejected,
                "A no-OpenMP build silently accepted multiple threads");
            }
          const int maximum = aprl::computeThreadingAvailable() ? 4 : 1;
          for (int threads = 1; threads <= maximum; threads *= 2)
            {
              for (unsigned int haloMask = 0; haloMask < 4; ++haloMask)
                {
                  // The minimum grid has no complete interior stencil;
                  // larger irregular grids exercise both interior and edges.
                  compareCowanReference(provided, threads, 2, 2, 2, haloMask);
                  compareCowanReference(provided, threads, 3, 4, 5, haloMask);
                  compareCowanReference(provided, threads, 17, 15, 2,
                                        haloMask);
                  compareCowanReference(provided, threads, 17, 15, 12,
                                        haloMask);
                }
              for (unsigned int solver = 0; solver < 2; ++solver)
                for (unsigned int boundary = 0; boundary < 2; ++boundary)
                  {
                    compare(provided, threads, solver == 0 ?
                      aprl::EBMaxwellSolver::Yee : aprl::EBMaxwellSolver::CowanZ,
                      boundary != 0, false);
                    if (solver == 1)
                      compare(provided, threads, aprl::EBMaxwellSolver::CowanZ,
                        boundary != 0, true);
                  }
            }
#ifdef DIRECT_EB_BOUNDS_CHECK
          aprl::YeeComponent checked(2, 3, 4);
          bool boundsRejected = false;
          try { checked(2, 0, 0) = 1.0; }
          catch (const std::out_of_range&) { boundsRejected = true; }
          require(boundsRejected, "Inline accessor lost its bounds check");
#endif
          std::cout << "serial_threaded_fields=bitwise_identical "
                       "cowan_general_reference=bitwise_identical "
                       "yee_cowan_cpml_halo=passed max_threads=" << maximum
                    << std::endl;
        }
    }
  catch (const std::exception& error)
    {
      std::cerr << error.what() << std::endl;
      status = 1;
    }
  MPI_Finalize();
  return status;
}
