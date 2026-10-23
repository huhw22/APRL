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
#ifdef DIRECT_EB_BOUNDS_CHECK
          aprl::YeeComponent checked(2, 3, 4);
          bool boundsRejected = false;
          try { checked(2, 0, 0) = 1.0; }
          catch (const std::out_of_range&) { boundsRejected = true; }
          require(boundsRejected, "Inline accessor lost its bounds check");
#endif
          std::cout << "serial_threaded_fields=bitwise_identical "
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
