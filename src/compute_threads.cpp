#include "compute_threads.h"

#include <cerrno>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include <mpi.h>
#ifdef APRL_HAS_OPENMP
#include <omp.h>
#endif

namespace aprl
{
  namespace
  {
    int configuredThreads = 1;

    int environmentThreadCount()
    {
      const char* value = std::getenv("OMP_NUM_THREADS");
      if (!value) return 1;
      errno = 0;
      char* end = 0;
      const long count = std::strtol(value, &end, 10);
      const bool hadDigits = end != value;
      while (end && std::isspace(static_cast<unsigned char>(*end))) ++end;
      if (!hadDigits || errno == ERANGE || !end || *end != '\0' ||
          count < 1 || count > INT_MAX)
        throw std::invalid_argument(
          "APRL requires OMP_NUM_THREADS to be one positive integer; "
          "nested thread-count lists are not supported");
      return static_cast<int>(count);
    }
  }

  bool computeThreadingAvailable()
  {
#ifdef APRL_HAS_OPENMP
    return true;
#else
    return false;
#endif
  }

  int computeThreadCount() { return configuredThreads; }

  void configureComputeThreads(int providedMpiThreadLevel,
                               int requestedThreads)
  {
    if (requestedThreads < 0)
      throw std::invalid_argument("Compute thread count cannot be negative");
    const int requested = requestedThreads == 0 ?
      environmentThreadCount() : requestedThreads;
    if (requested > 1 && providedMpiThreadLevel < MPI_THREAD_FUNNELED)
      throw std::runtime_error(
        "Multiple compute threads require MPI_THREAD_FUNNELED support; "
        "this MPI runtime supplied a lower thread level");
#ifdef APRL_HAS_OPENMP
    if (omp_in_parallel())
      throw std::runtime_error(
        "Compute threads must be configured outside an OpenMP region");
    if (requested > omp_get_thread_limit())
      throw std::runtime_error(
        "Requested compute threads exceed the OpenMP thread limit");
    omp_set_dynamic(0);
    omp_set_num_threads(requested);
    int actual = 1;
    if (requested > 1)
      {
#pragma omp parallel num_threads(requested)
        {
#pragma omp single
          actual = omp_get_num_threads();
        }
      }
    if (actual != requested)
      throw std::runtime_error(
        "OpenMP could not supply the requested compute thread count (" +
        std::to_string(requested) + " requested, " +
        std::to_string(actual) + " available)");
#else
    if (requested > 1)
      throw std::runtime_error(
        "Multiple compute threads were requested, but APRL was built "
        "without OpenMP; rebuild with APRL_ENABLE_OPENMP=ON");
#endif
    configuredThreads = requested;
  }
}
