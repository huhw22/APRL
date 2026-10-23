#include <exception>
#include <iostream>
#include <string>

#include <mpi.h>

#include "config.h"
#include "compute_threads.h"
#include "runtime_control.h"
#include "runtime_util.h"
#include "shutdown_watchdog.h"
#include "simulation.h"
#include "teardown_trace.h"

int main(int argc, char** argv)
{
  int providedThreadLevel = MPI_THREAD_SINGLE;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,
                  &providedThreadLevel);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  aprl::initializeTeardownTrace(rank);
  aprl::initializeShutdownWatchdog(
    rank, providedThreadLevel >= MPI_THREAD_FUNNELED);

  if (argc != 2)
    {
      if (rank == 0)
        std::cerr << "Usage: aprl <configuration.yaml>" << std::endl;
      aprl::logTeardownTrace("mpi-finalize-enter");
      aprl::armShutdownWatchdogForFinalize("mpi-finalize");
      MPI_Finalize();
      aprl::disarmShutdownWatchdog();
      aprl::shutdownShutdownWatchdog();
      // Use only cached process data and the monotonic process clock here:
      // MPI calls are invalid after MPI_Finalize returns.
      aprl::logTeardownTrace("mpi-finalize-complete");
      return 1;
    }

  const double start = MPI_Wtime();
  try
    {
      aprl::configureComputeThreads(providedThreadLevel);
      aprl::SimulationConfig config =
        aprl::YamlConfigLoader::loadFile(argv[1]);
      aprl::configureShutdownWatchdog(
        config.runtime.shutdownWatchdogEnabled,
        config.runtime.shutdownWatchdogTimeoutSeconds,
        config.runtime.shutdownWatchdogFinalizeTimeoutSeconds);
      config.runtime.shutdownWatchdogEnabled =
        aprl::shutdownWatchdogOperational();
      if (config.runtime.interactive())
        aprl::RuntimeControl::installSignalHandlers();

      aprl::logRoot(MPI_COMM_WORLD,
        "APRL — Accelerator Particles and Radiation in Lorentz Frames");
      aprl::logRoot(MPI_COMM_WORLD,
        "Compute threads per MPI rank: " +
        std::to_string(aprl::computeThreadCount()) +
        "; MPI and HDF5 calls stay on the initialization thread.");
      aprl::logRoot(MPI_COMM_WORLD,
        std::string("Shutdown watchdog: enabled=") +
        (config.runtime.shutdownWatchdogEnabled ? "true" : "false") +
        ", timeout=" +
        std::to_string(config.runtime.shutdownWatchdogTimeoutSeconds) +
        " s, finalize_timeout_seconds=" +
        std::to_string(
          config.runtime.shutdownWatchdogFinalizeTimeoutSeconds) +
        " s; guards final MPI/runtime teardown only, not physical stepping "
        "or the normal explicit HDF5 output-close phase.");
      {
        aprl::Simulation simulation(config, MPI_COMM_WORLD);
        try
          {
            simulation.solve();

            if (rank == 0)
              std::cout << "Elapsed wall time [s]: "
                        << (MPI_Wtime() - start) << std::endl;
          }
        catch (...)
          {
            /* Stack unwinding destroys simulation before the outer catch can
             * call MPI_Abort, so guard that exceptional destructor too. */
            aprl::logTeardownTrace("simulation-destructor-enter");
            aprl::armShutdownWatchdog("simulation-destructor");
            throw;
          }
        aprl::logTeardownTrace("simulation-destructor-enter");
        aprl::armShutdownWatchdog("simulation-destructor");
      }
      aprl::disarmShutdownWatchdog();
      aprl::logTeardownTrace("simulation-destructor-complete");
    }
  catch (const std::exception& error)
    {
      aprl::armShutdownWatchdog("exception-report-and-mpi-abort");
      aprl::logTeardownTrace("exception-handler-enter");
      std::cerr << "Simulation failure on MPI rank " << rank
                << ": " << error.what() << std::endl;
      aprl::logTeardownTrace("mpi-abort-enter");
      aprl::armShutdownWatchdog("mpi-abort");
      MPI_Abort(MPI_COMM_WORLD, 1);
      aprl::disarmShutdownWatchdog();
      aprl::shutdownShutdownWatchdog();
      return 1;
    }

  /* Do not let a rank start the short MPI_Finalize timeout while a slower
   * peer is still closing output or destroying Simulation.  This readiness
   * collective runs only after every local scientific close has completed
   * and retains the longer general teardown timeout. */
  aprl::logTeardownTrace("finalize-readiness-sync-enter");
  aprl::armShutdownWatchdog("finalize-readiness-sync");
  MPI_Barrier(MPI_COMM_WORLD);
  aprl::disarmShutdownWatchdog();
  aprl::logTeardownTrace("finalize-readiness-sync-complete");

  aprl::logTeardownTrace("mpi-finalize-enter");
  aprl::armShutdownWatchdogForFinalize("mpi-finalize");
  MPI_Finalize();
  aprl::disarmShutdownWatchdog();
  aprl::shutdownShutdownWatchdog();
  // Do not query MPI state, rank, host identity, or MPI_Wtime after finalize.
  aprl::logTeardownTrace("mpi-finalize-complete");
  return 0;
}
