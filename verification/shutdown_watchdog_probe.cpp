#include <cstring>
#include <iostream>

#include <mpi.h>
#include <unistd.h>

#include "shutdown_watchdog.h"

namespace
{
  void waitForever()
  {
    for (;;) ::pause();
  }
}

int main(int argc, char** argv)
{
  if (argc != 2)
    {
      std::cerr << "usage: shutdown_watchdog_probe "
                   "<timeout|finalize-timeout|disarm|disabled|"
                   "mpi-collective|mpi-finalize>\n";
      return 2;
    }

  const char* mode = argv[1];
  if (std::strcmp(mode, "timeout") == 0)
    {
      aprl::initializeShutdownWatchdog(7, true);
      aprl::configureShutdownWatchdog(true, 1, 3);
      aprl::armShutdownWatchdog("probe-timeout");
      waitForever();
    }
  if (std::strcmp(mode, "finalize-timeout") == 0)
    {
      aprl::initializeShutdownWatchdog(7, true);
      aprl::configureShutdownWatchdog(true, 3, 1);
      aprl::armShutdownWatchdogForFinalize("probe-finalize-timeout");
      waitForever();
    }
  if (std::strcmp(mode, "disarm") == 0)
    {
      aprl::initializeShutdownWatchdog(7, true);
      aprl::configureShutdownWatchdog(true, 1, 1);
      aprl::armShutdownWatchdog("probe-disarm");
      aprl::disarmShutdownWatchdog();
      ::sleep(2);
      aprl::shutdownShutdownWatchdog();
      return 0;
    }
  if (std::strcmp(mode, "disabled") == 0)
    {
      aprl::initializeShutdownWatchdog(7, true);
      aprl::configureShutdownWatchdog(false, 1, 1);
      aprl::armShutdownWatchdog("probe-disabled");
      ::sleep(2);
      aprl::shutdownShutdownWatchdog();
      return 0;
    }
  int providedThreadLevel = MPI_THREAD_SINGLE;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,
                  &providedThreadLevel);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  aprl::initializeShutdownWatchdog(
    rank, providedThreadLevel >= MPI_THREAD_FUNNELED);
  aprl::configureShutdownWatchdog(true, 1, 1);

  if (std::strcmp(mode, "mpi-collective") == 0)
    {
      std::cout << "WATCHDOG_PROBE_READY rank=" << rank << std::endl;
      MPI_Barrier(MPI_COMM_WORLD);
      aprl::armShutdownWatchdog("probe-mpi-collective");
      if (rank == 0)
        MPI_Barrier(MPI_COMM_WORLD);
      else
        waitForever();
      return 3;
    }
  if (std::strcmp(mode, "mpi-finalize") == 0)
    {
      aprl::armShutdownWatchdogForFinalize("probe-mpi-finalize");
      MPI_Finalize();
      aprl::disarmShutdownWatchdog();
      aprl::shutdownShutdownWatchdog();
      return 0;
    }

  MPI_Finalize();
  aprl::shutdownShutdownWatchdog();
  std::cerr << "unknown shutdown watchdog probe mode: " << mode << "\n";
  return 2;
}
