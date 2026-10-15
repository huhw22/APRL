#ifndef DIRECT_EB_SHUTDOWN_WATCHDOG_H
#define DIRECT_EB_SHUTDOWN_WATCHDOG_H

namespace aprl
{
  enum
  {
    shutdownWatchdogExitCode = 124
  };

  /* Process-local last-resort protection for teardown calls that can block
   * inside MPI or the C++ runtime.  It deliberately uses no MPI calls and is
   * never armed during the physical time-stepping loop or the normal explicit
   * HDF5 output-close phase.  Exceptional stack unwinding is guarded as a
   * whole because its destructors must not occupy a batch allocation forever. */
  void initializeShutdownWatchdog(int rank, bool threadSupportAvailable);
  void configureShutdownWatchdog(bool enabled,
                                 unsigned int timeoutSeconds);
  void armShutdownWatchdog(const char* phase);
  void disarmShutdownWatchdog();
  bool shutdownWatchdogOperational();
  void shutdownShutdownWatchdog();
}

#endif
