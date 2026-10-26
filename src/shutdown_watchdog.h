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
   * never armed during ordinary physical time stepping or the normal explicit
   * HDF5 output-close phase. Checkpoint I/O may arm its separate deadline.
   * Exceptional stack unwinding is guarded as a
   * whole because its destructors must not occupy a batch allocation forever. */
  void initializeShutdownWatchdog(int rank, bool threadSupportAvailable);
  void configureShutdownWatchdog(bool enabled,
                                 unsigned int timeoutSeconds,
                                 unsigned int finalizeTimeoutSeconds);
  /* General teardown keeps the first limit; MPI_Finalize uses the shorter,
   * independent limit because all normal scientific file closes precede it. */
  void armShutdownWatchdog(const char* phase);
  void armShutdownWatchdogForFinalize(const char* phase);
  /* Uses the existing process-local worker with an independent I/O limit,
   * even when general teardown protection was disabled. Returns false if
   * the worker is unavailable or the limit is zero; callers must reject an
   * unprotected checkpoint operation. Disarm after checkpoint I/O completes. */
  bool armCheckpointWatchdog(unsigned int timeoutSeconds);
  void disarmShutdownWatchdog();
  bool shutdownWatchdogOperational();
  void shutdownShutdownWatchdog();
}

#endif
