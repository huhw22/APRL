// Test-only fault injection; this library is never linked into APRL.
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <sys/syscall.h>
#include <unistd.h>

#include <hdf5.h>
#include <mpi.h>

namespace
{
  bool ledgerClosed = false;
  std::atomic<bool> rootStalled(false);
  std::atomic<std::size_t> exitMessageLength(0);
  char exitMessage[128] = {};

  bool enabled()
  {
    const char* value = std::getenv("APRL_TEST_STALL_DURATION_BARRIER");
    return value && std::strcmp(value, "1") == 0;
  }

  void waitForever()
  {
    for (;;) ::pause();
  }
}

extern "C" herr_t H5Fclose(hid_t file)
{
  typedef herr_t (*CloseFunction)(hid_t);
  static CloseFunction closeFile =
    reinterpret_cast<CloseFunction>(dlsym(RTLD_NEXT, "H5Fclose"));
  if (!closeFile) std::abort();
  char filename[4096] = {};
  const ssize_t length = H5Fget_name(file, filename, sizeof(filename));
  const char* basename = std::strrchr(filename, '/');
  basename = basename ? basename + 1 : filename;
  const bool isLedger = enabled() && length >= 0 &&
    static_cast<std::size_t>(length) < sizeof(filename) &&
    std::strcmp(basename, "energy-ledger.h5") == 0;
  const herr_t result = closeFile(file);
  // In this throughput fixture only root owns this writer. Its first barrier
  // after the final ledger close is the duration-exhaustion readiness barrier.
  // The driver also checks the production phase traces so this cannot silently
  // turn into a test of an unrelated earlier collective.
  if (isLedger && result >= 0) ledgerClosed = true;
  return result;
}

extern "C" int MPI_Barrier(MPI_Comm communicator)
{
  if (enabled())
    {
      int rank = -1;
      PMPI_Comm_rank(communicator, &rank);
      if (exitMessageLength.load(std::memory_order_relaxed) == 0)
        {
          const int length = std::snprintf(exitMessage, sizeof(exitMessage),
            "DURATION_BARRIER_WATCHDOG_EXIT rank=%d status=124\n", rank);
          if (length <= 0 || static_cast<std::size_t>(length) >= sizeof(exitMessage))
            std::abort();
          exitMessageLength.store(static_cast<std::size_t>(length),
            std::memory_order_release);
        }
      if (rank == 0 && ledgerClosed)
        {
          rootStalled.store(true, std::memory_order_release);
          std::fprintf(stdout,
            "DURATION_BARRIER_FAULT_READY rank=0 ledger_closed=1 action=withhold-barrier\n");
          std::fflush(stdout);
          waitForever();
        }
    }
  return PMPI_Barrier(communicator);
}

extern "C" void _exit(int status)
{
  const std::size_t length = exitMessageLength.load(std::memory_order_acquire);
  if (status == 124 && length > 0)
    {
      if (rootStalled.load(std::memory_order_acquire))
        {
          // Freeze the injected rank's watchdog too: a healthy peer must
          // independently exit its real MPI wait. Otherwise root's watchdog
          // could kill the launcher first and conceal an unprotected peer.
          static const char held[] =
            "DURATION_BARRIER_FAULT_ROOT_EXIT_HELD status=124\n";
          (void)::write(STDOUT_FILENO, held, sizeof(held) - 1);
          waitForever();
        }
      // No MPI, formatting, allocation, or locks on the waiting peer's exit
      // path. Preserve the production watchdog's exact process exit status.
      (void)::write(STDOUT_FILENO, exitMessage, length);
    }
  ::syscall(SYS_exit_group, status);
  __builtin_unreachable();
}
