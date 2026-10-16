#include "shutdown_watchdog.h"

#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <exception>
#include <mutex>
#include <thread>

#include <unistd.h>

namespace aprl
{
  namespace
  {
    typedef std::chrono::steady_clock WatchdogClock;

    struct ShutdownWatchdogState
    {
      std::mutex mutex;
      std::condition_variable changed;
      std::thread worker;
      bool initialized;
      bool available;
      bool enabled;
      bool stop;
      bool armed;
      int rank;
      unsigned int timeoutSeconds;
      unsigned int finalizeTimeoutSeconds;
      unsigned long long generation;
      WatchdogClock::time_point deadline;

      ShutdownWatchdogState()
        : mutex(), changed(), worker(), initialized(false), available(false),
          enabled(true), stop(false), armed(false), rank(0),
          timeoutSeconds(600), finalizeTimeoutSeconds(60), generation(0),
          deadline()
      {}
    };

    ShutdownWatchdogState watchdogState;

    void armWithTimeout(bool finalize)
    {
      ShutdownWatchdogState& state = watchdogState;
      std::lock_guard<std::mutex> lock(state.mutex);
      if (!state.initialized || !state.available || !state.enabled) return;
      const unsigned int timeoutSeconds = finalize ?
        state.finalizeTimeoutSeconds : state.timeoutSeconds;
      state.armed = true;
      state.deadline = WatchdogClock::now() +
        std::chrono::seconds(timeoutSeconds);
      ++state.generation;
      state.changed.notify_all();
    }

    void watchdogLoop()
    {
      ShutdownWatchdogState& state = watchdogState;
      std::unique_lock<std::mutex> lock(state.mutex);
      while (!state.stop)
        {
          if (!state.armed)
            {
              state.changed.wait(lock, [&state]() {
                return state.stop || state.armed;
              });
              continue;
            }

          const unsigned long long generation = state.generation;
          const WatchdogClock::time_point deadline = state.deadline;
          const bool changed = state.changed.wait_until(
            lock, deadline, [&state, generation]() {
              return state.stop || !state.armed ||
                     state.generation != generation;
            });
          if (changed) continue;

          /* No MPI, HDF5, stream, allocation, or logging is attempted on the
           * timeout path.  _exit(124) makes cleanup failure unambiguously
           * non-successful even when scientific files were already closed. */
          lock.unlock();
          ::_exit(shutdownWatchdogExitCode);
        }
    }
  }

  void initializeShutdownWatchdog(int rank, bool threadSupportAvailable)
  {
    ShutdownWatchdogState& state = watchdogState;
    state.rank = rank;
    state.initialized = true;
    if (!threadSupportAvailable)
      {
        state.enabled = false;
        std::fprintf(stderr,
          "WARNING: shutdown watchdog unavailable on rank %d: MPI did not "
          "provide MPI_THREAD_FUNNELED support\n", rank);
        std::fflush(stderr);
        return;
      }

    try
      {
        state.worker = std::thread(watchdogLoop);
        state.available = true;
      }
    catch (const std::exception& error)
      {
        state.enabled = false;
        std::fprintf(stderr,
          "WARNING: shutdown watchdog unavailable on rank %d: could not "
          "start its local worker thread: %s\n", rank, error.what());
        std::fflush(stderr);
      }
  }

  void configureShutdownWatchdog(bool enabled,
                                 unsigned int timeoutSeconds,
                                 unsigned int finalizeTimeoutSeconds)
  {
    ShutdownWatchdogState& state = watchdogState;
    std::lock_guard<std::mutex> lock(state.mutex);
    state.armed = false;
    state.enabled = enabled && state.available && timeoutSeconds > 0 &&
      finalizeTimeoutSeconds > 0;
    state.timeoutSeconds = timeoutSeconds;
    state.finalizeTimeoutSeconds = finalizeTimeoutSeconds;
    ++state.generation;
    state.changed.notify_all();
  }

  void armShutdownWatchdog(const char* phase)
  {
    (void)phase;
    armWithTimeout(false);
  }

  void armShutdownWatchdogForFinalize(const char* phase)
  {
    (void)phase;
    armWithTimeout(true);
  }

  void disarmShutdownWatchdog()
  {
    ShutdownWatchdogState& state = watchdogState;
    std::lock_guard<std::mutex> lock(state.mutex);
    if (!state.initialized) return;
    state.armed = false;
    ++state.generation;
    state.changed.notify_all();
  }

  bool shutdownWatchdogOperational()
  {
    ShutdownWatchdogState& state = watchdogState;
    std::lock_guard<std::mutex> lock(state.mutex);
    return state.initialized && state.available && state.enabled;
  }

  void shutdownShutdownWatchdog()
  {
    ShutdownWatchdogState& state = watchdogState;
    {
      std::lock_guard<std::mutex> lock(state.mutex);
      if (!state.initialized) return;
      state.armed = false;
      state.stop = true;
      ++state.generation;
      state.changed.notify_all();
    }
    if (state.worker.joinable()) state.worker.join();
    state.available = false;
    state.enabled = false;
  }
}
