#include "teardown_trace.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <unistd.h>

namespace aprl
{
  namespace
  {
    typedef std::chrono::steady_clock TeardownClock;

    struct TeardownTraceState
    {
      bool initialized;
      int rank;
      bool allRanks;
      char hostname[256];
      long pid;
      TeardownClock::time_point start;

      TeardownTraceState()
        : initialized(false), rank(0), allRanks(false), hostname(), pid(0),
          start()
      {}
    };

    TeardownTraceState& traceState()
    {
      static TeardownTraceState state;
      return state;
    }
  }

  void initializeTeardownTrace(int rank)
  {
    TeardownTraceState& state = traceState();
    state.rank = rank;
    const char* value = std::getenv("APRL_TEARDOWN_TRACE");
    state.allRanks = value != 0 && std::strcmp(value, "all") == 0;
    std::memset(state.hostname, 0, sizeof(state.hostname));
    if (::gethostname(state.hostname, sizeof(state.hostname) - 1) != 0)
      std::strncpy(state.hostname, "unknown", sizeof(state.hostname) - 1);
    state.hostname[sizeof(state.hostname) - 1] = '\0';
    state.pid = static_cast<long>(::getpid());
    state.start = TeardownClock::now();
    state.initialized = true;
  }

  void logTeardownTrace(const char* phase)
  {
    const TeardownTraceState& state = traceState();
    if (!state.initialized || (state.rank != 0 && !state.allRanks)) return;
    const double elapsedSeconds = std::chrono::duration<double>(
      TeardownClock::now() - state.start).count();
    std::fprintf(stdout,
      "[teardown] rank=%d host=%s pid=%ld elapsed_s=%.9f phase=%s\n",
      state.rank, state.hostname, state.pid, elapsedSeconds,
      phase ? phase : "unknown");
    std::fflush(stdout);
  }
}
