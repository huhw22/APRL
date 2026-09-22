#include "runtime_control.h"

#include <csignal>

namespace fel
{
  namespace RuntimeControl
  {
    namespace
    {
      volatile std::sig_atomic_t stopSignal = 0;

      extern "C" void requestStop(int signalNumber)
      {
        stopSignal = signalNumber;
      }
    }

    void installSignalHandlers()
    {
      std::signal(SIGINT, requestStop);
      std::signal(SIGTERM, requestStop);
    }

    bool stopRequested()
    {
      return stopSignal != 0;
    }

    int requestedSignal()
    {
      return static_cast<int>(stopSignal);
    }
  }
}
