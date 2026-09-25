#ifndef DIRECT_EB_RUNTIME_CONTROL_H
#define DIRECT_EB_RUNTIME_CONTROL_H

namespace aprl
{
  namespace RuntimeControl
  {
    /* Signal handlers only set a flag.  All HDF5 and MPI cleanup remains in
     * the ordinary control flow at the next field-step boundary. */
    void installSignalHandlers();
    bool stopRequested();
    int requestedSignal();
  }
}

#endif
