#ifndef DIRECT_EB_TEARDOWN_TRACE_H
#define DIRECT_EB_TEARDOWN_TRACE_H

namespace aprl
{
  /* Cache process identity while MPI is available.  The logging function is
   * deliberately MPI-free, including when it is called after MPI_Finalize. */
  void initializeTeardownTrace(int rank);
  void logTeardownTrace(const char* phase);
}

#endif
