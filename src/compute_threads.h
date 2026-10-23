#ifndef APRL_COMPUTE_THREADS_H
#define APRL_COMPUTE_THREADS_H

namespace aprl
{
  /* Call on the MPI initialization thread, before numerical work.  Zero
   * selects OMP_NUM_THREADS when present and otherwise one compute thread.
   * Explicit positive counts are useful for controlled comparisons. */
  void configureComputeThreads(int providedMpiThreadLevel,
                               int requestedThreads = 0);
  int computeThreadCount();
  bool computeThreadingAvailable();
}

#endif
