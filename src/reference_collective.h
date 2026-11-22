#ifndef APRL_REFERENCE_COLLECTIVE_H
#define APRL_REFERENCE_COLLECTIVE_H

#include <mpi.h>
#include <algorithm>
#include <cstdio>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>

namespace aprl { namespace verification { namespace reference_detail {

// MPI stays on the calling thread (MPI_THREAD_FUNNELED). Propagate local I/O,
// validation and allocation failures before another rank enters the next phase.
template<class Action> void collectively(MPI_Comm comm, Action action) {
  int rank=0;
  MPI_Comm_rank(comm,&rank);
  char message[2048]={};
  int failed=std::numeric_limits<int>::max();
  try { action(); }
  catch(const std::exception& e) {
    failed=rank; std::snprintf(message,sizeof(message),"%s",e.what());
  }
  catch(...) {
    failed=rank; std::snprintf(message,sizeof(message),"Unknown reference preparation failure");
  }
  int first=failed;
  MPI_Allreduce(&failed,&first,1,MPI_INT,MPI_MIN,comm);
  if(first!=std::numeric_limits<int>::max()) {
    MPI_Bcast(message,sizeof(message),MPI_CHAR,first,comm);
    throw std::runtime_error("Collective reference rank "+std::to_string(first)+": "+message);
  }
}

// Bounded counts also support tables larger than MPI's signed-int count limit.
inline void broadcast(double* values,std::size_t count,MPI_Comm comm) {
  const std::size_t batch=4*1024*1024;
  for(std::size_t p=0;p<count;p+=batch)
    MPI_Bcast(values+p,static_cast<int>(std::min(batch,count-p)),MPI_DOUBLE,0,comm);
}

}}} // namespace aprl::verification::reference_detail
#endif
