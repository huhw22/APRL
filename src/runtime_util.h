#ifndef DIRECT_EB_RUNTIME_UTIL_H
#define DIRECT_EB_RUNTIME_UTIL_H

#include <string>

#include <mpi.h>

namespace fel
{
  void logRoot(MPI_Comm communicator, const std::string& message);
  void createDirectories(const std::string& directory,
                         MPI_Comm communicator);
  std::string joinPath(const std::string& directory,
                       const std::string& filename);
}

#endif
