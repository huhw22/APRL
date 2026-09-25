#ifndef DIRECT_EB_RUNTIME_UTIL_H
#define DIRECT_EB_RUNTIME_UTIL_H

#include <cstdint>
#include <string>

#include <mpi.h>

namespace aprl
{
  void logRoot(MPI_Comm communicator, const std::string& message);
  void createDirectories(const std::string& directory,
                         MPI_Comm communicator);
  std::string joinPath(const std::string& directory,
                       const std::string& filename);
  std::uint64_t currentResidentBytes();
  std::uint64_t peakResidentBytes();
}

#endif
