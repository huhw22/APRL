#include "runtime_util.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <system_error>
#include <thread>
#include <vector>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/stat.h>

namespace aprl
{
  void logRoot(MPI_Comm communicator, const std::string& message)
  {
    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    if (rank == 0) std::cout << message << std::endl;
  }

  namespace
  {
    bool directoryExists(const std::string& path)
    {
      struct stat status;
      return stat(path.c_str(), &status) == 0 && S_ISDIR(status.st_mode);
    }

    void createDirectoryTree(const std::string& directory)
    {
      if (directory.empty() || directory == "." || directoryExists(directory))
        return;
      std::string current;
      std::size_t offset = 0;
      if (directory[0] == '/')
        {
          current = "/";
          offset = 1;
        }
      while (offset <= directory.size())
        {
          const std::size_t separator = directory.find('/', offset);
          const std::string part = directory.substr(
            offset, separator == std::string::npos ? std::string::npos :
            separator - offset);
          if (!part.empty())
            {
              if (!current.empty() && current[current.size() - 1] != '/')
                current += '/';
              current += part;
              if (!directoryExists(current) &&
                  mkdir(current.c_str(), 0775) != 0 && errno != EEXIST)
                throw std::system_error(errno, std::generic_category(),
                  "Cannot create directory '" + current + "'");
            }
          if (separator == std::string::npos) break;
          offset = separator + 1;
        }
    }
  }

  void createDirectories(const std::string& directory,
                         MPI_Comm communicator)
  {
    int rank = 0, ranks = 0;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &ranks);
    int failedRank = ranks;
    std::string message;
    // A barrier after rank-zero mkdir does not invalidate another host's
    // cached negative directory lookup from output preflight. Each client
    // must establish visibility before its independent HDF5 file is opened.
    // Concurrent mkdir is idempotent; open(O_DIRECTORY) also rejects a regular
    // file which mkdir reports as EEXIST. No file is removed or overwritten.
    try
      {
        const auto deadline = std::chrono::steady_clock::now() +
          std::chrono::seconds(5);
        for (;;)
          {
            try
              {
                createDirectoryTree(directory);
                const char* path = directory.empty() ? "." : directory.c_str();
                const int fd = ::open(path, O_RDONLY | O_DIRECTORY | O_CLOEXEC);
                if (fd < 0)
                  throw std::system_error(errno, std::generic_category(),
                    "Cannot open output directory '" + directory + "'");
                ::close(fd);
                break;
              }
            catch (const std::system_error& error)
              {
                if ((error.code().value() != ENOENT &&
                     error.code().value() != ESTALE) ||
                    std::chrono::steady_clock::now() >= deadline)
                  throw;
                std::this_thread::sleep_for(std::chrono::milliseconds(20));
              }
          }
      }
    catch (const std::exception& error)
      {
        failedRank = rank;
        message = error.what();
        char cwd[4096];
        if (::getcwd(cwd, sizeof(cwd)))
          message += std::string(" (working directory '") + cwd + "')";
      }
    int firstFailure = ranks;
    MPI_Allreduce(&failedRank, &firstFailure, 1, MPI_INT, MPI_MIN,
                  communicator);
    if (firstFailure < ranks)
      {
        int length = rank == firstFailure ? static_cast<int>(message.size()) : 0;
        MPI_Bcast(&length, 1, MPI_INT, firstFailure, communicator);
        std::vector<char> bytes(static_cast<std::size_t>(length));
        if (rank == firstFailure)
          std::copy(message.begin(), message.end(), bytes.begin());
        MPI_Bcast(bytes.data(), length, MPI_CHAR, firstFailure, communicator);
        throw std::runtime_error(std::string(bytes.begin(), bytes.end()));
      }
  }

  std::string joinPath(const std::string& directory,
                       const std::string& filename)
  {
    if (directory.empty() || directory == ".") return filename;
    if (directory[directory.size() - 1] == '/') return directory + filename;
    return directory + "/" + filename;
  }

  std::uint64_t currentResidentBytes()
  {
    std::ifstream status("/proc/self/statm");
    unsigned long long totalPages = 0;
    unsigned long long residentPages = 0;
    if (!(status >> totalPages >> residentPages)) return 0;
    const long pageBytes = sysconf(_SC_PAGESIZE);
    if (pageBytes <= 0) return 0;
    if (residentPages >
        std::numeric_limits<std::uint64_t>::max() /
          static_cast<std::uint64_t>(pageBytes))
      return 0;
    return static_cast<std::uint64_t>(residentPages) *
      static_cast<std::uint64_t>(pageBytes);
  }

  std::uint64_t peakResidentBytes()
  {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0 ||
        usage.ru_maxrss < 0)
      return 0;
#ifdef __APPLE__
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    const std::uint64_t kib =
      static_cast<std::uint64_t>(usage.ru_maxrss);
    if (kib > std::numeric_limits<std::uint64_t>::max() / 1024)
      return 0;
    return kib * 1024;
#endif
  }
}
