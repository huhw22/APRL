#include "runtime_util.h"

#include <cerrno>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <sys/stat.h>

namespace fel
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
                throw std::runtime_error("Cannot create directory '" +
                  current + "': " + std::strerror(errno));
            }
          if (separator == std::string::npos) break;
          offset = separator + 1;
        }
    }
  }

  void createDirectories(const std::string& directory,
                         MPI_Comm communicator)
  {
    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    int failed = 0;
    std::string message;
    if (rank == 0)
      {
        try { createDirectoryTree(directory); }
        catch (const std::exception& error) {
          failed = 1;
          message = error.what();
        }
      }
    MPI_Bcast(&failed, 1, MPI_INT, 0, communicator);
    if (failed)
      {
        if (rank == 0) throw std::runtime_error(message);
        throw std::runtime_error("Output directory creation failed on rank zero");
      }
    MPI_Barrier(communicator);
  }

  std::string joinPath(const std::string& directory,
                       const std::string& filename)
  {
    if (directory.empty() || directory == ".") return filename;
    if (directory[directory.size() - 1] == '/') return directory + filename;
    return directory + "/" + filename;
  }
}
