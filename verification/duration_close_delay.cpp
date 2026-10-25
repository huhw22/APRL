// Test-only HDF5 interposer. It is never linked into the APRL executable.
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <time.h>

#include <hdf5.h>

extern "C" herr_t H5Fclose(hid_t file)
{
  typedef herr_t (*CloseFunction)(hid_t);
  static CloseFunction closeFile =
    reinterpret_cast<CloseFunction>(dlsym(RTLD_NEXT, "H5Fclose"));
  if (!closeFile)
    {
      std::fprintf(stderr, "DURATION_CLOSE_DELAY_ERROR cannot resolve H5Fclose\n");
      std::abort();
    }
  const char* enabled = std::getenv("APRL_TEST_DELAY_H5_CLOSE");
  char filename[4096] = {};
  const ssize_t length = H5Fget_name(file, filename, sizeof(filename));
  const char* basename = std::strrchr(filename, '/');
  basename = basename ? basename + 1 : filename;
  static bool delayed = false;
  if (!delayed && enabled && std::strcmp(enabled, "1") == 0 &&
      length >= 0 && static_cast<std::size_t>(length) < sizeof(filename) &&
      std::strcmp(basename, "duration-field.h5") == 0)
    {
      delayed = true;
      std::fprintf(stdout,
        "DURATION_CLOSE_DELAY_BEGIN file=duration-field.h5 seconds=3\n");
      std::fflush(stdout);
      struct timespec remaining = {3, 0};
      while (nanosleep(&remaining, &remaining) != 0 && errno == EINTR) {}
      std::fprintf(stdout, "DURATION_CLOSE_DELAY_COMPLETE file=duration-field.h5\n");
      std::fflush(stdout);
    }
  return closeFile(file);
}
