#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

#include <hdf5.h>
#include <mpi.h>
#include <yaml-cpp/yaml.h>

#ifdef APRL_PROBE_FFTW
#include <fftw3.h>
#endif

#ifdef APRL_PROBE_SDDS
#include <SDDS.h>
#endif

namespace
{
  void require(bool condition, const std::string& message)
  {
    if (!condition) throw std::runtime_error(message);
  }

  void requireHdf5(herr_t status, const std::string& message)
  {
    if (status < 0) throw std::runtime_error(message);
  }

  void requireHandle(hid_t handle, const std::string& message)
  {
    if (handle < 0) throw std::runtime_error(message);
  }

  void exerciseHdf5(const std::string& filename, int rank, int size)
  {
    unsigned int runtimeMajor = 0;
    unsigned int runtimeMinor = 0;
    unsigned int runtimeRelease = 0;
    requireHdf5(H5get_libversion(
      &runtimeMajor, &runtimeMinor, &runtimeRelease),
      "H5get_libversion failed");
    require(runtimeMajor == H5_VERS_MAJOR && runtimeMinor == H5_VERS_MINOR,
      "HDF5 headers and runtime library have different major/minor versions");

#ifdef H5_HAVE_PARALLEL
    const bool headerParallel = true;
#else
    const bool headerParallel = false;
#endif
    require(headerParallel == (APRL_PROBE_CMAKE_HDF5_PARALLEL != 0),
      "CMake and hdf5.h disagree about parallel-HDF5 support");
    if (!headerParallel && size > 1)
      throw std::runtime_error(
        "Serial HDF5 cannot perform the two-rank collective I/O probe");

    hid_t access = H5Pcreate(H5P_FILE_ACCESS);
    requireHandle(access, "Cannot create HDF5 file-access properties");
#ifdef H5_HAVE_PARALLEL
    requireHdf5(H5Pset_fapl_mpio(access, MPI_COMM_WORLD, MPI_INFO_NULL),
      "Cannot attach the active MPI communicator to HDF5");
#endif
    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,
      H5P_DEFAULT, access);
    H5Pclose(access);
    requireHandle(file, "Cannot collectively create the HDF5 probe file");

    const hsize_t globalCount[1] = {static_cast<hsize_t>(size)};
    hid_t fileSpace = H5Screate_simple(1, globalCount, nullptr);
    requireHandle(fileSpace, "Cannot create HDF5 probe dataspace");
    hid_t dataset = H5Dcreate2(file, "rank_values", H5T_IEEE_F64LE,
      fileSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(fileSpace);
    requireHandle(dataset, "Cannot create HDF5 probe dataset");

    fileSpace = H5Dget_space(dataset);
    requireHandle(fileSpace, "Cannot inspect HDF5 probe dataset");
    const hsize_t start[1] = {static_cast<hsize_t>(rank)};
    const hsize_t count[1] = {1};
    requireHdf5(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
      start, nullptr, count, nullptr), "Cannot select HDF5 rank hyperslab");
    hid_t memorySpace = H5Screate_simple(1, count, nullptr);
    requireHandle(memorySpace, "Cannot create HDF5 memory dataspace");
    hid_t transfer = H5Pcreate(H5P_DATASET_XFER);
    requireHandle(transfer, "Cannot create HDF5 transfer properties");
#ifdef H5_HAVE_PARALLEL
    requireHdf5(H5Pset_dxpl_mpio(transfer, H5FD_MPIO_COLLECTIVE),
      "Cannot select collective parallel-HDF5 transfer");
#endif
    const double value = static_cast<double>(rank + 1);
    requireHdf5(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memorySpace,
      fileSpace, transfer, &value), "HDF5 probe write failed");
    H5Pclose(transfer);
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
    H5Dclose(dataset);
    H5Fclose(file);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) std::remove(filename.c_str());

    if (rank == 0)
      std::cout << "hdf5_compile=" << H5_VERS_MAJOR << '.'
                << H5_VERS_MINOR << '.' << H5_VERS_RELEASE
                << " hdf5_runtime=" << runtimeMajor << '.'
                << runtimeMinor << '.' << runtimeRelease
                << " parallel=" << (headerParallel ? "yes" : "no")
                << '\n';
  }

  void exerciseFftw(int rank)
  {
#ifdef APRL_PROBE_FFTW
    require(fftw_init_threads() != 0, "fftw_init_threads failed");
    fftw_plan_with_nthreads(1);
    fftw_complex* values = static_cast<fftw_complex*>(
      fftw_malloc(8 * sizeof(fftw_complex)));
    require(values != nullptr, "fftw_malloc failed");
    for (int i = 0; i < 8; ++i)
      {
        values[i][0] = i == 0 ? 1.0 : 0.0;
        values[i][1] = 0.0;
      }
    fftw_plan plan = fftw_plan_dft_1d(
      8, values, values, FFTW_FORWARD, FFTW_ESTIMATE);
    require(plan != nullptr, "FFTW plan creation failed");
    fftw_execute(plan);
    require(values[0][0] == 1.0 && values[0][1] == 0.0,
      "FFTW execution returned an unexpected DC component");
    fftw_destroy_plan(plan);
    fftw_free(values);
    fftw_cleanup_threads();
    if (rank == 0) std::cout << "fftw_threads=ok\n";
#else
    (void)rank;
#endif
  }
}

int main(int argc, char** argv)
{
  int provided = MPI_THREAD_SINGLE;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int localStatus = EXIT_SUCCESS;
  try
    {
      require(provided >= MPI_THREAD_FUNNELED,
        "MPI does not provide MPI_THREAD_FUNNELED");
      const YAML::Node document = YAML::Load("probe: [1, 2, 3]");
      require(document["probe"].size() == 3,
        "yaml-cpp parse/link probe failed");
      const std::string filename = argc > 1 ? argv[1] :
        "fel-dependency-probe.h5";
      exerciseHdf5(filename, rank, size);
      exerciseFftw(rank);
#ifdef APRL_PROBE_SDDS
      int32_t (*initializeInput)(SDDS_DATASET*, char*) =
        &SDDS_InitializeInput;
      require(initializeInput != nullptr, "SDDS link probe failed");
      if (rank == 0) std::cout << "official_sdds_link=ok\n";
#endif
      if (rank == 0)
        {
          int mpiVersionLength = 0;
          char mpiVersion[MPI_MAX_LIBRARY_VERSION_STRING] = {0};
          MPI_Get_library_version(mpiVersion, &mpiVersionLength);
          while (mpiVersionLength > 0 &&
                 (mpiVersion[mpiVersionLength - 1] == '\0' ||
                  mpiVersion[mpiVersionLength - 1] == '\n' ||
                  mpiVersion[mpiVersionLength - 1] == '\r'))
            --mpiVersionLength;
          std::cout << "mpi_ranks=" << size
                    << " mpi_thread_funneled=yes yaml_cpp=ok\n"
                    << "mpi_library="
                    << std::string(mpiVersion,
                         static_cast<std::size_t>(mpiVersionLength)) << '\n';
        }
    }
  catch (const std::exception& error)
    {
      std::cerr << "dependency probe failed on rank " << rank << ": "
                << error.what() << '\n';
      localStatus = EXIT_FAILURE;
    }
  int globalStatus = EXIT_SUCCESS;
  MPI_Allreduce(&localStatus, &globalStatus, 1, MPI_INT, MPI_MAX,
    MPI_COMM_WORLD);
  MPI_Finalize();
  return globalStatus;
}
