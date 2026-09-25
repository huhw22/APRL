cmake_minimum_required(VERSION 3.16)

project(fel_dependency_probe LANGUAGES C CXX)

option(FEL_PROBE_FULL "Check FFTW post-processing dependencies" ON)
option(FEL_PROBE_REQUIRE_PARALLEL_HDF5
  "Require an MPI-enabled HDF5 C library" ON)
option(FEL_PROBE_SDDS "Check a supplied official SDDS build" OFF)
set(SDDS_ROOT "" CACHE PATH "Root of a built official SDDS source tree")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED COMPONENTS CXX)
set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 1.10 REQUIRED COMPONENTS C)
find_package(yaml-cpp 0.6 REQUIRED)
find_package(Threads REQUIRED)

if(FEL_PROBE_REQUIRE_PARALLEL_HDF5 AND NOT HDF5_IS_PARALLEL)
  message(FATAL_ERROR
    "The selected HDF5 is serial. A production MPI build needs parallel "
    "HDF5 compiled against the same MPI family.")
endif()

add_executable(fel_dependency_probe probe.cpp)
target_include_directories(fel_dependency_probe PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(fel_dependency_probe PRIVATE
  MPI::MPI_CXX
  ${HDF5_LIBRARIES}
  yaml-cpp
  Threads::Threads
)
target_compile_definitions(fel_dependency_probe PRIVATE
  OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX)
if(HDF5_IS_PARALLEL)
  target_compile_definitions(fel_dependency_probe PRIVATE
    FEL_PROBE_CMAKE_HDF5_PARALLEL=1)
else()
  target_compile_definitions(fel_dependency_probe PRIVATE
    FEL_PROBE_CMAKE_HDF5_PARALLEL=0)
endif()

if(FEL_PROBE_FULL)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
  find_library(FFTW3_THREADS_LIBRARY NAMES fftw3_threads REQUIRED)
  target_link_libraries(fel_dependency_probe PRIVATE
    PkgConfig::FFTW3 ${FFTW3_THREADS_LIBRARY})
  target_compile_definitions(fel_dependency_probe PRIVATE FEL_PROBE_FFTW=1)
endif()

if(FEL_PROBE_SDDS)
  if(NOT SDDS_ROOT)
    message(FATAL_ERROR "FEL_PROBE_SDDS requires -DSDDS_ROOT=/built/SDDS")
  endif()
  find_path(SDDS_INCLUDE_DIR SDDS.h
    HINTS ${SDDS_ROOT}
    PATH_SUFFIXES include)
  set(_SDDS_PLATFORM "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
  set(_SDDS_SUFFIXES
    lib/${_SDDS_PLATFORM} lib/Linux-x86_64
    SDDSlib/O.${_SDDS_PLATFORM} SDDSlib/O.Linux-x86_64
    rpns/code/O.${_SDDS_PLATFORM} rpns/code/O.Linux-x86_64
    mdbmth/O.${_SDDS_PLATFORM} mdbmth/O.Linux-x86_64
    mdblib/O.${_SDDS_PLATFORM} mdblib/O.Linux-x86_64 lib)
  find_library(SDDS_LIBRARY NAMES SDDS1 SDDS
    HINTS ${SDDS_ROOT} PATH_SUFFIXES ${_SDDS_SUFFIXES})
  find_library(SDDS_RPN_LIBRARY NAMES rpnlib
    HINTS ${SDDS_ROOT} PATH_SUFFIXES ${_SDDS_SUFFIXES})
  find_library(SDDS_MDBMTH_LIBRARY NAMES mdbmth
    HINTS ${SDDS_ROOT} PATH_SUFFIXES ${_SDDS_SUFFIXES})
  find_library(SDDS_MDB_LIBRARY NAMES mdblib
    HINTS ${SDDS_ROOT} PATH_SUFFIXES ${_SDDS_SUFFIXES})
  if(NOT SDDS_INCLUDE_DIR OR NOT SDDS_LIBRARY OR NOT SDDS_RPN_LIBRARY OR
     NOT SDDS_MDBMTH_LIBRARY OR NOT SDDS_MDB_LIBRARY)
    message(FATAL_ERROR
      "SDDS_ROOT does not contain the required SDDS1, rpnlib, mdbmth and "
      "mdblib build products for this platform")
  endif()
  find_package(ZLIB REQUIRED)
  find_package(LibLZMA REQUIRED)
  find_package(GSL REQUIRED)
  target_include_directories(fel_dependency_probe PRIVATE ${SDDS_INCLUDE_DIR})
  target_compile_definitions(fel_dependency_probe PRIVATE
    FEL_PROBE_SDDS=1 zLib)
  target_link_libraries(fel_dependency_probe PRIVATE
    ${SDDS_LIBRARY}
    ${SDDS_RPN_LIBRARY}
    ${SDDS_MDBMTH_LIBRARY}
    ${SDDS_MDB_LIBRARY}
    LibLZMA::LibLZMA ZLIB::ZLIB GSL::gsl GSL::gslcblas
    ${CMAKE_DL_LIBS})
  if(UNIX)
    target_link_libraries(fel_dependency_probe PRIVATE m)
  endif()
endif()

message(STATUS "Probe MPI compiler: ${MPI_CXX_COMPILER}")
message(STATUS "Probe HDF5: ${HDF5_VERSION}; parallel=${HDF5_IS_PARALLEL}")
message(STATUS "Probe full post-processing dependencies: ${FEL_PROBE_FULL}")
message(STATUS "Probe official SDDS build: ${FEL_PROBE_SDDS}")

