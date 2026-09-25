cmake_minimum_required(VERSION 3.16)

# Working project identifier only. The program's public name is intentionally
# left undecided; the executable name "simulator" is also temporary.
project(unnamed_fel_program LANGUAGES CXX)

include(CTest)

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(yaml-cpp REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB QUIET)
find_package(LibLZMA QUIET)
find_package(GSL QUIET)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(DIRECT_EB_SOURCES
  src/boostframe.cpp
  src/config.cpp
  src/eb_bunch.cpp
  src/eb_cpml.cpp
  src/eb_deposition.cpp
  src/eb_field.cpp
  src/eb_incident.cpp
  src/eb_initial_field.cpp
  src/eb_mpi.cpp
  src/eb_particles.cpp
  src/eb_sources.cpp
  src/energy_ledger.cpp
  src/lab_detectors.cpp
  src/particle_file.cpp
  src/particle_initializer.cpp
  src/particle_boundary.cpp
  src/runtime_control.cpp
  src/runtime_util.cpp
  src/simulation.cpp
  src/trajectory_writer.cpp
)

add_library(direct_eb_core STATIC ${DIRECT_EB_SOURCES})
target_include_directories(direct_eb_core
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(direct_eb_core
  PUBLIC
    MPI::MPI_CXX
    ${HDF5_LIBRARIES}
    yaml-cpp
)
target_compile_definitions(direct_eb_core
  PUBLIC OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX
)
target_compile_options(direct_eb_core PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)

add_executable(simulator src/main.cpp)
target_link_libraries(simulator PRIVATE direct_eb_core)

add_executable(particle_text_to_hdf5
  utilities/particle_text_to_hdf5.cpp
)
target_include_directories(particle_text_to_hdf5 PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(particle_text_to_hdf5 PRIVATE
  MPI::MPI_CXX
  ${HDF5_LIBRARIES}
)
target_compile_options(particle_text_to_hdf5 PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)

# Elegant's SDDS5 protocol includes binary layout, endian, page, compression,
# and 64-bit integer details.  Use the official SDDS library rather than a
# partial in-tree parser.  SDDS_ROOT may point at a built SDDS source tree;
# the sibling-tree hint matches the recommended Elegant/SDDS source layout.
set(SDDS_ROOT "" CACHE PATH "Root of a built official SDDS source tree")
find_path(SDDS_INCLUDE_DIR SDDS.h
  HINTS ${SDDS_ROOT} ENV SDDS_ROOT
  PATHS ${CMAKE_CURRENT_SOURCE_DIR}/../SDDS
  PATH_SUFFIXES include
)
set(_SDDS_LIBRARY_HINTS
  ${SDDS_ROOT}
  $ENV{SDDS_ROOT}
  ${CMAKE_CURRENT_SOURCE_DIR}/../SDDS
)
set(_SDDS_PLATFORM "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
find_library(SDDS_LIBRARY NAMES SDDS1 SDDS
  HINTS ${_SDDS_LIBRARY_HINTS}
  PATH_SUFFIXES lib/${_SDDS_PLATFORM} lib/Linux-x86_64
                SDDSlib/O.${_SDDS_PLATFORM} SDDSlib/O.Linux-x86_64 lib
)
find_library(SDDS_RPN_LIBRARY NAMES rpnlib
  HINTS ${_SDDS_LIBRARY_HINTS}
  PATH_SUFFIXES lib/${_SDDS_PLATFORM} lib/Linux-x86_64
                rpns/code/O.${_SDDS_PLATFORM} rpns/code/O.Linux-x86_64 lib
)
find_library(SDDS_MDBMTH_LIBRARY NAMES mdbmth
  HINTS ${_SDDS_LIBRARY_HINTS}
  PATH_SUFFIXES lib/${_SDDS_PLATFORM} lib/Linux-x86_64
                mdbmth/O.${_SDDS_PLATFORM} mdbmth/O.Linux-x86_64 lib
)
find_library(SDDS_MDB_LIBRARY NAMES mdblib
  HINTS ${_SDDS_LIBRARY_HINTS}
  PATH_SUFFIXES lib/${_SDDS_PLATFORM} lib/Linux-x86_64
                mdblib/O.${_SDDS_PLATFORM} mdblib/O.Linux-x86_64 lib
)
if(SDDS_INCLUDE_DIR AND SDDS_LIBRARY AND SDDS_RPN_LIBRARY AND
   SDDS_MDBMTH_LIBRARY AND SDDS_MDB_LIBRARY AND ZLIB_FOUND AND
   LIBLZMA_FOUND AND GSL_FOUND)
  add_executable(elegant_sdds_to_hdf5
    utilities/elegant_sdds_to_hdf5.cpp
  )
  target_include_directories(elegant_sdds_to_hdf5 PRIVATE
    ${HDF5_INCLUDE_DIRS}
    ${SDDS_INCLUDE_DIR}
  )
  target_compile_definitions(elegant_sdds_to_hdf5 PRIVATE zLib)
  target_link_libraries(elegant_sdds_to_hdf5 PRIVATE
    ${HDF5_LIBRARIES}
    ${SDDS_LIBRARY}
    ${SDDS_RPN_LIBRARY}
    ${SDDS_MDBMTH_LIBRARY}
    ${SDDS_MDB_LIBRARY}
    LibLZMA::LibLZMA
    ZLIB::ZLIB
    GSL::gsl
    GSL::gslcblas
    Threads::Threads
    ${CMAKE_DL_LIBS}
  )
  if(UNIX)
    target_link_libraries(elegant_sdds_to_hdf5 PRIVATE m)
  endif()
  target_compile_options(elegant_sdds_to_hdf5 PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
  )
  set(HAVE_ELEGANT_SDDS_CONVERTER TRUE)
else()
  message(STATUS
    "Official SDDS library not found; elegant_sdds_to_hdf5 will not be built. "
    "Set SDDS_ROOT to a built https://github.com/rtsoliday/SDDS tree.")
endif()

add_executable(undulator_resonance utilities/undulator_resonance.cpp)
target_link_libraries(undulator_resonance PRIVATE yaml-cpp)
target_compile_options(undulator_resonance PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)

add_executable(energy_ledger_report utilities/energy_ledger_report.cpp)
target_link_libraries(energy_ledger_report PRIVATE direct_eb_core)
target_compile_options(energy_ledger_report PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)

add_executable(lab_frame_energy_estimate
  utilities/lab_frame_energy_estimate.cpp
)
target_link_libraries(lab_frame_energy_estimate PRIVATE yaml-cpp)
target_compile_options(lab_frame_energy_estimate PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)

install(TARGETS simulator particle_text_to_hdf5 undulator_resonance
  energy_ledger_report lab_frame_energy_estimate
  RUNTIME DESTINATION bin)
if(HAVE_ELEGANT_SDDS_CONVERTER)
  install(TARGETS elegant_sdds_to_hdf5 RUNTIME DESTINATION bin)
endif()

if(BUILD_TESTING)
  add_subdirectory(verification)
endif()
