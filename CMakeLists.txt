cmake_minimum_required(VERSION 3.16)

# APRL: Accelerator Particles and Radiation in Lorentz Frames.
project(APRL VERSION 2.1.1 LANGUAGES C CXX)

execute_process(COMMAND git rev-parse --verify HEAD
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  OUTPUT_VARIABLE APRL_SOURCE_REVISION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT APRL_SOURCE_REVISION)
  set(APRL_SOURCE_REVISION "unknown")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/SOURCE_REVISION")
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/SOURCE_REVISION" APRL_ARCHIVE_REVISION LIMIT_COUNT 1)
    string(LENGTH "${APRL_ARCHIVE_REVISION}" APRL_ARCHIVE_REVISION_LENGTH)
    if(APRL_ARCHIVE_REVISION MATCHES "^[0-9a-f]+$" AND APRL_ARCHIVE_REVISION_LENGTH EQUAL 40)
      set(APRL_SOURCE_REVISION "${APRL_ARCHIVE_REVISION}")
    endif()
  endif()
else()
  execute_process(COMMAND git status --porcelain --untracked-files=normal
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE APRL_SOURCE_DIRTY
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  if(APRL_SOURCE_DIRTY)
    string(APPEND APRL_SOURCE_REVISION "-dirty")
  endif()
endif()
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/SOURCE_REVISION" "${APRL_SOURCE_REVISION}\n")

include(CTest)

option(APRL_REQUIRE_PARALLEL_HDF5
  "Reject serial HDF5; required for multi-rank particle input" ON)
option(APRL_BUILD_POSTPROCESSORS
  "Build the five production radiation post-processors" ON)
option(APRL_BUILD_SASE_BUNCH_GENERATOR
  "Build the standalone SASE shot-noise bunch generator" ON)
option(APRL_PORTABLE_INSTALL
  "Use relative install RPATHs for libraries bundled beside APRL" OFF)
option(APRL_ENABLE_OPENMP
  "Enable optional compute threading (runtime default remains one thread)" ON)
option(APRL_ENABLE_MOVING_WINDOW
  "Build the moving window required by both APRL 2.0 coupling algorithms" ON)
if(NOT APRL_ENABLE_MOVING_WINDOW)
  message(FATAL_ERROR "APRL 2.0 requires APRL_ENABLE_MOVING_WINDOW=ON for its two supported algorithms")
endif()
option(DIRECT_EB_BOUNDS_CHECK
  "Check Yee lattice indices consistently in the core and its consumers" OFF)

if(APRL_PORTABLE_INSTALL)
  if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(FATAL_ERROR
      "APRL_PORTABLE_INSTALL currently supports Linux runtime loaders")
  endif()
  # Every installed executable is placed in <prefix>/bin while bootstrapped
  # libraries are placed in <prefix>/lib (or lib64).  Do not bake build-host
  # module paths into the installed copy; a site MPI remains supplied by the
  # active scheduler/module environment.
  set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib;$ORIGIN/../lib64")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
endif()

find_package(MPI REQUIRED COMPONENTS CXX)
set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 1.10 REQUIRED COMPONENTS C)
find_package(yaml-cpp 0.6 REQUIRED)
find_package(Threads REQUIRED)
if(APRL_ENABLE_OPENMP)
  find_package(OpenMP QUIET COMPONENTS CXX)
  if(NOT OpenMP_CXX_FOUND)
    message(STATUS "OpenMP unavailable: APRL accepts one compute thread only")
  endif()
endif()
find_package(ZLIB QUIET)
find_package(LibLZMA QUIET)
find_package(GSL QUIET)

# yaml-cpp <= 0.7 commonly exported an un-namespaced target, while 0.8
# exports yaml-cpp::yaml-cpp.  A plain "yaml-cpp" item degrades to
# "-lyaml-cpp" when only the namespaced target exists and loses the private
# installation prefix's library search path.
if(TARGET yaml-cpp::yaml-cpp)
  set(APRL_YAML_CPP_TARGET yaml-cpp::yaml-cpp)
elseif(TARGET yaml-cpp)
  set(APRL_YAML_CPP_TARGET yaml-cpp)
else()
  message(FATAL_ERROR "yaml-cpp was found but exported no supported CMake target")
endif()

set(APRL_BUILD_DESCRIPTION
  "${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION},${CMAKE_BUILD_TYPE},HDF5-${HDF5_VERSION}")

if(APRL_REQUIRE_PARALLEL_HDF5 AND NOT HDF5_IS_PARALLEL)
  message(FATAL_ERROR
    "A serial HDF5 installation was selected. Multi-rank particle input "
    "requires parallel HDF5 built against the same MPI family used here. "
    "Load the matching MPI/HDF5 modules or set HDF5_ROOT explicitly. "
    "Use -DAPRL_REQUIRE_PARALLEL_HDF5=OFF only for single-rank local work.")
elseif(NOT HDF5_IS_PARALLEL)
  message(WARNING
    "Serial HDF5 accepted for a single-rank/local build; multi-rank HDF5 "
    "particle input will fail at runtime.")
endif()

message(STATUS "Selected MPI C++ compiler: ${MPI_CXX_COMPILER}")
message(STATUS "Selected HDF5 ${HDF5_VERSION}; parallel=${HDF5_IS_PARALLEL}")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(DIRECT_EB_SOURCES
  src/boostframe.cpp
  src/compute_threads.cpp
  src/config.cpp
  src/checkpoint_io.cpp
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
  src/finite_radiation_screen.cpp
  src/radiation_screen_output.cpp
  src/longitudinal_mesh.cpp
  src/moving_window_fields.cpp
  src/particle_field_halo.cpp
  src/particle_file.cpp
  src/particle_initializer.cpp
  src/pending_particle_queue.cpp
  src/pending_particle_statistics.cpp
  src/particle_boundary.cpp
  src/runtime_control.cpp
  src/run_metadata.cpp
  src/runtime_util.cpp
  src/shutdown_watchdog.cpp
  src/simulation.cpp
  src/simulation_checkpoint.cpp
  src/simulation_moving_window.cpp
  src/simulation_energy_work.cpp
  src/source_coupling.cpp
  src/simulation_source_coupling.cpp
  src/simulation_particle_balance.cpp
  src/teardown_trace.cpp
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
    ${APRL_YAML_CPP_TARGET}
  PRIVATE Threads::Threads
)
if(APRL_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(direct_eb_core PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(direct_eb_core PUBLIC APRL_HAS_OPENMP=1)
endif()
if(DIRECT_EB_BOUNDS_CHECK)
  # Accessors are inline: every translation unit must use the same definition.
  target_compile_definitions(direct_eb_core PUBLIC DIRECT_EB_BOUNDS_CHECK=1)
endif()
if(APRL_ENABLE_MOVING_WINDOW)
  find_path(APRL_VERIFY_FFTW_INCLUDE_DIR fftw3.h)
  find_library(APRL_VERIFY_FFTW_LIBRARY NAMES fftw3)
  if(NOT APRL_VERIFY_FFTW_INCLUDE_DIR OR NOT APRL_VERIFY_FFTW_LIBRARY)
    message(FATAL_ERROR "APRL moving window requires FFTW3 headers and library")
  endif()
  target_sources(direct_eb_core PRIVATE src/blocked_cowan_reference.cpp src/transverse_reference_current.cpp)
  target_include_directories(direct_eb_core PRIVATE ${APRL_VERIFY_FFTW_INCLUDE_DIR})
  target_link_libraries(direct_eb_core PUBLIC ${APRL_VERIFY_FFTW_LIBRARY})
  target_compile_definitions(direct_eb_core PUBLIC APRL_HAS_MOVING_WINDOW=1)
endif()
target_compile_definitions(direct_eb_core
  PUBLIC OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX
)
set_source_files_properties(src/run_metadata.cpp PROPERTIES COMPILE_DEFINITIONS
  "APRL_SOURCE_REVISION=\"${APRL_SOURCE_REVISION}\";APRL_VERSION=\"${PROJECT_VERSION}\";APRL_BUILD_DESCRIPTION=\"${APRL_BUILD_DESCRIPTION}\"")
target_compile_options(direct_eb_core PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)

add_executable(aprl src/main.cpp)
target_compile_definitions(aprl PRIVATE APRL_VERSION="${PROJECT_VERSION}")
target_link_libraries(aprl PRIVATE direct_eb_core)

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

add_executable(uniform_lab_plane_bunch_generator
  utilities/uniform_lab_plane_bunch_generator.cpp
)
target_include_directories(uniform_lab_plane_bunch_generator PRIVATE
  ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(uniform_lab_plane_bunch_generator PRIVATE
  MPI::MPI_CXX
  ${HDF5_LIBRARIES}
)
target_compile_options(uniform_lab_plane_bunch_generator PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)
if(BUILD_TESTING)
  file(MAKE_DIRECTORY
    "${CMAKE_BINARY_DIR}/verification/work/uniform-generator")
  set(APRL_UNIFORM_GENERATOR_TEST_FILE
    "${CMAKE_BINARY_DIR}/verification/work/uniform-generator/particles.h5")
  add_test(NAME required.uniform-generator.generate
    COMMAND $<TARGET_FILE:uniform_lab_plane_bunch_generator>
      --output "${APRL_UNIFORM_GENERATOR_TEST_FILE}"
      --macro-weight-electrons 1000 --seed 17 --overwrite)
  set_tests_properties(required.uniform-generator.generate PROPERTIES
    LABELS "required;functional")
  add_test(NAME required.uniform-generator.verify
    COMMAND $<TARGET_FILE:uniform_lab_plane_bunch_generator>
      --verify "${APRL_UNIFORM_GENERATOR_TEST_FILE}")
  set_tests_properties(required.uniform-generator.verify PROPERTIES
    DEPENDS required.uniform-generator.generate
    LABELS "required;functional")
endif()

if(APRL_BUILD_SASE_BUNCH_GENERATOR)
  add_executable(sase_bunch_generator
    utilities/sase_bunch_generator.cpp
  )
  target_include_directories(sase_bunch_generator PRIVATE
    ${HDF5_INCLUDE_DIRS}
  )
  target_link_libraries(sase_bunch_generator PRIVATE
    MPI::MPI_CXX
    ${HDF5_LIBRARIES}
  )
  target_compile_options(sase_bunch_generator PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
  )
endif()

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
    MPI::MPI_CXX
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
target_link_libraries(undulator_resonance PRIVATE ${APRL_YAML_CPP_TARGET})
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
target_link_libraries(lab_frame_energy_estimate PRIVATE ${APRL_YAML_CPP_TARGET})
target_compile_options(lab_frame_energy_estimate PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)

install(TARGETS aprl particle_text_to_hdf5
  uniform_lab_plane_bunch_generator undulator_resonance
  energy_ledger_report lab_frame_energy_estimate
  RUNTIME DESTINATION bin)
if(APRL_BUILD_SASE_BUNCH_GENERATOR)
  install(TARGETS sase_bunch_generator RUNTIME DESTINATION bin)
endif()
if(HAVE_ELEGANT_SDDS_CONVERTER)
  install(TARGETS elegant_sdds_to_hdf5 RUNTIME DESTINATION bin)
endif()
install(PROGRAMS utilities/crop_elegant_hdf5.py DESTINATION bin)

install(FILES README.md README.zh-CN.md VERSION LICENSE DESTINATION share/aprl)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/SOURCE_REVISION" DESTINATION share/aprl)
install(DIRECTORY docs config examples
  DESTINATION share/aprl)
install(DIRECTORY postprocess/
  DESTINATION share/aprl/postprocess
  FILES_MATCHING PATTERN "example.yaml")

if(APRL_BUILD_POSTPROCESSORS)
  install(PROGRAMS postprocess/finite_screen_spectrum.py DESTINATION bin)
  add_subdirectory(postprocess/trajectory_radiation)
  add_subdirectory(postprocess/field_reconstruction)
  add_subdirectory(postprocess/field_plane_analysis)
  add_subdirectory(postprocess/field_power_compare)
  add_subdirectory(postprocess/energy_closure)
endif()

if(BUILD_TESTING)
  add_subdirectory(verification)
endif()
