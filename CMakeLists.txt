cmake_minimum_required(VERSION 3.16)

# Working project identifier only. The program's public name is intentionally
# left undecided; the executable name "simulator" is also temporary.
project(unnamed_fel_program LANGUAGES CXX)

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(yaml-cpp REQUIRED)

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

add_executable(undulator_resonance utilities/undulator_resonance.cpp)
target_link_libraries(undulator_resonance PRIVATE yaml-cpp)
target_compile_options(undulator_resonance PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)

install(TARGETS simulator particle_text_to_hdf5 undulator_resonance
  RUNTIME DESTINATION bin)
