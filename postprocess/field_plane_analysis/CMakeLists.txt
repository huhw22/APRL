cmake_minimum_required(VERSION 3.16)

project(field_plane_analysis LANGUAGES CXX)

find_package(MPI REQUIRED COMPONENTS CXX)
set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 1.10 REQUIRED COMPONENTS C)
find_package(yaml-cpp 0.6 REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)
find_library(FFTW3_THREADS_LIBRARY NAMES fftw3_threads REQUIRED)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(field_plane_analysis
  configuration.cpp
  field_plane_analysis.cpp
  main.cpp
)

target_include_directories(field_plane_analysis PRIVATE
  ${HDF5_INCLUDE_DIRS}
  ${FFTW3_INCLUDE_DIRS}
)
target_link_libraries(field_plane_analysis PRIVATE
  MPI::MPI_CXX
  ${HDF5_LIBRARIES}
  yaml-cpp
  ${FFTW3_THREADS_LIBRARY}
  ${FFTW3_LIBRARIES}
  Threads::Threads
)
target_compile_definitions(field_plane_analysis PRIVATE
  OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX
)
target_compile_options(field_plane_analysis PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)

install(TARGETS field_plane_analysis RUNTIME DESTINATION bin)
