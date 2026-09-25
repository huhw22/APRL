cmake_minimum_required(VERSION 3.16)

project(field_power_compare LANGUAGES CXX)

find_package(MPI REQUIRED COMPONENTS CXX)
set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 1.10 REQUIRED COMPONENTS C)
find_package(yaml-cpp 0.6 REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(field_power_compare main.cpp)
target_include_directories(field_power_compare PRIVATE
  ${HDF5_INCLUDE_DIRS}
  ${FFTW3_INCLUDE_DIRS}
)
target_link_libraries(field_power_compare PRIVATE
  MPI::MPI_CXX
  ${HDF5_LIBRARIES}
  yaml-cpp
  ${FFTW3_LIBRARIES}
)
target_compile_definitions(field_power_compare PRIVATE
  OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX
)
target_compile_options(field_power_compare PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)

install(TARGETS field_power_compare RUNTIME DESTINATION bin)
