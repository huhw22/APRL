cmake_minimum_required(VERSION 3.16)

project(trajectory_radiation LANGUAGES CXX)

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(yaml-cpp REQUIRED)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(trajectory_radiation
  configuration.cpp
  main.cpp
  output_writer.cpp
  radiation_kernel.cpp
  trajectory_reader.cpp
)

target_include_directories(trajectory_radiation PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(trajectory_radiation PRIVATE
  MPI::MPI_CXX
  ${HDF5_LIBRARIES}
  yaml-cpp
)
target_compile_definitions(trajectory_radiation PRIVATE
  OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX
)
target_compile_options(trajectory_radiation PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)

install(TARGETS trajectory_radiation RUNTIME DESTINATION bin)
