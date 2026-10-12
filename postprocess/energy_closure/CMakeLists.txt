cmake_minimum_required(VERSION 3.16)

project(energy_closure LANGUAGES CXX)

find_package(MPI REQUIRED COMPONENTS CXX)
set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 1.10 REQUIRED COMPONENTS C)
find_package(yaml-cpp 0.6 REQUIRED)

if(NOT APRL_YAML_CPP_TARGET)
  if(TARGET yaml-cpp::yaml-cpp)
    set(APRL_YAML_CPP_TARGET yaml-cpp::yaml-cpp)
  elseif(TARGET yaml-cpp)
    set(APRL_YAML_CPP_TARGET yaml-cpp)
  else()
    message(FATAL_ERROR
      "yaml-cpp was found but exported no supported CMake target")
  endif()
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(energy_closure main.cpp)
target_include_directories(energy_closure PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(energy_closure PRIVATE
  MPI::MPI_CXX
  ${HDF5_LIBRARIES}
  ${APRL_YAML_CPP_TARGET}
)
target_compile_definitions(energy_closure PRIVATE
  OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX
)
target_compile_options(energy_closure PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O3;-Wall;-Wextra;-Wpedantic>
)

install(TARGETS energy_closure RUNTIME DESTINATION bin)
