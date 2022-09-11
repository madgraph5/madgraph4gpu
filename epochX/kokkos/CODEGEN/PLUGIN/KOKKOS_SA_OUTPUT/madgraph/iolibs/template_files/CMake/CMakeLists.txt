# Minimal CMake configuration to build a functional CPU version

cmake_minimum_required(VERSION 3.22)

project(Madgraph4GPU)

find_package(Kokkos REQUIRED)
# set cuda flags
if(DEFINED Kokkos_ENABLE_CUDA)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --expt-relaxed-constexpr")
endif(DEFINED Kokkos_ENABLE_CUDA)

include(${PROJECT_SOURCE_DIR}/CMake/Platforms.txt)
include(${PROJECT_SOURCE_DIR}/CMake/Compilers.txt)
include(${PROJECT_SOURCE_DIR}/CMake/Macros.txt)

find_package(Kokkos REQUIRED)

# set cuda flags
if(DEFINED Kokkos_ENABLE_CUDA)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --expt-relaxed-constexpr")
endif(DEFINED Kokkos_ENABLE_CUDA)


set(PROJECT_GITROOT_DIR ${PROJECT_SOURCE_DIR}/../../..)

add_subdirectory(src)
add_subdirectory(SubProcesses)
