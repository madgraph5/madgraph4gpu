# Minimal CMake configuration to build a functional CPU version

cmake_minimum_required(VERSION 3.22)

project(Madgraph4GPU)

include(${PROJECT_SOURCE_DIR}/CMake/Platforms.txt)
include(${PROJECT_SOURCE_DIR}/CMake/Compilers.txt)
include(${PROJECT_SOURCE_DIR}/CMake/Macros.txt)

set(PROJECT_GITROOT_DIR ${PROJECT_SOURCE_DIR}/../../..)

add_subdirectory(src)
add_subdirectory(SubProcesses)
