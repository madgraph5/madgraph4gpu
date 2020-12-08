cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_STANDARD 17)

project(check)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/SubProcesses/P1_Sigma_sm_epem_mupmum/)

set(ONEAPI_INSTALL_LOCATION /opt/intel/oneapi)

include_directories(
    src
    SubProcesses/P1_Sigma_sm_epem_mupmum
    ../../../tools
)

set (SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/src/read_slha.cc
    ${CMAKE_CURRENT_LIST_DIR}/src/Parameters_sm.cc
    ${CMAKE_CURRENT_LIST_DIR}/src/rambo.cc
    ${CMAKE_CURRENT_LIST_DIR}/SubProcesses/P1_Sigma_sm_epem_mupmum/CPPProcess.cc
    ${CMAKE_CURRENT_LIST_DIR}/SubProcesses/P1_Sigma_sm_epem_mupmum/check.cc
)
add_compile_definitions(MGONGPU_COMMONRAND_ONHOST)

set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

add_executable(check.exe ${SOURCES})
target_compile_options(check.exe PRIVATE -O0 -g -fsycl -fsycl-unnamed-lambda -ferror-limit=1 -Wall -pthread)
target_link_libraries(check.exe PRIVATE sycl -fsycl)


