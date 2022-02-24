get_filename_component(basename ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(TOLOWER ${basename} targadd)

file(GLOB_RECURSE HEADERS "*.h")

set(libname mg5amc_cxx_${targadd})
add_library(${libname} BridgeKernels.cc CPPProcess.cc CrossSectionKernels.cc MatrixElementKernels.cc
                       RamboSamplingKernels.cc RandomNumberKernels.cc ${HEADERS})
target_include_directories(${libname} PUBLIC "${PROJECT_SOURCE_DIR}/src"
                                                        "${PROJECT_GITROOT_DIR}/tools")

set(execname check_${targadd}.exe)
add_executable(${execname} check_sa.cc)
target_link_libraries(${execname} PUBLIC mg5amc_common ${libname})
target_include_directories(${execname} PRIVATE "${PROJECT_SOURCE_DIR}/src")
