file(GLOB_RECURSE HEADERS "*.h")
add_library(mg5amc_cxx BridgeKernels.cc CPPProcess.cc CrossSectionKernels.cc
                       MatrixElementKernels.cc RamboSamplingKernels.cc
                       RandomNumberKernels.cc ${HEADERS})
target_include_directories(mg5amc_cxx PUBLIC "${PROJECT_SOURCE_DIR}/src"
                                              "${PROJECT_GITROOT_DIR}/tools")

add_executable(check.exe check_sa.cc)
target_link_libraries(check.exe PUBLIC mg5amc_common mg5amc_cxx)
target_include_directories(check.exe PRIVATE "${PROJECT_SOURCE_DIR}/src")
