file(GLOB_RECURSE HEADERS "*.h")
add_library(mg5amc_common Parameters_sm.cc read_slha.cc ${HEADERS})
