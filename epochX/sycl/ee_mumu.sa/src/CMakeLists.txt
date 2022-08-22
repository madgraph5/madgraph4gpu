file(GLOB_RECURSE HEADERS "*.h")
add_library(mg5amc_common Parameters_sm.cc read_slha.cc ${HEADERS})

# some XCode specific stuff to make the executable run
set_property(TARGET mg5amc_common PROPERTY XCODE_GENERATE_SCHEME TRUE)
