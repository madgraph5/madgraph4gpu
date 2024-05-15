# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: S. Roiser (Feb 2022) for the MG5aMC CUDACPP plugin.
# Further modified by: S. Roiser (2022-2024) for the MG5aMC CUDACPP plugin.

file(GLOB_RECURSE HEADERS "*.h")
add_library(mg5amc_common Parameters_sm.cc read_slha.cc ${HEADERS})

# some XCode specific stuff to make the executable run
set_property(TARGET mg5amc_common PROPERTY XCODE_GENERATE_SCHEME TRUE)
