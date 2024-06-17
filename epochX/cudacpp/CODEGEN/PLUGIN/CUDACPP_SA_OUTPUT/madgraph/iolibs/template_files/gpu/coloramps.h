// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#ifndef COLORAMPS_H
#define COLORAMPS_H 1

#include <map>

namespace mgOnGpu
{
  // Summary of numbering and indexing conventions for the relevant concepts (see issue #826 and PR #852)
  // - Diagram number (no variable) in [1, N_diagrams]: all values are allowed (N_diagrams distinct values)
  //   => this number is displayed for information before each block of code in CPPProcess.cc
  // - Channel number (CHANNEL_ID) in [0, N_diagrams]: not all values are allowed (N_config <= N_diagrams distinct values)
  //   => this number (with indexing like ps/pdf output) is passed around as an API argument between cudacpp functions
  //   0 is allowed to fallback to no multi-channel mode.
  // - Channel number in C indexing: "IconfiC", this is the equivalent of the Fortran iconfig
  //   iconfigC = iconfig -1
  //   provides a continuous index [0, N_config-1] for array
  //  iconfigC = ChannelId_to_iconfigC[channelId]
  //NOTE: All those ordering are event by event specific (with the intent to have those fix within a vector size/wrap   
  
  // Map channelId to iconfigC
  // This array has N_diagrams+1 elements, but only N_config <= N_diagrams valid values
  // unvalid values are set to -1
  // The 0 entry is a fall back to still write events even if no multi-channel is setup (wrong color selected in that mode) 
    __device__ constexpr int channelId_to_iconfigC[%(nb_diag_plus_one)i] = {
     0, // channelId=0: This value means not multi-channel, color will be wrong anyway -> pick the first
%(diag_to_channel)s
  };

  // Map iconfigC (in C indexing, i.e. iconfig-1) to the set of allowed colors
  // This array has N_config <= N_diagrams elements    
    __device__ constexpr bool icolamp[%(nb_channel)s][%(nb_color)s] = {
%(is_LC)s
  };

}
#endif // COLORAMPS_H
