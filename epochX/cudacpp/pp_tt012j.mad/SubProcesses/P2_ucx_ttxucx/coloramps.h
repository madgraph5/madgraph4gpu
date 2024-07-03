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
    __device__ constexpr int channelId_to_iconfigC[8] = {
     0, // channelId=0: This value means not multi-channel, color will be wrong anyway -> pick the first
     0, // channelId=1 (diagram=1) --> iconfig=1 (f77 conv) and iconfigC=0 (c conv)
     1, // channelId=2 (diagram=2) --> iconfig=2 (f77 conv) and iconfigC=1 (c conv)
     2, // channelId=3 (diagram=3) --> iconfig=3 (f77 conv) and iconfigC=2 (c conv)
     3, // channelId=4 (diagram=4) --> iconfig=4 (f77 conv) and iconfigC=3 (c conv)
     4, // channelId=5 (diagram=5) --> iconfig=5 (f77 conv) and iconfigC=4 (c conv)
     5, // channelId=6 (diagram=6) --> iconfig=6 (f77 conv) and iconfigC=5 (c conv)
     6, // channelId=7 (diagram=7) --> iconfig=7 (f77 conv) and iconfigC=6 (c conv)
  };

  // Map iconfigC (in C indexing, i.e. iconfig-1) to the set of allowed colors
  // This array has N_config <= N_diagrams elements    
    __device__ constexpr bool icolamp[7][6] = {
    {false, false, true, false, false, false}, // iconfigC=0, diag=1
    {false, true, false, false, false, false}, // iconfigC=1, diag=2
    {false, true, true, false, false, false}, // iconfigC=2, diag=3
    {false, false, true, false, false, false}, // iconfigC=3, diag=4
    {false, true, false, false, false, false}, // iconfigC=4, diag=5
    {false, true, false, false, false, false}, // iconfigC=5, diag=6
    {false, false, true, false, false, false}, // iconfigC=6, diag=7
  };

}
#endif // COLORAMPS_H
