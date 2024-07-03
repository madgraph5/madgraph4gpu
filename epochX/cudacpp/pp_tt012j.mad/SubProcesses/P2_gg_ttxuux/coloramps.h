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
    __device__ constexpr int channelId_to_iconfigC[37] = {
     0, // channelId=0: This value means not multi-channel, color will be wrong anyway -> pick the first
     0, // channelId=1 (diagram=1) --> iconfig=1 (f77 conv) and iconfigC=0 (c conv)
     1, // channelId=2 (diagram=2) --> iconfig=2 (f77 conv) and iconfigC=1 (c conv)
     2, // channelId=3 (diagram=3) --> iconfig=3 (f77 conv) and iconfigC=2 (c conv)
     3, // channelId=4 (diagram=4) --> iconfig=4 (f77 conv) and iconfigC=3 (c conv)
     4, // channelId=5 (diagram=5) --> iconfig=5 (f77 conv) and iconfigC=4 (c conv)
     5, // channelId=6 (diagram=6) --> iconfig=6 (f77 conv) and iconfigC=5 (c conv)
     6, // channelId=7 (diagram=7) --> iconfig=7 (f77 conv) and iconfigC=6 (c conv)
     7, // channelId=8 (diagram=8) --> iconfig=8 (f77 conv) and iconfigC=7 (c conv)
     8, // channelId=9 (diagram=9) --> iconfig=9 (f77 conv) and iconfigC=8 (c conv)
     9, // channelId=10 (diagram=10) --> iconfig=10 (f77 conv) and iconfigC=9 (c conv)
     10, // channelId=11 (diagram=11) --> iconfig=11 (f77 conv) and iconfigC=10 (c conv)
     11, // channelId=12 (diagram=12) --> iconfig=12 (f77 conv) and iconfigC=11 (c conv)
     12, // channelId=13 (diagram=13) --> iconfig=13 (f77 conv) and iconfigC=12 (c conv)
     13, // channelId=14 (diagram=14) --> iconfig=14 (f77 conv) and iconfigC=13 (c conv)
     14, // channelId=15 (diagram=15) --> iconfig=15 (f77 conv) and iconfigC=14 (c conv)
     15, // channelId=16 (diagram=16) --> iconfig=16 (f77 conv) and iconfigC=15 (c conv)
     16, // channelId=17 (diagram=17) --> iconfig=17 (f77 conv) and iconfigC=16 (c conv)
     17, // channelId=18 (diagram=18) --> iconfig=18 (f77 conv) and iconfigC=17 (c conv)
     18, // channelId=19 (diagram=19) --> iconfig=19 (f77 conv) and iconfigC=18 (c conv)
     19, // channelId=20 (diagram=20) --> iconfig=20 (f77 conv) and iconfigC=19 (c conv)
     20, // channelId=21 (diagram=21) --> iconfig=21 (f77 conv) and iconfigC=20 (c conv)
     21, // channelId=22 (diagram=22) --> iconfig=22 (f77 conv) and iconfigC=21 (c conv)
     22, // channelId=23 (diagram=23) --> iconfig=23 (f77 conv) and iconfigC=22 (c conv)
     23, // channelId=24 (diagram=24) --> iconfig=24 (f77 conv) and iconfigC=23 (c conv)
     24, // channelId=25 (diagram=25) --> iconfig=25 (f77 conv) and iconfigC=24 (c conv)
     25, // channelId=26 (diagram=26) --> iconfig=26 (f77 conv) and iconfigC=25 (c conv)
     26, // channelId=27 (diagram=27) --> iconfig=27 (f77 conv) and iconfigC=26 (c conv)
     27, // channelId=28 (diagram=28) --> iconfig=28 (f77 conv) and iconfigC=27 (c conv)
     28, // channelId=29 (diagram=29) --> iconfig=29 (f77 conv) and iconfigC=28 (c conv)
     29, // channelId=30 (diagram=30) --> iconfig=30 (f77 conv) and iconfigC=29 (c conv)
     30, // channelId=31 (diagram=31) --> iconfig=31 (f77 conv) and iconfigC=30 (c conv)
     31, // channelId=32 (diagram=32) --> iconfig=32 (f77 conv) and iconfigC=31 (c conv)
     32, // channelId=33 (diagram=33) --> iconfig=33 (f77 conv) and iconfigC=32 (c conv)
     -1, // channelId=34 (diagram=34): Not consider as a channel of integration (presence of 4 point interaction?)
     33, // channelId=35 (diagram=35) --> iconfig=34 (f77 conv) and iconfigC=33 (c conv)
     34, // channelId=36 (diagram=36) --> iconfig=35 (f77 conv) and iconfigC=34 (c conv)
  };

  // Map iconfigC (in C indexing, i.e. iconfig-1) to the set of allowed colors
  // This array has N_config <= N_diagrams elements    
    __device__ constexpr bool icolamp[35][12] = {
    {false, false, true, false, false, false, false, false, false, false, true, false}, // iconfigC=0, diag=1
    {false, true, false, false, false, false, false, false, false, true, false, false}, // iconfigC=1, diag=2
    {false, true, true, false, false, false, false, false, false, true, true, false}, // iconfigC=2, diag=3
    {false, true, false, false, false, false, false, false, false, true, false, false}, // iconfigC=3, diag=4
    {false, false, true, false, false, false, false, false, false, false, true, false}, // iconfigC=4, diag=5
    {false, false, false, false, false, true, false, false, false, false, false, false}, // iconfigC=5, diag=6
    {false, false, false, false, false, true, false, false, false, false, false, false}, // iconfigC=6, diag=7
    {false, true, false, false, false, false, false, false, false, false, false, false}, // iconfigC=7, diag=8
    {false, true, false, false, false, false, false, false, false, false, false, false}, // iconfigC=8, diag=9
    {false, true, false, false, false, true, false, false, false, false, false, false}, // iconfigC=9, diag=10
    {false, false, false, false, false, false, true, false, false, false, false, false}, // iconfigC=10, diag=11
    {false, false, false, false, false, false, false, false, false, false, true, false}, // iconfigC=11, diag=12
    {false, false, false, false, false, false, true, false, false, false, false, false}, // iconfigC=12, diag=13
    {false, false, false, false, false, false, false, false, false, false, true, false}, // iconfigC=13, diag=14
    {false, false, false, false, false, false, true, false, false, false, true, false}, // iconfigC=14, diag=15
    {false, false, false, false, false, false, true, false, false, false, false, false}, // iconfigC=15, diag=16
    {false, false, true, false, false, false, false, false, false, false, false, false}, // iconfigC=16, diag=17
    {false, false, false, false, false, false, true, false, false, false, false, false}, // iconfigC=17, diag=18
    {false, false, true, false, false, false, false, false, false, false, false, false}, // iconfigC=18, diag=19
    {false, false, true, false, false, false, true, false, false, false, false, false}, // iconfigC=19, diag=20
    {false, false, false, false, false, false, false, false, false, true, false, false}, // iconfigC=20, diag=21
    {false, false, false, false, false, true, false, false, false, false, false, false}, // iconfigC=21, diag=22
    {false, false, false, false, false, true, false, false, false, false, false, false}, // iconfigC=22, diag=23
    {false, false, false, false, false, false, false, false, false, true, false, false}, // iconfigC=23, diag=24
    {false, false, false, false, false, true, false, false, false, true, false, false}, // iconfigC=24, diag=25
    {false, false, false, false, false, false, false, false, false, true, false, false}, // iconfigC=25, diag=26
    {false, false, false, false, false, false, true, false, false, true, false, false}, // iconfigC=26, diag=27
    {false, false, true, false, false, false, false, false, false, false, false, false}, // iconfigC=27, diag=28
    {false, false, true, false, false, true, false, false, false, false, false, false}, // iconfigC=28, diag=29
    {false, false, false, false, false, false, false, false, false, false, true, false}, // iconfigC=29, diag=30
    {false, false, false, false, false, true, false, false, false, false, true, false}, // iconfigC=30, diag=31
    {false, true, false, false, false, false, false, false, false, false, false, false}, // iconfigC=31, diag=32
    {false, true, false, false, false, false, true, false, false, false, false, false}, // iconfigC=32, diag=33
    {false, true, true, false, false, true, true, false, false, true, true, false}, // iconfigC=33, diag=35
    {false, true, false, false, false, true, true, false, false, false, true, false}, // iconfigC=34, diag=36
  };

}
#endif // COLORAMPS_H
