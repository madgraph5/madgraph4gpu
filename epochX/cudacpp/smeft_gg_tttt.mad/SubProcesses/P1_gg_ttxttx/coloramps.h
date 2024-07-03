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
    __device__ constexpr int channelId_to_iconfigC[73] = {
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
     33, // channelId=34 (diagram=34) --> iconfig=34 (f77 conv) and iconfigC=33 (c conv)
     34, // channelId=35 (diagram=35) --> iconfig=35 (f77 conv) and iconfigC=34 (c conv)
     35, // channelId=36 (diagram=36) --> iconfig=36 (f77 conv) and iconfigC=35 (c conv)
     36, // channelId=37 (diagram=37) --> iconfig=37 (f77 conv) and iconfigC=36 (c conv)
     37, // channelId=38 (diagram=38) --> iconfig=38 (f77 conv) and iconfigC=37 (c conv)
     38, // channelId=39 (diagram=39) --> iconfig=39 (f77 conv) and iconfigC=38 (c conv)
     39, // channelId=40 (diagram=40) --> iconfig=40 (f77 conv) and iconfigC=39 (c conv)
     40, // channelId=41 (diagram=41) --> iconfig=41 (f77 conv) and iconfigC=40 (c conv)
     41, // channelId=42 (diagram=42) --> iconfig=42 (f77 conv) and iconfigC=41 (c conv)
     42, // channelId=43 (diagram=43) --> iconfig=43 (f77 conv) and iconfigC=42 (c conv)
     43, // channelId=44 (diagram=44) --> iconfig=44 (f77 conv) and iconfigC=43 (c conv)
     44, // channelId=45 (diagram=45) --> iconfig=45 (f77 conv) and iconfigC=44 (c conv)
     45, // channelId=46 (diagram=46) --> iconfig=46 (f77 conv) and iconfigC=45 (c conv)
     46, // channelId=47 (diagram=47) --> iconfig=47 (f77 conv) and iconfigC=46 (c conv)
     47, // channelId=48 (diagram=48) --> iconfig=48 (f77 conv) and iconfigC=47 (c conv)
     48, // channelId=49 (diagram=49) --> iconfig=49 (f77 conv) and iconfigC=48 (c conv)
     49, // channelId=50 (diagram=50) --> iconfig=50 (f77 conv) and iconfigC=49 (c conv)
     50, // channelId=51 (diagram=51) --> iconfig=51 (f77 conv) and iconfigC=50 (c conv)
     51, // channelId=52 (diagram=52) --> iconfig=52 (f77 conv) and iconfigC=51 (c conv)
     52, // channelId=53 (diagram=53) --> iconfig=53 (f77 conv) and iconfigC=52 (c conv)
     53, // channelId=54 (diagram=54) --> iconfig=54 (f77 conv) and iconfigC=53 (c conv)
     54, // channelId=55 (diagram=55) --> iconfig=55 (f77 conv) and iconfigC=54 (c conv)
     55, // channelId=56 (diagram=56) --> iconfig=56 (f77 conv) and iconfigC=55 (c conv)
     56, // channelId=57 (diagram=57) --> iconfig=57 (f77 conv) and iconfigC=56 (c conv)
     57, // channelId=58 (diagram=58) --> iconfig=58 (f77 conv) and iconfigC=57 (c conv)
     58, // channelId=59 (diagram=59) --> iconfig=59 (f77 conv) and iconfigC=58 (c conv)
     59, // channelId=60 (diagram=60) --> iconfig=60 (f77 conv) and iconfigC=59 (c conv)
     60, // channelId=61 (diagram=61) --> iconfig=61 (f77 conv) and iconfigC=60 (c conv)
     61, // channelId=62 (diagram=62) --> iconfig=62 (f77 conv) and iconfigC=61 (c conv)
     62, // channelId=63 (diagram=63) --> iconfig=63 (f77 conv) and iconfigC=62 (c conv)
     63, // channelId=64 (diagram=64) --> iconfig=64 (f77 conv) and iconfigC=63 (c conv)
     64, // channelId=65 (diagram=65) --> iconfig=65 (f77 conv) and iconfigC=64 (c conv)
     65, // channelId=66 (diagram=66) --> iconfig=66 (f77 conv) and iconfigC=65 (c conv)
     -1, // channelId=67 (diagram=67): Not consider as a channel of integration (presence of 4 point interaction?)
     66, // channelId=68 (diagram=68) --> iconfig=67 (f77 conv) and iconfigC=66 (c conv)
     67, // channelId=69 (diagram=69) --> iconfig=68 (f77 conv) and iconfigC=67 (c conv)
     -1, // channelId=70 (diagram=70): Not consider as a channel of integration (presence of 4 point interaction?)
     68, // channelId=71 (diagram=71) --> iconfig=69 (f77 conv) and iconfigC=68 (c conv)
     69, // channelId=72 (diagram=72) --> iconfig=70 (f77 conv) and iconfigC=69 (c conv)
  };

  // Map iconfigC (in C indexing, i.e. iconfig-1) to the set of allowed colors
  // This array has N_config <= N_diagrams elements    
    __device__ constexpr bool icolamp[70][12] = {
    {false, false, true, false, false, false, false, false, false, false, true, false}, // iconfigC=0, diag=1
    {false, true, false, false, false, false, false, false, false, true, false, false}, // iconfigC=1, diag=2
    {false, true, true, false, false, false, false, false, false, true, true, false}, // iconfigC=2, diag=3
    {true, false, false, false, false, false, false, false, true, false, false, false}, // iconfigC=3, diag=4
    {false, false, false, true, false, false, false, false, false, false, false, true}, // iconfigC=4, diag=5
    {true, false, false, true, false, false, false, false, true, false, false, true}, // iconfigC=5, diag=6
    {true, false, false, false, false, false, false, false, true, false, false, false}, // iconfigC=6, diag=7
    {false, false, false, true, false, false, false, false, false, false, false, true}, // iconfigC=7, diag=8
    {false, true, false, false, false, false, false, false, false, true, false, false}, // iconfigC=8, diag=9
    {false, false, true, false, false, false, false, false, false, false, true, false}, // iconfigC=9, diag=10
    {true, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=10, diag=11
    {false, false, false, false, false, true, false, false, false, false, false, false}, // iconfigC=11, diag=12
    {false, false, false, false, false, true, false, false, false, false, false, false}, // iconfigC=12, diag=13
    {false, false, false, false, true, false, false, false, false, false, false, false}, // iconfigC=13, diag=14
    {false, true, false, false, false, false, false, false, false, false, false, false}, // iconfigC=14, diag=15
    {false, false, false, false, true, false, false, false, false, false, false, false}, // iconfigC=15, diag=16
    {true, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=16, diag=17
    {true, false, false, false, true, false, false, false, false, false, false, false}, // iconfigC=17, diag=18
    {false, true, false, false, false, false, false, false, false, false, false, false}, // iconfigC=18, diag=19
    {false, true, false, false, false, true, false, false, false, false, false, false}, // iconfigC=19, diag=20
    {false, false, false, false, false, false, false, false, true, false, false, false}, // iconfigC=20, diag=21
    {false, false, false, false, false, false, true, false, false, false, false, false}, // iconfigC=21, diag=22
    {false, false, false, false, false, false, false, false, false, false, true, false}, // iconfigC=22, diag=23
    {false, false, false, false, true, false, false, false, false, false, false, false}, // iconfigC=23, diag=24
    {false, false, false, false, false, false, true, false, false, false, false, false}, // iconfigC=24, diag=25
    {false, false, false, false, true, false, false, false, false, false, false, false}, // iconfigC=25, diag=26
    {false, false, false, false, false, false, false, false, true, false, false, false}, // iconfigC=26, diag=27
    {false, false, false, false, true, false, false, false, true, false, false, false}, // iconfigC=27, diag=28
    {false, false, false, false, false, false, false, false, false, false, true, false}, // iconfigC=28, diag=29
    {false, false, false, false, false, false, true, false, false, false, true, false}, // iconfigC=29, diag=30
    {false, false, false, false, false, false, false, true, false, false, false, false}, // iconfigC=30, diag=31
    {false, false, false, false, false, false, true, false, false, false, false, false}, // iconfigC=31, diag=32
    {false, false, true, false, false, false, false, false, false, false, false, false}, // iconfigC=32, diag=33
    {false, false, false, false, false, false, false, true, false, false, false, false}, // iconfigC=33, diag=34
    {false, false, false, true, false, false, false, false, false, false, false, false}, // iconfigC=34, diag=35
    {false, false, false, false, false, false, true, false, false, false, false, false}, // iconfigC=35, diag=36
    {false, false, true, false, false, false, false, false, false, false, false, false}, // iconfigC=36, diag=37
    {false, false, true, false, false, false, true, false, false, false, false, false}, // iconfigC=37, diag=38
    {false, false, false, true, false, false, false, false, false, false, false, false}, // iconfigC=38, diag=39
    {false, false, false, true, false, false, false, true, false, false, false, false}, // iconfigC=39, diag=40
    {false, false, false, false, false, false, false, false, false, true, false, false}, // iconfigC=40, diag=41
    {false, false, false, false, false, false, false, true, false, false, false, false}, // iconfigC=41, diag=42
    {false, false, false, false, false, false, false, true, false, false, false, false}, // iconfigC=42, diag=43
    {false, false, false, false, false, true, false, false, false, false, false, false}, // iconfigC=43, diag=44
    {false, false, false, false, false, false, false, false, false, false, false, true}, // iconfigC=44, diag=45
    {false, false, false, false, false, true, false, false, false, false, false, false}, // iconfigC=45, diag=46
    {false, false, false, false, false, false, false, false, false, true, false, false}, // iconfigC=46, diag=47
    {false, false, false, false, false, true, false, false, false, true, false, false}, // iconfigC=47, diag=48
    {false, false, false, false, false, false, false, false, false, false, false, true}, // iconfigC=48, diag=49
    {false, false, false, false, false, false, false, true, false, false, false, true}, // iconfigC=49, diag=50
    {false, false, false, false, false, false, false, false, true, false, false, false}, // iconfigC=50, diag=51
    {false, false, false, false, false, false, false, true, true, false, false, false}, // iconfigC=51, diag=52
    {false, false, false, false, false, false, false, false, false, true, false, false}, // iconfigC=52, diag=53
    {false, false, false, false, false, false, true, false, false, true, false, false}, // iconfigC=53, diag=54
    {true, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=54, diag=55
    {true, false, false, false, false, false, false, true, false, false, false, false}, // iconfigC=55, diag=56
    {false, false, true, false, false, false, false, false, false, false, false, false}, // iconfigC=56, diag=57
    {false, false, true, false, false, true, false, false, false, false, false, false}, // iconfigC=57, diag=58
    {false, false, false, false, false, false, false, false, false, false, true, false}, // iconfigC=58, diag=59
    {false, false, false, false, false, true, false, false, false, false, true, false}, // iconfigC=59, diag=60
    {false, false, false, false, false, false, false, false, false, false, false, true}, // iconfigC=60, diag=61
    {false, false, false, false, true, false, false, false, false, false, false, true}, // iconfigC=61, diag=62
    {false, true, false, false, false, false, false, false, false, false, false, false}, // iconfigC=62, diag=63
    {false, true, false, false, false, false, true, false, false, false, false, false}, // iconfigC=63, diag=64
    {false, false, false, true, false, false, false, false, false, false, false, false}, // iconfigC=64, diag=65
    {false, false, false, true, true, false, false, false, false, false, false, false}, // iconfigC=65, diag=66
    {false, true, true, false, false, true, true, false, false, true, true, false}, // iconfigC=66, diag=68
    {false, true, false, false, false, true, true, false, false, false, true, false}, // iconfigC=67, diag=69
    {false, false, true, false, false, true, true, false, false, true, false, false}, // iconfigC=68, diag=71
    {true, false, false, true, true, false, false, true, true, false, false, true}, // iconfigC=69, diag=72
  };

}
#endif // COLORAMPS_H
