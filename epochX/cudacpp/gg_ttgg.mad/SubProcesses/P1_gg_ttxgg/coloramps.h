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
    __device__ constexpr int channelId_to_iconfigC[114] = {
     0, // channelId=0: This value means not multi-channel, color will be wrong anyway -> pick the first
     -1, // channelId=1 (diagram=1): Not consider as a channel of integration (presence of 4 point interaction?)
     0, // channelId=2 (diagram=2) --> iconfig=1 (f77 conv) and iconfigC=0 (c conv)
     1, // channelId=3 (diagram=3) --> iconfig=2 (f77 conv) and iconfigC=1 (c conv)
     2, // channelId=4 (diagram=4) --> iconfig=3 (f77 conv) and iconfigC=2 (c conv)
     3, // channelId=5 (diagram=5) --> iconfig=4 (f77 conv) and iconfigC=3 (c conv)
     4, // channelId=6 (diagram=6) --> iconfig=5 (f77 conv) and iconfigC=4 (c conv)
     5, // channelId=7 (diagram=7) --> iconfig=6 (f77 conv) and iconfigC=5 (c conv)
     6, // channelId=8 (diagram=8) --> iconfig=7 (f77 conv) and iconfigC=6 (c conv)
     7, // channelId=9 (diagram=9) --> iconfig=8 (f77 conv) and iconfigC=7 (c conv)
     8, // channelId=10 (diagram=10) --> iconfig=9 (f77 conv) and iconfigC=8 (c conv)
     9, // channelId=11 (diagram=11) --> iconfig=10 (f77 conv) and iconfigC=9 (c conv)
     10, // channelId=12 (diagram=12) --> iconfig=11 (f77 conv) and iconfigC=10 (c conv)
     11, // channelId=13 (diagram=13) --> iconfig=12 (f77 conv) and iconfigC=11 (c conv)
     12, // channelId=14 (diagram=14) --> iconfig=13 (f77 conv) and iconfigC=12 (c conv)
     13, // channelId=15 (diagram=15) --> iconfig=14 (f77 conv) and iconfigC=13 (c conv)
     14, // channelId=16 (diagram=16) --> iconfig=15 (f77 conv) and iconfigC=14 (c conv)
     15, // channelId=17 (diagram=17) --> iconfig=16 (f77 conv) and iconfigC=15 (c conv)
     16, // channelId=18 (diagram=18) --> iconfig=17 (f77 conv) and iconfigC=16 (c conv)
     17, // channelId=19 (diagram=19) --> iconfig=18 (f77 conv) and iconfigC=17 (c conv)
     18, // channelId=20 (diagram=20) --> iconfig=19 (f77 conv) and iconfigC=18 (c conv)
     19, // channelId=21 (diagram=21) --> iconfig=20 (f77 conv) and iconfigC=19 (c conv)
     20, // channelId=22 (diagram=22) --> iconfig=21 (f77 conv) and iconfigC=20 (c conv)
     21, // channelId=23 (diagram=23) --> iconfig=22 (f77 conv) and iconfigC=21 (c conv)
     22, // channelId=24 (diagram=24) --> iconfig=23 (f77 conv) and iconfigC=22 (c conv)
     23, // channelId=25 (diagram=25) --> iconfig=24 (f77 conv) and iconfigC=23 (c conv)
     24, // channelId=26 (diagram=26) --> iconfig=25 (f77 conv) and iconfigC=24 (c conv)
     25, // channelId=27 (diagram=27) --> iconfig=26 (f77 conv) and iconfigC=25 (c conv)
     26, // channelId=28 (diagram=28) --> iconfig=27 (f77 conv) and iconfigC=26 (c conv)
     27, // channelId=29 (diagram=29) --> iconfig=28 (f77 conv) and iconfigC=27 (c conv)
     28, // channelId=30 (diagram=30) --> iconfig=29 (f77 conv) and iconfigC=28 (c conv)
     29, // channelId=31 (diagram=31) --> iconfig=30 (f77 conv) and iconfigC=29 (c conv)
     -1, // channelId=32 (diagram=32): Not consider as a channel of integration (presence of 4 point interaction?)
     30, // channelId=33 (diagram=33) --> iconfig=31 (f77 conv) and iconfigC=30 (c conv)
     31, // channelId=34 (diagram=34) --> iconfig=32 (f77 conv) and iconfigC=31 (c conv)
     32, // channelId=35 (diagram=35) --> iconfig=33 (f77 conv) and iconfigC=32 (c conv)
     33, // channelId=36 (diagram=36) --> iconfig=34 (f77 conv) and iconfigC=33 (c conv)
     34, // channelId=37 (diagram=37) --> iconfig=35 (f77 conv) and iconfigC=34 (c conv)
     35, // channelId=38 (diagram=38) --> iconfig=36 (f77 conv) and iconfigC=35 (c conv)
     36, // channelId=39 (diagram=39) --> iconfig=37 (f77 conv) and iconfigC=36 (c conv)
     37, // channelId=40 (diagram=40) --> iconfig=38 (f77 conv) and iconfigC=37 (c conv)
     38, // channelId=41 (diagram=41) --> iconfig=39 (f77 conv) and iconfigC=38 (c conv)
     39, // channelId=42 (diagram=42) --> iconfig=40 (f77 conv) and iconfigC=39 (c conv)
     40, // channelId=43 (diagram=43) --> iconfig=41 (f77 conv) and iconfigC=40 (c conv)
     41, // channelId=44 (diagram=44) --> iconfig=42 (f77 conv) and iconfigC=41 (c conv)
     42, // channelId=45 (diagram=45) --> iconfig=43 (f77 conv) and iconfigC=42 (c conv)
     43, // channelId=46 (diagram=46) --> iconfig=44 (f77 conv) and iconfigC=43 (c conv)
     44, // channelId=47 (diagram=47) --> iconfig=45 (f77 conv) and iconfigC=44 (c conv)
     -1, // channelId=48 (diagram=48): Not consider as a channel of integration (presence of 4 point interaction?)
     45, // channelId=49 (diagram=49) --> iconfig=46 (f77 conv) and iconfigC=45 (c conv)
     46, // channelId=50 (diagram=50) --> iconfig=47 (f77 conv) and iconfigC=46 (c conv)
     47, // channelId=51 (diagram=51) --> iconfig=48 (f77 conv) and iconfigC=47 (c conv)
     48, // channelId=52 (diagram=52) --> iconfig=49 (f77 conv) and iconfigC=48 (c conv)
     49, // channelId=53 (diagram=53) --> iconfig=50 (f77 conv) and iconfigC=49 (c conv)
     50, // channelId=54 (diagram=54) --> iconfig=51 (f77 conv) and iconfigC=50 (c conv)
     51, // channelId=55 (diagram=55) --> iconfig=52 (f77 conv) and iconfigC=51 (c conv)
     52, // channelId=56 (diagram=56) --> iconfig=53 (f77 conv) and iconfigC=52 (c conv)
     53, // channelId=57 (diagram=57) --> iconfig=54 (f77 conv) and iconfigC=53 (c conv)
     -1, // channelId=58 (diagram=58): Not consider as a channel of integration (presence of 4 point interaction?)
     54, // channelId=59 (diagram=59) --> iconfig=55 (f77 conv) and iconfigC=54 (c conv)
     55, // channelId=60 (diagram=60) --> iconfig=56 (f77 conv) and iconfigC=55 (c conv)
     56, // channelId=61 (diagram=61) --> iconfig=57 (f77 conv) and iconfigC=56 (c conv)
     57, // channelId=62 (diagram=62) --> iconfig=58 (f77 conv) and iconfigC=57 (c conv)
     58, // channelId=63 (diagram=63) --> iconfig=59 (f77 conv) and iconfigC=58 (c conv)
     59, // channelId=64 (diagram=64) --> iconfig=60 (f77 conv) and iconfigC=59 (c conv)
     60, // channelId=65 (diagram=65) --> iconfig=61 (f77 conv) and iconfigC=60 (c conv)
     61, // channelId=66 (diagram=66) --> iconfig=62 (f77 conv) and iconfigC=61 (c conv)
     62, // channelId=67 (diagram=67) --> iconfig=63 (f77 conv) and iconfigC=62 (c conv)
     63, // channelId=68 (diagram=68) --> iconfig=64 (f77 conv) and iconfigC=63 (c conv)
     64, // channelId=69 (diagram=69) --> iconfig=65 (f77 conv) and iconfigC=64 (c conv)
     65, // channelId=70 (diagram=70) --> iconfig=66 (f77 conv) and iconfigC=65 (c conv)
     66, // channelId=71 (diagram=71) --> iconfig=67 (f77 conv) and iconfigC=66 (c conv)
     67, // channelId=72 (diagram=72) --> iconfig=68 (f77 conv) and iconfigC=67 (c conv)
     68, // channelId=73 (diagram=73) --> iconfig=69 (f77 conv) and iconfigC=68 (c conv)
     -1, // channelId=74 (diagram=74): Not consider as a channel of integration (presence of 4 point interaction?)
     69, // channelId=75 (diagram=75) --> iconfig=70 (f77 conv) and iconfigC=69 (c conv)
     70, // channelId=76 (diagram=76) --> iconfig=71 (f77 conv) and iconfigC=70 (c conv)
     71, // channelId=77 (diagram=77) --> iconfig=72 (f77 conv) and iconfigC=71 (c conv)
     72, // channelId=78 (diagram=78) --> iconfig=73 (f77 conv) and iconfigC=72 (c conv)
     73, // channelId=79 (diagram=79) --> iconfig=74 (f77 conv) and iconfigC=73 (c conv)
     74, // channelId=80 (diagram=80) --> iconfig=75 (f77 conv) and iconfigC=74 (c conv)
     75, // channelId=81 (diagram=81) --> iconfig=76 (f77 conv) and iconfigC=75 (c conv)
     76, // channelId=82 (diagram=82) --> iconfig=77 (f77 conv) and iconfigC=76 (c conv)
     77, // channelId=83 (diagram=83) --> iconfig=78 (f77 conv) and iconfigC=77 (c conv)
     78, // channelId=84 (diagram=84) --> iconfig=79 (f77 conv) and iconfigC=78 (c conv)
     79, // channelId=85 (diagram=85) --> iconfig=80 (f77 conv) and iconfigC=79 (c conv)
     80, // channelId=86 (diagram=86) --> iconfig=81 (f77 conv) and iconfigC=80 (c conv)
     81, // channelId=87 (diagram=87) --> iconfig=82 (f77 conv) and iconfigC=81 (c conv)
     82, // channelId=88 (diagram=88) --> iconfig=83 (f77 conv) and iconfigC=82 (c conv)
     83, // channelId=89 (diagram=89) --> iconfig=84 (f77 conv) and iconfigC=83 (c conv)
     84, // channelId=90 (diagram=90) --> iconfig=85 (f77 conv) and iconfigC=84 (c conv)
     85, // channelId=91 (diagram=91) --> iconfig=86 (f77 conv) and iconfigC=85 (c conv)
     86, // channelId=92 (diagram=92) --> iconfig=87 (f77 conv) and iconfigC=86 (c conv)
     -1, // channelId=93 (diagram=93): Not consider as a channel of integration (presence of 4 point interaction?)
     87, // channelId=94 (diagram=94) --> iconfig=88 (f77 conv) and iconfigC=87 (c conv)
     88, // channelId=95 (diagram=95) --> iconfig=89 (f77 conv) and iconfigC=88 (c conv)
     89, // channelId=96 (diagram=96) --> iconfig=90 (f77 conv) and iconfigC=89 (c conv)
     90, // channelId=97 (diagram=97) --> iconfig=91 (f77 conv) and iconfigC=90 (c conv)
     91, // channelId=98 (diagram=98) --> iconfig=92 (f77 conv) and iconfigC=91 (c conv)
     92, // channelId=99 (diagram=99) --> iconfig=93 (f77 conv) and iconfigC=92 (c conv)
     -1, // channelId=100 (diagram=100): Not consider as a channel of integration (presence of 4 point interaction?)
     93, // channelId=101 (diagram=101) --> iconfig=94 (f77 conv) and iconfigC=93 (c conv)
     94, // channelId=102 (diagram=102) --> iconfig=95 (f77 conv) and iconfigC=94 (c conv)
     95, // channelId=103 (diagram=103) --> iconfig=96 (f77 conv) and iconfigC=95 (c conv)
     96, // channelId=104 (diagram=104) --> iconfig=97 (f77 conv) and iconfigC=96 (c conv)
     97, // channelId=105 (diagram=105) --> iconfig=98 (f77 conv) and iconfigC=97 (c conv)
     98, // channelId=106 (diagram=106) --> iconfig=99 (f77 conv) and iconfigC=98 (c conv)
     -1, // channelId=107 (diagram=107): Not consider as a channel of integration (presence of 4 point interaction?)
     99, // channelId=108 (diagram=108) --> iconfig=100 (f77 conv) and iconfigC=99 (c conv)
     100, // channelId=109 (diagram=109) --> iconfig=101 (f77 conv) and iconfigC=100 (c conv)
     101, // channelId=110 (diagram=110) --> iconfig=102 (f77 conv) and iconfigC=101 (c conv)
     102, // channelId=111 (diagram=111) --> iconfig=103 (f77 conv) and iconfigC=102 (c conv)
     103, // channelId=112 (diagram=112) --> iconfig=104 (f77 conv) and iconfigC=103 (c conv)
     104, // channelId=113 (diagram=113) --> iconfig=105 (f77 conv) and iconfigC=104 (c conv)
  };

  // Map iconfigC (in C indexing, i.e. iconfig-1) to the set of allowed colors
  // This array has N_config <= N_diagrams elements    
    __device__ constexpr bool icolamp[105][24] = {
    {true, true, false, false, false, false, true, true, false, false, false, false, true, false, true, false, true, true, true, false, true, false, true, true}, // iconfigC=0, diag=2
    {true, false, false, false, false, false, true, false, false, false, false, false, true, false, true, false, false, false, true, false, true, false, true, true}, // iconfigC=1, diag=3
    {false, true, false, false, false, false, false, true, false, false, false, false, true, false, true, false, true, true, true, false, true, false, false, false}, // iconfigC=2, diag=4
    {true, true, false, false, false, false, true, true, false, false, false, false, false, false, false, false, true, true, false, false, false, false, true, true}, // iconfigC=3, diag=5
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, false, false, false, false, false}, // iconfigC=4, diag=6
    {false, false, false, false, false, false, false, false, false, false, false, false, true, false, true, false, true, true, false, false, false, false, false, false}, // iconfigC=5, diag=7
    {false, false, false, false, false, false, false, false, false, false, false, false, true, false, true, false, false, false, false, false, false, false, false, false}, // iconfigC=6, diag=8
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true}, // iconfigC=7, diag=9
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, true, false, true, true}, // iconfigC=8, diag=10
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, true, false, false, false}, // iconfigC=9, diag=11
    {false, true, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=10, diag=12
    {false, true, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, true, false, true, false, false, false}, // iconfigC=11, diag=13
    {true, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=12, diag=14
    {true, false, false, false, false, false, true, false, false, false, false, false, true, false, true, false, false, false, false, false, false, false, false, false}, // iconfigC=13, diag=15
    {true, true, false, false, false, false, true, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=14, diag=16
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, false, false, false, true, true}, // iconfigC=15, diag=17
    {false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=16, diag=18
    {false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=17, diag=19
    {false, false, false, true, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=18, diag=20
    {true, false, true, false, true, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=19, diag=21
    {false, false, false, false, true, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=20, diag=22
    {true, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=21, diag=23
    {false, true, true, true, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=22, diag=24
    {false, false, true, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=23, diag=25
    {false, true, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=24, diag=26
    {false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=25, diag=27
    {false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=26, diag=28
    {true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=27, diag=29
    {false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=28, diag=30
    {true, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=29, diag=31
    {true, true, false, true, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=30, diag=33
    {true, true, true, true, true, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=31, diag=34
    {false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=32, diag=35
    {false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=33, diag=36
    {false, false, false, false, false, false, false, false, false, true, false, true, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=34, diag=37
    {false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, true, false, false, false, false, false, true, false, true}, // iconfigC=35, diag=38
    {false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, true, false, false, false, false, false, false, false, false}, // iconfigC=36, diag=39
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, true}, // iconfigC=37, diag=40
    {false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, true, false, true, false, false, false, true, false, false}, // iconfigC=38, diag=41
    {false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, true, false, false}, // iconfigC=39, diag=42
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, true, false, false, false, false, false, false}, // iconfigC=40, diag=43
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false}, // iconfigC=41, diag=44
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false}, // iconfigC=42, diag=45
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true}, // iconfigC=43, diag=46
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false}, // iconfigC=44, diag=47
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, true}, // iconfigC=45, diag=49
    {false, false, false, false, false, false, false, false, false, true, false, true, false, false, false, false, false, true, false, false, false, false, false, true}, // iconfigC=46, diag=50
    {false, false, false, false, false, false, false, false, false, true, false, true, false, false, false, true, false, true, false, false, false, true, false, true}, // iconfigC=47, diag=51
    {false, false, false, false, false, false, false, false, false, false, true, true, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=48, diag=52
    {false, false, false, false, false, false, true, false, true, false, true, true, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=49, diag=53
    {false, false, false, false, false, false, true, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=50, diag=54
    {false, false, false, true, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false}, // iconfigC=51, diag=55
    {false, false, false, true, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, true, false, false, true, false}, // iconfigC=52, diag=56
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, true, false}, // iconfigC=53, diag=57
    {false, false, true, true, false, false, false, false, false, false, false, false, true, true, false, false, false, false, false, false, false, false, false, false}, // iconfigC=54, diag=59
    {false, false, false, false, false, false, false, false, false, false, true, true, false, false, false, false, false, false, false, false, true, true, false, false}, // iconfigC=55, diag=60
    {false, false, true, true, false, false, false, false, false, false, true, true, true, true, false, false, false, false, false, false, true, true, false, false}, // iconfigC=56, diag=61
    {false, false, true, true, false, false, true, false, true, false, true, true, true, true, false, false, false, false, false, true, true, true, true, false}, // iconfigC=57, diag=62
    {false, false, true, false, false, false, true, false, true, false, false, false, true, false, false, false, false, false, false, true, true, true, true, false}, // iconfigC=58, diag=63
    {false, false, false, true, false, false, true, false, true, false, true, true, false, true, false, false, false, false, false, true, false, false, true, false}, // iconfigC=59, diag=64
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, true, true, false}, // iconfigC=60, diag=65
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, false}, // iconfigC=61, diag=66
    {false, false, true, false, false, false, true, false, true, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=62, diag=67
    {false, false, true, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=63, diag=68
    {false, false, false, false, false, false, false, false, true, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=64, diag=69
    {false, false, false, false, false, false, false, true, true, true, true, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=65, diag=70
    {false, false, false, false, false, false, false, true, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=66, diag=71
    {false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false}, // iconfigC=67, diag=72
    {false, false, false, false, false, true, false, false, false, false, false, false, false, true, false, false, true, false, false, true, false, false, false, false}, // iconfigC=68, diag=73
    {false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, true, false, false, false, false, false, false, false}, // iconfigC=69, diag=75
    {false, false, false, false, true, true, false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, false, false, false}, // iconfigC=70, diag=76
    {false, false, false, false, false, false, false, false, true, true, false, false, false, false, true, true, false, false, false, false, false, false, false, false}, // iconfigC=71, diag=77
    {false, false, false, false, true, true, false, false, true, true, false, false, false, false, true, true, false, false, true, true, false, false, false, false}, // iconfigC=72, diag=78
    {false, false, false, false, true, true, false, true, true, true, true, false, false, true, true, true, true, false, true, true, false, false, false, false}, // iconfigC=73, diag=79
    {false, false, false, false, true, false, false, true, false, false, true, false, false, true, true, true, true, false, true, false, false, false, false, false}, // iconfigC=74, diag=80
    {false, false, false, false, false, true, false, true, true, true, true, false, false, true, false, false, true, false, false, true, false, false, false, false}, // iconfigC=75, diag=81
    {false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, true, true, false, false, false, false, false, false, false}, // iconfigC=76, diag=82
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, false, false, false, false, false, false, false}, // iconfigC=77, diag=83
    {false, false, false, false, true, false, false, true, false, false, true, false, false, false, false, false, false, false, true, false, false, false, false, false}, // iconfigC=78, diag=84
    {false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false}, // iconfigC=79, diag=85
    {false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=80, diag=86
    {false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=81, diag=87
    {false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=82, diag=88
    {false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=83, diag=89
    {false, false, false, false, false, false, true, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=84, diag=90
    {false, false, false, false, false, false, true, true, false, true, false, true, false, false, false, false, false, false, false, false, false, false, false, false}, // iconfigC=85, diag=91
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false}, // iconfigC=86, diag=92
    {false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false}, // iconfigC=87, diag=94
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false}, // iconfigC=88, diag=95
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false}, // iconfigC=89, diag=96
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, true, false}, // iconfigC=90, diag=97
    {false, false, false, true, false, true, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, true, false}, // iconfigC=91, diag=98
    {true, false, true, false, true, true, false, false, true, true, false, false, false, false, true, true, false, false, true, true, false, true, false, true}, // iconfigC=92, diag=99
    {true, false, true, false, false, false, false, false, true, false, false, false, false, false, true, false, false, false, true, true, false, true, false, true}, // iconfigC=93, diag=101
    {true, false, true, false, true, true, false, false, false, true, false, false, false, false, false, true, false, false, false, false, false, true, false, true}, // iconfigC=94, diag=102
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, true, false, true}, // iconfigC=95, diag=103
    {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, false, false, false}, // iconfigC=96, diag=104
    {true, false, true, false, false, false, false, false, true, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false}, // iconfigC=97, diag=105
    {false, false, false, false, false, false, false, false, true, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false}, // iconfigC=98, diag=106
    {false, true, true, true, true, false, false, false, false, false, true, true, true, true, false, true, false, true, false, false, true, true, false, false}, // iconfigC=99, diag=108
    {false, true, false, false, true, false, false, false, false, false, true, false, true, true, false, true, false, true, false, false, true, false, false, false}, // iconfigC=100, diag=109
    {false, true, true, true, true, false, false, false, false, false, false, true, false, false, false, true, false, true, false, false, false, true, false, false}, // iconfigC=101, diag=110
    {false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, true, false, true, false, false, false, false, false, false}, // iconfigC=102, diag=111
    {false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, false, false, false, false, false, false, false, false, false}, // iconfigC=103, diag=112
    {false, true, false, false, true, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, true, false, false, false}, // iconfigC=104, diag=113
  };

}
#endif // COLORAMPS_H
