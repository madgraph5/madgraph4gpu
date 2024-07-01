// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: O. Mattelaer, A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#ifndef COLORAMPS_H
#define COLORAMPS_H 1

namespace mgOnGpu /* clang-format off */
{
  // Summary of numbering and indexing conventions for the relevant concepts (see issue #826 and PR #852)
  // - Diagram number (no variable) in [1, N_diagrams]: all values are allowed (N_diagrams distinct values)
  //   => this number is displayed for information before each block of code in CPPProcess.cc
  // - Channel number ("channelId" in C, CHANNEL_ID in F) in [1, N_diagrams]: not all values are allowed (N_config <= N_diagrams distinct values)
  //   => this number (with F indexing as in ps/pdf output) is passed around as an API argument between cudacpp functions
  //   Note: the old API passes around a single CHANNEL_ID (and uses CHANNEL_ID=0 to indicate no-multichannel mode, but this is not used in coloramps.h),
  //   while the new API passes around an array of CHANNEL_ID's (and uses a NULL array pointer to indicate no-multichannel mode)
  // - Channel number in C indexing: "channelIdC" = channelID - 1
  //   => this number (with C indexing) is used as the index of the channelIdC_to_iconfig array below
  // - Config number ("iconfig" in C, ICONFIG in F) in [1, N_config]: all values are allowed (N_config <= N_diagrams distinct values)
  // - Config number in C indexing: "iconfigC" = iconfig - 1
  //   => this number (with C indexing) is used as the index of the icolamp array below

  // Map channelIdC (in C indexing, i.e. channelId-1) to iconfig (in F indexing)
  // Note: iconfig=0 indicates invalid values, i.e. channels/diagrams with no single-diagram enhancement in the MadEvent sampling algorithm (presence of 4-point interaction?)
  // This array has N_diagrams elements, but only N_config <= N_diagrams valid (non-zero) values
  __device__ constexpr int channelIdC_to_iconfig[72] = { // note: a trailing comma in the initializer list is allowed
     1, // CHANNEL_ID=1  i.e. DIAGRAM=1  --> ICONFIG=1
     2, // CHANNEL_ID=2  i.e. DIAGRAM=2  --> ICONFIG=2
     3, // CHANNEL_ID=3  i.e. DIAGRAM=3  --> ICONFIG=3
     4, // CHANNEL_ID=4  i.e. DIAGRAM=4  --> ICONFIG=4
     5, // CHANNEL_ID=5  i.e. DIAGRAM=5  --> ICONFIG=5
     6, // CHANNEL_ID=6  i.e. DIAGRAM=6  --> ICONFIG=6
     7, // CHANNEL_ID=7  i.e. DIAGRAM=7  --> ICONFIG=7
     8, // CHANNEL_ID=8  i.e. DIAGRAM=8  --> ICONFIG=8
     9, // CHANNEL_ID=9  i.e. DIAGRAM=9  --> ICONFIG=9
    10, // CHANNEL_ID=10 i.e. DIAGRAM=10 --> ICONFIG=10
    11, // CHANNEL_ID=11 i.e. DIAGRAM=11 --> ICONFIG=11
    12, // CHANNEL_ID=12 i.e. DIAGRAM=12 --> ICONFIG=12
    13, // CHANNEL_ID=13 i.e. DIAGRAM=13 --> ICONFIG=13
    14, // CHANNEL_ID=14 i.e. DIAGRAM=14 --> ICONFIG=14
    15, // CHANNEL_ID=15 i.e. DIAGRAM=15 --> ICONFIG=15
    16, // CHANNEL_ID=16 i.e. DIAGRAM=16 --> ICONFIG=16
    17, // CHANNEL_ID=17 i.e. DIAGRAM=17 --> ICONFIG=17
    18, // CHANNEL_ID=18 i.e. DIAGRAM=18 --> ICONFIG=18
    19, // CHANNEL_ID=19 i.e. DIAGRAM=19 --> ICONFIG=19
    20, // CHANNEL_ID=20 i.e. DIAGRAM=20 --> ICONFIG=20
    21, // CHANNEL_ID=21 i.e. DIAGRAM=21 --> ICONFIG=21
    22, // CHANNEL_ID=22 i.e. DIAGRAM=22 --> ICONFIG=22
    23, // CHANNEL_ID=23 i.e. DIAGRAM=23 --> ICONFIG=23
    24, // CHANNEL_ID=24 i.e. DIAGRAM=24 --> ICONFIG=24
    25, // CHANNEL_ID=25 i.e. DIAGRAM=25 --> ICONFIG=25
    26, // CHANNEL_ID=26 i.e. DIAGRAM=26 --> ICONFIG=26
    27, // CHANNEL_ID=27 i.e. DIAGRAM=27 --> ICONFIG=27
    28, // CHANNEL_ID=28 i.e. DIAGRAM=28 --> ICONFIG=28
    29, // CHANNEL_ID=29 i.e. DIAGRAM=29 --> ICONFIG=29
    30, // CHANNEL_ID=30 i.e. DIAGRAM=30 --> ICONFIG=30
    31, // CHANNEL_ID=31 i.e. DIAGRAM=31 --> ICONFIG=31
    32, // CHANNEL_ID=32 i.e. DIAGRAM=32 --> ICONFIG=32
    33, // CHANNEL_ID=33 i.e. DIAGRAM=33 --> ICONFIG=33
    34, // CHANNEL_ID=34 i.e. DIAGRAM=34 --> ICONFIG=34
    35, // CHANNEL_ID=35 i.e. DIAGRAM=35 --> ICONFIG=35
    36, // CHANNEL_ID=36 i.e. DIAGRAM=36 --> ICONFIG=36
    37, // CHANNEL_ID=37 i.e. DIAGRAM=37 --> ICONFIG=37
    38, // CHANNEL_ID=38 i.e. DIAGRAM=38 --> ICONFIG=38
    39, // CHANNEL_ID=39 i.e. DIAGRAM=39 --> ICONFIG=39
    40, // CHANNEL_ID=40 i.e. DIAGRAM=40 --> ICONFIG=40
    41, // CHANNEL_ID=41 i.e. DIAGRAM=41 --> ICONFIG=41
    42, // CHANNEL_ID=42 i.e. DIAGRAM=42 --> ICONFIG=42
    43, // CHANNEL_ID=43 i.e. DIAGRAM=43 --> ICONFIG=43
    44, // CHANNEL_ID=44 i.e. DIAGRAM=44 --> ICONFIG=44
    45, // CHANNEL_ID=45 i.e. DIAGRAM=45 --> ICONFIG=45
    46, // CHANNEL_ID=46 i.e. DIAGRAM=46 --> ICONFIG=46
    47, // CHANNEL_ID=47 i.e. DIAGRAM=47 --> ICONFIG=47
    48, // CHANNEL_ID=48 i.e. DIAGRAM=48 --> ICONFIG=48
    49, // CHANNEL_ID=49 i.e. DIAGRAM=49 --> ICONFIG=49
    50, // CHANNEL_ID=50 i.e. DIAGRAM=50 --> ICONFIG=50
    51, // CHANNEL_ID=51 i.e. DIAGRAM=51 --> ICONFIG=51
    52, // CHANNEL_ID=52 i.e. DIAGRAM=52 --> ICONFIG=52
    53, // CHANNEL_ID=53 i.e. DIAGRAM=53 --> ICONFIG=53
    54, // CHANNEL_ID=54 i.e. DIAGRAM=54 --> ICONFIG=54
    55, // CHANNEL_ID=55 i.e. DIAGRAM=55 --> ICONFIG=55
    56, // CHANNEL_ID=56 i.e. DIAGRAM=56 --> ICONFIG=56
    57, // CHANNEL_ID=57 i.e. DIAGRAM=57 --> ICONFIG=57
    58, // CHANNEL_ID=58 i.e. DIAGRAM=58 --> ICONFIG=58
    59, // CHANNEL_ID=59 i.e. DIAGRAM=59 --> ICONFIG=59
    60, // CHANNEL_ID=60 i.e. DIAGRAM=60 --> ICONFIG=60
    61, // CHANNEL_ID=61 i.e. DIAGRAM=61 --> ICONFIG=61
    62, // CHANNEL_ID=62 i.e. DIAGRAM=62 --> ICONFIG=62
    63, // CHANNEL_ID=63 i.e. DIAGRAM=63 --> ICONFIG=63
    64, // CHANNEL_ID=64 i.e. DIAGRAM=64 --> ICONFIG=64
    65, // CHANNEL_ID=65 i.e. DIAGRAM=65 --> ICONFIG=65
    66, // CHANNEL_ID=66 i.e. DIAGRAM=66 --> ICONFIG=66
     0, // CHANNEL_ID=67 i.e. DIAGRAM=67 --> ICONFIG=None
    67, // CHANNEL_ID=68 i.e. DIAGRAM=68 --> ICONFIG=67
    68, // CHANNEL_ID=69 i.e. DIAGRAM=69 --> ICONFIG=68
     0, // CHANNEL_ID=70 i.e. DIAGRAM=70 --> ICONFIG=None
    69, // CHANNEL_ID=71 i.e. DIAGRAM=71 --> ICONFIG=69
    70, // CHANNEL_ID=72 i.e. DIAGRAM=72 --> ICONFIG=70
  };

  // Map iconfigC (in C indexing, i.e. iconfig-1) to the set of allowed colors
  // This array has N_config <= N_diagrams elements
  __device__ constexpr bool icolamp[70][12] = { // note: a trailing comma in the initializer list is allowed
    { false, false,  true, false, false, false, false, false, false, false,  true, false }, // ICONFIG=1  <-- CHANNEL_ID=1
    { false,  true, false, false, false, false, false, false, false,  true, false, false }, // ICONFIG=2  <-- CHANNEL_ID=2
    { false,  true,  true, false, false, false, false, false, false,  true,  true, false }, // ICONFIG=3  <-- CHANNEL_ID=3
    {  true, false, false, false, false, false, false, false,  true, false, false, false }, // ICONFIG=4  <-- CHANNEL_ID=4
    { false, false, false,  true, false, false, false, false, false, false, false,  true }, // ICONFIG=5  <-- CHANNEL_ID=5
    {  true, false, false,  true, false, false, false, false,  true, false, false,  true }, // ICONFIG=6  <-- CHANNEL_ID=6
    {  true, false, false, false, false, false, false, false,  true, false, false, false }, // ICONFIG=7  <-- CHANNEL_ID=7
    { false, false, false,  true, false, false, false, false, false, false, false,  true }, // ICONFIG=8  <-- CHANNEL_ID=8
    { false,  true, false, false, false, false, false, false, false,  true, false, false }, // ICONFIG=9  <-- CHANNEL_ID=9
    { false, false,  true, false, false, false, false, false, false, false,  true, false }, // ICONFIG=10 <-- CHANNEL_ID=10
    {  true, false, false, false, false, false, false, false, false, false, false, false }, // ICONFIG=11 <-- CHANNEL_ID=11
    { false, false, false, false, false,  true, false, false, false, false, false, false }, // ICONFIG=12 <-- CHANNEL_ID=12
    { false, false, false, false, false,  true, false, false, false, false, false, false }, // ICONFIG=13 <-- CHANNEL_ID=13
    { false, false, false, false,  true, false, false, false, false, false, false, false }, // ICONFIG=14 <-- CHANNEL_ID=14
    { false,  true, false, false, false, false, false, false, false, false, false, false }, // ICONFIG=15 <-- CHANNEL_ID=15
    { false, false, false, false,  true, false, false, false, false, false, false, false }, // ICONFIG=16 <-- CHANNEL_ID=16
    {  true, false, false, false, false, false, false, false, false, false, false, false }, // ICONFIG=17 <-- CHANNEL_ID=17
    {  true, false, false, false,  true, false, false, false, false, false, false, false }, // ICONFIG=18 <-- CHANNEL_ID=18
    { false,  true, false, false, false, false, false, false, false, false, false, false }, // ICONFIG=19 <-- CHANNEL_ID=19
    { false,  true, false, false, false,  true, false, false, false, false, false, false }, // ICONFIG=20 <-- CHANNEL_ID=20
    { false, false, false, false, false, false, false, false,  true, false, false, false }, // ICONFIG=21 <-- CHANNEL_ID=21
    { false, false, false, false, false, false,  true, false, false, false, false, false }, // ICONFIG=22 <-- CHANNEL_ID=22
    { false, false, false, false, false, false, false, false, false, false,  true, false }, // ICONFIG=23 <-- CHANNEL_ID=23
    { false, false, false, false,  true, false, false, false, false, false, false, false }, // ICONFIG=24 <-- CHANNEL_ID=24
    { false, false, false, false, false, false,  true, false, false, false, false, false }, // ICONFIG=25 <-- CHANNEL_ID=25
    { false, false, false, false,  true, false, false, false, false, false, false, false }, // ICONFIG=26 <-- CHANNEL_ID=26
    { false, false, false, false, false, false, false, false,  true, false, false, false }, // ICONFIG=27 <-- CHANNEL_ID=27
    { false, false, false, false,  true, false, false, false,  true, false, false, false }, // ICONFIG=28 <-- CHANNEL_ID=28
    { false, false, false, false, false, false, false, false, false, false,  true, false }, // ICONFIG=29 <-- CHANNEL_ID=29
    { false, false, false, false, false, false,  true, false, false, false,  true, false }, // ICONFIG=30 <-- CHANNEL_ID=30
    { false, false, false, false, false, false, false,  true, false, false, false, false }, // ICONFIG=31 <-- CHANNEL_ID=31
    { false, false, false, false, false, false,  true, false, false, false, false, false }, // ICONFIG=32 <-- CHANNEL_ID=32
    { false, false,  true, false, false, false, false, false, false, false, false, false }, // ICONFIG=33 <-- CHANNEL_ID=33
    { false, false, false, false, false, false, false,  true, false, false, false, false }, // ICONFIG=34 <-- CHANNEL_ID=34
    { false, false, false,  true, false, false, false, false, false, false, false, false }, // ICONFIG=35 <-- CHANNEL_ID=35
    { false, false, false, false, false, false,  true, false, false, false, false, false }, // ICONFIG=36 <-- CHANNEL_ID=36
    { false, false,  true, false, false, false, false, false, false, false, false, false }, // ICONFIG=37 <-- CHANNEL_ID=37
    { false, false,  true, false, false, false,  true, false, false, false, false, false }, // ICONFIG=38 <-- CHANNEL_ID=38
    { false, false, false,  true, false, false, false, false, false, false, false, false }, // ICONFIG=39 <-- CHANNEL_ID=39
    { false, false, false,  true, false, false, false,  true, false, false, false, false }, // ICONFIG=40 <-- CHANNEL_ID=40
    { false, false, false, false, false, false, false, false, false,  true, false, false }, // ICONFIG=41 <-- CHANNEL_ID=41
    { false, false, false, false, false, false, false,  true, false, false, false, false }, // ICONFIG=42 <-- CHANNEL_ID=42
    { false, false, false, false, false, false, false,  true, false, false, false, false }, // ICONFIG=43 <-- CHANNEL_ID=43
    { false, false, false, false, false,  true, false, false, false, false, false, false }, // ICONFIG=44 <-- CHANNEL_ID=44
    { false, false, false, false, false, false, false, false, false, false, false,  true }, // ICONFIG=45 <-- CHANNEL_ID=45
    { false, false, false, false, false,  true, false, false, false, false, false, false }, // ICONFIG=46 <-- CHANNEL_ID=46
    { false, false, false, false, false, false, false, false, false,  true, false, false }, // ICONFIG=47 <-- CHANNEL_ID=47
    { false, false, false, false, false,  true, false, false, false,  true, false, false }, // ICONFIG=48 <-- CHANNEL_ID=48
    { false, false, false, false, false, false, false, false, false, false, false,  true }, // ICONFIG=49 <-- CHANNEL_ID=49
    { false, false, false, false, false, false, false,  true, false, false, false,  true }, // ICONFIG=50 <-- CHANNEL_ID=50
    { false, false, false, false, false, false, false, false,  true, false, false, false }, // ICONFIG=51 <-- CHANNEL_ID=51
    { false, false, false, false, false, false, false,  true,  true, false, false, false }, // ICONFIG=52 <-- CHANNEL_ID=52
    { false, false, false, false, false, false, false, false, false,  true, false, false }, // ICONFIG=53 <-- CHANNEL_ID=53
    { false, false, false, false, false, false,  true, false, false,  true, false, false }, // ICONFIG=54 <-- CHANNEL_ID=54
    {  true, false, false, false, false, false, false, false, false, false, false, false }, // ICONFIG=55 <-- CHANNEL_ID=55
    {  true, false, false, false, false, false, false,  true, false, false, false, false }, // ICONFIG=56 <-- CHANNEL_ID=56
    { false, false,  true, false, false, false, false, false, false, false, false, false }, // ICONFIG=57 <-- CHANNEL_ID=57
    { false, false,  true, false, false,  true, false, false, false, false, false, false }, // ICONFIG=58 <-- CHANNEL_ID=58
    { false, false, false, false, false, false, false, false, false, false,  true, false }, // ICONFIG=59 <-- CHANNEL_ID=59
    { false, false, false, false, false,  true, false, false, false, false,  true, false }, // ICONFIG=60 <-- CHANNEL_ID=60
    { false, false, false, false, false, false, false, false, false, false, false,  true }, // ICONFIG=61 <-- CHANNEL_ID=61
    { false, false, false, false,  true, false, false, false, false, false, false,  true }, // ICONFIG=62 <-- CHANNEL_ID=62
    { false,  true, false, false, false, false, false, false, false, false, false, false }, // ICONFIG=63 <-- CHANNEL_ID=63
    { false,  true, false, false, false, false,  true, false, false, false, false, false }, // ICONFIG=64 <-- CHANNEL_ID=64
    { false, false, false,  true, false, false, false, false, false, false, false, false }, // ICONFIG=65 <-- CHANNEL_ID=65
    { false, false, false,  true,  true, false, false, false, false, false, false, false }, // ICONFIG=66 <-- CHANNEL_ID=66
    { false,  true, false, false, false,  true,  true, false, false, false,  true, false }, // ICONFIG=67 <-- CHANNEL_ID=68
    { false, false,  true, false, false,  true,  true, false, false,  true, false, false }, // ICONFIG=68 <-- CHANNEL_ID=69
    {  true, false, false, false,  true, false, false,  true, false, false, false,  true }, // ICONFIG=69 <-- CHANNEL_ID=71
    { false, false, false,  true,  true, false, false,  true,  true, false, false, false }, // ICONFIG=70 <-- CHANNEL_ID=72
  }; /* clang-format on */

}
#endif // COLORAMPS_H
