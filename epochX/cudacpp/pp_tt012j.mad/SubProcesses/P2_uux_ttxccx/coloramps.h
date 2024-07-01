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
  __device__ constexpr int channelIdC_to_iconfig[7] = { // note: a trailing comma in the initializer list is allowed
    1, // CHANNEL_ID=1 i.e. DIAGRAM=1 --> ICONFIG=1
    2, // CHANNEL_ID=2 i.e. DIAGRAM=2 --> ICONFIG=2
    3, // CHANNEL_ID=3 i.e. DIAGRAM=3 --> ICONFIG=3
    4, // CHANNEL_ID=4 i.e. DIAGRAM=4 --> ICONFIG=4
    5, // CHANNEL_ID=5 i.e. DIAGRAM=5 --> ICONFIG=5
    6, // CHANNEL_ID=6 i.e. DIAGRAM=6 --> ICONFIG=6
    7, // CHANNEL_ID=7 i.e. DIAGRAM=7 --> ICONFIG=7
  };

  // Map iconfigC (in C indexing, i.e. iconfig-1) to the set of allowed colors
  // This array has N_config <= N_diagrams elements
  __device__ constexpr bool icolamp[7][6] = { // note: a trailing comma in the initializer list is allowed
    { false, false, false, false,  true, false }, // ICONFIG=1 <-- CHANNEL_ID=1
    { false, false, false,  true, false, false }, // ICONFIG=2 <-- CHANNEL_ID=2
    { false, false, false,  true,  true, false }, // ICONFIG=3 <-- CHANNEL_ID=3
    { false, false, false, false,  true, false }, // ICONFIG=4 <-- CHANNEL_ID=4
    { false, false, false,  true, false, false }, // ICONFIG=5 <-- CHANNEL_ID=5
    { false, false, false,  true, false, false }, // ICONFIG=6 <-- CHANNEL_ID=6
    { false, false, false, false,  true, false }, // ICONFIG=7 <-- CHANNEL_ID=7
  }; /* clang-format on */

}
#endif // COLORAMPS_H
