// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: O. Mattelaer, A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#ifndef COLORAMPS_H
#define COLORAMPS_H 1

namespace mgOnGpu
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
  __device__ constexpr int channelIdC_to_iconfig[6] = { // note: a trailing comma in the initializer list is allowed
    0, // channelIdC=0 i.e. channelId=1 (diagram=1) --> iconfig=None
    1, // channelIdC=1 i.e. channelId=2 (diagram=2) --> iconfig=1
    2, // channelIdC=2 i.e. channelId=3 (diagram=3) --> iconfig=2
    3, // channelIdC=3 i.e. channelId=4 (diagram=4) --> iconfig=3
    4, // channelIdC=4 i.e. channelId=5 (diagram=5) --> iconfig=4
    5, // channelIdC=5 i.e. channelId=6 (diagram=6) --> iconfig=5
  };

  // Map iconfigC (in C indexing, i.e. iconfig-1) to the set of allowed colors
  // This array has N_config <= N_diagrams elements
  __device__ constexpr bool icolamp[5][2] = { // note: a trailing comma in the initializer list is allowed
    {  true,  true }, // iconfigC=0 i.e. iconfig=1 <-- diagram=2
    {  true,  true }, // iconfigC=1 i.e. iconfig=2 <-- diagram=3
    {  true, false }, // iconfigC=2 i.e. iconfig=3 <-- diagram=4
    {  true, false }, // iconfigC=3 i.e. iconfig=4 <-- diagram=5
    { false,  true }, // iconfigC=4 i.e. iconfig=5 <-- diagram=6
  };

}
#endif // COLORAMPS_H
