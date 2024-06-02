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
  //   => this number (with F indexing) is passed around as an API argument between cudacpp functions
  // - Channel number in C indexing: "channelIdC" = channelID - 1
  //   => this number (with C indexing) is used as the index of the channelIdC_to_iconfig array below
  // - Config number ("iconfig" in C, ICONFIG in F) in [1, N_config]: all values are allowed (N_config <= N_diagrams distinct values)
  // - Config number in C indexing: "iconfigC" = iconfig - 1
  
  // Map channelIdC (in C indexing, i.e. channelId-1) to iconfig (in F indexing)
  // This array has N_diagrams elements, but only N_config <= N_diagrams valid (non-zero) values
  __device__ constexpr int channelIdC_to_iconfig[6] = {
    0, // channelId=1 (diagram=1) i.e. channelIdC=0 --> iconfig=None
    1, // channelId=2 (diagram=2) i.e. channelIdC=1 --> iconfig=1
    2, // channelId=3 (diagram=3) i.e. channelIdC=2 --> iconfig=2
    2, // channelId=4 (diagram=4) i.e. channelIdC=3 --> iconfig=3
    3, // channelId=5 (diagram=5) i.e. channelIdC=4 --> iconfig=4
    4  // channelId=6 (diagram=6) i.e. channelIdC=5 --> iconfig=5
  };

  // Map iconfigC (in C indexing, i.e. iconfig-1) to the set of allowed colors
  // This array has N_config <= N_diagrams elements
  __device__ constexpr bool icolamp[5][2] = {
    {  true,  true }, // iconfig=1 i.e. iconfigC=0
    {  true,  true }, // iconfig=2 i.e. iconfigC=1
    {  true, false }, // iconfig=3 i.e. iconfigC=2
    {  true, false }, // iconfig=4 i.e. iconfigC=3
    { false,  true }  // iconfig=5 i.e. iconfigC=4
  };

} /* clang-format off */
#endif // COLORAMPS_H
