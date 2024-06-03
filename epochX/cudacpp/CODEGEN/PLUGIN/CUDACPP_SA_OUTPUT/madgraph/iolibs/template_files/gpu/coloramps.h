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
  // - Channel number ("channelId" in C, CHANNEL_ID in F) in [1, N_diagrams]: not all values are allowed (N_config <= N_diagrams distinct values)
  //   0 is allowed to fallback to no multi-channel mode.
  //   => this number (with indexing like ps/pdf output) is passed around as an API argument between cudacpp functions
  // - Channel number in C indexing: "channelIdC" = channelIdC_to_iconfig[channelId]
  //   => this number (with C indexing) is used as the index of the channelIdC_to_iconfig array below
  //   This correspond to iconfig on the fortran side (with iconfig = channelIdC + 1)
  //NOTE: All those ordering are event by event specific (with the intent to have those fix within a vector size/wrap   
  
  // Map channelId to channelIdC
  // This array has N_diagrams+1 elements, but only N_config <= N_diagrams valid (non-zero) values
  // The 0 entry is a fall back to still write events even if no multi-channel is setup (wrong color selected in that mode) 
    __device__ constexpr int channelId_to_channelIdC[%(nb_diag_plus_one)i] = {
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
