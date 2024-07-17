// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: O. Mattelaer, A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#ifndef COLORAMPS_H
#define COLORAMPS_H 1

#include "CPPProcess.h"

// Note: strictly speaking the check '#ifdef MGONGPU_SUPPORTS_MULTICHANNEL' is not needed here,
// because coloramps.h is not included otherwise, but adding it does not harm and makes the code clearer

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL /* clang-format off */

namespace mgOnGpu
{
  // Summary of numbering and indexing conventions for the relevant concepts (see issue #826 and PR #852)
  // - Diagram number (no variable) in [1, N_diagrams]: all values are allowed (N_diagrams distinct values)
  //   => this number is displayed for information before each block of code in CPPProcess.cc
  // - Channel number ("channelId" in C, CHANNEL_ID in F) in [1, N_diagrams]: not all values are allowed (N_config <= N_diagrams distinct values)
  //   => this number (with F indexing as in ps/pdf output) is passed around as an API argument between cudacpp functions
  //   Note: the old API passes around a single CHANNEL_ID (and uses CHANNEL_ID=0 to indicate no-multichannel mode, but this is not used in coloramps.h),
  //   while the new API passes around an array of CHANNEL_ID's (and uses a NULL array pointer to indicate no-multichannel mode)
  // - Channel number in C indexing: "channelID - 1"
  //   => this number (with C indexing) is used as the index of the channel2iconfig array below
  // - Config number ("iconfig" in C, ICONFIG in F) in [1, N_config]: all values are allowed (N_config <= N_diagrams distinct values)
  // - Config number in C indexing: "iconfig - 1"
  //   => this number (with C indexing) is used as the index of the icolamp array below

  // The number of channels in the channel2iconfig array below
  // !!!FIXME!!! This should be N_diagrams, but it is now different from CPPProcess::ndiagrams (see #919 and #910)
  constexpr unsigned int nchannels = 6;
#ifdef MGONGPUCPP_GPUIMPL
  //static_assert( nchannels == mg5amcGpu::CPPProcess::ndiagrams, "mismatch between nchannels and ndiagrams?!" ); // sanity check #910 (fails #919)
#else
  //static_assert( nchannels == mg5amcCpu::CPPProcess::ndiagrams, "mismatch between nchannels and ndiagrams?!" ); // sanity check #910 (fails #919)
#endif
  
  // Map channel to iconfig (e.g. "iconfig = channel2iconfig[channelId - 1]": input index uses C indexing, output index uses F indexing)
  // Note: iconfig=-1 indicates channels/diagrams with no associated iconfig for single-diagram enhancement in the MadEvent sampling algorithm (presence of 4-point interaction?)
  // This array has N_diagrams elements, but only N_config <= N_diagrams valid values (iconfig>0)
  // (NB: this array is created on the host in C++ code and on the device in GPU code, but a host copy is also needed in runTest #917)
  __device__ constexpr int channel2iconfig[6] = { // note: a trailing comma in the initializer list is allowed
    -1, // CHANNEL_ID=1  i.e. DIAGRAM=1  --> ICONFIG=-1 (diagram with no associated iconfig for single-diagram enhancement)
     1, // CHANNEL_ID=2  i.e. DIAGRAM=2  --> ICONFIG=1
     2, // CHANNEL_ID=3  i.e. DIAGRAM=3  --> ICONFIG=2
     3, // CHANNEL_ID=4  i.e. DIAGRAM=4  --> ICONFIG=3
     4, // CHANNEL_ID=5  i.e. DIAGRAM=5  --> ICONFIG=4
     5, // CHANNEL_ID=6  i.e. DIAGRAM=6  --> ICONFIG=5
  };

  // Host copy of the channel2iconfig array (this is needed in runTest #917)
#ifndef MGONGPUCPP_GPUIMPL
  constexpr const int* hostChannel2iconfig = channel2iconfig;
#else
  constexpr int hostChannel2iconfig[6] = { // note: a trailing comma in the initializer list is allowed
    -1, // CHANNEL_ID=1  i.e. DIAGRAM=1  --> ICONFIG=-1 (diagram with no associated iconfig for single-diagram enhancement)
     1, // CHANNEL_ID=2  i.e. DIAGRAM=2  --> ICONFIG=1
     2, // CHANNEL_ID=3  i.e. DIAGRAM=3  --> ICONFIG=2
     3, // CHANNEL_ID=4  i.e. DIAGRAM=4  --> ICONFIG=3
     4, // CHANNEL_ID=5  i.e. DIAGRAM=5  --> ICONFIG=4
     5, // CHANNEL_ID=6  i.e. DIAGRAM=6  --> ICONFIG=5
  };
#endif

  // The number N_config of channels/diagrams with an associated iconfig for single-diagram enhancement in the MadEvent sampling algorithm (#917)
  constexpr unsigned int nconfigSDE = 5;

  // Map iconfig to the mask of allowed colors (e.g. "colormask = icolamp[iconfig - 1]": input index uses C indexing)
  // This array has N_config <= N_diagrams elements
  // (NB: this array is created on the host in C++ code and on the device in GPU code)
  __device__ constexpr bool icolamp[5][2] = { // note: a trailing comma in the initializer list is allowed
    {  true,  true }, // ICONFIG=1  <-- CHANNEL_ID=2
    {  true, false }, // ICONFIG=2  <-- CHANNEL_ID=3
    {  true, false }, // ICONFIG=3  <-- CHANNEL_ID=4
    { false,  true }, // ICONFIG=4  <-- CHANNEL_ID=5
    { false,  true }, // ICONFIG=5  <-- CHANNEL_ID=6
  };

}
#endif /* clang-format on */

#endif // COLORAMPS_H
