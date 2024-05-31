// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#ifndef COLORAMPS_H
#define COLORAMPS_H 1

#include <map>

namespace mgOnGpu
{

  __device__ constexpr int diag_to_channel[7] = {
    -1, // 0 --> None
    -1, // 1 --> None
    +0, // 2 --> 0
    +1, // 3 --> 1
    +2, // 4 --> 2
    +3, // 5 --> 3
    +4  // 6 --> 4
  };

  __device__ constexpr bool icolamp[5][2] = {
    { true, true },
    { true, true },
    { true, false },
    { true, false },
    { false, true } };

}
#endif // COLORAMPS_H
