// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#ifndef COLORAMPS_H
#define COLORAMPS_H 1

#include <map>

namespace mgOnGpu
{

  __device__ std::map<int, int> diag_to_channel = {
    { 2, 0 },
    { 3, 1 },
    { 4, 2 },
    { 5, 3 },
    { 6, 4 }, // note: a trailing comma in the initializer list is allowed
  };

  __device__ constexpr bool icolamp[5][2] = {
    { true, true },
    { true, true },
    { true, false },
    { true, false },
    { false, true } };

}
#endif // COLORAMPS_H
