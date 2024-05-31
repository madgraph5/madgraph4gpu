// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#ifndef COLORAMPS_H
#define COLORAMPS_H 1

#include <map>

namespace mgOnGpu
{

  __device__ constexpr int diag_to_channel[%(nb_diagmax)s] = {
%(diag_to_channel)s
  };

  __device__ constexpr bool icolamp[%(nb_channel)s][%(nb_color)s] = {
%(is_LC)s };

}
#endif // COLORAMPS_H
