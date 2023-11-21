// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2023) for the MG5aMC CUDACPP plugin.
//==========================================================================
// Copyright (C) 2021-2023 Argonne National Laboratory.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: N. Nichols (2022-2023) for the MG5aMC SYCL plugin.
//==========================================================================

#ifndef COLORAMPS_H
#define COLORAMPS_H 1

namespace mgOnGpu
{

  template <typename T>
  constexpr T icolamp[] = { // FIXME: assume process.nprocesses == 1 for the moment
     true, true,
     true, false,
     false, true
  };

}
#endif // COLORAMPS_H
