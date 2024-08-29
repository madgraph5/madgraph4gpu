// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jan 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2023) for the MG5aMC CUDACPP plugin.
//
// Copyright (C) 2021-2023 Argonne National Laboratory.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.

#ifndef MGONGPUFPTYPES_H
#define MGONGPUFPTYPES_H 1

#include "mgOnGpuConfig.h"

#include <CL/sycl.hpp>
#include <algorithm>
#include <cmath>

//==========================================================================

//------------------------------
// Floating point types - SYCL
//------------------------------

/*
SYCL_EXTERNAL inline
fptype fpmax( const fptype& a, const fptype& b )
{
  return max( a, b );
}

SYCL_EXTERNAL inline
fptype fpmin( const fptype& a, const fptype& b )
{
  return min( a, b );
}
*/

SYCL_EXTERNAL inline
const fptype& fpmax( const fptype& a, const fptype& b )
{
  return ( ( b < a ) ? a : b );
}

SYCL_EXTERNAL inline
const fptype& fpmin( const fptype& a, const fptype& b )
{
  return ( ( a < b ) ? a : b );
}

SYCL_EXTERNAL inline
fptype fpsqrt( const fptype& f )
{
#if defined MGONGPU_FPTYPE_FLOAT
  //FIXME wait for sqrtf support
  return sycl::sqrt( f );
#else
  //FIXME wait for std::sqrt support
  return sycl::sqrt( f );
#endif
}

//==========================================================================

#endif // MGONGPUFPTYPES_H
