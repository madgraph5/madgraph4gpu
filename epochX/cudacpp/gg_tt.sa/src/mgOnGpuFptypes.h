// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jan 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2023) for the MG5aMC CUDACPP plugin.

#ifndef MGONGPUFPTYPES_H
#define MGONGPUFPTYPES_H 1

#include "mgOnGpuConfig.h"

#include <algorithm>
#include <cmath>

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //==========================================================================

#ifdef __CUDACC__ // cuda

  //------------------------------
  // Floating point types - Cuda
  //------------------------------

  /*
  inline __host__ __device__ fptype
  fpmax( const fptype& a, const fptype& b )
  {
    return max( a, b );
  }

  inline __host__ __device__ fptype
  fpmin( const fptype& a, const fptype& b )
  {
    return min( a, b );
  }
  */

  inline __host__ __device__ const fptype&
  fpmax( const fptype& a, const fptype& b )
  {
    return ( ( b < a ) ? a : b );
  }

  inline __host__ __device__ const fptype&
  fpmin( const fptype& a, const fptype& b )
  {
    return ( ( a < b ) ? a : b );
  }

  inline __host__ __device__ fptype
  fpsqrt( const fptype& f )
  {
#if defined MGONGPU_FPTYPE_FLOAT
    // See https://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__SINGLE.html
    return sqrtf( f );
#else
    // See https://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__DOUBLE.html
    return sqrt( f );
#endif
  }

#endif // #ifdef __CUDACC__

  //==========================================================================

#ifndef __CUDACC__

  //------------------------------
  // Floating point types - C++
  //------------------------------

  inline const fptype&
  fpmax( const fptype& a, const fptype& b )
  {
    return std::max( a, b );
  }

  inline const fptype&
  fpmin( const fptype& a, const fptype& b )
  {
    return std::min( a, b );
  }

  inline fptype
  fpsqrt( const fptype& f )
  {
    return std::sqrt( f );
  }

#endif // #ifndef __CUDACC__

  //==========================================================================

} // end namespace mg5amcGpu/mg5amcCpu

#endif // MGONGPUFPTYPES_H
