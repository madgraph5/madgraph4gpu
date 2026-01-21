// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jan 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: J. Teig, A. Valassi, Z. Wettersten (2022-2025) for the MG5aMC CUDACPP plugin.

#ifndef MGONGPUFPTYPES_H
#define MGONGPUFPTYPES_H 1

#include "mgOnGpuConfig.h"

#include <algorithm>
#include <cmath>
#ifdef MGONGPU_FPTYPE_QUAD
#include <quadmath.h>
#include <ostream> // for operator<< overloading
#endif

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef MGONGPUCPP_GPUIMPL // cuda
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //==========================================================================

#ifdef MGONGPUCPP_GPUIMPL // cuda

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

#endif // #ifdef MGONGPUCPP_GPUIMPL

  //==========================================================================

#ifndef MGONGPUCPP_GPUIMPL

#ifdef MGONGPU_FPTYPE_QUAD

  //------------------------------
  // Quad precision types - C++
  //------------------------------

  // quad max returns by value, not by reference
  inline fptype
  fpmax( const fptype& a, const fptype& b )
  {
    return fmaxq( a, b );
  }

  // quad min returns by value, not by reference
  inline fptype
  fpmin( const fptype& a, const fptype& b )
  {
    return fminq( a, b );
  }

  inline fptype
  fpsqrt( const fptype& f )
  {
    return sqrtq( f );
  }

  inline fptype
  fpabs( const fptype& f )
  {
    return fabsq( f );
  }

  inline fptype
  fppow( const fptype& base, const int& exp )
  {
    // Special case for positive integer exponents
    if ( exp >= 0 )
    {
      fptype result = static_cast<__float128>( 1.0 );
      for ( int i = 0; i < exp; ++i )
      {
        result *= base;
      }
      return result;
    }
    // General case
    return powq( base, static_cast<__float128>( exp ) );
  }

  inline fptype
  fppow( const fptype& base, const fptype& exp )
  {
    return powq( base, exp );
  }

  // Overload operator<< for quad precision (need to convert to string first)
  inline std::ostream&
  operator<<( std::ostream& os, const fptype& f )
  {
    char buffer[128];
    // Convert __float128 to string with 36 digits of precision
    int n = quadmath_snprintf( buffer, sizeof( buffer ), "%.36Qg", f );
    if ( n > 0 && n < static_cast<int>( sizeof( buffer ) ) )
    {
      os << buffer;
    }
    else
    {
      os << "ConversionError";
    }
    return os;
  }

  inline bool
  fp_is_infinite( const fptype& fp )
  {
    return isinfq( fp );
  }

  inline bool
  fp_is_nan( const fptype& fp )
  {
    return isnanq( fp );
  }

#define UNUSED(x) (void)(x)
  inline bool
  fp_is_normal( const fptype& fp )
  {
    UNUSED(fp);
    return true; // no isnormalq in quadmath.h
  }

#else

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

  inline fptype
  fpabs( const fptype& f )
  {
    return std::abs( f );
  }

  inline fptype
  fppow( const fptype& base, const int& exp )
  {
    return std::pow( base, static_cast<fptype>( exp ) );
  }

  inline bool
  fp_is_infinite( const fptype& fp )
  {
    return std::isinf( fp );
  }

  inline bool
  fp_is_nan( const fptype& fp )
  {
    //#pragma clang diagnostic push
    //#pragma clang diagnostic ignored "-Wtautological-compare" // for icpx2021/clang13 (https://stackoverflow.com/a/15864661)
    return std::isnan( fp ); // always false for clang in fast math mode (tautological compare)?
    //#pragma clang diagnostic pop
  }

  inline bool
  fp_is_normal( const fptype& fp )
  {
    return std::isnormal( fp );
  }

#endif // #ifdef MGONGPU_FPTYPE_QUAD
#endif // #ifndef MGONGPUCPP_GPUIMPL

  //==========================================================================

} // end namespace mg5amcGpu/mg5amcCpu

#endif // MGONGPUFPTYPES_H
