#ifndef MGONGPUFPTYPES_H
#define MGONGPUFPTYPES_H 1

#include "mgOnGpuConfig.h"

//#include <CL/sycl.hpp>
#include <algorithm>
#include <cmath>

//==========================================================================

//------------------------------
// Floating point types - SYCL
//------------------------------


/* KOKKOS_INLINE_FUNCTION */
/* fptype fpmax( const fptype& a, const fptype& b ) */
/* { */
/* 	return std::max( a, b ); */
/* } */

/* KOKKOS_INLINE_FUNCTION */
/* fptype fpmin( const fptype& a, const fptype& b ) */
/* { */
/* 	return std::min( a, b ); */
/* } */


KOKKOS_INLINE_FUNCTION
const fptype& fpmax( const fptype& a, const fptype& b )
{
  return ( ( b < a ) ? a : b );
}

KOKKOS_INLINE_FUNCTION
const fptype& fpmin( const fptype& a, const fptype& b )
{
  return ( ( a < b ) ? a : b );
}

KOKKOS_INLINE_FUNCTION
fptype fpsqrt( const fptype& f )
{
#if defined MGONGPU_FPTYPE_FLOAT
  //FIXME wait for sqrtf support
	return std::sqrt( f );
#else
  //FIXME wait for std::sqrt support
  return std::sqrt( f );
#endif
}

//==========================================================================

#endif // MGONGPUFPTYPES_H
