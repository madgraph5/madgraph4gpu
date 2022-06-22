#ifndef MGONGPUVECTORS_H
#define MGONGPUVECTORS_H 1

#include "mgOnGpuFptypes.h"
#include "mgOnGpuCxtypes.h"

#include <iostream>

//------------------------------
// Vector types - C++
//------------------------------

namespace mgOnGpu
{
	const int neppV = 1;
}

//--------------------------------------------------------------------------

// Expose typedefs outside the namespace
using mgOnGpu::neppV;

//--------------------------------------------------------------------------

//==========================================================================


//------------------------------
// Vector types - SYCL
//------------------------------

// Printout to std::cout for user defined types
inline void print( const fptype& f ){ printf( "%f\n", f ); }
inline void print( const cxtype& c ){ printf( "[%f, %f]\n", cxreal(c), cximag(c) ); }


KOKKOS_INLINE_FUNCTION
fptype fpternary( const bool& mask, const fptype& a, const fptype& b )
{
  return ( mask ? a : b );
}

KOKKOS_INLINE_FUNCTION
cxtype cxternary( const bool& mask, const cxtype& a, const cxtype& b )
{
  return ( mask ? a : b );
}


//==========================================================================

// Scalar-or-vector types: scalar in SYCL
typedef bool bool_sv;
typedef fptype fptype_sv;
typedef cxtype cxtype_sv;


// Scalar-or-vector zeros: scalar in CUDA, vector or scalar in C++
KOKKOS_INLINE_FUNCTION cxtype cxzero_sv(){ return cxtype{ fptype{0}, fptype{0} }; }

//--------------------------------------------------------------------------

#endif // MGONGPUVECTORS_H
