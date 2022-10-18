#ifndef MGONGPUTYPES_H
#define MGONGPUTYPES_H 1

#include "mgOnGpuConfig.h"

#ifdef KOKKOS_ENABLE_CUDA
#include "nvToolsExt.h"
#else
#include <cmath>
inline void nvtxRangePush(char* temp){return;}
inline void nvtxRangePop(void){return;};
#endif


#ifdef MGONGPU_CXTYPE_THRUST
  #include <thrust/complex.h>
  typedef thrust::complex<fptype> cxtype; // two doubles: RI
  #define COMPLEX_TYPE_NAME "THRUST::COMPLEX"
#else
  #include "Kokkos_Complex.hpp"
  typedef Kokkos::complex<fptype> cxtype; // two doubles: RI
  #define COMPLEX_TYPE_NAME "KOKKOS::COMPLEX"
#endif

// --- Functions and operators for floating point types


// use predefined cuda math library functions for double or float

fptype fpmax(const fptype&,const fptype&);
fptype fpmin(const fptype& ,const fptype&);
fptype fpsqrt(const fptype&);
fptype fpabs(const fptype&);
int iabs(const int&);


// --- Functions and operators for complex types

//------------------------------
// CUDA
//------------------------------

#ifdef KOKKOS_ENABLE_CUDA // cuda

//------------------------------
// CUDA - using thrust::complex
//------------------------------

#if defined MGONGPU_CXTYPE_THRUST // cuda + thrust

//+++++++++++++++++++++++++
// thrust::complex<double>
// thrust::complex<float>
//+++++++++++++++++++++++++

KOKKOS_INLINE_FUNCTION
cxtype cxmake( const fptype& r, const fptype& i )
{
  return cxtype( r, i ); // thrust::complex<fptype> constructor
}

KOKKOS_INLINE_FUNCTION
fptype cxreal( const cxtype& c )
{
  return c.real(); // thrust::complex<fptype>::real()
}

KOKKOS_INLINE_FUNCTION
fptype cximag( const cxtype& c )
{
  return c.imag(); // thrust::complex<fptype>::imag()
}

KOKKOS_INLINE_FUNCTION
cxtype cxconj( const cxtype& c )
{
  return conj( c ); // conj( thrust::complex<fptype> )
}

KOKKOS_INLINE_FUNCTION
const cxtype& cxmake( const cxtype& c )
{
  return c;
}

//------------------------------
// CUDA - using cuComplex
//------------------------------

#elif defined MGONGPU_CXTYPE_CUCOMPLEX // cuda + cucomplex

//+++++++++++++++++++++++++
// cuDoubleComplex
//+++++++++++++++++++++++++

#if defined MGONGPU_FPTYPE_DOUBLE  // cuda + cucomplex + double

KOKKOS_INLINE_FUNCTION
cxtype cxmake( const fptype& r, const fptype& i )
{
  return make_cuDoubleComplex( r, i );
}

KOKKOS_INLINE_FUNCTION
fptype cxreal( const cxtype& c )
{
  return cuCreal(c); // returns by value
}

KOKKOS_INLINE_FUNCTION
fptype cximag( const cxtype& c )
{
  return cuCimag(c); // returns by value
}

KOKKOS_INLINE_FUNCTION
cxtype operator+( const cxtype& a, const cxtype& b )
{
  return cuCadd( a, b );
}

KOKKOS_INLINE_FUNCTION
cxtype& operator+=( cxtype& a, const cxtype& b )
{
  a = cuCadd( a, b ); return a;
}

KOKKOS_INLINE_FUNCTION
cxtype operator-( const cxtype& a, const cxtype& b )
{
  return cuCsub( a, b );
}

KOKKOS_INLINE_FUNCTION
cxtype& operator-=( cxtype& a, const cxtype& b )
{
  a = cuCsub( a, b ); return a;
}

KOKKOS_INLINE_FUNCTION
cxtype operator*( const cxtype& a, const cxtype& b )
{
  return cuCmul( a, b );
}

KOKKOS_INLINE_FUNCTION
cxtype operator/( const cxtype& a, const cxtype& b )
{
  return cuCdiv( a, b );
}

//+++++++++++++++++++++++++
// cuFloatComplex
//+++++++++++++++++++++++++

#elif defined MGONGPU_FPTYPE_FLOAT  // cuda + cucomplex + float

KOKKOS_INLINE_FUNCTION
cxtype cxmake( const fptype& r, const fptype& i )
{
  return make_cuFloatComplex( r, i );
}

KOKKOS_INLINE_FUNCTION
fptype cxreal( const cxtype& c )
{
  return cuCrealf(c); // returns by value
}

KOKKOS_INLINE_FUNCTION
fptype cximag( const cxtype& c )
{
  return cuCimagf(c); // returns by value
}

KOKKOS_INLINE_FUNCTION
cxtype operator+( const cxtype& a, const cxtype& b )
{
  return cuCaddf( a, b );
}

KOKKOS_INLINE_FUNCTION
cxtype& operator+=( cxtype& a, const cxtype& b )
{
  a = cuCaddf( a, b ); return a;
}

KOKKOS_INLINE_FUNCTION
cxtype operator-( const cxtype& a, const cxtype& b )
{
  return cuCsubf( a, b );
}

KOKKOS_INLINE_FUNCTION
cxtype& operator-=( cxtype& a, const cxtype& b )
{
  a = cuCsubf( a, b ); return a;
}

KOKKOS_INLINE_FUNCTION
cxtype operator*( const cxtype& a, const cxtype& b )
{
  return cuCmulf( a, b );
}

KOKKOS_INLINE_FUNCTION
cxtype operator/( const cxtype& a, const cxtype& b )
{
  return cuCdivf( a, b );
}

inline __host__ // NOT __device__
cxtype cxmake( const std::complex<double>& c ) // std::complex to cucomplex (cast double-to-float)
{
  return cxmake( (fptype)c.real(), (fptype)c.imag() );
}

//+++++++++++++++++++++++++
// cuDoubleComplex
// cuFloatComplex
//+++++++++++++++++++++++++

#endif  // END cuda + cucomplex + double/float

KOKKOS_INLINE_FUNCTION
cxtype operator+( const cxtype a )
{
  return a;
}

KOKKOS_INLINE_FUNCTION
cxtype operator-( const cxtype& a )
{
  return cxmake( -cxreal(a), -cximag(a) );
}

KOKKOS_INLINE_FUNCTION
cxtype operator+( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) + b;
}

KOKKOS_INLINE_FUNCTION
cxtype operator-( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) - b;
}

KOKKOS_INLINE_FUNCTION
cxtype operator*( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) * b;
}

KOKKOS_INLINE_FUNCTION
cxtype operator/( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) / b;
}

KOKKOS_INLINE_FUNCTION
cxtype operator+( const cxtype& a, const fptype& b )
{
  return a + cxmake( b, 0 );
}

KOKKOS_INLINE_FUNCTION
cxtype operator-( const cxtype& a, const fptype& b )
{
  return a - cxmake( b, 0 );
}

KOKKOS_INLINE_FUNCTION
cxtype operator*( const cxtype& a, const fptype& b )
{
  return a * cxmake( b, 0 );
}

KOKKOS_INLINE_FUNCTION
cxtype operator/( const cxtype& a, const fptype& b )
{
  return a / cxmake( b, 0 );
}

KOKKOS_INLINE_FUNCTION
cxtype cxconj( const cxtype& c )
{
  return cxmake( cxreal( c ), -cximag( c ) );
}

KOKKOS_INLINE_FUNCTION
cxtype cxmake( const std::complex<fptype>& c ) // std::complex to cucomplex (float-to-float or double-to-double)
{
  return cxmake( c.real(), c.imag() );
}

#endif  // END cuda + thrust/cucomplex

//------------------------------
// C++ - using std::complex
//------------------------------

#else  // c++

//+++++++++++++++++++++++++
// std::complex<float>
// std::complex<double>
//+++++++++++++++++++++++++

KOKKOS_INLINE_FUNCTION
cxtype cxmake( const fptype& r, const fptype& i )
{
  return cxtype( r, i ); // std::complex<fptype> constructor
}

KOKKOS_INLINE_FUNCTION
fptype cxreal( const cxtype& c )
{
  return c.real(); // std::complex<fptype>::real()
}

KOKKOS_INLINE_FUNCTION
fptype cximag( const cxtype& c )
{
  return c.imag(); // std::complex<fptype>::imag()
}

KOKKOS_INLINE_FUNCTION
cxtype cxconj( const cxtype& c )
{
  return conj( c ); // conj( std::complex<fptype> )
}

KOKKOS_INLINE_FUNCTION
const cxtype& cxmake( const cxtype& c ) // std::complex to std::complex (float-to-float or double-to-double)
{
  return c;
}

#if defined MGONGPU_FPTYPE_FLOAT
KOKKOS_INLINE_FUNCTION
cxtype cxmake( const std::complex<double>& c ) // std::complex to std::complex (cast double-to-float)
{
  return cxmake( (fptype)c.real(), (fptype)c.imag() );
}
#endif

#endif  // END cuda/c++

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

#endif // MGONGPUTYPES_H
