#ifndef MGONGPUCOMPLEX_H
#define MGONGPUCOMPLEX_H 1

#include <cuda_to_cupla.hpp>
#include "mgOnGpuConfig.h"

#if defined MGONGPU_CXTYPE_THRUST
#include <thrust/complex.h>
#elif defined MGONGPU_CXTYPE_CUCOMPLEX
#include <complex>
#include <cuComplex.h>
#elif defined MGONGPU_CXTYPE_ALSIMPLE
#include <complex>
#include <alsimple.h>
#endif

namespace mgOnGpu
{

  // --- Type definitions

  // Floating point type: fptype
#if defined MGONGPU_FPTYPE_DOUBLE
  typedef double fptype; // double precision (8 bytes, fp64)
#elif defined MGONGPU_FPTYPE_FLOAT
  typedef float fptype; // single precision (4 bytes, fp32)
#endif

  // Complex type: cxtype
#if defined MGONGPU_CXTYPE_THRUST
  typedef thrust::complex<fptype> cxtype; // two doubles: RI
#elif defined MGONGPU_CXTYPE_CUCOMPLEX
#if defined MGONGPU_FPTYPE_DOUBLE
  typedef cuDoubleComplex cxtype;
#elif defined MGONGPU_FPTYPE_FLOAT
  typedef cuFloatComplex cxtype;
#endif
#elif defined MGONGPU_CXTYPE_ALSIMPLE
  typedef alsimple::complex<fptype> cxtype;
#else
#error unrecognised MGONGPU_CXTYPE selection
#endif

}

// Expose typedefs and operators outside the namespace
using mgOnGpu::fptype;
using mgOnGpu::cxtype;

// --- Functions and operators for complex types

//------------------------------
// using thrust::complex
//------------------------------

#if defined MGONGPU_CXTYPE_THRUST // thrust

//+++++++++++++++++++++++++
// thrust::complex<double>
// thrust::complex<float>
//+++++++++++++++++++++++++

inline ALPAKA_FN_ACC
cxtype cxmake( const fptype& r, const fptype& i )
{
  return cxtype( r, i ); // thrust::complex<fptype> constructor
}

inline ALPAKA_FN_ACC
fptype cxreal( const cxtype& c )
{
  return c.real(); // thrust::complex<fptype>::real()
}

inline ALPAKA_FN_ACC
fptype cximag( const cxtype& c )
{
  return c.imag(); // thrust::complex<fptype>::imag()
}

inline ALPAKA_FN_ACC
const cxtype& cxmake( const cxtype& c )
{
  return c;
}

//------------------------------
// using cuComplex
//------------------------------

#elif defined MGONGPU_CXTYPE_CUCOMPLEX // cucomplex

//+++++++++++++++++++++++++
// cuDoubleComplex
//+++++++++++++++++++++++++

#if defined MGONGPU_FPTYPE_DOUBLE  // cucomplex + double

inline ALPAKA_FN_ACC
cxtype cxmake( const fptype& r, const fptype& i )
{
  return make_cuDoubleComplex( r, i );
}

inline ALPAKA_FN_ACC
fptype cxreal( const cxtype& c )
{
  return cuCreal(c); // returns by value
}

inline ALPAKA_FN_ACC
fptype cximag( const cxtype& c )
{
  return cuCimag(c); // returns by value
}

inline ALPAKA_FN_ACC
cxtype operator+( const cxtype& a, const cxtype& b )
{
  return cuCadd( a, b );
}

inline ALPAKA_FN_ACC
cxtype operator-( const cxtype& a, const cxtype& b )
{
  return cuCsub( a, b );
}

inline ALPAKA_FN_ACC
cxtype operator*( const cxtype& a, const cxtype& b )
{
  return cuCmul( a, b );
}

inline ALPAKA_FN_ACC
cxtype operator/( const cxtype& a, const cxtype& b )
{
  return cuCdiv( a, b );
}

//+++++++++++++++++++++++++
// cuFloatComplex
//+++++++++++++++++++++++++

#elif defined MGONGPU_FPTYPE_FLOAT  // cucomplex + float

inline ALPAKA_FN_ACC
cxtype cxmake( const fptype& r, const fptype& i )
{
  return make_cuFloatComplex( r, i );
}

inline ALPAKA_FN_ACC
fptype cxreal( const cxtype& c )
{
  return cuCrealf(c); // returns by value
}

inline ALPAKA_FN_ACC
fptype cximag( const cxtype& c )
{
  return cuCimagf(c); // returns by value
}

inline ALPAKA_FN_ACC
cxtype operator+( const cxtype& a, const cxtype& b )
{
  return cuCaddf( a, b );
}

inline ALPAKA_FN_ACC
cxtype operator-( const cxtype& a, const cxtype& b )
{
  return cuCsubf( a, b );
}

inline ALPAKA_FN_ACC
cxtype operator*( const cxtype& a, const cxtype& b )
{
  return cuCmulf( a, b );
}

inline ALPAKA_FN_ACC
cxtype operator/( const cxtype& a, const cxtype& b )
{
  return cuCdivf( a, b );
}

inline ALPAKA_FN_HOST // NOT __device__
cxtype cxmake( const std::complex<double>& c ) // std::complex to cucomplex (cast double-to-float)
{
  return cxmake( (fptype)c.real(), (fptype)c.imag() );
}

#endif  // END cucomplex + double/float

//------------------------------
// using alsimple::complex
//------------------------------

#elif defined MGONGPU_CXTYPE_ALSIMPLE

//+++++++++++++++++++++++++
// alsimple::complex<double>
// alsimple::complex<float>
//+++++++++++++++++++++++++

inline ALPAKA_FN_ACC
cxtype cxmake( const fptype& r, const fptype& i )
{
  return cxtype( r, i ); // alsimple::complex<fptype> constructor
}

inline ALPAKA_FN_ACC
fptype cxreal( const cxtype& c )
{
  return c.real(); // alsime::complex<fptype>::real()
}

inline ALPAKA_FN_ACC
fptype cximag( const cxtype& c )
{
  return c.imag(); // alsime::complex<fptype>::imag()
}

inline ALPAKA_FN_ACC
const cxtype& cxmake( const cxtype& c )
{
  return c;
}

#endif  // END cuda + thrust/cucomplex

//+++++++++++++++++++++++++
// required non-member operator overloads for ALSIMPLE or CUCOMPLEX
//+++++++++++++++++++++++++

#if defined MGONGPU_CXTYPE_ALSIMPLE || defined MGONGPU_CXTYPE_CUCOMPLEX

inline ALPAKA_FN_ACC
cxtype operator+( const cxtype a )
{
  return a;
}

inline ALPAKA_FN_ACC
cxtype operator-( const cxtype& a )
{
  return cxmake( -cxreal(a), -cximag(a) );
}

inline ALPAKA_FN_ACC
cxtype operator+( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) + b;
}

inline ALPAKA_FN_ACC
cxtype operator-( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) - b;
}

inline ALPAKA_FN_ACC
cxtype operator*( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) * b;
}

inline ALPAKA_FN_ACC
cxtype operator/( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) / b;
}

inline ALPAKA_FN_ACC
cxtype operator+( const cxtype& a, const fptype& b )
{
  return a + cxmake( b, 0 );
}

inline ALPAKA_FN_ACC
cxtype operator-( const cxtype& a, const fptype& b )
{
  return a - cxmake( b, 0 );
}

inline ALPAKA_FN_ACC
cxtype operator*( const cxtype& a, const fptype& b )
{
  return a * cxmake( b, 0 );
}

inline ALPAKA_FN_ACC
cxtype operator/( const cxtype& a, const fptype& b )
{
  return a / cxmake( b, 0 );
}

inline ALPAKA_FN_ACC
cxtype conj( const cxtype& c )
{
  return cxmake( cxreal( c ), -cximag( c ) );
}

inline ALPAKA_FN_HOST // NOT __device__
cxtype cxmake( const std::complex<fptype>& c ) // std::complex to cucomplex (float-to-float or double-to-double)
{
  return cxmake( c.real(), c.imag() );
}

#endif

#endif // MGONGPUCOMPLEX_H
