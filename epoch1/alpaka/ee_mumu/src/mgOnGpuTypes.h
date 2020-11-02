#ifndef MGONGPUCOMPLEX_H
#define MGONGPUCOMPLEX_H 1

#include "mgOnGpuConfig.h"

#if defined MGONGPU_CXTYPE_THRUST
#include <thrust/complex.h>
#elif defined MGONGPU_CXTYPE_CUCOMPLEX
#include <complex>
#include <cuComplex.h>
#elif defined MGONGPU_CXTYPE_EXTRAS
#include <complex>
#include "extras.h"
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
  typedef thrust::complex<fptype> cxtypeparam; // two doubles: RI

#elif defined MGONGPU_CXTYPE_EXTRAS
#ifdef MGONGPU_FPTYPE_DOUBLE
  // use class directly from extras (see 'using' below)
#else
  #error not handled comples type for MGONGPU_CXTYPE_EXTRAS
#endif

#elif defined MGONGPU_FPTYPE_DOUBLE
  typedef cuDoubleComplex cxtype;
  typedef cuDoubleComplex cxtypeparam;
#elif defined MGONGPU_FPTYPE_FLOAT
  typedef cuFloatComplex cxtype;
  typedef cuFloatComplex cxtypeparam;
#endif

  // Vector types: <type>_v is a <type>[256]
  //typedef double double_v[ntpbMAX];
  //typedef cxtype cxtype_v[ntpbMAX]; // RIRIRIRI: eventually move to RRRRIIII?

}

// Expose typedefs and operators outside the namespace
using mgOnGpu::fptype;

#if MGONGPU_CXTYPE_EXTRAS
using cxtypeparam = extras::cxtypeparam<fptype>;
using extras::cxtyperef;
using extras::cxtypeconstref;
using extras::cxtypeptr;
using extras::cxtypeop;
#else
using mgOnGpu::cxtype;
// fixme: wip
#define cxtyperef(a,b) \
    auto &a = b;
#define cxtypeptr(a,b) \
    auto &a = b;
#endif

// --- Functions and operators for complex types

//------------------------------
// CUPLA - using thrust::complex
//------------------------------

#if defined MGONGPU_CXTYPE_THRUST // cupla + thrust

//+++++++++++++++++++++++++
// thrust::complex<double>
// thrust::complex<float>
//+++++++++++++++++++++++++

inline __host__ __device__
cxtype cxmake( const fptype& r, const fptype& i )
{
  return cxtype( r, i ); // thrust::complex<fptype> constructor
}

inline __host__ __device__
fptype cxreal( const cxtype& c )
{
  return c.real(); // thrust::complex<fptype>::real()
}

inline __host__ __device__
fptype cximag( const cxtype& c )
{
  return c.imag(); // thrust::complex<fptype>::imag()
}

inline __host__ __device__
const cxtype& cxmake( const cxtype& c )
{
  return c;
}

//------------------------------
// CUPLA - using extras::compelx
//------------------------------

#elif defined MGONGPU_CXTYPE_EXTRAS // cupla + "extras" complex impl

inline ALPAKA_FN_ACC
cxtypeparam cxmake( const fptype& r, const fptype& i )
{
  cxtypeparam p = { .r = r, .i = i };
  return p;
}

template< typename T_Acc, typename FPType >
inline ALPAKA_FN_ACC
FPType cxreal( const cxtyperef<T_Acc,FPType>& c )
{
  return c.real();
}

template< typename T_Acc, typename FPType >
inline ALPAKA_FN_ACC
FPType cxreal( const cxtypeop<T_Acc,FPType>& c )
{
  return c.real();
}

template< typename T_Acc, typename FPType >
inline ALPAKA_FN_ACC
FPType cximag( const cxtyperef<T_Acc,FPType>& c )
{
  return c.imag();
}

template< typename T_Acc, typename FPType >
inline ALPAKA_FN_ACC
FPType cximag( const cxtypeop<T_Acc,FPType>& c )
{
  return c.imag();
}

template< typename T_Acc, typename FPType >
inline ALPAKA_FN_ACC
const cxtypeparam& cxmake( const cxtyperef<T_Acc,FPType>& c )
{
  return cxtypeparam(c.real(), c.imag());
}

inline ALPAKA_FN_HOST
cxtypeparam cxmake( const std::complex<fptype>& c )
{
  return cxmake( (fptype)c.real(), (fptype)c.imag() );
}


//------------------------------
// CUPLA - using cuComplex: separate functions for fptype double or float
//------------------------------

#elif defined MGONGPU_CXTYPE_CUCOMPLEX // cupla + cucomplex

//+++++++++++++++++++++++++
// cuDoubleComplex
//+++++++++++++++++++++++++

#if defined MGONGPU_FPTYPE_DOUBLE  // cupla + cucomplex + double

inline __host__ __device__
cxtype cxmake( const fptype& r, const fptype& i )
{
  return make_cuDoubleComplex( r, i );
}

inline __host__ __device__
fptype cxreal( const cxtype& c )
{
  return cuCreal(c); // returns by value
}

inline __host__ __device__
fptype cximag( const cxtype& c )
{
  return cuCimag(c); // returns by value
}

inline __host__ __device__
cxtype operator+( const cxtype& a, const cxtype& b )
{
  return cuCadd( a, b );
}

inline __host__ __device__
cxtype operator-( const cxtype& a, const cxtype& b )
{
  return cuCsub( a, b );
}

inline __host__ __device__
cxtype operator*( const cxtype& a, const cxtype& b )
{
  return cuCmul( a, b );
}

inline __host__ __device__
cxtype operator/( const cxtype& a, const cxtype& b )
{
  return cuCdiv( a, b );
}

//+++++++++++++++++++++++++
// cuFloatComplex
//+++++++++++++++++++++++++

#elif defined MGONGPU_FPTYPE_FLOAT  // cupla + cucomplex + float

inline __host__ __device__
cxtype cxmake( const fptype& r, const fptype& i )
{
  return make_cuFloatComplex( r, i );
}

inline __host__ __device__
fptype cxreal( const cxtype& c )
{
  return cuCrealf(c); // returns by value
}

inline __host__ __device__
fptype cximag( const cxtype& c )
{
  return cuCimagf(c); // returns by value
}

inline __host__ __device__
cxtype operator+( const cxtype& a, const cxtype& b )
{
  return cuCaddf( a, b );
}

inline __host__ __device__
cxtype operator-( const cxtype& a, const cxtype& b )
{
  return cuCsubf( a, b );
}

inline __host__ __device__
cxtype operator*( const cxtype& a, const cxtype& b )
{
  return cuCmulf( a, b );
}

inline __host__ __device__
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

#endif  // END cupla + cucomplex + double/float

inline __host__ __device__
cxtype operator+( const cxtype a )
{
  return a;
}

inline __host__ __device__
cxtype operator-( const cxtype& a )
{
  return cxmake( -cxreal(a), -cximag(a) );
}

inline __host__ __device__
cxtype operator+( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) + b;
}

inline __host__ __device__
cxtype operator-( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) - b;
}

inline __host__ __device__
cxtype operator*( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) * b;
}

inline __host__ __device__
cxtype operator/( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) / b;
}

inline __host__ __device__
cxtype operator+( const cxtype& a, const fptype& b )
{
  return a + cxmake( b, 0 );
}

inline __host__ __device__
cxtype operator-( const cxtype& a, const fptype& b )
{
  return a - cxmake( b, 0 );
}

inline __host__ __device__
cxtype operator*( const cxtype& a, const fptype& b )
{
  return a * cxmake( b, 0 );
}

inline __host__ __device__
cxtype operator/( const cxtype& a, const fptype& b )
{
  return a / cxmake( b, 0 );
}

inline __host__ __device__
cxtype conj( const cxtype& c )
{
  return cxmake( cxreal( c ), -cximag( c ) );
}

inline __host__ // NOT __device__
cxtype cxmake( const std::complex<fptype>& c ) // std::complex to cucomplex (float-to-float or double-to-double)
{
  return cxmake( c.real(), c.imag() );
}

#endif  // END cupla + cucomplex

#endif // MGONGPUCOMPLEX_H
