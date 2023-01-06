#ifndef MGONGPUTYPES_H
#define MGONGPUTYPES_H 1

#include <CL/sycl.hpp>
#ifdef MGONGPU_COMPLEX_CXSMPL
    #include "mgOnGpuCxtypes.h"
#endif

#ifdef MGONGPU_COMPLEX_EXTRAS
    #include "extras.h"
#endif

#ifdef MGONGPU_COMPLEX_STD
    #include <complex>
#endif

#ifdef MGONGPU_COMPLEX_ONEAPI
    #define SYCL_EXT_ONEAPI_COMPLEX 1
    #include <sycl/ext/oneapi/experimental/sycl_complex.hpp>
#endif

#ifdef MGONGPU_COMPLEX_CUTHRUST
    #include <thrust/complex.h>
#endif

#ifdef MGONGPU_COMPLEX_CUCOMPLEX
    #include <cuComplex.h>
#endif

#include "mgOnGpuConfig.h"
#include "mgOnGpuCxtypes.h"

namespace mgOnGpu
{

  // --- Type definitions (complex type: cxtype)
  #ifdef MGONGPU_COMPLEX_CXSMPL
      typedef mgOnGpu::cxsmpl<fptype> cxtype;
  #endif

  #ifdef MGONGPU_COMPLEX_EXTRAS
      typedef extras::complex<fptype> cxtype;
  #endif

  #ifdef MGONGPU_COMPLEX_STD
      typedef std::complex<fptype> cxtype;
  #endif

  #ifdef MGONGPU_COMPLEX_ONEAPI
      typedef sycl::ext::oneapi::experimental::complex<fptype> cxtype;
  #endif

  #ifdef MGONGPU_COMPLEX_CUTHRUST
      typedef thrust::complex<fptype> cxtype;
  #endif
  
  #ifdef MGONGPU_COMPLEX_CUCOMPLEX
      #if defined MGONGPU_FPTYPE_DOUBLE
          typedef cuDoubleComplex cxtype;
      #elif defined MGONGPU_FPTYPE_FLOAT
          typedef cuFloatComplex cxtype;
      #endif
  #endif

  // The number of floating point types in a complex type (real, imaginary)
  constexpr int nx2 = 2;
  
  // SANITY CHECK: memory access may be based on casts of fptype[2] to cxtype (e.g. for wavefunctions)
  static_assert( sizeof(cxtype) == nx2 * sizeof(fptype), "sizeof(cxtype) is not 2*sizeof(fptype)" );
}

// Expose typedefs and operators outside the namespace
using mgOnGpu::cxtype;

//==========================================================================
// COMPLEX TYPES: (PLATFORM-SPECIFIC) FUNCTIONS AND OPERATORS
//==========================================================================
//------------------------------
// SYCL - extras::complex
//------------------------------
SYCL_EXTERNAL inline
cxtype cxmake( const fptype& r, const fptype& i )
{
#if defined MGONGPU_COMPLEX_CUCOMPLEX && MGONGPU_FPTYPE_DOUBLE
  return make_cuDoubleComplex( r, i );
#elif defined MGONGPU_COMPLEX_CUCOMPLEX && MGONGPU_FPTYPE_FLOAT
  return make_cuFloatComplex( r, i );
#else
  return cxtype( r, i ); // thrust::complex<fptype> constructor
#endif
}

SYCL_EXTERNAL inline
fptype cxreal( const cxtype& c )
{
#if defined MGONGPU_COMPLEX_CUCOMPLEX && MGONGPU_FPTYPE_DOUBLE
  return cuCreal( c ); // returns by value
#elif defined MGONGPU_COMPLEX_CUCOMPLEX && MGONGPU_FPTYPE_FLOAT
  return cuCrealf( c ); // returns by value
#else
  return c.real(); // thrust::complex<fptype>::real()
#endif
}

SYCL_EXTERNAL inline
fptype cximag( const cxtype& c )
{
#if defined MGONGPU_COMPLEX_CUCOMPLEX && MGONGPU_FPTYPE_DOUBLE
  return cuCimag( c ); // returns by value
#elif defined MGONGPU_COMPLEX_CUCOMPLEX && MGONGPU_FPTYPE_FLOAT
  return cuCimagf( c ); // returns by value
#else
  return c.imag(); // thrust::complex<fptype>::imag()
#endif
}

#if defined MGONGPU_COMPLEX_CUCOMPLEX
SYCL_EXTERNAL inline cxtype
operator+( const cxtype& a, const cxtype& b )
{
#if MGONGPU_FPTYPE_DOUBLE
  return cuCadd( a, b );
#elif defined MGONGPU_FPTYPE_FLOAT
  return cuCaddf( a, b );
#endif
}

SYCL_EXTERNAL inline cxtype&
operator+=( cxtype& a, const cxtype& b )
{
#if MGONGPU_FPTYPE_DOUBLE
  a = cuCadd( a, b );
#elif defined MGONGPU_FPTYPE_FLOAT
  a = cuCaddf( a, b );
#endif
  return a;
}

SYCL_EXTERNAL inline cxtype
operator-( const cxtype& a, const cxtype& b )
{
#if MGONGPU_FPTYPE_DOUBLE
  return cuCsub( a, b );
#elif defined MGONGPU_FPTYPE_FLOAT
  return cuCsubf( a, b );
#endif
}

SYCL_EXTERNAL inline cxtype&
operator-=( cxtype& a, const cxtype& b )
{
#if MGONGPU_FPTYPE_DOUBLE
  a = cuCsub( a, b );
#elif defined MGONGPU_FPTYPE_FLOAT
  a = cuCsubf( a, b );
#endif
  return a;
}

SYCL_EXTERNAL inline cxtype
operator*( const cxtype& a, const cxtype& b )
{
#if MGONGPU_FPTYPE_DOUBLE
  return cuCmul( a, b );
#elif defined MGONGPU_FPTYPE_FLOAT
  return cuCmulf( a, b );
#endif
}

SYCL_EXTERNAL inline cxtype
operator/( const cxtype& a, const cxtype& b )
{
#if MGONGPU_FPTYPE_DOUBLE
  return cuCdiv( a, b );
#elif defined MGONGPU_FPTYPE_FLOAT
  return cuCdivf( a, b );
#endif
}

SYCL_EXTERNAL inline cxtype
operator+( const cxtype a )
{
  return a;
}

SYCL_EXTERNAL inline cxtype
operator-( const cxtype& a )
{
  return cxmake( -cxreal( a ), -cximag( a ) );
}

SYCL_EXTERNAL inline cxtype
operator+( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) + b;
}

SYCL_EXTERNAL inline cxtype
operator-( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) - b;
}

SYCL_EXTERNAL inline cxtype
operator*( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) * b;
}

SYCL_EXTERNAL inline cxtype
operator/( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) / b;
}

SYCL_EXTERNAL inline cxtype
operator+( const cxtype& a, const fptype& b )
{
  return a + cxmake( b, 0 );
}

SYCL_EXTERNAL inline cxtype
operator-( const cxtype& a, const fptype& b )
{
  return a - cxmake( b, 0 );
}

SYCL_EXTERNAL inline cxtype
operator*( const cxtype& a, const fptype& b )
{
  return a * cxmake( b, 0 );
}

SYCL_EXTERNAL inline cxtype
operator/( const cxtype& a, const fptype& b )
{
  return a / cxmake( b, 0 );
}
#endif

SYCL_EXTERNAL inline
cxtype cxconj( const cxtype& c )
{
  #ifdef MGONGPU_COMPLEX_EXTRAS
      return extras::conj( c ); // extras::conj( extras::complex<fptype> )
  #endif

  #ifdef MGONGPU_COMPLEX_CXSMPL
      return mgOnGpu::conj( c );
  #endif

  #ifdef MGONGPU_COMPLEX_STD
      return std::conj( c ); // std::conj( std::complex<fptype> )
  #endif

  #ifdef MGONGPU_COMPLEX_ONEAPI
      return sycl::ext::oneapi::experimental::conj( c ); // sycl::ext::oneapi::experimental::conj( sycl::ext::oneapi::experimental::complex<fptype> )
  #endif

  #ifdef MGONGPU_COMPLEX_CUTHRUST
      return thrust::conj( c );
  #endif
}

SYCL_EXTERNAL inline
const cxtype& cxmake( const cxtype& c )
{
  return c;
}

#ifndef MGONGPU_COMPLEX_CXSMPL
inline // NOT __device__
cxtype cxmake( const mgOnGpu::cxsmpl<fptype>& c ) // mgOnGpu::cxsmpl to extras::complex (float-to-float or double-to-double)
{
  return cxmake( c.real(), c.imag() );
}
#endif

#if defined MGONGPU_FPTYPE_FLOAT
inline
cxtype cxmake( const mgOnGpu::cxsmpl<double>& c ) // mgOnGpu::cxsmpl to mgOnGpu::cxsmpl (cast double-to-float)
{
  return cxmake( (fptype)c.real(), (fptype)c.imag() );
}
#endif

//==========================================================================

//==========================================================================
// COMPLEX TYPES: WRAPPER OVER RI FLOATING POINT PAIR (cxtype_ref)
//==========================================================================

namespace mgOnGpu
{
  // The cxtype_ref class (a non-const reference to two fp variables) was originally designed for cxtype_v::operator[]
  // It used to be included in the code only when MGONGPU_HAS_CPPCXTYPEV_BRK (originally MGONGPU_HAS_CPPCXTYPE_REF) is defined
  // It is now always included in the code because it is needed also to access an fptype wavefunction buffer as a cxtype
  class cxtype_ref
  {
  public:
    cxtype_ref() = delete;
    cxtype_ref( const cxtype_ref& ) = delete;
    cxtype_ref( cxtype_ref&& ) = default;
    cxtype_ref( fptype& r, fptype& i ) : m_real( r ), m_imag( i ) {}
    cxtype_ref& operator=( const cxtype_ref& ) = delete;
    cxtype_ref& operator=( cxtype_ref&& c ) { m_real = cxreal( c ); m_imag = cximag( c ); return *this; } // for cxternary
    cxtype_ref& operator=( const cxtype& c ) { m_real = cxreal( c ); m_imag = cximag( c ); return *this; }
    SYCL_EXTERNAL operator cxtype() const { return cxmake( m_real, m_imag ); }
  private:
    fptype &m_real, &m_imag; // RI
  };
}

//// Printout to stream for user defined types
//SYCL_EXTERNAL inline
//std::ostream& operator<<( std::ostream& out, const cxtype& c ) {
//    out << "[" << cxreal(c) << "," << cximag(c) << "]";
//    return out;
//}
//
//SYCL_EXTERNAL inline
//std::ostream& operator<<( std::ostream& out, const mgOnGpu::cxtype_ref& c ){ out << (cxtype)c; return out; }

//==========================================================================

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
#endif // MGONGPUTYPES_H
