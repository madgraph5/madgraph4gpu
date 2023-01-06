#ifndef MGONGPUTYPES_H
#define MGONGPUTYPES_H 1

#include <CL/sycl.hpp>
#if defined MGONGPU_COMPLEX_EXTRAS || defined MGONGPU_COMPLEX_EXTRASVEC
    #include "extras.h"
#endif

#if defined MGONGPU_COMPLEX_STD || defined MGONGPU_COMPLEX_STDVEC
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
  #if defined MGONGPU_COMPLEX_CXSMPL || defined MGONGPU_COMPLEX_CXSMPLVEC
      typedef mgOnGpu::cxsmpl<fptype> cxtype;
  #endif

  #if defined MGONGPU_COMPLEX_EXTRAS || defined MGONGPU_COMPLEX_EXTRASVEC
      typedef extras::complex<fptype> cxtype;
  #endif

  #if defined MGONGPU_COMPLEX_STD || defined MGONGPU_COMPLEX_STDVEC
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

//------------------------------
// Floating point types - SYCL
//------------------------------
// sqrt(2)
#define SQRT2_F 0x1.6a09e6p+0
#define SQRT2_D 0x1.6a09e667f3bcdp+0
#if defined MGONGPU_FPTYPE_DOUBLE
    #define SQRT2 SQRT2_D
#elif defined MGONGPU_FPTYPE_FLOAT
    #define SQRT2 SQRT2_F
#endif

// sqrt(1/2)
#define SQRTH_F 0x1.6a09e6p-1
#define SQRTH_D 0x1.6a09e667f3bcdp-1
#if defined MGONGPU_FPTYPE_DOUBLE
    #define SQRTH SQRTH_D
#elif defined MGONGPU_FPTYPE_FLOAT
    #define SQRTH SQRTH_F
#endif

// fpzero
#define FPZERO fptype(0.0)
#define FPZERO_SV fptype_sv(0.0)

// fpone
#define FPONE fptype(1.0)
#define FPONE_SV fptype_sv(1.0)

// fpmax FIXME add description
#if defined MGONGPU_COMPLEX_CXSMPLVEC || defined MGONGPU_COMPLEX_EXTRASVEC || defined MGONGPU_COMPLEX_STDVEC
    #define FPMAX(a, b) sycl::fmax(a, b)
#else
    #define FPMAX(a, b) (b < a) ? a : b
#endif

// fpmin FIXME add description
#if defined MGONGPU_COMPLEX_CXSMPLVEC || defined MGONGPU_COMPLEX_EXTRASVEC || defined MGONGPU_COMPLEX_STDVEC
    #define FPMIN(a, b) sycl::fmin(a, b)
#else
    #define FPMIN(a, b) (a < b) ? a : b
#endif

// fpsqrt FIXME add description
#define FPSQRT(a) sycl::sqrt(a)

// fpternary or fpconditional FIXME add better description
// fpconditional(a, b, c) = c ? b : a
#define FPCONDITIONAL(a, b, c) c ? b : a
#if defined MGONGPU_COMPLEX_CXSMPLVEC || defined MGONGPU_COMPLEX_EXTRASVEC || defined MGONGPU_COMPLEX_STDVEC
    #define FPCONDITIONAL_SV(a, b, c) sycl::select(a, b, c)
#else
    #define FPCONDITIONAL_SV(a, b, c) FPCONDITIONAL(a, b, c)
#endif

// fpany
#define FPANY(a) a
#if defined MGONGPU_COMPLEX_CXSMPLVEC || defined MGONGPU_COMPLEX_EXTRASVEC || defined MGONGPU_COMPLEX_STDVEC
    #if MGONGPU_MARRAY_DIM == 1
    //    #define FPANY_SV(a) FPANY(a)
        #define FPANY_SV(a) sycl::any((a)[0])
    #else
        #define FPANY_SV(a) sycl::any(a)
    #endif
#else
    #define FPANY_SV(a) FPANY(a)
#endif

//==========================================================================

//==========================================================================
// COMPLEX TYPES: (PLATFORM-SPECIFIC) FUNCTIONS AND OPERATORS
//==========================================================================

// cxmake FIXME add description
#if defined MGONGPU_COMPLEX_CUCOMPLEX && MGONGPU_FPTYPE_DOUBLE
    #define CXMAKE_1ARG(r) make_cuDoubleComplex(r, 0.0)
    #define CXMAKE_2ARG(r, i) make_cuDoubleComplex(r, i)
#elif defined MGONGPU_COMPLEX_CUCOMPLEX && MGONGPU_FPTYPE_FLOAT
    #define CXMAKE_1ARG(r) make_cuFloatComplex(r, 0.0)
    #define CXMAKE_2ARG(r, i) make_cuFloatComplex(r, i)
#else
    #define CXMAKE_1ARG(r) cxtype(r)
    #define CXMAKE_2ARG(r, i) cxtype(r, i)
#endif

// cxreal FIXME add description
#if defined MGONGPU_COMPLEX_CUCOMPLEX && MGONGPU_FPTYPE_DOUBLE
  #define CXREAL(c) cuCreal(c)
#elif defined MGONGPU_COMPLEX_CUCOMPLEX && MGONGPU_FPTYPE_FLOAT
  #define CXREAL(c) cuCrealf(c)
#else
  #define CXREAL(c) c.real()
#endif

#if defined MGONGPU_COMPLEX_CUCOMPLEX && MGONGPU_FPTYPE_DOUBLE
  #define CXIMAG(c) cuCimag(c)
#elif defined MGONGPU_COMPLEX_CUCOMPLEX && MGONGPU_FPTYPE_FLOAT
  #define CXIMAG(c) cuCimagf(c)
#else
  #define CXIMAG(c) c.imag()
#endif

#if defined MGONGPU_COMPLEX_CUCOMPLEX
SYCL_EXTERNAL inline
cxtype operator + (const cxtype& a,const cxtype& b) {
    #if MGONGPU_FPTYPE_DOUBLE
        return cuCadd(a, b);
    #elif defined MGONGPU_FPTYPE_FLOAT
        return cuCaddf(a, b);
    #endif
}

SYCL_EXTERNAL inline
cxtype& operator += (cxtype& a, const cxtype& b) {
    #if MGONGPU_FPTYPE_DOUBLE
        a = cuCadd(a, b);
    #elif defined MGONGPU_FPTYPE_FLOAT
        a = cuCaddf(a, b);
    #endif
    return a;
}

SYCL_EXTERNAL inline
cxtype operator - (const cxtype& a, const cxtype& b) {
    #if MGONGPU_FPTYPE_DOUBLE
        return cuCsub(a, b);
    #elif defined MGONGPU_FPTYPE_FLOAT
        return cuCsubf(a, b);
    #endif
}

SYCL_EXTERNAL inline
cxtype& operator -= (cxtype& a, const cxtype& b) {
    #if MGONGPU_FPTYPE_DOUBLE
        a = cuCsub(a, b);
    #elif defined MGONGPU_FPTYPE_FLOAT
        a = cuCsubf(a, b);
    #endif
    return a;
}

SYCL_EXTERNAL inline
cxtype operator * (const cxtype& a, const cxtype& b) {
    #if MGONGPU_FPTYPE_DOUBLE
        return cuCmul(a, b);
    #elif defined MGONGPU_FPTYPE_FLOAT
        return cuCmulf(a, b);
    #endif
}

SYCL_EXTERNAL inline
cxtype operator / (const cxtype& a, const cxtype& b) {
    #if MGONGPU_FPTYPE_DOUBLE
        return cuCdiv(a, b);
    #elif defined MGONGPU_FPTYPE_FLOAT
        return cuCdivf(a, b);
    #endif
}

SYCL_EXTERNAL inline
cxtype operator + (const cxtype a) {
    return a;
}

SYCL_EXTERNAL inline
cxtype operator - (const cxtype& a) {
    return CXMAKE_2ARG(-CXREAL(a), -CXIMAG(a));
}

SYCL_EXTERNAL inline
cxtype operator + (const fptype& a, const cxtype& b) {
    return CXMAKE_2ARG(a, 0.0) + b;
}

SYCL_EXTERNAL inline
cxtype operator - (const fptype& a, const cxtype& b) {
    return CXMAKE_2ARG(a, 0.0) - b;
}

SYCL_EXTERNAL inline
cxtype operator * (const fptype& a, const cxtype& b) {
    return CXMAKE_2ARG(a, 0.0)*b;
}

SYCL_EXTERNAL inline
cxtype operator / (const fptype& a, const cxtype& b) {
    return CXMAKE_2ARG(a, 0.0)/b;
}

SYCL_EXTERNAL inline
cxtype operator + (const cxtype& a, const fptype& b) {
    return a + CXMAKE_2ARG(b, 0.0);
}

SYCL_EXTERNAL inline
cxtype operator - (const cxtype& a, const fptype& b) {
    return a - CXMAKE_2ARG(b, 0.0);
}

SYCL_EXTERNAL inline
cxtype operator * (const cxtype& a, const fptype& b) {
    return a*CXMAKE_2ARG(b, 0.0);
}

SYCL_EXTERNAL inline
cxtype operator / (const cxtype& a, const fptype& b) {
    return a/CXMAKE_2ARG(b, 0.0);
}
#endif

// cxconj FIXME better description
#if defined MGONGPU_COMPLEX_EXTRAS || defined MGONGPU_COMPLEX_EXTRASVEC
    #define CXCONJ(c) extras::conj(c)
#endif

#if defined MGONGPU_COMPLEX_CXSMPL || defined MGONGPU_COMPLEX_CXSMPLVEC
    #define CXCONJ(c) mgOnGpu::conj(c)
#endif

#ifdef MGONGPU_COMPLEX_STD || defined MGONGPU_COMPLEX_STDVEC
    #define CXCONJ(c) std::conj(c)
#endif

#ifdef MGONGPU_COMPLEX_ONEAPI
    #define CXCONJ(c) sycl::ext::oneapi::experimental::conj(c)
#endif

#ifdef MGONGPU_COMPLEX_CUTHRUST
    #define CXCONJ(c) thrust::conj(c)
#endif

// cxmake_sv FIXME better description
#define CXMAKE_SV_1ARG(a) cxtype_sv(a)
#define CXMAKE_SV_2ARG(a, b) cxtype_sv(a, b)

// cxzero
#define CXZERO CXMAKE_2ARG(FPZERO, FPZERO)
#define CXZERO_SV CXMAKE_SV_2ARG(FPZERO_SV, FPZERO_SV)

// cI imaginary i
#define CXIMAGINARYI CXMAKE_2ARG(FPZERO, FPONE)
#define CXIMAGINARYI_SV CXMAKE_SV_2ARG(FPZERO_SV, FPONE_SV)

// cxternary or cxconditional FIXME better description 
// cxconditional(a, b, c) = c ? b : a
#define CXCONDITIONAL(a, b, c) c ? b : a
#if defined MGONGPU_COMPLEX_CXSMPLVEC || defined MGONGPU_COMPLEX_EXTRASVEC || defined MGONGPU_COMPLEX_STDVEC
    #define CXCONDITIONAL_SV(a, b, c) CXMAKE_SV_2ARG(FPCONDITIONAL_SV(CXREAL(a), CXREAL(b), c), FPCONDITIONAL_SV(CXIMAG(a), CXIMAG(b), c))
#else
    #define CXCONDITIONAL_SV(a, b, c) CXCONDITIONAL(a, b, c)
#endif

// cxabs2
#define CXABS2(c) CXREAL(c)*CXREAL(c) + CXIMAG(c)*CXIMAG(c)

#if not defined MGONGPU_COMPLEX_CXSMPL or not defined MGONGPU_COMPLEX_CXSMPLVEC
inline // NOT __device__
cxtype cxmake( const mgOnGpu::cxsmpl<fptype>& c ) { // mgOnGpu::cxsmpl to cxtype (float-to-float or double-to-double)
    return CXMAKE_2ARG( c.real(), c.imag() );
}
#endif

#if defined MGONGPU_FPTYPE_FLOAT
inline
cxtype cxmake( const mgOnGpu::cxsmpl<double>& c ) { // mgOnGpu::cxsmpl to mgOnGpu::cxsmpl (cast double-to-float)
  return CXMAKE_2ARG( (fptype)c.real(), (fptype)c.imag() );
}
#endif

//==========================================================================

#endif // MGONGPUTYPES_H
