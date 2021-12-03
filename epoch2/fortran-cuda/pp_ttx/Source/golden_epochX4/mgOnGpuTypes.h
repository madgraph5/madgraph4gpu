#ifndef MGONGPUTYPES_H
#define MGONGPUTYPES_H 1

#include "mgOnGpuConfig.h"

#ifdef __CUDACC__
#if defined MGONGPU_CXTYPE_THRUST
#include <thrust/complex.h>
#elif defined MGONGPU_CXTYPE_CUCOMPLEX
#include <cuComplex.h>

#include <complex>
#endif
#else
#include <cmath>
#include <complex>
#endif

namespace mgOnGpu {
// --- Type definitions

// Complex type: cxtype
#ifdef __CUDACC__ // cuda
#if defined MGONGPU_CXTYPE_THRUST
typedef thrust::complex<fptype> cxtype; // two doubles: RI
#elif defined MGONGPU_FPTYPE_DOUBLE
typedef cuDoubleComplex cxtype;
#elif defined MGONGPU_FPTYPE_FLOAT
typedef cuFloatComplex cxtype;
#endif
#else // c++
typedef std::complex<fptype> cxtype; // two doubles: RI
#endif

} // namespace mgOnGpu

// Expose typedefs and operators outside the namespace
using mgOnGpu::cxtype;

// --- Functions and operators for floating point types

#ifdef __CUDACC__ // cuda

/*
inline __host__ __device__
fptype fpmax( const fptype& a, const fptype& b )
{
  return max( a, b );
}

inline __host__ __device__
fptype fpmin( const fptype& a, const fptype& b )
{
  return min( a, b );
}
*/

inline __host__ __device__ const fptype &fpmax(const fptype &a,
                                               const fptype &b) {
  return ((b < a) ? a : b);
}

inline __host__ __device__ const fptype &fpmin(const fptype &a,
                                               const fptype &b) {
  return ((a < b) ? a : b);
}

inline __host__ __device__ fptype fpsqrt(const fptype &f) {
#if defined MGONGPU_FPTYPE_FLOAT
  // See
  // https://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__SINGLE.html
  return sqrtf(f);
#else
  // See
  // https://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__DOUBLE.html
  return sqrt(f);
#endif
}

#else // c++

inline const fptype &fpmax(const fptype &a, const fptype &b) {
  return std::max(a, b);
}

inline const fptype &fpmin(const fptype &a, const fptype &b) {
  return std::min(a, b);
}

inline fptype fpsqrt(const fptype &f) { return std::sqrt(f); }

#endif

// --- Functions and operators for complex types

//------------------------------
// CUDA
//------------------------------

#ifdef __CUDACC__ // cuda

//------------------------------
// CUDA - using thrust::complex
//------------------------------

#if defined MGONGPU_CXTYPE_THRUST // cuda + thrust

//+++++++++++++++++++++++++
// thrust::complex<double>
// thrust::complex<float>
//+++++++++++++++++++++++++

inline __host__ __device__ cxtype cxmake(const fptype &r, const fptype &i) {
  return cxtype(r, i); // thrust::complex<fptype> constructor
}

inline __host__ __device__ fptype cxreal(const cxtype &c) {
  return c.real(); // thrust::complex<fptype>::real()
}

inline __host__ __device__ fptype cximag(const cxtype &c) {
  return c.imag(); // thrust::complex<fptype>::imag()
}

inline __host__ __device__ cxtype cxconj(const cxtype &c) {
  return conj(c); // conj( thrust::complex<fptype> )
}

inline __host__ __device__ const cxtype &cxmake(const cxtype &c) { return c; }

//------------------------------
// CUDA - using cuComplex
//------------------------------

#elif defined MGONGPU_CXTYPE_CUCOMPLEX // cuda + cucomplex

//+++++++++++++++++++++++++
// cuDoubleComplex
//+++++++++++++++++++++++++

#if defined MGONGPU_FPTYPE_DOUBLE // cuda + cucomplex + double

inline __host__ __device__ cxtype cxmake(const fptype &r, const fptype &i) {
  return make_cuDoubleComplex(r, i);
}

inline __host__ __device__ fptype cxreal(const cxtype &c) {
  return cuCreal(c); // returns by value
}

inline __host__ __device__ fptype cximag(const cxtype &c) {
  return cuCimag(c); // returns by value
}

inline __host__ __device__ cxtype operator+(const cxtype &a, const cxtype &b) {
  return cuCadd(a, b);
}

inline __host__ __device__ cxtype &operator+=(cxtype &a, const cxtype &b) {
  a = cuCadd(a, b);
  return a;
}

inline __host__ __device__ cxtype operator-(const cxtype &a, const cxtype &b) {
  return cuCsub(a, b);
}

inline __host__ __device__ cxtype &operator-=(cxtype &a, const cxtype &b) {
  a = cuCsub(a, b);
  return a;
}

inline __host__ __device__ cxtype operator*(const cxtype &a, const cxtype &b) {
  return cuCmul(a, b);
}

inline __host__ __device__ cxtype operator/(const cxtype &a, const cxtype &b) {
  return cuCdiv(a, b);
}

//+++++++++++++++++++++++++
// cuFloatComplex
//+++++++++++++++++++++++++

#elif defined MGONGPU_FPTYPE_FLOAT // cuda + cucomplex + float

inline __host__ __device__ cxtype cxmake(const fptype &r, const fptype &i) {
  return make_cuFloatComplex(r, i);
}

inline __host__ __device__ fptype cxreal(const cxtype &c) {
  return cuCrealf(c); // returns by value
}

inline __host__ __device__ fptype cximag(const cxtype &c) {
  return cuCimagf(c); // returns by value
}

inline __host__ __device__ cxtype operator+(const cxtype &a, const cxtype &b) {
  return cuCaddf(a, b);
}

inline __host__ __device__ cxtype &operator+=(cxtype &a, const cxtype &b) {
  a = cuCaddf(a, b);
  return a;
}

inline __host__ __device__ cxtype operator-(const cxtype &a, const cxtype &b) {
  return cuCsubf(a, b);
}

inline __host__ __device__ cxtype &operator-=(cxtype &a, const cxtype &b) {
  a = cuCsubf(a, b);
  return a;
}

inline __host__ __device__ cxtype operator*(const cxtype &a, const cxtype &b) {
  return cuCmulf(a, b);
}

inline __host__ __device__ cxtype operator/(const cxtype &a, const cxtype &b) {
  return cuCdivf(a, b);
}

inline __host__ // NOT __device__
    cxtype
    cxmake(const std::complex<double>
               &c) // std::complex to cucomplex (cast double-to-float)
{
  return cxmake((fptype)c.real(), (fptype)c.imag());
}

//+++++++++++++++++++++++++
// cuDoubleComplex
// cuFloatComplex
//+++++++++++++++++++++++++

#endif // END cuda + cucomplex + double/float

inline __host__ __device__ cxtype operator+(const cxtype a) { return a; }

inline __host__ __device__ cxtype operator-(const cxtype &a) {
  return cxmake(-cxreal(a), -cximag(a));
}

inline __host__ __device__ cxtype operator+(const fptype &a, const cxtype &b) {
  return cxmake(a, 0) + b;
}

inline __host__ __device__ cxtype operator-(const fptype &a, const cxtype &b) {
  return cxmake(a, 0) - b;
}

inline __host__ __device__ cxtype operator*(const fptype &a, const cxtype &b) {
  return cxmake(a, 0) * b;
}

inline __host__ __device__ cxtype operator/(const fptype &a, const cxtype &b) {
  return cxmake(a, 0) / b;
}

inline __host__ __device__ cxtype operator+(const cxtype &a, const fptype &b) {
  return a + cxmake(b, 0);
}

inline __host__ __device__ cxtype operator-(const cxtype &a, const fptype &b) {
  return a - cxmake(b, 0);
}

inline __host__ __device__ cxtype operator*(const cxtype &a, const fptype &b) {
  return a * cxmake(b, 0);
}

inline __host__ __device__ cxtype operator/(const cxtype &a, const fptype &b) {
  return a / cxmake(b, 0);
}

inline __host__ __device__ cxtype cxconj(const cxtype &c) {
  return cxmake(cxreal(c), -cximag(c));
}

inline __host__ // NOT __device__
    cxtype
    cxmake(
        const std::complex<fptype> &
            c) // std::complex to cucomplex (float-to-float or double-to-double)
{
  return cxmake(c.real(), c.imag());
}

#endif // END cuda + thrust/cucomplex

//------------------------------
// C++ - using std::complex
//------------------------------

#else // c++

//+++++++++++++++++++++++++
// std::complex<float>
// std::complex<double>
//+++++++++++++++++++++++++

inline cxtype cxmake(const fptype &r, const fptype &i) {
  return cxtype(r, i); // std::complex<fptype> constructor
}

inline fptype cxreal(const cxtype &c) {
  return c.real(); // std::complex<fptype>::real()
}

inline fptype cximag(const cxtype &c) {
  return c.imag(); // std::complex<fptype>::imag()
}

inline cxtype cxconj(const cxtype &c) {
  return conj(c); // conj( std::complex<fptype> )
}

inline const cxtype &
cxmake(const cxtype &c) // std::complex to std::complex (float-to-float or
                        // double-to-double)
{
  return c;
}

#if defined MGONGPU_FPTYPE_FLOAT
inline cxtype
cxmake(const std::complex<double>
           &c) // std::complex to std::complex (cast double-to-float)
{
  return cxmake((fptype)c.real(), (fptype)c.imag());
}
#endif

#endif // END cuda/c++

#endif // MGONGPUTYPES_H
