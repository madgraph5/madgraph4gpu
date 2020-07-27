#ifndef mg5Complex_h
#define mg5Complex_h

#ifdef MG5_COMPLEX_THRUST

#include <thrust/complex.h>

#ifdef MG5_CMPLX_SINGLE_PREC
#define mg5Complex thrust::complex<float>
#else
#define mg5Complex thrust::complex<double>
#endif // MG5_CMPLX_SINGLE_PREC

#else // MG5_COMPLEX_THRUST

#include <complex>
#include <cuComplex.h>

#ifdef MG5_CMPLX_SINGLE_PREC

#define mg5Complex cuFloatComplex

inline __host__ __device__ mg5Complex mg5Complex(double r, double i) {
  return make_cuFloatComplex(r, i);
}

inline __host__ __device__ cuFloatComplex operator+(cuFloatComplex a,
                                                    cuFloatComplex b) {
  return cuCaddf(a, b);
}

inline __host__ __device__ cuFloatComplex operator-(cuFloatComplex a,
                                                    cuFloatComplex b) {
  return cuCsubf(a, b);
}

inline __host__ __device__ cuFloatComplex operator*(cuFloatComplex a,
                                                    cuFloatComplex b) {
  return cuCmulf(a, b);
}

inline __host__ __device__ cuFloatComplex operator/(cuFloatComplex a,
                                                    cuFloatComplex b) {
  return cuCdivf(a, b);
}

#else // MG5_CMPLX_SINGLE_PREC

#define mg5Complex cuDoubleComplex

inline __host__ __device__ cuDoubleComplex make_mg5Complex(double r, double i) {
  return make_cuDoubleComplex(r, i);
}

inline __host__ __device__ cuDoubleComplex
make_mg5Complex2(std::complex<double> c) {
  return make_cuDoubleComplex(c.real(), c.imag());
}

inline __host__ __device__ cuDoubleComplex operator-(cuDoubleComplex c) {
  return cuCmul(make_cuDoubleComplex(-1, 0), c);
}

inline __host__ __device__ cuDoubleComplex operator+(cuDoubleComplex c) {
  return c;
}

inline __host__ __device__ cuDoubleComplex operator+(cuDoubleComplex a,
                                                     cuDoubleComplex b) {
  return cuCadd(a, b);
}

inline __host__ __device__ cuDoubleComplex operator+(double a,
                                                     cuDoubleComplex b) {
  return cuCadd(make_cuDoubleComplex(a, 0), b);
}

inline __host__ __device__ cuDoubleComplex operator+(cuDoubleComplex a,
                                                     double b) {
  return cuCadd(a, make_cuDoubleComplex(b, 0));
}

inline __host__ __device__ cuDoubleComplex operator-(cuDoubleComplex a,
                                                     cuDoubleComplex b) {
  return cuCsub(a, b);
}

inline __host__ __device__ cuDoubleComplex operator-(double a,
                                                     cuDoubleComplex b) {
  return cuCsub(make_cuDoubleComplex(a, 0), b);
}

inline __host__ __device__ cuDoubleComplex operator-(cuDoubleComplex a,
                                                     double b) {
  return cuCsub(a, make_cuDoubleComplex(b, 0));
}

inline __host__ __device__ cuDoubleComplex operator*(cuDoubleComplex a,
                                                     cuDoubleComplex b) {
  return cuCmul(a, b);
}

inline __host__ __device__ cuDoubleComplex operator*(double a,
                                                     cuDoubleComplex b) {
  return cuCmul(make_cuDoubleComplex(a, 0), b);
}

inline __host__ __device__ cuDoubleComplex operator*(cuDoubleComplex a,
                                                     double b) {
  return cuCmul(a, make_cuDoubleComplex(b, 0));
}

inline __host__ __device__ cuDoubleComplex operator/(cuDoubleComplex a,
                                                     cuDoubleComplex b) {
  return cuCdiv(a, b);
}

inline __host__ __device__ cuDoubleComplex operator/(double a,
                                                     cuDoubleComplex b) {
  return cuCdiv(make_cuDoubleComplex(a, 0), b);
}

inline __host__ __device__ cuDoubleComplex operator/(cuDoubleComplex a,
                                                     double b) {
  return cuCdiv(a, make_cuDoubleComplex(b, 0));
}

inline __host__ __device__ cuDoubleComplex conj(cuDoubleComplex c) {
  return make_cuDoubleComplex(cuCreal(c), -cuCimag(c));
}

inline __host__ __device__ double real(cuDoubleComplex c) { return cuCreal(c); }

inline __host__ __device__ double imag(cuDoubleComplex c) { return cuCimag(c); }

#endif // MG5_CMPLX_SINGLE_PREC

#endif // MG5_COMPLEX_THRUST

#endif // mg5Complex_h
