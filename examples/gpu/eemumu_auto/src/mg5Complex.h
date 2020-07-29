#ifndef mg5Complex_h
#define mg5Complex_h

#define MG5_COMPLEX_THRUST

#ifdef MG5_COMPLEX_THRUST

#include <thrust/complex.h>

#ifdef MG5_CMPLX_SINGLE_PREC

#define mg5Complex thrust::complex<float>

inline __host__ __device__ thrust::complex<float> make_mg5Complex(float r,
                                                                  float i) {
  return thrust::complex<float>(r, i);
}

inline __host__ thrust::complex<float> make_mg5Complex(std::complex<float> c) {
  return thrust::complex<float>(c.real(), c.imag());
}

inline __host__ __device__ float real(thrust::complex<float> c) {
  return c.real();
}

inline __host__ __device__ float imag(thrust::complex<float> c) {
  return c.imag();
}

#else

#define mg5Complex thrust::complex<double>

inline __host__ __device__ thrust::complex<double> make_mg5Complex(float r,
                                                                   float i) {
  return thrust::complex<double>(r, i);
}

inline __host__ thrust::complex<double>
make_mg5Complex(std::complex<double> c) {
  return thrust::complex<double>(c.real(), c.imag());
}

inline __host__ __device__ float real(thrust::complex<double> c) {
  return c.real();
}

inline __host__ __device__ float imag(thrust::complex<double> c) {
  return c.imag();
}

#endif // MG5_CMPLX_SINGLE_PREC

#else // MG5_COMPLEX_THRUST

#include <complex>
#include <cuComplex.h>

#ifdef MG5_CMPLX_SINGLE_PREC

#define mg5Complex cuFloatComplex

inline __host__ __device__ cuFloatComplex make_mg5Complex(float r, float i) {
  return make_cuFloatComplex(r, i);
}

inline __host__ cuFloatComplex make_mg5Complex(std::complex<float> c) {
  return make_cuFloatComplex(c.real(), c.imag());
}

inline __host__ __device__ cuFloatComplex operator-(cuFloatComplex c) {
  return cuCmulf(make_cuFloatComplex(-1, 0), c);
}

inline __host__ __device__ cuFloatComplex operator+(cuFloatComplex c) {
  return c;
}

inline __host__ __device__ cuFloatComplex operator+(cuFloatComplex a,
                                                    cuFloatComplex b) {
  return cuCaddf(a, b);
}

inline __host__ __device__ cuFloatComplex operator+(float a, cuFloatComplex b) {
  return cuCaddf(make_cuFloatComplex(a, 0), b);
}

inline __host__ __device__ cuFloatComplex operator+(cuFloatComplex a, float b) {
  return cuCaddf(a, make_cuFloatComplex(b, 0));
}

inline __host__ __device__ cuFloatComplex operator-(cuFloatComplex a,
                                                    cuFloatComplex b) {
  return cuCsubf(a, b);
}

inline __host__ __device__ cuFloatComplex operator-(float a, cuFloatComplex b) {
  return cuCsubf(make_cuFloatComplex(a, 0), b);
}

inline __host__ __device__ cuFloatComplex operator-(cuFloatComplex a, float b) {
  return cuCsubf(a, make_cuFloatComplex(b, 0));
}

inline __host__ __device__ cuFloatComplex operator*(cuFloatComplex a,
                                                    cuFloatComplex b) {
  return cuCmulf(a, b);
}

inline __host__ __device__ cuFloatComplex operator*(float a, cuFloatComplex b) {
  return cuCmulf(make_cuFloatComplex(a, 0), b);
}

inline __host__ __device__ cuFloatComplex operator*(cuFloatComplex a, float b) {
  return cuCmulf(a, make_cuFloatComplex(b, 0));
}

inline __host__ __device__ cuFloatComplex operator/(cuFloatComplex a,
                                                    cuFloatComplex b) {
  return cuCdivf(a, b);
}

inline __host__ __device__ cuFloatComplex operator/(float a, cuFloatComplex b) {
  return cuCdivf(make_cuFloatComplex(a, 0), b);
}

inline __host__ __device__ cuFloatComplex operator/(cuFloatComplex a, float b) {
  return cuCdivf(a, make_cuFloatComplex(b, 0));
}

inline __host__ __device__ cuFloatComplex conj(cuFloatComplex c) {
  return make_cuFloatComplex(cuCrealf(c), -cuCimagf(c));
}

inline __host__ __device__ float real(cuFloatComplex c) { return cuCrealf(c); }

inline __host__ __device__ float imag(cuFloatComplex c) { return cuCimagf(c); }

#else // MG5_CMPLX_SINGLE_PREC

#define mg5Complex cuDoubleComplex

inline __host__ __device__ cuDoubleComplex make_mg5Complex(double r, double i) {
  return make_cuDoubleComplex(r, i);
}

inline __host__ cuDoubleComplex make_mg5Complex(std::complex<double> c) {
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
