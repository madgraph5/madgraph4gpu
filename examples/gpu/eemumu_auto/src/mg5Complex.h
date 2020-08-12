#ifndef mg5Complex_h
#define mg5Complex_h

// #define MG5_COMPLEX_THRUST
// #define MG5_CMPLX_SINGLE_PREC

#ifdef MG5_CMPLX_SINGLE_PREC
#define floa_t float
#define DECL_NUM(X) X##f
#else
#define floa_t double
#define DECL_NUM(X) X
#endif


#ifdef MG5_COMPLEX_THRUST

#include <thrust/complex.h>

#define mg5Complex thrust::complex<floa_t>

inline __host__ __device__ thrust::complex<floa_t> make_mg5Complex(floa_t r,
                                                                   floa_t i) {
  return thrust::complex<floa_t>(r, i);
}

inline __host__ thrust::complex<floa_t> make_mg5Complex(std::complex<floa_t> c) {
  return thrust::complex<floa_t>(c.real(), c.imag());
}

inline __host__ __device__ floa_t real(thrust::complex<floa_t> c) {
  return c.real();
}

inline __host__ __device__ floa_t imag(thrust::complex<floa_t> c) {
  return c.imag();
}

#else // MG5_COMPLEX_THRUST

#include <complex>
#include <cuComplex.h>

#ifdef MG5_CMPLX_SINGLE_PREC
#define mg5Complex cuFloatComplex
#else
#define mg5Complex cuDoubleComplex
#endif

inline __host__ __device__ mg5Complex make_mg5Complex(floa_t r, floa_t i) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return make_cuFloatComplex(r, i);
#else
  return make_cuDoubleComplex(r, i);
#endif
}

inline __host__ mg5Complex make_mg5Complex(std::complex<floa_t> c) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return make_cuFloatComplex(c.real(), c.imag());
#else
  return make_cuDoubleComplex(c.real(), c.imag());
#endif
}

inline __host__ __device__ mg5Complex operator-(mg5Complex c) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCmulf(make_cuFloatComplex(-1, 0), c);
#else
  return cuCmul(make_cuDoubleComplex(-1, 0), c);
#endif
}

inline __host__ __device__ mg5Complex operator+(mg5Complex c) {
  return c;
}

inline __host__ __device__ mg5Complex operator+(mg5Complex a,
                                                mg5Complex b) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCaddf(a, b);
#else
  return cuCadd(a, b);
#endif
}

inline __host__ __device__ mg5Complex operator+(floa_t a, mg5Complex b) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCaddf(make_cuFloatComplex(a, 0), b);
#else
  return cuCadd(make_cuDoubleComplex(a, 0), b);
#endif
}

inline __host__ __device__ mg5Complex operator+(mg5Complex a, floa_t b) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCaddf(a, make_cuFloatComplex(b, 0));
#else
  return cuCadd(a, make_cuDoubleComplex(b, 0));
#endif
}

inline __host__ __device__ mg5Complex operator-(mg5Complex a,
                                                mg5Complex b) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCsubf(a, b);
#else
  return cuCsub(a, b);
#endif
}

inline __host__ __device__ mg5Complex operator-(floa_t a, mg5Complex b) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCsubf(make_cuFloatComplex(a, 0), b);
#else
  return cuCsub(make_cuDoubleComplex(a, 0), b);
#endif
}

inline __host__ __device__ mg5Complex operator-(mg5Complex a, floa_t b) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCsubf(a, make_cuFloatComplex(b, 0));
#else
  return cuCsub(a, make_cuDoubleComplex(b, 0));
#endif
}

inline __host__ __device__ mg5Complex operator*(mg5Complex a,
                                                mg5Complex b) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCmulf(a, b);
#else
  return cuCmul(a, b);
#endif
}

inline __host__ __device__ mg5Complex operator*(floa_t a, mg5Complex b) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCmulf(make_cuFloatComplex(a, 0), b);
#else
  return cuCmul(make_cuDoubleComplex(a, 0), b);
#endif
}

inline __host__ __device__ mg5Complex operator*(mg5Complex a, floa_t b) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCmulf(a, make_cuFloatComplex(b, 0));
#else
  return cuCmul(a, make_cuDoubleComplex(b, 0));
#endif
}

inline __host__ __device__ mg5Complex operator/(mg5Complex a,
                                                mg5Complex b) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCdivf(a, b);
#else
  return cuCdiv(a, b);
#endif
}

inline __host__ __device__ mg5Complex operator/(floa_t a, mg5Complex b) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCdivf(make_cuFloatComplex(a, 0), b);
#else
  return cuCdiv(make_cuDoubleComplex(a, 0), b);
#endif
}

inline __host__ __device__ mg5Complex operator/(mg5Complex a, floa_t b) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCdivf(a, make_cuFloatComplex(b, 0));
#else
  return cuCdiv(a, make_cuDoubleComplex(b, 0));
#endif
}

inline __host__ __device__ mg5Complex conj(mg5Complex c) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return make_cuFloatComplex(cuCrealf(c), -cuCimagf(c));
#else
  return make_cuDoubleComplex(cuCreal(c), -cuCimag(c));
#endif
}

inline __host__ __device__ floa_t real(mg5Complex c) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCrealf(c);
#else
  return cuCreal(c);
#endif
}

inline __host__ __device__ floa_t imag(mg5Complex c) {
#ifdef MG5_CMPLX_SINGLE_PREC
  return cuCimagf(c);
#else
  return cuCimag(c);
#endif
}

#endif // MG5_COMPLEX_THRUST

#endif // mg5Complex_h
