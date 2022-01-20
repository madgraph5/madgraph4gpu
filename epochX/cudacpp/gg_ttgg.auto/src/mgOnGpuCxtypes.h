#ifndef MGONGPUCXTYPES_H
#define MGONGPUCXTYPES_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuFptypes.h"

//==========================================================================
// COMPLEX TYPES: (PLATFORM-SPECIFIC) HEADERS
//==========================================================================

#include <complex>

// Complex type in cuda: thrust or cucomplex or cxsmpl
#ifdef __CUDACC__
#if defined MGONGPU_CUCXTYPE_THRUST
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtautological-compare" // for icpx2021/clang13 (https://stackoverflow.com/a/15864661)
#include <thrust/complex.h>
#pragma clang diagnostic pop
#elif defined MGONGPU_CUCXTYPE_CUCOMPLEX
#include <cuComplex.h>
#elif not defined MGONGPU_CUCXTYPE_CXSMPL
#error You must CHOOSE (ONE AND) ONLY ONE of MGONGPU_CUCXTYPE_THRUST or MGONGPU_CUCXTYPE_CUCOMPLEX or MGONGPU_CUCXTYPE_CXSMPL
#endif
#else
// Complex type in c++: std::complex or cxsmpl
#if defined MGONGPU_CPPCXTYPE_STDCOMPLEX
#include <cmath>
#elif not defined MGONGPU_CPPCXTYPE_CXSMPL
#error You must CHOOSE (ONE AND) ONLY ONE of MGONGPU_CPPCXTYPE_STDCOMPLEX or MGONGPU_CPPCXTYPE_CXSMPL
#endif
#endif

//==========================================================================
// COMPLEX TYPES: SIMPLE COMPLEX CLASS (cxsmpl)
//==========================================================================

#if defined MGONGPU_CUCXTYPE_CXSMPL or defined MGONGPU_CPPCXTYPE_CXSMPL

namespace mgOnGpu
{
  // --- Type definition (simple complex type derived from cxtype_v)
  class cxsmpl
  {
  public:
    __host__ __device__ constexpr cxsmpl() : m_real{0}, m_imag{0} {}
    cxsmpl( const cxsmpl&  ) = default;
    cxsmpl( cxsmpl&&  ) = default;
    __host__ __device__ constexpr cxsmpl( const fptype& r, const fptype& i = 0 ) : m_real{r}, m_imag{i} {}
    cxsmpl& operator=( const cxsmpl& ) = default;
    cxsmpl& operator=( cxsmpl&& ) = default;
    __host__ __device__ constexpr cxsmpl& operator+=( const cxsmpl& c ){ m_real += c.real(); m_imag += c.imag(); return *this; }
    __host__ __device__ constexpr cxsmpl& operator-=( const cxsmpl& c ){ m_real -= c.real(); m_imag -= c.imag(); return *this; }
    __host__ __device__ constexpr const fptype& real() const { return m_real; }
    __host__ __device__ constexpr const fptype& imag() const { return m_imag; }
  private:
    fptype m_real, m_imag; // RI
  };
  inline __host__ __device__ //constexpr (NB: a constexpr function cannot have a nonliteral return type "mgOnGpu::cxsmpl")
  cxsmpl conj( const cxsmpl& c ){ return cxsmpl( c.real(), -c.imag() ); }
}
#endif

//==========================================================================
// COMPLEX TYPES: (PLATFORM-SPECIFIC) TYPEDEFS
//==========================================================================

namespace mgOnGpu
{

  // --- Type definitions (complex type: cxtype)
#ifdef __CUDACC__ // cuda
#if defined MGONGPU_CUCXTYPE_THRUST
  typedef thrust::complex<fptype> cxtype;
#elif defined MGONGPU_CUCXTYPE_CUCOMPLEX
#if defined MGONGPU_FPTYPE_DOUBLE
  typedef cuDoubleComplex cxtype;
#elif defined MGONGPU_FPTYPE_FLOAT
  typedef cuFloatComplex cxtype;
#endif
#else
  typedef cxsmpl cxtype;
#endif
#else // c++
#if defined MGONGPU_CPPCXTYPE_STDCOMPLEX
  typedef std::complex<fptype> cxtype;
#else
  typedef cxsmpl cxtype;
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

#if defined MGONGPU_CUCXTYPE_CXSMPL or defined MGONGPU_CPPCXTYPE_CXSMPL

//------------------------------
// CUDA or C++ - using cxsmpl
//------------------------------

inline __host__ __device__
cxtype cxmake( const fptype& r, const fptype& i )
{
  return cxtype( r, i ); // cxsmpl constructor
}

inline __host__ __device__
fptype cxreal( const cxtype& c )
{
  return c.real(); // cxsmpl::real()
}

inline __host__ __device__
fptype cximag( const cxtype& c )
{
  return c.imag(); // cxsmpl::imag()
}

inline __host__ __device__
cxtype cxconj( const cxtype& c )
{
  return conj( c ); // conj( cxsmpl )
}

inline __host__ __device__
const cxtype& cxmake( const cxtype& c )
{
  return c;
}

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
cxtype operator+( const cxtype& a, const cxtype& b )
{
  return cxmake( a.real() + b.real(), a.imag() + b.imag() );
}

inline __host__ __device__
cxtype operator+( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) + b;
}

inline __host__ __device__
cxtype operator-( const cxtype& a, const cxtype& b )
{
  return cxmake( a.real() - b.real(), a.imag() - b.imag() );
}

inline __host__ __device__
cxtype operator-( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) - b;
}

inline __host__ __device__
cxtype operator*( const cxtype& a, const cxtype& b )
{
  return cxmake( a.real() * b.real() - a.imag() * b.imag(), a.imag() * b.real() + a.real() * b.imag() );
}

inline __host__ __device__
cxtype operator*( const fptype& a, const cxtype& b )
{
  return cxmake( a, 0 ) * b;
}

inline __host__ __device__
cxtype operator/( const cxtype& a, const cxtype& b )
{
  fptype bnorm = b.real()*b.real() + b.imag()*b.imag();
  return cxmake( ( a.real() * b.real() + a.imag() * b.imag() ) / bnorm,
                 ( a.imag() * b.real() - a.real() * b.imag() ) / bnorm );
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

inline __host__ // NOT __device__
cxtype cxmake( const std::complex<float>& c ) // std::complex to cxsmpl (float-to-float or float-to-double)
{
  return cxmake( c.real(), c.imag() );
}

inline __host__ // NOT __device__
cxtype cxmake( const std::complex<double>& c ) // std::complex to cxsmpl (double-to-float or double-to-double)
{
  return cxmake( c.real(), c.imag() );
}

#endif // #if defined MGONGPU_CUCXTYPE_CXSMPL or defined MGONGPU_CPPCXTYPE_CXSMPL

//==========================================================================

#if defined __CUDACC__ and defined MGONGPU_CUCXTYPE_THRUST // cuda + thrust

//------------------------------
// CUDA - using thrust::complex
//------------------------------

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
cxtype cxconj( const cxtype& c )
{
  return conj( c ); // conj( thrust::complex<fptype> )
}

inline __host__ __device__
const cxtype& cxmake( const cxtype& c )
{
  return c;
}

#endif // #if defined __CUDACC__ and defined MGONGPU_CUCXTYPE_THRUST

//==========================================================================

#if defined __CUDACC__ and defined MGONGPU_CUCXTYPE_CUCOMPLEX // cuda + cucomplex

//------------------------------
// CUDA - using cuComplex
//------------------------------

#if defined MGONGPU_FPTYPE_DOUBLE  // cuda + cucomplex + double

//+++++++++++++++++++++++++
// cuDoubleComplex ONLY
//+++++++++++++++++++++++++

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
cxtype& operator+=( cxtype& a, const cxtype& b )
{
  a = cuCadd( a, b ); return a;
}

inline __host__ __device__
cxtype operator-( const cxtype& a, const cxtype& b )
{
  return cuCsub( a, b );
}

inline __host__ __device__
cxtype& operator-=( cxtype& a, const cxtype& b )
{
  a = cuCsub( a, b ); return a;
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

#elif defined MGONGPU_FPTYPE_FLOAT  // cuda + cucomplex + float

//+++++++++++++++++++++++++
// cuFloatComplex ONLY
//+++++++++++++++++++++++++

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
cxtype& operator+=( cxtype& a, const cxtype& b )
{
  a = cuCaddf( a, b ); return a;
}

inline __host__ __device__
cxtype operator-( const cxtype& a, const cxtype& b )
{
  return cuCsubf( a, b );
}

inline __host__ __device__
cxtype& operator-=( cxtype& a, const cxtype& b )
{
  a = cuCsubf( a, b ); return a;
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

#endif

//+++++++++++++++++++++++++
// cuDoubleComplex OR
// cuFloatComplex
//+++++++++++++++++++++++++

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
cxtype cxconj( const cxtype& c )
{
  return cxmake( cxreal( c ), -cximag( c ) );
}

inline __host__ // NOT __device__
cxtype cxmake( const std::complex<fptype>& c ) // std::complex to cucomplex (float-to-float or double-to-double)
{
  return cxmake( c.real(), c.imag() );
}

#endif // #if defined __CUDACC__ and defined MGONGPU_CUCXTYPE_CUCOMPLEX

//==========================================================================

#if not defined __CUDACC__ and defined MGONGPU_CPPCXTYPE_STDCOMPLEX // c++ + stdcomplex

//------------------------------
// C++ - using std::complex
//------------------------------

inline
cxtype cxmake( const fptype& r, const fptype& i )
{
  return cxtype( r, i ); // std::complex<fptype> constructor
}

inline
fptype cxreal( const cxtype& c )
{
  return c.real(); // std::complex<fptype>::real()
}

inline
fptype cximag( const cxtype& c )
{
  return c.imag(); // std::complex<fptype>::imag()
}

inline
cxtype cxconj( const cxtype& c )
{
  return conj( c ); // conj( std::complex<fptype> )
}

inline
const cxtype& cxmake( const cxtype& c ) // std::complex to std::complex (float-to-float or double-to-double)
{
  return c;
}

#if defined MGONGPU_FPTYPE_FLOAT
inline
cxtype cxmake( const std::complex<double>& c ) // std::complex to std::complex (cast double-to-float)
{
  return cxmake( (fptype)c.real(), (fptype)c.imag() );
}
#endif

#endif // #if not defined __CUDACC__ and defined MGONGPU_CPPCXTYPE_STDCOMPLEX

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
    cxtype_ref( fptype& r, fptype& i ) : m_real{r}, m_imag{i} {}
    cxtype_ref& operator=( const cxtype_ref& ) = delete;
    cxtype_ref& operator=( cxtype_ref&& c ) { m_real = cxreal( c ); m_imag = cximag( c ); return *this; } // for cxternary
    cxtype_ref& operator=( const cxtype& c ) { m_real = cxreal( c ); m_imag = cximag( c ); return *this; }
    operator cxtype() const { return cxmake( m_real, m_imag ); }
  private:
    fptype &m_real, &m_imag; // RI
  };
}

//==========================================================================

#endif // MGONGPUCXTYPES_H
