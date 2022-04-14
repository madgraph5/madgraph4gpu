#ifndef MGONGPUVECTORS_H
#define MGONGPUVECTORS_H 1

#include <iostream>

//------------------------------
// Vector types - C++
//------------------------------

#if !defined(ALPAKA) && !defined(__CUDACC__)

#include "mgOnGpuTypes.h"

namespace mgOnGpu
{
#ifdef MGONGPU_CPPSIMD

  const int neppV = neppM;

  // --- Type definition (using vector compiler extensions: need -march=...)
#ifdef __clang__ // https://clang.llvm.org/docs/LanguageExtensions.html#vectors-and-extended-vectors
  typedef fptype fptype_v __attribute__ ((ext_vector_type(neppV))); // RRRR
#else
  typedef fptype fptype_v __attribute__ ((vector_size (neppV*sizeof(fptype)))); // RRRR
#endif

#ifdef __clang__
  // If set: return a pair of (fptype&, fptype&) by non-const reference in cxtype_v::operator[]
  // This is forbidden in clang ("non-const reference cannot bind to vector element")
  // See also https://stackoverflow.com/questions/26554829
  //#define MGONGPU_HAS_CXTYPE_REF 1 // clang test (compilation fails also on clang 12.0, issue #182)
#undef MGONGPU_HAS_CXTYPE_REF // clang default
#elif defined __INTEL_COMPILER
  //#define MGONGPU_HAS_CXTYPE_REF 1 // icc default?
#undef MGONGPU_HAS_CXTYPE_REF // icc test
#else
#define MGONGPU_HAS_CXTYPE_REF 1 // gcc default
  //#undef MGONGPU_HAS_CXTYPE_REF // gcc test (very slightly slower? issue #172)
#endif

#ifdef MGONGPU_HAS_CXTYPE_REF
  // NB: the alternative "clang" implementation is simpler: it simply does not have class cxtype_ref
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
#endif

  // --- Type definition (using vector compiler extensions: need -march=...)
  class cxtype_v // no need for "class alignas(2*sizeof(fptype_v)) cxtype_v"
  {
  public:
    cxtype_v() : m_real{0}, m_imag{0} {}
    cxtype_v( const cxtype_v&  ) = default;
    cxtype_v( cxtype_v&&  ) = default;
    cxtype_v( const fptype_v& r, const fptype_v& i ) : m_real{r}, m_imag{i} {}
    cxtype_v& operator=( const cxtype_v& ) = default;
    cxtype_v& operator=( cxtype_v&& ) = default;
    cxtype_v& operator+=( const cxtype_v& c ){ m_real += c.real(); m_imag += c.imag(); return *this; }
    cxtype_v& operator-=( const cxtype_v& c ){ m_real -= c.real(); m_imag -= c.imag(); return *this; }
#ifdef MGONGPU_HAS_CXTYPE_REF
    // NB: the alternative "clang" implementation is simpler: it simply does not have any operator[]
    // NB: ** do NOT implement operator[] to return a value: it does not fail the build (why?) and gives unexpected results! **
    cxtype_ref operator[]( size_t i ) const { return cxtype_ref( m_real[i], m_imag[i] ); }
#endif
    const fptype_v& real() const { return m_real; }
    const fptype_v& imag() const { return m_imag; }
  private:
    fptype_v m_real, m_imag; // RRRRIIII
  };

  // --- Type definition (using vector compiler extensions: need -march=...)
#ifdef __clang__ // https://clang.llvm.org/docs/LanguageExtensions.html#vectors-and-extended-vectors
#if defined MGONGPU_FPTYPE_DOUBLE
  typedef long int bool_v __attribute__ ((ext_vector_type(neppV))); // bbbb
#elif defined MGONGPU_FPTYPE_FLOAT
  typedef int bool_v __attribute__ ((ext_vector_type(neppV))); // bbbb
#endif
#else // gcc
#if defined MGONGPU_FPTYPE_DOUBLE
  typedef long int bool_v __attribute__ ((vector_size (neppV*sizeof(long int)))); // bbbb
#elif defined MGONGPU_FPTYPE_FLOAT
  typedef int bool_v __attribute__ ((vector_size (neppV*sizeof(int)))); // bbbb
#endif
#endif

#else // MGONGPU_CPPSIMD not defined

  const int neppV = 1; // Note: also neppM is equal to 1

#endif
}

// Expose typedefs outside the namespace
using mgOnGpu::neppV;
#ifdef MGONGPU_CPPSIMD
using mgOnGpu::fptype_v;
using mgOnGpu::cxtype_v;
using mgOnGpu::bool_v;
#endif

// Printout to stream for user defined types
#ifdef MGONGPU_CPPSIMD
inline std::ostream& operator<<( std::ostream& out, const bool_v& v )
{
  out << "{ " << v[0];
  for ( int i=1; i<neppV; i++ ) std::cout << ", " << v[i];
  out << " }";
  return out;
}

inline std::ostream& operator<<( std::ostream& out, const fptype_v& v )
{
  out << "{ " << v[0];
  for ( int i=1; i<neppV; i++ ) std::cout << ", " << v[i];
  out << " }";
  return out;
}
#endif

inline std::ostream& operator<<( std::ostream& out, const cxtype& c )
{
  out << "[" << cxreal(c) << "," << cximag(c) << "]";
  //out << cxreal(c) << "+i" << cximag(c);
  return out;
}

#ifdef MGONGPU_CPPSIMD
inline std::ostream& operator<<( std::ostream& out, const cxtype_v& v )
{
#ifdef MGONGPU_HAS_CXTYPE_REF
  out << "{ " << v[0];
  for ( int i=1; i<neppV; i++ ) std::cout << ", " << v[i];
#else
  out << "{ " << cxmake( v.real()[0], v.imag()[0] );
  for ( int i=1; i<neppV; i++ ) std::cout << ", " << cxmake( v.real()[i], v.imag()[i] );
#endif
  out << " }";
  return out;
}
#endif

// Printout to std::cout for user defined types
inline void print( const fptype& f ) { std::cout << f << std::endl; }

#ifdef MGONGPU_CPPSIMD
inline void print( const fptype_v& v ) { std::cout << v << std::endl; }
#endif

inline void print( const cxtype& c ) { std::cout << c << std::endl; }

#ifdef MGONGPU_CPPSIMD
inline void print( const cxtype_v& v ) { std::cout << v << std::endl; }
#endif

// Operators for fptype_v
#ifdef MGONGPU_CPPSIMD
inline
fptype_v fpsqrt( const fptype_v& v )
{
  // See https://stackoverflow.com/questions/18921049/gcc-vector-extensions-sqrt
  fptype_v out;
  for ( int i=0; i<neppV; i++ ) out[i]=fpsqrt(v[i]);
  return out;
}

/*
inline
fptype_v fpvmake( const fptype v[neppV] )
{
  fptype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = v[i];
  return out;
}
*/
#endif

// Operators for cxtype_v
#ifdef MGONGPU_CPPSIMD
/*
inline
cxtype_v cxvmake( const cxtype c )
{
  cxtype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = c;
  return out;
}
*/

inline
cxtype_v cxmake( const fptype_v& r, const fptype_v& i )
{
  return cxtype_v{ r, i };
}

inline
cxtype_v cxmake( const fptype_v& r, const fptype& i )
{
  return cxtype_v{ r, fptype_v{i} };
}

inline
cxtype_v cxmake( const fptype& r, const fptype_v& i )
{
  return cxtype_v{ fptype_v{r}, i };
}

inline
const fptype_v& cxreal( const cxtype_v& c )
{
  return c.real(); // returns by reference
}

inline
const fptype_v& cximag( const cxtype_v& c )
{
  return c.imag(); // returns by reference
}

inline
const cxtype_v cxconj( const cxtype_v& c )
{
  return cxmake( c.real(), -c.imag() );
}

inline
cxtype_v operator+( const cxtype_v& a, const cxtype_v& b )
{
  return cxmake( a.real() + b.real(), a.imag() + b.imag() );
}

inline
cxtype_v operator+( const fptype_v& a, const cxtype_v& b )
{
  return cxmake( a + b.real(), b.imag() );
}

inline
cxtype_v operator+( const cxtype_v& a, const fptype_v& b )
{
  return cxmake( a.real() + b, a.imag() );
}

inline
const cxtype_v& operator+( const cxtype_v& a )
{
  return a;
}

inline
cxtype_v operator-( const cxtype_v& a, const cxtype_v& b )
{
  return cxmake( a.real() - b.real(), a.imag() - b.imag() );
}

inline
cxtype_v operator-( const fptype& a, const cxtype_v& b )
{
  return cxmake( a - b.real(), - b.imag() );
}

inline
cxtype_v operator-( const cxtype_v& a )
{
  return 0 - a;
}

inline
cxtype_v operator-( const cxtype_v& a, const fptype& b )
{
  return cxmake( a.real() - b, a.imag() );
}

inline
cxtype_v operator-( const fptype_v& a, const cxtype_v& b )
{
  return cxmake( a - b.real(), - b.imag() );
}

inline
cxtype_v operator-( const cxtype_v& a, const fptype_v& b )
{
  return cxmake( a.real() - b, a.imag() );
}

inline
cxtype_v operator-( const fptype_v& a, const cxtype& b )
{
  return cxmake( a - b.real(), fptype_v{0} - b.imag() );
}

inline
cxtype_v operator*( const cxtype_v& a, const cxtype_v& b )
{
  return cxmake( a.real() * b.real() - a.imag() * b.imag(), a.imag() * b.real() + a.real() * b.imag() );
}

inline
cxtype_v operator*( const cxtype& a, const cxtype_v& b )
{
  return cxmake( a.real() * b.real() - a.imag() * b.imag(), a.imag() * b.real() + a.real() * b.imag() );
}

inline
cxtype_v operator*( const cxtype_v& a, const cxtype& b )
{
  return cxmake( a.real() * b.real() - a.imag() * b.imag(), a.imag() * b.real() + a.real() * b.imag() );
}

inline
cxtype_v operator*( const fptype& a, const cxtype_v& b )
{
  return cxmake( a * b.real(), a * b.imag() );
}

inline
cxtype_v operator*( const cxtype_v& a, const fptype& b )
{
  return cxmake( a.real() * b, a.imag() * b );
}

inline
cxtype_v operator*( const fptype_v& a, const cxtype_v& b )
{
  return cxmake( a * b.real(), a * b.imag() );
}

inline
cxtype_v operator*( const cxtype_v& a, const fptype_v& b )
{
  return cxmake( a.real() * b, a.imag() * b );
}

inline
cxtype_v operator*( const fptype_v& a, const cxtype& b )
{
  return cxmake( a * b.real(), a * b.imag() );
}

inline
cxtype_v operator*( const cxtype& a, const fptype_v& b )
{
  return cxmake( a.real() * b, a.imag() * b );
}

inline
cxtype_v operator/( const cxtype_v& a, const cxtype_v& b )
{
  fptype_v bnorm = b.real()*b.real() + b.imag()*b.imag();
  return cxmake( ( a.real() * b.real() + a.imag() * b.imag() ) / bnorm,
                 ( a.imag() * b.real() - a.real() * b.imag() ) / bnorm );
}

inline
cxtype_v operator/( const cxtype& a, const cxtype_v& b )
{
  fptype_v bnorm = b.real()*b.real() + b.imag()*b.imag();
  return cxmake( ( cxreal( a ) * b.real() + cximag( a ) * b.imag() ) / bnorm,
                 ( cximag( a ) * b.real() - cxreal( a ) * b.imag() ) / bnorm );
}

/*
inline
cxtype_v operator/( const fptype& a, const cxtype_v& b )
{
  fptype_v bnorm = b.real()*b.real() + b.imag()*b.imag();
  return cxmake( ( a * b.real() ) / bnorm, ( - a * b.imag() ) / bnorm );
}
*/

inline
cxtype_v operator/( const cxtype_v& a, const fptype_v& b )
{
  return cxmake( a.real() / b, a.imag() / b );
}

inline
cxtype_v operator/( const cxtype_v& a, const fptype& b )
{
  return cxmake( a.real() / b, a.imag() / b );
}
#endif

// Operators for bool_v
#ifdef MGONGPU_CPPSIMD

inline
fptype_v fpternary( const bool_v& mask, const fptype_v& a, const fptype_v& b )
{
  fptype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = ( mask[i] ? a[i] : b[i] );
  return out;
}

inline
fptype_v fpternary( const bool_v& mask, const fptype_v& a, const fptype& b )
{
  fptype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = ( mask[i] ? a[i] : b );
  return out;
}

inline
fptype_v fpternary( const bool_v& mask, const fptype& a, const fptype_v& b )
{
  fptype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = ( mask[i] ? a : b[i] );
  return out;
}

inline
fptype_v fpternary( const bool_v& mask, const fptype& a, const fptype& b )
{
  fptype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = ( mask[i] ? a : b );
  return out;
}

inline
cxtype_v cxternary( const bool_v& mask, const cxtype_v& a, const cxtype_v& b )
{
#ifdef MGONGPU_HAS_CXTYPE_REF
  cxtype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = ( mask[i] ? a[i] : b[i] );
  return out;
#else
  fptype_v outr, outi;
  for ( int i=0; i<neppV; i++ )
  {
    outr[i] = ( mask[i] ? a.real()[i] : b.real()[i] );
    outi[i] = ( mask[i] ? a.imag()[i] : b.imag()[i] );
  }
  return cxmake( outr, outi );
#endif
}

inline
cxtype_v cxternary( const bool_v& mask, const cxtype_v& a, const cxtype& b )
{
#ifdef MGONGPU_HAS_CXTYPE_REF
  cxtype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = ( mask[i] ? a[i] : b );
  return out;
#else
  fptype_v outr, outi;
  for ( int i=0; i<neppV; i++ )
  {
    outr[i] = ( mask[i] ? a.real()[i] : b.real() );
    outi[i] = ( mask[i] ? a.imag()[i] : b.imag() );
  }
  return cxmake( outr, outi );
#endif
}

inline
cxtype_v cxternary( const bool_v& mask, const cxtype& a, const cxtype_v& b )
{
#ifdef MGONGPU_HAS_CXTYPE_REF
  cxtype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = ( mask[i] ? a : b[i] );
  return out;
#else
  fptype_v outr, outi;
  for ( int i=0; i<neppV; i++ )
  {
    outr[i] = ( mask[i] ? a.real() : b.real()[i] );
    outi[i] = ( mask[i] ? a.imag() : b.imag()[i] );
  }
  return cxmake( outr, outi );
#endif
}

inline
cxtype_v cxternary( const bool_v& mask, const cxtype& a, const cxtype& b )
{
#ifdef MGONGPU_HAS_CXTYPE_REF
  cxtype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = ( mask[i] ? a : b );
  return out;
#else
  fptype_v outr, outi;
  for ( int i=0; i<neppV; i++ )
  {
    outr[i] = ( mask[i] ? a.real() : b.real() );
    outi[i] = ( mask[i] ? a.imag() : b.imag() );
  }
  return cxmake( outr, outi );
#endif
}

inline
fptype_v fpmax( const fptype_v& a, const fptype_v& b )
{
  return fpternary( ( b < a ), a, b );
}

inline
fptype_v fpmax( const fptype_v& a, const fptype& b )
{
  return fpternary( ( b < a ), a, b );
}

/*
inline
fptype_v fpmax( const fptype& a, const fptype_v& b )
{
  return fpternary( ( b < a ), a, b );
}
*/

inline
fptype_v fpmin( const fptype_v& a, const fptype_v& b )
{
  return fpternary( ( a < b ), a, b );
}

/*
inline
fptype_v fpmin( const fptype_v& a, const fptype& b )
{
  return fpternary( ( a < b ), a, b );
}

inline
fptype_v fpmin( const fptype& a, const fptype_v& b )
{
  return fpternary( ( a < b ), a, b );
}
*/

inline
bool maskor( const bool_v& mask )
{
  bool out = false;
  for ( int i=0; i<neppV; i++ ) out = out || mask[i];
  return out;
}

#else

inline
fptype fpternary( const bool& mask, const fptype& a, const fptype& b )
{
  return ( mask ? a : b );
}

inline
cxtype cxternary( const bool& mask, const cxtype& a, const cxtype& b )
{
  return ( mask ? a : b );
}

inline
bool maskor( const bool& mask )
{
  return mask;
}

#endif

//------------------------------
// Vector types - CUDA or ALPAKA
//------------------------------

#else

// Printout to std::cout for user defined types
#ifdef ALPAKA
inline ALPAKA_FN_ACC void print( const fptype& f ){ printf( "%f\n", f ); }
inline ALPAKA_FN_ACC void print( const cxtype& c ){ printf( "[%f, %f]\n", cxreal(c), cximag(c) ); }
#else
inline __device__ void print( const fptype& f ){ printf( "%f\n", f ); }
inline __device__ void print( const cxtype& c ){ printf( "[%f, %f]\n", cxreal(c), cximag(c) ); }
#endif

/*
inline __device__
const cxtype& cxvmake( const cxtype& c )
{
  return c;
}
*/

#ifdef ALPAKA
inline ALPAKA_FN_HOST_ACC
#else
inline __host__ __device__
#endif
fptype fpternary( const bool& mask, const fptype& a, const fptype& b )
{
  return ( mask ? a : b );
}

#ifdef ALPAKA
inline ALPAKA_FN_HOST_ACC
#else
inline __host__ __device__
#endif
cxtype cxternary( const bool& mask, const cxtype& a, const cxtype& b )
{
  return ( mask ? a : b );
}

#endif

//--------------------------------------------------------------------------

// Scalar-or-vector types: scalar in CUDA or ALPAKA, vector or scalar in C++
#if defined(ALPAKA) || defined(__CUDACC__)
typedef bool bool_sv;
typedef fptype fptype_sv;
typedef cxtype cxtype_sv;
#elif defined MGONGPU_CPPSIMD
typedef bool_v bool_sv;
typedef fptype_v fptype_sv;
typedef cxtype_v cxtype_sv;
#else
typedef bool bool_sv;
typedef fptype fptype_sv;
typedef cxtype cxtype_sv;
#endif

// Scalar-or-vector zeros: scalar in CUDA or ALPAKA, vector or scalar in C++
#ifdef ALPAKA
inline ALPAKA_FN_ACC cxtype cxzero_sv(){ return cxmake( 0, 0 ); }
#elif defined(__CUDACC__)
inline __device__ cxtype cxzero_sv(){ return cxmake( 0, 0 ); }
#elif defined MGONGPU_CPPSIMD
inline cxtype_v cxzero_sv(){ return cxtype_v{ fptype_v{0}, fptype_v{0} }; }
#else
inline cxtype cxzero_sv(){ return cxtype{ fptype{0}, fptype{0} }; }
#endif

//--------------------------------------------------------------------------

#endif // MGONGPUVECTORS_H
