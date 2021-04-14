#ifndef MGONGPUVECTORS_H
#define MGONGPUVECTORS_H 1

#include <iostream>

//------------------------------
// Vector types - C++
//------------------------------

#ifndef __CUDACC__

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
    //cxtype_v& operator+=( const cxtype_v& c ){ m_real += c.real(); m_imag += c.imag(); return *this; }
    cxtype_v& operator-=( const cxtype_v& c ){ m_real -= c.real(); m_imag -= c.imag(); return *this; }

#ifdef __clang__
    // ERROR! clang build fails: [] is a value not a ref?
    //cxtype_ref operator[]( size_t i ) const { return cxtype_ref( m_real[i], m_imag[i] ); }
    // ERROR! clang build crashes (probably because [] is a value not a ref)
    cxtype_ref operator[]( size_t i ) const { return cxtype_ref( (fptype&)(m_real[i]), (fptype&)(m_imag[i]) ); }
#else
    cxtype_ref operator[]( size_t i ) const { return cxtype_ref( m_real[i], m_imag[i] ); }
#endif

    const fptype_v& real() const { return m_real; }
    const fptype_v& imag() const { return m_imag; }
  private:
    fptype_v m_real, m_imag; // RRRRIIII
  };

  // --- Type definition (using vector compiler extensions: need -march=...)
#ifdef __clang__ // https://clang.llvm.org/docs/LanguageExtensions.html#vectors-and-extended-vectors
  typedef long int bool_v __attribute__ ((ext_vector_type(neppV))); // bbbb
#else
  typedef long int bool_v __attribute__ ((vector_size (neppV*sizeof(long int)))); // bbbb
#endif

#else

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
  out << "{ " << v[0];
  for ( int i=1; i<neppV; i++ ) std::cout << ", " << v[i];
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
fptype_v sqrt( const fptype_v& v )
{
  fptype_v out;
  for ( int i=0; i<neppV; i++ ) out[i]=sqrt(v[i]);
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
  cxtype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = ( mask[i] ? a[i] : b[i] );
  return out;
}

inline
cxtype_v cxternary( const bool_v& mask, const cxtype_v& a, const cxtype& b )
{
  cxtype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = ( mask[i] ? a[i] : b );
  return out;
}

inline
cxtype_v cxternary( const bool_v& mask, const cxtype& a, const cxtype_v& b )
{
  cxtype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = ( mask[i] ? a : b[i] );
  return out;
}

inline
cxtype_v cxternary( const bool_v& mask, const cxtype& a, const cxtype& b )
{
  cxtype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = ( mask[i] ? a : b );
  return out;
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

#endif

//------------------------------
// Vector types - CUDA
//------------------------------

#else

// Printout to std::cout for user defined types
inline __device__ void print( const fptype& f ){ printf( "%f\n", f ); }
inline __device__ void print( const cxtype& c ){ printf( "[%f, %f]\n", cxreal(c), cximag(c) ); }

/*
inline __device__
const cxtype& cxvmake( const cxtype& c )
{
  return c;
}
*/

inline __host__ __device__
fptype fpternary( const bool& mask, const fptype& a, const fptype& b )
{
  return ( mask ? a : b );
}

inline __host__ __device__
cxtype cxternary( const bool& mask, const cxtype& a, const cxtype& b )
{
  return ( mask ? a : b );
}

#endif

//--------------------------------------------------------------------------

// Scalar-or-vector types: scalar in CUDA, vector or scalar in C++
#ifdef __CUDACC__
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

// Scalar-or-vector zeros: scalar in CUDA, vector or scalar in C++
#ifdef __CUDACC__
inline __device__ cxtype cxzero_sv(){ return cxmake( 0, 0 ); }
#elif defined MGONGPU_CPPSIMD
inline cxtype_v cxzero_sv(){ return cxtype_v{ fptype_v{0}, fptype_v{0} }; }
#else
inline cxtype cxzero_sv(){ return cxtype{ fptype{0}, fptype{0} }; }
#endif

//--------------------------------------------------------------------------

#endif // MGONGPUVECTORS_H
