#ifndef MGONGPUVECTORS_H
#define MGONGPUVECTORS_H 1

//------------------------------
// Vector types - C++
//------------------------------

#ifndef __CUDACC__

#include "mgOnGpuTypes.h"

namespace mgOnGpu
{
  const int neppV = neppM;

  // --- Type definition (using vector compiler extensions: need -march=...)
  typedef fptype fptype_v alignas(neppV*sizeof(fptype)) __attribute__ ((vector_size (neppV*sizeof(fptype)))); // RRRR

  class cxtype_ref
  {
  public:
    cxtype_ref() = delete;
    cxtype_ref( const cxtype_ref& ) = delete;
    cxtype_ref( cxtype_ref&& ) = default;
    cxtype_ref( fptype& r, fptype& i ) : m_real{r}, m_imag{i} {}
    cxtype_ref& operator=( const cxtype_ref& ) = delete;
    cxtype_ref& operator=( cxtype_ref&& ) = delete;
    cxtype_ref& operator=( const cxtype& c ) { m_real = cxreal( c ); m_imag = cximag( c ); return *this; }
    operator cxtype() const { return cxmake( m_real, m_imag ); }
  private:
    fptype &m_real, &m_imag; // RI
  };

  // --- Type definition (using vector compiler extensions: need -march=...)
  class alignas(2*sizeof(fptype_v)) cxtype_v
  {
  public:
    cxtype_v() : m_real{0}, m_imag{0} {}
    cxtype_v( const cxtype_v&  ) = default;
    cxtype_v( cxtype_v&&  ) = default;
    cxtype_v( const fptype_v& r, const fptype_v& i ) : m_real{r}, m_imag{i} {}
    cxtype_v& operator=( const cxtype_v& ) = default;
    cxtype_v& operator=( cxtype_v&& ) = default;
    cxtype_ref operator[]( size_t i ) const { return cxtype_ref( m_real[i], m_imag[i] ); }
    const fptype_v& real() const { return m_real; }
    const fptype_v& imag() const { return m_imag; }
  private:
    fptype_v m_real, m_imag; // RRRRIIII
  };

}

// Expose typedefs outside the namespace
using mgOnGpu::neppV;
using mgOnGpu::fptype_v;
using mgOnGpu::cxtype_v;

// Printout to stream for user defined types
inline std::ostream& operator<<( std::ostream& out, const fptype_v& v )
{
  out << "{ " << v[0];
  for ( int i=1; i<neppV; i++ ) std::cout << ", " << v[i];
  out << " }";
  return out;
}

inline std::ostream& operator<<( std::ostream& out, const cxtype& c )
{
  out << "[" << cxreal(c) << "," << cximag(c) << "]";
  //out << cxreal(c) << "+i" << cximag(c);
  return out;
}

inline std::ostream& operator<<( std::ostream& out, const cxtype_v& v )
{
  out << "{ " << v[0];
  for ( int i=1; i<neppV; i++ ) std::cout << ", " << v[i];
  out << " }";
  return out;
}

// Printout to std::cout for user defined types
inline void print( const fptype& f ) { std::cout << f << std::endl; }
inline void print( const fptype_v& v ) { std::cout << v << std::endl; }
inline void print( const cxtype& c ) { std::cout << c << std::endl; }
inline void print( const cxtype_v& v ) { std::cout << v << std::endl; }

// Operators for fptype_v
inline
fptype_v sqrt( const fptype_v& v )
{
  fptype_v out;
  for ( int i=0; i<neppV; i++ ) out[i]=sqrt(v[i]);
  return out;
}

inline
fptype_v fpmake_v( const fptype v[neppV] )
{
  fptype_v out;
  for ( int i=0; i<neppV; i++ ) out[i] = v[i];
  return out;
}

// Operators for cxtype_v
/*
inline
cxtype_v cxmake_v( const cxtype c )
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
cxtype_v cxmake00()
{
  return cxtype_v{ fptype_v{0}, fptype_v{0} };
}

inline
cxtype_v cxmaker0( const fptype_v& r )
{
  return cxtype_v{ r, fptype_v{0} };
}

/*
inline
cxtype_v cxmake0i( const fptype_v& i )
{
  return cxtype_v{ fptype_v{0}, i };
}
*/

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
cxtype_v operator/( const cxtype_v& a, const fptype& b )
{
  return cxmake( a.real() / b, a.imag() / b );
}

//------------------------------
// Vector types - CUDA
//------------------------------

#else

// Printout to std::cout for user defined types
inline __device__ void print( const fptype& f ){ printf( "%f\n", f ); }
inline __device__ void print( const cxtype& c ){ printf( "[%f, %f]\n", cxreal(c), cximag(c) ); }

inline __device__
const cxtype& cxvmake( const cxtype& c )
{
  return c;
}

inline __device__
cxtype cxmake00()
{
  return cxmake( 0, 0 );
}

inline __device__
cxtype cxmaker0( const fptype& r )
{
  return cxtype( r, 0 );
}

/*
inline
cxtype cxmake0i( const fptype& i )
{
  return cxmake( 0, i );
}
*/

#endif

#endif // MGONGPUVECTORS_H
