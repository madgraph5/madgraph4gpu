#ifndef MGONGPUVECTORS_H
#define MGONGPUVECTORS_H 1

#ifndef __CUDACC__

#include "mgOnGpuTypes.h"

namespace mgOnGpu
{
  const int neppV = neppM;

  // --- Type definition (using vector compiler extensions: need -march=...)
  typedef fptype fptype_v __attribute__ ((vector_size (neppV * sizeof(fptype)))); // RRRR

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
  class cxtype_v
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

// DEBUG - START
void print( const fptype& f )
{
  std::cout << f << std::endl;
}

void print( const fptype_v& v )
{
  std::cout << "{ " << v[0];
  for ( int i=1; i<neppV; i++ ) std::cout << ", " << v[i];
  std::cout << " }" << std::endl;
}

std::string str( const cxtype& c )
{
  std::stringstream ss;
  ss << "[" << cxreal(c) << "," << cximag(c) << "]";
  return ss.str();
}

void print( const cxtype& c )
{
  std::cout << str(c) << std::endl;
}

void print( const cxtype_v& v )
{
  std::cout << "{ " << str( v[0] );
  for ( int i=1; i<neppV; i++ ) std::cout << ", " << str( v[i] );
  std::cout << " }" << std::endl;
}
// DEBUG - END

// Operators for fptype_v
inline
fptype_v sqrt( const fptype_v& v )
{
  fptype_v out;
  for ( int i=0; i<neppV; i++ ) out[i]=sqrt(v[i]);
  return out;
}

inline
fptype_v fpvmake( const fptype v[neppV] )
{
  fptype_v out;
  for ( int i=0; i<neppV; i++ ) out[i]=v[i];
  return out;
}

// Operators for cxtype_v
inline
cxtype_v cxvmake( const cxtype c )
{
  print( c );
  cxtype_v out;
  for ( int i=0; i<neppV; i++ )
  {
    std::cout << "START" << std::endl;
    print( c );
    out[i] = c;
    print( out[i] );
    std::cout << "END" << std::endl;
  }
  print( out );
  return out;
}

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
cxtype_v& operator-( const cxtype_v& a )
{
  return -a;
}

inline
cxtype_v operator-( const fptype& a, const cxtype_v& b )
{
  return cxmake( a - b.real(), - b.imag() );
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
cxtype_v operator*( const cxtype_v& a, const cxtype_v& b )
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

#else

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
