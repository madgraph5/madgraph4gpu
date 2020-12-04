#ifndef MGONGPUVECTORS_H
#define MGONGPUVECTORS_H 1

#ifndef __CUDACC__

#include "mgOnGpuTypes.h"

namespace mgOnGpu
{
  const int neppV = neppM;

  // --- Type definitions
  //typedef fptype fptype_v[neppV]; // RRRR
  //typedef cxtype cxtype_v[neppV]; // RIRIRIRI (not RRRRIIII)

  // --- Type definitions (using vector compiler extensions: need -march=native)
  typedef fptype fptype_v __attribute__ ((vector_size (neppV * sizeof(fptype)))); // RRRR
  struct cxtype_v
  {
    fptype_v real, imag; // RRRRIIII
    cxtype operator[]( size_t i ) const { return cxmake( real[i], imag[i] ); }
  };

}

// Expose typedefs outside the namespace
using mgOnGpu::neppV;
using mgOnGpu::fptype_v;
using mgOnGpu::cxtype_v;

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
  for ( int i=1; i<neppV; i++ ) std::cout << ", " << str( v[i] ) << " }";
  std::cout << std::endl;
}
// DEBUG - END

#else

inline __device__
cxtype cxmake00()
{
  return cxmake( 0, 0 );
}

/*
inline
cxtype cxmaker0( const fptype& r )
{
  return cxtype( r, 0 );
}

inline
cxtype cxmake0i( const fptype& i )
{
  return cxmake( 0, i );
}
*/

#endif

#endif // MGONGPUVECTORS_H
