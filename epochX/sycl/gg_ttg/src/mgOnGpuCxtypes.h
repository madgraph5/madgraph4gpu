#ifndef MGONGPUCXTYPES_H
#define MGONGPUCXTYPES_H 1

#include "extras.h"
#include "mgOnGpuConfig.h"
#include "mgOnGpuFptypes.h"

//==========================================================================
// COMPLEX TYPES: (PLATFORM-SPECIFIC) HEADERS
//==========================================================================

#include <complex>
#include <cmath>

// Complex type in SYCL

//==========================================================================
// COMPLEX TYPES: (PLATFORM-SPECIFIC) TYPEDEFS
//==========================================================================

namespace mgOnGpu
{

  // --- Type definitions (complex type: cxtype)
  typedef extras::complex<fptype> cxtype;

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

//------------------------------
// SYCL - extras::complex
//------------------------------
SYCL_EXTERNAL inline
cxtype cxmake( const fptype& r, const fptype& i )
{
  return cxtype( r, i ); // thrust::complex<fptype> constructor
}

SYCL_EXTERNAL inline
fptype cxreal( const cxtype& c )
{
  return c.real(); // thrust::complex<fptype>::real()
}

SYCL_EXTERNAL inline
fptype cximag( const cxtype& c )
{
  return c.imag(); // thrust::complex<fptype>::imag()
}

SYCL_EXTERNAL inline
cxtype cxconj( const cxtype& c )
{
  return extras::conj( c ); // extras::conj( extras::complex<fptype> )
}

SYCL_EXTERNAL inline
const cxtype& cxmake( const cxtype& c )
{
  return c;
}

inline // NOT __device__
cxtype cxmake( const std::complex<fptype>& c ) // std::complex to extras::complex (float-to-float or double-to-double)
{
  return cxmake( c.real(), c.imag() );
}

#if defined MGONGPU_FPTYPE_FLOAT
inline
cxtype cxmake( const std::complex<double>& c ) // std::complex to std::complex (cast double-to-float)
{
  return cxmake( (fptype)c.real(), (fptype)c.imag() );
}
#endif

//==========================================================================

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
    cxtype_ref( fptype& r, fptype& i ) : m_real( r ), m_imag( i ) {}
    cxtype_ref& operator=( const cxtype_ref& ) = delete;
    cxtype_ref& operator=( cxtype_ref&& c ) { m_real = cxreal( c ); m_imag = cximag( c ); return *this; } // for cxternary
    cxtype_ref& operator=( const cxtype& c ) { m_real = cxreal( c ); m_imag = cximag( c ); return *this; }
    SYCL_EXTERNAL operator cxtype() const { return cxmake( m_real, m_imag ); }
  private:
    fptype &m_real, &m_imag; // RI
  };
}

//// Printout to stream for user defined types
//SYCL_EXTERNAL inline
//std::ostream& operator<<( std::ostream& out, const cxtype& c ) {
//    out << "[" << cxreal(c) << "," << cximag(c) << "]";
//    return out;
//}
//
//SYCL_EXTERNAL inline
//std::ostream& operator<<( std::ostream& out, const mgOnGpu::cxtype_ref& c ){ out << (cxtype)c; return out; }

//==========================================================================

#endif // MGONGPUCXTYPES_H
