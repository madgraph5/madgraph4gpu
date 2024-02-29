#ifndef MGONGPUCXTYPES_H
#define MGONGPUCXTYPES_H 1

#include "mgOnGpuConfig.h"

#include <iostream>

namespace mgOnGpu /* clang-format off */
{
  // --- Type definition (simple complex type derived from cxtype_v)
  template<typename FP>
  class cxsmpl
  {
  public:
    constexpr cxsmpl( const FP& r = FP(), const FP& i = FP() ) : m_real( r ), m_imag( i ) {}
    cxsmpl( const cxsmpl& ) = default;
    cxsmpl( cxsmpl&& ) = default;
    cxsmpl& operator=( const cxsmpl& ) = default;
    cxsmpl& operator=( cxsmpl&& ) = default;
    constexpr cxsmpl& operator+=( const cxsmpl& c ) { m_real += c.real(); m_imag += c.imag(); return *this; }
    constexpr cxsmpl& operator-=( const cxsmpl& c ) { m_real -= c.real(); m_imag -= c.imag(); return *this; }
    constexpr const FP& real() const { return m_real; }
    constexpr const FP& imag() const { return m_imag; }
  private:
    FP m_real, m_imag; // RI
  };

  template<typename FP>
  inline cxsmpl<FP> // (NB: cannot be constexpr as a constexpr function cannot have a nonliteral return type "mgOnGpu::cxsmpl")
  conj( const cxsmpl<FP>& c )
  {
    return cxsmpl<FP>( c.real(), -c.imag() );
  }
} /* clang-format on */

// Expose the cxsmpl class outside the namespace
using mgOnGpu::cxsmpl;

// Printout to stream for user defined types
template<typename FP>
inline std::ostream&
operator<<( std::ostream& out, const cxsmpl<FP>& c )
{
  //out << std::complex( c.real(), c.imag() );
  out << "[" << c.real() << "," << c.imag() << "]";
  return out;
}

// Operators for cxsmpl
template<typename FP>
inline constexpr cxsmpl<FP>
operator+( const cxsmpl<FP> a )
{
  return a;
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator-( const cxsmpl<FP>& a )
{
  return cxsmpl<FP>( -a.real(), -a.imag() );
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator+( const cxsmpl<FP>& a, const cxsmpl<FP>& b )
{
  return cxsmpl<FP>( a.real() + b.real(), a.imag() + b.imag() );
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator+( const FP& a, const cxsmpl<FP>& b )
{
  return cxsmpl<FP>( a, (FP)0 ) + b;
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator-( const cxsmpl<FP>& a, const cxsmpl<FP>& b )
{
  return cxsmpl<FP>( a.real() - b.real(), a.imag() - b.imag() );
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator-( const FP& a, const cxsmpl<FP>& b )
{
  return cxsmpl<FP>( a, (FP)0 ) - b;
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator*( const cxsmpl<FP>& a, const cxsmpl<FP>& b )
{
  return cxsmpl<FP>( a.real() * b.real() - a.imag() * b.imag(), a.imag() * b.real() + a.real() * b.imag() );
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator*( const FP& a, const cxsmpl<FP>& b )
{
  return cxsmpl<FP>( a, (FP)0 ) * b;
}

inline constexpr cxsmpl<float>
operator*( const double& a, const cxsmpl<float>& b )
{
  return cxsmpl<float>( a, (float)0 ) * b;
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator/( const cxsmpl<FP>& a, const cxsmpl<FP>& b )
{
  FP bnorm = b.real() * b.real() + b.imag() * b.imag();
  return cxsmpl<FP>( ( a.real() * b.real() + a.imag() * b.imag() ) / bnorm,
                     ( a.imag() * b.real() - a.real() * b.imag() ) / bnorm );
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator/( const FP& a, const cxsmpl<FP>& b )
{
  return cxsmpl<FP>( a, (FP)0 ) / b;
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator+( const cxsmpl<FP>& a, const FP& b )
{
  return a + cxsmpl<FP>( b, (FP)0 );
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator-( const cxsmpl<FP>& a, const FP& b )
{
  return a - cxsmpl<FP>( b, (FP)0 );
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator*( const cxsmpl<FP>& a, const FP& b )
{
  return a * cxsmpl<FP>( b, (FP)0 );
}

template<typename FP>
inline constexpr cxsmpl<FP>
operator/( const cxsmpl<FP>& a, const FP& b )
{
  return a / cxsmpl<FP>( b, (FP)0 );
}

#endif
