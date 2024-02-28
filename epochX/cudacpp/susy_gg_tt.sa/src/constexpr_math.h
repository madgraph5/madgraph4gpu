// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Feb 2024) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

#ifndef constexpr_math_h
#define constexpr_math_h 1

#include "mgOnGpuConfig.h"

#include <cassert>
#include <cmath>
#include <limits>

// FOR DEBUGGING!
#undef CONSTEXPR_MATH_DEBUG // no-debug
//#define CONSTEXPR_MATH_DEBUG 1 // debug
#ifdef CONSTEXPR_MATH_DEBUG
#define constexpr const
#endif

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{  
  // Constexpr implementation of sqrt (see https://stackoverflow.com/a/34134071)
  constexpr long double sqrtNewtonRaphson( const long double x, const long double curr, const long double prev )
  {
    return curr == prev ? curr : sqrtNewtonRaphson( x, 0.5 * ( curr + x / curr ), curr );
  }
  constexpr long double constexpr_sqrt( const long double x )
  {
    return x >= 0 // && x < std::numeric_limits<double>::infinity() // avoid -Wtautological-constant-compare warning in fast math
      ? sqrtNewtonRaphson( x, x, 0 )
      : std::numeric_limits<long double>::quiet_NaN();
  }

  // Constexpr implementation of floor (see https://stackoverflow.com/a/66146159)
  constexpr int constexpr_floor( const long double d )
  {
    const int i = static_cast<int>( d );
    return d < i ? i - 1 : i;
  }

  // Constexpr implementation of pow
  constexpr long double constexpr_pow( const long double base, const long double exp )
  {
    // NB(1): this implementation of constexpr_pow requires exponent >= 0
    assert( exp >= 0 ); // NB would fail at compile time with "error: call to non-‘constexpr’ function ‘void __assert_fail'"
    // NB(2): this implementation of constexpr_pow requires an integer exponent
    const int iexp = constexpr_floor( exp );
    assert( static_cast<long double>( iexp ) == exp ); // NB would fail at compile time with "error: call to non-‘constexpr’ function ‘void __assert_fail'"
    // Iterative implementation of pow if exp is a non negative integer
    return iexp == 0 ? 1 : base * constexpr_pow( base, iexp - 1 );
  }

  // PI from cmath
  constexpr long double constexpr_pi = M_PIl; // pi
  constexpr long double constexpr_pi_by_2 = M_PI_2l; // pi/2
  constexpr long double constexpr_pi_by_4 = M_PI_4l; // pi/4

  // Constexpr implementation of sin for 0<x<pi/4 (long double signature)
  // Taylor expansion : x - x**3/3! + x**5/5!
  constexpr long double sinTaylor( const long double x )
  {
    assert( x >= 0 && "The argument of sinTaylor is assumed to be in [0,pi/4)" );
    assert( x < constexpr_pi_by_4 && "The argument of sinTaylor is assumed to be in [0,pi/4)" );
    long double sinx = 0;
    int ipow = 1;
    long double delta = x;
    while( true )
    {
      long double sinxlast = sinx;
      sinx += delta;
#ifdef CONSTEXPR_MATH_DEBUG
      std::cout << "ipow=" << ipow << ", delta=" << delta << ", sinx=" << sinx <<std::endl; // for debugging (not constexpr)
#endif
      if ( sinx == sinxlast ) break;
      // Next iteration
      ipow += 2;
      delta *= -x * x / ( ipow - 1 ) / ipow;
    }
    return sinx;
  }

  // Constexpr implementation of cos (long double signature)
  constexpr long double constexpr_cos_quad( const long double x, const bool assume0to2Pi = false )
  {
    if( assume0to2Pi )
    {
      assert( x >= 0 && "The argument of constexpr_cos_quad is assumed to be in [0,2*pi)" );
      assert( x < 2 * constexpr_pi && "The argument of constexpr_cos_quad is assumed to be in [0,2*pi)" );
    }
    if( x < 0 )
      return constexpr_cos_quad( x + ( constexpr_floor( -x / constexpr_pi ) + 1 ) * constexpr_pi, true );
    else if( x < constexpr_pi_by_4 ) // [0/4*pi, 1/4*pi)
      return constexpr_sqrt( 1 - constexpr_pow( sinTaylor( (long double)x ), 2 ) );
    else if( x < constexpr_pi_by_2 ) // [1/4*pi, 2/4*pi)
      return sinTaylor( constexpr_pi_by_2 - (long double)x );
    else if( x < 3 * constexpr_pi_by_4 ) // [2/4*pi, 3/4*pi)
      return -sinTaylor( (long double)x - constexpr_pi_by_2 );
    else if( x < constexpr_pi ) // [3/4*pi, 4/4*pi)
      return -constexpr_sqrt( 1 - constexpr_pow( sinTaylor( constexpr_pi - (long double)x ), 2 ) );
    else if( x < 2 * constexpr_pi ) // [4/4*pi, 8/4*pi)
      return constexpr_cos_quad( 2 * constexpr_pi - (long double)x, true );
    else // [8/4*pi, +inf)
      return constexpr_cos_quad( x - constexpr_floor( x / constexpr_pi ) * constexpr_pi, true );
  }

  // Constexpr implementation of cos (double signature, internally implemented as long double)
  constexpr double constexpr_cos( const double x )
  {
    return constexpr_cos_quad( x );
  }

  // Constexpr implementation of sin (long double signature)
  constexpr long double constexpr_sin_quad( const long double x, const bool assume0toPi = false )
  {
    if( assume0toPi )
    {
      assert( x >= 0 && "The argument of constexpr_sin_quad is assumed to be in [0,pi)" );
      assert( x < constexpr_pi && "The argument of constexpr_sin_quad is assumed to be in [0,pi)" );
    }
    if( x < 0 )
      return constexpr_sin_quad( x + ( constexpr_floor( -x / constexpr_pi ) + 1 ) * constexpr_pi, true );
    else if( x < constexpr_pi_by_4 ) // [0/4*pi, 1/4*pi)
      return sinTaylor( (long double)x );
    else if( x < constexpr_pi_by_2 ) // [1/4*pi, 2/4*pi)
      return constexpr_sqrt( 1 - constexpr_pow( sinTaylor( constexpr_pi_by_2 - (long double)x ), 2 ) );
    else if( x < 3 * constexpr_pi_by_4 ) // [2/4*pi, 3/4*pi)
      return constexpr_sqrt( 1 - constexpr_pow( sinTaylor( (long double)x - constexpr_pi_by_2 ), 2 ) );
    else if( x < constexpr_pi ) // [3/4*pi, 4/4*pi)
      return sinTaylor( constexpr_pi - (long double)x );
    else if( x < 2 * constexpr_pi ) // [4/4*pi, 8/4*pi)
      return -constexpr_sin_quad( 2 * constexpr_pi - (long double)x, true );
    else // [8/4*pi, +inf)
      return constexpr_sin_quad( x - constexpr_floor( x / constexpr_pi ) * constexpr_pi, true );
  }

  // Constexpr implementation of sin (double signature, internally implemented as long double)
  constexpr double constexpr_sin( const double x )
  {
    return constexpr_sin_quad( x );
  }

  // Constexpr implementation of tan (long double signature)
  constexpr long double constexpr_tan_quad( const long double x, const bool assume0to2Pi = false )
  {
    if( assume0to2Pi )
    {
      assert( x >= 0 && "The argument of constexpr_sin_quad is assumed to be in [0,2*pi)" );
      assert( x < 2 * constexpr_pi && "The argument of constexpr_sin_quad is assumed to be in [0,2*pi)" );
    }
    if( x < 0 )
      return constexpr_tan_quad( x + ( constexpr_floor( -x / constexpr_pi ) + 1 ) * constexpr_pi, true );
    else if( x < 2 * constexpr_pi ) // [0, 2*pi)
      return constexpr_sin_quad( x, assume0to2Pi ) / constexpr_cos_quad( x, assume0to2Pi );
    else // [8/4*pi, +inf)
      return constexpr_tan_quad( x - constexpr_floor( x / constexpr_pi ) * constexpr_pi, true );
  }

  // Constexpr implementation of tan (double signature, internally implemented as long double)
  constexpr double constexpr_tan( const double x )
  {
    return constexpr_tan_quad( x );
  }
  
  // Constexpr implementation of atan for -1<x<1 (long double signature)
  // Taylor expansion : x - x**3/3 + x**5/5...
  constexpr long double atanTaylor( const long double x )
  {
    assert( x >= -1 && "The argument of atanTaylor is assumed to be in (-1,+1)" );
    assert( x < 1 && "The argument of atanTaylor is assumed to be in (-1,+1)" );
    long double atanx = 0;
    int ipow = 1;
    long double xpow = x;
    while( true )
    {
      long double atanxlast = atanx;
      atanx += xpow / ipow;
#ifdef CONSTEXPR_MATH_DEBUG
      std::cout << "ipow=" << ipow << ", xpow=" << xpow << ", atanx=" << atanx <<std::endl; // for debugging (not constexpr)
#endif
      if ( atanx == atanxlast ) break;
      // Next iteration
      ipow += 2;
      xpow *= -x * x;
    }
    return atanx;
  }

  // Constexpr implementation of atan (long double signature)
  constexpr long double constexpr_atan_quad( const long double x )
  {
    if( x > 1 )
      return atanTaylor( 1 / (long double)x );
    else if( x == 1 )
      return constexpr_pi_by_4;
    else if( x > -1 )
      return atanTaylor( (long double)x );
    else if( x == -1 )
      return -constexpr_pi_by_4;
    else if( x < -1 )
      return atanTaylor( 1 / (long double)x );
  }

  // Constexpr implementation of atan (double signature, internally implemented as long double)
  constexpr double constexpr_atan( const double x )
  {
    return constexpr_atan_quad( x );
  }
}

#endif // constexpr_math_h
