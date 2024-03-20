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
  constexpr long double sqrtNewtonRaphson( const long double xx, const long double curr, const long double prev )
  {
    return curr == prev ? curr : sqrtNewtonRaphson( xx, 0.5 * ( curr + xx / curr ), curr );
  }
  constexpr long double constexpr_sqrt( const long double xx )
  {
    return xx >= 0 // && x < std::numeric_limits<double>::infinity() // avoid -Wtautological-constant-compare warning in fast math
      ? sqrtNewtonRaphson( xx, xx, 0 )
      : std::numeric_limits<long double>::quiet_NaN();
  }

  // Constexpr implementation of floor (see https://stackoverflow.com/a/66146159)
  constexpr int constexpr_floor( const long double xx )
  {
    const int i = static_cast<int>( xx );
    return xx < i ? i - 1 : i;
  }

  // Constexpr implementation of pow
  constexpr long double constexpr_pow( const long double base, const long double exp, const bool requireExpGe0 = false )
  {
    // NB(1): this iterative integer implementation of constexpr_pow requires exponent >= 0
    if ( requireExpGe0 ) assert( exp >= 0 ); // NB would fail at compile time with "error: call to non-‘constexpr’ function ‘void __assert_fail'"
    if ( exp < 0 ) return 1. / constexpr_pow( base, -exp, true );
    // NB(2): this iterative integer implementation of constexpr_pow requires an integer exponent
    const int iexp = constexpr_floor( exp );
    assert( static_cast<long double>( iexp ) == exp ); // NB would fail at compile time with "error: call to non-‘constexpr’ function ‘void __assert_fail'"
    // Iterative implementation of pow if exp is a non negative integer
    return iexp == 0 ? 1 : base * constexpr_pow( base, iexp - 1 );
  }

  // PI constants
  // NB1: M_PIl from from cmath is not defined on Mac
  // NB2: std::numbers::pi needs c++20 but we are still using c++17
  // NB3: I could use my constexpr_atan(1)*4... but a literal is better?
  //constexpr long double constexpr_pi = M_PIl;        // pi
  //constexpr long double constexpr_pi_by_2 = M_PI_2l; // pi/2
  //constexpr long double constexpr_pi_by_4 = M_PI_4l; // pi/4
  constexpr long double constexpr_pi = 3.141592653589793238462643383279502884L;      // same as M_PIl in gcc
  constexpr long double constexpr_pi_by_2 = 1.570796326794896619231321691639751442L; // same as M_PI_2l in gcc
  constexpr long double constexpr_pi_by_4 = 0.785398163397448309615660845819875721L; // same as M_PI_4l in gcc
  static_assert( constexpr_pi_by_4 * 4 == constexpr_pi );
  static_assert( constexpr_pi_by_4 * 2 == constexpr_pi_by_2 );
  static_assert( constexpr_pi_by_2 * 2 == constexpr_pi );

  // Constexpr implementation of sin for 0<x<pi/4 (long double signature)
  // Taylor expansion : x - x**3/3! + x**5/5!
  constexpr long double sinTaylor( const long double xx )
  {
    assert( xx >= 0 && "The argument of sinTaylor is assumed to be in [0,pi/4)" );
    assert( xx < constexpr_pi_by_4 && "The argument of sinTaylor is assumed to be in [0,pi/4)" );
    long double sinx = 0;
    int ipow = 1;
    long double delta = xx;
    while( true )
    {
      long double sinxlast = sinx;
      sinx += delta;
#ifdef CONSTEXPR_MATH_DEBUG
      std::cout << "ipow=" << ipow << ", delta=" << delta << ", sinx=" << sinx << std::endl; // for debugging (not constexpr)
#endif
      if( sinx == sinxlast ) break;
      // Next iteration
      ipow += 2;
      delta *= -xx * xx / ( ipow - 1 ) / ipow;
    }
    return sinx;
  }

  // Mapping to [0,2*pi) range (long double signature)
  constexpr long double mapIn0to2Pi( const long double xx )
  {
    return xx - constexpr_floor( xx / 2 / constexpr_pi ) * 2 * constexpr_pi;
  }

  // Constexpr implementation of cos (long double signature)
  constexpr long double constexpr_cos_quad( const long double xx, const bool assume0to2Pi = false )
  {
    if( assume0to2Pi )
    {
      assert( xx >= 0 && "The argument of constexpr_cos_quad is assumed to be in [0,2*pi)" );
      assert( xx < 2 * constexpr_pi && "The argument of constexpr_cos_quad is assumed to be in [0,2*pi)" );
    }
    if( xx < 0 )
      return constexpr_cos_quad( mapIn0to2Pi( xx ), true );
    else if( xx < constexpr_pi_by_4 ) // [0/4*pi, 1/4*pi)
      return constexpr_sqrt( 1 - constexpr_pow( sinTaylor( xx ), 2 ) );
    else if( xx < constexpr_pi_by_2 ) // [1/4*pi, 2/4*pi)
      return sinTaylor( constexpr_pi_by_2 - xx );
    else if( xx < 3 * constexpr_pi_by_4 ) // [2/4*pi, 3/4*pi)
      return -sinTaylor( xx - constexpr_pi_by_2 );
    else if( xx < constexpr_pi ) // [3/4*pi, 4/4*pi)
      return -constexpr_sqrt( 1 - constexpr_pow( sinTaylor( constexpr_pi - xx ), 2 ) );
    else if( xx < 2 * constexpr_pi ) // [4/4*pi, 8/4*pi)
      return constexpr_cos_quad( 2 * constexpr_pi - xx, true );
    else // [8/4*pi, +inf)
      return constexpr_cos_quad( mapIn0to2Pi( xx ), true );
  }

  // Constexpr implementation of cos (double signature, internally implemented as long double)
  constexpr double constexpr_cos( const double x )
  {
    return constexpr_cos_quad( x );
  }

  // Constexpr implementation of sin (long double signature)
  constexpr long double constexpr_sin_quad( const long double xx, const bool assume0to2Pi = false )
  {
    if( assume0to2Pi )
    {
      assert( xx >= 0 && "The argument of constexpr_sin_quad is assumed to be in [0,2*pi)" );
      assert( xx < 2 * constexpr_pi && "The argument of constexpr_sin_quad is assumed to be in [0,2*pi)" );
    }
    if( xx < 0 )
      return constexpr_sin_quad( mapIn0to2Pi( xx ), true );
    else if( xx < constexpr_pi_by_4 ) // [0/4*pi, 1/4*pi)
      return sinTaylor( xx );
    else if( xx < constexpr_pi_by_2 ) // [1/4*pi, 2/4*pi)
      return constexpr_sqrt( 1 - constexpr_pow( sinTaylor( constexpr_pi_by_2 - xx ), 2 ) );
    else if( xx < 3 * constexpr_pi_by_4 ) // [2/4*pi, 3/4*pi)
      return constexpr_sqrt( 1 - constexpr_pow( sinTaylor( xx - constexpr_pi_by_2 ), 2 ) );
    else if( xx < constexpr_pi ) // [3/4*pi, 4/4*pi)
      return sinTaylor( constexpr_pi - xx );
    else if( xx < 2 * constexpr_pi ) // [4/4*pi, 8/4*pi)
      return -constexpr_sin_quad( 2 * constexpr_pi - xx, true );
    else // [8/4*pi, +inf)
      return constexpr_sin_quad( mapIn0to2Pi( xx ), true );
  }

  // Constexpr implementation of sin (double signature, internally implemented as long double)
  constexpr double constexpr_sin( const double x )
  {
    return constexpr_sin_quad( x );
  }

  // Constexpr implementation of tan (long double signature)
  constexpr long double constexpr_tan_quad( const long double xx, const bool assume0to2Pi = false )
  {
    if( assume0to2Pi )
    {
      assert( xx >= 0 && "The argument of constexpr_sin_quad is assumed to be in [0,2*pi)" );
      assert( xx < 2 * constexpr_pi && "The argument of constexpr_sin_quad is assumed to be in [0,2*pi)" );
    }
    if( xx < 0 )
      return constexpr_tan_quad( mapIn0to2Pi( xx ), true );
    else if( xx < 2 * constexpr_pi ) // [0, 2*pi)
      return constexpr_sin_quad( xx, assume0to2Pi ) / constexpr_cos_quad( xx, assume0to2Pi );
    else // [8/4*pi, +inf)
      return constexpr_tan_quad( mapIn0to2Pi( xx ), true );
  }

  // Constexpr implementation of tan (double signature, internally implemented as long double)
  constexpr double constexpr_tan( const double x )
  {
    return constexpr_tan_quad( x );
  }

  // Constexpr implementation of atan for -1<x<1 (long double signature)
  // Taylor expansion : x - x**3/3 + x**5/5...
  constexpr long double atanTaylor( const long double xx )
  {
    assert( xx >= -1 && "The argument of atanTaylor is assumed to be in (-1,+1)" );
    assert( xx < 1 && "The argument of atanTaylor is assumed to be in (-1,+1)" );
    long double atanx = 0;
    int ipow = 1;
    long double xpow = xx;
    while( true )
    {
      long double atanxlast = atanx;
      atanx += xpow / ipow;
#ifdef CONSTEXPR_MATH_DEBUG
      std::cout << "ipow=" << ipow << ", xpow=" << xpow << ", atanx=" << atanx << std::endl; // for debugging (not constexpr)
#endif
      if( atanx == atanxlast ) break;
      // Next iteration
      ipow += 2;
      xpow *= -xx * xx;
    }
    return atanx;
  }

  // Constexpr implementation of atan (long double signature)
  constexpr long double constexpr_atan_quad( const long double xx )
  {
    if( xx > 1 )
      return constexpr_pi_by_2 - atanTaylor( 1 / xx );
    else if( xx == 1 )
      return constexpr_pi_by_4;
    else if( xx > -1 )
      return atanTaylor( xx );
    else if( xx == -1 )
      return -constexpr_pi_by_4;
    else // if( xx < -1 )
      return -constexpr_pi_by_2 - atanTaylor( 1 / xx );
  }

  // Constexpr implementation of atan (double signature, internally implemented as long double)
  constexpr double constexpr_atan( const double x )
  {
    return constexpr_atan_quad( x );
  }
}

#endif // constexpr_math_h
