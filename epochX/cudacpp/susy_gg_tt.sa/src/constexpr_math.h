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

  // Constexpr implementation of cos
  // Taylor expansion : 1 - x**2/2! + x**4/4!
  constexpr long double cosTaylor( const long double x )
  {
    long double cosx = 0;
    int ipow = 0;
    long double delta = 1;
    while( true )
    {
      long double cosxlast = cosx;
      cosx += delta;
      //std::cout << "ipow=" << ipow << ", delta=" << delta << ", cosx=" << cosx <<std::endl; // for debugging (not constexpr)
      if ( cosx == cosxlast ) break;
      // Next iteration
      ipow += 2;
      delta *= -x * x / ( ipow - 1 ) / ipow;
    }
    return cosx;
  }

  // Constexpr implementation of sin
  // Taylor expansion : x - x**3/3! + x**5/5!
  constexpr long double sinTaylor( const long double x )
  {
    long double sinx = 0;
    int ipow = 1;
    long double delta = x;
    while( true )
    {
      long double sinxlast = sinx;
      sinx += delta;
      //std::cout << "ipow=" << ipow << ", delta=" << delta << ", sinx=" << sinx <<std::endl; // for debugging (not constexpr)
      if ( sinx == sinxlast ) break;
      // Next iteration
      ipow += 2;
      delta *= -x * x / ( ipow - 1 ) / ipow;
    }
    return sinx;
  }

  // Constexpr implementation of cos (double signature, internally implemented as long double)
  constexpr double constexpr_cos( const long double x )
  {
    assert( x >= 0 && "The argument of constexpr_cos must be between 0 and pi/2" );
    assert( x <= constexpr_pi_by_2 && "The argument of constexpr_cos must be between 0 and pi/2" );
    //if ( x < constexpr_pi_by_4 ) return cosTaylor( (long double)x );
    //else return constexpr_sqrt( 1 - constexpr_pow( cosTaylor( constexpr_pi_by_2 - (long double)x ), 2 ) );
    if ( x < constexpr_pi_by_4 ) return constexpr_sqrt( 1 - constexpr_pow( sinTaylor( (long double)x ), 2 ) );
    else return sinTaylor( constexpr_pi_by_2 - (long double)x );
  }

  // Constexpr implementation of sin (double signature, internally implemented as long double)
  constexpr double constexpr_sin( const double x )
  {
    assert( x >= 0 && "The argument of constexpr_sin must be between 0 and pi/2" );
    assert( x <= constexpr_pi_by_2 && "The argument of constexpr_sin must be between 0 and pi/2" );
    if ( x < constexpr_pi_by_4 ) return sinTaylor( (long double)x );
    else return constexpr_sqrt( 1 - constexpr_pow( sinTaylor( constexpr_pi_by_2 - (long double)x ), 2 ) );
  }

  // Constexpr implementation of sin (double signature, internally implemented as long double)
  constexpr double constexpr_tan( const double x )
  {
    assert( x >= 0 && "The argument of constexpr_sin must be between 0 and pi/2" );
    assert( x <= constexpr_pi_by_2 && "The argument of constexpr_sin must be between 0 and pi/2" );
    if ( x < constexpr_pi_by_4 )
      return sinTaylor( (long double)x ) / constexpr_sqrt( 1 - constexpr_pow( sinTaylor( (long double)x ), 2 ) );
    else
      return constexpr_sqrt( 1 - constexpr_pow( sinTaylor( constexpr_pi_by_2 - (long double)x ), 2 ) ) / sinTaylor( constexpr_pi_by_2 - (long double)x );
  }
}

#endif // constexpr_math_h
