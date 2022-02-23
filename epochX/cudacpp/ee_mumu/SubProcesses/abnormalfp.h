#ifndef ABNORMALFP_H 
#define ABNORMALFP_H 1

#include "mgOnGpuFptypes.h"

#include <cmath>
#include <iostream>

//--------------------------------------------------------------------------

// Disabling fast math is essential here, otherwise results are undefined
// *** NB: these __attributes__ are now superseded by -fno-fast-math in the Makefile, but keep them anyway ***
// See https://stackoverflow.com/a/40702790 about __attribute__ on gcc
// See https://stackoverflow.com/a/32292725 about __attribute__ on clang
// (probably this does not work on clang? see https://groups.google.com/g/llvm-dev/c/Ys0hpgTFMH8)
#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
inline bool fp_is_nan( const fptype& fp )
{
  //#pragma clang diagnostic push
  //#pragma clang diagnostic ignored "-Wtautological-compare" // for icpx2021/clang13 (https://stackoverflow.com/a/15864661)
  return std::isnan( fp ); // always false for clang in fast math mode (tautological compare)?
  //#pragma clang diagnostic pop
}

//--------------------------------------------------------------------------

#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
inline bool fp_is_abnormal( const fptype& fp )
{
  if ( fp_is_nan( fp ) ) return true;
  if ( fp != fp ) return true;
  return false;
}

//--------------------------------------------------------------------------

#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
inline bool fp_is_zero( const fptype& fp )
{
  if ( fp == 0 ) return true;
  return false;
}

//--------------------------------------------------------------------------

// See https://en.cppreference.com/w/cpp/numeric/math/FP_categories
#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
inline const char* fp_show_class( const fptype& fp )
{
  switch( std::fpclassify( fp ) ) {
  case FP_INFINITE:  return "Inf";
  case FP_NAN:       return "NaN";
  case FP_NORMAL:    return "normal";
  case FP_SUBNORMAL: return "subnormal";
  case FP_ZERO:      return "zero";
  default:           return "unknown";
  }
}

//--------------------------------------------------------------------------

#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
inline void debug_me_is_abnormal( const fptype& me, size_t ievt )
{
  std::cout << "DEBUG[" << ievt << "]"
            << " ME=" << me
            << " fpisabnormal=" << fp_is_abnormal( me )
            << " fpclass=" << fp_show_class( me )
            << " (me==me)=" << ( me == me )
            << " (me==me+1)=" << ( me == me+1 )
            << " isnan=" << fp_is_nan( me )
            << " isfinite=" << std::isfinite( me )
            << " isnormal=" << std::isnormal( me )
            << " is0=" << ( me == 0 )
            << " is1=" << ( me == 1 )
            << " abs(ME)=" << std::abs( me )
            << " isnan=" << fp_is_nan( std::abs( me ) )
            << std::endl;
}

//--------------------------------------------------------------------------

#endif
