// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Nov 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#ifndef MGONGPUVECTORSSPLITMERGE_H
#define MGONGPUVECTORSSPLITMERGE_H 1

#include "mgOnGpuVectors.h"

// Disable all implementations
#undef MGONGPU_FPVFUN_INTRINSICS
#undef MGONGPU_FPVFUN_SCALAR
#undef MGONGPU_FPVFUN_INITLIST

// Non-default implementation of fpvmerge using intrinsics (only on x86-64)
#ifdef __x86_64__
//#define MGONGPU_FPVFUN_INTRINSICS 1 // NON-DEFAULT FOR TESTS
#endif

// Non-default scalar implementation of fpvmerge for tests
#define MGONGPU_FPVFUN_SCALAR 1 // NON-DEFAULT FOR TESTS

// Default implementation of fpvmerge using initializer lists
//#define MGONGPU_FPVFUN_INITLIST 1 // DEFAULT

// SANITY CHECKS
#if defined MGONGPU_FPVFUN_INTRINSICS and ( defined MGONGPU_FPVFUN_SCALAR or defined MGONGPU_FPVFUN_INITLIST )
#error You must CHOOSE AT MOST ONE of MGONGPU_FPVFUN_INTRINSICS or MGONGPU_FPVFUN_SCALAR or MGONGPU_FPVFUN_INITLIST
#elif defined MGONGPU_FPVFUN_SCALAR and defined MGONGPU_FPVFUN_INITLIST
#error You must CHOOSE AT MOST ONE of MGONGPU_FPVFUN_SCALAR or MGONGPU_FPVFUN_INITLIST
#endif

// Headers for intrinsics 
#ifdef MGONGPU_FPVFUN_INTRINSICS
#include <x86intrin.h>
#endif

//==========================================================================

#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
namespace mg5amcCpu
{
  //--------------------------------------------------------------------------

  inline fptype2_v
  fpvmerge_scalar( const fptype_v& v1, const fptype_v& v2 )
  {
    // Scalar implementation for sanity checks (slower? auto-vectorized?)
    fptype2_v out;
    for( int ieppV = 0; ieppV < neppV; ieppV++ )
    {
      out[ieppV] = v1[ieppV];
      out[ieppV+neppV] = v2[ieppV];
    }
    return out;
  }  

  //--------------------------------------------------------------------------

  inline fptype2_v
  fpvmerge_initializerlist( const fptype_v& v1, const fptype_v& v2 )
  {
    // AV's original implementation with initializer lists
    // I initially thought that this was inefficient as it seemed as slow as double (#537)
    // Later tests show that this is as fast as intrinsics and faster than experimental SIMD
#if MGONGPU_CPPSIMD == 2
    // --- CUDACPP "sse4" ---
    fptype2_v out =
      { (fptype2)v1[0], (fptype2)v1[1], (fptype2)v2[0], (fptype2)v2[1] };
#elif MGONGPU_CPPSIMD == 4
    // --- CUDACPP "avx2" or "512y" ---
    fptype2_v out =
      { (fptype2)v1[0], (fptype2)v1[1], (fptype2)v1[2], (fptype2)v1[3], (fptype2)v2[0], (fptype2)v2[1], (fptype2)v2[2], (fptype2)v2[3] };
#elif MGONGPU_CPPSIMD == 8
    // --- CUDACPP "512z" ---
    fptype2_v out =
      { (fptype2)v1[0], (fptype2)v1[1], (fptype2)v1[2], (fptype2)v1[3], (fptype2)v1[4], (fptype2)v1[5], (fptype2)v1[6], (fptype2)v1[7], (fptype2)v2[0], (fptype2)v2[1], (fptype2)v2[2], (fptype2)v2[3], (fptype2)v2[4], (fptype2)v2[5], (fptype2)v2[6], (fptype2)v2[7] };
#endif
    return out;
  }

  //--------------------------------------------------------------------------

#ifdef MGONGPU_FPVFUN_INTRINSICS
  inline fptype2_v
  fpvmerge_intrinsics( const fptype_v& v1, const fptype_v& v2 )
  {
    // AV's implementation with x86-64 intrinsics
#if MGONGPU_CPPSIMD == 2
    // --- CUDACPP "sse4" ---
    union { fptype_v v; __m128d i; } u1, u2; // bitcast fptype_v to __m128d
    u1.v = v1; u2.v = v2;
    __m128 f10 = _mm_cvtpd_ps( u1.i );      // converts 2 doubles to 2 floats into lower 64 bits of of __m128
    __m128 f20 = _mm_cvtpd_ps( u2.i );      // converts 2 doubles to 2 floats into lower 64 bits of of __m128
    __m128 f12 = _mm_movelh_ps( f10, f20 ); // places lower half of f10 then lower half of f20 into f12
    union { __m128 i; fptype2_v v; } u12;
    u12.i = f12;
    fptype2_v out = u12.v;
#elif MGONGPU_CPPSIMD == 4
    // --- CUDACPP "avx2" or "512y" ---
    union { fptype_v v; __m256d i; } u1, u2; // bitcast fptype_v to __m256d
    u1.v = v1; u2.v = v2;
    __m128 f1 = _mm256_cvtpd_ps( u1.i );             // converts 2 doubles to 2 floats
    __m128 f2 = _mm256_cvtpd_ps( u2.i );             // converts 2 doubles to 2 floats
    __m256 f10 = _mm256_castps128_ps256( f1 );       // insert f1 into lower 128 bits of f12
    __m256 f12 = _mm256_insertf128_ps( f10, f2, 1 ); // copy f10 to f12 and insert f2 into higher 128 bits
    union { __m256 i; fptype2_v v; } u12;
    u12.i = f12;
    fptype2_v out = u12.v;
#elif MGONGPU_CPPSIMD == 8
    // --- CUDACPP "512z" ---
    union { fptype_v v; __m512d i; } u1, u2; // bitcast fptype_v to __512d
    u1.v = v1; u2.v = v2;
    __m256 f1 = _mm512_cvtpd_ps( u1.i );           // converts 2 doubles to 2 floats
    __m256 f2 = _mm512_cvtpd_ps( u2.i );           // converts 2 doubles to 2 floats
    __m512 f10 = _mm512_castps256_ps512( f1 );     // insert f1 into lower 256 bits of f12
    __m512 f12 = _mm512_insertf32x8( f10, f2, 1 ); // copy f10 to f12 and insert f2 into higher 256 bits
    union { __m512 i; fptype2_v v; } u12;
    u12.i = f12;
    fptype2_v out = u12.v;
#endif
    return out;
  }
#endif

  //--------------------------------------------------------------------------

  inline fptype2_v
  fpvmerge( const fptype_v& v1, const fptype_v& v2 )
  {
#ifdef MGONGPU_FPVFUN_SCALAR
    return fpvmerge_scalar( v1, v2 );
#elif defined MGONGPU_FPVFUN_INTRINSICS
    return fpvmerge_intrinsics( v1, v2 );
#elif defined MGONGPU_FPVFUN_INITLIST
    return fpvmerge_initializerlist( v1, v2 );
#else
#error No implementation found for fpvmerge
#endif
  }

  //--------------------------------------------------------------------------

  inline fptype_v
  fpvsplit0_scalar( const fptype2_v& v )
  {
    fptype_v out = {};
    for( int ieppV = 0; ieppV < neppV; ieppV++ )
    {
      out[ieppV] = v[ieppV];
    }
    return out;
  }

  //--------------------------------------------------------------------------

  inline fptype_v
  fpvsplit0_initializerlist( const fptype2_v& v )
  {
#if MGONGPU_CPPSIMD == 2
    fptype_v out =
      { (fptype)v[0], (fptype)v[1] };
#elif MGONGPU_CPPSIMD == 4
    fptype_v out =
      { (fptype)v[0], (fptype)v[1], (fptype)v[2], (fptype)v[3] };
#elif MGONGPU_CPPSIMD == 8
    fptype_v out =
      { (fptype)v[0], (fptype)v[1], (fptype)v[2], (fptype)v[3], (fptype)v[4], (fptype)v[5], (fptype)v[6], (fptype)v[7] };
#endif
    return out;
  }

  //--------------------------------------------------------------------------

  inline fptype_v
  fpvsplit0( const fptype2_v& v )
  {
#ifdef MGONGPU_FPVFUN_SCALAR
    return fpvsplit0_scalar( v );
#elif defined MGONGPU_FPVFUN_INTRINSICS
    //return fpvsplit0_intrinsics( v );
    return fpvsplit0_initializerlist( v );
#elif defined MGONGPU_FPVFUN_INITLIST
    return fpvsplit0_initializerlist( v );
#else
#error No implementation found for fpvsplit0
#endif
  }

  //--------------------------------------------------------------------------

  inline fptype_v
  fpvsplit1_scalar( const fptype2_v& v )
  {
    fptype_v out = {};
    for( int ieppV = 0; ieppV < neppV; ieppV++ )
    {
      out[ieppV] = v[ieppV+neppV];
    }
    return out;
  }

  //--------------------------------------------------------------------------

  inline fptype_v
  fpvsplit1_initializerlist( const fptype2_v& v )
  {
#if MGONGPU_CPPSIMD == 2
    fptype_v out =
      { (fptype)v[2], (fptype)v[3] };
#elif MGONGPU_CPPSIMD == 4
    fptype_v out =
      { (fptype)v[4], (fptype)v[5], (fptype)v[6], (fptype)v[7] };
#elif MGONGPU_CPPSIMD == 8
    fptype_v out =
      { (fptype)v[8], (fptype)v[9], (fptype)v[10], (fptype)v[11], (fptype)v[12], (fptype)v[13], (fptype)v[14], (fptype)v[15] };
#endif
    return out;
  }

  //--------------------------------------------------------------------------

  inline fptype_v
  fpvsplit1( const fptype2_v& v )
  {
#ifdef MGONGPU_FPVFUN_SCALAR
    return fpvsplit1_scalar( v );
#elif defined MGONGPU_FPVFUN_INTRINSICS
    //return fpvsplit1_intrinsics( v );
    return fpvsplit1_initializerlist( v );
#elif defined MGONGPU_FPVFUN_INITLIST
    return fpvsplit1_initializerlist( v );
#else
#error No implementation found for fpvsplit1
#endif
  }

  //--------------------------------------------------------------------------
}
#endif

#endif // MGONGPUVECTORSSPLITMERGE_H
