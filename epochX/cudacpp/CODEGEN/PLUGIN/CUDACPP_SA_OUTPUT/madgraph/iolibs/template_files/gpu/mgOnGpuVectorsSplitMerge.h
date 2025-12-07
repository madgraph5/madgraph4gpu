// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Nov 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#ifndef MGONGPUVECTORSSPLITMERGE_H
#define MGONGPUVECTORSSPLITMERGE_H 1

#include "mgOnGpuVectors.h"

// Disable all implementations
#undef MGONGPU_FPVFUN_EXPSIMD
#undef MGONGPU_FPVFUN_INTRINSICS
#undef MGONGPU_FPVFUN_SCALAR
#undef MGONGPU_FPVFUN_INITLIST

// Non-default implementation of fpvmerge using experimental simd (tested with gcc11)
//#define MGONGPU_FPVFUN_EXPSIMD 1 // NON-DEFAULT FOR TESTS

// Non-default implementation of fpvmerge using intrinsics (only on x86-64)
#ifdef __x86_64__
//#define MGONGPU_FPVFUN_INTRINSICS 1 // NON-DEFAULT FOR TESTS
#endif

// Non-default scalar implementation of fpvmerge for tests
//#define MGONGPU_FPVFUN_SCALAR 1 // NON-DEFAULT FOR TESTS

// Default implementation of fpvmerge using initializer lists
#define MGONGPU_FPVFUN_INITLIST 1 // DEFAULT

// SANITY CHECKS
#if defined MGONGPU_FPVFUN_EXPSIMD and ( defined MGONGPU_FPVFUN_INTRINSICS or defined MGONGPU_FPVFUN_SCALAR or defined MGONGPU_FPVFUN_INITLIST )
#error You must CHOOSE AT MOST ONE of MGONGPU_FPVFUN_EXPSIMD or MGONGPU_FPVFUN_INTRINSICS or MGONGPU_FPVFUN_SCALAR or MGONGPU_FPVFUN_INITLIST
#elif defined MGONGPU_FPVFUN_INTRINSICS and ( defined MGONGPU_FPVFUN_SCALAR or defined MGONGPU_FPVFUN_INITLIST )
#error You must CHOOSE AT MOST ONE of MGONGPU_FPVFUN_INTRINSICS or MGONGPU_FPVFUN_SCALAR or MGONGPU_FPVFUN_INITLIST
#elif defined MGONGPU_FPVFUN_SCALAR and defined MGONGPU_FPVFUN_INITLIST
#error You must CHOOSE AT MOST ONE of MGONGPU_FPVFUN_SCALAR or MGONGPU_FPVFUN_INITLIST
#endif

// Headers for intrinsics 
#ifdef MGONGPU_FPVFUN_INTRINSICS
#include <x86intrin.h>
#endif

// Headers for experimental simd
#ifdef MGONGPU_FPVFUN_EXPSIMD
#include <experimental/simd>
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
    // AV's original implementation with initializer lists (Oct 2022)
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
    // AV's implementation with x86-64 intrinsics (Nov 2025)
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
    __m128 f1 = _mm256_cvtpd_ps( u1.i );             // converts 4 doubles to 4 floats into __m128
    __m128 f2 = _mm256_cvtpd_ps( u2.i );             // converts 4 doubles to 4 floats into __m128
    __m256 f10 = _mm256_castps128_ps256( f1 );       // insert f1 into lower 128 bits of f12
    __m256 f12 = _mm256_insertf128_ps( f10, f2, 1 ); // copy f10 to f12 and insert f2 into higher 128 bits
    union { __m256 i; fptype2_v v; } u12;
    u12.i = f12;
    fptype2_v out = u12.v;
#elif MGONGPU_CPPSIMD == 8
    // --- CUDACPP "512z" ---
    union { fptype_v v; __m512d i; } u1, u2; // bitcast fptype_v to __512d
    u1.v = v1; u2.v = v2;
    __m256 f1 = _mm512_cvtpd_ps( u1.i );           // converts 8 doubles to 8 floats into __m256
    __m256 f2 = _mm512_cvtpd_ps( u2.i );           // converts 8 doubles to 8 floats into __m256
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

#ifdef MGONGPU_FPVFUN_EXPSIMD
  inline fptype2_v
  fpvmerge_expsimd( const fptype_v& v1, const fptype_v& v2 )
  {
    // AV's implementation with experimental simd (Nov 2025)
    namespace stdx = std::experimental;
    // Convert each fptype_v into a stdx::fixed_size_simd<fptype, n_d>
    constexpr size_t n_d = sizeof( fptype_v ) / sizeof( fptype ); // MGONGPU_CPPSIMD
    stdx::fixed_size_simd<fptype, n_d> sd1( reinterpret_cast<const fptype*>( &v1 ), stdx::element_aligned );
    stdx::fixed_size_simd<fptype, n_d> sd2( reinterpret_cast<const fptype*>( &v2 ), stdx::element_aligned );
    // Cast each stdx::fixed_size_simd<fptype, n_d> into a stdx::fixed_size_simd<fptype2, n_d>
    // (use static_simd_cast for vectorized double-to-float narrowing: simd_cast can only be used for non-narrowing casts)
    stdx::fixed_size_simd<fptype2, n_d> sf1 = stdx::static_simd_cast<stdx::fixed_size_simd<fptype2, n_d> >( sd1 );
    stdx::fixed_size_simd<fptype2, n_d> sf2 = stdx::static_simd_cast<stdx::fixed_size_simd<fptype2, n_d> >( sd2 );
    // Now concatenate sf1 (low half) and sf2 (high half) into one stdx::fixed_size_simd<fptype2, n_d*2>
    // Many TS implementations provide stdx::simd_cat, but some do not: do a safe copy to buffer instead
    fptype2_v out;
    sf1.copy_to( reinterpret_cast<fptype2*>( &out ), stdx::element_aligned );
    sf2.copy_to( reinterpret_cast<fptype2*>( &out ) + n_d, stdx::element_aligned );
    return out;
  }
#endif

  //--------------------------------------------------------------------------

  inline fptype2_v
  fpvmerge( const fptype_v& v1, const fptype_v& v2 )
  {
#ifdef MGONGPU_FPVFUN_SCALAR
    return fpvmerge_scalar( v1, v2 );
#elif defined MGONGPU_FPVFUN_EXPSIMD
    return fpvmerge_expsimd( v1, v2 );
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
#elif defined MGONGPU_FPVFUN_EXPSIMD
    //return fpvsplit0_expsimd( v );
    return fpvsplit0_initializerlist( v );
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
#elif defined MGONGPU_FPVFUN_EXPSIMD
    //return fpvsplit1_expsimd( v );
    return fpvsplit1_initializerlist( v );
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
