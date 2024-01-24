// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jul 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, O. Mattelaer, S. Roiser, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.

#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

// HARDCODED AT CODE GENERATION TIME: DO NOT MODIFY (#473)
// There are two different code bases for standalone_cudacpp (without multichannel) and madevent+cudacpp (with multichannel)
#define MGONGPU_SUPPORTS_MULTICHANNEL 1

// ** NB1 Throughputs (e.g. 6.8E8) are events/sec for "./gcheck.exe -p 65536 128 12"
// ** NB2 Baseline on b7g47n0004 fluctuates (probably depends on load on other VMs)

// Choose if curand is supported for generating random numbers
// For both CUDA and C++, by default, do not inline, but allow this macro to be set from outside with e.g. -DMGONGPU_HAS_NO_CURAND
// (there exist CUDA installations, e.g. using the HPC package, which do not include curand - see PR #784)
//#ifdef __CUDACC__
//#undef MGONGPU_HAS_NO_CURAND // default
////#define MGONGPU_HAS_NO_CURAND 1
//#else
//#undef MGONGPU_HAS_NO_CURAND // default
////#define MGONGPU_HAS_NO_CURAND 1
//#endif

// Choose floating point precision (for everything but color algebra #537)
// If one of these macros has been set from outside with e.g. -DMGONGPU_FPTYPE_FLOAT, nothing happens (issue #167)
#if not defined MGONGPU_FPTYPE_DOUBLE and not defined MGONGPU_FPTYPE_FLOAT
// Floating point precision (CHOOSE ONLY ONE)
#define MGONGPU_FPTYPE_DOUBLE 1 // default
//#define MGONGPU_FPTYPE_FLOAT 1 // 2x faster
#endif

// Choose floating point precision (for color algebra alone #537)
// If one of these macros has been set from outside with e.g. -DMGONGPU_FPTYPE2_FLOAT, nothing happens (issue #167)
#if not defined MGONGPU_FPTYPE2_DOUBLE and not defined MGONGPU_FPTYPE2_FLOAT
// Floating point precision (CHOOSE ONLY ONE)
#define MGONGPU_FPTYPE2_DOUBLE 1 // default
//#define MGONGPU_FPTYPE2_FLOAT 1 // 2x faster
#endif

// Choose whether to inline all HelAmps functions
// This optimization can gain almost a factor 4 in C++, similar to -flto (issue #229)
// By default, do not inline, but allow this macro to be set from outside with e.g. -DMGONGPU_INLINE_HELAMPS
//#undef MGONGPU_INLINE_HELAMPS // default
////#define MGONGPU_INLINE_HELAMPS 1

// Choose whether to hardcode the cIPD physics parameters rather than reading them from user cards
// This optimization can gain 20% in CUDA in eemumu (issue #39)
// By default, do not hardcode, but allow this macro to be set from outside with e.g. -DMGONGPU_HARDCODE_PARAM
// ** NB: The option to use hardcoded cIPD physics parameters is supported again even now when alphas is running (#373)
// ** NB: Note however that it now only refers to cIPD parameters (cIPC parameters are always accessed through global memory)
//#undef MGONGPU_HARDCODE_PARAM // default
////#define MGONGPU_HARDCODE_PARAM 1

// Complex type in c++: std::complex or cxsmpl (CHOOSE ONLY ONE)
#ifndef __CUDACC__
//#define MGONGPU_CPPCXTYPE_STDCOMPLEX 1 // ~8 percent slower on float, same on double (5.1E6/double, 9.4E6/float)
#define MGONGPU_CPPCXTYPE_CXSMPL 1 // new default (5.1E6/double, 10.2E6/float)
#endif

// Complex type in cuda: thrust or cucomplex or cxsmpl (CHOOSE ONLY ONE)
#ifdef __CUDACC__
#define MGONGPU_CUCXTYPE_THRUST 1 // default (~1.15E9/double, ~3.2E9/float)
//#define MGONGPU_CUCXTYPE_CUCOMPLEX 1 // ~10 percent slower (1.03E9/double, ~2.8E9/float)
//#define MGONGPU_CUCXTYPE_CXSMPL 1 // ~10 percent slower (1.00E9/double, ~2.9E9/float)
#endif

// Cuda nsight compute (ncu) debug: add dummy lines to ease SASS program flow navigation
#ifdef __CUDACC__
#undef MGONGPU_NSIGHT_DEBUG // default
//#define MGONGPU_NSIGHT_DEBUG 1
#endif

// SANITY CHECKS (floating point precision for everything but color algebra #537)
#if defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE_FLOAT
#error You must CHOOSE (ONE AND) ONLY ONE of MGONGPU_FPTYPE_DOUBLE or defined MGONGPU_FPTYPE_FLOAT
#endif

// SANITY CHECKS (floating point precision for color algebra alone #537)
#if defined MGONGPU_FPTYPE2_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
#error You must CHOOSE (ONE AND) ONLY ONE of MGONGPU_FPTYPE2_DOUBLE or defined MGONGPU_FPTYPE2_FLOAT
#endif
#if defined MGONGPU_FPTYPE2_DOUBLE and defined MGONGPU_FPTYPE_FLOAT
#error You cannot use double precision for color algebra and single precision elsewhere
#endif

// SANITY CHECKS (c++ complex number implementation)
#ifndef __CUDACC__
#if defined MGONGPU_CPPCXTYPE_STDCOMPLEX and defined MGONGPU_CPPCXTYPE_CXSMPL
#error You must CHOOSE (ONE AND) ONLY ONE of MGONGPU_CPPCXTYPE_STDCOMPLEX or MGONGPU_CPPCXTYPE_CXSMPL
#endif
#endif

// SANITY CHECKS (cuda complex number implementation)
#ifdef __CUDACC__
#if defined MGONGPU_CUCXTYPE_THRUST and defined MGONGPU_CUCXTYPE_CUCOMPLEX and defined MGONGPU_CUCXTYPE_CXSMPL
#error You must CHOOSE (ONE AND) ONLY ONE of MGONGPU_CUCXTYPE_THRUST or MGONGPU_CUCXTYPE_CUCOMPLEX or MGONGPU_CUCXTYPE_CXSMPL
#endif
#endif

// NB: namespace mgOnGpu includes types which are defined in exactly the same way for CPU and GPU builds (see #318 and #725)
namespace mgOnGpu
{

  // --- Type definitions

  // Floating point type (for everything but color algebra #537): fptype
#if defined MGONGPU_FPTYPE_DOUBLE
  typedef double fptype; // double precision (8 bytes, fp64)
#elif defined MGONGPU_FPTYPE_FLOAT
  typedef float fptype;  // single precision (4 bytes, fp32)
#endif

  // Floating point type (for color algebra alone #537): fptype2
#if defined MGONGPU_FPTYPE2_DOUBLE
  typedef double fptype2; // double precision (8 bytes, fp64)
#elif defined MGONGPU_FPTYPE2_FLOAT
  typedef float fptype2; // single precision (4 bytes, fp32)
#endif

  // --- Platform-specific software implementation details

  // Maximum number of blocks per grid
  // ** NB Some arrays of pointers will be allocated statically to fit all these blocks
  // ** (the actual memory for each block will then be allocated dynamically only for existing blocks)
  //const int nbpgMAX = 2048;

  // Maximum number of threads per block
  //const int ntpbMAX = 256; // AV Apr2021: why had I set this to 256?
  const int ntpbMAX = 1024; // NB: 512 is ok, but 1024 does fail with "too many resources requested for launch"

  // Alignment requirement for using reinterpret_cast with SIMD vectorized code
  // (using reinterpret_cast with non aligned memory may lead to segmentation faults!)
  // Only needed for C++ code but can be enforced also in NVCC builds of C++ code using CUDA>=11.2 and C++17 (#318, #319, #333)
#ifndef __CUDACC__
  constexpr int cppAlign = 64; // alignment requirement for SIMD vectorization (64-byte i.e. 512-bit)
#endif

}

// Expose typedefs and operators outside the namespace
using mgOnGpu::fptype;
using mgOnGpu::fptype2;

// C++ SIMD vectorization width (this will be used to set neppV)
#ifdef __CUDACC__ // CUDA implementation has no SIMD
#undef MGONGPU_CPPSIMD
#elif defined __AVX512VL__ && defined MGONGPU_PVW512 // C++ "512z" AVX512 with 512 width (512-bit ie 64-byte): 8 (DOUBLE) or 16 (FLOAT)
#ifdef MGONGPU_FPTYPE_DOUBLE
#define MGONGPU_CPPSIMD 8
#else
#define MGONGPU_CPPSIMD 16
#endif
#elif defined __AVX512VL__ // C++ "512y" AVX512 with 256 width (256-bit ie 32-byte): 4 (DOUBLE) or 8 (FLOAT) [gcc DEFAULT]
#ifdef MGONGPU_FPTYPE_DOUBLE
#define MGONGPU_CPPSIMD 4
#else
#define MGONGPU_CPPSIMD 8
#endif
#elif defined __AVX2__ // C++ "avx2" AVX2 (256-bit ie 32-byte): 4 (DOUBLE) or 8 (FLOAT) [clang DEFAULT]
#ifdef MGONGPU_FPTYPE_DOUBLE
#define MGONGPU_CPPSIMD 4
#else
#define MGONGPU_CPPSIMD 8
#endif
#elif defined __SSE4_2__ // C++ "sse4" SSE4.2 (128-bit ie 16-byte): 2 (DOUBLE) or 4 (FLOAT) [Power9 and ARM default]
#ifdef MGONGPU_FPTYPE_DOUBLE
#define MGONGPU_CPPSIMD 2
#else
#define MGONGPU_CPPSIMD 4
#endif
#else // C++ "none" i.e. no SIMD
#undef MGONGPU_CPPSIMD
#endif

// Cuda nsight compute (ncu) debug: add dummy lines to ease SASS program flow navigation
// Arguments (not used so far): text is __FUNCTION__, code is 0 (start) or 1 (end)
#if defined __CUDACC__ && defined MGONGPU_NSIGHT_DEBUG /* clang-format off */
#define mgDebugDeclare() __shared__ float mgDebugCounter[mgOnGpu::ntpbMAX];
#define mgDebugInitialise() { mgDebugCounter[threadIdx.x] = 0; }
#define mgDebug( code, text ) { mgDebugCounter[threadIdx.x] += 1; }
#define mgDebugFinalise() { if ( blockIdx.x == 0 && threadIdx.x == 0 ) printf( "MGDEBUG: counter=%f\n", mgDebugCounter[threadIdx.x] ); }
#else
#define mgDebugDeclare() /*noop*/
#define mgDebugInitialise() { /*noop*/ }
#define mgDebug( code, text ) { /*noop*/ }
#define mgDebugFinalise() { /*noop*/ }
#endif /* clang-format on */

// Define empty CUDA declaration specifiers for C++
#ifndef __CUDACC__
#define __global__
#define __host__
#define __device__
#endif

// For SANITY CHECKS: check that neppR, neppM, neppV... are powers of two (https://stackoverflow.com/a/108360)
inline constexpr bool
ispoweroftwo( int n )
{
  return ( n > 0 ) && !( n & ( n - 1 ) );
}

// Compiler version support (#96): require nvcc from CUDA >= 11.2, e.g. to use C++17 (see #333)
#ifdef __NVCC__
#if( __CUDACC_VER_MAJOR__ < 11 ) || ( __CUDACC_VER_MAJOR__ == 11 && __CUDACC_VER_MINOR__ < 2 )
#error Unsupported CUDA version: please use CUDA >= 11.2
#endif
#endif

// Compiler version support (#96): require clang >= 11
#if defined __clang__
#if( __clang_major__ < 11 )
#error Unsupported clang version: please use clang >= 11
#endif
// Compiler version support (#96): require gcc >= 9.3, e.g. for some OMP issues (see #269)
// [NB skip this check for the gcc toolchain below clang or icx (TEMPORARY? #355)]
#elif defined __GNUC__
#if( __GNUC__ < 9 ) || ( __GNUC__ == 9 && __GNUC_MINOR__ < 3 )
#error Unsupported gcc version: please gcc >= 9.3
#endif
#endif

#endif // MGONGPUCONFIG_H
