// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jan 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2023) for the MG5aMC CUDACPP plugin.

#include "MatrixElementKernels.h"

#include "CPPProcess.h"
#include "CudaRuntime.h"
#include "MemoryAccessMomenta.h"
#include "MemoryBuffers.h"

#include <sstream>

//============================================================================

#ifndef __CUDACC__
namespace mg5amcCpu
{

  //--------------------------------------------------------------------------

  MatrixElementKernelHost::MatrixElementKernelHost( const BufferMomenta& momenta,         // input: momenta
                                                    const BufferGs& gs,                   // input: gs for alphaS
                                                    const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                                                    const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                                                    BufferMatrixElements& matrixElements, // output: matrix elements
                                                    BufferSelectedHelicity& selhel,       // output: helicity selection
                                                    BufferSelectedColor& selcol,          // output: color selection
                                                    const size_t nevt )
    : MatrixElementKernelBase( momenta, gs, rndhel, rndcol, matrixElements, selhel, selcol )
    , NumberOfEvents( nevt )
    , m_couplings( nevt )
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    , m_numerators( nevt )
    , m_denominators( nevt )
#endif
  {
    if( m_momenta.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelHost: momenta must be a host array" );
    if( m_matrixElements.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelHost: matrixElements must be a host array" );
    if( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "MatrixElementKernelHost: nevt mismatch with momenta" );
    if( this->nevt() != m_matrixElements.nevt() ) throw std::runtime_error( "MatrixElementKernelHost: nevt mismatch with matrixElements" );
    // Sanity checks for memory access (momenta buffer)
    constexpr int neppM = MemoryAccessMomenta::neppM; // AOSOA layout
    static_assert( ispoweroftwo( neppM ), "neppM is not a power of 2" );
    if( nevt % neppM != 0 )
    {
      std::ostringstream sstr;
      sstr << "MatrixElementKernelHost: nevt should be a multiple of neppM=" << neppM;
      throw std::runtime_error( sstr.str() );
    }
    // Fail gently and avoid "Illegal instruction (core dumped)" if the host does not support the SIMD used in the ME calculation
    // Note: this prevents a crash on pmpe04 but not on some github CI nodes?
    // [NB: SIMD vectorization in mg5amc C++ code is only used in the ME calculation below MatrixElementKernelHost!]
    if( !MatrixElementKernelHost::hostSupportsSIMD() )
      throw std::runtime_error( "Host does not support the SIMD implementation of MatrixElementKernelsHost" );
  }

  //--------------------------------------------------------------------------

  int MatrixElementKernelHost::computeGoodHelicities()
  {
    constexpr int ncomb = CPPProcess::ncomb; // the number of helicity combinations
    HostBufferHelicityMask hstIsGoodHel( ncomb );
    // ... 0d1. Compute good helicity mask on the host
    computeDependentCouplings( m_gs.data(), m_couplings.data(), m_gs.size() );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    sigmaKin_getGoodHel( m_momenta.data(), m_couplings.data(), m_matrixElements.data(), m_numerators.data(), m_denominators.data(), hstIsGoodHel.data(), nevt() );
#else
    sigmaKin_getGoodHel( m_momenta.data(), m_couplings.data(), m_matrixElements.data(), hstIsGoodHel.data(), nevt() );
#endif
    // ... 0d2. Copy back good helicity list to static memory on the host
    // [FIXME! REMOVE THIS STATIC THAT BREAKS MULTITHREADING?]
    return sigmaKin_setGoodHel( hstIsGoodHel.data() );
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelHost::computeMatrixElements( const unsigned int channelId )
  {
    computeDependentCouplings( m_gs.data(), m_couplings.data(), m_gs.size() );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    sigmaKin( m_momenta.data(), m_couplings.data(), m_rndhel.data(), m_rndcol.data(), m_matrixElements.data(), channelId, m_numerators.data(), m_denominators.data(), m_selhel.data(), m_selcol.data(), nevt() );
#else
    sigmaKin( m_momenta.data(), m_couplings.data(), m_rndhel.data(), m_rndcol.data(), m_matrixElements.data(), m_selhel.data(), m_selcol.data(), nevt() );
#endif
  }

  //--------------------------------------------------------------------------

  // Does this host system support the SIMD used in the matrix element calculation?
  bool MatrixElementKernelHost::hostSupportsSIMD( const bool verbose )
  {
#if defined __AVX512VL__
    bool known = true;
    bool ok = __builtin_cpu_supports( "avx512vl" );
    const std::string tag = "skylake-avx512 (AVX512VL)";
#elif defined __AVX2__
    bool known = true;
    bool ok = __builtin_cpu_supports( "avx2" );
    const std::string tag = "haswell (AVX2)";
#elif defined __SSE4_2__
#ifdef __PPC__
    // See https://gcc.gnu.org/onlinedocs/gcc/Basic-PowerPC-Built-in-Functions-Available-on-all-Configurations.html
    bool known = true;
    bool ok = __builtin_cpu_supports( "vsx" );
    const std::string tag = "powerpc vsx (128bit as in SSE4.2)";
#elif defined __ARM_NEON__ // consider using __BUILTIN_CPU_SUPPORTS__
    bool known = false; // __builtin_cpu_supports is not supported
    // See https://gcc.gnu.org/onlinedocs/gcc/Basic-PowerPC-Built-in-Functions-Available-on-all-Configurations.html
    // See https://stackoverflow.com/q/62783908
    // See https://community.arm.com/arm-community-blogs/b/operating-systems-blog/posts/runtime-detection-of-cpu-features-on-an-armv8-a-cpu
    bool ok = true; // this is just an assumption!
    const std::string tag = "arm neon (128bit as in SSE4.2)";
#elif defined( __x86_64__ ) || defined( __i386__ )
    bool known = true;
    bool ok = __builtin_cpu_supports( "sse4.2" );
    const std::string tag = "nehalem (SSE4.2)";
#else // AV FIXME! Added by OM for Mac, should identify the correct __xxx__ flag that should be targeted
    bool known = false; // __builtin_cpu_supports is not supported
    // See https://gcc.gnu.org/onlinedocs/gcc/Basic-PowerPC-Built-in-Functions-Available-on-all-Configurations.html
    // See https://stackoverflow.com/q/62783908
    // See https://community.arm.com/arm-community-blogs/b/operating-systems-blog/posts/runtime-detection-of-cpu-features-on-an-armv8-a-cpu
    bool ok = true; // this is just an assumption!
    const std::string tag = "arm neon (128bit as in SSE4.2)";
#endif
#else
    bool known = true;
    bool ok = true;
    const std::string tag = "none";
#endif
    if( verbose )
    {
      if( tag == "none" )
        std::cout << "INFO: The application does not require the host to support any AVX feature" << std::endl;
      else if( ok && known )
        std::cout << "INFO: The application is built for " << tag << " and the host supports it" << std::endl;
      else if( ok )
        std::cout << "WARNING: The application is built for " << tag << " but it is unknown if the host supports it" << std::endl;
      else
        std::cout << "ERROR! The application is built for " << tag << " but the host does not support it" << std::endl;
    }
    return ok;
  }

  //--------------------------------------------------------------------------

}
#endif

//============================================================================

#ifdef __CUDACC__
namespace mg5amcGpu
{

  //--------------------------------------------------------------------------

  MatrixElementKernelDevice::MatrixElementKernelDevice( const BufferMomenta& momenta,         // input: momenta
                                                        const BufferGs& gs,                   // input: gs for alphaS
                                                        const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                                                        const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                                                        BufferMatrixElements& matrixElements, // output: matrix elements
                                                        BufferSelectedHelicity& selhel,       // output: helicity selection
                                                        BufferSelectedColor& selcol,          // output: color selection
                                                        const size_t gpublocks,
                                                        const size_t gputhreads )
    : MatrixElementKernelBase( momenta, gs, rndhel, rndcol, matrixElements, selhel, selcol )
    , NumberOfEvents( gpublocks * gputhreads )
    , m_couplings( this->nevt() )
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    , m_numerators( this->nevt() )
    , m_denominators( this->nevt() )
#endif
    , m_gpublocks( gpublocks )
    , m_gputhreads( gputhreads )
  {
    if( !m_momenta.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelDevice: momenta must be a device array" );
    if( !m_matrixElements.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelDevice: matrixElements must be a device array" );
    if( m_gpublocks == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gpublocks must be > 0" );
    if( m_gputhreads == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gputhreads must be > 0" );
    if( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "MatrixElementKernelDevice: nevt mismatch with momenta" );
    if( this->nevt() != m_matrixElements.nevt() ) throw std::runtime_error( "MatrixElementKernelDevice: nevt mismatch with matrixElements" );
    // Sanity checks for memory access (momenta buffer)
    constexpr int neppM = MemoryAccessMomenta::neppM; // AOSOA layout
    static_assert( ispoweroftwo( neppM ), "neppM is not a power of 2" );
    if( m_gputhreads % neppM != 0 )
    {
      std::ostringstream sstr;
      sstr << "MatrixElementKernelHost: gputhreads should be a multiple of neppM=" << neppM;
      throw std::runtime_error( sstr.str() );
    }
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelDevice::setGrid( const int gpublocks, const int gputhreads )
  {
    if( m_gpublocks == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gpublocks must be > 0 in setGrid" );
    if( m_gputhreads == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gputhreads must be > 0 in setGrid" );
    if( this->nevt() != m_gpublocks * m_gputhreads ) throw std::runtime_error( "MatrixElementKernelDevice: nevt mismatch in setGrid" );
  }

  //--------------------------------------------------------------------------

  int MatrixElementKernelDevice::computeGoodHelicities()
  {
    constexpr int ncomb = CPPProcess::ncomb; // the number of helicity combinations
    PinnedHostBufferHelicityMask hstIsGoodHel( ncomb );
    DeviceBufferHelicityMask devIsGoodHel( ncomb );
    // ... 0d1. Compute good helicity mask on the device
    computeDependentCouplings<<<m_gpublocks, m_gputhreads>>>( m_gs.data(), m_couplings.data() );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    sigmaKin_getGoodHel<<<m_gpublocks, m_gputhreads>>>( m_momenta.data(), m_couplings.data(), m_matrixElements.data(), m_numerators.data(), m_denominators.data(), devIsGoodHel.data() );
#else
    sigmaKin_getGoodHel<<<m_gpublocks, m_gputhreads>>>( m_momenta.data(), m_couplings.data(), m_matrixElements.data(), devIsGoodHel.data() );
#endif
    checkCuda( cudaPeekAtLastError() );
    // ... 0d2. Copy back good helicity mask to the host
    copyHostFromDevice( hstIsGoodHel, devIsGoodHel );
    // ... 0d3. Copy back good helicity list to constant memory on the device
    return sigmaKin_setGoodHel( hstIsGoodHel.data() );
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelDevice::computeMatrixElements( const unsigned int channelId )
  {
    computeDependentCouplings<<<m_gpublocks, m_gputhreads>>>( m_gs.data(), m_couplings.data() );
#ifndef MGONGPU_NSIGHT_DEBUG
    constexpr unsigned int sharedMemSize = 0;
#else
    constexpr unsigned int sharedMemSize = ntpbMAX * sizeof( float );
#endif
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    sigmaKin<<<m_gpublocks, m_gputhreads, sharedMemSize>>>( m_momenta.data(), m_couplings.data(), m_rndhel.data(), m_rndcol.data(), m_matrixElements.data(), channelId, m_numerators.data(), m_denominators.data(), m_selhel.data(), m_selcol.data() );
#else
    sigmaKin<<<m_gpublocks, m_gputhreads, sharedMemSize>>>( m_momenta.data(), m_couplings.data(), m_rndhel.data(), m_rndcol.data(), m_matrixElements.data(), m_selhel.data(), m_selcol.data() );
#endif
    checkCuda( cudaPeekAtLastError() );
    checkCuda( cudaDeviceSynchronize() );
  }

  //--------------------------------------------------------------------------

}
#endif

//============================================================================
