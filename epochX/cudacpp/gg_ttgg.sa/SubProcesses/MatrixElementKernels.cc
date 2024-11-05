// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jan 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: J. Teig, A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#include "MatrixElementKernels.h"

#include "CPPProcess.h"
#include "GpuRuntime.h" // Includes the abstraction for Nvidia/AMD compilation
#include "MemoryAccessMomenta.h"
#include "MemoryBuffers.h"

#include <cfenv> // for fetestexcept
#include <iostream>
#include <sstream>

//============================================================================

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

  MatrixElementKernelBase::MatrixElementKernelBase( const BufferMomenta& momenta,         // input: momenta
                                                    const BufferGs& gs,                   // input: gs for alphaS
                                                    const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                                                    const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                                                    const BufferChannelIds& channelIds,   // input: channel ids for single-diagram enhancement
                                                    BufferMatrixElements& matrixElements, // output: matrix elements
                                                    BufferSelectedHelicity& selhel,       // output: helicity selection
                                                    BufferSelectedColor& selcol )         // output: color selection
    : m_momenta( momenta )
    , m_gs( gs )
    , m_rndhel( rndhel )
    , m_rndcol( rndcol )
    , m_channelIds( channelIds )
    , m_matrixElements( matrixElements )
    , m_selhel( selhel )
    , m_selcol( selcol )
#ifdef MGONGPU_CHANNELID_DEBUG
    , m_nevtProcessedByChannel()
    , m_tag()
#endif
  {
    //std::cout << "DEBUG: MatrixElementKernelBase ctor " << this << std::endl;
#ifdef MGONGPU_CHANNELID_DEBUG
    for( size_t channelId = 0; channelId < CPPProcess::ndiagrams + 1; channelId++ ) // [0...ndiagrams] (TEMPORARY: 0=multichannel)
      m_nevtProcessedByChannel[channelId] = 0;
#endif
  }

  //--------------------------------------------------------------------------

  MatrixElementKernelBase::~MatrixElementKernelBase()
  {
    //std::cout << "DEBUG: MatrixElementKernelBase dtor " << this << std::endl;
#ifdef MGONGPU_CHANNELID_DEBUG
    MatrixElementKernelBase::dumpNevtProcessedByChannel();
#endif
    MatrixElementKernelBase::dumpSignallingFPEs();
  }

  //--------------------------------------------------------------------------

#ifdef MGONGPU_CHANNELID_DEBUG
  void MatrixElementKernelBase::updateNevtProcessedByChannel( const unsigned int* pHstChannelIds, const size_t nevt )
  {
    if( pHstChannelIds != nullptr )
    {
      //std::cout << "DEBUG " << this << ": not nullptr " << nevt << std::endl;
      for( unsigned int ievt = 0; ievt < nevt; ievt++ )
      {
        const size_t channelId = pHstChannelIds[ievt]; // Fortran indexing
        //assert( channelId > 0 );
        //assert( channelId < CPPProcess::ndiagrams );
        m_nevtProcessedByChannel[channelId]++;
      }
    }
    else
    {
      //std::cout << "DEBUG " << this << ": nullptr " << std::endl;
      m_nevtProcessedByChannel[0] += nevt;
    }
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPU_CHANNELID_DEBUG
  void MatrixElementKernelBase::dumpNevtProcessedByChannel()
  {
    size_t nevtProcessed = 0;
    for( size_t channelId = 0; channelId < CPPProcess::ndiagrams + 1; channelId++ ) // [0...ndiagrams] (TEMPORARY: 0=multichannel)
      nevtProcessed += m_nevtProcessedByChannel[channelId];
    std::ostringstream sstr;
    sstr << " {";
    for( size_t channelId = 0; channelId < CPPProcess::ndiagrams + 1; channelId++ ) // [0...ndiagrams] (TEMPORARY: 0=multichannel)
    {
      if( m_nevtProcessedByChannel[channelId] > 0 )
      {
        if( sstr.str() != " {" ) sstr << ",";
        if( channelId == 0 )
          sstr << " no-multichannel";
        else
          sstr << " " << channelId;
        sstr << " : " << m_nevtProcessedByChannel[channelId];
      }
    }
    sstr << " }";
    std::cout << "DEBUG: MEK " << this;
    if( m_tag != "" ) std::cout << " " << m_tag;
    std::cout << " processed " << nevtProcessed << " events across " << CPPProcess::ndiagrams << " channels" << sstr.str() << std::endl;
  }
#endif

  //--------------------------------------------------------------------------

  void MatrixElementKernelBase::dumpSignallingFPEs()
  {
    // New strategy for issue #831: add a final report of FPEs
    // Note: normally only underflow will be reported here (inexact is switched off because it would almost always signal;
    // divbyzero, invalid and overflow are configured by feenablexcept to send a SIGFPE signal, and are normally fixed in the code)
    // Note: this is now called in the individual destructors of MEK classes rather than in that of MatrixElementKernelBase(#837)
    std::string fpes;
    if( std::fetestexcept( FE_DIVBYZERO ) ) fpes += " FE_DIVBYZERO";
    if( std::fetestexcept( FE_INVALID ) ) fpes += " FE_INVALID";
    if( std::fetestexcept( FE_OVERFLOW ) ) fpes += " FE_OVERFLOW";
    if( std::fetestexcept( FE_UNDERFLOW ) ) fpes += " FE_UNDERFLOW";
    //if( std::fetestexcept( FE_INEXACT ) ) fpes += " FE_INEXACT"; // do not print this out: this would almost always signal!
    if( fpes == "" )
      std::cout << "INFO: No Floating Point Exceptions have been reported" << std::endl;
    else
      std::cerr << "INFO: The following Floating Point Exceptions have been reported:" << fpes << std::endl;
  }

  //--------------------------------------------------------------------------
}

//============================================================================

#ifndef MGONGPUCPP_GPUIMPL
namespace mg5amcCpu
{

  //--------------------------------------------------------------------------

  MatrixElementKernelHost::MatrixElementKernelHost( const BufferMomenta& momenta,         // input: momenta
                                                    const BufferGs& gs,                   // input: gs for alphaS
                                                    const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                                                    const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                                                    const BufferChannelIds& channelIds,   // input: channel ids for single-diagram enhancement
                                                    BufferMatrixElements& matrixElements, // output: matrix elements
                                                    BufferSelectedHelicity& selhel,       // output: helicity selection
                                                    BufferSelectedColor& selcol,          // output: color selection
                                                    const size_t nevt )
    : MatrixElementKernelBase( momenta, gs, rndhel, rndcol, channelIds, matrixElements, selhel, selcol )
    , NumberOfEvents( nevt )
    , m_couplings( nevt )
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    , m_numerators( nevt )
    , m_denominators( nevt )
#endif
  {
    //std::cout << "DEBUG: MatrixElementKernelHost::ctor " << this << std::endl;
    if( m_momenta.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelHost: momenta must be a host array" );
    if( m_matrixElements.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelHost: matrixElements must be a host array" );
    if( m_channelIds.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelHost: channelIds must be a device array" );
    if( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "MatrixElementKernelHost: nevt mismatch with momenta" );
    if( this->nevt() != m_matrixElements.nevt() ) throw std::runtime_error( "MatrixElementKernelHost: nevt mismatch with matrixElements" );
    if( this->nevt() != m_channelIds.nevt() ) throw std::runtime_error( "MatrixElementKernelHost: nevt mismatch with channelIds" );
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

  MatrixElementKernelHost::~MatrixElementKernelHost()
  {
    //std::cout << "DEBUG: MatrixElementKernelBase::dtor " << this << std::endl;
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
    // ... 0d2. Copy good helicity list to static memory on the host
    // [FIXME! REMOVE THIS STATIC THAT BREAKS MULTITHREADING?]
    return sigmaKin_setGoodHel( hstIsGoodHel.data() );
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelHost::computeMatrixElements( const bool useChannelIds )
  {
    computeDependentCouplings( m_gs.data(), m_couplings.data(), m_gs.size() );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    const unsigned int* pChannelIds = ( useChannelIds ? m_channelIds.data() : nullptr );
    sigmaKin( m_momenta.data(), m_couplings.data(), m_rndhel.data(), m_rndcol.data(), m_matrixElements.data(), pChannelIds, m_numerators.data(), m_denominators.data(), m_selhel.data(), m_selcol.data(), nevt() );
#else
    assert( useChannelIds == false );
    sigmaKin( m_momenta.data(), m_couplings.data(), m_rndhel.data(), m_matrixElements.data(), m_selhel.data(), nevt() );
#endif
#ifdef MGONGPU_CHANNELID_DEBUG
    //std::cout << "DEBUG: MatrixElementKernelHost::computeMatrixElements " << this << " " << ( useChannelIds ? "T" : "F" ) << " " << nevt() << std::endl;
    MatrixElementKernelBase::updateNevtProcessedByChannel( pChannelIds, nevt() );
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

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
{

  //--------------------------------------------------------------------------

  MatrixElementKernelDevice::MatrixElementKernelDevice( const BufferMomenta& momenta,         // input: momenta
                                                        const BufferGs& gs,                   // input: gs for alphaS
                                                        const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                                                        const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                                                        const BufferChannelIds& channelIds,   // input: channel ids for single-diagram enhancement
                                                        BufferMatrixElements& matrixElements, // output: matrix elements
                                                        BufferSelectedHelicity& selhel,       // output: helicity selection
                                                        BufferSelectedColor& selcol,          // output: color selection
                                                        const size_t gpublocks,
                                                        const size_t gputhreads )
    : MatrixElementKernelBase( momenta, gs, rndhel, rndcol, channelIds, matrixElements, selhel, selcol )
    , NumberOfEvents( gpublocks * gputhreads )
    , m_couplings( this->nevt() )
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    , m_numerators( this->nevt() )
    , m_denominators( this->nevt() )
#endif
    , m_helStreams()
    , m_pHelSelAux()
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    , m_pColSelAux()
#endif
#ifdef MGONGPU_CHANNELID_DEBUG
    , m_hstChannelIds( this->nevt() )
#endif
    , m_gpublocks( gpublocks )
    , m_gputhreads( gputhreads )
  {
    //std::cout << "DEBUG: MatrixElementKernelDevice::ctor " << this << std::endl;
    if( !m_momenta.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelDevice: momenta must be a device array" );
    if( !m_matrixElements.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelDevice: matrixElements must be a device array" );
    if( !m_channelIds.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelDevice: channelIds must be a device array" ); // FIXME?!
    if( m_gpublocks == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gpublocks must be > 0" );
    if( m_gputhreads == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gputhreads must be > 0" );
    if( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "MatrixElementKernelDevice: nevt mismatch with momenta" );
    if( this->nevt() != m_matrixElements.nevt() ) throw std::runtime_error( "MatrixElementKernelDevice: nevt mismatch with matrixElements" );
    if( this->nevt() != m_channelIds.nevt() ) throw std::runtime_error( "MatrixElementKernelDevice: nevt mismatch with channelIds" );
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

  MatrixElementKernelDevice::~MatrixElementKernelDevice()
  {
    //std::cout << "DEBUG: MatrixElementKernelDevice::dtor " << this << std::endl;
    for( int ihel = 0; ihel < CPPProcess::ncomb; ihel++ )
      if( m_helStreams[ihel] ) cudaStreamDestroy( m_helStreams[ihel] ); // do not destroy if nullptr
  }

  //--------------------------------------------------------------------------

  // FIXME! The relevance of this function should be reassessed (#543 and #902)
  void MatrixElementKernelDevice::setGrid( const int /*gpublocks*/, const int /*gputhreads*/ )
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
    // ... 0d1. Compute good helicity mask (a host variable) on the device
    gpuLaunchKernel( computeDependentCouplings, m_gpublocks, m_gputhreads, m_gs.data(), m_couplings.data() );
    const int nevt = m_gpublocks * m_gputhreads;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    sigmaKin_getGoodHel( m_momenta.data(), m_couplings.data(), m_matrixElements.data(), m_numerators.data(), m_denominators.data(), hstIsGoodHel.data(), nevt );
#else
    sigmaKin_getGoodHel( m_momenta.data(), m_couplings.data(), m_matrixElements.data(), hstIsGoodHel.data(), nevt );
#endif
    // ... 0d3. Set good helicity list in host static memory
    int nGoodHel = sigmaKin_setGoodHel( hstIsGoodHel.data() );
    // Create one GPU stream for each good helicity
    for( int ighel = 0; ighel < nGoodHel; ighel++ )
      cudaStreamCreate( &m_helStreams[ighel] );
    // ... Create the auxiliary buffer for helicity selection (for each event: matrix element sum up to each good helicity)
    m_pHelSelAux.reset( new DeviceBufferSimple( nevt * nGoodHel ) );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // ... Create the auxiliary buffer for color selection (for each event: matrix element sum up to each color)
    constexpr int ncolor = CPPProcess::ncolor; // the number of leading colors
    m_pColSelAux.reset( new DeviceBufferSimple( nevt * ncolor ) );
#endif
    // Return the number of good helicities
    return nGoodHel;
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelDevice::computeMatrixElements( const bool useChannelIds )
  {
    gpuLaunchKernel( computeDependentCouplings, m_gpublocks, m_gputhreads, m_gs.data(), m_couplings.data() );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    const unsigned int* pChannelIds = ( useChannelIds ? m_channelIds.data() : nullptr );
    sigmaKin( m_momenta.data(), m_couplings.data(), m_rndhel.data(), m_rndcol.data(), m_matrixElements.data(), pChannelIds, m_numerators.data(), m_denominators.data(), m_selhel.data(), m_selcol.data(), m_gpublocks, m_gputhreads, m_helStreams, m_pHelSelAux->data(), m_pColSelAux->data() );
#else
    assert( useChannelIds == false );
    sigmaKin( m_momenta.data(), m_couplings.data(), m_rndhel.data(), m_matrixElements.data(), m_selhel.data(), m_gpublocks, m_gputhreads, m_helStreams, m_pHelSelAux->data() );
#endif
#ifdef MGONGPU_CHANNELID_DEBUG
    //std::cout << "DEBUG: MatrixElementKernelDevice::computeMatrixElements " << this << " " << ( useChannelIds ? "T" : "F" ) << " " << nevt() << std::endl;
    copyHostFromDevice( m_hstChannelIds, m_channelIds ); // FIXME?!
    const unsigned int* pHstChannelIds = ( useChannelIds ? m_hstChannelIds.data() : nullptr );
    MatrixElementKernelBase::updateNevtProcessedByChannel( pHstChannelIds, nevt() );
#endif
    checkGpu( gpuPeekAtLastError() );   // is this needed?
    checkGpu( gpuDeviceSynchronize() ); // probably not needed? but it avoids errors in sigmaKin above from appearing later on in random places...
  }

  //--------------------------------------------------------------------------

}
#endif

//============================================================================
