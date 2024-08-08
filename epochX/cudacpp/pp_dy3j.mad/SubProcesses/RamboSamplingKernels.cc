// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: J. Teig, A. Valassi (2021-2024) for the MG5aMC CUDACPP plugin.

#include "RamboSamplingKernels.h"

#include "GpuRuntime.h"
#include "MemoryAccessMomenta.h"
#include "MemoryAccessRandomNumbers.h"
#include "MemoryAccessWeights.h"
#include "MemoryBuffers.h"
#include "rambo.h" // inline implementation of RAMBO algorithms and kernels

#include <sstream>

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

  RamboSamplingKernelHost::RamboSamplingKernelHost( const fptype energy,               // input: energy
                                                    const BufferRndNumMomenta& rndmom, // input: random numbers in [0,1]
                                                    BufferMomenta& momenta,            // output: momenta
                                                    BufferWeights& weights,            // output: weights
                                                    const size_t nevt )
    : SamplingKernelBase( energy, rndmom, momenta, weights )
    , NumberOfEvents( nevt )
  {
    if( m_rndmom.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelHost: rndmom must be a host array" );
    if( m_momenta.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelHost: momenta must be a host array" );
    if( m_weights.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelHost: weights must be a host array" );
    if( this->nevt() != m_rndmom.nevt() ) throw std::runtime_error( "RamboSamplingKernelHost: nevt mismatch with rndmom" );
    if( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "RamboSamplingKernelHost: nevt mismatch with momenta" );
    if( this->nevt() != m_weights.nevt() ) throw std::runtime_error( "RamboSamplingKernelHost: nevt mismatch with weights" );
    // Sanity checks for memory access (momenta buffer)
    constexpr int neppM = MemoryAccessMomenta::neppM; // AOSOA layout
    static_assert( ispoweroftwo( neppM ), "neppM is not a power of 2" );
    if( nevt % neppM != 0 )
    {
      std::ostringstream sstr;
      sstr << "RamboSamplingKernelHost: nevt should be a multiple of neppM=" << neppM;
      throw std::runtime_error( sstr.str() );
    }
    // Sanity checks for memory access (random number buffer)
    constexpr int neppR = MemoryAccessRandomNumbers::neppR; // AOSOA layout
    static_assert( ispoweroftwo( neppR ), "neppR is not a power of 2" );
    if( nevt % neppR != 0 )
    {
      std::ostringstream sstr;
      sstr << "RamboSamplingKernelHost: nevt should be a multiple of neppR=" << neppR;
      throw std::runtime_error( sstr.str() );
    }
  }

  //--------------------------------------------------------------------------

  void
  RamboSamplingKernelHost::getMomentaInitial()
  {
    constexpr auto getMomentaInitial = ramboGetMomentaInitial<HostAccessMomenta>;
    // ** START LOOP ON IEVT **
    for( size_t ievt = 0; ievt < nevt(); ++ievt )
    {
      // NB all KernelLaunchers assume that memory access can be decomposed as "accessField = decodeRecord( accessRecord )"
      fptype* ievtMomenta = MemoryAccessMomenta::ieventAccessRecord( m_momenta.data(), ievt );
      getMomentaInitial( m_energy, ievtMomenta );
    }
    // ** END LOOP ON IEVT **
  }

  //--------------------------------------------------------------------------

  void
  RamboSamplingKernelHost::getMomentaFinal()
  {
    constexpr auto getMomentaFinal = ramboGetMomentaFinal<HostAccessRandomNumbers, HostAccessMomenta, HostAccessWeights>;
    // ** START LOOP ON IEVT **
    for( size_t ievt = 0; ievt < nevt(); ++ievt )
    {
      // NB all KernelLaunchers assume that memory access can be decomposed as "accessField = decodeRecord( accessRecord )"
      const fptype* ievtRndmom = MemoryAccessRandomNumbers::ieventAccessRecordConst( m_rndmom.data(), ievt );
      fptype* ievtMomenta = MemoryAccessMomenta::ieventAccessRecord( m_momenta.data(), ievt );
      fptype* ievtWeights = MemoryAccessWeights::ieventAccessRecord( m_weights.data(), ievt );
      getMomentaFinal( m_energy, ievtRndmom, ievtMomenta, ievtWeights );
    }
    // ** END LOOP ON IEVT **
  }

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  RamboSamplingKernelDevice::RamboSamplingKernelDevice( const fptype energy,               // input: energy
                                                        const BufferRndNumMomenta& rndmom, // input: random numbers in [0,1]
                                                        BufferMomenta& momenta,            // output: momenta
                                                        BufferWeights& weights,            // output: weights
                                                        const size_t gpublocks,
                                                        const size_t gputhreads )
    : SamplingKernelBase( energy, rndmom, momenta, weights )
    , NumberOfEvents( gpublocks * gputhreads )
    , m_gpublocks( gpublocks )
    , m_gputhreads( gputhreads )
  {
    if( !m_rndmom.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelDevice: rndmom must be a device array" );
    if( !m_momenta.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelDevice: momenta must be a device array" );
    if( !m_weights.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelDevice: weights must be a device array" );
    if( m_gpublocks == 0 ) throw std::runtime_error( "RamboSamplingKernelDevice: gpublocks must be > 0" );
    if( m_gputhreads == 0 ) throw std::runtime_error( "RamboSamplingKernelDevice: gputhreads must be > 0" );
    if( this->nevt() != m_rndmom.nevt() ) throw std::runtime_error( "RamboSamplingKernelDevice: nevt mismatch with rndmom" );
    if( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "RamboSamplingKernelDevice: nevt mismatch with momenta" );
    if( this->nevt() != m_weights.nevt() ) throw std::runtime_error( "RamboSamplingKernelDevice: nevt mismatch with weights" );
    // Sanity checks for memory access (momenta buffer)
    constexpr int neppM = MemoryAccessMomenta::neppM; // AOSOA layout
    static_assert( ispoweroftwo( neppM ), "neppM is not a power of 2" );
    if( m_gputhreads % neppM != 0 )
    {
      std::ostringstream sstr;
      sstr << "RamboSamplingKernelHost: gputhreads should be a multiple of neppM=" << neppM;
      throw std::runtime_error( sstr.str() );
    }
    // Sanity checks for memory access (random number buffer)
    constexpr int neppR = MemoryAccessRandomNumbers::neppR; // AOSOA layout
    static_assert( ispoweroftwo( neppR ), "neppR is not a power of 2" );
    if( m_gputhreads % neppR != 0 )
    {
      std::ostringstream sstr;
      sstr << "RamboSamplingKernelDevice: gputhreads should be a multiple of neppR=" << neppR;
      throw std::runtime_error( sstr.str() );
    }
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  __global__ void
  getMomentaInitialDevice( const fptype energy,
                           fptype* momenta )
  {
    constexpr auto getMomentaInitial = ramboGetMomentaInitial<DeviceAccessMomenta>;
    return getMomentaInitial( energy, momenta );
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  void
  RamboSamplingKernelDevice::getMomentaInitial()
  {
    gpuLaunchKernel( getMomentaInitialDevice, m_gpublocks, m_gputhreads, m_energy, m_momenta.data() );
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  __global__ void
  getMomentaFinalDevice( const fptype energy,
                         const fptype* rndmom,
                         fptype* momenta,
                         fptype* wgts )
  {
    constexpr auto getMomentaFinal = ramboGetMomentaFinal<DeviceAccessRandomNumbers, DeviceAccessMomenta, DeviceAccessWeights>;
    return getMomentaFinal( energy, rndmom, momenta, wgts );
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  void
  RamboSamplingKernelDevice::getMomentaFinal()
  {
    gpuLaunchKernel( getMomentaFinalDevice, m_gpublocks, m_gputhreads, m_energy, m_rndmom.data(), m_momenta.data(), m_weights.data() );
  }
#endif

  //--------------------------------------------------------------------------
}
