#include "RamboSamplingKernels.h"

#include "checkCuda.h"
#include "MemoryAccessMomenta.h"
#include "MemoryAccessRandomNumbers.h"
#include "MemoryAccessWeights.h"
#include "MemoryBuffers.h"
#include "rambo.h" // inline implementation of RAMBO algorithms and kernels

#include <sstream>

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

  RamboSamplingKernelHost::RamboSamplingKernelHost( const fptype energy,                // input: energy
                                                    const BufferRandomNumbers& rnarray, // input: random numbers in [0,1]
                                                    BufferMomenta& momenta,             // output: momenta
                                                    BufferWeights& weights,             // output: weights
                                                    const size_t nevt )
    : SamplingKernelBase( energy, rnarray, momenta, weights )
    , NumberOfEvents( nevt )
  {
    if ( m_rnarray.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelHost: rnarray must be a host array" );
    if ( m_momenta.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelHost: momenta must be a host array" );
    if ( m_weights.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelHost: weights must be a host array" );
    // Sanity checks for memory access (are these really strictly needed?)
    constexpr int neppR = MemoryAccessRandomNumbers::neppR; // AOSOA layout
#ifndef __CUDACC__
    auto ispoweroftwo = []( int n ) { return ( n > 0 ) && !( n & ( n - 1 ) ); }; // https://stackoverflow.com/a/108360
    static_assert( ispoweroftwo( neppR ), "neppR is not a power of 2" );
#else
    static_assert( ( neppR > 0 ) && !( neppR & ( neppR - 1 ) ), "neppR is not a power of 2" ); // implementation without c++17 lambda
#endif
  }

  //--------------------------------------------------------------------------

  void RamboSamplingKernelHost::getMomentaInitial()
  {
    constexpr auto getMomentaInitial = ramboGetMomentaInitial<HostAccessMomenta>;
    // ** START LOOP ON IEVT **
    for ( size_t ievt = 0; ievt < nevt(); ++ievt )
    {
      // NB all KernelLaunchers assume that memory access can be decomposed as "accessField = decodeRecord( accessRecord )"
      fptype* ievtMomenta = MemoryAccessMomenta::ieventAccessRecord( m_momenta.data(), ievt );
      getMomentaInitial( m_energy, ievtMomenta );
    }
    // ** END LOOP ON IEVT **
  }

  //--------------------------------------------------------------------------

  void RamboSamplingKernelHost::getMomentaFinal()
  {
    constexpr auto getMomentaFinal = ramboGetMomentaFinal<HostAccessRandomNumbers, HostAccessMomenta, HostAccessWeights>;
    // ** START LOOP ON IEVT **
    for ( size_t ievt = 0; ievt < nevt(); ++ievt )
    {
      // NB all KernelLaunchers assume that memory access can be decomposed as "accessField = decodeRecord( accessRecord )"
      const fptype* ievtRnarray = MemoryAccessRandomNumbers::ieventAccessRecordConst( m_rnarray.data(), ievt );
      fptype* ievtMomenta = MemoryAccessMomenta::ieventAccessRecord( m_momenta.data(), ievt );
      fptype* ievtWeights = MemoryAccessWeights::ieventAccessRecord( m_weights.data(), ievt );
      getMomentaFinal( m_energy, ievtRnarray, ievtMomenta, ievtWeights );
    }
    // ** END LOOP ON IEVT **
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  RamboSamplingKernelDevice::RamboSamplingKernelDevice( const fptype energy,                // input: energy
                                                        const BufferRandomNumbers& rnarray, // input: random numbers in [0,1]
                                                        BufferMomenta& momenta,             // output: momenta
                                                        BufferWeights& weights,             // output: weights
                                                        const size_t gpublocks,
                                                        const size_t gputhreads )
    : SamplingKernelBase( energy, rnarray, momenta, weights )
    , m_gpublocks( gpublocks )
    , m_gputhreads( gputhreads )
  {
    if ( ! m_rnarray.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelDevice: rnarray must be a device array" );
    if ( ! m_momenta.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelDevice: momenta must be a device array" );
    if ( ! m_weights.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelDevice: weights must be a device array" );
    if ( m_gpublocks == 0 ) throw std::runtime_error( "RamboSamplingKernelDevice: gpublocks must be > 0" );
    if ( m_gputhreads == 0 ) throw std::runtime_error( "RamboSamplingKernelDevice: gputhreads must be > 0" );
    // Sanity checks for memory access (are these really strictly needed?)
    constexpr int neppR = MemoryAccessRandomNumbers::neppR; // AOSOA layout
    static_assert( ( neppR > 0 ) && !( neppR & ( neppR - 1 ) ), "neppR is not a power of 2" ); // implementation without c++17 lambda
    if ( m_gputhreads%neppR != 0 )
    {
      std::ostringstream sstr;
      sstr << "RamboSamplingKernelDevice: gputhreads should be a multiple of neppR=" << neppR;
      throw std::runtime_error( sstr.str() );
    }
  }
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __global__
  void getMomentaInitialDevice( const fptype energy,
                                fptype* momenta )
  {
    constexpr auto getMomentaInitial = ramboGetMomentaInitial<DeviceAccessMomenta>;
    return getMomentaInitial( energy, momenta );
  }
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void RamboSamplingKernelDevice::getMomentaInitial()
  {
    getMomentaInitialDevice<<<m_gpublocks, m_gputhreads>>>( m_energy, m_momenta.data() );
  }
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __global__
  void getMomentaFinalDevice( const fptype energy,
                              const fptype* rnarray,
                              fptype* momenta,
                              fptype* wgts )
  {
    constexpr auto getMomentaFinal = ramboGetMomentaFinal<DeviceAccessRandomNumbers, DeviceAccessMomenta, DeviceAccessWeights>;
    return getMomentaFinal( energy, rnarray, momenta, wgts );
  }
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void RamboSamplingKernelDevice::getMomentaFinal()
  {
    getMomentaFinalDevice<<<m_gpublocks, m_gputhreads>>>( m_energy, m_rnarray.data(), m_momenta.data(), m_weights.data() );
  }
#endif

  //--------------------------------------------------------------------------

}
