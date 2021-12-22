#include "RamboSamplingKernels.h"

#include "checkCuda.h"
#include "MemoryBuffers.h"

#include "rambo.h" // inline implementation of RAMBO algorithms and kernels

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
#ifdef __CUDACC__
    throw std::runtime_error( "RamboSamplingKernelHost is not yet implemented in CUDA" ); // FIXME!
#endif
    if ( m_rnarray.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelHost: rnarray must be a host array" );
    if ( m_momenta.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelHost: momenta must be a host array" );
    if ( m_weights.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelHost: weights must be a host array" );
  }

  //--------------------------------------------------------------------------

  void RamboSamplingKernelHost::getMomentaInitial()
  {
#ifdef __CUDACC__
    throw std::runtime_error( "RamboSamplingKernelHost is not yet implemented in CUDA" ); // FIXME!
#else
    // ** START LOOP ON IEVT **    
    for ( size_t ievt = 0; ievt < nevt(); ++ievt )
    {
      fptype* ievtMomenta = &( MemoryAccessMomenta::ieventAccessIp4Ipar( m_momenta.data(), ievt, 0, 0 ) ); // FIXME: document constraints on functions
      ramboGetMomentaInitial( m_energy, ievtMomenta ); // FIXME!
    }
    // ** END LOOP ON IEVT **    
#endif
  }

  //--------------------------------------------------------------------------

  void RamboSamplingKernelHost::getMomentaFinal()
  {
#ifdef __CUDACC__
    throw std::runtime_error( "RamboSamplingKernelHost is not yet implemented in CUDA" ); // FIXME!
#else
    ramboGetMomentaFinal( m_energy, m_rnarray.data(), m_momenta.data(), m_weights.data(), nevt() );
#endif
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
  }
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void RamboSamplingKernelDevice::getMomentaInitial()
  {
    ramboGetMomentaInitial<<<m_gpublocks, m_gputhreads>>>( m_energy, m_momenta.data() );
  }
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void RamboSamplingKernelDevice::getMomentaFinal()
  {
    ramboGetMomentaFinal<<<m_gpublocks, m_gputhreads>>>( m_energy,
                                                         m_rnarray.data(),
                                                         m_momenta.data(),
                                                         m_weights.data() );
  }
#endif

  //--------------------------------------------------------------------------

}
