#include "RamboSamplingKernels.h"

#include "MemoryBuffers.h"
#include "rambo.h" // inline implementation of RAMBO algorithms and kernels

#include <sstream>

namespace mg5amcGpu
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
    if ( this->nevt() != m_rnarray.nevt() ) throw std::runtime_error( "RamboSamplingKernelHost: nevt mismatch with rnarray" );
    if ( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "RamboSamplingKernelHost: nevt mismatch with momenta" );
    if ( this->nevt() != m_weights.nevt() ) throw std::runtime_error( "RamboSamplingKernelHost: nevt mismatch with weights" );
    // Sanity checks for memory access (momenta buffer)
    constexpr int neppM = mgOnGpu::neppM; // AOSOA layout
    static_assert( ispoweroftwo( neppM ), "neppM is not a power of 2" );
    if ( nevt%neppM != 0 )
    {
      std::ostringstream sstr;
      sstr << "RamboSamplingKernelHost: nevt should be a multiple of neppM=" << neppM;
      throw std::runtime_error( sstr.str() );
    }
    // Sanity checks for memory access (random number buffer)
    constexpr int neppR = 8; // FIXME HARDCODED for consistant physics
    static_assert( ispoweroftwo( neppR ), "neppR is not a power of 2" );
    if ( nevt%neppR != 0 )
    {
      std::ostringstream sstr;
      sstr << "RamboSamplingKernelHost: nevt should be a multiple of neppR=" << neppR;
      throw std::runtime_error( sstr.str() );
    }
  }

  //--------------------------------------------------------------------------

  void RamboSamplingKernelHost::getMomentaInitial()
  {
    constexpr auto getMomentaInitial = ramboGetMomentaInitial;
    // ** START LOOP ON IEVT **
    for ( size_t ievt = 0; ievt < nevt(); ++ievt )
    {
      getMomentaInitial( m_energy, m_momenta.data(), ievt );
    }
    // ** END LOOP ON IEVT **
  }

  //--------------------------------------------------------------------------

  void RamboSamplingKernelHost::getMomentaFinal()
  {
    constexpr auto getMomentaFinal = ramboGetMomentaFinal;
    // ** START LOOP ON IEVT **
    for ( size_t ievt = 0; ievt < nevt(); ++ievt )
    {
      getMomentaFinal( m_energy, m_rnarray.data(), m_momenta.data(), m_weights.data(), ievt );
    }
    // ** END LOOP ON IEVT **
  }

  //--------------------------------------------------------------------------

  RamboSamplingKernelDevice::RamboSamplingKernelDevice( const fptype energy,                // input: energy
                                                        const BufferRandomNumbers& rnarray, // input: random numbers in [0,1]
                                                        BufferMomenta& momenta,             // output: momenta
                                                        BufferWeights& weights,             // output: weights
                                                        sycl::queue q,
                                                        const size_t gpublocks,
                                                        const size_t gputhreads )
    : SamplingKernelBase( energy, rnarray, momenta, weights )
    , NumberOfEvents( gpublocks*gputhreads )
    , m_gpublocks( gpublocks )
    , m_gputhreads( gputhreads )
    , m_q( q )
  {
    if ( ! m_rnarray.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelDevice: rnarray must be a device array" );
    if ( ! m_momenta.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelDevice: momenta must be a device array" );
    if ( ! m_weights.isOnDevice() ) throw std::runtime_error( "RamboSamplingKernelDevice: weights must be a device array" );
    if ( m_gpublocks == 0 ) throw std::runtime_error( "RamboSamplingKernelDevice: gpublocks must be > 0" );
    if ( m_gputhreads == 0 ) throw std::runtime_error( "RamboSamplingKernelDevice: gputhreads must be > 0" );
    if ( this->nevt() != m_rnarray.nevt() ) throw std::runtime_error( "RamboSamplingKernelDevice: nevt mismatch with rnarray" );
    if ( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "RamboSamplingKernelDevice: nevt mismatch with momenta" );
    if ( this->nevt() != m_weights.nevt() ) throw std::runtime_error( "RamboSamplingKernelDevice: nevt mismatch with weights" );
    // Sanity checks for memory access (momenta buffer)
    constexpr int neppM = mgOnGpu::neppM; // AOSOA layout
    static_assert( ispoweroftwo( neppM ), "neppM is not a power of 2" );
    if ( m_gputhreads%neppM != 0 )
    {
      std::ostringstream sstr;
      sstr << "RamboSamplingKernelHost: gputhreads should be a multiple of neppM=" << neppM;
      throw std::runtime_error( sstr.str() );
    }
    // Sanity checks for memory access (random number buffer)
    constexpr int neppR = 8; // FIXME HARDCODED
    static_assert( ispoweroftwo( neppR ), "neppR is not a power of 2" );
    if ( m_gputhreads%neppR != 0 )
    {
      std::ostringstream sstr;
      sstr << "RamboSamplingKernelDevice: gputhreads should be a multiple of neppR=" << neppR;
      throw std::runtime_error( sstr.str() );
    }
  }

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void getMomentaInitialDevice( const fptype energy,
                                fptype* momenta,
                                const size_t ievt )
  {
    return ramboGetMomentaInitial( energy, momenta, ievt );
  }

  //--------------------------------------------------------------------------

  void RamboSamplingKernelDevice::getMomentaInitial()
  {
    auto energy = m_energy;
    auto momenta = m_momenta.data();
    m_q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for_work_group(sycl::range<1>{m_gpublocks}, sycl::range<1>{m_gputhreads}, ([=](sycl::group<1> wGroup) {
            wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                size_t ievt = index.get_global_id(0);
                getMomentaInitialDevice( energy, momenta, ievt );
            });
        }));
    });
    m_q.wait();
  }

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void getMomentaFinalDevice( const fptype energy,
                              const fptype* rnarray,
                              fptype* momenta,
                              fptype* wgts,
                              const size_t ievt )
  {
    return ramboGetMomentaFinal( energy, rnarray, momenta, wgts, ievt );
  }

  //--------------------------------------------------------------------------

  void RamboSamplingKernelDevice::getMomentaFinal()
  {
    auto energy = m_energy;
    auto rnarray = m_rnarray.data();
    auto momenta = m_momenta.data();
    auto weights = m_weights.data();
    m_q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for_work_group(sycl::range<1>{m_gpublocks}, sycl::range<1>{m_gputhreads}, ([=](sycl::group<1> wGroup) {
            wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                size_t ievt = index.get_global_id(0);
                getMomentaFinalDevice( energy, rnarray, momenta, weights, ievt );
            });
        }));
    });
    m_q.wait();
  }

  //--------------------------------------------------------------------------

}
