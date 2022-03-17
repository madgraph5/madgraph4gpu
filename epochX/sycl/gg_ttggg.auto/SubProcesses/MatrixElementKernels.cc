#include "MatrixElementKernels.h"

#include "CPPProcess.h"
#include "MemoryBuffers.h"

#include <sstream>

//============================================================================

namespace mg5amcGpu
{

  //--------------------------------------------------------------------------

  MatrixElementKernelDevice::MatrixElementKernelDevice( const BufferMomenta& momenta,         // input: momenta
                                                        BufferMatrixElements& matrixElements, // output: matrix elements
                                                        sycl::queue q,
                                                        const size_t gpublocks,
                                                        const size_t gputhreads )
    : MatrixElementKernelBase( momenta, matrixElements, q )
    , NumberOfEvents( gpublocks*gputhreads )
    , m_gpublocks( gpublocks )
    , m_gputhreads( gputhreads )
    , m_cHel( mgOnGpu::ncomb*mgOnGpu::npar, q )
    , m_cIPC( mgOnGpu::ncouplingstimes2, q )
    , m_cIPD( mgOnGpu::nparams, q )
    , m_cNGoodHel( 1, q )
    , m_cGoodHel( mgOnGpu::ncomb, q )
  {
    if ( ! m_momenta.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelDevice: momenta must be a device array" );
    if ( ! m_matrixElements.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelDevice: matrixElements must be a device array" );
    if ( m_gpublocks == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gpublocks must be > 0" );
    if ( m_gputhreads == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gputhreads must be > 0" );
    if ( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "MatrixElementKernelDevice: nevt mismatch with momenta" );
    if ( this->nevt() != m_matrixElements.nevt() ) throw std::runtime_error( "MatrixElementKernelDevice: nevt mismatch with matrixElements" );
    // Sanity checks for memory access (momenta buffer)
    constexpr int neppM = mgOnGpu::neppM; // AOSOA layout
    static_assert( ispoweroftwo( neppM ), "neppM is not a power of 2" );
    if ( m_gputhreads%neppM != 0 )
    {
      std::ostringstream sstr;
      sstr << "MatrixElementKernelHost: gputhreads should be a multiple of neppM=" << neppM;
      throw std::runtime_error( sstr.str() );
    }

  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelDevice::setGrid( const int gpublocks, const int gputhreads )
  {
    m_gpublocks = gpublocks;
    m_gputhreads = gputhreads;
    if ( m_gpublocks == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gpublocks must be > 0 in setGrid" );
    if ( m_gputhreads == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gputhreads must be > 0 in setGrid" );
    if ( this->nevt() != m_gpublocks * m_gputhreads ) throw std::runtime_error( "MatrixElementKernelDevice: nevt mismatch in setGrid" );
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelDevice::setDeviceArrays( const short* tHel, const cxtype* tIPC, const fptype* tIPD )
  {
    m_q.memcpy( m_cHel.data(), tHel, m_cHel.bytes() ).wait();
    m_q.memcpy( m_cIPC.data(), tIPC, m_cIPC.bytes() ).wait();
    m_q.memcpy( m_cIPD.data(), tIPD, m_cIPD.bytes() ).wait();
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelDevice::computeGoodHelicities()
  {
    PinnedHostBufferHelicityMask l_hstIsGoodHel( mgOnGpu::ncomb, m_q );
    DeviceBufferHelicityMask l_devIsGoodHel( mgOnGpu::ncomb, m_q );
    auto momenta = m_momenta.data();
    auto matrixElements = m_matrixElements.data();
    auto cHel = m_cHel.data();
    auto cIPC = m_cIPC.data();
    auto cIPD = m_cIPD.data();
    auto devIsGoodHel = l_devIsGoodHel.data();
    auto hstIsGoodHel = l_hstIsGoodHel.data();

    // ... 0d1. Compute good helicity mask on the device
    m_q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for_work_group(sycl::range<1>{m_gpublocks}, sycl::range<1>{m_gputhreads}, ([=](sycl::group<1> wGroup) {
            wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                size_t ievt = index.get_global_id(0);
                sigmaKin_getGoodHel( momenta, matrixElements, devIsGoodHel, ievt, cHel, cIPC, cIPD );
            });
        }));
    });
    m_q.wait();

    // ... 0d2. Copy back good helicity mask to the host
    copyHostFromDevice( l_hstIsGoodHel, l_devIsGoodHel );
    // ... 0d3. Copy back good helicity list to constant memory on the device
    int goodHel[mgOnGpu::ncomb] = {0}; // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes
    int nGoodHel = sigmaKin_setGoodHel( hstIsGoodHel, goodHel );

    m_q.memcpy( m_cNGoodHel.data(), &nGoodHel, m_cNGoodHel.bytes() ).wait();
    m_q.memcpy( m_cGoodHel.data(), goodHel, m_cGoodHel.bytes() ).wait();
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelDevice::computeMatrixElements()
  {
    auto momenta = m_momenta.data();
    auto matrixElements = m_matrixElements.data();
    auto cHel = m_cHel.data();
    auto cIPC = m_cIPC.data();
    auto cIPD = m_cIPD.data();
    auto cNGoodHel = m_cNGoodHel.data();
    auto cGoodHel = m_cGoodHel.data();

    m_q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for_work_group(sycl::range<1>{m_gpublocks}, sycl::range<1>{m_gputhreads}, ([=](sycl::group<1> wGroup) {
            wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                size_t ievt = index.get_global_id(0);
                sigmaKin( momenta, matrixElements, ievt, cHel, cIPC, cIPD, cNGoodHel, cGoodHel );
            });
        }));
    });
    m_q.wait();
  }

  //--------------------------------------------------------------------------

}

//============================================================================
