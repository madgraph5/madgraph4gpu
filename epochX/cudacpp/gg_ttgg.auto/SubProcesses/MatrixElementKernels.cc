#include "MatrixElementKernels.h"

#include "checkCuda.h"
#include "CPPProcess.h"
#include "MemoryAccessMomenta.h"
#include "MemoryBuffers.h"

#include <sstream>

//============================================================================

#ifndef __CUDACC__
namespace mg5amcCpu
{  

  //--------------------------------------------------------------------------

  MatrixElementKernelHost::MatrixElementKernelHost( const BufferMomenta& momenta,         // input: momenta
                                                    BufferMatrixElements& matrixElements, // output: matrix elements
                                                    const size_t nevt )
    : MatrixElementKernelBase( momenta, matrixElements )
    , NumberOfEvents( nevt )
  {
    if ( m_momenta.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelHost: momenta must be a host array" );
    if ( m_matrixElements.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelHost: matrixElements must be a host array" );
    if ( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "MatrixElementKernelHost: nevt mismatch with momenta" );
    if ( this->nevt() != m_matrixElements.nevt() ) throw std::runtime_error( "MatrixElementKernelHost: nevt mismatch with matrixElements" );
    // Sanity checks for memory access (momenta buffer)
    constexpr int neppM = MemoryAccessMomenta::neppM; // AOSOA layout
    static_assert( ispoweroftwo( neppM ), "neppM is not a power of 2" );
    if ( nevt%neppM != 0 )
    {
      std::ostringstream sstr;
      sstr << "MatrixElementKernelHost: nevt should be a multiple of neppM=" << neppM;
      throw std::runtime_error( sstr.str() );
    }
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelHost::computeGoodHelicities()
  {
    using mgOnGpu::ncomb; // the number of helicity combinations
    HostBufferHelicityMask hstIsGoodHel( ncomb );
    // ... 0d1. Compute good helicity mask on the host
    sigmaKin_getGoodHel( m_momenta.data(), m_matrixElements.data(), hstIsGoodHel.data(), nevt() );
    // ... 0d2. Copy back good helicity list to static memory on the host
    // [FIXME! REMOVE THIS STATIC THAT BREAKS MULTITHREADING?]
    sigmaKin_setGoodHel( hstIsGoodHel.data() );
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelHost::computeMatrixElements()
  {
    sigmaKin( m_momenta.data(), m_matrixElements.data(), nevt() );
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
                                                        BufferMatrixElements& matrixElements, // output: matrix elements
                                                        const size_t gpublocks,
                                                        const size_t gputhreads )
    : MatrixElementKernelBase( momenta, matrixElements )
    , NumberOfEvents( gpublocks*gputhreads )
    , m_gpublocks( gpublocks )
    , m_gputhreads( gputhreads )
  {
    if ( ! m_momenta.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelDevice: momenta must be a device array" );
    if ( ! m_matrixElements.isOnDevice() ) throw std::runtime_error( "MatrixElementKernelDevice: matrixElements must be a device array" );
    if ( m_gpublocks == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gpublocks must be > 0" );
    if ( m_gputhreads == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gputhreads must be > 0" );
    if ( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "MatrixElementKernelDevice: nevt mismatch with momenta" );
    if ( this->nevt() != m_matrixElements.nevt() ) throw std::runtime_error( "MatrixElementKernelDevice: nevt mismatch with matrixElements" );
    // Sanity checks for memory access (momenta buffer)
    constexpr int neppM = MemoryAccessMomenta::neppM; // AOSOA layout
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
    if ( m_gpublocks == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gpublocks must be > 0 in setGrid" );
    if ( m_gputhreads == 0 ) throw std::runtime_error( "MatrixElementKernelDevice: gputhreads must be > 0 in setGrid" );
    if ( this->nevt() != m_gpublocks * m_gputhreads ) throw std::runtime_error( "MatrixElementKernelDevice: nevt mismatch in setGrid" );
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelDevice::computeGoodHelicities()
  {
    using mgOnGpu::ncomb; // the number of helicity combinations
    PinnedHostBufferHelicityMask hstIsGoodHel( ncomb );
    DeviceBufferHelicityMask devIsGoodHel( ncomb );
    // ... 0d1. Compute good helicity mask on the device
    sigmaKin_getGoodHel<<<m_gpublocks, m_gputhreads>>>( m_momenta.data(), m_matrixElements.data(), devIsGoodHel.data() );
    checkCuda( cudaPeekAtLastError() );
    // ... 0d2. Copy back good helicity mask to the host
    copyHostFromDevice( hstIsGoodHel, devIsGoodHel );
    // ... 0d3. Copy back good helicity list to constant memory on the device
    sigmaKin_setGoodHel( hstIsGoodHel.data() );
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelDevice::computeMatrixElements()
  {
#ifndef MGONGPU_NSIGHT_DEBUG
    sigmaKin<<<m_gpublocks, m_gputhreads>>>( m_momenta.data(), m_matrixElements.data() );
#else
    sigmaKin<<<m_gpublocks, m_gputhreads, ntpbMAX*sizeof(float)>>>( m_momenta.data(), m_matrixElements.data() );
#endif
    checkCuda( cudaPeekAtLastError() );
    checkCuda( cudaDeviceSynchronize() );
  }

  //--------------------------------------------------------------------------

}
#endif

//============================================================================
