#include "MatrixElementKernels.h"

#include "checkCuda.h"
#include "MemoryAccessMomenta.h"
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
  }

  //--------------------------------------------------------------------------

  void MatrixElementKernelHost::computeMatrixElements()
  {
    // ** START LOOP ON IEVT **
    for ( size_t ievt = 0; ievt < nevt(); ++ievt )
    {
    }
    // ** END LOOP ON IEVT **
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
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
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void MatrixElementKernelDevice::computeGoodHelicities()
  {
  }
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void MatrixElementKernelDevice::computeMatrixElements()
  {
  }
#endif

  //--------------------------------------------------------------------------

}
