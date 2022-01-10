#include "CrossSectionKernels.h"

#include "MemoryAccessMatrixElements.h"
#include "MemoryAccessWeights.h"
#include "MemoryBuffers.h"

#include <sstream>

//============================================================================

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  //--------------------------------------------------------------------------

  CrossSectionKernelHost::CrossSectionKernelHost( const BufferWeights& samplingWeights,       // input: sampling weights
                                                  const BufferMatrixElements& matrixElements, // input: matrix elements
                                                  EventStatistics& stats,                     // output: event statistics
                                                  const size_t nevt )
    : CrossSectionKernelBase( samplingWeights, matrixElements, stats )
    , NumberOfEvents( nevt )
  {
    if ( m_samplingWeights.isOnDevice() ) throw std::runtime_error( "CrossSectionKernelHost: samplingWeights must be a host array" );
    if ( m_matrixElements.isOnDevice() ) throw std::runtime_error( "CrossSectionKernelHost: matrixElements must be a host array" );
    if ( this->nevt() != m_samplingWeights.nevt() ) throw std::runtime_error( "CrossSectionKernelHost: nevt mismatch with samplingWeights" );
    if ( this->nevt() != m_matrixElements.nevt() ) throw std::runtime_error( "CrossSectionKernelHost: nevt mismatch with matrixElements" );
  }

  //--------------------------------------------------------------------------

  void CrossSectionKernelHost::updateEventStatistics()
  {
  }

  //--------------------------------------------------------------------------

}

//============================================================================

#ifdef __CUDACC__
namespace mg5amcGpu
{

  //--------------------------------------------------------------------------

  CrossSectionKernelDevice::CrossSectionKernelDevice( const BufferWeights& samplingWeights,       // input: sampling weights
                                                      const BufferMatrixElements& matrixElements, // input: matrix elements
                                                      EventStatistics& stats,                     // output: event statistics
                                                      const size_t gpublocks,
                                                      const size_t gputhreads )
    : CrossSectionKernelBase( samplingWeights, matrixElements, stats )
    , NumberOfEvents( gpublocks*gputhreads )
    , m_gpublocks( gpublocks )
    , m_gputhreads( gputhreads )
  {
    if ( ! m_samplingWeights.isOnDevice() ) throw std::runtime_error( "CrossSectionKernelDevice: samplingWeights must be a device array" );
    if ( ! m_matrixElements.isOnDevice() ) throw std::runtime_error( "CrossSectionKernelDevice: matrixElements must be a device array" );
    if ( m_gpublocks == 0 ) throw std::runtime_error( "CrossSectionKernelDevice: gpublocks must be > 0" );
    if ( m_gputhreads == 0 ) throw std::runtime_error( "CrossSectionKernelDevice: gputhreads must be > 0" );
    if ( this->nevt() != m_samplingWeights.nevt() ) throw std::runtime_error( "CrossSectionKernelDevice: nevt mismatch with samplingWeights" );
    if ( this->nevt() != m_matrixElements.nevt() ) throw std::runtime_error( "CrossSectionKernelDevice: nevt mismatch with matrixElements" );
  }

  //--------------------------------------------------------------------------

  void CrossSectionKernelDevice::setGrid( const int gpublocks, const int gputhreads )
  {
    if ( m_gpublocks == 0 ) throw std::runtime_error( "CrossSectionKernelDevice: gpublocks must be > 0 in setGrid" );
    if ( m_gputhreads == 0 ) throw std::runtime_error( "CrossSectionKernelDevice: gputhreads must be > 0 in setGrid" );
    if ( this->nevt() != m_gpublocks * m_gputhreads ) throw std::runtime_error( "CrossSectionKernelDevice: nevt mismatch in setGrid" );
  }

  //--------------------------------------------------------------------------

  void CrossSectionKernelDevice::updateEventStatistics()
  {
  }

  //--------------------------------------------------------------------------

}
#endif

//============================================================================
