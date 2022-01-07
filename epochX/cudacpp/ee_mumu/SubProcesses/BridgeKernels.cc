#include "BridgeKernels.h"
#include "MemoryAccessMomenta.h"

#include <sstream>

using mgOnGpu::npar;  // the number of particles (external = initial + final)
using mgOnGpu::np4;   // the number of dimensions of 4-momenta (E,px,py,pz)
using mgOnGpu::ncomb; // the number of helicity combinations

//============================================================================

#ifndef __CUDACC__
namespace mg5amcCpu
{  

  //--------------------------------------------------------------------------

  BridgeKernelHost::BridgeKernelHost( const BufferMomenta& momenta,         // input: momenta
                                      BufferMatrixElements& matrixElements, // output: matrix elements
                                      const size_t nevt )
    : MatrixElementKernelBase( momenta, matrixElements )
    , NumberOfEvents( nevt )
    , m_bridge( nevt, npar, np4, MemoryAccessMomenta::neppM, ncomb )
  {
    if ( m_momenta.isOnDevice() ) throw std::runtime_error( "BridgeKernelHost: momenta must be a host array" );
    if ( m_matrixElements.isOnDevice() ) throw std::runtime_error( "BridgeKernelHost: matrixElements must be a host array" );
    if ( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "BridgeKernelHost: nevt mismatch with momenta" );
    if ( this->nevt() != m_matrixElements.nevt() ) throw std::runtime_error( "BridgeKernelHost: nevt mismatch with matrixElements" );
  }

  //--------------------------------------------------------------------------

  void BridgeKernelHost::computeGoodHelicities()
  {
    constexpr bool goodHelOnly=true;
    m_bridge.cpu_sequence( m_momenta.data(), m_matrixElements.data(), goodHelOnly );
  }

  //--------------------------------------------------------------------------

  void BridgeKernelHost::computeMatrixElements()
  {
    constexpr bool goodHelOnly=false;
    m_bridge.cpu_sequence( m_momenta.data(), m_matrixElements.data(), goodHelOnly );
  }

  //--------------------------------------------------------------------------

}
#endif

//============================================================================

#ifdef __CUDACC__
namespace mg5amcGpu
{

  //--------------------------------------------------------------------------

  BridgeKernelDevice::BridgeKernelDevice( const BufferMomenta& momenta,         // input: momenta
                                          BufferMatrixElements& matrixElements, // output: matrix elements
                                          const size_t gpublocks,
                                          const size_t gputhreads )
    : MatrixElementKernelBase( momenta, matrixElements )
    , NumberOfEvents( gpublocks*gputhreads )
    , m_bridge( nevt(), npar, np4, MemoryAccessMomenta::neppM, ncomb )
    , m_gpublocks( gpublocks )
    , m_gputhreads( gputhreads )
  {
    if ( m_momenta.isOnDevice() ) throw std::runtime_error( "BridgeKernelDevice: momenta must be a host array" );
    if ( m_matrixElements.isOnDevice() ) throw std::runtime_error( "BridgeKernelDevice: matrixElements must be a host array" );
    if ( m_gpublocks == 0 ) throw std::runtime_error( "BridgeKernelDevice: gpublocks must be > 0" );
    if ( m_gputhreads == 0 ) throw std::runtime_error( "BridgeKernelDevice: gputhreads must be > 0" );
    if ( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "BridgeKernelDevice: nevt mismatch with momenta" );
    if ( this->nevt() != m_matrixElements.nevt() ) throw std::runtime_error( "BridgeKernelDevice: nevt mismatch with matrixElements" );
  }

  //--------------------------------------------------------------------------

  void BridgeKernelDevice::computeGoodHelicities()
  {
    constexpr bool goodHelOnly=true;
    m_bridge.gpu_sequence( m_momenta.data(), m_matrixElements.data(), goodHelOnly );
  }

  //--------------------------------------------------------------------------

  void BridgeKernelDevice::computeMatrixElements()
  {
    constexpr bool goodHelOnly=false;
    m_bridge.gpu_sequence( m_momenta.data(), m_matrixElements.data(), goodHelOnly );
  }

  //--------------------------------------------------------------------------

}
#endif

//============================================================================
