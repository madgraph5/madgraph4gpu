#include "BridgeKernels.h"

#include "MemoryAccessMomenta.h"

#include <sstream>

using mgOnGpu::npar; // the number of particles (external = initial + final)
using mgOnGpu::np4;  // the number of dimensions of 4-momenta (E,px,py,pz)

//============================================================================

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

  BridgeKernelBase::BridgeKernelBase( const BufferMomenta& momenta,         // input: momenta
                                      const BufferGs& gs,                   // input: gs for alphaS
                                      const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                                      const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                                      BufferMatrixElements& matrixElements, // output: matrix elements
                                      const size_t nevt )
    : MatrixElementKernelBase( momenta, gs, rndhel, rndcol, matrixElements )
    , NumberOfEvents( nevt )
    , m_bridge( nevt, npar, np4 )
  {
    if( m_momenta.isOnDevice() ) throw std::runtime_error( "BridgeKernelBase: momenta must be a host array" );
    if( m_matrixElements.isOnDevice() ) throw std::runtime_error( "BridgeKernelBase: matrixElements must be a host array" );
    if( this->nevt() != m_momenta.nevt() ) throw std::runtime_error( "BridgeKernelBase: nevt mismatch with momenta" );
    if( this->nevt() != m_matrixElements.nevt() ) throw std::runtime_error( "BridgeKernelBase: nevt mismatch with matrixElements" );
  }

  //--------------------------------------------------------------------------
}

//============================================================================

#ifndef __CUDACC__
namespace mg5amcCpu
{

  //--------------------------------------------------------------------------

  BridgeKernelHost::BridgeKernelHost( const BufferMomenta& momenta,         // input: momenta
                                      const BufferGs& gs,                   // input: Gs for alphaS
                                      const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                                      const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                                      BufferMatrixElements& matrixElements, // output: matrix elements
                                      const size_t nevt )
    : BridgeKernelBase( momenta, gs, rndhel, rndcol, matrixElements, nevt )
    , m_fortranMomenta( nevt )
  {
  }

  //--------------------------------------------------------------------------

  void BridgeKernelHost::transposeInputMomentaC2F()
  {
    hst_transposeMomentaC2F( m_momenta.data(), m_fortranMomenta.data(), nevt() );
  }

  //--------------------------------------------------------------------------

  int BridgeKernelHost::computeGoodHelicities()
  {
    constexpr bool goodHelOnly = true;
    constexpr unsigned int channelId = 0; // disable multi-channel for helicity filtering
    m_bridge.cpu_sequence( m_fortranMomenta.data(), m_gs.data(), m_rndhel.data(), m_rndcol.data(), channelId, m_matrixElements.data(), goodHelOnly );
    return m_bridge.nGoodHel();
  }

  //--------------------------------------------------------------------------

  void BridgeKernelHost::computeMatrixElements( const unsigned int channelId )
  {
    constexpr bool goodHelOnly = false;
    m_bridge.cpu_sequence( m_fortranMomenta.data(), m_gs.data(), m_rndhel.data(), m_rndcol.data(), channelId, m_matrixElements.data(), goodHelOnly );
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
                                          const BufferGs& gs,                   // input: Gs for alphaS
                                          const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                                          const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                                          BufferMatrixElements& matrixElements, // output: matrix elements
                                          const size_t gpublocks,
                                          const size_t gputhreads )
    : BridgeKernelBase( momenta, gs, rndhel, rndcol, matrixElements, gpublocks * gputhreads )
    , m_fortranMomenta( nevt() )
    , m_gpublocks( gpublocks )
    , m_gputhreads( gputhreads )
  {
    if( m_gpublocks == 0 ) throw std::runtime_error( "BridgeKernelDevice: gpublocks must be > 0" );
    if( m_gputhreads == 0 ) throw std::runtime_error( "BridgeKernelDevice: gputhreads must be > 0" );
    m_bridge.set_gpugrid( gpublocks, gputhreads );
  }

  //--------------------------------------------------------------------------

  void BridgeKernelDevice::transposeInputMomentaC2F()
  {
    hst_transposeMomentaC2F( m_momenta.data(), m_fortranMomenta.data(), nevt() );
  }

  //--------------------------------------------------------------------------

  int BridgeKernelDevice::computeGoodHelicities()
  {
    constexpr bool goodHelOnly = true;
    constexpr unsigned int channelId = 0; // disable multi-channel for helicity filtering
    m_bridge.gpu_sequence( m_fortranMomenta.data(), m_gs.data(), m_rndhel.data(), m_rndcol.data(), channelId, m_matrixElements.data(), goodHelOnly );
    return m_bridge.nGoodHel();
  }

  //--------------------------------------------------------------------------

  void BridgeKernelDevice::computeMatrixElements( const unsigned int channelId )
  {
    constexpr bool goodHelOnly = false;
    m_bridge.gpu_sequence( m_fortranMomenta.data(), m_gs.data(), m_rndhel.data(), m_rndcol.data(), channelId, m_matrixElements.data(), goodHelOnly );
  }

  //--------------------------------------------------------------------------

}
#endif

//============================================================================
