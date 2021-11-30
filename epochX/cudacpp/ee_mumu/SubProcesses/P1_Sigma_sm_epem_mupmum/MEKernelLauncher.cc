//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 2.9.5, 2021-08-22
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "MEKernelLauncher.h"

#include "checkCuda.h"

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  MEKernelLauncher::MEKernelLauncher( int ngpublocks, int ngputhreads )
    : m_ngpublocks( ngpublocks )
    , m_ngputhreads( ngputhreads )
    , m_nevt( ngpublocks * ngputhreads )
    , m_momenta( nullptr )
    , m_MEs( nullptr )
  {
    // Allocate the buffer for the input momenta
    checkCuda( cudaMalloc( &m_momenta, nbytesMomenta() * sizeof(fptype) ) );

    // Allocate the buffer for the input MEs
    checkCuda( cudaMalloc( &m_MEs, nbytesMEs() * sizeof(fptype) ) );
  }
#else
  MEKernelLauncher::MEKernelLauncher( int nevt )
    : m_nevt( nevt )
    , m_momenta( nullptr )
    , m_MEs( nullptr )
  {
    // Allocate the buffer for the input momenta
    m_momenta = new( std::align_val_t{ cppAlign } ) fptype[ nbytesMomenta() ]();

    // Allocate the buffer for the output MEs
    m_MEs = new( std::align_val_t{ cppAlign } ) fptype[ nbytesMEs() ]();
  }
#endif

  //--------------------------------------------------------------------------

  MEKernelLauncher::~MEKernelLauncher()
  {
#ifndef __CUDACC__
    // Deallocate the host buffer for the input momenta
    ::operator delete( m_momenta, std::align_val_t{ cppAlign } );

    // Deallocate the host buffer for the output MEs
    ::operator delete( m_MEs, std::align_val_t{ cppAlign } );
#else
    // Deallocate the device buffer for the input momenta
    checkCuda( cudaFree( m_momenta ) );

    // Deallocate the device buffer for the output MEs
    checkCuda( cudaFree( m_MEs ) );
#endif
  }

  //--------------------------------------------------------------------------

  void MEKernelLauncher::computeMEs() const
  {
  }

  //--------------------------------------------------------------------------

}
