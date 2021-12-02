//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 2.9.5, 2021-08-22
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "MEKernelLauncher.h"

#include "checkCuda.h"
#include "CPPProcess.h"

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
    , m_devMomenta( nullptr )
    , m_devMEs( nullptr )
  {
    // Allocate the device buffer for the input momenta
    checkCuda( cudaMalloc( &m_devMomenta, nbytesMomenta() * sizeof(fptype) ) );

    // Allocate the device buffer for the output MEs
    checkCuda( cudaMalloc( &m_devMEs, nbytesMEs() * sizeof(fptype) ) );
  }
#else
  MEKernelLauncher::MEKernelLauncher( int nevt )
    : m_nevt( nevt )
    , m_hstMomenta( nullptr )
    , m_hstMEs( nullptr )
  {
    // Allocate the host buffer for the input momenta
    m_hstMomenta = new( std::align_val_t{ cppAlign } ) fptype[ nbytesMomenta() ]();

    // Allocate the host buffer for the output MEs
    m_hstMEs = new( std::align_val_t{ cppAlign } ) fptype[ nbytesMEs() ]();
  }
#endif

  //--------------------------------------------------------------------------

  MEKernelLauncher::~MEKernelLauncher()
  {
#ifdef __CUDACC__
    // Deallocate the device buffer for the input momenta
    checkCuda( cudaFree( m_devMomenta ) );

    // Deallocate the device buffer for the output MEs
    checkCuda( cudaFree( m_devMEs ) );
#else
    // Deallocate the host buffer for the input momenta
    ::operator delete( m_hstMomenta, std::align_val_t{ cppAlign } );

    // Deallocate the host buffer for the output MEs
    ::operator delete( m_hstMEs, std::align_val_t{ cppAlign } );
#endif
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void MEKernelLauncher::computeDevMEs() const
  {
#ifndef MGONGPU_NSIGHT_DEBUG
    gProc::sigmaKin<<<m_ngpublocks, m_ngputhreads>>>( m_devMomenta, m_devMEs );
#else
    gProc::sigmaKin<<<m_ngpublocks, m_gputhreads, mgOnGpu::ntpbMAX * sizeof(float)>>>( m_devMomenta, m_devMEs );
#endif
  }
#else
  void MEKernelLauncher::computeHstMEs() const
  {
  Proc::sigmaKin( m_hstMomenta, m_hstMEs, m_nevt );
  }
#endif

  //--------------------------------------------------------------------------

}
