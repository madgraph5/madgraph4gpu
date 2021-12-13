//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 2.9.5, 2021-08-22
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "MEKernelLauncher.h"

#include "checkCuda.h"
#include "CPPProcess.h"

#include <stdexcept>

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  MEKernelLauncher::MEKernelLauncher( int ngpublocks, int ngputhreads, UseHstBuffers useHstBuffers )
    : m_ngpublocks( ngpublocks )
    , m_ngputhreads( ngputhreads )
    , m_nevt( ngpublocks * ngputhreads )
    , m_hstMomenta( nullptr )
    , m_devMomenta( nullptr )
    , m_hstMEs( nullptr )
    , m_devMEs( nullptr )
    , m_hstIsGoodHel( nullptr )
  {
    // Allocate the host buffer for the input momenta
    if ( useHstBuffers == UseHstMomenta || useHstBuffers == UseBoth )
      checkCuda( cudaMallocHost( &m_hstMomenta, nMomenta() * sizeof(fptype) ) );

    // Allocate the device buffer for the input momenta
    checkCuda( cudaMalloc( &m_devMomenta, nMomenta() * sizeof(fptype) ) );

    // Allocate the host buffer for the output MEs
    if ( useHstBuffers == UseHstMEs || useHstBuffers == UseBoth )
      checkCuda( cudaMallocHost( &m_hstMEs, nMEs() * sizeof(fptype) ) );

    // Allocate the device buffer for the output MEs
    checkCuda( cudaMalloc( &m_devMEs, nMEs() * sizeof(fptype) ) );
  }
#else
  MEKernelLauncher::MEKernelLauncher( int nevt )
    : m_nevt( nevt )
    , m_hstMomenta( nullptr )
    , m_hstMEs( nullptr )
    , m_hstIsGoodHel( nullptr )
  {
    // Allocate the host buffer for the input momenta
    m_hstMomenta = new( std::align_val_t{ cppAlign } ) fptype[ nMomenta() ]();

    // Allocate the host buffer for the output MEs
    m_hstMEs = new( std::align_val_t{ cppAlign } ) fptype[ nMEs() ]();
  }
#endif

  //--------------------------------------------------------------------------

  MEKernelLauncher::~MEKernelLauncher()
  {
#ifdef __CUDACC__
    // Deallocate the host buffer for the input momenta
    if ( m_hstMomenta ) checkCuda( cudaFreeHost( m_hstMomenta ) );

    // Deallocate the device buffer for the input momenta
    checkCuda( cudaFree( m_devMomenta ) );

    // Deallocate the host buffer for the output MEs
    if ( m_hstMEs ) checkCuda( cudaFreeHost( m_hstMEs ) );

    // Deallocate the device buffer for the output MEs
    checkCuda( cudaFree( m_devMEs ) );

    // Deallocate the host buffer for the helicity mask
    if ( m_hstIsGoodHel ) checkCuda( cudaFreeHost( m_hstIsGoodHel ) );
#else
    // Deallocate the host buffer for the input momenta
    ::operator delete( m_hstMomenta, std::align_val_t{ cppAlign } );

    // Deallocate the host buffer for the output MEs
    ::operator delete( m_hstMEs, std::align_val_t{ cppAlign } );

    // Deallocate the host buffer for the helicity mask
    if ( m_hstIsGoodHel ) delete( m_hstIsGoodHel );
#endif
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void MEKernelLauncher::copyDevMomentaToHstMomenta() const
  {
    if ( !m_hstMomenta ) throw std::logic_error( "Cannot copyDevMomentaToHstMomenta unless UseHstMomenta or UseBoth are specified" );
    checkCuda( cudaMemcpy( m_hstMomenta, m_devMomenta, nMomenta() * sizeof(fptype), cudaMemcpyDeviceToHost ) );
  }
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void MEKernelLauncher::copyHstMomentaToDevMomenta() const
  {
    if ( !m_hstMomenta ) throw std::logic_error( "Cannot copyHstMomentaToDevMomenta unless UseHstMomenta or UseBoth are specified" );
    checkCuda( cudaMemcpy( m_devMomenta, m_hstMomenta, nMomenta() * sizeof(fptype), cudaMemcpyHostToDevice ) );
  }
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void MEKernelLauncher::computeDevMEsFromDevMomenta() const
  {
#ifndef MGONGPU_NSIGHT_DEBUG
    gProc::sigmaKin<<<m_ngpublocks, m_ngputhreads>>>( m_devMomenta, m_devMEs );
#else
    gProc::sigmaKin<<<m_ngpublocks, m_gputhreads, mgOnGpu::ntpbMAX * sizeof(float)>>>( m_devMomenta, m_devMEs );
#endif
    checkCuda( cudaPeekAtLastError() );
  }
#else
  void MEKernelLauncher::computeHstMEsFromHstMomenta() const
  {
    Proc::sigmaKin( m_hstMomenta, m_hstMEs, m_nevt );
  }
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void MEKernelLauncher::copyDevMEsToHstMEs() const
  {
    if ( !m_hstMEs ) throw std::logic_error( "Cannot copyDevMEsToHstMEs unless UseHstMEs or UseBoth are specified" );
    checkCuda( cudaMemcpy( m_hstMEs, m_devMEs, nMEs() * sizeof(fptype), cudaMemcpyDeviceToHost ) );
  }
#endif

  //--------------------------------------------------------------------------

#ifndef __CUDACC__
  void MEKernelLauncher::computeGoodHelFromHstMomenta()
  {
    if ( !m_hstIsGoodHel )
    {
      // Allocate the host buffer for the helicity mask
      m_hstIsGoodHel = new bool[ ncomb ]();
      // Compute the good helicity mask on the host
      Proc::sigmaKin_getGoodHel( m_hstMomenta, m_hstMEs, m_hstIsGoodHel, m_nevt );
      // Copy the helicity mask to static memory on the host
      // [FIXME! REMOVE THIS STATIC THAT BREAKS MULTITHREADING?]
      Proc::sigmaKin_setGoodHel( m_hstIsGoodHel );
    }
  }
#else
  void MEKernelLauncher::computeGoodHelFromDevMomenta()
  {
    if ( !m_hstIsGoodHel )
    {
      // Allocate the host buffer for the helicity mask
      checkCuda( cudaMallocHost( &m_hstIsGoodHel, ncomb * sizeof(bool) ) );
      // Allocate the device buffer for the helicity mask
      bool* devIsGoodHel( nullptr );
      checkCuda( cudaMalloc( &devIsGoodHel, ncomb * sizeof(bool) ) );
      // Compute the good helicity mask on the device
      gProc::sigmaKin_getGoodHel<<<m_ngpublocks, m_ngputhreads>>>( m_devMomenta, m_devMEs, devIsGoodHel );
      checkCuda( cudaPeekAtLastError() );
      // Copy the helicity mask from the device to the host
      checkCuda( cudaMemcpy( m_hstIsGoodHel, devIsGoodHel, ncomb * sizeof(bool), cudaMemcpyDeviceToHost ) );
      // Deallocate the device buffer for the helicity mask
      checkCuda( cudaFree( devIsGoodHel ) );
      // Copy the helicity mask to constant memory on the device
      gProc::sigmaKin_setGoodHel( m_hstIsGoodHel );
    }
  }
#endif

  //--------------------------------------------------------------------------

  /*
  void MEKernelLauncher::setGoodHel( const bool* isGoodHel )
  {
    if ( !m_hstIsGoodHel )
    {
#ifdef __CUDACC__
      // Allocate the host buffer for the helicity mask
      checkCuda( cudaMallocHost( &m_hstIsGoodHel, ncomb * sizeof(bool) ) );
#else
      // Allocate the host buffer for the helicity mask
      m_hstIsGoodHel = new bool[ ncomb ]();
#endif
      // Set the helicity mask from user input
      for ( int ihel = 0; ihel < ncomb; ihel++ ) m_hstIsGoodHel[ihel] = isGoodHel[ihel];
#ifdef __CUDACC__
      // Copy the helicity mask to constant memory on the device
      gProc::sigmaKin_setGoodHel( m_hstIsGoodHel );
#else
      // Copy the helicity mask to static memory on the host
      // [FIXME! REMOVE THIS STATIC THAT BREAKS MULTITHREADING?]
      Proc::sigmaKin_setGoodHel( m_hstIsGoodHel );
#endif
    }
    else throw std::logic_error( "Cannot setGoodHel: helicity mask has already been set or computed" );
  }
  */
  
  //--------------------------------------------------------------------------

}
