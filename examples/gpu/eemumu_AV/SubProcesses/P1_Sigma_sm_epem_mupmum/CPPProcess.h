//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_epem_mupmum_H
#define MG5_Sigma_sm_epem_mupmum_H

#include <cassert>
#include <complex>
#include <vector>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "Parameters_sm.h"

//--------------------------------------------------------------------------

#ifdef __CUDACC__

#define checkCuda( code )                       \
  { assertCuda( code, __FILE__, __LINE__ ); }

inline void assertCuda( cudaError_t code, const char *file, int line, bool abort = true )
{
  if ( code != cudaSuccess )
  {
    printf( "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line );
    if ( abort ) assert( code == cudaSuccess );
  }
}

#endif

//--------------------------------------------------------------------------

#ifdef __CUDACC__
namespace gProc
#else
namespace Proc
#endif
{

  //==========================================================================
  // A class for calculating the matrix elements for
  // Process: e+ e- > mu+ mu- WEIGHTED<=4 @1
  //--------------------------------------------------------------------------

  class CPPProcess
  {
  public:

    CPPProcess( int numiterations, int gpublocks, int gputhreads, bool verbose = false );

    ~CPPProcess();

    // Initialize process.
    virtual void initProc(std::string param_card_name);


    virtual int code() const {return 1;}

    const std::vector<fptype> &getMasses() const;

    void setInitial(int inid1, int inid2)
    {
      id1 = inid1;
      id2 = inid2;
    }

    int getDim() const {return dim;}

    int getNIOParticles() const {return nexternal;}


    // Constants for array limits
    static const int ninitial = mgOnGpu::npari;
    static const int nexternal = mgOnGpu::npar;
    static const int nprocesses = 1;

  private:
    int m_numiterations;
    // gpu variables
    int gpu_nblocks;
    int gpu_nthreads;
    int dim;  // gpu_nblocks * gpu_nthreads;

    // print verbose info
    bool m_verbose;

    static const int nwavefuncs = 6;
    static const int namplitudes = 2;
    static const int ncomb = 16;
    static const int wrows = 6;
    // static const int nioparticles = 4;

    cxtype * * amp;


    // Pointer to the model parameters
    Parameters_sm * pars;

    // vector with external particle masses
    std::vector<fptype> mME;

    // Initial particle ids
    int id1, id2;

  };

  //--------------------------------------------------------------------------

#if defined MGONGPU_LAYOUT_ASA
  void sigmakin_setNeppM( const int neppM ); // input: n_events_per_page for momenta AOSOA (nevt=npagM*neppM)
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __global__
#if defined MGONGPU_WFMEM_GLOBAL
  void sigmaKin( const fptype* allmomenta, fptype* output, cxtype* tmpWFs );
#else
  void sigmaKin( const fptype* allmomenta, fptype* output );
#endif
#else
  void sigmaKin( const fptype* allmomenta, fptype* output, const int nevt );
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
#if defined MGONGPU_WFMEM_SHARED
  int sigmakin_sharedmem_nbytes( const int ntpb ); // input: #threads per block
#endif
#endif

  //--------------------------------------------------------------------------

}

#endif  // MG5_Sigma_sm_epem_mupmum_H
