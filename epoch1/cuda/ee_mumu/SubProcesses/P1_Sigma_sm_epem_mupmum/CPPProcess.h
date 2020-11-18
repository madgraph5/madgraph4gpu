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
#include <fstream>
#include <iostream>

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
    printf( "GPUassert: %s %s:%d\n", cudaGetErrorString(code), file, line );
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
    //static const int nprocesses = 1; // FIXME: assume process.nprocesses == 1

  private:

    int m_numiterations;

    // gpu variables
    int gpu_nblocks;
    int gpu_nthreads;
    int dim; // gpu_nblocks * gpu_nthreads;

    bool m_verbose; // print verbose info

    static const int nwavefuncs = 6;
    static const int namplitudes = 2;
    static const int ncomb = 16;
    static const int wrows = 6;

    cxtype** amp;

    // Pointer to the model parameters
    Parameters_sm * pars;

    // vector with external particle masses
    std::vector<fptype> mME;

    // Initial particle ids
    int id1, id2;

  };

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __global__
  void sigmaKin_getGoodHel( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            bool* isGoodHel );        // output: isGoodHel[ncomb] - device array
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void sigmaKin_setGoodHel( const bool* isGoodHel ); // input: isGoodHel[ncomb] - host array
#endif

  //--------------------------------------------------------------------------

  __global__
  void sigmaKin( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                 fptype* allMEs            // output: allMEs[nevt], final |M|^2 averaged over all helicities
#ifndef __CUDACC__
                 , const int nevt          // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                 );

  //--------------------------------------------------------------------------

}

#endif  // MG5_Sigma_sm_epem_mupmum_H
