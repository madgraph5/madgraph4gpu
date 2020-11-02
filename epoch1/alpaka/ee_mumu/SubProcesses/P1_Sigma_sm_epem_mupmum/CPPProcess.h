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

#include <cuda_to_cupla.hpp>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "Parameters_sm.h"

//--------------------------------------------------------------------------

#define checkCupla( code )                       \
  { assertCupla( code, __FILE__, __LINE__ ); }

inline void assertCupla( cuplaError_t code, const char *file, int line, bool abort = true )
{
  if ( code != cuplaSuccess )
  {
    printf( "GPUassert: %s %s %d\n", cuplaGetErrorString(code), file, line );
    if ( abort ) assert( code == cuplaSuccess );
  }
}

//--------------------------------------------------------------------------

namespace gProc
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

//    cxtype * * amp;


    // Pointer to the model parameters
    Parameters_sm * pars;

    // vector with external particle masses
    std::vector<fptype> mME;

    // Initial particle ids
    int id1, id2;

  };

  //--------------------------------------------------------------------------

struct sigmaKin
{
  template< typename T_Acc >
  ALPAKA_FN_ACC
#if defined MGONGPU_WFMEM_GLOBAL
  void operator()( T_Acc const &acc, const fptype* allmomenta, fptype* output, cxtypeparam* tmpWFs ) const;
#else
  void operator()( T_Acc const &acc, const fptype* allmomenta, fptype* output ) const;
#endif
};

  //--------------------------------------------------------------------------

#if defined MGONGPU_WFMEM_SHARED
  int sigmakin_sharedmem_nbytes( const int ntpb ); // input: #threads per block
#endif

  //--------------------------------------------------------------------------

}

#include "CPPProcess.cc"

#endif  // MG5_Sigma_sm_epem_mupmum_H
