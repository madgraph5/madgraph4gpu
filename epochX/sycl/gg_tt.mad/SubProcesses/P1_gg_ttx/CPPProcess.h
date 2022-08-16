//==========================================================================
// This file has been automatically generated for SYCL standalone by
// MadGraph5_aMC@NLO v. 3.4.0_lo_vect, 2022-05-06
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_gg_ttx_H
#define MG5_Sigma_sm_gg_ttx_H

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include "Parameters_sm.h"

#include <vector>

//--------------------------------------------------------------------------

namespace Proc
{

  //==========================================================================
  // A class for calculating the matrix elements for
  // Process: g g > t t~ WEIGHTED<=2 @1
  //--------------------------------------------------------------------------

  class CPPProcess
  {
  public:

    // Constructor (from command line arguments)
    CPPProcess( int numiterations, int gpublocks, int gputhreads, bool verbose = false, bool debug = false );

    // Destructor
    ~CPPProcess();

    // Initialize process (read model parameters from file)
    virtual void initProc( const std::string& param_card_name );

    // Pointer accessors
    const short* get_tHel_ptr() const;

    cxtype* get_tIPC_ptr();
    const cxtype* get_tIPC_ptr() const;

    fptype* get_tIPD_ptr();
    const fptype* get_tIPD_ptr() const;

    // Other methods of this instance (???)
    //const std::vector<fptype>& getMasses() const { return m_masses; }
    //virtual int code() const{ return 1; }
    //void setInitial( int inid1, int inid2 ){ id1 = inid1; id2 = inid2; }
    //int getDim() const { return dim; }
    //int getNIOParticles() const { return nexternal; } // nexternal was nioparticles

    // Accessors (unused so far: add four of them only to fix a clang build warning)
    int numiterations() const { return m_numiterations; }
    int gpublocks() const { return m_ngpublocks; }
    int gputhreads() const { return m_ngputhreads; }
    //bool verbose() const { return m_verbose; }
    bool debug() const { return m_debug; }

  public:

    // Hardcoded parameters for this process (constant class variables)
    //static const int ninitial = mgOnGpu::npari;
    //static const int nexternal = 4; // mgOnGpu::npar (nexternal was nioparticles)
    //static const int nprocesses = 1; // FIXME: assume process.nprocesses == 1
    //static const int nwavefuncs = 6; // mgOnGpu::nwf
    //static const int namplitudes = 3;
    //static const int ncomb = 16; // mgOnGpu::ncomb
    //static const int wrows = 5; // mgOnGpu::nw6;

  private:

    // Command line arguments (constructor)
    int m_numiterations; // number of iterations (each iteration has nblocks*nthreads events)
    int m_ngpublocks; // number of GPU blocks in one grid (i.e. one iteration)
    int m_ngputhreads; // number of GPU threads in a block
    bool m_verbose;
    bool m_debug;

    // Physics model parameters to be read from file (initProc function)
    Parameters_sm* m_pars;
    std::vector<fptype> m_masses; // external particle masses

    // Helicities for the process
    const short m_tHel[mgOnGpu::ncomb][mgOnGpu::npar];

    // Physics parameters (masses, coupling, etc...)
    cxtype m_tIPC[mgOnGpu::ncouplings];
    fptype m_tIPD[mgOnGpu::nparams];

    // Other variables of this instance (???)
    //int id1, id2; // initial particle ids
    //cxtype** amp; // ???

  };

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin_getGoodHel( const fptype* __restrict__ allmomenta, // input: momenta[nevt*npar*4]
                            bool* isGoodHel,          // output: isGoodHel[ncomb] - device array
                            const short* __restrict__ cHel,
                            const fptype* __restrict__ cIPC,
                            const fptype* __restrict__ cIPD
                            );

  //--------------------------------------------------------------------------

  int sigmaKin_setGoodHel( const bool* isGoodHel, int* goodHel ); // input: isGoodHel[ncomb] - host array

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  fptype sigmaKin( const fptype* __restrict__ allmomenta, // input: momenta[nevt*npar*4]
                   const short* __restrict__ cHel,
                   const fptype* __restrict__ cIPC,
                   const fptype* __restrict__ cIPD,
                   const int* __restrict__ cNGoodHel,
                   const int* __restrict__ cGoodHel
                 );

  //--------------------------------------------------------------------------
}

#endif // MG5_Sigma_sm_gg_ttx_H
