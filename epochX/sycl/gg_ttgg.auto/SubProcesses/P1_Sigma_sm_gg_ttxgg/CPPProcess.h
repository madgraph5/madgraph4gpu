//==========================================================================
// This file has been automatically generated for SYCL/C++ standalone by
// MadGraph5_aMC@NLO v. 2.9.5, 2021-08-22
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_gg_ttxgg_H
#define MG5_Sigma_sm_gg_ttxgg_H

#include <cassert>
#include <complex>
#include <iostream>
#include <vector>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include "Parameters_sm.h"

//--------------------------------------------------------------------------

namespace Proc
{

  //==========================================================================
  // A class for calculating the matrix elements for
  // Process: g g > t t~ g g WEIGHTED<=4 @1
  //--------------------------------------------------------------------------

  class CPPProcess
  {
  public:

    // Constructor (from command line arguments)
    CPPProcess( int numiterations, int gpublocks, int gputhreads, bool verbose = false, bool debug = false );

    // Destructor
    ~CPPProcess();

    // Device Pointer Accessors
    short * get_cHel_ptr() const;
    fptype * get_cIPC_ptr();
    fptype * get_cIPD_ptr();

    // Initialize process (read model parameters from file)
    virtual void initProc( const std::string& param_card_name );

    // Retrieve the compiler that was used to build this module
    static const std::string getCompiler();

    // Other methods of this instance (???)
    //const std::vector<fptype>& getMasses() const { return m_masses; }
    //virtual int code() const{ return 1; }
    //void setInitial( int inid1, int inid2 ){ id1 = inid1; id2 = inid2; }
    //int getDim() const { return dim; }
    //int getNIOParticles() const { return nexternal; } // nexternal was nioparticles

    // Accessors (unused so far: add them to fix a clang build warning)
    //int numiterations() const { return m_numiterations; }
    //int gpublocks() const { return m_ngpublocks; }
    //int gputhreads() const { return m_ngputhreads; }
    //bool verbose() const { return m_verbose; }
    //bool debug() const { return m_debug; }

  public:

    // Hardcoded parameters for this process (constant class variables)
    //static const int ninitial = mgOnGpu::npari;
    static const int nexternal = 6; // mgOnGpu::npar (nexternal was nioparticles)
    //static const int nprocesses = 1; // FIXME: assume process.nprocesses == 1
    //static const int nwavefuncs = 6; // mgOnGpu::nwf
    //static const int namplitudes = 159;
    //static const int ncomb = 64; // mgOnGpu::ncomb
    //static const int wrows = 63; // mgOnGpu::nw6;

  private:

    // Command line arguments (constructor)
    int m_numiterations; // number of iterations (each iteration has nblocks*nthreads events)
    int m_ngpublocks; // number of GPU blocks in one grid (i.e. one iteration)
    int m_ngputhreads; // number of GPU threads in a block
    bool m_verbose;
    bool m_debug;

    // Physics model parameters to be read from file (initProc function)
#ifndef MGONGPU_HARDCODE_CIPC
    Parameters_sm* m_pars;
#endif
    std::vector<fptype> m_masses; // external particle masses

    // Other variables of this instance (???)
    //int id1, id2; // initial particle ids
    //cxtype** amp; // ???

  };

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin_getGoodHel(
          const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
          fptype* allMEs,           // output: allMEs[nevt], |M|^2 final_avg_over_helicities
          bool * isGoodHel,         // output: isGoodHel[ncomb] - device array
          size_t ievt,
          short *cHel,
          const fptype *cIPC,
          const fptype *cIPD
          );

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin_setGoodHel(
          const bool * isGoodHel,  // input: isGoodHel[ncomb] - host array
          int * cNGoodHel_ptr,
          int* cGoodHel_ptr
          );

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin(
          const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
          fptype* allMEs,           // output: allMEs[nevt], |M|^2 final_avg_over_helicities
          size_t ievt,
          short *cHel,
          const fptype *cIPC,
          const fptype *cIPD,
          int *cNGoodHel,
          int *cGoodHel
          );

  //--------------------------------------------------------------------------
}

#endif // MG5_Sigma_sm_gg_ttxgg_H
