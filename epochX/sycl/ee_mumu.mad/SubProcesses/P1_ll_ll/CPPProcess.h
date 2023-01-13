//==========================================================================
// This file has been automatically generated for SYCL standalone by
// MadGraph5_aMC@NLO v. 3.4.0_lo_vect, 2022-05-06
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_epem_mupmum_H
#define MG5_Sigma_sm_epem_mupmum_H

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include "Parameters_sm.h"

#include <vector>

//--------------------------------------------------------------------------

namespace Proc
{
  namespace dependentCouplings = Parameters_sm_dependentCouplings;
  namespace independentCouplings = Parameters_sm_independentCouplings;

  template <typename T>
  constexpr T helicities[] { 
    -1, 1, 1, -1,
    -1, 1, 1, 1,
    -1, 1, -1, -1,
    -1, 1, -1, 1,
    -1, -1, 1, -1,
    -1, -1, 1, 1,
    -1, -1, -1, -1,
    -1, -1, -1, 1,
    1, 1, 1, -1,
    1, 1, 1, 1,
    1, 1, -1, -1,
    1, 1, -1, 1,
    1, -1, 1, -1,
    1, -1, 1, 1,
    1, -1, -1, -1,
    1, -1, -1, 1
  };

#ifdef MGONGPU_HARDCODE_PARAM
  template <typename FPType>
  constexpr FPType independent_parameters[] { (FPType)Parameters_sm::mdl_MZ, (FPType)Parameters_sm::mdl_WZ };
#endif

  //==========================================================================
  // A class for calculating the matrix elements for
  // Process: e+ e- > mu+ mu- WEIGHTED<=4 @1
  //--------------------------------------------------------------------------

  class CPPProcess
  {
  public:

    // Constructor (from command line arguments)
    CPPProcess( size_t numiterations, size_t gpublocks, size_t gputhreads, bool verbose = false, bool debug = false );

    // Destructor
    ~CPPProcess();

    // Initialize process (read model parameters from file)
    virtual void initProc( const std::string& param_card_name );

#ifndef MGONGPU_HARDCODE_PARAM
    // Pointer accessors
    cxtype* get_tIPC_ptr();
    const cxtype* get_tIPC_ptr() const;

    fptype* get_tIPD_ptr();
    const fptype* get_tIPD_ptr() const;
#endif

    // Other methods of this instance (???)
    //const std::vector<fptype>& getMasses() const { return m_masses; }
    //virtual int code() const{ return 1; }
    //void setInitial( int inid1, int inid2 ){ id1 = inid1; id2 = inid2; }
    //int getDim() const { return dim; }
    //int getNIOParticles() const { return nexternal; } // nexternal was nioparticles

    // Accessors (unused so far: add four of them only to fix a clang build warning)
    size_t numiterations() const { return m_numiterations; }
    size_t gpublocks() const { return m_ngpublocks; }
    size_t gputhreads() const { return m_ngputhreads; }
    //bool verbose() const { return m_verbose; }
    bool debug() const { return m_debug; }

  public:

    // Hardcoded parameters for this process (constant class variables)
    //static const int ninitial = mgOnGpu::npari;
    //static const int nexternal = 4; // mgOnGpu::npar (nexternal was nioparticles)
    //static const int nprocesses = 1; // FIXME: assume process.nprocesses == 1
    //static const int nwavefuncs = 6; // mgOnGpu::nwf
    //static const int namplitudes = 2;
    //static const int ncomb = 16; // mgOnGpu::ncomb
    //static const int wrows = 5; // mgOnGpu::nw6;

  private:

    // Command line arguments (constructor)
    size_t m_numiterations; // number of iterations (each iteration has nblocks*nthreads events)
    size_t m_ngpublocks; // number of GPU blocks in one grid (i.e. one iteration)
    size_t m_ngputhreads; // number of GPU threads in a block
    bool m_verbose;
    bool m_debug;

#ifndef MGONGPU_HARDCODE_PARAM
    // Physics model parameters to be read from file (initProc function)
    Parameters_sm* m_pars;
    std::vector<fptype> m_masses; // external particle masses

    // Physics parameters (masses, coupling, etc...)
    cxtype m_tIPC[independentCouplings::nicoup];
    fptype m_tIPD[mgOnGpu::nparams];
#endif

    // Other variables of this instance (???)
    //int id1, id2; // initial particle ids
    //cxtype** amp; // ???

  };

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin_getGoodHel( const vector4* __restrict__ allmomenta, // input: momenta[nevt*npar*4]
                            bool* isGoodHel,                        // output: isGoodHel[ncomb] - device array
                            const signed char* __restrict__ cHel,
                            const cxtype_sv* __restrict__ COUPs,
                            const fptype* __restrict__ cIPD
                            );

  //--------------------------------------------------------------------------

  size_t sigmaKin_setGoodHel( const bool* isGoodHel, size_t* goodHel ); // input: isGoodHel[ncomb] - host array

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  fptype_sv sigmaKin( const vector4* __restrict__ allmomenta, // input: momenta[]
                      const fptype_sv* __restrict__ rndhel,   // input: random numbers[] for helicity selection
                      const fptype_sv* __restrict__ rndcol,   // input: random numbers[] for color selection
                      int_sv* __restrict__ selhel,            // output: helicity selection[]
                      int_sv* __restrict__ selcol,            // output: color selection[]
                      #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                          const size_t channelId,             // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
                      #endif
                      const signed char* __restrict__ cHel,
                      const cxtype_sv* __restrict__ COUPs,
                      const fptype* __restrict__ cIPD,
                      const size_t* __restrict__ cNGoodHel,
                      const size_t* __restrict__ cGoodHel
                    );

  //--------------------------------------------------------------------------
}

#endif // MG5_Sigma_sm_epem_mupmum_H
