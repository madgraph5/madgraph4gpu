// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: O. Mattelaer, S. Roiser, J. Teig, A. Valassi (2020-2025) for the MG5aMC CUDACPP plugin.
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.6.3, 2025-06-12
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_uux_ttx_H
#define MG5_Sigma_sm_uux_ttx_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuVectors.h"

#include "GpuAbstraction.h"
#include "Parameters_sm.h"

#include <vector>

//--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //==========================================================================
  // A class for calculating the matrix elements for
  // Process: u u~ > t t~ WEIGHTED<=2
  // Process: c c~ > t t~ WEIGHTED<=2
  // Process: d d~ > t t~ WEIGHTED<=2
  // Process: s s~ > t t~ WEIGHTED<=2
  //--------------------------------------------------------------------------

  class CPPProcess
  {
  public: /* clang-format off */

    // Constructor (from command line arguments)
    CPPProcess( bool verbose = false, bool debug = false );

    // Destructor
    ~CPPProcess();

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

    // Accessors (unused so far: add four of them only to fix a clang build warning)
    //bool verbose() const { return m_verbose; }
    bool debug() const { return m_debug; }

  public:

    // Process-independent compile-time constants
    static constexpr int np4 = 4; // dimensions of 4-momenta (E,px,py,pz)
    static constexpr int nw6 = 6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)

    // Process-dependent compile-time constants
    static constexpr int npari = 2; // #particles in the initial state (incoming): e.g. 2 (e+ e-) for e+ e- -> mu+ mu-
    static constexpr int nparf = 2; // #particles in the final state (outgoing): e.g. 2 (mu+ mu-) for e+ e- -> mu+ mu-
    static constexpr int npar = npari + nparf; // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
    static constexpr int ncomb = 16; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)
    static constexpr int ndiagrams = 1; // #Feynman diagrams: e.g. 3 for e+ e- -> mu+ mu-
    static constexpr int ncolor = 2; // the number of leading colors: e.g. 1 for e+ e- -> mu+ mu-

    // Process-dependent (and generation-choice-dependent) compile-time constants
    static constexpr int ndiagramgroups = 1; // #groups of Feynman diagrams (with at most 100 diagrams per group)

    // Hardcoded parameters for this process (constant class variables)
    // [NB: this class assumes nprocesses==1 i.e. a single DSIG1 and no DSIG2 in Fortran (#272 and #343)]
    // [NB: these parameters (e.g. nwf) are P1-specific, i.e. they are different for different P1 subdirectories (#644)]
    // [NB: I was unable to get the right value of nwf in CPPProcess.h directly, so I added it with a hack after generating CPPProcess.cc (#644)]
    static const int nwf = 5; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)

    // Process-dependent (but event-independent) parameters and couplings
    // Note: in the Python code generator, nIPD == nparam, while nIPC <= nicoup, because (see #823)
    // nIPC may vary from one P*/CPPProcess.cc to another, while nicoup is defined in src/Param.h and is common to all P*
    static const int nIPD = 2; // SM independent parameters (FIXME? rename as sm_IndepParam?)
    static const int nIPC = 0; // SM independent couplings (FIXME? rename as sm_IndepCoupl?)

    // Other variables of this instance (???)
    //static const int ninitial = CPPProcess::npari;
    //static const int nexternal = 4; // CPPProcess::npar (nexternal was nioparticles)
    //static const int namplitudes = 1;
    //static const int ncomb = 16; // CPPProcess::ncomb

  private: /* clang-format on */

    // Command line arguments (constructor)
    bool m_verbose;
    bool m_debug;

    // Physics model parameters to be read from file (initProc function)
#ifndef MGONGPU_HARDCODE_PARAM
    Parameters_sm* m_pars;
#endif
    std::vector<fptype> m_masses; // external particle masses

    // Other variables of this instance (???)
    //int id1, id2; // initial particle ids
    //cxtype** amp; // ???
  };

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  __global__ void
  computeDependentCouplings( const fptype* allgs,    // input: Gs[nevt]
                             fptype* allcouplings ); // output: couplings[nevt*ndcoup*2]
#else
  __global__ void
  computeDependentCouplings( const fptype* allgs,  // input: Gs[nevt]
                             fptype* allcouplings, // output: couplings[nevt*ndcoup*2]
                             const int nevt );     // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL /* clang-format off */
  void
  sigmaKin_getGoodHel( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
                       const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
                       fptype* allMEs,             // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                       fptype* allNumerators,      // output: numerators[nevt], running_sum_over_helicities
                       fptype* allDenominators,    // output: denominators[nevt], running_sum_over_helicities
#endif
                       fptype* allJamps,           // output: jamp[ncolor*2*nevt]
                       fptype* allWfs,             // output: wf[nwf*nw6*2*nevt]
                       bool* isGoodHel,            // output: isGoodHel[ncomb] - device array (GPU device implementation)
                       const int nevt );           // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#else
  void
  sigmaKin_getGoodHel( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
                       const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
                       fptype* allMEs,             // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                       fptype* allNumerators,      // output: numerators[nevt], running_sum_over_helicities
                       fptype* allDenominators,    // output: denominators[nevt], running_sum_over_helicities
#endif
                       bool* isGoodHel,            // output: isGoodHel[ncomb] - host array (C++ implementation)
                       const int nevt );           // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif /* clang-format on */

  //--------------------------------------------------------------------------

  int                                           // output: nGoodHel (the number of good helicity combinations out of ncomb)
  sigmaKin_setGoodHel( const bool* isGoodHel ); // input: isGoodHel[ncomb] - host array

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL /* clang-format off */
  void
  sigmaKin( const fptype* allmomenta,           // input: momenta[nevt*npar*4]
            const fptype* allcouplings,         // input: couplings[nevt*ndcoup*2]
            const fptype* allrndhel,            // input: random numbers[nevt] for helicity selection
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            const fptype* allrndcol,            // input: random numbers[nevt] for color selection
            const unsigned int* allChannelIds,  // input: channelIds[nevt] (1 to #diagrams); nullptr to disable single-diagram enhancement (fix #899/#911)
#endif
            fptype* allMEs,                     // output: allMEs[nevt], |M|^2 final_avg_over_helicities
            int* allselhel,                     // output: helicity selection[nevt]
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            int* allselcol,                     // output: helicity selection[nevt]
            fptype* colAllJamp2s,               // tmp: allJamp2s super-buffer for ncolor individual colors, running sum over colors and helicities
            fptype* ghelAllNumerators,          // tmp: allNumerators super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            fptype* ghelAllDenominators,        // tmp: allDenominators super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
#endif
            fptype* ghelAllMEs,                 // tmp: allMEs super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            fptype* ghelAllJamps,               // tmp: allJamps super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            fptype* ghelAllWfs,                 // tmp: allWfs super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            fptype2* ghelAllBlasTmp,            // tmp: allBlasTmp super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            gpuBlasHandle_t* ghelBlasHandles,   // input: cuBLAS/hipBLAS handles (index is ighel: only the first nGoodHel <= ncomb are non-null)
            gpuStream_t* ghelStreams,           // input: cuda streams (index is ighel: only the first nGoodHel <= ncomb are non-null)
            const int gpublocks,                // input: cuda gpublocks
            const int gputhreads );             // input: cuda gputhreads
#else
  void
  sigmaKin( const fptype* allmomenta,           // input: momenta[nevt*npar*4]
            const fptype* allcouplings,         // input: couplings[nevt*ndcoup*2]
            const fptype* allrndhel,            // input: random numbers[nevt] for helicity selection
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            const fptype* allrndcol,            // input: random numbers[nevt] for color selection
            const unsigned int* allChannelIds,  // input: channelIds[nevt] (1 to #diagrams); nullptr to disable single-diagram enhancement (fix #899)
#endif
            fptype* allMEs,                     // output: allMEs[nevt], |M|^2 final_avg_over_helicities
            int* allselhel,                     // output: helicity selection[nevt]
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            int* allselcol,                     // output: helicity selection[nevt]
            fptype* allNumerators,              // tmp: numerators[nevt], running_sum_over_helicities
            fptype* allDenominators,            // tmp: denominators[nevt], running_sum_over_helicities
#endif
            const int nevt );                   // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif /* clang-format on */

  //--------------------------------------------------------------------------
}

#endif // MG5_Sigma_sm_uux_ttx_H
