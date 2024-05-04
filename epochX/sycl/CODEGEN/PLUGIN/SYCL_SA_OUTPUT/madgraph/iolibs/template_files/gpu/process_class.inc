! Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
! Created by: J. Alwall (Jul 2010) for the MG5aMC CPP backend.
!==========================================================================
! Copyright (C) 2020-2023 CERN and UCLouvain.
! Licensed under the GNU Lesser General Public License (version 3 or later).
! Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
! Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
!==========================================================================
! Copyright (C) 2021-2023 Argonne National Laboratory.
! Licensed under the GNU Lesser General Public License (version 3 or later).
! Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.
!==========================================================================

// --- Physics process-specific constants that are best declared at compile time
// dimensions of 4-momenta (E,px,py,pz)
#define CPPPROCESS_NP4 4

// #particles in the initial state (incoming): e.g. 2 (e+ e-) for e+ e- -> mu+ mu-
#define CPPPROCESS_NPARI %(nincoming)d

// #particles in the final state (outgoing): e.g. 2 (mu+ mu-) for e+ e- -> mu+ mu-
#define CPPPROCESS_NPARF %(noutcoming)d

// #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
#define CPPPROCESS_NPAR %(nexternal)d

// #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)
#define CPPPROCESS_NCOMB %(nbhel)d

// dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
#define CPPPROCESS_NW6 %(wavefuncsize)d

// #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
#define CPPPROCESS_NWF %(nwavefunc)d
  

namespace Proc
{
  //==========================================================================
  // A class for calculating the matrix elements for
  %(process_lines)s
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

  private:

    // Command line arguments (constructor)
    size_t m_numiterations; // number of iterations (each iteration has nblocks*nthreads events)
    size_t m_ngpublocks; // number of GPU blocks in one grid (i.e. one iteration)
    size_t m_ngputhreads; // number of GPU threads in a block
    bool m_verbose;
    bool m_debug;

#ifndef MGONGPU_HARDCODE_PARAM
    // Physics model parameters to be read from file (initProc function)
    Parameters_%(model_name)s* m_pars;
    std::vector<fptype> m_masses; // external particle masses

    // Physics parameters (masses, coupling, etc...)
    cxtype m_tIPC[independentCouplings::nicoup];
    fptype m_tIPD[mgOnGpu::nparams];
#endif

    // Other variables of this instance (???)
    //int id1, id2; // initial particle ids
    //cxtype** amp; // ???

  };
