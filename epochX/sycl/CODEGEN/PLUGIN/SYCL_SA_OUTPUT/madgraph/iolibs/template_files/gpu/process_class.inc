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
    //static const int ninitial = mgOnGpu::npari;
    //static const int nexternal = %(nexternal)d; // mgOnGpu::npar (nexternal was nioparticles)
    //static const int nprocesses = %(nprocesses)d; // FIXME: assume process.nprocesses == 1
    //static const int nwavefuncs = %(nwavefuncs)d; // mgOnGpu::nwf
    //static const int namplitudes = %(namp)d;
    //static const int ncomb = %(ncomb)d; // mgOnGpu::ncomb
    //static const int wrows = %(sizew)d; // mgOnGpu::nw6;

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
