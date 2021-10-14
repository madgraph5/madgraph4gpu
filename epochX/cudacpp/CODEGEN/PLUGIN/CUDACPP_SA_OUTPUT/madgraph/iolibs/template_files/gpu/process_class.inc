  //==========================================================================
  // A class for calculating the matrix elements for
  %(process_lines)s
  //--------------------------------------------------------------------------

  class CPPProcess
  {
  public:

    CPPProcess( int numiterations,
                int gpublocks,
                int gputhreads,
                bool verbose = false,
                bool debug = false );

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

    // Accessors (unused so far: add them to fix a clang build warning)
    //int numiterations() const { return m_numiterations; }
    //int gpublocks() const { return m_ngpublocks; }
    //int gputhreads() const { return m_ngputhreads; }
    //bool verbose() const { return m_verbose; }
    //bool debug() const { return m_debug; }

    // Hardcoded parameters for this process (constant class variables)
    static const int ninitial = mgOnGpu::npari;
    static const int nexternal = mgOnGpu::npar;
    //static const int nprocesses = %(nprocesses)d; // FIXME: assume process.nprocesses == 1

  private:

    int m_numiterations;
    // gpu variables
    int gpu_nblocks;
    int gpu_nthreads;
    int dim; // gpu_nblocks * gpu_nthreads;

    // print verbose info
    bool m_verbose;

    // print debug info
    bool m_debug;

    static const int nwavefuncs = %(nwavefuncs)d;
    static const int namplitudes = %(namp)d;
    static const int ncomb = %(ncomb)d;
    static const int wrows = %(sizew)d;
    //static const int nioparticles = %(nexternal)d;

    cxtype** amp;

    // Physics model parameters to be read from file (initProc function)
    Parameters_sm* m_pars;
    std::vector<fptype> m_masses; // external particle masses

    // Initial particle ids
    int id1, id2;

  };
