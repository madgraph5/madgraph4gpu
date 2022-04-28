  //==========================================================================
  // A class for calculating the matrix elements for
  %(process_lines)s
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

  public: /* clang-format on */

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
    bool m_verbose;
    bool m_debug;

    // Physics model parameters to be read from file (initProc function)
#ifndef MGONGPU_HARDCODE_CIPD
    Parameters_%(model_name)s* m_pars;
#endif
    std::vector<fptype> m_masses; // external particle masses

    // Other variables of this instance (???)
    //int id1, id2; // initial particle ids
    //cxtype** amp; // ???
  };
