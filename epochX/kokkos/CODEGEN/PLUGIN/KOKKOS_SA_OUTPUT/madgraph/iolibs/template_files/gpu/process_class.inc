//==========================================================================
// A class for calculating the matrix elements for
%(process_lines)s
//--------------------------------------------------------------------------

class CPPProcess
{
public:

  // Constructor (from command line arguments)
  CPPProcess( bool verbose = false, bool debug = false );


  // Destructor
  ~CPPProcess() = default;

  // Initialize process (read model parameters from file)
  virtual void initProc(const std::string& param_card_name);
  
  // Pointer accessors
  cxtype* get_tIPC_ptr();
  const cxtype* get_tIPC_ptr() const;
  
  fptype* get_tIPD_ptr();
  const fptype* get_tIPD_ptr() const;

private:

  // Command line arguments (constructor)
  bool m_verbose;
  bool m_debug;

  // Physics model parameters to be read from file (initProc function)
  Parameters_%(model_name)s* m_pars;
  std::vector<fptype> m_masses; // external particle masses

  // Physics parameters (masses, coupling, etc...)
  cxtype m_tIPC[independentCouplings::nicoup];
  fptype m_tIPD[mgOnGpu::nparams];

};
