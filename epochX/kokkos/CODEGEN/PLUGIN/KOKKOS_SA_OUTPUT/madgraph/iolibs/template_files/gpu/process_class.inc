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

  Kokkos::View<cxtype*> m_cIPC;
  typename Kokkos::View<cxtype*,Kokkos::HostSpace> m_hIPC;

  Kokkos::View<fptype*> m_cIPD;
  typename Kokkos::View<fptype*,Kokkos::HostSpace> m_hIPD;

private:

  // Command line arguments (constructor)
  bool m_verbose;
  bool m_debug;

    // Physics model parameters to be read from file (initProc function)
  Parameters_%(model_name)s* m_pars;
  std::vector<fptype> m_masses; // external particle masses

};
