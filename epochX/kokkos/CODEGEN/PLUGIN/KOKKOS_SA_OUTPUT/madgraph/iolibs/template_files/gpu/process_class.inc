///////////////////////// -*- C++ -*- /////////////////////////////

//==========================================================================
// A class for calculating the matrix elements for
%(process_lines)s
//--------------------------------------------------------------------------

template <class ExecSpace>
class CPPProcess
{
public:

  CPPProcess(int numiterations, int leaguesize, int teamsize,
             bool verbose = false, bool debug = false);


  ~CPPProcess() = default;

  // Initialize process.
	virtual void initProc(const std::string& param_card_name);
  const std::string getCompiler();

  virtual int code() const { return 1; }

  const typename Kokkos::View<fptype*,ExecSpace>::HostMirror &getMasses() const { return hmME; }

  void setInitial(int inid1, int inid2) {
    id1 = inid1;
    id2 = inid2;
  }

  int getDim() const { return dim; }

  int getNIOParticles() const { return mgOnGpu::npar; }


  // Constants for array limits
  static const int ninitial = mgOnGpu::npari;
  static const int nexternal = mgOnGpu::npar;
  static const int nprocesses = %(nprocesses)d;
  static const int ncomb = mgOnGpu::ncomb;

  Kokkos::View<int**,ExecSpace> cHel;
  typename Kokkos::View<int**,ExecSpace>::HostMirror hHel;

  Kokkos::View<fptype*,ExecSpace> cmME;
  typename Kokkos::View<fptype*,ExecSpace>::HostMirror hmME;

  Kokkos::View<fptype*,ExecSpace> cIPC;
  typename Kokkos::View<fptype*,ExecSpace>::HostMirror hIPC;

  Kokkos::View<fptype*,ExecSpace> cIPD;
  typename Kokkos::View<fptype*,ExecSpace>::HostMirror hIPD;

  const short* get_tHel_ptr() const;

  cxtype* get_tIPC_ptr();
  const cxtype* get_tIPC_ptr() const;

  fptype* get_tIPD_ptr();
  const fptype* get_tIPD_ptr() const;
private:
  int m_numiterations;
  // gpu variables
  int league_size;
  int team_size;
  int dim; // league_size * team_size;

  // print verbose info
  bool m_verbose;

  // print debug info
  bool m_debug;

  static const int nwavefuncs = %(nwavefuncs)d;
  static const int namplitudes = %(namp)d;
  static const int wrows = %(sizew)d;
  //static const int nioparticles = %(nexternal)d;

  // Pointer to the model parameters
  Parameters_%(model)s *pars;

  // Initial particle ids
  int id1, id2;

};

// Constants for array limits
template <typename ExecSpace> const int CPPProcess<ExecSpace>::ninitial;
template <typename ExecSpace> const int CPPProcess<ExecSpace>::nprocesses;
template <typename ExecSpace> const int CPPProcess<ExecSpace>::ncomb;
template <typename ExecSpace> const int CPPProcess<ExecSpace>::nexternal;
