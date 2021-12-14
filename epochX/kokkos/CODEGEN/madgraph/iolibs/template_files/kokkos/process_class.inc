//==========================================================================
// A class for calculating the matrix elements for
%(process_lines)s
//--------------------------------------------------------------------------

template <class ExecSpace>
class CPPProcess
{
public:

  CPPProcess(int numiterations, int leaguesize, int teamsize,
             bool verbose = false, bool debug = false): 
        m_numiterations(numiterations), league_size(leaguesize), 
        team_size(teamsize), 
        dim(league_size * team_size),
        cHel("cHel",ncomb,nexternal), hHel("hHel",ncomb,nexternal), 
        cmME("cmME",nexternal), hmME("hmME",nexternal),
        cIPC("cIPC",%(ncouplings)i), hIPC("hIPC",%(ncouplings)i), 
        cIPD("cIPD",%(nparams)i), hIPD("hIPD",%(nparams)i) {

    // Helicities for the process - nodim
    %(all_helicities)s

    for(int i=0;i<ncomb;++i)
        for(int j=0;j<nexternal;++j){
            hHel(i,j) = tHel[i][j];
        }
    Kokkos::deep_copy(cHel,hHel);
  }


  ~CPPProcess() = default;

  // Initialize process.
  virtual void initProc(std::string param_card_name){
    // Instantiate the model class and set parameters that stay fixed during run
    pars = Parameters_%(model_name)s::getInstance();
    SLHAReader slha(param_card_name, m_verbose);
    pars->setIndependentParameters(slha);
    pars->setIndependentCouplings();
    if (m_verbose) {
        pars->printIndependentParameters();
      pars->printIndependentCouplings();
    }
    pars->setDependentParameters();
    pars->setDependentCouplings();
    %(initProc_lines)s
    
    %(assign_coupling)s

  }


  virtual int code() const { return 1; }

  const typename Kokkos::View<double*,ExecSpace>::HostMirror &getMasses() const { return hmME; }

  void setInitial(int inid1, int inid2) {
    id1 = inid1;
    id2 = inid2;
  }

  int getDim() const { return dim; }

  int getNIOParticles() const { return nexternal; }


  // Constants for array limits
  static const int ninitial = mgKokkos::npari;
  static const int nexternal = mgKokkos::npar;
  static const int nprocesses = %(nprocesses)d;
  static const int ncomb = mgKokkos::ncomb;

  Kokkos::View<int**,ExecSpace> cHel;
  typename Kokkos::View<int**,ExecSpace>::HostMirror hHel;

  Kokkos::View<double*,ExecSpace> cmME;
  typename Kokkos::View<double*,ExecSpace>::HostMirror hmME;

  Kokkos::View<complex_t<double>*,ExecSpace> cIPC;
  typename Kokkos::View<complex_t<double>*,ExecSpace>::HostMirror hIPC;

  Kokkos::View<double*,ExecSpace> cIPD;
  typename Kokkos::View<double*,ExecSpace>::HostMirror hIPD;

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
