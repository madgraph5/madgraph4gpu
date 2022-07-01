///////////////////////// -*- C++ -*- /////////////////////////////
//==========================================================================
// Class member functions for calculating the matrix elements for
%(process_lines)s

using mgOnGpu::ncomb; // number of helicity combinations
//--------------------------------------------------------------------------

// Evaluate |M|^2 for each subprocess
// NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event


%(all_sigmaKin)s

//--------------------------------------------------------------------------

template <typename ExecSpace>
CPPProcess<ExecSpace>::CPPProcess(
    int numiterations, int leaguesize, int teamsize,
    bool verbose, bool debug): 
      m_numiterations(numiterations), league_size(leaguesize), 
      team_size(teamsize), 
      dim(league_size * team_size),
      cHel("cHel",mgOnGpu::ncomb,mgOnGpu::npar), hHel("hHel",mgOnGpu::ncomb,mgOnGpu::npar), 
      cmME("cmME",mgOnGpu::npar), hmME("hmME",mgOnGpu::npar),
      cIPC("cIPC",%(ncouplingstimes2)s), hIPC("hIPC",%(ncouplingstimes2)s), 
      cIPD("cIPD",%(nparams)s), hIPD("hIPD",%(nparams)s)
{

  // Helicities for the process - nodim
  %(all_helicities)s

  for(int i=0;i<mgOnGpu::ncomb;++i)
      for(int j=0;j<mgOnGpu::npar;++j){
          hHel(i,j) = tHel[i][j];
      }
  Kokkos::deep_copy(cHel,hHel);
}

//--------------------------------------------------------------------------

// Initialize process (with parameters read from user cards)
template <typename ExecSpace>
void CPPProcess<ExecSpace>::initProc( const std::string& param_card_name )
{
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


// Retrieve the compiler that was used to build this module
template<typename ExecSpace>
const std::string CPPProcess<ExecSpace>::getCompiler()
{
  std::stringstream out;
  // CUDA version (NVCC)
#ifdef __CUDACC__
#if defined __CUDACC_VER_MAJOR__ && defined __CUDACC_VER_MINOR__ && defined __CUDACC_VER_BUILD__
  out << "nvcc " << __CUDACC_VER_MAJOR__ << "." << __CUDACC_VER_MINOR__ << "." << __CUDACC_VER_BUILD__;
#else
  out << "nvcc UNKNOWN";
#endif
  out << " (";
#endif
  // ICX version (either as CXX or as host compiler inside NVCC)
#if defined __INTEL_COMPILER
#error "icc is no longer supported: please use icx"
#elif defined __INTEL_LLVM_COMPILER // alternative: __INTEL_CLANG_COMPILER
  out << "icx " << __INTEL_LLVM_COMPILER << " (";
#endif
  // CLANG version (either as CXX or as host compiler inside NVCC or inside ICX)
#if defined __clang__
#if defined __clang_major__ && defined __clang_minor__ && defined __clang_patchlevel__
  out << "clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
  // GCC toolchain version inside CLANG
  std::string tchainout;
  std::string tchaincmd = "readelf -p .comment $(${CXX} -print-libgcc-file-name) |& grep 'GCC: (GNU)' | grep -v Warning | sort -u | awk '{print $5}'";
  std::unique_ptr<FILE, decltype(&pclose)> tchainpipe( popen( tchaincmd.c_str(), "r" ), pclose );
  if ( !tchainpipe ) throw std::runtime_error( "`readelf ...` failed?" );
  std::array<char, 128> tchainbuf;
  while ( fgets( tchainbuf.data(), tchainbuf.size(), tchainpipe.get() ) != nullptr ) tchainout += tchainbuf.data();
  if(tchainout.size() > 0) tchainout.pop_back(); // remove trailing newline
#if defined __CUDACC__ or defined __INTEL_LLVM_COMPILER
  out << ", gcc " << tchainout;
#else
  out << " (gcc " << tchainout << ")";
#endif
#else
  out << "clang UNKNOWKN";
#endif
#else
  // GCC version (either as CXX or as host compiler inside NVCC)
#if defined __GNUC__ && defined __GNUC_MINOR__ && defined __GNUC_PATCHLEVEL__
  out << "gcc " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#else
  out << "gcc UNKNOWKN";
#endif
#endif
#if defined __CUDACC__ or defined __INTEL_LLVM_COMPILER
  out << ")";
#endif
  return out.str();
}

//--------------------------------------------------------------------------
// Define pointer accessors
template <typename ExecSpace>
const short* CPPProcess<ExecSpace>::get_tHel_ptr() const {return &(**hHel);}

template <typename ExecSpace>
cxtype* CPPProcess<ExecSpace>::get_tIPC_ptr() {return hIPC;}
template <typename ExecSpace>
const cxtype* CPPProcess<ExecSpace>::get_tIPC_ptr() const {return hIPC;}

template <typename ExecSpace>
fptype* CPPProcess<ExecSpace>::get_tIPD_ptr() {return hIPD;}
template <typename ExecSpace>
const fptype* CPPProcess<ExecSpace>::get_tIPD_ptr() const {return hIPD;}

//--------------------------------------------------------------------------
template <typename mom_t, typename out_t, typename hel_t, typename ipd_t, 
          typename ipc_t, typename igh_t, typename ngh_t>
void sigmaKin_setup(
    const mom_t& momenta,
    out_t& allMEs,
    const hel_t& cHel,
    const ipd_t& cIPD,
    const ipc_t& cIPC,
    igh_t& iGoodHel,
    ngh_t& nGoodHel,
    const int& ncomb,
    const int& league_size,
    const int& team_size) 
{
  Kokkos::View<int*,Kokkos::DefaultExecutionSpace> isGoodHel("isGoodHel",ncomb); // has to be constant index, but should be `ncomb`

  using member_type = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy( league_size, team_size );
  Kokkos::parallel_for(__func__,policy, 
  KOKKOS_LAMBDA(const member_type& team_member){
    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();
    fptype allMEsLast = 0;

    auto local_mom = Kokkos::subview(momenta,ievt,Kokkos::ALL,Kokkos::ALL);
    for (int ihel = 0;ihel < ncomb;++ihel)
    {
      auto local_cHel = Kokkos::subview(cHel,ihel,Kokkos::ALL);
      calculate_wavefunctions(local_mom, local_cHel, cIPD, cIPC, allMEs[ievt]);
      if (allMEs[ievt] != allMEsLast)
      {
        isGoodHel(ihel) = true;
      }
      allMEsLast = allMEs[ievt];
    }
  });

  Kokkos::parallel_for(__func__,Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0,1),
  KOKKOS_LAMBDA(const int& i){
    for(int ihel=0;ihel < ncomb;++ihel){

      if(isGoodHel(ihel)){
        iGoodHel(nGoodHel(0)) = ihel;
        nGoodHel(0)++;
      }

    }
  });
}


//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.
template <typename mom_t, typename out_t, typename hel_t, typename ipd_t, 
          typename ipc_t, typename igh_t, typename ngh_t>
void sigmaKin(const mom_t& momenta, out_t& allMEs, const hel_t& cHel,
    const ipd_t& cIPD, const ipc_t& cIPC, const igh_t& iGoodHel,
    const ngh_t& nGoodHel, const int& ncomb, const int& league_size,
    const int& team_size)
{
  using member_type = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy( league_size, team_size );
  Kokkos::parallel_for(__func__,policy, 
  KOKKOS_LAMBDA(const member_type& team_member){

    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();
    
    // Denominators: spins, colors and identical particles
    constexpr int denominators = 512;
    allMEs[ievt] = 0;

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    auto local_mom = Kokkos::subview(momenta,ievt,Kokkos::ALL,Kokkos::ALL);
    for (int ighel = 0;ighel < nGoodHel(0);++ighel)
    {
      auto local_cHel = Kokkos::subview(cHel,iGoodHel(ighel),Kokkos::ALL);
      calculate_wavefunctions(local_mom, local_cHel, cIPD, cIPC, allMEs[ievt]);
    }
    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    allMEs[ievt] /= (fptype)denominators;

  });// end parallel for

}


//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
