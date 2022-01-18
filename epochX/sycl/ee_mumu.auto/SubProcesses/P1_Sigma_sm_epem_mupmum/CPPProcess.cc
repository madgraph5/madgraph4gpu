//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 2.9.5, 2021-08-22
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <CL/sycl.hpp>
#include "../../src/HelAmps_sm.cc"

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <memory>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "CPPProcess.h"
#include "extras.h"

#undef MGONGPU_TEST_DIVERGENCE
//#define MGONGPU_TEST_DIVERGENCE 1

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1

namespace Proc
{
  using mgOnGpu::np4; // dimensions of 4-momenta (E,px,py,pz)
  using mgOnGpu::npar; // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  using mgOnGpu::ncomb; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  using mgOnGpu::nwf; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  using mgOnGpu::nw6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)

  // Physics parameters (masses, coupling, etc...)
#ifdef MGONGPU_HARDCODE_CIPC
  __device__ const fptype cIPC[6] = {
    Parameters_sm::GC_3.real(), Parameters_sm::GC_3.imag(),
    Parameters_sm::GC_50.real(), Parameters_sm::GC_50.imag(),
    Parameters_sm::GC_59.real(), Parameters_sm::GC_59.imag() };
  __device__ const fptype cIPD[2] = {
    Parameters_sm::mdl_MZ,
    Parameters_sm::mdl_WZ };
#else
  static fptype cIPC[6];
  static fptype cIPD[2];
#endif

  // Helicity combinations (and filtering of "good" helicity combinations)
  static short cHel[ncomb][npar];

  //--------------------------------------------------------------------------

    // Evaluate |M|^2 for each subprocess
    // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
    void calculate_wavefunctions( int ihel,
                                  const fptype *allmomenta,
                                  fptype *allMEs,
                                  size_t ievt,
                                  short *cHel,
                                  const fptype *cIPC,
                                  const fptype *cIPD
                                ) {
        using namespace MG5_sm;
        mgDebug( 0, __FUNCTION__ );

    // The number of colors
    constexpr int ncolor = 1;

    // Local TEMPORARY variables for a subset of Feynman diagrams in the given CUDA event (ievt) or C++ event page (ipagV)
    // [NB these variables are reused several times (and re-initialised each time) within the same event or event page]
    cxtype_sv w_sv[nwf][nw6]; // particle wavefunctions within Feynman diagrams (nw6 is often 6, the dimension of spin 1/2 or spin 1 particles)
    cxtype_sv amp_sv[1]; // invariant amplitude for one given Feynman diagram

    // Local variables for the given CUDA event (ievt) or C++ event page (ipagV)
    cxtype_sv jamp_sv[ncolor] = {}; // sum of the invariant amplitudes for all Feynman diagrams in the event or event page

    // === Calculate wavefunctions and amplitudes for all diagrams in all processes - Loop over nevt events ===

    {
      // Reset color flows (reset jamp_sv) at the beginning of a new event or event page
      for( int i=0; i<ncolor; i++ ){ jamp_sv[i] = cxzero_sv(); }

      // *** DIAGRAM 1 OF 2 ***

      // Wavefunction(s) for diagram number 1
#ifdef SYCL_LANGUAGE_VERSION
#ifndef MGONGPU_TEST_DIVERGENCE
      opzxxx( allmomenta, cHel[ihel][0], -1, w_sv[0], 0 ); // NB: opzxxx only uses pz
#else
      if ( ( blockDim.x * blockIdx.x + threadIdx.x ) % 2 == 0 )
        opzxxx( allmomenta, cHel[ihel][0], -1, w_sv[0], 0 ); // NB: opzxxx only uses pz
      else
        oxxxxx( allmomenta, 0, cHel[ihel][0], -1, w_sv[0], 0 );
#endif
#else
      opzxxx( p4IparIpagV( allmomenta, 0, ipagV ), cHel[ihel][0], -1, w_sv[0] ); // NB: opzxxx only uses pz
#endif

#ifdef SYCL_LANGUAGE_VERSION
      imzxxx( allmomenta, cHel[ihel][1], +1, w_sv[1], 1 ); // NB: imzxxx only uses pz
#else
      imzxxx( p4IparIpagV( allmomenta, 1, ipagV ), cHel[ihel][1], +1, w_sv[1] ); // NB: imzxxx only uses pz
#endif

#ifdef SYCL_LANGUAGE_VERSION
      ixzxxx( allmomenta, cHel[ihel][2], -1, w_sv[2], 2 );
#else
      ixzxxx( p4IparIpagV( allmomenta, 2, ipagV ), cHel[ihel][2], -1, w_sv[2] );
#endif

#ifdef SYCL_LANGUAGE_VERSION
      oxzxxx( allmomenta, cHel[ihel][3], +1, w_sv[3], 3 );
#else
      oxzxxx( p4IparIpagV( allmomenta, 3, ipagV ), cHel[ihel][3], +1, w_sv[3] );
#endif

      FFV1P0_3( w_sv[1], w_sv[0], cxmake( cIPC[0], cIPC[1] ), 0., 0., w_sv[4] );

      // Amplitude(s) for diagram number 1
      FFV1_0( w_sv[2], w_sv[3], w_sv[4], cxmake( cIPC[0], cIPC[1] ), &amp_sv[0] );
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 2 OF 2 ***

      // Wavefunction(s) for diagram number 2
      FFV2_4_3( w_sv[1], w_sv[0], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), cIPD[0], cIPD[1], w_sv[4] );

      // Amplitude(s) for diagram number 2
      FFV2_4_0( w_sv[2], w_sv[3], w_sv[4], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), &amp_sv[0] );
      jamp_sv[0] -= amp_sv[0];

      // *** COLOR ALGEBRA BELOW ***
      // (This method used to be called CPPProcess::matrix_1_epem_mupmum()?)

      // The color matrix [NB do keep 'static' for these constexpr arrays, see issue #283]
      static constexpr fptype denom[ncolor] = {1};
      static constexpr fptype cf[ncolor][ncolor] = {{1}};

      // Sum and square the color flows to get the matrix element
      // (compute |M|^2 by squaring |M|, taking into account colours)
      fptype_sv deltaMEs = { 0 }; // all zeros
      for( int icol = 0; icol < ncolor; icol++ )
      {
        cxtype_sv ztemp_sv = cxzero_sv();
        for( int jcol = 0; jcol < ncolor; jcol++ )
          ztemp_sv += cf[icol][jcol] * jamp_sv[jcol];
        deltaMEs += cxreal( ztemp_sv * cxconj( jamp_sv[icol] ) ) / denom[icol];
      }

      // *** STORE THE RESULTS ***

      // Store the leading color flows for choice of color
      // (NB: jamp2_sv must be an array of fptype_sv)
      // for( int icol = 0; icol < ncolor; icol++ )
      // jamp2_sv[0][icol] += cxreal( jamp_sv[icol]*cxconj( jamp_sv[icol] ) );

      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ievt] += deltaMEs;
      //if ( cNGoodHel > 0 ) printf( "calculate_wavefunction: %6d %2d %f\n", ievt, ihel, allMEs[ievt] );
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess(
          int numiterations,
          int ngpublocks,
          int ngputhreads,
          bool verbose,
          bool debug )
    : m_numiterations( numiterations )
    , m_ngpublocks( ngpublocks )
    , m_ngputhreads( ngputhreads )
    , m_verbose( verbose )
    , m_debug( debug )
#ifndef MGONGPU_HARDCODE_CIPC
    , m_pars( 0 )
#endif
    , m_masses()
  {
    // Helicities for the process [NB do keep 'static' for this constexpr array, see issue #283]
    static constexpr short tHel[ncomb][mgOnGpu::npar] = {
      {-1, -1, -1, -1},
      {-1, -1, -1, 1},
      {-1, -1, 1, -1},
      {-1, -1, 1, 1},
      {-1, 1, -1, -1},
      {-1, 1, -1, 1},
      {-1, 1, 1, -1},
      {-1, 1, 1, 1},
      {1, -1, -1, -1},
      {1, -1, -1, 1},
      {1, -1, 1, -1},
      {1, -1, 1, 1},
      {1, 1, -1, -1},
      {1, 1, -1, 1},
      {1, 1, 1, -1},
      {1, 1, 1, 1}};
    memcpy( *cHel, *tHel, ncomb * mgOnGpu::npar * sizeof(short) );

    // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
    static_assert( sizeof(cxtype) == 2 * sizeof(fptype) );
    // SANITY CHECK: check that neppR, neppM and neppV are powers of two (https://stackoverflow.com/a/108360)
    auto ispoweroftwo = []( int n ) { return ( n > 0 ) && !( n & ( n - 1 ) ); };
    static_assert( ispoweroftwo( mgOnGpu::neppR ) );
    static_assert( ispoweroftwo( mgOnGpu::neppM ) );
    static_assert( ispoweroftwo( neppV ) );
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

  //const std::vector<fptype>& CPPProcess::getMasses() const { return mME; }

  //--------------------------------------------------------------------------

  short * CPPProcess::get_cHel_ptr() const {return *cHel;}

  //--------------------------------------------------------------------------

  fptype* CPPProcess::get_cIPC_ptr() {return cIPC;}

  //--------------------------------------------------------------------------

  fptype* CPPProcess::get_cIPD_ptr() {return cIPD;}

  //--------------------------------------------------------------------------

#ifndef MGONGPU_HARDCODE_CIPC
  // Initialize process (with parameters read from user cards)
  void CPPProcess::initProc( const std::string& param_card_name )
  {
    // Instantiate the model class and set parameters that stay fixed during run
    m_pars = Parameters_sm::getInstance();
    SLHAReader slha( param_card_name, m_verbose );
    m_pars->setIndependentParameters( slha );
    m_pars->setIndependentCouplings();
    m_pars->setDependentParameters();
    m_pars->setDependentCouplings();
    if ( m_verbose )
    {
      m_pars->printIndependentParameters();
      m_pars->printIndependentCouplings();
      m_pars->printDependentParameters();
      m_pars->printDependentCouplings();
    }
    // Set external particle masses for this matrix element
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    // Then copy them to CUDA constant memory (issue #39) or its C++ emulation in file-scope static memory
    const cxtype tIPC[3] = { cxmake( m_pars->GC_3 ), cxmake( m_pars->GC_50 ), cxmake( m_pars->GC_59 ) };
    const fptype tIPD[2] = { (fptype)m_pars->mdl_MZ, (fptype)m_pars->mdl_WZ };
    memcpy( cIPC, tIPC, 3 * sizeof(cxtype) );
    memcpy( cIPD, tIPD, 2 * sizeof(fptype) );
  }
#else
  // Initialize process (with hardcoded parameters)
  void CPPProcess::initProc( const std::string& /*param_card_name*/ )
  {
    // Use hardcoded physics parameters
    if ( m_verbose )
    {
      Parameters_sm::printIndependentParameters();
      Parameters_sm::printIndependentCouplings();
      Parameters_sm::printDependentParameters();
      Parameters_sm::printDependentCouplings();
    }
    // Set external particle masses for this matrix element
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
  }
#endif

  //--------------------------------------------------------------------------

  // Retrieve the compiler that was used to build this module
  const std::string CPPProcess::getCompiler() {
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
    tchainout.pop_back(); // remove trailing newline
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

  SYCL_EXTERNAL
  void sigmaKin_getGoodHel( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype* allMEs,           // output: allMEs[nevt], |M|^2 final_avg_over_helicities
                            bool * isGoodHel,  // output: isGoodHel[ncomb] - device array
                            size_t ievt,
                            short *cHel,
                            const fptype *cIPC,
                            const fptype *cIPD
                            ) {
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      fptype allMEsLast = 0;
      for ( int ihel = 0; ihel < ncomb; ihel++ )
      {
        // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
        calculate_wavefunctions( ihel, allmomenta, allMEs, ievt, cHel, cIPC, cIPD );
        if ( allMEs[ievt] != allMEsLast )
        {
          //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
          isGoodHel[ihel] = true;
        }
        allMEsLast = allMEs[ievt]; // running sum up to helicity ihel for event ievt
      }
  }

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin_setGoodHel( const bool* isGoodHel, int* cNGoodHel, int* cGoodHel) {
    *cNGoodHel = 0;
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      if ( isGoodHel[ihel] ) {
        cGoodHel[*cNGoodHel] = ihel;
        (*cNGoodHel)++;
      }
    }
  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour
  // FIXME: assume process.nprocesses == 1 (eventually: allMEs[nevt] -> allMEs[nevt*nprocesses]?)

  SYCL_EXTERNAL
  void sigmaKin(
          const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
          fptype* allMEs,            // output: allMEs[nevt], |M|^2 final_avg_over_helicities
          size_t ievt,
          short *cHel,
          const fptype *cIPC,
          const fptype *cIPD,
          int *cNGoodHel,
          int *cGoodHel
          )  {
    mgDebugInitialise();

    // Denominators: spins, colors and identical particles
    const int denominators = 4; // FIXME: assume process.nprocesses == 1 for the moment (eventually denominators[nprocesses]?)

    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    //m_pars->setDependentParameters();
    //m_pars->setDependentCouplings();

    // Start sigmaKin_lines
    // PART 0 - INITIALISATION (before calculate_wavefunctions)
    // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    allMEs[ievt] = 0;

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    for ( int ighel = 0; ighel < *cNGoodHel; ighel++ ) {
        const int ihel = cGoodHel[ighel];
        calculate_wavefunctions( ihel, allmomenta, allMEs, ievt, cHel, cIPC, cIPD );
    }

    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    allMEs[ievt] /= denominators;
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------

} // end namespace


