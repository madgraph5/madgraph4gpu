//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 2.9.5, 2021-08-22
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <algorithm>
#include <cstring>
#include <iostream>
#include <memory>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"
#include "HelAmps_sm.h"

#include "CPPProcess.h"

// Test ncu metrics for CUDA thread divergence
#undef MGONGPU_TEST_DIVERGENCE
//#define MGONGPU_TEST_DIVERGENCE 1

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1

#ifdef __CUDACC__
namespace gProc
#else
namespace Proc
#endif
{
  using mgOnGpu::np4; // dimensions of 4-momenta (E,px,py,pz)
  using mgOnGpu::npar; // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  using mgOnGpu::ncomb; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  using mgOnGpu::nwf; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  using mgOnGpu::nw6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)

#ifdef __CUDACC__
  __device__ __constant__ int cHel[ncomb][npar];
  __device__ __constant__ fptype cIPC[6];
  __device__ __constant__ fptype cIPD[2];
  __device__ __constant__ int cNGoodHel[1]; // FIXME: assume process.nprocesses == 1 for the moment
  __device__ __constant__ int cGoodHel[ncomb];
#else
  static int cHel[ncomb][npar];
  static fptype cIPC[6];
  static fptype cIPD[2];
  //static int cNGoodHel[1]; // FIXME: assume process.nprocesses == 1 for the moment
  static int cNGoodHel;
  static int cGoodHel[ncomb];
#endif

  //--------------------------------------------------------------------------

  // Evaluate |M|^2 for each subprocess
  // NB: calculate_wavefunctions ADDS |M|^2 for given ihel
  // to running sum of |M|^2 over helicities for given event(s)
  __device__ void calculate_wavefunctions( int ihel,
                                           const fptype* allmomenta,
                                           fptype& meHelSum
#ifndef __CUDACC__
                                           , const int ievt
#endif
                                           )
  {
    using namespace MG5_sm;
    mgDebug( 0, __FUNCTION__ );
    cxtype amp[1]; // was 2
    const int ncolor = 1;
    cxtype jamp[ncolor];

    // Calculate wavefunctions for all processes
    cxtype w[nwf][nw6];
    for( int i=0; i<1; i++ ){ jamp[i] = cxtype( 0., 0. ); } // reset jamp (reset color flows)

    // *** DIAGRAM 1 OF 2 ***

    // Wavefunction(s) for diagram number 1
#ifdef __CUDACC__
    opzxxx( allmomenta, cHel[ihel][0], -1, w[0], 0 );
#else
    opzxxx( allmomenta, cHel[ihel][0], -1, w[0], ievt, 0 );
#endif

#ifdef __CUDACC__
    imzxxx( allmomenta, cHel[ihel][1], +1, w[1], 1 );
#else
    imzxxx( allmomenta, cHel[ihel][1], +1, w[1], ievt, 1 );
#endif

#ifdef __CUDACC__
    ixzxxx( allmomenta, cHel[ihel][2], -1, w[2], 2 );
#else
    ixzxxx( allmomenta, cHel[ihel][2], -1, w[2], ievt, 2 );
#endif

#ifdef __CUDACC__
    oxzxxx( allmomenta, cHel[ihel][3], +1, w[3], 3 );
#else
    oxzxxx( allmomenta, cHel[ihel][3], +1, w[3], ievt, 3 );
#endif

    FFV1P0_3( w[1], w[0], cxtype( cIPC[0], cIPC[1] ), 0., 0., w[4] );

    // Amplitude(s) for diagram number 1
    FFV1_0( w[2], w[3], w[4], cxtype( cIPC[0], cIPC[1] ), &amp[0] );
    jamp[0] += -amp[0];

    // *** DIAGRAM 2 OF 2 ***

    // Wavefunction(s) for diagram number 2
    FFV2_4_3( w[1], w[0], cxtype( cIPC[2], cIPC[3] ), cxtype( cIPC[4], cIPC[5] ), cIPD[0], cIPD[1], w[4] );

    // Amplitude(s) for diagram number 2
    FFV2_4_0( w[2], w[3], w[4], cxtype( cIPC[2], cIPC[3] ), cxtype( cIPC[4], cIPC[5] ), &amp[0] );
    jamp[0] += -amp[0];

    // *** COLOR ALGEBRA BELOW ***
    // (This method used to be called CPPProcess::matrix_1_epem_mupmum()?)

    // The color matrix;
    static const fptype denom[ncolor] = {1};
    static const fptype cf[ncolor][ncolor] = {
      {1}};

    // Sum and square the color flows to get the matrix element
    for( int icol=0; icol < ncolor; icol++ )
    {
      cxtype ztemp = cxmake( 0, 0 );
      for( int jcol = 0; jcol < ncolor; jcol++ )
        ztemp = ztemp + cf[icol][jcol] * jamp[jcol];
      meHelSum = meHelSum + cxreal(ztemp*conj(jamp[icol])) / denom[icol];
    }

    // Store the leading color flows for choice of color
    //for( i=0; i < ncolor; i++ )
    //  jamp2[0][i] += real(jamp[i]*conj(jamp[i]));

    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess( int numiterations,
                          int ngpublocks,
                          int ngputhreads,
                          bool verbose,
                          bool debug )
    : m_numiterations( numiterations )
    , gpu_nblocks( ngpublocks )
    , gpu_nthreads( ngputhreads )
    , m_verbose( verbose )
    , m_debug( debug )
    , m_pars( 0 )
    , m_masses()
    , dim( gpu_nblocks * gpu_nthreads )
  {
    // Helicities for the process - nodim
    static const int tHel[ncomb][nexternal] = {
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
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cHel, tHel, ncomb * nexternal * sizeof(int) ) );
#else
    memcpy( cHel, tHel, ncomb * nexternal * sizeof(int) );
#endif
    // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
    assert( sizeof(cxtype) == 2 * sizeof(fptype) );
#ifndef __CUDACC__
    // SANITY CHECK: momenta AOSOA uses vectors with the same size as fptype_v
    //assert( neppV == mgOnGpu::neppM );
#endif
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

  // Initialize process
  void CPPProcess::initProc( const std::string& param_card_name )
  {
    // Instantiate the model class and set parameters that stay fixed during run
    m_pars = Parameters_sm::getInstance();
    SLHAReader slha( param_card_name, m_verbose );
    m_pars->setIndependentParameters( slha );
    m_pars->setIndependentCouplings();
    if ( m_verbose )
    {
      m_pars->printIndependentParameters();
      m_pars->printIndependentCouplings();
    }
    m_pars->setDependentParameters();
    m_pars->setDependentCouplings();

    // Set external particle masses for this matrix element
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );

    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    // Then copy them to CUDA constant memory (issue #39) or its C++ emulation in file-scope static memory
    const cxtype tIPC[3] = { cxmake(m_pars->GC_3), cxmake(m_pars->GC_50), cxmake(m_pars->GC_59) };
    const fptype tIPD[2] = { (fptype)m_pars->mdl_MZ, (fptype)m_pars->mdl_WZ };
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cIPC, tIPC, 3 * sizeof(cxtype) ) );
    checkCuda( cudaMemcpyToSymbol( cIPD, tIPD, 2 * sizeof(fptype) ) );
#else
    memcpy( cIPC, tIPC, 3 * sizeof(cxtype) );
    memcpy( cIPD, tIPD, 2 * sizeof(fptype) );
#endif

    //std::cout << std::setprecision(17) << "tIPC[0] = " << tIPC[0] << std::endl;
    //std::cout << std::setprecision(17) << "tIPC[1] = " << tIPC[1] << std::endl;
    //std::cout << std::setprecision(17) << "tIPC[2] = " << tIPC[2] << std::endl;
    //std::cout << std::setprecision(17) << "tIPD[0] = " << tIPD[0] << std::endl;
    //std::cout << std::setprecision(17) << "tIPD[1] = " << tIPD[1] << std::endl;
  }

  //--------------------------------------------------------------------------

  // Retrieve the compiler that was used to build this module
  const std::string CPPProcess::getCompiler()
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

#ifdef __CUDACC__
  __global__
  void sigmaKin_getGoodHel( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            bool* isGoodHel )         // output: isGoodHel[ncomb] - device array
  {
    const int nprocesses = 1; // FIXME: assume process.nprocesses == 1
    fptype meHelSum[nprocesses] = { 0 }; // all zeros
    fptype meHelSumLast = 0;
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running
      // sum of |M|^2 over helicities for the given event
      calculate_wavefunctions( ihel, allmomenta, meHelSum[0] );
      if ( meHelSum[0] != meHelSumLast )
      {
        isGoodHel[ihel] = true;
        meHelSumLast = meHelSum[0];
      }
    }
  }
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void sigmaKin_setGoodHel( const bool* isGoodHel ) // input: isGoodHel[ncomb] - host array
  {
    int nGoodHel[1] = { 0 }; // FIXME: assume process.nprocesses == 1 for the moment
    int goodHel[ncomb] = { 0 };
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      if ( isGoodHel[ihel] )
      {
        goodHel[nGoodHel[0]] = ihel;
        nGoodHel[0]++;
      }
    }
    checkCuda( cudaMemcpyToSymbol( cNGoodHel, nGoodHel, sizeof(int) ) );
    checkCuda( cudaMemcpyToSymbol( cGoodHel, goodHel, ncomb * sizeof(int) ) );
  }
#endif

  //--------------------------------------------------------------------------

  // Evaluate |M|^2, part independent of incoming flavour
  __global__
  void sigmaKin( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                 fptype_sv* allMEs            // output: allMEs[nevt], final |M|^2 averaged over all helicities
#ifndef __CUDACC__
                 , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                 )
  {
    mgDebugInitialise();

    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    //m_pars->setDependentParameters();
    //m_pars->setDependentCouplings();

#ifndef __CUDACC__
    const int maxtry = 10;
    static unsigned long long sigmakin_itry = 0; // first iteration over nevt events
    static bool sigmakin_goodhel[ncomb] = { false };
#endif

    // Start sigmaKin_lines
#ifndef __CUDACC__
    // ** START LOOP ON IEVT **
    for( int ievt = 0; ievt < nevt; ++ievt )
#endif
    {
#ifdef __CUDACC__
      const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
      const int ievt = idim;
      //printf( "sigmakin: ievt %d\n", ievt );
#endif

      // Denominators: spins, colors and identical particles
      const int nprocesses = 1; // FIXME: assume process.nprocesses == 1
      const int denominators[1] = { 4 };

      // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
      fptype meHelSum[nprocesses] = { 0 }; // all zeros

#ifdef __CUDACC__
      // CUDA - using precomputed good helicities
      for ( int ighel = 0; ighel < cNGoodHel[0]; ighel++ )
      {
        const int ihel = cGoodHel[ighel];
        calculate_wavefunctions( ihel, allmomenta, meHelSum[0] );
      }
#else
      // C++ - compute good helicities within this loop
      fptype meHelSumLast = 0; // check for good helicities
      for ( int ihel = 0; ihel < ncomb; ihel++ )
      {
        if ( sigmakin_itry > maxtry && !sigmakin_goodhel[ihel] ) continue;
        // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event
        calculate_wavefunctions( ihel, allmomenta, meHelSum[0], ievt );
        if ( sigmakin_itry <= maxtry )
        {
          if ( !sigmakin_goodhel[ihel] && meHelSum[0] > meHelSumLast ) sigmakin_goodhel[ihel] = true;
          meHelSumLast = meHelSum[0];
        }
      }
#endif

      // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
      // [NB 'sum over final spins, average over initial spins', eg see
      // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
      for( int iproc = 0; iproc < nprocesses; ++iproc )
      {
        meHelSum[iproc] /= denominators[iproc];
      }

      // Set the final average |M|^2 for this event in the output array for all events
      for( int iproc = 0; iproc < nprocesses; ++iproc )
      {
        allMEs[iproc*nprocesses + ievt] = meHelSum[iproc];
      }

#ifndef __CUDACC__
      if ( sigmakin_itry <= maxtry )
        sigmakin_itry++;
      //if ( sigmakin_itry == maxtry )
      //  for (int ihel = 0; ihel < ncomb; ihel++ )
      //    printf( "sigmakin: ihelgood %2d %d\n", ihel, sigmakin_goodhel[ihel] );
#endif
      // ** END LOOP ON IEVT **
    }
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------

} // end namespace

//==========================================================================

// This was initially added to both C++ and CUDA in order to avoid RDC in CUDA (issue #51)
// This is now also needed by C++ LTO-like optimizations via inlining (issue #229)
#include "HelAmps_sm.cc"

//==========================================================================

