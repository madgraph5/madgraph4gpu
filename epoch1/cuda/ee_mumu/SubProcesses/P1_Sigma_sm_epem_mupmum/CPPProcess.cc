//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "../../src/HelAmps_sm.cc"

#include <algorithm>
#include <iostream>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "CPPProcess.h"

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1

#ifdef __CUDACC__
namespace gProc
#else
namespace Proc
#endif
{
  using mgOnGpu::np4; // 4: the dimension of 4-momenta (E,px,py,pz)
  using mgOnGpu::npar; // number of particles in total (initial + final)
  using mgOnGpu::ncomb; // number of helicity combinations

  const int nwf = 5; // #wavefunctions: npar (4 external) + 1 (internal, reused for gamma and Z)
  const int nw6 = 6; // dimension of each wavefunction (see KEK 91-11)

#ifdef __CUDACC__
  __device__ __constant__ short cHel[ncomb][npar];
  __device__ __constant__ fptype cIPC[6];
  __device__ __constant__ fptype cIPD[2];
  //__device__ __constant__ int cNGoodHel[1]; // FIXME: assume process.nprocesses == 1 for the moment
  __device__ __constant__ int cNGoodHel;
  __device__ __constant__ int cGoodHel[ncomb];
#else
  static short cHel[ncomb][npar];
  static fptype cIPC[6];
  static fptype cIPD[2];
  //static int cNGoodHel[1]; // FIXME: assume process.nprocesses == 1 for the moment
  static int cNGoodHel;
  static int cGoodHel[ncomb];
#endif

  //--------------------------------------------------------------------------

  // Evaluate |M|^2 for each subprocess
  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum
  // of |M|^2 over helicities for the given event

  __device__
  void calculate_wavefunctions( int ihel,
                                const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                                fptype* allMEs               // output: allMEs[nevt], final |M|^2 averaged over all helicities
#ifndef __CUDACC__
                                , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                                )
  {
    using namespace MG5_sm;
    mgDebug( 0, __FUNCTION__ );
#ifndef __CUDACC__
    //printf( "calculate_wavefunctions: nevt %d\n", nevt );
#endif

#ifdef __CUDACC__
    const fptype cIPC[6] = { 0, -0.30795376724436879, 0, -0.28804415396362731, 0, 0.082309883272248419 };
    const fptype cIPD[2] = { 91.188000000000002, 2.4414039999999999 };
#endif

    // The number of colors
    const int ncolor = 1;

    // The color matrix
    const fptype denom[ncolor] = {1};
    const fptype cf[ncolor][ncolor] = {{1}};

#ifdef __CUDACC__
    // Local variables for the given event (ievt)
    cxtype w[nwf][nw6]; // w[5][6]
    cxtype amp[1]; // was 2
#else
    // Local variables for the given event page (ipagV)
    cxtype_sv w_v[nwf][nw6]; // w_v[5][6]
    cxtype_sv amp_v[1]; // was 2
#endif

    // For CUDA performance, this is ~better: fewer registers, even if no throughput increase (issue #39)
    // However, physics parameters like masses and couplings must be read from user parameter files
    //const fptype cIPC[6] = { 0, -0.30795376724436879, 0, -0.28804415396362731, 0, 0.082309883272248419 };
    //const fptype cIPD[2] = { 91.188000000000002, 2.4414039999999999 };

#ifndef __CUDACC__
    const int npagV = nevt / neppV;
    // ** START LOOP ON IPAGV **
    // - default(none): no variables are shared by default
    // - shared: as the name says
    // - private: give each thread its own copy, without initialising
    // - firstprivate: give each thread its own copy, and initialise with value from outside
#pragma omp parallel for default(none) shared(allmomenta,allMEs,cf,cHel,cIPC,cIPD,denom,ihel,npagV) private (amp_v,w_v)
    for ( int ipagV = 0; ipagV < npagV; ++ipagV )
#endif
    {
#ifdef __CUDACC__
      opzxxx( allmomenta, cHel[ihel][0], -1, w[0], 0 );
      //oxxxxx( allmomenta, 0, cHel[ihel][0], -1, w[0], 0 ); // tested ok (much slower)
#else
      opzxxx( allmomenta, cHel[ihel][0], -1, w_v[0], ipagV, 0 );
      //oxxxxx( allmomenta, 0, cHel[ihel][0], -1, w[0], ievt, 0 ); // tested ok (slower)
#endif

#ifdef __CUDACC__
      imzxxx( allmomenta, cHel[ihel][1], +1, w[1], 1 );
      //ixxxxx( allmomenta, 0, cHel[ihel][1], +1, w[1], 1 ); // tested ok (slower)
#else
      imzxxx( allmomenta, cHel[ihel][1], +1, w_v[1], ipagV, 1 );
      //ixxxxx( allmomenta, 0, cHel[ihel][1], +1, w[1], ievt, 1 ); // tested ok (a bit slower)
#endif

#ifdef __CUDACC__
      ixzxxx( allmomenta, cHel[ihel][2], -1, w[2], 2 );
      //ixxxxx( allmomenta, 0, cHel[ihel][2], -1, w[2], 2 ); // tested ok (a bit slower)
#else
      ixzxxx( allmomenta, cHel[ihel][2], -1, w_v[2], ipagV, 2 );
      //ixxxxx( allmomenta, 0, cHel[ihel][2], -1, w[2], ievt, 2 ); // tested ok (a bit slower)
#endif

#ifdef __CUDACC__
      oxzxxx( allmomenta, cHel[ihel][3], +1, w[3], 3 );
      //oxxxxx( allmomenta, 0, cHel[ihel][3], +1, w[3], 3 ); // tested ok (a bit slower)
#else
      oxzxxx( allmomenta, cHel[ihel][3], +1, w_v[3], ipagV, 3 );
      //oxxxxx( allmomenta, 0, cHel[ihel][3], +1, w[3], ievt, 3 ); // tested ok (a bit slower)
#endif

#ifdef __CUDACC__
      // Local variables for the given event (ievt)
      cxtype jamp[ncolor] = {}; // sum of the invariant amplitudes for all Feynman diagrams
#else
      // Local variables for the given event page (ipagV)
      cxtype_sv jamp_v[ncolor] = {}; // sum of the invariant amplitudes for all Feynman diagrams
#endif

#ifndef __CUDACC__
      FFV1P0_3( w_v[1], w_v[0], cxmake( cIPC[0], cIPC[1] ), 0., 0., w_v[4] );
      // Amplitude(s) for diagram number 1
      FFV1_0( w_v[2], w_v[3], w_v[4], cxmake( cIPC[0], cIPC[1] ), &amp_v[0] );
      jamp_v[0] -= amp_v[0];

      FFV2_4_3( w_v[1], w_v[0], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), cIPD[0], cIPD[1], w_v[4] );
      // Amplitude(s) for diagram number 2
      FFV2_4_0( w_v[2], w_v[3], w_v[4], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), &amp_v[0] );
      jamp_v[0] -= amp_v[0];

      // ** START LOOP ON IEPPV **
      for ( int ieppV = 0; ieppV < neppV; ++ieppV )
#endif
      {
#ifdef __CUDACC__
        FFV1P0_3( w[1], w[0], cxmake( cIPC[0], cIPC[1] ), 0., 0., w[4] );
        // Amplitude(s) for diagram number 1
        FFV1_0( w[2], w[3], w[4], cxmake( cIPC[0], cIPC[1] ), &amp[0] );
        jamp[0] -= amp[0];

        FFV2_4_3( w[1], w[0], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), cIPD[0], cIPD[1], w[4] );
        // Amplitude(s) for diagram number 2
        FFV2_4_0( w[2], w[3], w[4], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), &amp[0] );
        jamp[0] -= amp[0];
#endif

#ifdef __CUDACC__
        const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
        const int ievt = idim;
        //printf( "calculate_wavefunctions: ievt %d\n", ievt );
#else
        const int ievt = ipagV*neppV + ieppV;
        //printf( "calculate_wavefunctions: ievt %d\n", ievt );
#endif

#ifndef __CUDACC__
        // Local variables for the given event (ievt)
#ifdef MGONGPU_CPPSIMD
        cxtype jamp[ncolor];
        jamp[0] = jamp_v[0][ieppV];
#else
        cxtype* jamp = jamp_v;
#endif
#endif

        // Sum and square the color flows to get the matrix element
        // (compute |M|^2 by squaring |M|, taking into account colours)
        for( int icol = 0; icol < ncolor; icol++ )
        {
          cxtype ztemp = cxmake( 0, 0 );
          for( int jcol = 0; jcol < ncolor; jcol++ )
            ztemp += cf[icol][jcol] * jamp[jcol];
          // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event
          // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
          allMEs[ievt] += cxreal( ztemp * conj( jamp[icol] ) ) / denom[icol];
        }

        // Store the leading color flows for choice of color
        // for(i=0;i < ncolor; i++)
        // jamp2[0][i] += cxreal( jamp[i]*conj( jamp[i] ) );

        //printf( "calculate_wavefunction: %6d %2d %f\n", ievt, ihel, allMEs[ievt] );
      }
    }

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
    , m_ngpublocks( ngpublocks )
    , m_ngputhreads( ngputhreads )
    , m_verbose( verbose )
    , m_debug( debug )
    , m_pars( 0 )
    , m_masses()
  {
    // Helicities for the process - nodim
    const short tHel[ncomb][nexternal] =
      { {-1, -1, -1, -1}, {-1, -1, -1, +1}, {-1, -1, +1, -1}, {-1, -1, +1, +1},
        {-1, +1, -1, -1}, {-1, +1, -1, +1}, {-1, +1, +1, -1}, {-1, +1, +1, +1},
        {+1, -1, -1, -1}, {+1, -1, -1, +1}, {+1, -1, +1, -1}, {+1, -1, +1, +1},
        {+1, +1, -1, -1}, {+1, +1, -1, +1}, {+1, +1, +1, -1}, {+1, +1, +1, +1} };
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cHel, tHel, ncomb * nexternal * sizeof(short) ) );
#else
    memcpy( cHel, tHel, ncomb * nexternal * sizeof(short) );
#endif
    // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
    assert( sizeof(cxtype) == 2 * sizeof(fptype) );
#ifndef __CUDACC__
    // SANITY CHECK: momenta AOSOA uses vectors with the same size as fptype_v
    assert( neppV == mgOnGpu::neppM );
#endif
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

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
    const cxtype tIPC[3] = { cxmake( m_pars->GC_3 ), cxmake( m_pars->GC_50 ), cxmake( m_pars->GC_59 ) };
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
#ifdef __CUDACC__
#if defined __CUDACC_VER_MAJOR__ && defined __CUDACC_VER_MINOR__ && defined __CUDACC_VER_BUILD__
    out << "nvcc " << __CUDACC_VER_MAJOR__ << "." << __CUDACC_VER_MINOR__ << "." << __CUDACC_VER_BUILD__;
#else
    out << "nvcc UNKNOWN";
#endif
#elif defined __clang__
#if defined __clang_major__ && defined __clang_minor__ && defined __clang_patchlevel__
    out << "clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
#else
    out << "gcc UNKNOWKN";
#endif

#else
#if defined __GNUC__ && defined __GNUC_MINOR__ && defined __GNUC_PATCHLEVEL__
    out << "gcc (GCC) " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#else
    out << "gcc UNKNOWKN";
#endif
#endif
    return out.str();
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __global__
  void sigmaKin_getGoodHel( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype* allMEs,              // output: allMEs[nevt], final |M|^2 averaged over all helicities
                            bool* isGoodHel )            // output: isGoodHel[ncomb] - device array
  {
    const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
    const int ievt = idim;
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    fptype allMEsLast = 0;
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to running sum of |M|^2 over helicities for the given event(s)
      calculate_wavefunctions( ihel, allmomenta, allMEs );
      if ( allMEs[ievt] != allMEsLast )
      {
        //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
        isGoodHel[ihel] = true;
      }
      allMEsLast = allMEs[ievt]; // running sum up to helicity ihel for event ievt
    }
  }
#else
  void sigmaKin_getGoodHel( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype* allMEs,              // output: allMEs[nevt], final |M|^2 averaged over all helicities
                            bool* isGoodHel              // output: isGoodHel[ncomb] - device array
                            , const int nevt )           // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    const int maxtry0 = ( neppV > 10 ? neppV : 10 ); // 10, but at least neppV (otherwise the npagV loop does not even start)
    fptype allMEsLast[maxtry0] = { 0 };
    const int maxtry = std::min( maxtry0, nevt ); // 10, but at most nevt (avoid invalid memory access if nevt<maxtry0)
    for ( int ievt = 0; ievt < maxtry; ++ievt )
    {
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ievt] = 0; // all zeros
    }
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      calculate_wavefunctions( ihel, allmomenta, allMEs, maxtry );
      for ( int ievt = 0; ievt < maxtry; ++ievt )
      {
        // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
        if ( allMEs[ievt] != allMEsLast[ievt] )
        {
          //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
          isGoodHel[ihel] = true;
        }
        allMEsLast[ievt] = allMEs[ievt]; // running sum up to helicity ihel
      }
    }
  }
#endif

  //--------------------------------------------------------------------------

  void sigmaKin_setGoodHel( const bool* isGoodHel ) // input: isGoodHel[ncomb] - host array
  {
    // FIXME: assume process.nprocesses == 1 for the moment
    //int nGoodHel[1] = { 0 };
    int nGoodHel = 0;
    int goodHel[ncomb] = { 0 };
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      if ( isGoodHel[ihel] )
      {
        goodHel[nGoodHel] = ihel;
        nGoodHel++;
      }
    }
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cNGoodHel, &nGoodHel, sizeof(int) ) );
    checkCuda( cudaMemcpyToSymbol( cGoodHel, goodHel, ncomb*sizeof(int) ) );
#else
    cNGoodHel = nGoodHel;
    for ( int ihel = 0; ihel < ncomb; ihel++ ) cGoodHel[ihel] = goodHel[ihel];
#endif
  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour
  // FIXME: assume process.nprocesses == 1 (eventually: allMEs[nevt] -> allMEs[nevt*nprocesses]?)

  __global__
  void sigmaKin( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                 fptype* allMEs               // output: allMEs[nevt], final |M|^2 averaged over all helicities
#ifndef __CUDACC__
                 , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                 )
  {
    mgDebugInitialise();
    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    // pars->setDependentParameters();
    // pars->setDependentCouplings();
    // Reset color flows
    // start sigmakin_lines

    // Denominators: spins, colors and identical particles
    //const int nprocesses = 1;
    //const int denominators[nprocesses] = { 4 };
    const int denominators = 4;

#ifdef __CUDACC__
    // Remember: in CUDA this is a kernel for one event, in c++ this processes n events
    const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
    const int ievt = idim;
    //printf( "sigmakin: ievt %d\n", ievt );
#endif

    // PART 0 - INITIALISATION (before calculate_wavefunctions)

#ifndef __CUDACC__
    for ( int ievt = 0; ievt < nevt; ++ievt )
#endif
    {
      // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ievt] = 0; // all zeros
    }

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    for ( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      const int ihel = cGoodHel[ighel];
#ifdef __CUDACC__
      calculate_wavefunctions( ihel, allmomenta, allMEs );
#else
      calculate_wavefunctions( ihel, allmomenta, allMEs, nevt );
#endif
    }

    // PART 2 - FINALISATION (after calculate_wavefunctions)
#ifndef __CUDACC__
    for ( int ievt = 0; ievt < nevt; ++ievt )
#endif
    {
      // Get the final |M|^2 as an average over helicities/colors of running sum of |M|^2 over helicities for the given event
      // [NB 'sum over final spins, average over initial spins', eg see
      // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ievt] /= denominators;
    }
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------

} // end namespace
