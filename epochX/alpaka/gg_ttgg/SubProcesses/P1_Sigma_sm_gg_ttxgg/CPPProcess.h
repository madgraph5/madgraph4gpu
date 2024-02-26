//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 2.9.5, 2021-08-22
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_gg_ttxgg_H
#define MG5_Sigma_sm_gg_ttxgg_H

#include <cassert>
#include <complex>
#include <iostream>
#include <vector>

#ifdef ALPAKA
#include <cuda_to_cupla.hpp>
#endif

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include "Parameters_sm.h"

//--------------------------------------------------------------------------

#if defined(ALPAKA) || defined(__CUDACC__)

#define checkCuda( code ) { assertCuda( code, __FILE__, __LINE__ ); }

inline void assertCuda( cudaError_t code, const char* file, int line, bool abort = true )
{
  if ( code != cudaSuccess )
  {
    printf( "GPUassert: %s %s:%d\n", cudaGetErrorString(code), file, line );
    if ( abort ) assert( code == cudaSuccess );
  }
}

#endif

//--------------------------------------------------------------------------

#ifdef ALPAKA
namespace gProc
{
  //==========================================================================
  // A class for calculating the matrix elements for
  // Process: g g > t t~ g g WEIGHTED<=4 @1
  //--------------------------------------------------------------------------

  class CPPProcess
  {
  public:

    // Constructor (from command line arguments)
    CPPProcess( int numiterations, int gpublocks, int gputhreads, bool verbose = false, bool debug = false );

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

    // Accessors (unused so far: add them to fix a clang build warning)
    //int numiterations() const { return m_numiterations; }
    //int gpublocks() const { return m_ngpublocks; }
    //int gputhreads() const { return m_ngputhreads; }
    //bool verbose() const { return m_verbose; }
    //bool debug() const { return m_debug; }

  public:

    // Hardcoded parameters for this process (constant class variables)
    //static const int ninitial = mgOnGpu::npari;
    //static const int nexternal = 6; // mgOnGpu::npar (nexternal was nioparticles)
    //static const int nprocesses = 1; // FIXME: assume process.nprocesses == 1
    //static const int nwavefuncs = 6; // mgOnGpu::nwf
    //static const int namplitudes = 159;
    //static const int ncomb = 64; // mgOnGpu::ncomb
    //static const int wrows = 63; // mgOnGpu::nw6;

  private:

    // Command line arguments (constructor)
    int m_numiterations; // number of iterations (each iteration has nblocks*nthreads events)
    int m_ngpublocks; // number of GPU blocks in one grid (i.e. one iteration)
    int m_ngputhreads; // number of GPU threads in a block
    bool m_verbose;
    bool m_debug;

    // Physics model parameters to be read from file (initProc function)
    Parameters_sm* m_pars;
    std::vector<fptype> m_masses; // external particle masses

    // Other variables of this instance (???)
    //int id1, id2; // initial particle ids
    //cxtype** amp; // ???

  };

  //--------------------------------------------------------------------------

  struct sigmaKin_getGoodHel
  {
    template< typename T_Acc >
    ALPAKA_FN_ACC
    void operator()( T_Acc const &acc,
                     const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                     fptype_sv* allMEs,           // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                     bool* isGoodHel              // output: isGoodHel[ncomb] - device array
                   ) const;
  };

  //--------------------------------------------------------------------------

  void sigmaKin_setGoodHel( const bool* isGoodHel ); // input: isGoodHel[ncomb] - host array

  //--------------------------------------------------------------------------

  struct sigmaKin
  {
    template< typename T_Acc >
    ALPAKA_FN_ACC
    void operator()( T_Acc const &acc,
                     const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                     fptype_sv* allMEs            // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                   ) const;
  };

  //--------------------------------------------------------------------------
}
#else
#ifdef __CUDACC__
namespace gProc
#else
namespace Proc
#endif
{

  //==========================================================================
  // A class for calculating the matrix elements for
  // Process: g g > t t~ g g WEIGHTED<=4 @1
  //--------------------------------------------------------------------------

  class CPPProcess
  {
  public:

    // Constructor (from command line arguments)
    CPPProcess( int numiterations, int gpublocks, int gputhreads, bool verbose = false, bool debug = false );

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

    // Accessors (unused so far: add them to fix a clang build warning)
    //int numiterations() const { return m_numiterations; }
    //int gpublocks() const { return m_ngpublocks; }
    //int gputhreads() const { return m_ngputhreads; }
    //bool verbose() const { return m_verbose; }
    //bool debug() const { return m_debug; }

  public:

    // Hardcoded parameters for this process (constant class variables)
    //static const int ninitial = mgOnGpu::npari;
    //static const int nexternal = 6; // mgOnGpu::npar (nexternal was nioparticles)
    //static const int nprocesses = 1; // FIXME: assume process.nprocesses == 1
    //static const int nwavefuncs = 6; // mgOnGpu::nwf
    //static const int namplitudes = 159;
    //static const int ncomb = 64; // mgOnGpu::ncomb
    //static const int wrows = 63; // mgOnGpu::nw6;

  private:

    // Command line arguments (constructor)
    int m_numiterations; // number of iterations (each iteration has nblocks*nthreads events)
    int m_ngpublocks; // number of GPU blocks in one grid (i.e. one iteration)
    int m_ngputhreads; // number of GPU threads in a block
    bool m_verbose;
    bool m_debug;

    // Physics model parameters to be read from file (initProc function)
    Parameters_sm* m_pars;
    std::vector<fptype> m_masses; // external particle masses

    // Other variables of this instance (???)
    //int id1, id2; // initial particle ids
    //cxtype** amp; // ???

  };

  //--------------------------------------------------------------------------

  __global__
  void sigmaKin_getGoodHel( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype_sv* allMEs,           // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                            bool* isGoodHel              // output: isGoodHel[ncomb] - device array
#ifndef __CUDACC__
                            , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                            );

  //--------------------------------------------------------------------------

  void sigmaKin_setGoodHel( const bool* isGoodHel ); // input: isGoodHel[ncomb] - host array

  //--------------------------------------------------------------------------

  __global__
  void sigmaKin( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                 fptype_sv* allMEs            // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
#ifndef __CUDACC__
                 , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                 );

  //--------------------------------------------------------------------------
}
#endif

#endif // MG5_Sigma_sm_gg_ttxgg_H
