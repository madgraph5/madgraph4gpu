//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_epem_mupmum_H
#define MG5_Sigma_sm_epem_mupmum_H

#define DPCT_USM_LEVEL_NONE
#include <CL/sycl.hpp>
#include <dpct/dpct.hpp>
#include <cassert>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "Parameters_sm.h"

//--------------------------------------------------------------------------

#ifdef CL_SYCL_LANGUAGE_VERSION

#define checkCuda( code )                       \
  { assertCuda( code, __FILE__, __LINE__ ); }

inline void assertCuda(int code, const char *file, int line, bool abort = true)
{
}

#endif

//--------------------------------------------------------------------------

#ifdef CL_SYCL_LANGUAGE_VERSION
namespace gProc
#else
namespace Proc
#endif
{

  //==========================================================================
  // A class for calculating the matrix elements for
  // Process: e+ e- > mu+ mu- WEIGHTED<=4 @1
  //--------------------------------------------------------------------------

  class CPPProcess
  {
  public:

    CPPProcess( int numiterations, int gpublocks, int gputhreads, bool verbose = false );

    ~CPPProcess();

    // Initialize process.
    virtual void initProc(std::string param_card_name);

    virtual int code() const {return 1;}

    const std::vector<fptype> &getMasses() const;

    void setInitial(int inid1, int inid2)
    {
      id1 = inid1;
      id2 = inid2;
    }

    int getDim() const {return dim;}

    int getNIOParticles() const {return nexternal;}

    // Constants for array limits
    static const int ninitial = mgOnGpu::npari;
    static const int nexternal = mgOnGpu::npar;
    //static const int nprocesses = 1; // FIXME: assume process.nprocesses == 1

  private:

    int m_numiterations;

    // gpu variables
    int gpu_nblocks;
    int gpu_nthreads;
    int dim; // gpu_nblocks * gpu_nthreads;

    bool m_verbose; // print verbose info

    static const int nwavefuncs = 6;
    static const int namplitudes = 2;
    static const int ncomb = 16;
    static const int wrows = 6;

    cxtype** amp;

    // Pointer to the model parameters
    Parameters_sm * pars;

    // vector with external particle masses
    std::vector<fptype> mME;

    // Initial particle ids
    int id1, id2;

  };

  //--------------------------------------------------------------------------

#ifdef CL_SYCL_LANGUAGE_VERSION
  SYCL_EXTERNAL
  void sigmaKin_getGoodHel( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            bool* isGoodHel,
                            sycl::nd_item<3> item_ct1,
                            dpct::accessor<int, dpct::device, 2> cHel );        // output: isGoodHel[ncomb] - device array
#endif

  //--------------------------------------------------------------------------

#ifdef CL_SYCL_LANGUAGE_VERSION
  void sigmaKin_setGoodHel( const bool* isGoodHel ); // input: isGoodHel[ncomb] - host array
#endif

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin(const fptype *allmomenta, // input: momenta as
                                          // AOSOA[npagM][npar][4][neppM] with
                                          // nevt=npagM*neppM
                fptype *allMEs, sycl::nd_item<3> item_ct1,
                dpct::accessor<int, dpct::device, 2> cHel, int *cNGoodHel,
                int *cGoodHel // output: allMEs[nevt], final |M|^2 averaged over
                              // all helicities
#ifndef CL_SYCL_LANGUAGE_VERSION
                ,
                const int nevt // input: #events (for cuda: nevt == ndim ==
                               // gpublocks*gputhreads)
#endif
  );

  //--------------------------------------------------------------------------

}

#endif  // MG5_Sigma_sm_epem_mupmum_H
