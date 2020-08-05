#include "mgOnGpuConfig.h"
using mgOnGpu::dcomplex;

//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

/*
#ifndef HelAmps_sm_H
#define HelAmps_sm_H

#include <cmath>

namespace MG5_sm
{
  __device__ void oxxxxx(const double p[4],
                         const double fmass,
                         const int nhel,
                         const int nsf,
                         dcomplex fo[6]);

  __device__ void sxxxxx(const double p[4],
                         const int nss,
                         dcomplex sc[3]);

  __device__ void ixxxxx(const double p[4],
                         const double fmass,
                         const int nhel,
                         const int nsf,
                         dcomplex fi[6]);

  __device__ void txxxxx(const double p[4],
                         const double tmass,
                         const int nhel,
                         const int nst,
                         dcomplex fi[18]);

  __device__ void vxxxxx(const double p[4],
                         const double vmass,
                         const int nhel,
                         const int nsv,
                         dcomplex v[6]);

  __device__ void FFV2_0(const dcomplex F1[],
                         const dcomplex F2[],
                         const dcomplex V3[],
                         const dcomplex COUP,
                         dcomplex * vertex);

  __device__ void FFV2_3(const dcomplex F1[],
                         const dcomplex F2[],
                         const dcomplex COUP,
                         const double M3,
                         const double W3,
                         dcomplex V3[]);

  __device__ void FFV4_0(const dcomplex F1[],
                         const dcomplex F2[],
                         const dcomplex V3[],
                         const dcomplex COUP,
                         dcomplex * vertex);

  __device__ void FFV4_3(const dcomplex F1[],
                         const dcomplex F2[],
                         const dcomplex COUP,
                         const double M3,
                         const double W3,
                         dcomplex V3[]);

  __device__ void FFV1_0(const dcomplex F1[],
                         const dcomplex F2[],
                         const dcomplex V3[],
                         const dcomplex COUP,
                         dcomplex * vertex);

  __device__ void FFV1P0_3(const dcomplex F1[],
                           const dcomplex F2[],
                           const dcomplex COUP,
                           const double M3,
                           const double W3,
                           dcomplex V3[]);

  __device__ void FFV2_4_0(const dcomplex F1[],
                           const dcomplex F2[],
                           const dcomplex V3[],
                           const dcomplex COUP1,
                           const dcomplex COUP2,
                           dcomplex * vertex);

  __device__ void FFV2_4_3(const dcomplex F1[],
                           const dcomplex F2[],
                           const dcomplex COUP1,
                           const dcomplex COUP2,
                           const double M3,
                           const double W3,
                           dcomplex V3[]);

}  // end namespace MG5_sm

#endif  // HelAmps_sm_H
*/

//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_epem_mupmum_H
#define MG5_Sigma_sm_epem_mupmum_H

#include <cassert>
#include <complex>
#include <vector>

#include "Parameters_sm.h"

#ifdef __CUDACC__
#define gpuErrchk3( code )                      \
  { gpuAssert3( code, __FILE__, __LINE__ ); }

inline void gpuAssert3( cudaError_t code, const char *file, int line, bool abort = true )
{
  if ( code != cudaSuccess )
  {
    printf( "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line );
    if ( abort ) assert( code == cudaSuccess );
  }
}
#endif

#ifdef __CUDACC__
__global__ 
#endif
void sigmaKin( const double* allmomenta, double* output );

//==========================================================================
// A class for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1
//--------------------------------------------------------------------------

class CPPProcess
{
public:

  CPPProcess(int numiterations, int gpublocks, int gputhreads,
             bool verbose = false, bool debug = false);

  ~CPPProcess();

  // Initialize process.
  virtual void initProc(std::string param_card_name);


  virtual int code() const {return 1;}

  const std::vector<double> &getMasses() const;

  void setInitial(int inid1, int inid2)
  {
    id1 = inid1;
    id2 = inid2;
  }

  int getDim() const {return dim;}

  int getNIOParticles() const {return nexternal;}


  // Constants for array limits
  static const int ninitial = 2;
  static const int nexternal = 4;
  static const int nprocesses = 1;

private:
  int m_numiterations;
  // gpu variables
  int gpu_nblocks;
  int gpu_nthreads;
  int dim;  // gpu_nblocks * gpu_nthreads;

  // print verbose info
  bool m_verbose;

  // print debug info
  bool m_debug;

  static const int nwavefuncs = 6;
  static const int namplitudes = 2;
  static const int ncomb = 16;
  static const int wrows = 6;
  // static const int nioparticles = 4;

  dcomplex * * amp;


  // Pointer to the model parameters
  Parameters_sm * pars;

  // vector with external particle masses
  std::vector<double> mME;

  // Initial particle ids
  int id1, id2;

};


#endif  // MG5_Sigma_sm_epem_mupmum_H
