//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>

#include "mgOnGpuConfig.h"
using mgOnGpu::dcomplex;
using mgOnGpu::double_v;
using mgOnGpu::dcomplex_v;

namespace MG5_sm
{

#ifndef __CUDACC__
  // Quick and dirty way to share nevt across all computational kernels
  int nevt;
#endif


#ifdef __CUDACC__
  __device__
#endif
  inline const double& pIparIp4Ievt( const double* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                                     const int ipar,
                                     const int ip4,
                                     const int ievt )
  {
    using mgOnGpu::np4;
    using mgOnGpu::npar;
#if defined MGONGPU_LAYOUT_ASA
    using mgOnGpu::nepp;
    const int ipag = ievt/nepp; // #eventpage in this iteration
    const int iepp = ievt%nepp; // #event in the current eventpage in this iteration
    // ASA: allmomenta[npag][npar][np4][nepp]
    return allmomenta[ipag*npar*np4*nepp + ipar*nepp*np4 + ip4*nepp + iepp]; // AOSOA[ipag][ipar][ip4][iepp]
#elif defined MGONGPU_LAYOUT_SOA
#ifdef __CUDACC__
    const int ndim = blockDim.x * gridDim.x; // (previously was: DIM)
    const int nevt = ndim;
#else
    using MG5_sm::nevt;
#endif
    // SOA: allmomenta[npar][np4][ndim]
    return allmomenta[ipar*np4*nevt + ip4*nevt + ievt]; // SOA[ipar][ip4][ievt]
#elif defined MGONGPU_LAYOUT_AOS
    // AOS: allmomenta[ndim][npar][np4]
    return allmomenta[ievt*npar*np4 + ipar*np4 + ip4]; // AOS[ievt][ipar][ip4]
#endif
  }


#ifdef __CUDACC__
  __device__
#endif
  void imzxxxM0( const double* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                 //const double fmass,
                 const int nhel,
                 const int nsf,
#ifndef __CUDACC__
                 dcomplex fi[6],
                 const int ievt,
#else
                 dcomplex_v fiv[6],
#endif
                 const int ipar )          // input: particle# out of npar
  {
#ifndef __CUDACC__
    // ** START LOOP ON IEVT **
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef __CUDACC__
      const int ieib = threadIdx.x; // event in block
      const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
      const int ievt = idim;
      //printf( "imzxxxM0: ievt=%d ieib=%d\n", ievt, ieib );
#endif
      const double& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const double& pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const double& pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const double& pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#ifndef __CUDACC__
      fi[0] = dcomplex (-pvec0 * nsf, -pvec3 * nsf);
      fi[1] = dcomplex (-pvec1 * nsf, -pvec2 * nsf);
#else
      fiv[0][ieib] = dcomplex (-pvec0 * nsf, -pvec3 * nsf);
      fiv[1][ieib] = dcomplex (-pvec1 * nsf, -pvec2 * nsf);
#endif
      const int nh = nhel * nsf;
      // ASSUMPTIONS FMASS = 0 and
      // (PX = PY = 0 and E = -P3 > 0)
      {
        const dcomplex chi0( 0, 0 );
        const dcomplex chi1( -nhel * sqrt(2 * pvec0), 0 );
        if (nh == 1)
        {
#ifndef __CUDACC__
          fi[2] = dcomplex (0, 0);
          fi[3] = dcomplex (0, 0);
          fi[4] = chi0;
          fi[5] = chi1;
#else
          fiv[2][ieib] = dcomplex (0, 0);
          fiv[3][ieib] = dcomplex (0, 0);
          fiv[4][ieib] = chi0;
          fiv[5][ieib] = chi1;
#endif
        }
        else
        {
#ifndef __CUDACC__
          fi[2] = chi1;
          fi[3] = chi0;
          fi[4] = dcomplex (0, 0);
          fi[5] = dcomplex (0, 0);
#else
          fiv[2][ieib] = chi1;
          fiv[3][ieib] = chi0;
          fiv[4][ieib] = dcomplex (0, 0);
          fiv[5][ieib] = dcomplex (0, 0);
#endif
        }
      }
    }
    // ** END LOOP ON IEVT **
    return;
  }


#ifdef __CUDACC__
  __device__
#endif
  void ixzxxxM0(const double pvec[4],
                //const double fmass,
                const int nhel,
                const int nsf,
                dcomplex fi[6])
  {
    fi[0] = dcomplex (-pvec[0] * nsf, -pvec[3] * nsf);
    fi[1] = dcomplex (-pvec[1] * nsf, -pvec[2] * nsf);
    const int nh = nhel * nsf;
    // ASSUMPTIONS FMASS = 0 and
    // (PX and PY are not 0)
    {
      const double sqp0p3 = sqrt( pvec[0] + pvec[3] ) * nsf;
      const dcomplex chi0( sqp0p3, 0.0 );
      const dcomplex chi1( nh * pvec[1] / sqp0p3, pvec[2] / sqp0p3 );
      if (nh == 1)
      {
        fi[2] = dcomplex (0, 0);
        fi[3] = dcomplex (0, 0);
        fi[4] = chi0;
        fi[5] = chi1;
      }
      else
      {
        fi[2] = chi1;
        fi[3] = chi0;
        fi[4] = dcomplex (0, 0);
        fi[5] = dcomplex (0, 0);
      }
    }
    return;
  }


#ifdef __CUDACC__
  __device__
#endif
  void oxzxxxM0(const double pvec[4],
                //const double fmass,
                const int nhel,
                const int nsf,
                dcomplex fo[6])
  {
    fo[0] = dcomplex (pvec[0] * nsf, pvec[3] * nsf);
    fo[1] = dcomplex (pvec[1] * nsf, pvec[2] * nsf);
    const int nh = nhel * nsf;
    // ASSUMPTIONS FMASS = 0 and
    // EITHER (Px and Py are not zero)
    // OR (PX = PY = 0 and E = P3 > 0)
    {
      const double sqp0p3 = sqrt( pvec[0] + pvec[3] ) * nsf;
      const dcomplex chi0( sqp0p3, 0.0 );
      const dcomplex chi1( nh * pvec[1] / sqp0p3, -pvec[2] / sqp0p3 );
      if (nh == 1)
      {
        fo[2] = chi0;
        fo[3] = chi1;
        fo[4] = dcomplex ( 0, 0 );
        fo[5] = dcomplex ( 0, 0 );
      }
      else
      {
        fo[2] = dcomplex ( 0, 0 );
        fo[3] = dcomplex ( 0, 0 );
        fo[4] = chi1;
        fo[5] = chi0;
      }
    }
    return;
  }


#ifdef __CUDACC__
  __device__
#endif
  void FFV1_0(const dcomplex F1[],
              const dcomplex F2[],
              const dcomplex V3[],
              const dcomplex COUP,
              dcomplex * vertex)
  {
    const dcomplex cI = dcomplex (0., 1.);
    const dcomplex TMP4 =
      (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
       (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
        (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
         F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5])))));
    (*vertex) = COUP * - cI * TMP4;
  }


#ifdef __CUDACC__
  __device__
#endif
  void FFV1P0_3(const dcomplex F1[],
                const dcomplex F2[],
                const dcomplex COUP,
                const double M3,
                const double W3,
                dcomplex V3[])
  {
    const dcomplex cI = dcomplex (0., 1.);
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    const double P3[4] = { -V3[0].real(),
                           -V3[1].real(),
                           -V3[1].imag(),
                           -V3[0].imag() };
    const dcomplex denom =
      COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-cI) * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]);
    V3[3] = denom * (-cI) * (-F1[2] * F2[5] - F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]);
    V3[4] = denom * (-cI) * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3] * F2[4] + F1[4] * F2[3]));
    V3[5] = denom * (-cI) * (-F1[2] * F2[4] - F1[5] * F2[3] + F1[3] * F2[5] + F1[4] * F2[2]);
  }


#ifdef __CUDACC__
  __device__
#endif
  void FFV2_4_0(const dcomplex F1[],
                const dcomplex F2[],
                const dcomplex V3[],
                const dcomplex COUP1,
                const dcomplex COUP2,
                dcomplex * vertex)
  {
    const dcomplex cI = dcomplex (0., 1.);
    const dcomplex TMP2 =
      (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
       F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5])));
    const dcomplex TMP0 =
      (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
       F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
    (*vertex) = (-1.) * (COUP2 * (+cI * (TMP0) + 2. * cI * (TMP2)) + cI * (TMP0 * COUP1));
  }


#ifdef __CUDACC__
  __device__
#endif
  void FFV2_4_3(const dcomplex F1[],
                const dcomplex F2[],
                const dcomplex COUP1,
                const dcomplex COUP2,
                const double M3,
                const double W3,
                dcomplex V3[])
  {
    const dcomplex cI = dcomplex (0., 1.);
    double OM3 = 0.;
    if (M3 != 0.) OM3 = 1./(M3 * M3);
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    const double P3[4] = { -V3[0].real(),
                           -V3[1].real(),
                           -V3[1].imag(),
                           -V3[0].imag() };
    const dcomplex TMP1 =
      (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
       F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
    const dcomplex TMP3 =
      (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
       F1[5] * (F2[2] * (-P3[1] + cI * (P3[2])) + F2[3] * (P3[0] + P3[3])));
    const dcomplex denom =
      1./((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) -
          (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-2. * cI) *
      (COUP2 * (OM3 * - 1./2. * P3[0] * (TMP1 + 2. * (TMP3))
                + (+1./2. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] * F2[3]))
       + 1./2. * (COUP1 * (F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP1)));
    V3[3] = denom * (-2. * cI) *
      (COUP2 * (OM3 * - 1./2. * P3[1] * (TMP1 + 2. * (TMP3))
                + (-1./2. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] * F2[2]))
       - 1./2. * (COUP1 * (F1[2] * F2[5] + F1[3] * F2[4] + P3[1] * OM3 * TMP1)));
    V3[4] = denom * cI *
      (COUP2 * (OM3 * P3[2] * (TMP1 + 2. * (TMP3))
                + (+cI * (F1[2] * F2[5]) - cI * (F1[3] * F2[4])
                   - 2. * cI * (F1[4] * F2[3])
                   + 2. * cI * (F1[5] * F2[2])))
       + COUP1 * (+cI * (F1[2] * F2[5]) - cI * (F1[3] * F2[4]) + P3[2] * OM3 * TMP1));
    V3[5] = denom * 2. * cI *
      (COUP2 * (OM3 * 1./2. * P3[3] * (TMP1 + 2. * (TMP3)) +
                (+1./2. * (F1[2] * F2[4]) - 1./2. * (F1[3] * F2[5]) - F1[4] * F2[2] + F1[5] * F2[3]))
       + 1./2. * (COUP1 * (F1[2] * F2[4] + P3[3] * OM3 * TMP1 - F1[3] * F2[5])));
  }


}  // end namespace $(namespace)s_sm


//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <algorithm>
#include <iostream>

#include "mgOnGpuConfig.h"

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
  using mgOnGpu::np4;
  using mgOnGpu::npar;
  const int ncomb = 16; // #helicity combinations is hardcoded for this process (eemumu: ncomb=16)

#ifdef __CUDACC__
  __constant__ int cHel[ncomb][npar];
  __constant__ double cIPC[6];  // coupling ?
  __constant__ double cIPD[2];
#else
  static int cHel[ncomb][npar];
  static double cIPC[6];  // coupling ?
  static double cIPD[2];
#endif

#ifdef __CUDACC__
  __device__ unsigned long long sigmakin_itry = 0; // first iteration over nevt events
  __device__ bool sigmakin_goodhel[ncomb] = { false };
#endif

  // Evaluate |M|^2 for each subprocess
#ifdef __CUDACC__
  __device__
#endif
  // ** NB: allmomenta can have three different layouts
  // ASA: allmomenta[npag][npar][np4][nepp] where ndim=npag*nepp
  // SOA: allmomenta[npar][np4][ndim]
  // AOS: allmomenta[ndim][npar][np4]
  void calculate_wavefunctions( int ihel,
                                const double* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                                double &matrix
#ifndef __CUDACC__
                                , const int ievt
#endif
                                )
  {
#ifdef __CUDACC__
    const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
    const int ievt = idim;
    //printf( "sigmakin: ievt %d\n", ievt );
#endif

    double local_mom[npar][np4];
    for (int ipar = 0; ipar < npar; ipar++ )
      for (int ip4 = 0; ip4 < np4; ip4++ )
        local_mom[ipar][ip4] = MG5_sm::pIparIp4Ievt( allmomenta, ipar, ip4, ievt );

    dcomplex amp[2];
    dcomplex w[5][6];

#ifdef __CUDACC__
    //__shared__ double_v wv0[2*5*6]; // dcomplex_v wv[5][6] gives "dynamic initialization is not supported"
    //dcomplex_v (*wv)[6] = (dcomplex_v (*)[6]) wv0; // dcomplex_v wv[5][6] i.e. dcomplex[5][6][256]
    //__shared__ dcomplex_v wv1[6]; // dcomplex_v wv1[6] gives "dynamic initialization is not supported"
    //__shared__ double_v wv10[2*6]; // dcomplex_v wv1[6] gives "dynamic initialization is not supported"
    double_v wv10[2*6]; // dcomplex_v wv1[6] gives "dynamic initialization is not supported"
    dcomplex_v* wv1 = (dcomplex_v*) wv10; // dcomplex_v wv1[6] i.e. dcomplex[6][256]
#endif

    MG5_sm::oxzxxxM0(local_mom[0], cHel[ihel][0], -1, w[0]);
#ifdef __CUDACC__
    MG5_sm::imzxxxM0( allmomenta, cHel[ihel][1], +1, wv1, 1 );
    const int ieib = threadIdx.x; // event in block
    for (int i6=1; i6<6; i6++) w[1][i6] = wv1[i6][ieib];
#else
    MG5_sm::imzxxxM0( allmomenta, cHel[ihel][1], +1, w[1], ievt, 1 );
#endif
    MG5_sm::ixzxxxM0(local_mom[2], cHel[ihel][2], -1, w[2]);
    MG5_sm::oxzxxxM0(local_mom[3], cHel[ihel][3], +1, w[3]);

    // Diagram 1
    MG5_sm::FFV1P0_3(w[1], w[0], dcomplex (cIPC[0], cIPC[1]), 0., 0., w[4]);
    MG5_sm::FFV1_0(w[2], w[3], w[4], dcomplex (cIPC[0], cIPC[1]), &amp[0]);

    // Diagram 2
    MG5_sm::FFV2_4_3(w[1], w[0], dcomplex (cIPC[2], cIPC[3]), dcomplex (cIPC[4], cIPC[5]), cIPD[0], cIPD[1], w[4]);
    MG5_sm::FFV2_4_0(w[2], w[3], w[4], dcomplex (cIPC[2], cIPC[3]), dcomplex (cIPC[4], cIPC[5]), &amp[1]);

    const int ncolor = 1;
    dcomplex ztemp;
    dcomplex jamp[ncolor];

    // The color matrix;
    static const double denom[ncolor] = {1};
    static const double cf[ncolor][ncolor] = {{1}};

    // Calculate color flows
    jamp[0] = -amp[0] - amp[1];

    // Sum and square the color flows to get the matrix element
    for(int icol = 0; icol < ncolor; icol++ )
    {
      ztemp = 0.;
      for(int jcol = 0; jcol < ncolor; jcol++ )
        ztemp = ztemp + cf[icol][jcol] * jamp[jcol];
      matrix = matrix + (ztemp * conj(jamp[icol])).real()/denom[icol];
    }

    // Store the leading color flows for choice of color
    // for(i=0;i < ncolor; i++)
    // jamp2[0][i] += real(jamp[i]*conj(jamp[i]));

  }

  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess(int numiterations,
                         int gpublocks,
                         int gputhreads,
                         bool verbose,
                         bool debug)
    : m_numiterations(numiterations)
    , gpu_nblocks(gpublocks)
    , gpu_nthreads(gputhreads)
    , dim(gpu_nblocks * gpu_nthreads)
    , m_verbose(verbose)
    , m_debug(debug)
  {
    // Helicities for the process - nodim
    static const int tHel[ncomb][nexternal] =
      { {-1, -1, -1, -1}, {-1, -1, -1, +1}, {-1, -1, +1, -1}, {-1, -1, +1, +1},
        {-1, +1, -1, -1}, {-1, +1, -1, +1}, {-1, +1, +1, -1}, {-1, +1, +1, +1},
        {+1, -1, -1, -1}, {+1, -1, -1, +1}, {+1, -1, +1, -1}, {+1, -1, +1, +1},
        {+1, +1, -1, -1}, {+1, +1, -1, +1}, {+1, +1, +1, -1}, {+1, +1, +1, +1} };
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cHel, tHel, ncomb * nexternal * sizeof(int) ) );
#else
    memcpy( cHel, tHel, ncomb * nexternal * sizeof(int) );
#endif
    // SANITY CHECK: GPU shared memory usage is based on casts of double[2] to complex
    assert( sizeof(dcomplex) == 2*sizeof(double) );
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

  const std::vector<double> &CPPProcess::getMasses() const {return mME;}

  //--------------------------------------------------------------------------
  // Initialize process.

  void CPPProcess::initProc(std::string param_card_name)
  {
    // Instantiate the model class and set parameters that stay fixed during run
    pars = Parameters_sm::getInstance();
    SLHAReader slha(param_card_name, m_verbose);
    pars->setIndependentParameters(slha);
    pars->setIndependentCouplings();
    if (m_verbose) {
      pars->printIndependentParameters();
      pars->printIndependentCouplings();
    }
    pars->setDependentParameters();
    pars->setDependentCouplings();
    // Set external particle masses for this matrix element
    mME.push_back(pars->ZERO);
    mME.push_back(pars->ZERO);
    mME.push_back(pars->ZERO);
    mME.push_back(pars->ZERO);
    static dcomplex tIPC[3] = {pars->GC_3, pars->GC_50,
                               pars->GC_59};
    static double tIPD[2] = {pars->mdl_MZ, pars->mdl_WZ};

#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cIPC, tIPC, 3 * sizeof(dcomplex ) ) );
    checkCuda( cudaMemcpyToSymbol( cIPD, tIPD, 2 * sizeof(double) ) );
#else
    memcpy( cIPC, tIPC, 3 * sizeof(dcomplex ) );
    memcpy( cIPD, tIPD, 2 * sizeof(double) );
#endif

  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour.

  // ** NB: allmomenta can have three different layouts
  // ASA: allmomenta[npag][npar][np4][nepp] where ndim=npag*nepp
  // SOA: allmomenta[npar][np4][ndim]
  // AOS: allmomenta[ndim][npar][np4]
#ifdef __CUDACC__
  __global__
#endif
  void sigmaKin( const double* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                 double* output            // output[nevt]
#ifdef __CUDACC__
                 // NB: nevt == ndim=gpublocks*gputhreads in CUDA
#else
                 , const int nevt          // input: #events
#endif
                 )
  {
    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    // pars->setDependentParameters();
    // pars->setDependentCouplings();
    // Reset color flows

    const int maxtry = 10;
#ifndef __CUDACC__
    static unsigned long long sigmakin_itry = 0; // first iteration over nevt events
    static bool sigmakin_goodhel[ncomb] = { false };
#endif

#ifndef __CUDACC__
    MG5_sm::nevt = nevt;
    // ** START LOOP ON IEVT **
    for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef __CUDACC__
      const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
      const int ievt = idim;
      //printf( "sigmakin: ievt %d\n", ievt );
#endif

      // Denominators: spins, colors and identical particles
      const int nprocesses = 1;
      const int denominators[nprocesses] = {4};

      // Reset the matrix elements
      double matrix_element[nprocesses];
      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] = 0.;
      }

      double melast = matrix_element[0];
      for (int ihel = 0; ihel < ncomb; ihel++ )
      {
        if ( sigmakin_itry>maxtry && !sigmakin_goodhel[ihel] ) continue;
#ifdef __CUDACC__
        calculate_wavefunctions(ihel, allmomenta, matrix_element[0]); // adds ME for ihel to matrix_element[0]
#else
        calculate_wavefunctions(ihel, allmomenta, matrix_element[0], ievt); // adds ME for ihel to matrix_element[0]
#endif
        if ( sigmakin_itry<=maxtry )
        {
          if ( !sigmakin_goodhel[ihel] && matrix_element[0]>melast ) sigmakin_goodhel[ihel] = true;
          melast = matrix_element[0];
        }
      }

      for (int iproc = 0; iproc < nprocesses; ++iproc)
      {
        matrix_element[iproc] /= denominators[iproc];
      }

      for (int iproc = 0; iproc < nprocesses; ++iproc)
      {
        output[iproc*nprocesses + ievt] = matrix_element[iproc];
      }

#ifndef __CUDACC__
      //if ( sigmakin_itry == maxtry )
      //  for (int ihel = 0; ihel < ncomb; ihel++ )
      //    printf( "sigmakin: ihelgood %2d %d\n", ihel, sigmakin_goodhel[ihel] );
      if ( sigmakin_itry <= maxtry )
        sigmakin_itry++;
#else
      if ( sigmakin_itry <= maxtry )
        atomicAdd(&sigmakin_itry, 1);
#endif

    }
    // ** END LOOP ON IEVT **

  }

  //--------------------------------------------------------------------------

}
