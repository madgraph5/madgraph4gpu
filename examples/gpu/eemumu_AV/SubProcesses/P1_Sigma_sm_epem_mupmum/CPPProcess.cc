//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <complex> 
#include <cmath> 
#include <iostream> 
#include <cstdlib> 
#include <thrust/complex.h> 

#include "mgOnGpuConfig.h"

namespace MG5_sm 
{

__device__ void imzxxxM0(const double pvec[4], 
                         //const double fmass, 
                         const int nhel, 
                         const int nsf, 
                         thrust::complex<double> fi[6]) 
{
  fi[0] = thrust::complex<double> (-pvec[0] * nsf, -pvec[3] * nsf); 
  fi[1] = thrust::complex<double> (-pvec[1] * nsf, -pvec[2] * nsf); 
  const int nh = nhel * nsf; 
  // ASSUMPTIONS FMASS = 0 and
  // (PX = PY = 0 and E = -P3 > 0)
  {
    const thrust::complex<double> chi0( 0, 0 );
    const thrust::complex<double> chi1( -nhel * sqrt(2 * pvec[0]), 0 ); 
    if (nh == 1)
    {
      fi[2] = thrust::complex<double> (0, 0); 
      fi[3] = thrust::complex<double> (0, 0); 
      fi[4] = chi0; 
      fi[5] = chi1; 
    }
    else
    {
      fi[2] = chi1; 
      fi[3] = chi0; 
      fi[4] = thrust::complex<double> (0, 0); 
      fi[5] = thrust::complex<double> (0, 0); 
    }
  }
  return; 
}


__device__ void ixzxxxM0(const double pvec[4], 
                         //const double fmass, 
                         const int nhel, 
                         const int nsf, 
                         thrust::complex<double> fi[6]) 
{
  fi[0] = thrust::complex<double> (-pvec[0] * nsf, -pvec[3] * nsf); 
  fi[1] = thrust::complex<double> (-pvec[1] * nsf, -pvec[2] * nsf); 
  const int nh = nhel * nsf; 
  // ASSUMPTIONS FMASS = 0 and
  // (PX and PY are not 0)
  {
    const double sqp0p3 = sqrt( pvec[0] + pvec[3] ) * nsf; 
    const thrust::complex<double> chi0( sqp0p3, 0.0 );
    const thrust::complex<double> chi1( nh * pvec[1] / sqp0p3, pvec[2] / sqp0p3 );
    if (nh == 1)
    {
      fi[2] = thrust::complex<double> (0, 0); 
      fi[3] = thrust::complex<double> (0, 0); 
      fi[4] = chi0; 
      fi[5] = chi1; 
    }
    else
    {
      fi[2] = chi1; 
      fi[3] = chi0; 
      fi[4] = thrust::complex<double> (0, 0); 
      fi[5] = thrust::complex<double> (0, 0); 
    }
  }
  return; 
}



__device__ void oxzxxxM0(const double pvec[4], 
                         //const double fmass, 
                         const int nhel, 
                         const int nsf, 
                         thrust::complex<double> fo[6]) 
{
  fo[0] = thrust::complex<double> (pvec[0] * nsf, pvec[3] * nsf); 
  fo[1] = thrust::complex<double> (pvec[1] * nsf, pvec[2] * nsf); 
  const int nh = nhel * nsf; 
  // ASSUMPTIONS FMASS = 0 and
  // EITHER (Px and Py are not zero)
  // OR (PX = PY = 0 and E = P3 > 0)
  {
    const double sqp0p3 = sqrt( pvec[0] + pvec[3] ) * nsf;
    const thrust::complex<double> chi0( sqp0p3, 0.0 );
    const thrust::complex<double> chi1( nh * pvec[1] / sqp0p3, -pvec[2] / sqp0p3 );
    if (nh == 1)
    {
      fo[2] = chi0; 
      fo[3] = chi1; 
      fo[4] = thrust::complex<double> ( 0, 0 ); 
      fo[5] = thrust::complex<double> ( 0, 0 ); 
    }
    else
    {
      fo[2] = thrust::complex<double> ( 0, 0 ); 
      fo[3] = thrust::complex<double> ( 0, 0 ); 
      fo[4] = chi1; 
      fo[5] = chi0; 
    }
  }
  return; 
}



__device__ void FFV1_0(const thrust::complex<double> F1[], 
                       const thrust::complex<double> F2[], 
                       const thrust::complex<double> V3[], 
                       const thrust::complex<double> COUP, 
                       thrust::complex<double> * vertex)
{
  const thrust::complex<double> cI = thrust::complex<double> (0., 1.); 
  const thrust::complex<double> TMP4 = 
    (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) + 
     (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) + 
      (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) + 
       F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5]))))); 
  (*vertex) = COUP * - cI * TMP4; 
}


__device__ void FFV1P0_3(const thrust::complex<double> F1[], 
                         const thrust::complex<double> F2[], 
                         const thrust::complex<double> COUP, 
                         const double M3, 
                         const double W3, 
                         thrust::complex<double> V3[])
{
  const thrust::complex<double> cI = thrust::complex<double> (0., 1.); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  const double P3[4] = { -V3[0].real(),
                         -V3[1].real(),
                         -V3[1].imag(), 
                         -V3[0].imag() };
  const thrust::complex<double> denom = 
    COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * P3[3]) - M3 * (M3 - cI * W3)); 
  V3[2] = denom * (-cI) * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]); 
  V3[3] = denom * (-cI) * (-F1[2] * F2[5] - F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]); 
  V3[4] = denom * (-cI) * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3] * F2[4] + F1[4] * F2[3])); 
  V3[5] = denom * (-cI) * (-F1[2] * F2[4] - F1[5] * F2[3] + F1[3] * F2[5] + F1[4] * F2[2]); 
}


__device__ void FFV2_4_0(const thrust::complex<double> F1[], 
                         const thrust::complex<double> F2[], 
                         const thrust::complex<double> V3[], 
                         const thrust::complex<double> COUP1, 
                         const thrust::complex<double> COUP2, 
                         thrust::complex<double> * vertex)
{
  const thrust::complex<double> cI = thrust::complex<double> (0., 1.); 
  const thrust::complex<double> TMP2 = 
    (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) + 
     F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5]))); 
  const thrust::complex<double> TMP0 = 
    (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) + 
     F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5]))); 
  (*vertex) = (-1.) * (COUP2 * (+cI * (TMP0) + 2. * cI * (TMP2)) + cI * (TMP0 * COUP1)); 
}


__device__ void FFV2_4_3(const thrust::complex<double> F1[], 
                         const thrust::complex<double> F2[], 
                         const thrust::complex<double> COUP1, 
                         const thrust::complex<double> COUP2, 
                         const double M3, 
                         const double W3, 
                         thrust::complex<double> V3[])
{
  const thrust::complex<double> cI = thrust::complex<double> (0., 1.); 
  double OM3 = 0.; 
  if (M3 != 0.) OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1];
  const double P3[4] = { -V3[0].real(), 
                         -V3[1].real(), 
                         -V3[1].imag(),
                         -V3[0].imag() };  
  const thrust::complex<double> TMP1 = 
    (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) + 
     F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  const thrust::complex<double> TMP3 = 
    (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) + 
     F1[5] * (F2[2] * (-P3[1] + cI * (P3[2])) + F2[3] * (P3[0] + P3[3]))); 
  const thrust::complex<double> denom = 
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
#include <thrust/complex.h> 

#include "CPPProcess.h"

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1

__constant__ int cHel[16][4]; 

// __constant__ double cmME[4]; value hardcoded now
// extern __constant__ int cPerm[4];

__constant__ double cIPC[6];  // coupling ?
__constant__ double cIPD[2]; 


// Evaluate |M|^2 for each subprocess
__device__ void calculate_wavefunctions(int ihel, double local_mom[4][4],
    double &matrix)
{
  thrust::complex<double> amp[2]; 
  // Calculate wavefunctions for all processes
  thrust::complex<double> w[5][6]; 

  //MG5_sm::oxxxxx(local_mom[0], 0., cHel[ihel][0], -1, w[0]); 
  //MG5_sm::oxxxxxM0(local_mom[0], cHel[ihel][0], -1, w[0]); 
  MG5_sm::oxzxxxM0(local_mom[0], cHel[ihel][0], -1, w[0]); 

  //MG5_sm::ixxxxx(local_mom[1], 0., cHel[ihel][1], +1, w[1]); 
  //MG5_sm::ixxxxxM0(local_mom[1], cHel[ihel][1], +1, w[1]); 
  MG5_sm::imzxxxM0(local_mom[1], cHel[ihel][1], +1, w[1]); 

  //MG5_sm::ixxxxx(local_mom[2], 0., cHel[ihel][2], -1, w[2]);
  //MG5_sm::ixxxxxM0(local_mom[2], cHel[ihel][2], -1, w[2]);
  MG5_sm::ixzxxxM0(local_mom[2], cHel[ihel][2], -1, w[2]);

  //MG5_sm::oxxxxx(local_mom[3], 0., cHel[ihel][3], +1, w[3]); 
  //MG5_sm::oxxxxxM0(local_mom[3], cHel[ihel][3], +1, w[3]); 
  MG5_sm::oxzxxxM0(local_mom[3], cHel[ihel][3], +1, w[3]); 

  MG5_sm::FFV1P0_3(w[1], w[0], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0., w[4]);
  // Amplitude(s) for diagram number 1
  MG5_sm::FFV1_0(w[2], w[3], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), &amp[0]);
  MG5_sm::FFV2_4_3(w[1], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), 
                   thrust::complex<double> (cIPC[4], cIPC[5]), cIPD[0], cIPD[1], w[4]);
  // Amplitude(s) for diagram number 2
  MG5_sm::FFV2_4_0(w[2], w[3], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
                   thrust::complex<double> (cIPC[4], cIPC[5]), &amp[1]);

  const int ncolor = 1; 
  thrust::complex<double> ztemp; 
  thrust::complex<double> jamp[ncolor]; 

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



CPPProcess::CPPProcess(int numiterations, int gpublocks, int gputhreads, 
bool verbose, bool debug)
: m_numiterations(numiterations), gpu_nblocks(gpublocks), 
gpu_nthreads(gputhreads), m_verbose(verbose), m_debug(debug), 
dim(gpu_nblocks * gpu_nthreads) 
{


  // Helicities for the process - nodim
  static const int tHel[ncomb][nexternal] = {{-1, -1, -1, -1}, {-1, -1, -1, 1},
      {-1, -1, 1, -1}, {-1, -1, 1, 1}, {-1, 1, -1, -1}, {-1, 1, -1, 1}, {-1, 1,
      1, -1}, {-1, 1, 1, 1}, {1, -1, -1, -1}, {1, -1, -1, 1}, {1, -1, 1, -1},
      {1, -1, 1, 1}, {1, 1, -1, -1}, {1, 1, -1, 1}, {1, 1, 1, -1}, {1, 1, 1,
      1}};
  gpuErrchk3( cudaMemcpyToSymbol( cHel, tHel, ncomb * nexternal * sizeof(int) ) ); 
  // perm - nodim
  // static int perm[nexternal] = {0, 1, 2, 3};
}

CPPProcess::~CPPProcess() {}

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
  static thrust::complex<double> tIPC[3] = {pars->GC_3, pars->GC_50,
      pars->GC_59};
  static double tIPD[2] = {pars->mdl_MZ, pars->mdl_WZ}; 

  gpuErrchk3( cudaMemcpyToSymbol( cIPC, tIPC, 3 * sizeof(thrust::complex<double> ) ) );
  gpuErrchk3( cudaMemcpyToSymbol( cIPD, tIPD, 2 * sizeof(double) ) );
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

__global__ 
#if defined MGONGPU_LAYOUT_ASA
// AOSOA: allmomenta[npag][npar][np4][nepp] where ndim=npag*nepp
#elif defined MGONGPU_LAYOUT_SOA
// SOA: allmomenta[npar][np4][ndim]
#elif defined MGONGPU_LAYOUT_AOS
// AOS: allmomenta[ndim][npar][np4]
#endif
void sigmaKin( const double* allmomenta, // input[(npar=4)*(np4=4)*(ndim=gpublocks*gputhreads)]
               double* output ) // output[nevt]
{
  // Set the parameters which change event by event
  // Need to discuss this with Stefan
  // pars->setDependentParameters();
  // pars->setDependentCouplings();
  // Reset color flows

  const int npar = 4; // hardcoded for this process (eemumu): npar=4
  const int np4 = 4; // dimension of 4-momenta (E,px,py,pz): copy all of them from rambo

  const int idim = blockIdx.x * blockDim.x + threadIdx.x; // event# == threadid (previously was: tid)

  double local_m[npar][np4];
#if defined MGONGPU_LAYOUT_ASA
  using mgOnGpu::nepp;
  const int ipag = idim/nepp; // #eventpage in this iteration
  const int iepp = idim%nepp; // #event in the current eventpage in this iteration
  // ASA: allmomenta[npag][npar][np4][nepp]
  for (int ipar = 0; ipar < npar; ipar++ )
    for (int ip4 = 0; ip4 < np4; ip4++ )
      local_m[ipar][ip4] = allmomenta[ipag*npar*np4*nepp + ipar*nepp*np4 + ip4*nepp + iepp]; // AOSOA[ipag][ipar][ip4][iepp]
#elif defined MGONGPU_LAYOUT_SOA
  const int ndim = blockDim.x * gridDim.x; // (previously was: DIM)
  // SOA: allmomenta[npar][np4][ndim]
  for (int ipar = 0; ipar < npar; ipar++ )
    for (int ip4 = 0; ip4 < np4; ip4++ )
      local_m[ipar][ip4] = allmomenta[ipar*np4*ndim + ip4*ndim + idim]; // SOA[ipar][ip4][idim]
#elif defined MGONGPU_LAYOUT_AOS
  // AOS: allmomenta[ndim][npar][np4]
  for (int ipar = 0; ipar < npar; ipar++ )
    for (int ip4 = 0; ip4 < np4; ip4++ )
      local_m[ipar][ip4] = allmomenta[idim*npar*np4 + ipar*np4 + ip4]; // AOS[idim][ipar][ip4]
#endif

  // Helicity combinations
  const int ncomb = 16;

  // Denominators: spins, colors and identical particles
  const int nprocesses = 1;
  const int denominators[nprocesses] = {4};

  // Reset the matrix elements
  double matrix_element[nprocesses];
  for(int iproc = 0; iproc < nprocesses; iproc++ )
  {
    matrix_element[iproc] = 0.; 
  }

  for (int ihel = 0; ihel < ncomb; ihel++ )
  {
    calculate_wavefunctions(ihel, local_m, matrix_element[0]); 
  }

  for (int iproc = 0; iproc < nprocesses; ++iproc)
  {
    matrix_element[iproc] /= denominators[iproc]; 
  }

  for (int iproc = 0; iproc < nprocesses; ++iproc)
  {
    output[iproc*nprocesses + idim] = matrix_element[iproc]; 
  }

}

//--------------------------------------------------------------------------
