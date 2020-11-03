//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

#include <cmath> 
#include <thrust/complex.h> 


using namespace std; 

namespace MG5_sm 
{
__device__ void oxxxxx(double p[4], double fmass, int nhel, int nsf, 
thrust::complex<double> fo[6]); 

__device__ void sxxxxx(double p[4], int nss, thrust::complex<double> sc[3]); 

__device__ void ixxxxx(double p[4], double fmass, int nhel, int nsf, 
thrust::complex<double> fi[6]); 

__device__ void txxxxx(double p[4], double tmass, int nhel, int nst, 
thrust::complex<double> fi[18]); 

__device__ void vxxxxx(double p[4], double vmass, int nhel, int nsv, 
thrust::complex<double> v[6]); 

__device__ void VVV1P0_1(thrust::complex<double> V2[], const
thrust::complex<double> V3[], const thrust::complex<double> COUP, const
double M1, const double W1, thrust::complex<double> V1[]); 

__device__ void FFV1_0(thrust::complex<double> F1[], const
thrust::complex<double> F2[], const thrust::complex<double> V3[], const
thrust::complex<double> COUP, thrust::complex<double> * vertex); 

__device__ void FFV1_1(thrust::complex<double> F2[], const
thrust::complex<double> V3[], const thrust::complex<double> COUP, const
double M1, const double W1, thrust::complex<double> F1[]); 

__device__ void FFV1_2(thrust::complex<double> F1[], const
thrust::complex<double> V3[], const thrust::complex<double> COUP, const
double M2, const double W2, thrust::complex<double> F2[]); 

}  // end namespace MG5_sm

#endif  // HelAmps_sm_H

//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_gg_ttx_H
#define MG5_Sigma_sm_gg_ttx_H

#include <complex> 
#include <vector> 


#include "Parameters_sm.h"

#include <thrust/complex.h> 


__global__ void sigmaKin(double * allmomenta, double * output); 

//==========================================================================
// A class for calculating the matrix elements for
// Process: g g > t t~ WEIGHTED<=2 @1
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
    static const int namplitudes = 3; 
    static const int ncomb = 16; 
    static const int wrows = 7; 
    // static const int nioparticles = 4;

    thrust::complex<double> * * amp; 


    // Pointer to the model parameters
    Parameters_sm * pars; 

    // vector with external particle masses
    std::vector<double> mME; 

    // Initial particle ids
    int id1, id2; 

}; 


#endif  // MG5_Sigma_sm_gg_ttx_H
