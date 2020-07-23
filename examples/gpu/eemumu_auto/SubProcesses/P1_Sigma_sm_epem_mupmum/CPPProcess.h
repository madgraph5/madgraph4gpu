//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_epem_mupmum_H
#define MG5_Sigma_sm_epem_mupmum_H

#include <complex> 
#include <vector> 


#include "Parameters_sm.h"

#include <thrust/complex.h> 


__global__ void sigmaKin(double * allmomenta); 

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

    thrust::complex<double> * * amp; 


    // Pointer to the model parameters
    Parameters_sm * pars; 

    // vector with external particle masses
    std::vector<double> mME; 

    // Initial particle ids
    int id1, id2; 

}; 


#endif  // MG5_Sigma_sm_epem_mupmum_H
