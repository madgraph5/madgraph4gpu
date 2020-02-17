//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.0, 2020-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"
#include "HelAmps_sm_ckm_no_b_mass.h"

using namespace MG5_sm_ckm_no_b_mass; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: d s > e+ e- g g d s WEIGHTED<=8 / h
// Process: d s > mu+ mu- g g d s WEIGHTED<=8 / h
// Process: d b > e+ e- g g d b WEIGHTED<=8 / h
// Process: d b > mu+ mu- g g d b WEIGHTED<=8 / h
// Process: s b > e+ e- g g s b WEIGHTED<=8 / h
// Process: s b > mu+ mu- g g s b WEIGHTED<=8 / h

//--------------------------------------------------------------------------
// Initialize process.

void CPPProcess::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm_ckm_no_b_mass::getInstance(); 
  SLHAReader slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  pars->printIndependentParameters(); 
  pars->printIndependentCouplings(); 
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[12]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void CPPProcess::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 
  static bool firsttime = true; 
  if (firsttime)
  {
    pars->printDependentParameters(); 
    pars->printDependentCouplings(); 
    firsttime = false; 
  }

  // Reset color flows
  for(int i = 0; i < 12; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 256; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
  std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1, -1,
      -1}, {-1, -1, -1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, -1, -1, 1, -1},
      {-1, -1, -1, -1, -1, -1, 1, 1}, {-1, -1, -1, -1, -1, 1, -1, -1}, {-1, -1,
      -1, -1, -1, 1, -1, 1}, {-1, -1, -1, -1, -1, 1, 1, -1}, {-1, -1, -1, -1,
      -1, 1, 1, 1}, {-1, -1, -1, -1, 1, -1, -1, -1}, {-1, -1, -1, -1, 1, -1,
      -1, 1}, {-1, -1, -1, -1, 1, -1, 1, -1}, {-1, -1, -1, -1, 1, -1, 1, 1},
      {-1, -1, -1, -1, 1, 1, -1, -1}, {-1, -1, -1, -1, 1, 1, -1, 1}, {-1, -1,
      -1, -1, 1, 1, 1, -1}, {-1, -1, -1, -1, 1, 1, 1, 1}, {-1, -1, -1, 1, -1,
      -1, -1, -1}, {-1, -1, -1, 1, -1, -1, -1, 1}, {-1, -1, -1, 1, -1, -1, 1,
      -1}, {-1, -1, -1, 1, -1, -1, 1, 1}, {-1, -1, -1, 1, -1, 1, -1, -1}, {-1,
      -1, -1, 1, -1, 1, -1, 1}, {-1, -1, -1, 1, -1, 1, 1, -1}, {-1, -1, -1, 1,
      -1, 1, 1, 1}, {-1, -1, -1, 1, 1, -1, -1, -1}, {-1, -1, -1, 1, 1, -1, -1,
      1}, {-1, -1, -1, 1, 1, -1, 1, -1}, {-1, -1, -1, 1, 1, -1, 1, 1}, {-1, -1,
      -1, 1, 1, 1, -1, -1}, {-1, -1, -1, 1, 1, 1, -1, 1}, {-1, -1, -1, 1, 1, 1,
      1, -1}, {-1, -1, -1, 1, 1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1, -1, -1},
      {-1, -1, 1, -1, -1, -1, -1, 1}, {-1, -1, 1, -1, -1, -1, 1, -1}, {-1, -1,
      1, -1, -1, -1, 1, 1}, {-1, -1, 1, -1, -1, 1, -1, -1}, {-1, -1, 1, -1, -1,
      1, -1, 1}, {-1, -1, 1, -1, -1, 1, 1, -1}, {-1, -1, 1, -1, -1, 1, 1, 1},
      {-1, -1, 1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, 1, -1, -1, 1}, {-1, -1,
      1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, -1, 1, 1}, {-1, -1, 1, -1, 1, 1,
      -1, -1}, {-1, -1, 1, -1, 1, 1, -1, 1}, {-1, -1, 1, -1, 1, 1, 1, -1}, {-1,
      -1, 1, -1, 1, 1, 1, 1}, {-1, -1, 1, 1, -1, -1, -1, -1}, {-1, -1, 1, 1,
      -1, -1, -1, 1}, {-1, -1, 1, 1, -1, -1, 1, -1}, {-1, -1, 1, 1, -1, -1, 1,
      1}, {-1, -1, 1, 1, -1, 1, -1, -1}, {-1, -1, 1, 1, -1, 1, -1, 1}, {-1, -1,
      1, 1, -1, 1, 1, -1}, {-1, -1, 1, 1, -1, 1, 1, 1}, {-1, -1, 1, 1, 1, -1,
      -1, -1}, {-1, -1, 1, 1, 1, -1, -1, 1}, {-1, -1, 1, 1, 1, -1, 1, -1}, {-1,
      -1, 1, 1, 1, -1, 1, 1}, {-1, -1, 1, 1, 1, 1, -1, -1}, {-1, -1, 1, 1, 1,
      1, -1, 1}, {-1, -1, 1, 1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, -1, -1, 1}, {-1, 1, -1,
      -1, -1, -1, 1, -1}, {-1, 1, -1, -1, -1, -1, 1, 1}, {-1, 1, -1, -1, -1, 1,
      -1, -1}, {-1, 1, -1, -1, -1, 1, -1, 1}, {-1, 1, -1, -1, -1, 1, 1, -1},
      {-1, 1, -1, -1, -1, 1, 1, 1}, {-1, 1, -1, -1, 1, -1, -1, -1}, {-1, 1, -1,
      -1, 1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1, 1, -1}, {-1, 1, -1, -1, 1, -1,
      1, 1}, {-1, 1, -1, -1, 1, 1, -1, -1}, {-1, 1, -1, -1, 1, 1, -1, 1}, {-1,
      1, -1, -1, 1, 1, 1, -1}, {-1, 1, -1, -1, 1, 1, 1, 1}, {-1, 1, -1, 1, -1,
      -1, -1, -1}, {-1, 1, -1, 1, -1, -1, -1, 1}, {-1, 1, -1, 1, -1, -1, 1,
      -1}, {-1, 1, -1, 1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, 1, -1, -1}, {-1, 1,
      -1, 1, -1, 1, -1, 1}, {-1, 1, -1, 1, -1, 1, 1, -1}, {-1, 1, -1, 1, -1, 1,
      1, 1}, {-1, 1, -1, 1, 1, -1, -1, -1}, {-1, 1, -1, 1, 1, -1, -1, 1}, {-1,
      1, -1, 1, 1, -1, 1, -1}, {-1, 1, -1, 1, 1, -1, 1, 1}, {-1, 1, -1, 1, 1,
      1, -1, -1}, {-1, 1, -1, 1, 1, 1, -1, 1}, {-1, 1, -1, 1, 1, 1, 1, -1},
      {-1, 1, -1, 1, 1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, -1, -1, 1}, {-1, 1, 1, -1, -1, -1, 1, -1}, {-1, 1, 1, -1, -1, -1,
      1, 1}, {-1, 1, 1, -1, -1, 1, -1, -1}, {-1, 1, 1, -1, -1, 1, -1, 1}, {-1,
      1, 1, -1, -1, 1, 1, -1}, {-1, 1, 1, -1, -1, 1, 1, 1}, {-1, 1, 1, -1, 1,
      -1, -1, -1}, {-1, 1, 1, -1, 1, -1, -1, 1}, {-1, 1, 1, -1, 1, -1, 1, -1},
      {-1, 1, 1, -1, 1, -1, 1, 1}, {-1, 1, 1, -1, 1, 1, -1, -1}, {-1, 1, 1, -1,
      1, 1, -1, 1}, {-1, 1, 1, -1, 1, 1, 1, -1}, {-1, 1, 1, -1, 1, 1, 1, 1},
      {-1, 1, 1, 1, -1, -1, -1, -1}, {-1, 1, 1, 1, -1, -1, -1, 1}, {-1, 1, 1,
      1, -1, -1, 1, -1}, {-1, 1, 1, 1, -1, -1, 1, 1}, {-1, 1, 1, 1, -1, 1, -1,
      -1}, {-1, 1, 1, 1, -1, 1, -1, 1}, {-1, 1, 1, 1, -1, 1, 1, -1}, {-1, 1, 1,
      1, -1, 1, 1, 1}, {-1, 1, 1, 1, 1, -1, -1, -1}, {-1, 1, 1, 1, 1, -1, -1,
      1}, {-1, 1, 1, 1, 1, -1, 1, -1}, {-1, 1, 1, 1, 1, -1, 1, 1}, {-1, 1, 1,
      1, 1, 1, -1, -1}, {-1, 1, 1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, 1, 1, -1},
      {-1, 1, 1, 1, 1, 1, 1, 1}, {1, -1, -1, -1, -1, -1, -1, -1}, {1, -1, -1,
      -1, -1, -1, -1, 1}, {1, -1, -1, -1, -1, -1, 1, -1}, {1, -1, -1, -1, -1,
      -1, 1, 1}, {1, -1, -1, -1, -1, 1, -1, -1}, {1, -1, -1, -1, -1, 1, -1, 1},
      {1, -1, -1, -1, -1, 1, 1, -1}, {1, -1, -1, -1, -1, 1, 1, 1}, {1, -1, -1,
      -1, 1, -1, -1, -1}, {1, -1, -1, -1, 1, -1, -1, 1}, {1, -1, -1, -1, 1, -1,
      1, -1}, {1, -1, -1, -1, 1, -1, 1, 1}, {1, -1, -1, -1, 1, 1, -1, -1}, {1,
      -1, -1, -1, 1, 1, -1, 1}, {1, -1, -1, -1, 1, 1, 1, -1}, {1, -1, -1, -1,
      1, 1, 1, 1}, {1, -1, -1, 1, -1, -1, -1, -1}, {1, -1, -1, 1, -1, -1, -1,
      1}, {1, -1, -1, 1, -1, -1, 1, -1}, {1, -1, -1, 1, -1, -1, 1, 1}, {1, -1,
      -1, 1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1, -1, 1}, {1, -1, -1, 1, -1,
      1, 1, -1}, {1, -1, -1, 1, -1, 1, 1, 1}, {1, -1, -1, 1, 1, -1, -1, -1},
      {1, -1, -1, 1, 1, -1, -1, 1}, {1, -1, -1, 1, 1, -1, 1, -1}, {1, -1, -1,
      1, 1, -1, 1, 1}, {1, -1, -1, 1, 1, 1, -1, -1}, {1, -1, -1, 1, 1, 1, -1,
      1}, {1, -1, -1, 1, 1, 1, 1, -1}, {1, -1, -1, 1, 1, 1, 1, 1}, {1, -1, 1,
      -1, -1, -1, -1, -1}, {1, -1, 1, -1, -1, -1, -1, 1}, {1, -1, 1, -1, -1,
      -1, 1, -1}, {1, -1, 1, -1, -1, -1, 1, 1}, {1, -1, 1, -1, -1, 1, -1, -1},
      {1, -1, 1, -1, -1, 1, -1, 1}, {1, -1, 1, -1, -1, 1, 1, -1}, {1, -1, 1,
      -1, -1, 1, 1, 1}, {1, -1, 1, -1, 1, -1, -1, -1}, {1, -1, 1, -1, 1, -1,
      -1, 1}, {1, -1, 1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, -1, 1, 1}, {1,
      -1, 1, -1, 1, 1, -1, -1}, {1, -1, 1, -1, 1, 1, -1, 1}, {1, -1, 1, -1, 1,
      1, 1, -1}, {1, -1, 1, -1, 1, 1, 1, 1}, {1, -1, 1, 1, -1, -1, -1, -1}, {1,
      -1, 1, 1, -1, -1, -1, 1}, {1, -1, 1, 1, -1, -1, 1, -1}, {1, -1, 1, 1, -1,
      -1, 1, 1}, {1, -1, 1, 1, -1, 1, -1, -1}, {1, -1, 1, 1, -1, 1, -1, 1}, {1,
      -1, 1, 1, -1, 1, 1, -1}, {1, -1, 1, 1, -1, 1, 1, 1}, {1, -1, 1, 1, 1, -1,
      -1, -1}, {1, -1, 1, 1, 1, -1, -1, 1}, {1, -1, 1, 1, 1, -1, 1, -1}, {1,
      -1, 1, 1, 1, -1, 1, 1}, {1, -1, 1, 1, 1, 1, -1, -1}, {1, -1, 1, 1, 1, 1,
      -1, 1}, {1, -1, 1, 1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1, 1, 1}, {1, 1, -1,
      -1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, -1, -1, 1}, {1, 1, -1, -1, -1,
      -1, 1, -1}, {1, 1, -1, -1, -1, -1, 1, 1}, {1, 1, -1, -1, -1, 1, -1, -1},
      {1, 1, -1, -1, -1, 1, -1, 1}, {1, 1, -1, -1, -1, 1, 1, -1}, {1, 1, -1,
      -1, -1, 1, 1, 1}, {1, 1, -1, -1, 1, -1, -1, -1}, {1, 1, -1, -1, 1, -1,
      -1, 1}, {1, 1, -1, -1, 1, -1, 1, -1}, {1, 1, -1, -1, 1, -1, 1, 1}, {1, 1,
      -1, -1, 1, 1, -1, -1}, {1, 1, -1, -1, 1, 1, -1, 1}, {1, 1, -1, -1, 1, 1,
      1, -1}, {1, 1, -1, -1, 1, 1, 1, 1}, {1, 1, -1, 1, -1, -1, -1, -1}, {1, 1,
      -1, 1, -1, -1, -1, 1}, {1, 1, -1, 1, -1, -1, 1, -1}, {1, 1, -1, 1, -1,
      -1, 1, 1}, {1, 1, -1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1, -1, 1}, {1,
      1, -1, 1, -1, 1, 1, -1}, {1, 1, -1, 1, -1, 1, 1, 1}, {1, 1, -1, 1, 1, -1,
      -1, -1}, {1, 1, -1, 1, 1, -1, -1, 1}, {1, 1, -1, 1, 1, -1, 1, -1}, {1, 1,
      -1, 1, 1, -1, 1, 1}, {1, 1, -1, 1, 1, 1, -1, -1}, {1, 1, -1, 1, 1, 1, -1,
      1}, {1, 1, -1, 1, 1, 1, 1, -1}, {1, 1, -1, 1, 1, 1, 1, 1}, {1, 1, 1, -1,
      -1, -1, -1, -1}, {1, 1, 1, -1, -1, -1, -1, 1}, {1, 1, 1, -1, -1, -1, 1,
      -1}, {1, 1, 1, -1, -1, -1, 1, 1}, {1, 1, 1, -1, -1, 1, -1, -1}, {1, 1, 1,
      -1, -1, 1, -1, 1}, {1, 1, 1, -1, -1, 1, 1, -1}, {1, 1, 1, -1, -1, 1, 1,
      1}, {1, 1, 1, -1, 1, -1, -1, -1}, {1, 1, 1, -1, 1, -1, -1, 1}, {1, 1, 1,
      -1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, -1, 1, 1}, {1, 1, 1, -1, 1, 1, -1,
      -1}, {1, 1, 1, -1, 1, 1, -1, 1}, {1, 1, 1, -1, 1, 1, 1, -1}, {1, 1, 1,
      -1, 1, 1, 1, 1}, {1, 1, 1, 1, -1, -1, -1, -1}, {1, 1, 1, 1, -1, -1, -1,
      1}, {1, 1, 1, 1, -1, -1, 1, -1}, {1, 1, 1, 1, -1, -1, 1, 1}, {1, 1, 1, 1,
      -1, 1, -1, -1}, {1, 1, 1, 1, -1, 1, -1, 1}, {1, 1, 1, 1, -1, 1, 1, -1},
      {1, 1, 1, 1, -1, 1, 1, 1}, {1, 1, 1, 1, 1, -1, -1, -1}, {1, 1, 1, 1, 1,
      -1, -1, 1}, {1, 1, 1, 1, 1, -1, 1, -1}, {1, 1, 1, 1, 1, -1, 1, 1}, {1, 1,
      1, 1, 1, 1, -1, -1}, {1, 1, 1, 1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, 1, 1,
      -1}, {1, 1, 1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {72, 72}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_ds_epemggds_no_h(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_ds_epemggds_no_h(); 
        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc]; 
          tsum += t[iproc]; 
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true; 
          ngood++; 
          igood[ngood] = ihel; 
        }
      }
    }
    jhel = 0; 
    sum_hel = min(sum_hel, ngood); 
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_ds_epemggds_no_h(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_ds_epemggds_no_h(); 
      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double CPPProcess::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 1 && id2 == 3)
  {
    // Add matrix elements for processes with beams (1, 3)
    return matrix_element[0] * 2; 
  }
  else if(id1 == 3 && id2 == 1)
  {
    // Add matrix elements for processes with beams (3, 1)
    return matrix_element[1] * 2; 
  }
  else if(id1 == 1 && id2 == 5)
  {
    // Add matrix elements for processes with beams (1, 5)
    return matrix_element[0] * 2; 
  }
  else if(id1 == 5 && id2 == 1)
  {
    // Add matrix elements for processes with beams (5, 1)
    return matrix_element[1] * 2; 
  }
  else if(id1 == 3 && id2 == 5)
  {
    // Add matrix elements for processes with beams (3, 5)
    return matrix_element[0] * 2; 
  }
  else if(id1 == 5 && id2 == 3)
  {
    // Add matrix elements for processes with beams (5, 3)
    return matrix_element[1] * 2; 
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.; 
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void CPPProcess::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j; 

  // Calculate all wavefunctions
  ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
  ixxxxx(p[perm[1]], mME[1], hel[1], +1, w[1]); 
  ixxxxx(p[perm[2]], mME[2], hel[2], -1, w[2]); 
  oxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]); 
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  vxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]); 
  oxxxxx(p[perm[6]], mME[6], hel[6], +1, w[6]); 
  oxxxxx(p[perm[7]], mME[7], hel[7], +1, w[7]); 
  VVV1P0_1(w[4], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[8]); 
  FFV1P0_3(w[2], w[3], pars->GC_3, pars->ZERO, pars->ZERO, w[9]); 
  FFV1_1(w[6], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[10]); 
  FFV1_1(w[7], w[9], pars->GC_1, pars->ZERO, pars->ZERO, w[11]); 
  FFV1P0_3(w[0], w[10], pars->GC_11, pars->ZERO, pars->ZERO, w[12]); 
  FFV1_2(w[1], w[9], pars->GC_1, pars->ZERO, pars->ZERO, w[13]); 
  FFV1_1(w[7], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[14]); 
  FFV1_1(w[6], w[9], pars->GC_1, pars->ZERO, pars->ZERO, w[15]); 
  FFV1P0_3(w[1], w[14], pars->GC_11, pars->ZERO, pars->ZERO, w[16]); 
  FFV1_2(w[0], w[9], pars->GC_1, pars->ZERO, pars->ZERO, w[17]); 
  FFV1_2(w[0], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[18]); 
  FFV1P0_3(w[18], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[19]); 
  FFV1_2(w[1], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[20]); 
  FFV1P0_3(w[20], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[21]); 
  FFV2_4_3(w[2], w[3], -pars->GC_51, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ,
      w[22]);
  FFV2_3_1(w[7], w[22], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[23]);
  FFV2_3_2(w[1], w[22], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[24]);
  FFV2_3_1(w[6], w[22], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[25]);
  FFV2_3_2(w[0], w[22], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[26]);
  FFV1P0_3(w[0], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[27]); 
  VVV1P0_1(w[8], w[27], pars->GC_10, pars->ZERO, pars->ZERO, w[28]); 
  FFV1_2(w[1], w[27], pars->GC_11, pars->ZERO, pars->ZERO, w[29]); 
  FFV1_1(w[7], w[27], pars->GC_11, pars->ZERO, pars->ZERO, w[30]); 
  FFV1P0_3(w[1], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[31]); 
  FFV1_2(w[0], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[32]); 
  VVV1P0_1(w[8], w[31], pars->GC_10, pars->ZERO, pars->ZERO, w[33]); 
  FFV1_1(w[6], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[34]); 
  FFV1_1(w[6], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[35]); 
  FFV1_1(w[7], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[36]); 
  FFV1_1(w[35], w[9], pars->GC_1, pars->ZERO, pars->ZERO, w[37]); 
  FFV1P0_3(w[1], w[36], pars->GC_11, pars->ZERO, pars->ZERO, w[38]); 
  FFV1P0_3(w[0], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[39]); 
  FFV1_1(w[36], w[9], pars->GC_1, pars->ZERO, pars->ZERO, w[40]); 
  FFV2_3_1(w[35], w[22], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[41]);
  FFV2_3_1(w[36], w[22], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[42]);
  FFV1_2(w[0], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[43]); 
  FFV1P0_3(w[43], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[44]); 
  FFV1_1(w[35], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[45]); 
  FFV1_2(w[1], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[46]); 
  FFV1P0_3(w[46], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[47]); 
  FFV1_2(w[46], w[9], pars->GC_1, pars->ZERO, pars->ZERO, w[48]); 
  FFV2_3_2(w[46], w[22], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[49]);
  FFV1_1(w[35], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[50]); 
  FFV1P0_3(w[0], w[50], pars->GC_11, pars->ZERO, pars->ZERO, w[51]); 
  VVV1P0_1(w[39], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[52]); 
  FFV1_2(w[1], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[53]); 
  FFV1_1(w[7], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[54]); 
  VVV1P0_1(w[5], w[31], pars->GC_10, pars->ZERO, pars->ZERO, w[55]); 
  FFV1_1(w[7], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[56]); 
  FFV1_1(w[6], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[57]); 
  FFV1_1(w[56], w[9], pars->GC_1, pars->ZERO, pars->ZERO, w[58]); 
  FFV1P0_3(w[0], w[57], pars->GC_11, pars->ZERO, pars->ZERO, w[59]); 
  FFV1P0_3(w[1], w[56], pars->GC_11, pars->ZERO, pars->ZERO, w[60]); 
  FFV1_1(w[57], w[9], pars->GC_1, pars->ZERO, pars->ZERO, w[61]); 
  FFV2_3_1(w[56], w[22], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[62]);
  FFV2_3_1(w[57], w[22], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[63]);
  FFV1P0_3(w[43], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[64]); 
  FFV1_2(w[43], w[9], pars->GC_1, pars->ZERO, pars->ZERO, w[65]); 
  FFV2_3_2(w[43], w[22], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[66]);
  FFV1P0_3(w[46], w[56], pars->GC_11, pars->ZERO, pars->ZERO, w[67]); 
  FFV1_1(w[56], w[27], pars->GC_11, pars->ZERO, pars->ZERO, w[68]); 
  FFV1_1(w[56], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[69]); 
  FFV1P0_3(w[1], w[69], pars->GC_11, pars->ZERO, pars->ZERO, w[70]); 
  VVV1P0_1(w[60], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[71]); 
  FFV1_2(w[0], w[60], pars->GC_11, pars->ZERO, pars->ZERO, w[72]); 
  FFV1_1(w[6], w[60], pars->GC_11, pars->ZERO, pars->ZERO, w[73]); 
  VVV1P0_1(w[5], w[27], pars->GC_10, pars->ZERO, pars->ZERO, w[74]); 
  FFV1_2(w[0], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[75]); 
  FFV1P0_3(w[75], w[57], pars->GC_11, pars->ZERO, pars->ZERO, w[76]); 
  FFV1_2(w[75], w[9], pars->GC_1, pars->ZERO, pars->ZERO, w[77]); 
  FFV1_2(w[75], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[78]); 
  FFV2_3_2(w[75], w[22], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[79]);
  FFV1P0_3(w[75], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[80]); 
  FFV1_2(w[75], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[81]); 
  FFV1P0_3(w[81], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[82]); 
  VVV1P0_1(w[80], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[83]); 
  FFV1_2(w[1], w[80], pars->GC_11, pars->ZERO, pars->ZERO, w[84]); 
  FFV1_1(w[7], w[80], pars->GC_11, pars->ZERO, pars->ZERO, w[85]); 
  FFV1_2(w[1], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[86]); 
  FFV1_2(w[86], w[9], pars->GC_1, pars->ZERO, pars->ZERO, w[87]); 
  FFV1P0_3(w[86], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[88]); 
  FFV2_3_2(w[86], w[22], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[89]);
  FFV1P0_3(w[86], w[36], pars->GC_11, pars->ZERO, pars->ZERO, w[90]); 
  FFV1_2(w[86], w[27], pars->GC_11, pars->ZERO, pars->ZERO, w[91]); 
  FFV1_2(w[86], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[92]); 
  FFV1P0_3(w[92], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[93]); 
  VVV1P0_1(w[88], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[94]); 
  FFV1_2(w[0], w[88], pars->GC_11, pars->ZERO, pars->ZERO, w[95]); 
  FFV1_1(w[6], w[88], pars->GC_11, pars->ZERO, pars->ZERO, w[96]); 
  FFV1_1(w[57], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[97]); 
  FFV1P0_3(w[0], w[97], pars->GC_11, pars->ZERO, pars->ZERO, w[98]); 
  VVV1P0_1(w[4], w[59], pars->GC_10, pars->ZERO, pars->ZERO, w[99]); 
  FFV1_1(w[11], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[100]); 
  FFV1_2(w[13], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[101]); 
  FFV1_1(w[23], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[102]); 
  FFV1_2(w[24], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[103]); 
  VVV1P0_1(w[4], w[31], pars->GC_10, pars->ZERO, pars->ZERO, w[104]); 
  FFV1_1(w[57], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[105]); 
  FFV1_1(w[36], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[106]); 
  FFV1P0_3(w[1], w[106], pars->GC_11, pars->ZERO, pars->ZERO, w[107]); 
  VVV1P0_1(w[4], w[38], pars->GC_10, pars->ZERO, pars->ZERO, w[108]); 
  FFV1_1(w[15], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[109]); 
  FFV1_2(w[17], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[110]); 
  FFV1_1(w[25], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[111]); 
  FFV1_2(w[26], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[112]); 
  VVV1P0_1(w[4], w[27], pars->GC_10, pars->ZERO, pars->ZERO, w[113]); 
  FFV1_1(w[36], w[27], pars->GC_11, pars->ZERO, pars->ZERO, w[114]); 
  FFV1_2(w[43], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[115]); 
  FFV1P0_3(w[115], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[116]); 
  VVV1P0_1(w[4], w[64], pars->GC_10, pars->ZERO, pars->ZERO, w[117]); 
  FFV1_2(w[43], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[118]); 
  FFV1_2(w[46], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[119]); 
  FFV1P0_3(w[119], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[120]); 
  VVV1P0_1(w[4], w[47], pars->GC_10, pars->ZERO, pars->ZERO, w[121]); 
  FFV1_2(w[46], w[27], pars->GC_11, pars->ZERO, pars->ZERO, w[122]); 
  VVV1P0_1(w[113], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[123]); 
  FFV1_2(w[1], w[113], pars->GC_11, pars->ZERO, pars->ZERO, w[124]); 
  FFV1_1(w[7], w[113], pars->GC_11, pars->ZERO, pars->ZERO, w[125]); 
  VVV1P0_1(w[4], w[74], pars->GC_10, pars->ZERO, pars->ZERO, w[126]); 
  FFV1_2(w[29], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[127]); 
  FFV1_1(w[30], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[128]); 
  VVVV1P0_1(w[4], w[5], w[27], pars->GC_12, pars->ZERO, pars->ZERO, w[129]); 
  VVVV3P0_1(w[4], w[5], w[27], pars->GC_12, pars->ZERO, pars->ZERO, w[130]); 
  VVVV4P0_1(w[4], w[5], w[27], pars->GC_12, pars->ZERO, pars->ZERO, w[131]); 
  VVV1P0_1(w[104], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[132]); 
  FFV1_2(w[0], w[104], pars->GC_11, pars->ZERO, pars->ZERO, w[133]); 
  FFV1_1(w[6], w[104], pars->GC_11, pars->ZERO, pars->ZERO, w[134]); 
  VVV1P0_1(w[4], w[55], pars->GC_10, pars->ZERO, pars->ZERO, w[135]); 
  FFV1_2(w[32], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[136]); 
  FFV1_1(w[34], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[137]); 
  VVVV1P0_1(w[4], w[5], w[31], pars->GC_12, pars->ZERO, pars->ZERO, w[138]); 
  VVVV3P0_1(w[4], w[5], w[31], pars->GC_12, pars->ZERO, pars->ZERO, w[139]); 
  VVVV4P0_1(w[4], w[5], w[31], pars->GC_12, pars->ZERO, pars->ZERO, w[140]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[1], w[11], w[12], pars->GC_11, amp[0]); 
  FFV1_0(w[13], w[7], w[12], pars->GC_11, amp[1]); 
  FFV1_0(w[0], w[15], w[16], pars->GC_11, amp[2]); 
  FFV1_0(w[17], w[6], w[16], pars->GC_11, amp[3]); 
  FFV1_0(w[1], w[11], w[19], pars->GC_11, amp[4]); 
  FFV1_0(w[13], w[7], w[19], pars->GC_11, amp[5]); 
  FFV1_0(w[0], w[15], w[21], pars->GC_11, amp[6]); 
  FFV1_0(w[17], w[6], w[21], pars->GC_11, amp[7]); 
  FFV1_0(w[1], w[23], w[12], pars->GC_11, amp[8]); 
  FFV1_0(w[24], w[7], w[12], pars->GC_11, amp[9]); 
  FFV1_0(w[0], w[25], w[16], pars->GC_11, amp[10]); 
  FFV1_0(w[26], w[6], w[16], pars->GC_11, amp[11]); 
  FFV1_0(w[1], w[23], w[19], pars->GC_11, amp[12]); 
  FFV1_0(w[24], w[7], w[19], pars->GC_11, amp[13]); 
  FFV1_0(w[0], w[25], w[21], pars->GC_11, amp[14]); 
  FFV1_0(w[26], w[6], w[21], pars->GC_11, amp[15]); 
  FFV1_0(w[1], w[11], w[28], pars->GC_11, amp[16]); 
  FFV1_0(w[13], w[7], w[28], pars->GC_11, amp[17]); 
  FFV1_0(w[13], w[14], w[27], pars->GC_11, amp[18]); 
  FFV1_0(w[29], w[14], w[9], pars->GC_1, amp[19]); 
  FFV1_0(w[20], w[11], w[27], pars->GC_11, amp[20]); 
  FFV1_0(w[20], w[30], w[9], pars->GC_1, amp[21]); 
  FFV1_0(w[29], w[11], w[8], pars->GC_11, amp[22]); 
  FFV1_0(w[13], w[30], w[8], pars->GC_11, amp[23]); 
  FFV1_0(w[1], w[23], w[28], pars->GC_11, amp[24]); 
  FFV1_0(w[24], w[7], w[28], pars->GC_11, amp[25]); 
  FFV1_0(w[24], w[14], w[27], pars->GC_11, amp[26]); 
  FFV2_3_0(w[29], w[14], w[22], -pars->GC_51, pars->GC_58, amp[27]); 
  FFV1_0(w[20], w[23], w[27], pars->GC_11, amp[28]); 
  FFV2_3_0(w[20], w[30], w[22], -pars->GC_51, pars->GC_58, amp[29]); 
  FFV1_0(w[29], w[23], w[8], pars->GC_11, amp[30]); 
  FFV1_0(w[24], w[30], w[8], pars->GC_11, amp[31]); 
  FFV1_0(w[17], w[10], w[31], pars->GC_11, amp[32]); 
  FFV1_0(w[32], w[10], w[9], pars->GC_1, amp[33]); 
  FFV1_0(w[0], w[15], w[33], pars->GC_11, amp[34]); 
  FFV1_0(w[17], w[6], w[33], pars->GC_11, amp[35]); 
  FFV1_0(w[18], w[15], w[31], pars->GC_11, amp[36]); 
  FFV1_0(w[18], w[34], w[9], pars->GC_1, amp[37]); 
  FFV1_0(w[32], w[15], w[8], pars->GC_11, amp[38]); 
  FFV1_0(w[17], w[34], w[8], pars->GC_11, amp[39]); 
  FFV1_0(w[26], w[10], w[31], pars->GC_11, amp[40]); 
  FFV2_3_0(w[32], w[10], w[22], -pars->GC_51, pars->GC_58, amp[41]); 
  FFV1_0(w[0], w[25], w[33], pars->GC_11, amp[42]); 
  FFV1_0(w[26], w[6], w[33], pars->GC_11, amp[43]); 
  FFV1_0(w[18], w[25], w[31], pars->GC_11, amp[44]); 
  FFV2_3_0(w[18], w[34], w[22], -pars->GC_51, pars->GC_58, amp[45]); 
  FFV1_0(w[32], w[25], w[8], pars->GC_11, amp[46]); 
  FFV1_0(w[26], w[34], w[8], pars->GC_11, amp[47]); 
  FFV1_0(w[0], w[37], w[38], pars->GC_11, amp[48]); 
  FFV1_0(w[1], w[40], w[39], pars->GC_11, amp[49]); 
  FFV1_0(w[13], w[36], w[39], pars->GC_11, amp[50]); 
  FFV1_0(w[17], w[35], w[38], pars->GC_11, amp[51]); 
  FFV1_0(w[0], w[41], w[38], pars->GC_11, amp[52]); 
  FFV1_0(w[1], w[42], w[39], pars->GC_11, amp[53]); 
  FFV1_0(w[24], w[36], w[39], pars->GC_11, amp[54]); 
  FFV1_0(w[26], w[35], w[38], pars->GC_11, amp[55]); 
  FFV1_0(w[1], w[11], w[44], pars->GC_11, amp[56]); 
  FFV1_0(w[13], w[7], w[44], pars->GC_11, amp[57]); 
  FFV1_0(w[1], w[23], w[44], pars->GC_11, amp[58]); 
  FFV1_0(w[24], w[7], w[44], pars->GC_11, amp[59]); 
  FFV1_0(w[43], w[37], w[31], pars->GC_11, amp[60]); 
  FFV1_0(w[43], w[45], w[9], pars->GC_1, amp[61]); 
  FFV1_0(w[43], w[41], w[31], pars->GC_11, amp[62]); 
  FFV2_3_0(w[43], w[45], w[22], -pars->GC_51, pars->GC_58, amp[63]); 
  FFV1_0(w[0], w[37], w[47], pars->GC_11, amp[64]); 
  FFV1_0(w[48], w[7], w[39], pars->GC_11, amp[65]); 
  FFV1_0(w[46], w[11], w[39], pars->GC_11, amp[66]); 
  FFV1_0(w[17], w[35], w[47], pars->GC_11, amp[67]); 
  FFV1_0(w[0], w[41], w[47], pars->GC_11, amp[68]); 
  FFV1_0(w[49], w[7], w[39], pars->GC_11, amp[69]); 
  FFV1_0(w[46], w[23], w[39], pars->GC_11, amp[70]); 
  FFV1_0(w[26], w[35], w[47], pars->GC_11, amp[71]); 
  FFV1_0(w[1], w[11], w[51], pars->GC_11, amp[72]); 
  FFV1_0(w[13], w[7], w[51], pars->GC_11, amp[73]); 
  FFV1_0(w[1], w[11], w[52], pars->GC_11, amp[74]); 
  FFV1_0(w[53], w[11], w[5], pars->GC_11, amp[75]); 
  FFV1_0(w[13], w[7], w[52], pars->GC_11, amp[76]); 
  FFV1_0(w[13], w[54], w[5], pars->GC_11, amp[77]); 
  FFV1_0(w[1], w[23], w[51], pars->GC_11, amp[78]); 
  FFV1_0(w[24], w[7], w[51], pars->GC_11, amp[79]); 
  FFV1_0(w[1], w[23], w[52], pars->GC_11, amp[80]); 
  FFV1_0(w[53], w[23], w[5], pars->GC_11, amp[81]); 
  FFV1_0(w[24], w[7], w[52], pars->GC_11, amp[82]); 
  FFV1_0(w[24], w[54], w[5], pars->GC_11, amp[83]); 
  FFV1_0(w[17], w[50], w[31], pars->GC_11, amp[84]); 
  FFV1_0(w[32], w[50], w[9], pars->GC_1, amp[85]); 
  FFV1_0(w[0], w[37], w[55], pars->GC_11, amp[86]); 
  FFV1_0(w[32], w[37], w[5], pars->GC_11, amp[87]); 
  FFV1_0(w[17], w[45], w[5], pars->GC_11, amp[88]); 
  FFV1_0(w[17], w[35], w[55], pars->GC_11, amp[89]); 
  FFV1_0(w[26], w[50], w[31], pars->GC_11, amp[90]); 
  FFV2_3_0(w[32], w[50], w[22], -pars->GC_51, pars->GC_58, amp[91]); 
  FFV1_0(w[0], w[41], w[55], pars->GC_11, amp[92]); 
  FFV1_0(w[32], w[41], w[5], pars->GC_11, amp[93]); 
  FFV1_0(w[26], w[45], w[5], pars->GC_11, amp[94]); 
  FFV1_0(w[26], w[35], w[55], pars->GC_11, amp[95]); 
  FFV1_0(w[1], w[58], w[59], pars->GC_11, amp[96]); 
  FFV1_0(w[0], w[61], w[60], pars->GC_11, amp[97]); 
  FFV1_0(w[17], w[57], w[60], pars->GC_11, amp[98]); 
  FFV1_0(w[13], w[56], w[59], pars->GC_11, amp[99]); 
  FFV1_0(w[1], w[62], w[59], pars->GC_11, amp[100]); 
  FFV1_0(w[0], w[63], w[60], pars->GC_11, amp[101]); 
  FFV1_0(w[26], w[57], w[60], pars->GC_11, amp[102]); 
  FFV1_0(w[24], w[56], w[59], pars->GC_11, amp[103]); 
  FFV1_0(w[1], w[58], w[64], pars->GC_11, amp[104]); 
  FFV1_0(w[65], w[6], w[60], pars->GC_11, amp[105]); 
  FFV1_0(w[43], w[15], w[60], pars->GC_11, amp[106]); 
  FFV1_0(w[13], w[56], w[64], pars->GC_11, amp[107]); 
  FFV1_0(w[1], w[62], w[64], pars->GC_11, amp[108]); 
  FFV1_0(w[66], w[6], w[60], pars->GC_11, amp[109]); 
  FFV1_0(w[43], w[25], w[60], pars->GC_11, amp[110]); 
  FFV1_0(w[24], w[56], w[64], pars->GC_11, amp[111]); 
  FFV1_0(w[0], w[15], w[67], pars->GC_11, amp[112]); 
  FFV1_0(w[17], w[6], w[67], pars->GC_11, amp[113]); 
  FFV1_0(w[0], w[25], w[67], pars->GC_11, amp[114]); 
  FFV1_0(w[26], w[6], w[67], pars->GC_11, amp[115]); 
  FFV1_0(w[46], w[58], w[27], pars->GC_11, amp[116]); 
  FFV1_0(w[46], w[68], w[9], pars->GC_1, amp[117]); 
  FFV1_0(w[46], w[62], w[27], pars->GC_11, amp[118]); 
  FFV2_3_0(w[46], w[68], w[22], -pars->GC_51, pars->GC_58, amp[119]); 
  FFV1_0(w[0], w[15], w[70], pars->GC_11, amp[120]); 
  FFV1_0(w[17], w[6], w[70], pars->GC_11, amp[121]); 
  FFV1_0(w[0], w[15], w[71], pars->GC_11, amp[122]); 
  FFV1_0(w[72], w[15], w[5], pars->GC_11, amp[123]); 
  FFV1_0(w[17], w[6], w[71], pars->GC_11, amp[124]); 
  FFV1_0(w[17], w[73], w[5], pars->GC_11, amp[125]); 
  FFV1_0(w[0], w[25], w[70], pars->GC_11, amp[126]); 
  FFV1_0(w[26], w[6], w[70], pars->GC_11, amp[127]); 
  FFV1_0(w[0], w[25], w[71], pars->GC_11, amp[128]); 
  FFV1_0(w[72], w[25], w[5], pars->GC_11, amp[129]); 
  FFV1_0(w[26], w[6], w[71], pars->GC_11, amp[130]); 
  FFV1_0(w[26], w[73], w[5], pars->GC_11, amp[131]); 
  FFV1_0(w[13], w[69], w[27], pars->GC_11, amp[132]); 
  FFV1_0(w[29], w[69], w[9], pars->GC_1, amp[133]); 
  FFV1_0(w[1], w[58], w[74], pars->GC_11, amp[134]); 
  FFV1_0(w[29], w[58], w[5], pars->GC_11, amp[135]); 
  FFV1_0(w[13], w[68], w[5], pars->GC_11, amp[136]); 
  FFV1_0(w[13], w[56], w[74], pars->GC_11, amp[137]); 
  FFV1_0(w[24], w[69], w[27], pars->GC_11, amp[138]); 
  FFV2_3_0(w[29], w[69], w[22], -pars->GC_51, pars->GC_58, amp[139]); 
  FFV1_0(w[1], w[62], w[74], pars->GC_11, amp[140]); 
  FFV1_0(w[29], w[62], w[5], pars->GC_11, amp[141]); 
  FFV1_0(w[24], w[68], w[5], pars->GC_11, amp[142]); 
  FFV1_0(w[24], w[56], w[74], pars->GC_11, amp[143]); 
  FFV1_0(w[1], w[11], w[76], pars->GC_11, amp[144]); 
  FFV1_0(w[13], w[7], w[76], pars->GC_11, amp[145]); 
  FFV1_0(w[1], w[23], w[76], pars->GC_11, amp[146]); 
  FFV1_0(w[24], w[7], w[76], pars->GC_11, amp[147]); 
  FFV1_0(w[77], w[57], w[31], pars->GC_11, amp[148]); 
  FFV1_0(w[78], w[57], w[9], pars->GC_1, amp[149]); 
  FFV1_0(w[79], w[57], w[31], pars->GC_11, amp[150]); 
  FFV2_3_0(w[78], w[57], w[22], -pars->GC_51, pars->GC_58, amp[151]); 
  FFV1_0(w[77], w[6], w[38], pars->GC_11, amp[152]); 
  FFV1_0(w[1], w[40], w[80], pars->GC_11, amp[153]); 
  FFV1_0(w[13], w[36], w[80], pars->GC_11, amp[154]); 
  FFV1_0(w[75], w[15], w[38], pars->GC_11, amp[155]); 
  FFV1_0(w[79], w[6], w[38], pars->GC_11, amp[156]); 
  FFV1_0(w[1], w[42], w[80], pars->GC_11, amp[157]); 
  FFV1_0(w[24], w[36], w[80], pars->GC_11, amp[158]); 
  FFV1_0(w[75], w[25], w[38], pars->GC_11, amp[159]); 
  FFV1_0(w[77], w[6], w[47], pars->GC_11, amp[160]); 
  FFV1_0(w[48], w[7], w[80], pars->GC_11, amp[161]); 
  FFV1_0(w[46], w[11], w[80], pars->GC_11, amp[162]); 
  FFV1_0(w[75], w[15], w[47], pars->GC_11, amp[163]); 
  FFV1_0(w[79], w[6], w[47], pars->GC_11, amp[164]); 
  FFV1_0(w[49], w[7], w[80], pars->GC_11, amp[165]); 
  FFV1_0(w[46], w[23], w[80], pars->GC_11, amp[166]); 
  FFV1_0(w[75], w[25], w[47], pars->GC_11, amp[167]); 
  FFV1_0(w[1], w[11], w[82], pars->GC_11, amp[168]); 
  FFV1_0(w[13], w[7], w[82], pars->GC_11, amp[169]); 
  FFV1_0(w[1], w[11], w[83], pars->GC_11, amp[170]); 
  FFV1_0(w[84], w[11], w[5], pars->GC_11, amp[171]); 
  FFV1_0(w[13], w[7], w[83], pars->GC_11, amp[172]); 
  FFV1_0(w[13], w[85], w[5], pars->GC_11, amp[173]); 
  FFV1_0(w[1], w[23], w[82], pars->GC_11, amp[174]); 
  FFV1_0(w[24], w[7], w[82], pars->GC_11, amp[175]); 
  FFV1_0(w[1], w[23], w[83], pars->GC_11, amp[176]); 
  FFV1_0(w[84], w[23], w[5], pars->GC_11, amp[177]); 
  FFV1_0(w[24], w[7], w[83], pars->GC_11, amp[178]); 
  FFV1_0(w[24], w[85], w[5], pars->GC_11, amp[179]); 
  FFV1_0(w[81], w[15], w[31], pars->GC_11, amp[180]); 
  FFV1_0(w[81], w[34], w[9], pars->GC_1, amp[181]); 
  FFV1_0(w[77], w[6], w[55], pars->GC_11, amp[182]); 
  FFV1_0(w[77], w[34], w[5], pars->GC_11, amp[183]); 
  FFV1_0(w[78], w[15], w[5], pars->GC_11, amp[184]); 
  FFV1_0(w[75], w[15], w[55], pars->GC_11, amp[185]); 
  FFV1_0(w[81], w[25], w[31], pars->GC_11, amp[186]); 
  FFV2_3_0(w[81], w[34], w[22], -pars->GC_51, pars->GC_58, amp[187]); 
  FFV1_0(w[79], w[6], w[55], pars->GC_11, amp[188]); 
  FFV1_0(w[79], w[34], w[5], pars->GC_11, amp[189]); 
  FFV1_0(w[78], w[25], w[5], pars->GC_11, amp[190]); 
  FFV1_0(w[75], w[25], w[55], pars->GC_11, amp[191]); 
  FFV1_0(w[87], w[7], w[59], pars->GC_11, amp[192]); 
  FFV1_0(w[0], w[61], w[88], pars->GC_11, amp[193]); 
  FFV1_0(w[17], w[57], w[88], pars->GC_11, amp[194]); 
  FFV1_0(w[86], w[11], w[59], pars->GC_11, amp[195]); 
  FFV1_0(w[89], w[7], w[59], pars->GC_11, amp[196]); 
  FFV1_0(w[0], w[63], w[88], pars->GC_11, amp[197]); 
  FFV1_0(w[26], w[57], w[88], pars->GC_11, amp[198]); 
  FFV1_0(w[86], w[23], w[59], pars->GC_11, amp[199]); 
  FFV1_0(w[0], w[15], w[90], pars->GC_11, amp[200]); 
  FFV1_0(w[17], w[6], w[90], pars->GC_11, amp[201]); 
  FFV1_0(w[0], w[25], w[90], pars->GC_11, amp[202]); 
  FFV1_0(w[26], w[6], w[90], pars->GC_11, amp[203]); 
  FFV1_0(w[87], w[36], w[27], pars->GC_11, amp[204]); 
  FFV1_0(w[91], w[36], w[9], pars->GC_1, amp[205]); 
  FFV1_0(w[89], w[36], w[27], pars->GC_11, amp[206]); 
  FFV2_3_0(w[91], w[36], w[22], -pars->GC_51, pars->GC_58, amp[207]); 
  FFV1_0(w[87], w[7], w[64], pars->GC_11, amp[208]); 
  FFV1_0(w[65], w[6], w[88], pars->GC_11, amp[209]); 
  FFV1_0(w[43], w[15], w[88], pars->GC_11, amp[210]); 
  FFV1_0(w[86], w[11], w[64], pars->GC_11, amp[211]); 
  FFV1_0(w[89], w[7], w[64], pars->GC_11, amp[212]); 
  FFV1_0(w[66], w[6], w[88], pars->GC_11, amp[213]); 
  FFV1_0(w[43], w[25], w[88], pars->GC_11, amp[214]); 
  FFV1_0(w[86], w[23], w[64], pars->GC_11, amp[215]); 
  FFV1_0(w[0], w[15], w[93], pars->GC_11, amp[216]); 
  FFV1_0(w[17], w[6], w[93], pars->GC_11, amp[217]); 
  FFV1_0(w[0], w[15], w[94], pars->GC_11, amp[218]); 
  FFV1_0(w[95], w[15], w[5], pars->GC_11, amp[219]); 
  FFV1_0(w[17], w[6], w[94], pars->GC_11, amp[220]); 
  FFV1_0(w[17], w[96], w[5], pars->GC_11, amp[221]); 
  FFV1_0(w[0], w[25], w[93], pars->GC_11, amp[222]); 
  FFV1_0(w[26], w[6], w[93], pars->GC_11, amp[223]); 
  FFV1_0(w[0], w[25], w[94], pars->GC_11, amp[224]); 
  FFV1_0(w[95], w[25], w[5], pars->GC_11, amp[225]); 
  FFV1_0(w[26], w[6], w[94], pars->GC_11, amp[226]); 
  FFV1_0(w[26], w[96], w[5], pars->GC_11, amp[227]); 
  FFV1_0(w[92], w[11], w[27], pars->GC_11, amp[228]); 
  FFV1_0(w[92], w[30], w[9], pars->GC_1, amp[229]); 
  FFV1_0(w[87], w[7], w[74], pars->GC_11, amp[230]); 
  FFV1_0(w[87], w[30], w[5], pars->GC_11, amp[231]); 
  FFV1_0(w[91], w[11], w[5], pars->GC_11, amp[232]); 
  FFV1_0(w[86], w[11], w[74], pars->GC_11, amp[233]); 
  FFV1_0(w[92], w[23], w[27], pars->GC_11, amp[234]); 
  FFV2_3_0(w[92], w[30], w[22], -pars->GC_51, pars->GC_58, amp[235]); 
  FFV1_0(w[89], w[7], w[74], pars->GC_11, amp[236]); 
  FFV1_0(w[89], w[30], w[5], pars->GC_11, amp[237]); 
  FFV1_0(w[91], w[23], w[5], pars->GC_11, amp[238]); 
  FFV1_0(w[86], w[23], w[74], pars->GC_11, amp[239]); 
  FFV1_0(w[1], w[11], w[98], pars->GC_11, amp[240]); 
  FFV1_0(w[13], w[7], w[98], pars->GC_11, amp[241]); 
  FFV1_0(w[1], w[11], w[99], pars->GC_11, amp[242]); 
  FFV1_0(w[1], w[100], w[59], pars->GC_11, amp[243]); 
  FFV1_0(w[13], w[7], w[99], pars->GC_11, amp[244]); 
  FFV1_0(w[101], w[7], w[59], pars->GC_11, amp[245]); 
  FFV1_0(w[1], w[23], w[98], pars->GC_11, amp[246]); 
  FFV1_0(w[24], w[7], w[98], pars->GC_11, amp[247]); 
  FFV1_0(w[1], w[23], w[99], pars->GC_11, amp[248]); 
  FFV1_0(w[1], w[102], w[59], pars->GC_11, amp[249]); 
  FFV1_0(w[24], w[7], w[99], pars->GC_11, amp[250]); 
  FFV1_0(w[103], w[7], w[59], pars->GC_11, amp[251]); 
  FFV1_0(w[17], w[97], w[31], pars->GC_11, amp[252]); 
  FFV1_0(w[32], w[97], w[9], pars->GC_1, amp[253]); 
  FFV1_0(w[0], w[61], w[104], pars->GC_11, amp[254]); 
  FFV1_0(w[17], w[57], w[104], pars->GC_11, amp[255]); 
  FFV1_0(w[32], w[61], w[4], pars->GC_11, amp[256]); 
  FFV1_0(w[17], w[105], w[4], pars->GC_11, amp[257]); 
  FFV1_0(w[26], w[97], w[31], pars->GC_11, amp[258]); 
  FFV2_3_0(w[32], w[97], w[22], -pars->GC_51, pars->GC_58, amp[259]); 
  FFV1_0(w[0], w[63], w[104], pars->GC_11, amp[260]); 
  FFV1_0(w[26], w[57], w[104], pars->GC_11, amp[261]); 
  FFV1_0(w[32], w[63], w[4], pars->GC_11, amp[262]); 
  FFV1_0(w[26], w[105], w[4], pars->GC_11, amp[263]); 
  FFV1_0(w[0], w[15], w[107], pars->GC_11, amp[264]); 
  FFV1_0(w[17], w[6], w[107], pars->GC_11, amp[265]); 
  FFV1_0(w[0], w[15], w[108], pars->GC_11, amp[266]); 
  FFV1_0(w[0], w[109], w[38], pars->GC_11, amp[267]); 
  FFV1_0(w[17], w[6], w[108], pars->GC_11, amp[268]); 
  FFV1_0(w[110], w[6], w[38], pars->GC_11, amp[269]); 
  FFV1_0(w[0], w[25], w[107], pars->GC_11, amp[270]); 
  FFV1_0(w[26], w[6], w[107], pars->GC_11, amp[271]); 
  FFV1_0(w[0], w[25], w[108], pars->GC_11, amp[272]); 
  FFV1_0(w[0], w[111], w[38], pars->GC_11, amp[273]); 
  FFV1_0(w[26], w[6], w[108], pars->GC_11, amp[274]); 
  FFV1_0(w[112], w[6], w[38], pars->GC_11, amp[275]); 
  FFV1_0(w[13], w[106], w[27], pars->GC_11, amp[276]); 
  FFV1_0(w[29], w[106], w[9], pars->GC_1, amp[277]); 
  FFV1_0(w[1], w[40], w[113], pars->GC_11, amp[278]); 
  FFV1_0(w[13], w[36], w[113], pars->GC_11, amp[279]); 
  FFV1_0(w[29], w[40], w[4], pars->GC_11, amp[280]); 
  FFV1_0(w[13], w[114], w[4], pars->GC_11, amp[281]); 
  FFV1_0(w[24], w[106], w[27], pars->GC_11, amp[282]); 
  FFV2_3_0(w[29], w[106], w[22], -pars->GC_51, pars->GC_58, amp[283]); 
  FFV1_0(w[1], w[42], w[113], pars->GC_11, amp[284]); 
  FFV1_0(w[24], w[36], w[113], pars->GC_11, amp[285]); 
  FFV1_0(w[29], w[42], w[4], pars->GC_11, amp[286]); 
  FFV1_0(w[24], w[114], w[4], pars->GC_11, amp[287]); 
  FFV1_0(w[1], w[11], w[116], pars->GC_11, amp[288]); 
  FFV1_0(w[13], w[7], w[116], pars->GC_11, amp[289]); 
  FFV1_0(w[1], w[11], w[117], pars->GC_11, amp[290]); 
  FFV1_0(w[1], w[100], w[64], pars->GC_11, amp[291]); 
  FFV1_0(w[13], w[7], w[117], pars->GC_11, amp[292]); 
  FFV1_0(w[101], w[7], w[64], pars->GC_11, amp[293]); 
  FFV1_0(w[1], w[23], w[116], pars->GC_11, amp[294]); 
  FFV1_0(w[24], w[7], w[116], pars->GC_11, amp[295]); 
  FFV1_0(w[1], w[23], w[117], pars->GC_11, amp[296]); 
  FFV1_0(w[1], w[102], w[64], pars->GC_11, amp[297]); 
  FFV1_0(w[24], w[7], w[117], pars->GC_11, amp[298]); 
  FFV1_0(w[103], w[7], w[64], pars->GC_11, amp[299]); 
  FFV1_0(w[115], w[15], w[31], pars->GC_11, amp[300]); 
  FFV1_0(w[115], w[34], w[9], pars->GC_1, amp[301]); 
  FFV1_0(w[65], w[6], w[104], pars->GC_11, amp[302]); 
  FFV1_0(w[43], w[15], w[104], pars->GC_11, amp[303]); 
  FFV1_0(w[65], w[34], w[4], pars->GC_11, amp[304]); 
  FFV1_0(w[118], w[15], w[4], pars->GC_11, amp[305]); 
  FFV1_0(w[115], w[25], w[31], pars->GC_11, amp[306]); 
  FFV2_3_0(w[115], w[34], w[22], -pars->GC_51, pars->GC_58, amp[307]); 
  FFV1_0(w[66], w[6], w[104], pars->GC_11, amp[308]); 
  FFV1_0(w[43], w[25], w[104], pars->GC_11, amp[309]); 
  FFV1_0(w[66], w[34], w[4], pars->GC_11, amp[310]); 
  FFV1_0(w[118], w[25], w[4], pars->GC_11, amp[311]); 
  FFV1_0(w[0], w[15], w[120], pars->GC_11, amp[312]); 
  FFV1_0(w[17], w[6], w[120], pars->GC_11, amp[313]); 
  FFV1_0(w[0], w[15], w[121], pars->GC_11, amp[314]); 
  FFV1_0(w[0], w[109], w[47], pars->GC_11, amp[315]); 
  FFV1_0(w[17], w[6], w[121], pars->GC_11, amp[316]); 
  FFV1_0(w[110], w[6], w[47], pars->GC_11, amp[317]); 
  FFV1_0(w[0], w[25], w[120], pars->GC_11, amp[318]); 
  FFV1_0(w[26], w[6], w[120], pars->GC_11, amp[319]); 
  FFV1_0(w[0], w[25], w[121], pars->GC_11, amp[320]); 
  FFV1_0(w[0], w[111], w[47], pars->GC_11, amp[321]); 
  FFV1_0(w[26], w[6], w[121], pars->GC_11, amp[322]); 
  FFV1_0(w[112], w[6], w[47], pars->GC_11, amp[323]); 
  FFV1_0(w[119], w[11], w[27], pars->GC_11, amp[324]); 
  FFV1_0(w[119], w[30], w[9], pars->GC_1, amp[325]); 
  FFV1_0(w[48], w[7], w[113], pars->GC_11, amp[326]); 
  FFV1_0(w[46], w[11], w[113], pars->GC_11, amp[327]); 
  FFV1_0(w[48], w[30], w[4], pars->GC_11, amp[328]); 
  FFV1_0(w[122], w[11], w[4], pars->GC_11, amp[329]); 
  FFV1_0(w[119], w[23], w[27], pars->GC_11, amp[330]); 
  FFV2_3_0(w[119], w[30], w[22], -pars->GC_51, pars->GC_58, amp[331]); 
  FFV1_0(w[49], w[7], w[113], pars->GC_11, amp[332]); 
  FFV1_0(w[46], w[23], w[113], pars->GC_11, amp[333]); 
  FFV1_0(w[49], w[30], w[4], pars->GC_11, amp[334]); 
  FFV1_0(w[122], w[23], w[4], pars->GC_11, amp[335]); 
  FFV1_0(w[1], w[11], w[123], pars->GC_11, amp[336]); 
  FFV1_0(w[124], w[11], w[5], pars->GC_11, amp[337]); 
  FFV1_0(w[13], w[7], w[123], pars->GC_11, amp[338]); 
  FFV1_0(w[13], w[125], w[5], pars->GC_11, amp[339]); 
  FFV1_0(w[1], w[11], w[126], pars->GC_11, amp[340]); 
  FFV1_0(w[1], w[100], w[74], pars->GC_11, amp[341]); 
  FFV1_0(w[13], w[7], w[126], pars->GC_11, amp[342]); 
  FFV1_0(w[101], w[7], w[74], pars->GC_11, amp[343]); 
  FFV1_0(w[29], w[100], w[5], pars->GC_11, amp[344]); 
  FFV1_0(w[127], w[11], w[5], pars->GC_11, amp[345]); 
  FFV1_0(w[101], w[30], w[5], pars->GC_11, amp[346]); 
  FFV1_0(w[13], w[128], w[5], pars->GC_11, amp[347]); 
  FFV1_0(w[1], w[11], w[129], pars->GC_11, amp[348]); 
  FFV1_0(w[1], w[11], w[130], pars->GC_11, amp[349]); 
  FFV1_0(w[1], w[11], w[131], pars->GC_11, amp[350]); 
  FFV1_0(w[13], w[7], w[129], pars->GC_11, amp[351]); 
  FFV1_0(w[13], w[7], w[130], pars->GC_11, amp[352]); 
  FFV1_0(w[13], w[7], w[131], pars->GC_11, amp[353]); 
  FFV1_0(w[1], w[23], w[123], pars->GC_11, amp[354]); 
  FFV1_0(w[124], w[23], w[5], pars->GC_11, amp[355]); 
  FFV1_0(w[24], w[7], w[123], pars->GC_11, amp[356]); 
  FFV1_0(w[24], w[125], w[5], pars->GC_11, amp[357]); 
  FFV1_0(w[1], w[23], w[126], pars->GC_11, amp[358]); 
  FFV1_0(w[1], w[102], w[74], pars->GC_11, amp[359]); 
  FFV1_0(w[24], w[7], w[126], pars->GC_11, amp[360]); 
  FFV1_0(w[103], w[7], w[74], pars->GC_11, amp[361]); 
  FFV1_0(w[29], w[102], w[5], pars->GC_11, amp[362]); 
  FFV1_0(w[127], w[23], w[5], pars->GC_11, amp[363]); 
  FFV1_0(w[103], w[30], w[5], pars->GC_11, amp[364]); 
  FFV1_0(w[24], w[128], w[5], pars->GC_11, amp[365]); 
  FFV1_0(w[1], w[23], w[129], pars->GC_11, amp[366]); 
  FFV1_0(w[1], w[23], w[130], pars->GC_11, amp[367]); 
  FFV1_0(w[1], w[23], w[131], pars->GC_11, amp[368]); 
  FFV1_0(w[24], w[7], w[129], pars->GC_11, amp[369]); 
  FFV1_0(w[24], w[7], w[130], pars->GC_11, amp[370]); 
  FFV1_0(w[24], w[7], w[131], pars->GC_11, amp[371]); 
  FFV1_0(w[0], w[15], w[132], pars->GC_11, amp[372]); 
  FFV1_0(w[133], w[15], w[5], pars->GC_11, amp[373]); 
  FFV1_0(w[17], w[6], w[132], pars->GC_11, amp[374]); 
  FFV1_0(w[17], w[134], w[5], pars->GC_11, amp[375]); 
  FFV1_0(w[0], w[15], w[135], pars->GC_11, amp[376]); 
  FFV1_0(w[0], w[109], w[55], pars->GC_11, amp[377]); 
  FFV1_0(w[17], w[6], w[135], pars->GC_11, amp[378]); 
  FFV1_0(w[110], w[6], w[55], pars->GC_11, amp[379]); 
  FFV1_0(w[32], w[109], w[5], pars->GC_11, amp[380]); 
  FFV1_0(w[136], w[15], w[5], pars->GC_11, amp[381]); 
  FFV1_0(w[110], w[34], w[5], pars->GC_11, amp[382]); 
  FFV1_0(w[17], w[137], w[5], pars->GC_11, amp[383]); 
  FFV1_0(w[0], w[15], w[138], pars->GC_11, amp[384]); 
  FFV1_0(w[0], w[15], w[139], pars->GC_11, amp[385]); 
  FFV1_0(w[0], w[15], w[140], pars->GC_11, amp[386]); 
  FFV1_0(w[17], w[6], w[138], pars->GC_11, amp[387]); 
  FFV1_0(w[17], w[6], w[139], pars->GC_11, amp[388]); 
  FFV1_0(w[17], w[6], w[140], pars->GC_11, amp[389]); 
  FFV1_0(w[0], w[25], w[132], pars->GC_11, amp[390]); 
  FFV1_0(w[133], w[25], w[5], pars->GC_11, amp[391]); 
  FFV1_0(w[26], w[6], w[132], pars->GC_11, amp[392]); 
  FFV1_0(w[26], w[134], w[5], pars->GC_11, amp[393]); 
  FFV1_0(w[0], w[25], w[135], pars->GC_11, amp[394]); 
  FFV1_0(w[0], w[111], w[55], pars->GC_11, amp[395]); 
  FFV1_0(w[26], w[6], w[135], pars->GC_11, amp[396]); 
  FFV1_0(w[112], w[6], w[55], pars->GC_11, amp[397]); 
  FFV1_0(w[32], w[111], w[5], pars->GC_11, amp[398]); 
  FFV1_0(w[136], w[25], w[5], pars->GC_11, amp[399]); 
  FFV1_0(w[112], w[34], w[5], pars->GC_11, amp[400]); 
  FFV1_0(w[26], w[137], w[5], pars->GC_11, amp[401]); 
  FFV1_0(w[0], w[25], w[138], pars->GC_11, amp[402]); 
  FFV1_0(w[0], w[25], w[139], pars->GC_11, amp[403]); 
  FFV1_0(w[0], w[25], w[140], pars->GC_11, amp[404]); 
  FFV1_0(w[26], w[6], w[138], pars->GC_11, amp[405]); 
  FFV1_0(w[26], w[6], w[139], pars->GC_11, amp[406]); 
  FFV1_0(w[26], w[6], w[140], pars->GC_11, amp[407]); 

}
double CPPProcess::matrix_ds_epemggds_no_h() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 408; 
  const int ncolor = 12; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}; 
  static const double cf[ncolor][ncolor] = {{48, 16, 16, 6, 0, 16, -2, 0, -6,
      -2, -2, 6}, {16, 48, 6, 16, 16, 0, 0, -2, -2, -6, 6, -2}, {16, 6, 48, 16,
      -2, 0, 0, 16, -2, 6, -6, -2}, {6, 16, 16, 48, 0, -2, 16, 0, 6, -2, -2,
      -6}, {0, 16, -2, 0, 48, 16, 16, 6, 0, -2, 16, 0}, {16, 0, 0, -2, 16, 48,
      6, 16, -2, 0, 0, 16}, {-2, 0, 0, 16, 16, 6, 48, 16, 16, 0, 0, -2}, {0,
      -2, 16, 0, 6, 16, 16, 48, 0, 16, -2, 0}, {-6, -2, -2, 6, 0, -2, 16, 0,
      48, 16, 16, 6}, {-2, -6, 6, -2, -2, 0, 0, 16, 16, 48, 6, 16}, {-2, 6, -6,
      -2, 16, 0, 0, -2, 16, 6, 48, 16}, {6, -2, -2, -6, 0, 16, -2, 0, 6, 16,
      16, 48}};

  // Calculate color flows
  jamp[0] = +1./2. * (-1./3. * std::complex<double> (0, 1) * amp[0] - 1./3. *
      std::complex<double> (0, 1) * amp[1] - 1./3. * std::complex<double> (0,
      1) * amp[4] - 1./3. * std::complex<double> (0, 1) * amp[5] - 1./3. *
      std::complex<double> (0, 1) * amp[8] - 1./3. * std::complex<double> (0,
      1) * amp[9] - 1./3. * std::complex<double> (0, 1) * amp[12] - 1./3. *
      std::complex<double> (0, 1) * amp[13] - 1./3. * std::complex<double> (0,
      1) * amp[32] - 1./3. * std::complex<double> (0, 1) * amp[33] - 1./3. *
      std::complex<double> (0, 1) * amp[36] - 1./3. * std::complex<double> (0,
      1) * amp[37] - 1./3. * std::complex<double> (0, 1) * amp[38] - 1./3. *
      std::complex<double> (0, 1) * amp[39] - 1./3. * std::complex<double> (0,
      1) * amp[40] - 1./3. * std::complex<double> (0, 1) * amp[41] - 1./3. *
      std::complex<double> (0, 1) * amp[44] - 1./3. * std::complex<double> (0,
      1) * amp[45] - 1./3. * std::complex<double> (0, 1) * amp[46] - 1./3. *
      std::complex<double> (0, 1) * amp[47] + 1./3. * amp[56] + 1./3. * amp[57]
      + 1./3. * amp[58] + 1./3. * amp[59] + 1./3. * amp[60] + 1./3. * amp[61] +
      1./3. * amp[62] + 1./3. * amp[63] + 1./3. * amp[72] + 1./3. * amp[73] +
      1./3. * amp[78] + 1./3. * amp[79] + 1./3. * amp[84] + 1./3. * amp[85] +
      1./3. * amp[87] + 1./3. * amp[88] + 1./3. * amp[90] + 1./3. * amp[91] +
      1./3. * amp[93] + 1./3. * amp[94] + 1./3. * amp[288] + 1./3. * amp[289] +
      1./3. * amp[294] + 1./3. * amp[295] + 1./3. * amp[300] + 1./3. * amp[301]
      + 1./3. * amp[304] + 1./3. * amp[305] + 1./3. * amp[306] + 1./3. *
      amp[307] + 1./3. * amp[310] + 1./3. * amp[311] + 1./3. * amp[380] + 1./3.
      * amp[383] + 1./3. * amp[398] + 1./3. * amp[401]);
  jamp[1] = +1./2. * (+std::complex<double> (0, 1) * amp[0] +
      std::complex<double> (0, 1) * amp[1] + std::complex<double> (0, 1) *
      amp[6] + std::complex<double> (0, 1) * amp[7] + std::complex<double> (0,
      1) * amp[8] + std::complex<double> (0, 1) * amp[9] + std::complex<double>
      (0, 1) * amp[14] + std::complex<double> (0, 1) * amp[15] - amp[16] -
      amp[17] + std::complex<double> (0, 1) * amp[20] + std::complex<double>
      (0, 1) * amp[21] + std::complex<double> (0, 1) * amp[23] - amp[24] -
      amp[25] + std::complex<double> (0, 1) * amp[28] + std::complex<double>
      (0, 1) * amp[29] + std::complex<double> (0, 1) * amp[31] +
      std::complex<double> (0, 1) * amp[32] + std::complex<double> (0, 1) *
      amp[33] + amp[34] + amp[35] + std::complex<double> (0, 1) * amp[38] +
      std::complex<double> (0, 1) * amp[40] + std::complex<double> (0, 1) *
      amp[41] + amp[42] + amp[43] + std::complex<double> (0, 1) * amp[46] -
      amp[64] - amp[65] - amp[66] - amp[67] - amp[68] - amp[69] - amp[70] -
      amp[71] - amp[72] - amp[73] + std::complex<double> (0, 1) * amp[74] +
      std::complex<double> (0, 1) * amp[76] - amp[77] - amp[78] - amp[79] +
      std::complex<double> (0, 1) * amp[80] + std::complex<double> (0, 1) *
      amp[82] - amp[83] - amp[84] - amp[85] + std::complex<double> (0, 1) *
      amp[86] - amp[87] + std::complex<double> (0, 1) * amp[89] - amp[90] -
      amp[91] + std::complex<double> (0, 1) * amp[92] - amp[93] +
      std::complex<double> (0, 1) * amp[95] - amp[312] - amp[313] +
      std::complex<double> (0, 1) * amp[314] - amp[315] + std::complex<double>
      (0, 1) * amp[316] - amp[318] - amp[319] + std::complex<double> (0, 1) *
      amp[320] - amp[321] + std::complex<double> (0, 1) * amp[322] - amp[324] -
      amp[325] - std::complex<double> (0, 1) * amp[326] - std::complex<double>
      (0, 1) * amp[327] - amp[328] - amp[330] - amp[331] - std::complex<double>
      (0, 1) * amp[332] - std::complex<double> (0, 1) * amp[333] - amp[334] -
      amp[336] - amp[338] - std::complex<double> (0, 1) * amp[339] - amp[347] +
      amp[350] + amp[349] + amp[353] + amp[352] - amp[354] - amp[356] -
      std::complex<double> (0, 1) * amp[357] - amp[365] + amp[368] + amp[367] +
      amp[371] + amp[370] + amp[376] + std::complex<double> (0, 1) * amp[377] +
      amp[378] - amp[380] - amp[386] + amp[384] - amp[389] + amp[387] +
      amp[394] + std::complex<double> (0, 1) * amp[395] + amp[396] - amp[398] -
      amp[404] + amp[402] - amp[407] + amp[405]);
  jamp[2] = +1./2. * (+std::complex<double> (0, 1) * amp[2] +
      std::complex<double> (0, 1) * amp[3] + std::complex<double> (0, 1) *
      amp[4] + std::complex<double> (0, 1) * amp[5] + std::complex<double> (0,
      1) * amp[10] + std::complex<double> (0, 1) * amp[11] +
      std::complex<double> (0, 1) * amp[12] + std::complex<double> (0, 1) *
      amp[13] + amp[16] + amp[17] + std::complex<double> (0, 1) * amp[18] +
      std::complex<double> (0, 1) * amp[19] + std::complex<double> (0, 1) *
      amp[22] + amp[24] + amp[25] + std::complex<double> (0, 1) * amp[26] +
      std::complex<double> (0, 1) * amp[27] + std::complex<double> (0, 1) *
      amp[30] - amp[34] - amp[35] + std::complex<double> (0, 1) * amp[36] +
      std::complex<double> (0, 1) * amp[37] + std::complex<double> (0, 1) *
      amp[39] - amp[42] - amp[43] + std::complex<double> (0, 1) * amp[44] +
      std::complex<double> (0, 1) * amp[45] + std::complex<double> (0, 1) *
      amp[47] - amp[104] - amp[105] - amp[106] - amp[107] - amp[108] - amp[109]
      - amp[110] - amp[111] - amp[120] - amp[121] + std::complex<double> (0, 1)
      * amp[122] + std::complex<double> (0, 1) * amp[124] - amp[125] - amp[126]
      - amp[127] + std::complex<double> (0, 1) * amp[128] +
      std::complex<double> (0, 1) * amp[130] - amp[131] - amp[132] - amp[133] +
      std::complex<double> (0, 1) * amp[134] - amp[135] + std::complex<double>
      (0, 1) * amp[137] - amp[138] - amp[139] + std::complex<double> (0, 1) *
      amp[140] - amp[141] + std::complex<double> (0, 1) * amp[143] - amp[288] -
      amp[289] + std::complex<double> (0, 1) * amp[290] - amp[291] +
      std::complex<double> (0, 1) * amp[292] - amp[294] - amp[295] +
      std::complex<double> (0, 1) * amp[296] - amp[297] + std::complex<double>
      (0, 1) * amp[298] - amp[300] - amp[301] - std::complex<double> (0, 1) *
      amp[302] - std::complex<double> (0, 1) * amp[303] - amp[304] - amp[306] -
      amp[307] - std::complex<double> (0, 1) * amp[308] - std::complex<double>
      (0, 1) * amp[309] - amp[310] + amp[340] + std::complex<double> (0, 1) *
      amp[341] + amp[342] - amp[344] - amp[350] + amp[348] - amp[353] +
      amp[351] + amp[358] + std::complex<double> (0, 1) * amp[359] + amp[360] -
      amp[362] - amp[368] + amp[366] - amp[371] + amp[369] - amp[372] -
      amp[374] - std::complex<double> (0, 1) * amp[375] - amp[383] + amp[386] +
      amp[385] + amp[389] + amp[388] - amp[390] - amp[392] -
      std::complex<double> (0, 1) * amp[393] - amp[401] + amp[404] + amp[403] +
      amp[407] + amp[406]);
  jamp[3] = +1./2. * (-1./3. * std::complex<double> (0, 1) * amp[2] - 1./3. *
      std::complex<double> (0, 1) * amp[3] - 1./3. * std::complex<double> (0,
      1) * amp[6] - 1./3. * std::complex<double> (0, 1) * amp[7] - 1./3. *
      std::complex<double> (0, 1) * amp[10] - 1./3. * std::complex<double> (0,
      1) * amp[11] - 1./3. * std::complex<double> (0, 1) * amp[14] - 1./3. *
      std::complex<double> (0, 1) * amp[15] - 1./3. * std::complex<double> (0,
      1) * amp[18] - 1./3. * std::complex<double> (0, 1) * amp[19] - 1./3. *
      std::complex<double> (0, 1) * amp[20] - 1./3. * std::complex<double> (0,
      1) * amp[21] - 1./3. * std::complex<double> (0, 1) * amp[22] - 1./3. *
      std::complex<double> (0, 1) * amp[23] - 1./3. * std::complex<double> (0,
      1) * amp[26] - 1./3. * std::complex<double> (0, 1) * amp[27] - 1./3. *
      std::complex<double> (0, 1) * amp[28] - 1./3. * std::complex<double> (0,
      1) * amp[29] - 1./3. * std::complex<double> (0, 1) * amp[30] - 1./3. *
      std::complex<double> (0, 1) * amp[31] + 1./3. * amp[112] + 1./3. *
      amp[113] + 1./3. * amp[114] + 1./3. * amp[115] + 1./3. * amp[116] + 1./3.
      * amp[117] + 1./3. * amp[118] + 1./3. * amp[119] + 1./3. * amp[120] +
      1./3. * amp[121] + 1./3. * amp[126] + 1./3. * amp[127] + 1./3. * amp[132]
      + 1./3. * amp[133] + 1./3. * amp[135] + 1./3. * amp[136] + 1./3. *
      amp[138] + 1./3. * amp[139] + 1./3. * amp[141] + 1./3. * amp[142] + 1./3.
      * amp[312] + 1./3. * amp[313] + 1./3. * amp[318] + 1./3. * amp[319] +
      1./3. * amp[324] + 1./3. * amp[325] + 1./3. * amp[328] + 1./3. * amp[329]
      + 1./3. * amp[330] + 1./3. * amp[331] + 1./3. * amp[334] + 1./3. *
      amp[335] + 1./3. * amp[344] + 1./3. * amp[347] + 1./3. * amp[362] + 1./3.
      * amp[365]);
  jamp[4] = +1./2. * (+1./3. * amp[48] + 1./3. * amp[49] + 1./3. * amp[50] +
      1./3. * amp[51] + 1./3. * amp[52] + 1./3. * amp[53] + 1./3. * amp[54] +
      1./3. * amp[55] + 1./3. * amp[64] + 1./3. * amp[65] + 1./3. * amp[66] +
      1./3. * amp[67] + 1./3. * amp[68] + 1./3. * amp[69] + 1./3. * amp[70] +
      1./3. * amp[71] + 1./3. * amp[75] + 1./3. * amp[77] + 1./3. * amp[81] +
      1./3. * amp[83] + 1./3. * amp[152] + 1./3. * amp[153] + 1./3. * amp[154]
      + 1./3. * amp[155] + 1./3. * amp[156] + 1./3. * amp[157] + 1./3. *
      amp[158] + 1./3. * amp[159] + 1./3. * amp[160] + 1./3. * amp[161] + 1./3.
      * amp[162] + 1./3. * amp[163] + 1./3. * amp[164] + 1./3. * amp[165] +
      1./3. * amp[166] + 1./3. * amp[167] + 1./3. * amp[171] + 1./3. * amp[173]
      + 1./3. * amp[177] + 1./3. * amp[179] + 1./3. * amp[267] + 1./3. *
      amp[269] + 1./3. * amp[273] + 1./3. * amp[275] + 1./3. * amp[315] + 1./3.
      * amp[317] + 1./3. * amp[321] + 1./3. * amp[323]);
  jamp[5] = +1./2. * (-amp[48] - amp[49] - amp[50] - amp[51] - amp[52] -
      amp[53] - amp[54] - amp[55] - amp[56] - amp[57] - amp[58] - amp[59] -
      amp[60] - amp[61] - amp[62] - amp[63] - std::complex<double> (0, 1) *
      amp[74] - amp[75] - std::complex<double> (0, 1) * amp[76] -
      std::complex<double> (0, 1) * amp[80] - amp[81] - std::complex<double>
      (0, 1) * amp[82] - std::complex<double> (0, 1) * amp[86] - amp[88] -
      std::complex<double> (0, 1) * amp[89] - std::complex<double> (0, 1) *
      amp[92] - amp[94] - std::complex<double> (0, 1) * amp[95] - amp[200] -
      amp[201] - amp[202] - amp[203] - amp[204] - amp[205] - amp[206] -
      amp[207] - amp[208] - amp[209] - amp[210] - amp[211] - amp[212] -
      amp[213] - amp[214] - amp[215] + std::complex<double> (0, 1) * amp[218] +
      std::complex<double> (0, 1) * amp[220] - amp[221] + std::complex<double>
      (0, 1) * amp[224] + std::complex<double> (0, 1) * amp[226] - amp[227] +
      std::complex<double> (0, 1) * amp[230] - amp[232] + std::complex<double>
      (0, 1) * amp[233] + std::complex<double> (0, 1) * amp[236] - amp[238] +
      std::complex<double> (0, 1) * amp[239] + std::complex<double> (0, 1) *
      amp[266] - amp[267] + std::complex<double> (0, 1) * amp[268] +
      std::complex<double> (0, 1) * amp[272] - amp[273] + std::complex<double>
      (0, 1) * amp[274] - std::complex<double> (0, 1) * amp[278] -
      std::complex<double> (0, 1) * amp[279] - amp[281] - std::complex<double>
      (0, 1) * amp[284] - std::complex<double> (0, 1) * amp[285] - amp[287] -
      std::complex<double> (0, 1) * amp[290] - std::complex<double> (0, 1) *
      amp[292] - amp[293] - std::complex<double> (0, 1) * amp[296] -
      std::complex<double> (0, 1) * amp[298] - amp[299] + std::complex<double>
      (0, 1) * amp[302] + std::complex<double> (0, 1) * amp[303] - amp[305] +
      std::complex<double> (0, 1) * amp[308] + std::complex<double> (0, 1) *
      amp[309] - amp[311] + amp[336] - std::complex<double> (0, 1) * amp[337] +
      amp[338] - amp[340] - amp[342] + std::complex<double> (0, 1) * amp[343] -
      amp[348] - amp[349] - amp[351] - amp[352] + amp[354] -
      std::complex<double> (0, 1) * amp[355] + amp[356] - amp[358] - amp[360] +
      std::complex<double> (0, 1) * amp[361] - amp[366] - amp[367] - amp[369] -
      amp[370] + amp[372] + amp[374] + std::complex<double> (0, 1) * amp[375] -
      amp[376] - std::complex<double> (0, 1) * amp[377] - amp[378] - amp[384] -
      amp[385] - amp[387] - amp[388] + amp[390] + amp[392] +
      std::complex<double> (0, 1) * amp[393] - amp[394] - std::complex<double>
      (0, 1) * amp[395] - amp[396] - amp[402] - amp[403] - amp[405] - amp[406]);
  jamp[6] = +1./2. * (-amp[96] - amp[97] - amp[98] - amp[99] - amp[100] -
      amp[101] - amp[102] - amp[103] - amp[112] - amp[113] - amp[114] -
      amp[115] - amp[116] - amp[117] - amp[118] - amp[119] -
      std::complex<double> (0, 1) * amp[122] - amp[123] - std::complex<double>
      (0, 1) * amp[124] - std::complex<double> (0, 1) * amp[128] - amp[129] -
      std::complex<double> (0, 1) * amp[130] - std::complex<double> (0, 1) *
      amp[134] - amp[136] - std::complex<double> (0, 1) * amp[137] -
      std::complex<double> (0, 1) * amp[140] - amp[142] - std::complex<double>
      (0, 1) * amp[143] - amp[144] - amp[145] - amp[146] - amp[147] - amp[148]
      - amp[149] - amp[150] - amp[151] - amp[160] - amp[161] - amp[162] -
      amp[163] - amp[164] - amp[165] - amp[166] - amp[167] +
      std::complex<double> (0, 1) * amp[170] + std::complex<double> (0, 1) *
      amp[172] - amp[173] + std::complex<double> (0, 1) * amp[176] +
      std::complex<double> (0, 1) * amp[178] - amp[179] + std::complex<double>
      (0, 1) * amp[182] - amp[184] + std::complex<double> (0, 1) * amp[185] +
      std::complex<double> (0, 1) * amp[188] - amp[190] + std::complex<double>
      (0, 1) * amp[191] + std::complex<double> (0, 1) * amp[242] - amp[243] +
      std::complex<double> (0, 1) * amp[244] + std::complex<double> (0, 1) *
      amp[248] - amp[249] + std::complex<double> (0, 1) * amp[250] -
      std::complex<double> (0, 1) * amp[254] - std::complex<double> (0, 1) *
      amp[255] - amp[257] - std::complex<double> (0, 1) * amp[260] -
      std::complex<double> (0, 1) * amp[261] - amp[263] - std::complex<double>
      (0, 1) * amp[314] - std::complex<double> (0, 1) * amp[316] - amp[317] -
      std::complex<double> (0, 1) * amp[320] - std::complex<double> (0, 1) *
      amp[322] - amp[323] + std::complex<double> (0, 1) * amp[326] +
      std::complex<double> (0, 1) * amp[327] - amp[329] + std::complex<double>
      (0, 1) * amp[332] + std::complex<double> (0, 1) * amp[333] - amp[335] +
      amp[336] + amp[338] + std::complex<double> (0, 1) * amp[339] - amp[340] -
      std::complex<double> (0, 1) * amp[341] - amp[342] - amp[348] - amp[349] -
      amp[351] - amp[352] + amp[354] + amp[356] + std::complex<double> (0, 1) *
      amp[357] - amp[358] - std::complex<double> (0, 1) * amp[359] - amp[360] -
      amp[366] - amp[367] - amp[369] - amp[370] + amp[372] -
      std::complex<double> (0, 1) * amp[373] + amp[374] - amp[376] - amp[378] +
      std::complex<double> (0, 1) * amp[379] - amp[384] - amp[385] - amp[387] -
      amp[388] + amp[390] - std::complex<double> (0, 1) * amp[391] + amp[392] -
      amp[394] - amp[396] + std::complex<double> (0, 1) * amp[397] - amp[402] -
      amp[403] - amp[405] - amp[406]);
  jamp[7] = +1./2. * (+1./3. * amp[96] + 1./3. * amp[97] + 1./3. * amp[98] +
      1./3. * amp[99] + 1./3. * amp[100] + 1./3. * amp[101] + 1./3. * amp[102]
      + 1./3. * amp[103] + 1./3. * amp[104] + 1./3. * amp[105] + 1./3. *
      amp[106] + 1./3. * amp[107] + 1./3. * amp[108] + 1./3. * amp[109] + 1./3.
      * amp[110] + 1./3. * amp[111] + 1./3. * amp[123] + 1./3. * amp[125] +
      1./3. * amp[129] + 1./3. * amp[131] + 1./3. * amp[192] + 1./3. * amp[193]
      + 1./3. * amp[194] + 1./3. * amp[195] + 1./3. * amp[196] + 1./3. *
      amp[197] + 1./3. * amp[198] + 1./3. * amp[199] + 1./3. * amp[208] + 1./3.
      * amp[209] + 1./3. * amp[210] + 1./3. * amp[211] + 1./3. * amp[212] +
      1./3. * amp[213] + 1./3. * amp[214] + 1./3. * amp[215] + 1./3. * amp[219]
      + 1./3. * amp[221] + 1./3. * amp[225] + 1./3. * amp[227] + 1./3. *
      amp[243] + 1./3. * amp[245] + 1./3. * amp[249] + 1./3. * amp[251] + 1./3.
      * amp[291] + 1./3. * amp[293] + 1./3. * amp[297] + 1./3. * amp[299]);
  jamp[8] = +1./2. * (+1./3. * std::complex<double> (0, 1) * amp[0] + 1./3. *
      std::complex<double> (0, 1) * amp[1] + 1./3. * std::complex<double> (0,
      1) * amp[4] + 1./3. * std::complex<double> (0, 1) * amp[5] + 1./3. *
      std::complex<double> (0, 1) * amp[8] + 1./3. * std::complex<double> (0,
      1) * amp[9] + 1./3. * std::complex<double> (0, 1) * amp[12] + 1./3. *
      std::complex<double> (0, 1) * amp[13] + 1./3. * std::complex<double> (0,
      1) * amp[32] + 1./3. * std::complex<double> (0, 1) * amp[33] + 1./3. *
      std::complex<double> (0, 1) * amp[36] + 1./3. * std::complex<double> (0,
      1) * amp[37] + 1./3. * std::complex<double> (0, 1) * amp[38] + 1./3. *
      std::complex<double> (0, 1) * amp[39] + 1./3. * std::complex<double> (0,
      1) * amp[40] + 1./3. * std::complex<double> (0, 1) * amp[41] + 1./3. *
      std::complex<double> (0, 1) * amp[44] + 1./3. * std::complex<double> (0,
      1) * amp[45] + 1./3. * std::complex<double> (0, 1) * amp[46] + 1./3. *
      std::complex<double> (0, 1) * amp[47] + 1./3. * amp[144] + 1./3. *
      amp[145] + 1./3. * amp[146] + 1./3. * amp[147] + 1./3. * amp[148] + 1./3.
      * amp[149] + 1./3. * amp[150] + 1./3. * amp[151] + 1./3. * amp[168] +
      1./3. * amp[169] + 1./3. * amp[174] + 1./3. * amp[175] + 1./3. * amp[180]
      + 1./3. * amp[181] + 1./3. * amp[183] + 1./3. * amp[184] + 1./3. *
      amp[186] + 1./3. * amp[187] + 1./3. * amp[189] + 1./3. * amp[190] + 1./3.
      * amp[240] + 1./3. * amp[241] + 1./3. * amp[246] + 1./3. * amp[247] +
      1./3. * amp[252] + 1./3. * amp[253] + 1./3. * amp[256] + 1./3. * amp[257]
      + 1./3. * amp[258] + 1./3. * amp[259] + 1./3. * amp[262] + 1./3. *
      amp[263] + 1./3. * amp[381] + 1./3. * amp[382] + 1./3. * amp[399] + 1./3.
      * amp[400]);
  jamp[9] = +1./2. * (-std::complex<double> (0, 1) * amp[0] -
      std::complex<double> (0, 1) * amp[1] - std::complex<double> (0, 1) *
      amp[6] - std::complex<double> (0, 1) * amp[7] - std::complex<double> (0,
      1) * amp[8] - std::complex<double> (0, 1) * amp[9] - std::complex<double>
      (0, 1) * amp[14] - std::complex<double> (0, 1) * amp[15] + amp[16] +
      amp[17] - std::complex<double> (0, 1) * amp[20] - std::complex<double>
      (0, 1) * amp[21] - std::complex<double> (0, 1) * amp[23] + amp[24] +
      amp[25] - std::complex<double> (0, 1) * amp[28] - std::complex<double>
      (0, 1) * amp[29] - std::complex<double> (0, 1) * amp[31] -
      std::complex<double> (0, 1) * amp[32] - std::complex<double> (0, 1) *
      amp[33] - amp[34] - amp[35] - std::complex<double> (0, 1) * amp[38] -
      std::complex<double> (0, 1) * amp[40] - std::complex<double> (0, 1) *
      amp[41] - amp[42] - amp[43] - std::complex<double> (0, 1) * amp[46] -
      amp[192] - amp[193] - amp[194] - amp[195] - amp[196] - amp[197] -
      amp[198] - amp[199] - amp[216] - amp[217] - std::complex<double> (0, 1) *
      amp[218] - amp[219] - std::complex<double> (0, 1) * amp[220] - amp[222] -
      amp[223] - std::complex<double> (0, 1) * amp[224] - amp[225] -
      std::complex<double> (0, 1) * amp[226] - amp[228] - amp[229] -
      std::complex<double> (0, 1) * amp[230] - amp[231] - std::complex<double>
      (0, 1) * amp[233] - amp[234] - amp[235] - std::complex<double> (0, 1) *
      amp[236] - amp[237] - std::complex<double> (0, 1) * amp[239] - amp[240] -
      amp[241] - std::complex<double> (0, 1) * amp[242] - std::complex<double>
      (0, 1) * amp[244] - amp[245] - amp[246] - amp[247] - std::complex<double>
      (0, 1) * amp[248] - std::complex<double> (0, 1) * amp[250] - amp[251] -
      amp[252] - amp[253] + std::complex<double> (0, 1) * amp[254] +
      std::complex<double> (0, 1) * amp[255] - amp[256] - amp[258] - amp[259] +
      std::complex<double> (0, 1) * amp[260] + std::complex<double> (0, 1) *
      amp[261] - amp[262] + amp[340] + amp[342] - std::complex<double> (0, 1) *
      amp[343] - amp[346] - amp[350] + amp[348] - amp[353] + amp[351] +
      amp[358] + amp[360] - std::complex<double> (0, 1) * amp[361] - amp[364] -
      amp[368] + amp[366] - amp[371] + amp[369] - amp[372] +
      std::complex<double> (0, 1) * amp[373] - amp[374] - amp[381] + amp[386] +
      amp[385] + amp[389] + amp[388] - amp[390] + std::complex<double> (0, 1) *
      amp[391] - amp[392] - amp[399] + amp[404] + amp[403] + amp[407] +
      amp[406]);
  jamp[10] = +1./2. * (-std::complex<double> (0, 1) * amp[2] -
      std::complex<double> (0, 1) * amp[3] - std::complex<double> (0, 1) *
      amp[4] - std::complex<double> (0, 1) * amp[5] - std::complex<double> (0,
      1) * amp[10] - std::complex<double> (0, 1) * amp[11] -
      std::complex<double> (0, 1) * amp[12] - std::complex<double> (0, 1) *
      amp[13] - amp[16] - amp[17] - std::complex<double> (0, 1) * amp[18] -
      std::complex<double> (0, 1) * amp[19] - std::complex<double> (0, 1) *
      amp[22] - amp[24] - amp[25] - std::complex<double> (0, 1) * amp[26] -
      std::complex<double> (0, 1) * amp[27] - std::complex<double> (0, 1) *
      amp[30] + amp[34] + amp[35] - std::complex<double> (0, 1) * amp[36] -
      std::complex<double> (0, 1) * amp[37] - std::complex<double> (0, 1) *
      amp[39] + amp[42] + amp[43] - std::complex<double> (0, 1) * amp[44] -
      std::complex<double> (0, 1) * amp[45] - std::complex<double> (0, 1) *
      amp[47] - amp[152] - amp[153] - amp[154] - amp[155] - amp[156] - amp[157]
      - amp[158] - amp[159] - amp[168] - amp[169] - std::complex<double> (0, 1)
      * amp[170] - amp[171] - std::complex<double> (0, 1) * amp[172] - amp[174]
      - amp[175] - std::complex<double> (0, 1) * amp[176] - amp[177] -
      std::complex<double> (0, 1) * amp[178] - amp[180] - amp[181] -
      std::complex<double> (0, 1) * amp[182] - amp[183] - std::complex<double>
      (0, 1) * amp[185] - amp[186] - amp[187] - std::complex<double> (0, 1) *
      amp[188] - amp[189] - std::complex<double> (0, 1) * amp[191] - amp[264] -
      amp[265] - std::complex<double> (0, 1) * amp[266] - std::complex<double>
      (0, 1) * amp[268] - amp[269] - amp[270] - amp[271] - std::complex<double>
      (0, 1) * amp[272] - std::complex<double> (0, 1) * amp[274] - amp[275] -
      amp[276] - amp[277] + std::complex<double> (0, 1) * amp[278] +
      std::complex<double> (0, 1) * amp[279] - amp[280] - amp[282] - amp[283] +
      std::complex<double> (0, 1) * amp[284] + std::complex<double> (0, 1) *
      amp[285] - amp[286] - amp[336] + std::complex<double> (0, 1) * amp[337] -
      amp[338] - amp[345] + amp[350] + amp[349] + amp[353] + amp[352] -
      amp[354] + std::complex<double> (0, 1) * amp[355] - amp[356] - amp[363] +
      amp[368] + amp[367] + amp[371] + amp[370] + amp[376] + amp[378] -
      std::complex<double> (0, 1) * amp[379] - amp[382] - amp[386] + amp[384] -
      amp[389] + amp[387] + amp[394] + amp[396] - std::complex<double> (0, 1) *
      amp[397] - amp[400] - amp[404] + amp[402] - amp[407] + amp[405]);
  jamp[11] = +1./2. * (+1./3. * std::complex<double> (0, 1) * amp[2] + 1./3. *
      std::complex<double> (0, 1) * amp[3] + 1./3. * std::complex<double> (0,
      1) * amp[6] + 1./3. * std::complex<double> (0, 1) * amp[7] + 1./3. *
      std::complex<double> (0, 1) * amp[10] + 1./3. * std::complex<double> (0,
      1) * amp[11] + 1./3. * std::complex<double> (0, 1) * amp[14] + 1./3. *
      std::complex<double> (0, 1) * amp[15] + 1./3. * std::complex<double> (0,
      1) * amp[18] + 1./3. * std::complex<double> (0, 1) * amp[19] + 1./3. *
      std::complex<double> (0, 1) * amp[20] + 1./3. * std::complex<double> (0,
      1) * amp[21] + 1./3. * std::complex<double> (0, 1) * amp[22] + 1./3. *
      std::complex<double> (0, 1) * amp[23] + 1./3. * std::complex<double> (0,
      1) * amp[26] + 1./3. * std::complex<double> (0, 1) * amp[27] + 1./3. *
      std::complex<double> (0, 1) * amp[28] + 1./3. * std::complex<double> (0,
      1) * amp[29] + 1./3. * std::complex<double> (0, 1) * amp[30] + 1./3. *
      std::complex<double> (0, 1) * amp[31] + 1./3. * amp[200] + 1./3. *
      amp[201] + 1./3. * amp[202] + 1./3. * amp[203] + 1./3. * amp[204] + 1./3.
      * amp[205] + 1./3. * amp[206] + 1./3. * amp[207] + 1./3. * amp[216] +
      1./3. * amp[217] + 1./3. * amp[222] + 1./3. * amp[223] + 1./3. * amp[228]
      + 1./3. * amp[229] + 1./3. * amp[231] + 1./3. * amp[232] + 1./3. *
      amp[234] + 1./3. * amp[235] + 1./3. * amp[237] + 1./3. * amp[238] + 1./3.
      * amp[264] + 1./3. * amp[265] + 1./3. * amp[270] + 1./3. * amp[271] +
      1./3. * amp[276] + 1./3. * amp[277] + 1./3. * amp[280] + 1./3. * amp[281]
      + 1./3. * amp[282] + 1./3. * amp[283] + 1./3. * amp[286] + 1./3. *
      amp[287] + 1./3. * amp[345] + 1./3. * amp[346] + 1./3. * amp[363] + 1./3.
      * amp[364]);

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}



