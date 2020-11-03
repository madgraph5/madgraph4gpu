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
// Process: s s~ > ta+ ta- u d u~ d~ WEIGHTED<=8 / h
// Process: s s~ > ta+ ta- c d c~ d~ WEIGHTED<=8 / h
// Process: b b~ > ta+ ta- u d u~ d~ WEIGHTED<=8 / h
// Process: b b~ > ta+ ta- u s u~ s~ WEIGHTED<=8 / h
// Process: b b~ > ta+ ta- c d c~ d~ WEIGHTED<=8 / h
// Process: b b~ > ta+ ta- c s c~ s~ WEIGHTED<=8 / h

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
  mME.push_back(pars->mdl_MTA); 
  mME.push_back(pars->mdl_MTA); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[6]; 
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
  for(int i = 0; i < 6; i++ )
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
  const int denominators[nprocesses] = {36, 36}; 

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
        t[0] = matrix_ssx_taptamuduxdx_no_h(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_ssx_taptamuduxdx_no_h(); 
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
      t[0] = matrix_ssx_taptamuduxdx_no_h(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_ssx_taptamuduxdx_no_h(); 
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
  if(id1 == -3 && id2 == 3)
  {
    // Add matrix elements for processes with beams (-3, 3)
    return matrix_element[1] * 2; 
  }
  else if(id1 == 5 && id2 == -5)
  {
    // Add matrix elements for processes with beams (5, -5)
    return matrix_element[0] * 4; 
  }
  else if(id1 == -5 && id2 == 5)
  {
    // Add matrix elements for processes with beams (-5, 5)
    return matrix_element[1] * 4; 
  }
  else if(id1 == 3 && id2 == -3)
  {
    // Add matrix elements for processes with beams (3, -3)
    return matrix_element[0] * 2; 
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
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  ixxxxx(p[perm[2]], mME[2], hel[2], -1, w[2]); 
  oxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]); 
  oxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  oxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]); 
  ixxxxx(p[perm[6]], mME[6], hel[6], -1, w[6]); 
  ixxxxx(p[perm[7]], mME[7], hel[7], -1, w[7]); 
  FFV1P0_3(w[6], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[8]); 
  FFV1P0_3(w[7], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[9]); 
  FFV1P0_3(w[2], w[3], pars->GC_3, pars->ZERO, pars->ZERO, w[10]); 
  VVV1P0_1(w[8], w[9], pars->GC_10, pars->ZERO, pars->ZERO, w[11]); 
  FFV1_1(w[1], w[10], pars->GC_1, pars->ZERO, pars->ZERO, w[12]); 
  FFV1_2(w[0], w[10], pars->GC_1, pars->ZERO, pars->ZERO, w[13]); 
  FFV1_1(w[1], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[14]); 
  FFV1_2(w[0], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[15]); 
  FFV1_2(w[0], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[16]); 
  FFV1_1(w[1], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[17]); 
  FFV2_4_3(w[2], w[3], -pars->GC_51, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ,
      w[18]);
  FFV2_3_1(w[1], w[18], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[19]);
  FFV2_3_2(w[0], w[18], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[20]);
  FFV1_2(w[7], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[21]); 
  FFV1P0_3(w[21], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[22]); 
  FFV1_1(w[5], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[23]); 
  FFV1P0_3(w[7], w[23], pars->GC_11, pars->ZERO, pars->ZERO, w[24]); 
  FFV1_2(w[7], w[10], pars->GC_1, pars->ZERO, pars->ZERO, w[25]); 
  FFV1P0_3(w[0], w[14], pars->GC_11, pars->ZERO, pars->ZERO, w[26]); 
  FFV1_1(w[5], w[10], pars->GC_1, pars->ZERO, pars->ZERO, w[27]); 
  FFV1P0_3(w[16], w[1], pars->GC_11, pars->ZERO, pars->ZERO, w[28]); 
  FFV2_3_2(w[7], w[18], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[29]);
  FFV2_3_1(w[5], w[18], -pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[30]);
  FFV1P0_3(w[0], w[1], pars->GC_11, pars->ZERO, pars->ZERO, w[31]); 
  FFV1_1(w[5], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[32]); 
  FFV1_2(w[7], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[33]); 
  VVV1P0_1(w[8], w[31], pars->GC_10, pars->ZERO, pars->ZERO, w[34]); 
  FFV1_2(w[6], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[35]); 
  FFV1P0_3(w[35], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[36]); 
  FFV1_2(w[6], w[10], pars->GC_2, pars->ZERO, pars->ZERO, w[37]); 
  FFV1P0_3(w[37], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[38]); 
  FFV1_1(w[4], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[39]); 
  FFV1P0_3(w[6], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[40]); 
  FFV1_1(w[4], w[10], pars->GC_2, pars->ZERO, pars->ZERO, w[41]); 
  FFV1P0_3(w[6], w[41], pars->GC_11, pars->ZERO, pars->ZERO, w[42]); 
  FFV2_5_2(w[6], w[18], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[43]);
  FFV1P0_3(w[43], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[44]); 
  FFV2_5_1(w[4], w[18], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[45]);
  FFV1P0_3(w[6], w[45], pars->GC_11, pars->ZERO, pars->ZERO, w[46]); 
  FFV1_1(w[4], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[47]); 
  VVV1P0_1(w[9], w[31], pars->GC_10, pars->ZERO, pars->ZERO, w[48]); 
  FFV1_2(w[6], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[49]); 
  FFV1P0_3(w[49], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[50]); 
  FFV1P0_3(w[6], w[47], pars->GC_11, pars->ZERO, pars->ZERO, w[51]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[0], w[12], w[11], pars->GC_11, amp[0]); 
  FFV1_0(w[13], w[1], w[11], pars->GC_11, amp[1]); 
  FFV1_0(w[15], w[14], w[10], pars->GC_1, amp[2]); 
  FFV1_0(w[13], w[14], w[9], pars->GC_11, amp[3]); 
  FFV1_0(w[16], w[17], w[10], pars->GC_1, amp[4]); 
  FFV1_0(w[16], w[12], w[9], pars->GC_11, amp[5]); 
  FFV1_0(w[13], w[17], w[8], pars->GC_11, amp[6]); 
  FFV1_0(w[15], w[12], w[8], pars->GC_11, amp[7]); 
  FFV1_0(w[0], w[19], w[11], pars->GC_11, amp[8]); 
  FFV1_0(w[20], w[1], w[11], pars->GC_11, amp[9]); 
  FFV2_3_0(w[15], w[14], w[18], -pars->GC_51, pars->GC_58, amp[10]); 
  FFV1_0(w[20], w[14], w[9], pars->GC_11, amp[11]); 
  FFV2_3_0(w[16], w[17], w[18], -pars->GC_51, pars->GC_58, amp[12]); 
  FFV1_0(w[16], w[19], w[9], pars->GC_11, amp[13]); 
  FFV1_0(w[20], w[17], w[8], pars->GC_11, amp[14]); 
  FFV1_0(w[15], w[19], w[8], pars->GC_11, amp[15]); 
  FFV1_0(w[0], w[12], w[22], pars->GC_11, amp[16]); 
  FFV1_0(w[13], w[1], w[22], pars->GC_11, amp[17]); 
  FFV1_0(w[0], w[12], w[24], pars->GC_11, amp[18]); 
  FFV1_0(w[13], w[1], w[24], pars->GC_11, amp[19]); 
  FFV1_0(w[25], w[5], w[26], pars->GC_11, amp[20]); 
  FFV1_0(w[7], w[27], w[26], pars->GC_11, amp[21]); 
  FFV1_0(w[25], w[5], w[28], pars->GC_11, amp[22]); 
  FFV1_0(w[7], w[27], w[28], pars->GC_11, amp[23]); 
  FFV1_0(w[0], w[19], w[22], pars->GC_11, amp[24]); 
  FFV1_0(w[20], w[1], w[22], pars->GC_11, amp[25]); 
  FFV1_0(w[0], w[19], w[24], pars->GC_11, amp[26]); 
  FFV1_0(w[20], w[1], w[24], pars->GC_11, amp[27]); 
  FFV1_0(w[29], w[5], w[26], pars->GC_11, amp[28]); 
  FFV1_0(w[7], w[30], w[26], pars->GC_11, amp[29]); 
  FFV1_0(w[29], w[5], w[28], pars->GC_11, amp[30]); 
  FFV1_0(w[7], w[30], w[28], pars->GC_11, amp[31]); 
  FFV1_0(w[21], w[27], w[31], pars->GC_11, amp[32]); 
  FFV1_0(w[21], w[32], w[10], pars->GC_1, amp[33]); 
  FFV1_0(w[25], w[23], w[31], pars->GC_11, amp[34]); 
  FFV1_0(w[33], w[23], w[10], pars->GC_1, amp[35]); 
  FFV1_0(w[25], w[5], w[34], pars->GC_11, amp[36]); 
  FFV1_0(w[7], w[27], w[34], pars->GC_11, amp[37]); 
  FFV1_0(w[25], w[32], w[8], pars->GC_11, amp[38]); 
  FFV1_0(w[33], w[27], w[8], pars->GC_11, amp[39]); 
  FFV1_0(w[21], w[30], w[31], pars->GC_11, amp[40]); 
  FFV2_3_0(w[21], w[32], w[18], -pars->GC_51, pars->GC_58, amp[41]); 
  FFV1_0(w[29], w[23], w[31], pars->GC_11, amp[42]); 
  FFV2_3_0(w[33], w[23], w[18], -pars->GC_51, pars->GC_58, amp[43]); 
  FFV1_0(w[29], w[5], w[34], pars->GC_11, amp[44]); 
  FFV1_0(w[7], w[30], w[34], pars->GC_11, amp[45]); 
  FFV1_0(w[29], w[32], w[8], pars->GC_11, amp[46]); 
  FFV1_0(w[33], w[30], w[8], pars->GC_11, amp[47]); 
  FFV1_0(w[0], w[12], w[36], pars->GC_11, amp[48]); 
  FFV1_0(w[13], w[1], w[36], pars->GC_11, amp[49]); 
  FFV1_0(w[0], w[17], w[38], pars->GC_11, amp[50]); 
  FFV1_0(w[15], w[1], w[38], pars->GC_11, amp[51]); 
  FFV1_0(w[0], w[12], w[40], pars->GC_11, amp[52]); 
  FFV1_0(w[13], w[1], w[40], pars->GC_11, amp[53]); 
  FFV1_0(w[0], w[17], w[42], pars->GC_11, amp[54]); 
  FFV1_0(w[15], w[1], w[42], pars->GC_11, amp[55]); 
  FFV1_0(w[0], w[19], w[36], pars->GC_11, amp[56]); 
  FFV1_0(w[20], w[1], w[36], pars->GC_11, amp[57]); 
  FFV1_0(w[0], w[17], w[44], pars->GC_11, amp[58]); 
  FFV1_0(w[15], w[1], w[44], pars->GC_11, amp[59]); 
  FFV1_0(w[0], w[19], w[40], pars->GC_11, amp[60]); 
  FFV1_0(w[20], w[1], w[40], pars->GC_11, amp[61]); 
  FFV1_0(w[0], w[17], w[46], pars->GC_11, amp[62]); 
  FFV1_0(w[15], w[1], w[46], pars->GC_11, amp[63]); 
  FFV1_0(w[35], w[41], w[31], pars->GC_11, amp[64]); 
  FFV1_0(w[35], w[47], w[10], pars->GC_2, amp[65]); 
  FFV1_0(w[37], w[39], w[31], pars->GC_11, amp[66]); 
  FFV1_0(w[37], w[4], w[48], pars->GC_11, amp[67]); 
  FFV1_0(w[37], w[47], w[9], pars->GC_11, amp[68]); 
  FFV1_0(w[49], w[39], w[10], pars->GC_2, amp[69]); 
  FFV1_0(w[49], w[41], w[9], pars->GC_11, amp[70]); 
  FFV1_0(w[6], w[41], w[48], pars->GC_11, amp[71]); 
  FFV1_0(w[35], w[45], w[31], pars->GC_11, amp[72]); 
  FFV2_5_0(w[35], w[47], w[18], pars->GC_51, pars->GC_58, amp[73]); 
  FFV1_0(w[43], w[39], w[31], pars->GC_11, amp[74]); 
  FFV1_0(w[43], w[4], w[48], pars->GC_11, amp[75]); 
  FFV1_0(w[43], w[47], w[9], pars->GC_11, amp[76]); 
  FFV2_5_0(w[49], w[39], w[18], pars->GC_51, pars->GC_58, amp[77]); 
  FFV1_0(w[49], w[45], w[9], pars->GC_11, amp[78]); 
  FFV1_0(w[6], w[45], w[48], pars->GC_11, amp[79]); 
  FFV1_0(w[33], w[5], w[38], pars->GC_11, amp[80]); 
  FFV1_0(w[7], w[32], w[38], pars->GC_11, amp[81]); 
  FFV1_0(w[25], w[5], w[50], pars->GC_11, amp[82]); 
  FFV1_0(w[7], w[27], w[50], pars->GC_11, amp[83]); 
  FFV1_0(w[25], w[5], w[51], pars->GC_11, amp[84]); 
  FFV1_0(w[33], w[5], w[42], pars->GC_11, amp[85]); 
  FFV1_0(w[7], w[32], w[42], pars->GC_11, amp[86]); 
  FFV1_0(w[7], w[27], w[51], pars->GC_11, amp[87]); 
  FFV1_0(w[33], w[5], w[44], pars->GC_11, amp[88]); 
  FFV1_0(w[7], w[32], w[44], pars->GC_11, amp[89]); 
  FFV1_0(w[29], w[5], w[50], pars->GC_11, amp[90]); 
  FFV1_0(w[7], w[30], w[50], pars->GC_11, amp[91]); 
  FFV1_0(w[29], w[5], w[51], pars->GC_11, amp[92]); 
  FFV1_0(w[33], w[5], w[46], pars->GC_11, amp[93]); 
  FFV1_0(w[7], w[32], w[46], pars->GC_11, amp[94]); 
  FFV1_0(w[7], w[30], w[51], pars->GC_11, amp[95]); 

}
double CPPProcess::matrix_ssx_taptamuduxdx_no_h() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 96; 
  const int ncolor = 6; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1, 1, 1, 1, 1, 1}; 
  static const double cf[ncolor][ncolor] = {{27, 9, 9, 3, 3, 9}, {9, 27, 3, 9,
      9, 3}, {9, 3, 27, 9, 9, 3}, {3, 9, 9, 27, 3, 9}, {3, 9, 9, 3, 27, 9}, {9,
      3, 3, 9, 9, 27}};

  // Calculate color flows
  jamp[0] = +1./4. * (-1./9. * amp[2] - 1./9. * amp[3] - 1./9. * amp[4] - 1./9.
      * amp[5] - 1./9. * amp[6] - 1./9. * amp[7] - 1./9. * amp[10] - 1./9. *
      amp[11] - 1./9. * amp[12] - 1./9. * amp[13] - 1./9. * amp[14] - 1./9. *
      amp[15] - 1./9. * amp[16] - 1./9. * amp[17] - 1./9. * amp[18] - 1./9. *
      amp[19] - 1./9. * amp[20] - 1./9. * amp[21] - 1./9. * amp[22] - 1./9. *
      amp[23] - 1./9. * amp[24] - 1./9. * amp[25] - 1./9. * amp[26] - 1./9. *
      amp[27] - 1./9. * amp[28] - 1./9. * amp[29] - 1./9. * amp[30] - 1./9. *
      amp[31] - 1./9. * amp[32] - 1./9. * amp[33] - 1./9. * amp[34] - 1./9. *
      amp[35] - 1./9. * amp[38] - 1./9. * amp[39] - 1./9. * amp[40] - 1./9. *
      amp[41] - 1./9. * amp[42] - 1./9. * amp[43] - 1./9. * amp[46] - 1./9. *
      amp[47] - 1./9. * amp[48] - 1./9. * amp[49] - 1./9. * amp[50] - 1./9. *
      amp[51] - 1./9. * amp[52] - 1./9. * amp[53] - 1./9. * amp[54] - 1./9. *
      amp[55] - 1./9. * amp[56] - 1./9. * amp[57] - 1./9. * amp[58] - 1./9. *
      amp[59] - 1./9. * amp[60] - 1./9. * amp[61] - 1./9. * amp[62] - 1./9. *
      amp[63] - 1./9. * amp[64] - 1./9. * amp[65] - 1./9. * amp[66] - 1./9. *
      amp[68] - 1./9. * amp[69] - 1./9. * amp[70] - 1./9. * amp[72] - 1./9. *
      amp[73] - 1./9. * amp[74] - 1./9. * amp[76] - 1./9. * amp[77] - 1./9. *
      amp[78] - 1./9. * amp[80] - 1./9. * amp[81] - 1./9. * amp[82] - 1./9. *
      amp[83] - 1./9. * amp[84] - 1./9. * amp[85] - 1./9. * amp[86] - 1./9. *
      amp[87] - 1./9. * amp[88] - 1./9. * amp[89] - 1./9. * amp[90] - 1./9. *
      amp[91] - 1./9. * amp[92] - 1./9. * amp[93] - 1./9. * amp[94] - 1./9. *
      amp[95]);
  jamp[1] = +1./4. * (+1./3. * amp[16] + 1./3. * amp[17] + 1./3. * amp[18] +
      1./3. * amp[19] + 1./3. * amp[24] + 1./3. * amp[25] + 1./3. * amp[26] +
      1./3. * amp[27] + 1./3. * amp[32] + 1./3. * amp[33] + 1./3. * amp[34] +
      1./3. * amp[35] + 1./3. * amp[38] + 1./3. * amp[39] + 1./3. * amp[40] +
      1./3. * amp[41] + 1./3. * amp[42] + 1./3. * amp[43] + 1./3. * amp[46] +
      1./3. * amp[47] + 1./3. * amp[48] + 1./3. * amp[49] + 1./3. * amp[52] +
      1./3. * amp[53] + 1./3. * amp[56] + 1./3. * amp[57] + 1./3. * amp[60] +
      1./3. * amp[61] + 1./3. * amp[64] + 1./3. * amp[65] + 1./3. * amp[66] +
      1./3. * amp[68] + 1./3. * amp[69] + 1./3. * amp[70] + 1./3. * amp[72] +
      1./3. * amp[73] + 1./3. * amp[74] + 1./3. * amp[76] + 1./3. * amp[77] +
      1./3. * amp[78] + 1./3. * amp[80] + 1./3. * amp[81] + 1./3. * amp[82] +
      1./3. * amp[83] + 1./3. * amp[84] + 1./3. * amp[85] + 1./3. * amp[86] +
      1./3. * amp[87] + 1./3. * amp[88] + 1./3. * amp[89] + 1./3. * amp[90] +
      1./3. * amp[91] + 1./3. * amp[92] + 1./3. * amp[93] + 1./3. * amp[94] +
      1./3. * amp[95]);
  jamp[2] = +1./4. * (+1./3. * amp[2] + 1./3. * amp[3] + 1./3. * amp[4] + 1./3.
      * amp[5] + 1./3. * amp[6] + 1./3. * amp[7] + 1./3. * amp[10] + 1./3. *
      amp[11] + 1./3. * amp[12] + 1./3. * amp[13] + 1./3. * amp[14] + 1./3. *
      amp[15] + 1./3. * amp[20] + 1./3. * amp[21] + 1./3. * amp[22] + 1./3. *
      amp[23] + 1./3. * amp[28] + 1./3. * amp[29] + 1./3. * amp[30] + 1./3. *
      amp[31] + 1./3. * amp[48] + 1./3. * amp[49] + 1./3. * amp[50] + 1./3. *
      amp[51] + 1./3. * amp[52] + 1./3. * amp[53] + 1./3. * amp[54] + 1./3. *
      amp[55] + 1./3. * amp[56] + 1./3. * amp[57] + 1./3. * amp[58] + 1./3. *
      amp[59] + 1./3. * amp[60] + 1./3. * amp[61] + 1./3. * amp[62] + 1./3. *
      amp[63] + 1./3. * amp[64] + 1./3. * amp[65] + 1./3. * amp[66] + 1./3. *
      amp[68] + 1./3. * amp[69] + 1./3. * amp[70] + 1./3. * amp[72] + 1./3. *
      amp[73] + 1./3. * amp[74] + 1./3. * amp[76] + 1./3. * amp[77] + 1./3. *
      amp[78] + 1./3. * amp[82] + 1./3. * amp[83] + 1./3. * amp[84] + 1./3. *
      amp[87] + 1./3. * amp[90] + 1./3. * amp[91] + 1./3. * amp[92] + 1./3. *
      amp[95]);
  jamp[3] = +1./4. * (+std::complex<double> (0, 1) * amp[0] +
      std::complex<double> (0, 1) * amp[1] - amp[2] - amp[3] - amp[7] +
      std::complex<double> (0, 1) * amp[8] + std::complex<double> (0, 1) *
      amp[9] - amp[10] - amp[11] - amp[15] - amp[16] - amp[17] - amp[20] -
      amp[21] - amp[24] - amp[25] - amp[28] - amp[29] - amp[32] - amp[33] -
      std::complex<double> (0, 1) * amp[36] - std::complex<double> (0, 1) *
      amp[37] - amp[38] - amp[40] - amp[41] - std::complex<double> (0, 1) *
      amp[44] - std::complex<double> (0, 1) * amp[45] - amp[46] - amp[51] -
      amp[52] - amp[53] - amp[55] - amp[59] - amp[60] - amp[61] - amp[63] -
      amp[66] + std::complex<double> (0, 1) * amp[67] - amp[69] - amp[70] +
      std::complex<double> (0, 1) * amp[71] - amp[74] + std::complex<double>
      (0, 1) * amp[75] - amp[77] - amp[78] + std::complex<double> (0, 1) *
      amp[79] - amp[81] - amp[82] - amp[83] - amp[86] - amp[89] - amp[90] -
      amp[91] - amp[94]);
  jamp[4] = +1./4. * (-std::complex<double> (0, 1) * amp[0] -
      std::complex<double> (0, 1) * amp[1] - amp[4] - amp[5] - amp[6] -
      std::complex<double> (0, 1) * amp[8] - std::complex<double> (0, 1) *
      amp[9] - amp[12] - amp[13] - amp[14] - amp[18] - amp[19] - amp[22] -
      amp[23] - amp[26] - amp[27] - amp[30] - amp[31] - amp[34] - amp[35] +
      std::complex<double> (0, 1) * amp[36] + std::complex<double> (0, 1) *
      amp[37] - amp[39] - amp[42] - amp[43] + std::complex<double> (0, 1) *
      amp[44] + std::complex<double> (0, 1) * amp[45] - amp[47] - amp[48] -
      amp[49] - amp[50] - amp[54] - amp[56] - amp[57] - amp[58] - amp[62] -
      amp[64] - amp[65] - std::complex<double> (0, 1) * amp[67] - amp[68] -
      std::complex<double> (0, 1) * amp[71] - amp[72] - amp[73] -
      std::complex<double> (0, 1) * amp[75] - amp[76] - std::complex<double>
      (0, 1) * amp[79] - amp[80] - amp[84] - amp[85] - amp[87] - amp[88] -
      amp[92] - amp[93] - amp[95]);
  jamp[5] = +1./4. * (+1./3. * amp[2] + 1./3. * amp[3] + 1./3. * amp[4] + 1./3.
      * amp[5] + 1./3. * amp[6] + 1./3. * amp[7] + 1./3. * amp[10] + 1./3. *
      amp[11] + 1./3. * amp[12] + 1./3. * amp[13] + 1./3. * amp[14] + 1./3. *
      amp[15] + 1./3. * amp[16] + 1./3. * amp[17] + 1./3. * amp[18] + 1./3. *
      amp[19] + 1./3. * amp[20] + 1./3. * amp[21] + 1./3. * amp[22] + 1./3. *
      amp[23] + 1./3. * amp[24] + 1./3. * amp[25] + 1./3. * amp[26] + 1./3. *
      amp[27] + 1./3. * amp[28] + 1./3. * amp[29] + 1./3. * amp[30] + 1./3. *
      amp[31] + 1./3. * amp[32] + 1./3. * amp[33] + 1./3. * amp[34] + 1./3. *
      amp[35] + 1./3. * amp[38] + 1./3. * amp[39] + 1./3. * amp[40] + 1./3. *
      amp[41] + 1./3. * amp[42] + 1./3. * amp[43] + 1./3. * amp[46] + 1./3. *
      amp[47] + 1./3. * amp[50] + 1./3. * amp[51] + 1./3. * amp[54] + 1./3. *
      amp[55] + 1./3. * amp[58] + 1./3. * amp[59] + 1./3. * amp[62] + 1./3. *
      amp[63] + 1./3. * amp[80] + 1./3. * amp[81] + 1./3. * amp[85] + 1./3. *
      amp[86] + 1./3. * amp[88] + 1./3. * amp[89] + 1./3. * amp[93] + 1./3. *
      amp[94]);

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



