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
// Process: g u > ta+ ta- g g g u WEIGHTED<=8 / h
// Process: g c > ta+ ta- g g g c WEIGHTED<=8 / h

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
  jamp2[0] = new double[24]; 
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
  for(int i = 0; i < 24; i++ )
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
  const int denominators[nprocesses] = {576, 576}; 

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
        t[0] = matrix_gu_taptamgggu_no_h(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_gu_taptamgggu_no_h(); 
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
      t[0] = matrix_gu_taptamgggu_no_h(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_gu_taptamgggu_no_h(); 
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
  if(id1 == 4 && id2 == 21)
  {
    // Add matrix elements for processes with beams (4, 21)
    return matrix_element[1]; 
  }
  else if(id1 == 21 && id2 == 4)
  {
    // Add matrix elements for processes with beams (21, 4)
    return matrix_element[0]; 
  }
  else if(id1 == 21 && id2 == 2)
  {
    // Add matrix elements for processes with beams (21, 2)
    return matrix_element[0]; 
  }
  else if(id1 == 2 && id2 == 21)
  {
    // Add matrix elements for processes with beams (2, 21)
    return matrix_element[1]; 
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
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  ixxxxx(p[perm[1]], mME[1], hel[1], +1, w[1]); 
  ixxxxx(p[perm[2]], mME[2], hel[2], -1, w[2]); 
  oxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]); 
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  vxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]); 
  vxxxxx(p[perm[6]], mME[6], hel[6], +1, w[6]); 
  oxxxxx(p[perm[7]], mME[7], hel[7], +1, w[7]); 
  VVV1P0_1(w[0], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[8]); 
  FFV1P0_3(w[2], w[3], pars->GC_3, pars->ZERO, pars->ZERO, w[9]); 
  VVV1P0_1(w[8], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[10]); 
  FFV1_1(w[7], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[11]); 
  VVV1P0_1(w[10], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[12]); 
  FFV1_2(w[1], w[10], pars->GC_11, pars->ZERO, pars->ZERO, w[13]); 
  FFV1_2(w[1], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[14]); 
  FFV1_1(w[7], w[10], pars->GC_11, pars->ZERO, pars->ZERO, w[15]); 
  VVV1P0_1(w[8], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[16]); 
  VVV1P0_1(w[16], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[17]); 
  FFV1_2(w[1], w[16], pars->GC_11, pars->ZERO, pars->ZERO, w[18]); 
  FFV1_1(w[7], w[16], pars->GC_11, pars->ZERO, pars->ZERO, w[19]); 
  FFV1_1(w[7], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[20]); 
  FFV1_1(w[20], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[21]); 
  FFV1_1(w[20], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[22]); 
  FFV1_2(w[1], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[23]); 
  FFV1_2(w[23], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[24]); 
  FFV1_2(w[23], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[25]); 
  VVVV1P0_1(w[8], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[26]); 
  VVVV3P0_1(w[8], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[27]); 
  VVVV4P0_1(w[8], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[28]); 
  FFV2_4_3(w[2], w[3], -pars->GC_51, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ,
      w[29]);
  FFV2_5_1(w[7], w[29], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[30]);
  FFV2_5_2(w[1], w[29], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[31]);
  VVV1P0_1(w[5], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[32]); 
  VVV1P0_1(w[8], w[32], pars->GC_10, pars->ZERO, pars->ZERO, w[33]); 
  FFV1_2(w[1], w[32], pars->GC_11, pars->ZERO, pars->ZERO, w[34]); 
  FFV1_1(w[7], w[32], pars->GC_11, pars->ZERO, pars->ZERO, w[35]); 
  FFV1_1(w[7], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[36]); 
  FFV1_1(w[36], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[37]); 
  FFV1_1(w[36], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[38]); 
  FFV1_1(w[36], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[39]); 
  FFV2_5_1(w[36], w[29], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[40]);
  FFV1_2(w[1], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[41]); 
  FFV1_2(w[41], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[42]); 
  FFV1_2(w[1], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[43]); 
  FFV1_2(w[43], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[44]); 
  FFV1_2(w[43], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[45]); 
  FFV1_2(w[43], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[46]); 
  FFV2_5_2(w[43], w[29], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[47]);
  FFV1_1(w[7], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[48]); 
  FFV1_1(w[48], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[49]); 
  FFV1_1(w[48], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[50]); 
  FFV1_1(w[48], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[51]); 
  FFV2_5_1(w[48], w[29], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[52]);
  FFV1_2(w[41], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[53]); 
  FFV1_2(w[41], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[54]); 
  FFV2_5_2(w[41], w[29], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[55]);
  VVV1P0_1(w[0], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[56]); 
  VVV1P0_1(w[4], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[57]); 
  VVV1P0_1(w[56], w[57], pars->GC_10, pars->ZERO, pars->ZERO, w[58]); 
  FFV1_1(w[7], w[56], pars->GC_11, pars->ZERO, pars->ZERO, w[59]); 
  FFV1_2(w[1], w[57], pars->GC_11, pars->ZERO, pars->ZERO, w[60]); 
  FFV1_2(w[1], w[56], pars->GC_11, pars->ZERO, pars->ZERO, w[61]); 
  FFV1_1(w[7], w[57], pars->GC_11, pars->ZERO, pars->ZERO, w[62]); 
  FFV1_1(w[7], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[63]); 
  FFV1_1(w[63], w[56], pars->GC_11, pars->ZERO, pars->ZERO, w[64]); 
  VVV1P0_1(w[56], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[65]); 
  FFV1_1(w[63], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[66]); 
  FFV1_1(w[63], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[67]); 
  FFV2_5_1(w[63], w[29], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[68]);
  FFV1_2(w[41], w[56], pars->GC_11, pars->ZERO, pars->ZERO, w[69]); 
  FFV1_2(w[1], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[70]); 
  FFV1_2(w[70], w[56], pars->GC_11, pars->ZERO, pars->ZERO, w[71]); 
  FFV1_2(w[70], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[72]); 
  FFV1_2(w[70], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[73]); 
  FFV2_5_2(w[70], w[29], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[74]);
  FFV1_1(w[48], w[56], pars->GC_11, pars->ZERO, pars->ZERO, w[75]); 
  VVV1P0_1(w[56], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[76]); 
  VVV1P0_1(w[76], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[77]); 
  FFV1_2(w[1], w[76], pars->GC_11, pars->ZERO, pars->ZERO, w[78]); 
  FFV1_1(w[7], w[76], pars->GC_11, pars->ZERO, pars->ZERO, w[79]); 
  VVV1P0_1(w[65], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[80]); 
  FFV1_2(w[1], w[65], pars->GC_11, pars->ZERO, pars->ZERO, w[81]); 
  FFV1_1(w[7], w[65], pars->GC_11, pars->ZERO, pars->ZERO, w[82]); 
  FFV1_1(w[59], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[83]); 
  FFV1_1(w[59], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[84]); 
  FFV1_2(w[61], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[85]); 
  FFV1_2(w[61], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[86]); 
  VVVV1P0_1(w[56], w[4], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[87]); 
  VVVV3P0_1(w[56], w[4], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[88]); 
  VVVV4P0_1(w[56], w[4], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[89]); 
  FFV1_1(w[48], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[90]); 
  FFV1_2(w[41], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[91]); 
  VVV1P0_1(w[0], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[92]); 
  VVV1P0_1(w[4], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[93]); 
  VVV1P0_1(w[92], w[93], pars->GC_10, pars->ZERO, pars->ZERO, w[94]); 
  FFV1_1(w[7], w[92], pars->GC_11, pars->ZERO, pars->ZERO, w[95]); 
  FFV1_2(w[1], w[93], pars->GC_11, pars->ZERO, pars->ZERO, w[96]); 
  FFV1_2(w[1], w[92], pars->GC_11, pars->ZERO, pars->ZERO, w[97]); 
  FFV1_1(w[7], w[93], pars->GC_11, pars->ZERO, pars->ZERO, w[98]); 
  FFV1_1(w[63], w[92], pars->GC_11, pars->ZERO, pars->ZERO, w[99]); 
  VVV1P0_1(w[92], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[100]); 
  FFV1_1(w[63], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[101]); 
  FFV1_2(w[43], w[92], pars->GC_11, pars->ZERO, pars->ZERO, w[102]); 
  FFV1_2(w[70], w[92], pars->GC_11, pars->ZERO, pars->ZERO, w[103]); 
  FFV1_2(w[70], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[104]); 
  FFV1_1(w[36], w[92], pars->GC_11, pars->ZERO, pars->ZERO, w[105]); 
  VVV1P0_1(w[92], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[106]); 
  VVV1P0_1(w[106], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[107]); 
  FFV1_2(w[1], w[106], pars->GC_11, pars->ZERO, pars->ZERO, w[108]); 
  FFV1_1(w[7], w[106], pars->GC_11, pars->ZERO, pars->ZERO, w[109]); 
  VVV1P0_1(w[100], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[110]); 
  FFV1_2(w[1], w[100], pars->GC_11, pars->ZERO, pars->ZERO, w[111]); 
  FFV1_1(w[7], w[100], pars->GC_11, pars->ZERO, pars->ZERO, w[112]); 
  FFV1_1(w[95], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[113]); 
  FFV1_1(w[95], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[114]); 
  FFV1_2(w[97], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[115]); 
  FFV1_2(w[97], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[116]); 
  VVVV1P0_1(w[92], w[4], w[5], pars->GC_12, pars->ZERO, pars->ZERO, w[117]); 
  VVVV3P0_1(w[92], w[4], w[5], pars->GC_12, pars->ZERO, pars->ZERO, w[118]); 
  VVVV4P0_1(w[92], w[4], w[5], pars->GC_12, pars->ZERO, pars->ZERO, w[119]); 
  FFV1_1(w[36], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[120]); 
  FFV1_2(w[43], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[121]); 
  FFV1_1(w[7], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[122]); 
  FFV1_1(w[122], w[93], pars->GC_11, pars->ZERO, pars->ZERO, w[123]); 
  FFV1_1(w[122], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[124]); 
  VVV1P0_1(w[93], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[125]); 
  FFV1_1(w[122], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[126]); 
  FFV2_5_1(w[122], w[29], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[127]);
  FFV1_1(w[122], w[57], pars->GC_11, pars->ZERO, pars->ZERO, w[128]); 
  VVV1P0_1(w[57], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[129]); 
  FFV1_1(w[122], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[130]); 
  FFV1_1(w[122], w[32], pars->GC_11, pars->ZERO, pars->ZERO, w[131]); 
  FFV1_1(w[122], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[132]); 
  FFV1_1(w[132], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[133]); 
  FFV1_1(w[132], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[134]); 
  FFV1_1(w[130], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[135]); 
  FFV1_1(w[130], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[136]); 
  FFV1_1(w[126], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[137]); 
  FFV1_1(w[126], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[138]); 
  VVV1P0_1(w[4], w[32], pars->GC_10, pars->ZERO, pars->ZERO, w[139]); 
  VVVV1P0_1(w[4], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[140]); 
  VVVV3P0_1(w[4], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[141]); 
  VVVV4P0_1(w[4], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[142]); 
  FFV1_1(w[122], w[140], pars->GC_11, pars->ZERO, pars->ZERO, w[143]); 
  FFV1_1(w[122], w[141], pars->GC_11, pars->ZERO, pars->ZERO, w[144]); 
  FFV1_1(w[122], w[142], pars->GC_11, pars->ZERO, pars->ZERO, w[145]); 
  FFV1_2(w[1], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[146]); 
  FFV1_2(w[146], w[93], pars->GC_11, pars->ZERO, pars->ZERO, w[147]); 
  FFV1_2(w[146], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[148]); 
  FFV1_2(w[146], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[149]); 
  FFV2_5_2(w[146], w[29], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[150]);
  FFV1_2(w[146], w[57], pars->GC_11, pars->ZERO, pars->ZERO, w[151]); 
  FFV1_2(w[146], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[152]); 
  FFV1_2(w[146], w[32], pars->GC_11, pars->ZERO, pars->ZERO, w[153]); 
  FFV1_2(w[146], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[154]); 
  FFV1_2(w[154], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[155]); 
  FFV1_2(w[154], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[156]); 
  FFV1_2(w[152], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[157]); 
  FFV1_2(w[152], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[158]); 
  FFV1_2(w[149], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[159]); 
  FFV1_2(w[149], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[160]); 
  FFV1_2(w[146], w[140], pars->GC_11, pars->ZERO, pars->ZERO, w[161]); 
  FFV1_2(w[146], w[141], pars->GC_11, pars->ZERO, pars->ZERO, w[162]); 
  FFV1_2(w[146], w[142], pars->GC_11, pars->ZERO, pars->ZERO, w[163]); 
  VVV1P0_1(w[0], w[93], pars->GC_10, pars->ZERO, pars->ZERO, w[164]); 
  VVV1P0_1(w[164], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[165]); 
  FFV1_2(w[1], w[164], pars->GC_11, pars->ZERO, pars->ZERO, w[166]); 
  FFV1_1(w[7], w[164], pars->GC_11, pars->ZERO, pars->ZERO, w[167]); 
  VVV1P0_1(w[0], w[125], pars->GC_10, pars->ZERO, pars->ZERO, w[168]); 
  FFV1_1(w[11], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[169]); 
  FFV1_2(w[14], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[170]); 
  FFV1_1(w[98], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[171]); 
  FFV1_2(w[96], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[172]); 
  VVVV1P0_1(w[0], w[93], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[173]); 
  VVVV3P0_1(w[0], w[93], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[174]); 
  VVVV4P0_1(w[0], w[93], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[175]); 
  FFV1_1(w[30], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[176]); 
  FFV1_2(w[31], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[177]); 
  FFV1_1(w[48], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[178]); 
  FFV1_1(w[48], w[93], pars->GC_11, pars->ZERO, pars->ZERO, w[179]); 
  FFV1_2(w[41], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[180]); 
  FFV1_2(w[41], w[93], pars->GC_11, pars->ZERO, pars->ZERO, w[181]); 
  VVV1P0_1(w[0], w[57], pars->GC_10, pars->ZERO, pars->ZERO, w[182]); 
  VVV1P0_1(w[182], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[183]); 
  FFV1_2(w[1], w[182], pars->GC_11, pars->ZERO, pars->ZERO, w[184]); 
  FFV1_1(w[7], w[182], pars->GC_11, pars->ZERO, pars->ZERO, w[185]); 
  VVV1P0_1(w[0], w[129], pars->GC_10, pars->ZERO, pars->ZERO, w[186]); 
  FFV1_1(w[62], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[187]); 
  FFV1_2(w[60], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[188]); 
  VVVV1P0_1(w[0], w[57], w[5], pars->GC_12, pars->ZERO, pars->ZERO, w[189]); 
  VVVV3P0_1(w[0], w[57], w[5], pars->GC_12, pars->ZERO, pars->ZERO, w[190]); 
  VVVV4P0_1(w[0], w[57], w[5], pars->GC_12, pars->ZERO, pars->ZERO, w[191]); 
  FFV1_1(w[36], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[192]); 
  FFV1_1(w[36], w[57], pars->GC_11, pars->ZERO, pars->ZERO, w[193]); 
  FFV1_2(w[43], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[194]); 
  FFV1_2(w[43], w[57], pars->GC_11, pars->ZERO, pars->ZERO, w[195]); 
  FFV1_1(w[63], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[196]); 
  FFV1_1(w[196], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[197]); 
  FFV1_1(w[196], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[198]); 
  FFV1_1(w[101], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[199]); 
  FFV1_1(w[67], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[200]); 
  VVV1P0_1(w[0], w[32], pars->GC_10, pars->ZERO, pars->ZERO, w[201]); 
  FFV1_1(w[63], w[32], pars->GC_11, pars->ZERO, pars->ZERO, w[202]); 
  FFV1_2(w[70], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[203]); 
  FFV1_2(w[203], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[204]); 
  FFV1_2(w[203], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[205]); 
  FFV1_2(w[104], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[206]); 
  FFV1_2(w[73], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[207]); 
  FFV1_2(w[70], w[32], pars->GC_11, pars->ZERO, pars->ZERO, w[208]); 
  VVV1P0_1(w[201], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[209]); 
  FFV1_2(w[1], w[201], pars->GC_11, pars->ZERO, pars->ZERO, w[210]); 
  FFV1_1(w[7], w[201], pars->GC_11, pars->ZERO, pars->ZERO, w[211]); 
  VVV1P0_1(w[0], w[139], pars->GC_10, pars->ZERO, pars->ZERO, w[212]); 
  FFV1_2(w[34], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[213]); 
  FFV1_1(w[35], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[214]); 
  VVVV1P0_1(w[0], w[4], w[32], pars->GC_12, pars->ZERO, pars->ZERO, w[215]); 
  VVVV3P0_1(w[0], w[4], w[32], pars->GC_12, pars->ZERO, pars->ZERO, w[216]); 
  VVVV4P0_1(w[0], w[4], w[32], pars->GC_12, pars->ZERO, pars->ZERO, w[217]); 
  FFV1_1(w[192], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[218]); 
  FFV1_1(w[192], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[219]); 
  FFV1_1(w[120], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[220]); 
  FFV1_1(w[39], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[221]); 
  FFV1_2(w[194], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[222]); 
  FFV1_2(w[194], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[223]); 
  FFV1_2(w[121], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[224]); 
  FFV1_2(w[46], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[225]); 
  FFV1_1(w[178], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[226]); 
  FFV1_1(w[178], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[227]); 
  FFV1_1(w[90], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[228]); 
  FFV1_1(w[51], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[229]); 
  FFV1_2(w[180], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[230]); 
  FFV1_2(w[180], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[231]); 
  FFV1_2(w[91], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[232]); 
  FFV1_2(w[54], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[233]); 
  VVVV1P0_1(w[0], w[4], w[5], pars->GC_12, pars->ZERO, pars->ZERO, w[234]); 
  VVVV3P0_1(w[0], w[4], w[5], pars->GC_12, pars->ZERO, pars->ZERO, w[235]); 
  VVVV4P0_1(w[0], w[4], w[5], pars->GC_12, pars->ZERO, pars->ZERO, w[236]); 
  VVV1P0_1(w[234], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[237]); 
  VVV1P0_1(w[235], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[238]); 
  VVV1P0_1(w[236], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[239]); 
  FFV1_1(w[7], w[234], pars->GC_11, pars->ZERO, pars->ZERO, w[240]); 
  FFV1_1(w[7], w[235], pars->GC_11, pars->ZERO, pars->ZERO, w[241]); 
  FFV1_1(w[7], w[236], pars->GC_11, pars->ZERO, pars->ZERO, w[242]); 
  FFV1_2(w[1], w[234], pars->GC_11, pars->ZERO, pars->ZERO, w[243]); 
  FFV1_2(w[1], w[235], pars->GC_11, pars->ZERO, pars->ZERO, w[244]); 
  FFV1_2(w[1], w[236], pars->GC_11, pars->ZERO, pars->ZERO, w[245]); 
  FFV1_1(w[48], w[234], pars->GC_11, pars->ZERO, pars->ZERO, w[246]); 
  FFV1_1(w[48], w[235], pars->GC_11, pars->ZERO, pars->ZERO, w[247]); 
  FFV1_1(w[48], w[236], pars->GC_11, pars->ZERO, pars->ZERO, w[248]); 
  FFV1_2(w[41], w[234], pars->GC_11, pars->ZERO, pars->ZERO, w[249]); 
  FFV1_2(w[41], w[235], pars->GC_11, pars->ZERO, pars->ZERO, w[250]); 
  FFV1_2(w[41], w[236], pars->GC_11, pars->ZERO, pars->ZERO, w[251]); 
  VVVV1P0_1(w[0], w[4], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[252]); 
  VVVV3P0_1(w[0], w[4], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[253]); 
  VVVV4P0_1(w[0], w[4], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[254]); 
  VVV1P0_1(w[252], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[255]); 
  VVV1P0_1(w[253], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[256]); 
  VVV1P0_1(w[254], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[257]); 
  FFV1_1(w[7], w[252], pars->GC_11, pars->ZERO, pars->ZERO, w[258]); 
  FFV1_1(w[7], w[253], pars->GC_11, pars->ZERO, pars->ZERO, w[259]); 
  FFV1_1(w[7], w[254], pars->GC_11, pars->ZERO, pars->ZERO, w[260]); 
  FFV1_2(w[1], w[252], pars->GC_11, pars->ZERO, pars->ZERO, w[261]); 
  FFV1_2(w[1], w[253], pars->GC_11, pars->ZERO, pars->ZERO, w[262]); 
  FFV1_2(w[1], w[254], pars->GC_11, pars->ZERO, pars->ZERO, w[263]); 
  FFV1_1(w[36], w[252], pars->GC_11, pars->ZERO, pars->ZERO, w[264]); 
  FFV1_1(w[36], w[253], pars->GC_11, pars->ZERO, pars->ZERO, w[265]); 
  FFV1_1(w[36], w[254], pars->GC_11, pars->ZERO, pars->ZERO, w[266]); 
  FFV1_2(w[43], w[252], pars->GC_11, pars->ZERO, pars->ZERO, w[267]); 
  FFV1_2(w[43], w[253], pars->GC_11, pars->ZERO, pars->ZERO, w[268]); 
  FFV1_2(w[43], w[254], pars->GC_11, pars->ZERO, pars->ZERO, w[269]); 
  VVVV1P0_1(w[0], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[270]); 
  VVVV3P0_1(w[0], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[271]); 
  VVVV4P0_1(w[0], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[272]); 
  FFV1_1(w[63], w[270], pars->GC_11, pars->ZERO, pars->ZERO, w[273]); 
  FFV1_1(w[63], w[271], pars->GC_11, pars->ZERO, pars->ZERO, w[274]); 
  FFV1_1(w[63], w[272], pars->GC_11, pars->ZERO, pars->ZERO, w[275]); 
  FFV1_2(w[1], w[270], pars->GC_11, pars->ZERO, pars->ZERO, w[276]); 
  FFV1_2(w[1], w[271], pars->GC_11, pars->ZERO, pars->ZERO, w[277]); 
  FFV1_2(w[1], w[272], pars->GC_11, pars->ZERO, pars->ZERO, w[278]); 
  FFV1_2(w[70], w[270], pars->GC_11, pars->ZERO, pars->ZERO, w[279]); 
  FFV1_2(w[70], w[271], pars->GC_11, pars->ZERO, pars->ZERO, w[280]); 
  FFV1_2(w[70], w[272], pars->GC_11, pars->ZERO, pars->ZERO, w[281]); 
  FFV1_1(w[7], w[270], pars->GC_11, pars->ZERO, pars->ZERO, w[282]); 
  FFV1_1(w[7], w[271], pars->GC_11, pars->ZERO, pars->ZERO, w[283]); 
  FFV1_1(w[7], w[272], pars->GC_11, pars->ZERO, pars->ZERO, w[284]); 
  VVV1P0_1(w[270], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[285]); 
  VVV1P0_1(w[271], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[286]); 
  VVV1P0_1(w[272], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[287]); 
  VVV1P0_1(w[0], w[140], pars->GC_10, pars->ZERO, pars->ZERO, w[288]); 
  VVV1P0_1(w[0], w[141], pars->GC_10, pars->ZERO, pars->ZERO, w[289]); 
  VVV1P0_1(w[0], w[142], pars->GC_10, pars->ZERO, pars->ZERO, w[290]); 
  FFV1_1(w[7], w[140], pars->GC_11, pars->ZERO, pars->ZERO, w[291]); 
  FFV1_1(w[7], w[141], pars->GC_11, pars->ZERO, pars->ZERO, w[292]); 
  FFV1_1(w[7], w[142], pars->GC_11, pars->ZERO, pars->ZERO, w[293]); 
  FFV1_2(w[1], w[140], pars->GC_11, pars->ZERO, pars->ZERO, w[294]); 
  FFV1_2(w[1], w[141], pars->GC_11, pars->ZERO, pars->ZERO, w[295]); 
  FFV1_2(w[1], w[142], pars->GC_11, pars->ZERO, pars->ZERO, w[296]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[1], w[11], w[12], pars->GC_11, amp[0]); 
  FFV1_0(w[13], w[11], w[6], pars->GC_11, amp[1]); 
  FFV1_0(w[14], w[7], w[12], pars->GC_11, amp[2]); 
  FFV1_0(w[14], w[15], w[6], pars->GC_11, amp[3]); 
  FFV1_0(w[1], w[11], w[17], pars->GC_11, amp[4]); 
  FFV1_0(w[18], w[11], w[5], pars->GC_11, amp[5]); 
  FFV1_0(w[14], w[7], w[17], pars->GC_11, amp[6]); 
  FFV1_0(w[14], w[19], w[5], pars->GC_11, amp[7]); 
  FFV1_0(w[14], w[21], w[6], pars->GC_11, amp[8]); 
  FFV1_0(w[14], w[22], w[5], pars->GC_11, amp[9]); 
  FFV1_0(w[24], w[11], w[6], pars->GC_11, amp[10]); 
  FFV1_0(w[25], w[11], w[5], pars->GC_11, amp[11]); 
  FFV1_0(w[1], w[11], w[26], pars->GC_11, amp[12]); 
  FFV1_0(w[1], w[11], w[27], pars->GC_11, amp[13]); 
  FFV1_0(w[1], w[11], w[28], pars->GC_11, amp[14]); 
  FFV1_0(w[14], w[7], w[26], pars->GC_11, amp[15]); 
  FFV1_0(w[14], w[7], w[27], pars->GC_11, amp[16]); 
  FFV1_0(w[14], w[7], w[28], pars->GC_11, amp[17]); 
  FFV1_0(w[1], w[30], w[12], pars->GC_11, amp[18]); 
  FFV1_0(w[13], w[30], w[6], pars->GC_11, amp[19]); 
  FFV1_0(w[31], w[7], w[12], pars->GC_11, amp[20]); 
  FFV1_0(w[31], w[15], w[6], pars->GC_11, amp[21]); 
  FFV1_0(w[1], w[30], w[17], pars->GC_11, amp[22]); 
  FFV1_0(w[18], w[30], w[5], pars->GC_11, amp[23]); 
  FFV1_0(w[31], w[7], w[17], pars->GC_11, amp[24]); 
  FFV1_0(w[31], w[19], w[5], pars->GC_11, amp[25]); 
  FFV1_0(w[31], w[21], w[6], pars->GC_11, amp[26]); 
  FFV1_0(w[31], w[22], w[5], pars->GC_11, amp[27]); 
  FFV1_0(w[24], w[30], w[6], pars->GC_11, amp[28]); 
  FFV1_0(w[25], w[30], w[5], pars->GC_11, amp[29]); 
  FFV1_0(w[1], w[30], w[26], pars->GC_11, amp[30]); 
  FFV1_0(w[1], w[30], w[27], pars->GC_11, amp[31]); 
  FFV1_0(w[1], w[30], w[28], pars->GC_11, amp[32]); 
  FFV1_0(w[31], w[7], w[26], pars->GC_11, amp[33]); 
  FFV1_0(w[31], w[7], w[27], pars->GC_11, amp[34]); 
  FFV1_0(w[31], w[7], w[28], pars->GC_11, amp[35]); 
  FFV1_0(w[1], w[11], w[33], pars->GC_11, amp[36]); 
  FFV1_0(w[14], w[7], w[33], pars->GC_11, amp[37]); 
  FFV1_0(w[14], w[20], w[32], pars->GC_11, amp[38]); 
  FFV1_0(w[34], w[20], w[9], pars->GC_2, amp[39]); 
  FFV1_0(w[23], w[11], w[32], pars->GC_11, amp[40]); 
  FFV1_0(w[23], w[35], w[9], pars->GC_2, amp[41]); 
  FFV1_0(w[34], w[11], w[8], pars->GC_11, amp[42]); 
  FFV1_0(w[14], w[35], w[8], pars->GC_11, amp[43]); 
  FFV1_0(w[1], w[30], w[33], pars->GC_11, amp[44]); 
  FFV1_0(w[31], w[7], w[33], pars->GC_11, amp[45]); 
  FFV1_0(w[31], w[20], w[32], pars->GC_11, amp[46]); 
  FFV2_5_0(w[34], w[20], w[29], pars->GC_51, pars->GC_58, amp[47]); 
  FFV1_0(w[23], w[30], w[32], pars->GC_11, amp[48]); 
  FFV2_5_0(w[23], w[35], w[29], pars->GC_51, pars->GC_58, amp[49]); 
  FFV1_0(w[34], w[30], w[8], pars->GC_11, amp[50]); 
  FFV1_0(w[31], w[35], w[8], pars->GC_11, amp[51]); 
  FFV1_0(w[14], w[37], w[6], pars->GC_11, amp[52]); 
  FFV1_0(w[1], w[38], w[16], pars->GC_11, amp[53]); 
  FFV1_0(w[14], w[36], w[16], pars->GC_11, amp[54]); 
  FFV1_0(w[23], w[38], w[6], pars->GC_11, amp[55]); 
  FFV1_0(w[23], w[39], w[9], pars->GC_2, amp[56]); 
  FFV1_0(w[14], w[39], w[8], pars->GC_11, amp[57]); 
  FFV1_0(w[31], w[37], w[6], pars->GC_11, amp[58]); 
  FFV1_0(w[1], w[40], w[16], pars->GC_11, amp[59]); 
  FFV1_0(w[31], w[36], w[16], pars->GC_11, amp[60]); 
  FFV1_0(w[23], w[40], w[6], pars->GC_11, amp[61]); 
  FFV2_5_0(w[23], w[39], w[29], pars->GC_51, pars->GC_58, amp[62]); 
  FFV1_0(w[31], w[39], w[8], pars->GC_11, amp[63]); 
  FFV1_0(w[41], w[37], w[9], pars->GC_2, amp[64]); 
  FFV1_0(w[42], w[36], w[9], pars->GC_2, amp[65]); 
  FFV2_5_0(w[41], w[37], w[29], pars->GC_51, pars->GC_58, amp[66]); 
  FFV2_5_0(w[42], w[36], w[29], pars->GC_51, pars->GC_58, amp[67]); 
  FFV1_0(w[44], w[11], w[6], pars->GC_11, amp[68]); 
  FFV1_0(w[45], w[7], w[16], pars->GC_11, amp[69]); 
  FFV1_0(w[43], w[11], w[16], pars->GC_11, amp[70]); 
  FFV1_0(w[45], w[20], w[6], pars->GC_11, amp[71]); 
  FFV1_0(w[46], w[20], w[9], pars->GC_2, amp[72]); 
  FFV1_0(w[46], w[11], w[8], pars->GC_11, amp[73]); 
  FFV1_0(w[44], w[30], w[6], pars->GC_11, amp[74]); 
  FFV1_0(w[47], w[7], w[16], pars->GC_11, amp[75]); 
  FFV1_0(w[43], w[30], w[16], pars->GC_11, amp[76]); 
  FFV1_0(w[47], w[20], w[6], pars->GC_11, amp[77]); 
  FFV2_5_0(w[46], w[20], w[29], pars->GC_51, pars->GC_58, amp[78]); 
  FFV1_0(w[46], w[30], w[8], pars->GC_11, amp[79]); 
  FFV1_0(w[44], w[48], w[9], pars->GC_2, amp[80]); 
  FFV1_0(w[43], w[49], w[9], pars->GC_2, amp[81]); 
  FFV2_5_0(w[44], w[48], w[29], pars->GC_51, pars->GC_58, amp[82]); 
  FFV2_5_0(w[43], w[49], w[29], pars->GC_51, pars->GC_58, amp[83]); 
  FFV1_0(w[1], w[50], w[10], pars->GC_11, amp[84]); 
  FFV1_0(w[14], w[48], w[10], pars->GC_11, amp[85]); 
  FFV1_0(w[14], w[49], w[5], pars->GC_11, amp[86]); 
  FFV1_0(w[23], w[50], w[5], pars->GC_11, amp[87]); 
  FFV1_0(w[23], w[51], w[9], pars->GC_2, amp[88]); 
  FFV1_0(w[14], w[51], w[8], pars->GC_11, amp[89]); 
  FFV1_0(w[1], w[52], w[10], pars->GC_11, amp[90]); 
  FFV1_0(w[31], w[48], w[10], pars->GC_11, amp[91]); 
  FFV1_0(w[31], w[49], w[5], pars->GC_11, amp[92]); 
  FFV1_0(w[23], w[52], w[5], pars->GC_11, amp[93]); 
  FFV2_5_0(w[23], w[51], w[29], pars->GC_51, pars->GC_58, amp[94]); 
  FFV1_0(w[31], w[51], w[8], pars->GC_11, amp[95]); 
  FFV1_0(w[53], w[7], w[10], pars->GC_11, amp[96]); 
  FFV1_0(w[41], w[11], w[10], pars->GC_11, amp[97]); 
  FFV1_0(w[42], w[11], w[5], pars->GC_11, amp[98]); 
  FFV1_0(w[53], w[20], w[5], pars->GC_11, amp[99]); 
  FFV1_0(w[54], w[20], w[9], pars->GC_2, amp[100]); 
  FFV1_0(w[54], w[11], w[8], pars->GC_11, amp[101]); 
  FFV1_0(w[55], w[7], w[10], pars->GC_11, amp[102]); 
  FFV1_0(w[41], w[30], w[10], pars->GC_11, amp[103]); 
  FFV1_0(w[42], w[30], w[5], pars->GC_11, amp[104]); 
  FFV1_0(w[55], w[20], w[5], pars->GC_11, amp[105]); 
  FFV2_5_0(w[54], w[20], w[29], pars->GC_51, pars->GC_58, amp[106]); 
  FFV1_0(w[54], w[30], w[8], pars->GC_11, amp[107]); 
  FFV1_0(w[1], w[11], w[58], pars->GC_11, amp[108]); 
  FFV1_0(w[14], w[7], w[58], pars->GC_11, amp[109]); 
  FFV1_0(w[60], w[59], w[9], pars->GC_2, amp[110]); 
  FFV1_0(w[14], w[59], w[57], pars->GC_11, amp[111]); 
  FFV1_0(w[61], w[62], w[9], pars->GC_2, amp[112]); 
  FFV1_0(w[61], w[11], w[57], pars->GC_11, amp[113]); 
  FFV1_0(w[14], w[62], w[56], pars->GC_11, amp[114]); 
  FFV1_0(w[60], w[11], w[56], pars->GC_11, amp[115]); 
  FFV1_0(w[1], w[30], w[58], pars->GC_11, amp[116]); 
  FFV1_0(w[31], w[7], w[58], pars->GC_11, amp[117]); 
  FFV2_5_0(w[60], w[59], w[29], pars->GC_51, pars->GC_58, amp[118]); 
  FFV1_0(w[31], w[59], w[57], pars->GC_11, amp[119]); 
  FFV2_5_0(w[61], w[62], w[29], pars->GC_51, pars->GC_58, amp[120]); 
  FFV1_0(w[61], w[30], w[57], pars->GC_11, amp[121]); 
  FFV1_0(w[31], w[62], w[56], pars->GC_11, amp[122]); 
  FFV1_0(w[60], w[30], w[56], pars->GC_11, amp[123]); 
  FFV1_0(w[14], w[64], w[6], pars->GC_11, amp[124]); 
  FFV1_0(w[1], w[66], w[65], pars->GC_11, amp[125]); 
  FFV1_0(w[14], w[63], w[65], pars->GC_11, amp[126]); 
  FFV1_0(w[61], w[66], w[6], pars->GC_11, amp[127]); 
  FFV1_0(w[61], w[67], w[9], pars->GC_2, amp[128]); 
  FFV1_0(w[14], w[67], w[56], pars->GC_11, amp[129]); 
  FFV1_0(w[31], w[64], w[6], pars->GC_11, amp[130]); 
  FFV1_0(w[1], w[68], w[65], pars->GC_11, amp[131]); 
  FFV1_0(w[31], w[63], w[65], pars->GC_11, amp[132]); 
  FFV1_0(w[61], w[68], w[6], pars->GC_11, amp[133]); 
  FFV2_5_0(w[61], w[67], w[29], pars->GC_51, pars->GC_58, amp[134]); 
  FFV1_0(w[31], w[67], w[56], pars->GC_11, amp[135]); 
  FFV1_0(w[41], w[64], w[9], pars->GC_2, amp[136]); 
  FFV1_0(w[69], w[63], w[9], pars->GC_2, amp[137]); 
  FFV2_5_0(w[41], w[64], w[29], pars->GC_51, pars->GC_58, amp[138]); 
  FFV2_5_0(w[69], w[63], w[29], pars->GC_51, pars->GC_58, amp[139]); 
  FFV1_0(w[71], w[11], w[6], pars->GC_11, amp[140]); 
  FFV1_0(w[72], w[7], w[65], pars->GC_11, amp[141]); 
  FFV1_0(w[70], w[11], w[65], pars->GC_11, amp[142]); 
  FFV1_0(w[72], w[59], w[6], pars->GC_11, amp[143]); 
  FFV1_0(w[73], w[59], w[9], pars->GC_2, amp[144]); 
  FFV1_0(w[73], w[11], w[56], pars->GC_11, amp[145]); 
  FFV1_0(w[71], w[30], w[6], pars->GC_11, amp[146]); 
  FFV1_0(w[74], w[7], w[65], pars->GC_11, amp[147]); 
  FFV1_0(w[70], w[30], w[65], pars->GC_11, amp[148]); 
  FFV1_0(w[74], w[59], w[6], pars->GC_11, amp[149]); 
  FFV2_5_0(w[73], w[59], w[29], pars->GC_51, pars->GC_58, amp[150]); 
  FFV1_0(w[73], w[30], w[56], pars->GC_11, amp[151]); 
  FFV1_0(w[71], w[48], w[9], pars->GC_2, amp[152]); 
  FFV1_0(w[70], w[75], w[9], pars->GC_2, amp[153]); 
  FFV2_5_0(w[71], w[48], w[29], pars->GC_51, pars->GC_58, amp[154]); 
  FFV2_5_0(w[70], w[75], w[29], pars->GC_51, pars->GC_58, amp[155]); 
  FFV1_0(w[1], w[11], w[77], pars->GC_11, amp[156]); 
  FFV1_0(w[78], w[11], w[6], pars->GC_11, amp[157]); 
  FFV1_0(w[14], w[7], w[77], pars->GC_11, amp[158]); 
  FFV1_0(w[14], w[79], w[6], pars->GC_11, amp[159]); 
  FFV1_0(w[1], w[11], w[80], pars->GC_11, amp[160]); 
  FFV1_0(w[81], w[11], w[4], pars->GC_11, amp[161]); 
  FFV1_0(w[14], w[7], w[80], pars->GC_11, amp[162]); 
  FFV1_0(w[14], w[82], w[4], pars->GC_11, amp[163]); 
  FFV1_0(w[14], w[83], w[6], pars->GC_11, amp[164]); 
  FFV1_0(w[14], w[84], w[4], pars->GC_11, amp[165]); 
  FFV1_0(w[85], w[11], w[6], pars->GC_11, amp[166]); 
  FFV1_0(w[86], w[11], w[4], pars->GC_11, amp[167]); 
  FFV1_0(w[1], w[11], w[87], pars->GC_11, amp[168]); 
  FFV1_0(w[1], w[11], w[88], pars->GC_11, amp[169]); 
  FFV1_0(w[1], w[11], w[89], pars->GC_11, amp[170]); 
  FFV1_0(w[14], w[7], w[87], pars->GC_11, amp[171]); 
  FFV1_0(w[14], w[7], w[88], pars->GC_11, amp[172]); 
  FFV1_0(w[14], w[7], w[89], pars->GC_11, amp[173]); 
  FFV1_0(w[1], w[30], w[77], pars->GC_11, amp[174]); 
  FFV1_0(w[78], w[30], w[6], pars->GC_11, amp[175]); 
  FFV1_0(w[31], w[7], w[77], pars->GC_11, amp[176]); 
  FFV1_0(w[31], w[79], w[6], pars->GC_11, amp[177]); 
  FFV1_0(w[1], w[30], w[80], pars->GC_11, amp[178]); 
  FFV1_0(w[81], w[30], w[4], pars->GC_11, amp[179]); 
  FFV1_0(w[31], w[7], w[80], pars->GC_11, amp[180]); 
  FFV1_0(w[31], w[82], w[4], pars->GC_11, amp[181]); 
  FFV1_0(w[31], w[83], w[6], pars->GC_11, amp[182]); 
  FFV1_0(w[31], w[84], w[4], pars->GC_11, amp[183]); 
  FFV1_0(w[85], w[30], w[6], pars->GC_11, amp[184]); 
  FFV1_0(w[86], w[30], w[4], pars->GC_11, amp[185]); 
  FFV1_0(w[1], w[30], w[87], pars->GC_11, amp[186]); 
  FFV1_0(w[1], w[30], w[88], pars->GC_11, amp[187]); 
  FFV1_0(w[1], w[30], w[89], pars->GC_11, amp[188]); 
  FFV1_0(w[31], w[7], w[87], pars->GC_11, amp[189]); 
  FFV1_0(w[31], w[7], w[88], pars->GC_11, amp[190]); 
  FFV1_0(w[31], w[7], w[89], pars->GC_11, amp[191]); 
  FFV1_0(w[1], w[50], w[76], pars->GC_11, amp[192]); 
  FFV1_0(w[14], w[48], w[76], pars->GC_11, amp[193]); 
  FFV1_0(w[14], w[75], w[4], pars->GC_11, amp[194]); 
  FFV1_0(w[61], w[90], w[9], pars->GC_2, amp[195]); 
  FFV1_0(w[61], w[50], w[4], pars->GC_11, amp[196]); 
  FFV1_0(w[14], w[90], w[56], pars->GC_11, amp[197]); 
  FFV1_0(w[1], w[52], w[76], pars->GC_11, amp[198]); 
  FFV1_0(w[31], w[48], w[76], pars->GC_11, amp[199]); 
  FFV1_0(w[31], w[75], w[4], pars->GC_11, amp[200]); 
  FFV2_5_0(w[61], w[90], w[29], pars->GC_51, pars->GC_58, amp[201]); 
  FFV1_0(w[61], w[52], w[4], pars->GC_11, amp[202]); 
  FFV1_0(w[31], w[90], w[56], pars->GC_11, amp[203]); 
  FFV1_0(w[53], w[7], w[76], pars->GC_11, amp[204]); 
  FFV1_0(w[41], w[11], w[76], pars->GC_11, amp[205]); 
  FFV1_0(w[69], w[11], w[4], pars->GC_11, amp[206]); 
  FFV1_0(w[91], w[59], w[9], pars->GC_2, amp[207]); 
  FFV1_0(w[53], w[59], w[4], pars->GC_11, amp[208]); 
  FFV1_0(w[91], w[11], w[56], pars->GC_11, amp[209]); 
  FFV1_0(w[55], w[7], w[76], pars->GC_11, amp[210]); 
  FFV1_0(w[41], w[30], w[76], pars->GC_11, amp[211]); 
  FFV1_0(w[69], w[30], w[4], pars->GC_11, amp[212]); 
  FFV2_5_0(w[91], w[59], w[29], pars->GC_51, pars->GC_58, amp[213]); 
  FFV1_0(w[55], w[59], w[4], pars->GC_11, amp[214]); 
  FFV1_0(w[91], w[30], w[56], pars->GC_11, amp[215]); 
  FFV1_0(w[1], w[11], w[94], pars->GC_11, amp[216]); 
  FFV1_0(w[14], w[7], w[94], pars->GC_11, amp[217]); 
  FFV1_0(w[96], w[95], w[9], pars->GC_2, amp[218]); 
  FFV1_0(w[14], w[95], w[93], pars->GC_11, amp[219]); 
  FFV1_0(w[97], w[98], w[9], pars->GC_2, amp[220]); 
  FFV1_0(w[97], w[11], w[93], pars->GC_11, amp[221]); 
  FFV1_0(w[14], w[98], w[92], pars->GC_11, amp[222]); 
  FFV1_0(w[96], w[11], w[92], pars->GC_11, amp[223]); 
  FFV1_0(w[1], w[30], w[94], pars->GC_11, amp[224]); 
  FFV1_0(w[31], w[7], w[94], pars->GC_11, amp[225]); 
  FFV2_5_0(w[96], w[95], w[29], pars->GC_51, pars->GC_58, amp[226]); 
  FFV1_0(w[31], w[95], w[93], pars->GC_11, amp[227]); 
  FFV2_5_0(w[97], w[98], w[29], pars->GC_51, pars->GC_58, amp[228]); 
  FFV1_0(w[97], w[30], w[93], pars->GC_11, amp[229]); 
  FFV1_0(w[31], w[98], w[92], pars->GC_11, amp[230]); 
  FFV1_0(w[96], w[30], w[92], pars->GC_11, amp[231]); 
  FFV1_0(w[14], w[99], w[5], pars->GC_11, amp[232]); 
  FFV1_0(w[1], w[66], w[100], pars->GC_11, amp[233]); 
  FFV1_0(w[14], w[63], w[100], pars->GC_11, amp[234]); 
  FFV1_0(w[97], w[66], w[5], pars->GC_11, amp[235]); 
  FFV1_0(w[97], w[101], w[9], pars->GC_2, amp[236]); 
  FFV1_0(w[14], w[101], w[92], pars->GC_11, amp[237]); 
  FFV1_0(w[31], w[99], w[5], pars->GC_11, amp[238]); 
  FFV1_0(w[1], w[68], w[100], pars->GC_11, amp[239]); 
  FFV1_0(w[31], w[63], w[100], pars->GC_11, amp[240]); 
  FFV1_0(w[97], w[68], w[5], pars->GC_11, amp[241]); 
  FFV2_5_0(w[97], w[101], w[29], pars->GC_51, pars->GC_58, amp[242]); 
  FFV1_0(w[31], w[101], w[92], pars->GC_11, amp[243]); 
  FFV1_0(w[43], w[99], w[9], pars->GC_2, amp[244]); 
  FFV1_0(w[102], w[63], w[9], pars->GC_2, amp[245]); 
  FFV2_5_0(w[43], w[99], w[29], pars->GC_51, pars->GC_58, amp[246]); 
  FFV2_5_0(w[102], w[63], w[29], pars->GC_51, pars->GC_58, amp[247]); 
  FFV1_0(w[103], w[11], w[5], pars->GC_11, amp[248]); 
  FFV1_0(w[72], w[7], w[100], pars->GC_11, amp[249]); 
  FFV1_0(w[70], w[11], w[100], pars->GC_11, amp[250]); 
  FFV1_0(w[72], w[95], w[5], pars->GC_11, amp[251]); 
  FFV1_0(w[104], w[95], w[9], pars->GC_2, amp[252]); 
  FFV1_0(w[104], w[11], w[92], pars->GC_11, amp[253]); 
  FFV1_0(w[103], w[30], w[5], pars->GC_11, amp[254]); 
  FFV1_0(w[74], w[7], w[100], pars->GC_11, amp[255]); 
  FFV1_0(w[70], w[30], w[100], pars->GC_11, amp[256]); 
  FFV1_0(w[74], w[95], w[5], pars->GC_11, amp[257]); 
  FFV2_5_0(w[104], w[95], w[29], pars->GC_51, pars->GC_58, amp[258]); 
  FFV1_0(w[104], w[30], w[92], pars->GC_11, amp[259]); 
  FFV1_0(w[103], w[36], w[9], pars->GC_2, amp[260]); 
  FFV1_0(w[70], w[105], w[9], pars->GC_2, amp[261]); 
  FFV2_5_0(w[103], w[36], w[29], pars->GC_51, pars->GC_58, amp[262]); 
  FFV2_5_0(w[70], w[105], w[29], pars->GC_51, pars->GC_58, amp[263]); 
  FFV1_0(w[1], w[11], w[107], pars->GC_11, amp[264]); 
  FFV1_0(w[108], w[11], w[5], pars->GC_11, amp[265]); 
  FFV1_0(w[14], w[7], w[107], pars->GC_11, amp[266]); 
  FFV1_0(w[14], w[109], w[5], pars->GC_11, amp[267]); 
  FFV1_0(w[1], w[11], w[110], pars->GC_11, amp[268]); 
  FFV1_0(w[111], w[11], w[4], pars->GC_11, amp[269]); 
  FFV1_0(w[14], w[7], w[110], pars->GC_11, amp[270]); 
  FFV1_0(w[14], w[112], w[4], pars->GC_11, amp[271]); 
  FFV1_0(w[14], w[113], w[5], pars->GC_11, amp[272]); 
  FFV1_0(w[14], w[114], w[4], pars->GC_11, amp[273]); 
  FFV1_0(w[115], w[11], w[5], pars->GC_11, amp[274]); 
  FFV1_0(w[116], w[11], w[4], pars->GC_11, amp[275]); 
  FFV1_0(w[1], w[11], w[117], pars->GC_11, amp[276]); 
  FFV1_0(w[1], w[11], w[118], pars->GC_11, amp[277]); 
  FFV1_0(w[1], w[11], w[119], pars->GC_11, amp[278]); 
  FFV1_0(w[14], w[7], w[117], pars->GC_11, amp[279]); 
  FFV1_0(w[14], w[7], w[118], pars->GC_11, amp[280]); 
  FFV1_0(w[14], w[7], w[119], pars->GC_11, amp[281]); 
  FFV1_0(w[1], w[30], w[107], pars->GC_11, amp[282]); 
  FFV1_0(w[108], w[30], w[5], pars->GC_11, amp[283]); 
  FFV1_0(w[31], w[7], w[107], pars->GC_11, amp[284]); 
  FFV1_0(w[31], w[109], w[5], pars->GC_11, amp[285]); 
  FFV1_0(w[1], w[30], w[110], pars->GC_11, amp[286]); 
  FFV1_0(w[111], w[30], w[4], pars->GC_11, amp[287]); 
  FFV1_0(w[31], w[7], w[110], pars->GC_11, amp[288]); 
  FFV1_0(w[31], w[112], w[4], pars->GC_11, amp[289]); 
  FFV1_0(w[31], w[113], w[5], pars->GC_11, amp[290]); 
  FFV1_0(w[31], w[114], w[4], pars->GC_11, amp[291]); 
  FFV1_0(w[115], w[30], w[5], pars->GC_11, amp[292]); 
  FFV1_0(w[116], w[30], w[4], pars->GC_11, amp[293]); 
  FFV1_0(w[1], w[30], w[117], pars->GC_11, amp[294]); 
  FFV1_0(w[1], w[30], w[118], pars->GC_11, amp[295]); 
  FFV1_0(w[1], w[30], w[119], pars->GC_11, amp[296]); 
  FFV1_0(w[31], w[7], w[117], pars->GC_11, amp[297]); 
  FFV1_0(w[31], w[7], w[118], pars->GC_11, amp[298]); 
  FFV1_0(w[31], w[7], w[119], pars->GC_11, amp[299]); 
  FFV1_0(w[1], w[38], w[106], pars->GC_11, amp[300]); 
  FFV1_0(w[14], w[36], w[106], pars->GC_11, amp[301]); 
  FFV1_0(w[14], w[105], w[4], pars->GC_11, amp[302]); 
  FFV1_0(w[97], w[120], w[9], pars->GC_2, amp[303]); 
  FFV1_0(w[97], w[38], w[4], pars->GC_11, amp[304]); 
  FFV1_0(w[14], w[120], w[92], pars->GC_11, amp[305]); 
  FFV1_0(w[1], w[40], w[106], pars->GC_11, amp[306]); 
  FFV1_0(w[31], w[36], w[106], pars->GC_11, amp[307]); 
  FFV1_0(w[31], w[105], w[4], pars->GC_11, amp[308]); 
  FFV2_5_0(w[97], w[120], w[29], pars->GC_51, pars->GC_58, amp[309]); 
  FFV1_0(w[97], w[40], w[4], pars->GC_11, amp[310]); 
  FFV1_0(w[31], w[120], w[92], pars->GC_11, amp[311]); 
  FFV1_0(w[45], w[7], w[106], pars->GC_11, amp[312]); 
  FFV1_0(w[43], w[11], w[106], pars->GC_11, amp[313]); 
  FFV1_0(w[102], w[11], w[4], pars->GC_11, amp[314]); 
  FFV1_0(w[121], w[95], w[9], pars->GC_2, amp[315]); 
  FFV1_0(w[45], w[95], w[4], pars->GC_11, amp[316]); 
  FFV1_0(w[121], w[11], w[92], pars->GC_11, amp[317]); 
  FFV1_0(w[47], w[7], w[106], pars->GC_11, amp[318]); 
  FFV1_0(w[43], w[30], w[106], pars->GC_11, amp[319]); 
  FFV1_0(w[102], w[30], w[4], pars->GC_11, amp[320]); 
  FFV2_5_0(w[121], w[95], w[29], pars->GC_51, pars->GC_58, amp[321]); 
  FFV1_0(w[47], w[95], w[4], pars->GC_11, amp[322]); 
  FFV1_0(w[121], w[30], w[92], pars->GC_11, amp[323]); 
  FFV1_0(w[14], w[123], w[6], pars->GC_11, amp[324]); 
  FFV1_0(w[1], w[124], w[125], pars->GC_11, amp[325]); 
  FFV1_0(w[96], w[124], w[6], pars->GC_11, amp[326]); 
  FFV1_0(w[96], w[126], w[9], pars->GC_2, amp[327]); 
  FFV1_0(w[14], w[126], w[93], pars->GC_11, amp[328]); 
  FFV1_0(w[14], w[122], w[125], pars->GC_11, amp[329]); 
  FFV1_0(w[31], w[123], w[6], pars->GC_11, amp[330]); 
  FFV1_0(w[1], w[127], w[125], pars->GC_11, amp[331]); 
  FFV1_0(w[96], w[127], w[6], pars->GC_11, amp[332]); 
  FFV2_5_0(w[96], w[126], w[29], pars->GC_51, pars->GC_58, amp[333]); 
  FFV1_0(w[31], w[126], w[93], pars->GC_11, amp[334]); 
  FFV1_0(w[31], w[122], w[125], pars->GC_11, amp[335]); 
  FFV1_0(w[41], w[123], w[9], pars->GC_2, amp[336]); 
  FFV1_0(w[41], w[124], w[93], pars->GC_11, amp[337]); 
  FFV2_5_0(w[41], w[123], w[29], pars->GC_51, pars->GC_58, amp[338]); 
  FFV1_0(w[41], w[127], w[93], pars->GC_11, amp[339]); 
  FFV1_0(w[14], w[128], w[5], pars->GC_11, amp[340]); 
  FFV1_0(w[1], w[124], w[129], pars->GC_11, amp[341]); 
  FFV1_0(w[60], w[124], w[5], pars->GC_11, amp[342]); 
  FFV1_0(w[60], w[130], w[9], pars->GC_2, amp[343]); 
  FFV1_0(w[14], w[130], w[57], pars->GC_11, amp[344]); 
  FFV1_0(w[14], w[122], w[129], pars->GC_11, amp[345]); 
  FFV1_0(w[31], w[128], w[5], pars->GC_11, amp[346]); 
  FFV1_0(w[1], w[127], w[129], pars->GC_11, amp[347]); 
  FFV1_0(w[60], w[127], w[5], pars->GC_11, amp[348]); 
  FFV2_5_0(w[60], w[130], w[29], pars->GC_51, pars->GC_58, amp[349]); 
  FFV1_0(w[31], w[130], w[57], pars->GC_11, amp[350]); 
  FFV1_0(w[31], w[122], w[129], pars->GC_11, amp[351]); 
  FFV1_0(w[43], w[128], w[9], pars->GC_2, amp[352]); 
  FFV1_0(w[43], w[124], w[57], pars->GC_11, amp[353]); 
  FFV2_5_0(w[43], w[128], w[29], pars->GC_51, pars->GC_58, amp[354]); 
  FFV1_0(w[43], w[127], w[57], pars->GC_11, amp[355]); 
  FFV1_0(w[104], w[124], w[6], pars->GC_11, amp[356]); 
  FFV1_0(w[73], w[124], w[5], pars->GC_11, amp[357]); 
  FFV1_0(w[72], w[130], w[6], pars->GC_11, amp[358]); 
  FFV1_0(w[73], w[130], w[9], pars->GC_2, amp[359]); 
  FFV1_0(w[72], w[126], w[5], pars->GC_11, amp[360]); 
  FFV1_0(w[104], w[126], w[9], pars->GC_2, amp[361]); 
  FFV1_0(w[104], w[127], w[6], pars->GC_11, amp[362]); 
  FFV1_0(w[73], w[127], w[5], pars->GC_11, amp[363]); 
  FFV1_0(w[74], w[130], w[6], pars->GC_11, amp[364]); 
  FFV2_5_0(w[73], w[130], w[29], pars->GC_51, pars->GC_58, amp[365]); 
  FFV1_0(w[74], w[126], w[5], pars->GC_11, amp[366]); 
  FFV2_5_0(w[104], w[126], w[29], pars->GC_51, pars->GC_58, amp[367]); 
  FFV1_0(w[70], w[124], w[32], pars->GC_11, amp[368]); 
  FFV1_0(w[70], w[131], w[9], pars->GC_2, amp[369]); 
  FFV1_0(w[70], w[127], w[32], pars->GC_11, amp[370]); 
  FFV2_5_0(w[70], w[131], w[29], pars->GC_51, pars->GC_58, amp[371]); 
  FFV1_0(w[14], w[133], w[6], pars->GC_11, amp[372]); 
  FFV1_0(w[14], w[134], w[5], pars->GC_11, amp[373]); 
  FFV1_0(w[14], w[135], w[6], pars->GC_11, amp[374]); 
  FFV1_0(w[14], w[136], w[4], pars->GC_11, amp[375]); 
  FFV1_0(w[14], w[137], w[5], pars->GC_11, amp[376]); 
  FFV1_0(w[14], w[138], w[4], pars->GC_11, amp[377]); 
  FFV1_0(w[31], w[133], w[6], pars->GC_11, amp[378]); 
  FFV1_0(w[31], w[134], w[5], pars->GC_11, amp[379]); 
  FFV1_0(w[31], w[135], w[6], pars->GC_11, amp[380]); 
  FFV1_0(w[31], w[136], w[4], pars->GC_11, amp[381]); 
  FFV1_0(w[31], w[137], w[5], pars->GC_11, amp[382]); 
  FFV1_0(w[31], w[138], w[4], pars->GC_11, amp[383]); 
  FFV1_0(w[14], w[132], w[32], pars->GC_11, amp[384]); 
  FFV1_0(w[34], w[132], w[9], pars->GC_2, amp[385]); 
  FFV1_0(w[1], w[124], w[139], pars->GC_11, amp[386]); 
  FFV1_0(w[34], w[124], w[4], pars->GC_11, amp[387]); 
  FFV1_0(w[14], w[131], w[4], pars->GC_11, amp[388]); 
  FFV1_0(w[14], w[122], w[139], pars->GC_11, amp[389]); 
  FFV1_0(w[31], w[132], w[32], pars->GC_11, amp[390]); 
  FFV2_5_0(w[34], w[132], w[29], pars->GC_51, pars->GC_58, amp[391]); 
  FFV1_0(w[1], w[127], w[139], pars->GC_11, amp[392]); 
  FFV1_0(w[34], w[127], w[4], pars->GC_11, amp[393]); 
  FFV1_0(w[31], w[131], w[4], pars->GC_11, amp[394]); 
  FFV1_0(w[31], w[122], w[139], pars->GC_11, amp[395]); 
  FFV1_0(w[45], w[132], w[6], pars->GC_11, amp[396]); 
  FFV1_0(w[46], w[132], w[9], pars->GC_2, amp[397]); 
  FFV1_0(w[121], w[124], w[6], pars->GC_11, amp[398]); 
  FFV1_0(w[46], w[124], w[4], pars->GC_11, amp[399]); 
  FFV1_0(w[121], w[126], w[9], pars->GC_2, amp[400]); 
  FFV1_0(w[45], w[126], w[4], pars->GC_11, amp[401]); 
  FFV1_0(w[47], w[132], w[6], pars->GC_11, amp[402]); 
  FFV2_5_0(w[46], w[132], w[29], pars->GC_51, pars->GC_58, amp[403]); 
  FFV1_0(w[121], w[127], w[6], pars->GC_11, amp[404]); 
  FFV1_0(w[46], w[127], w[4], pars->GC_11, amp[405]); 
  FFV2_5_0(w[121], w[126], w[29], pars->GC_51, pars->GC_58, amp[406]); 
  FFV1_0(w[47], w[126], w[4], pars->GC_11, amp[407]); 
  FFV1_0(w[53], w[132], w[5], pars->GC_11, amp[408]); 
  FFV1_0(w[54], w[132], w[9], pars->GC_2, amp[409]); 
  FFV1_0(w[91], w[124], w[5], pars->GC_11, amp[410]); 
  FFV1_0(w[54], w[124], w[4], pars->GC_11, amp[411]); 
  FFV1_0(w[91], w[130], w[9], pars->GC_2, amp[412]); 
  FFV1_0(w[53], w[130], w[4], pars->GC_11, amp[413]); 
  FFV1_0(w[55], w[132], w[5], pars->GC_11, amp[414]); 
  FFV2_5_0(w[54], w[132], w[29], pars->GC_51, pars->GC_58, amp[415]); 
  FFV1_0(w[91], w[127], w[5], pars->GC_11, amp[416]); 
  FFV1_0(w[54], w[127], w[4], pars->GC_11, amp[417]); 
  FFV2_5_0(w[91], w[130], w[29], pars->GC_51, pars->GC_58, amp[418]); 
  FFV1_0(w[55], w[130], w[4], pars->GC_11, amp[419]); 
  FFV1_0(w[1], w[143], w[9], pars->GC_2, amp[420]); 
  FFV1_0(w[1], w[144], w[9], pars->GC_2, amp[421]); 
  FFV1_0(w[1], w[145], w[9], pars->GC_2, amp[422]); 
  FFV1_0(w[1], w[124], w[140], pars->GC_11, amp[423]); 
  FFV1_0(w[1], w[124], w[141], pars->GC_11, amp[424]); 
  FFV1_0(w[1], w[124], w[142], pars->GC_11, amp[425]); 
  FFV2_5_0(w[1], w[143], w[29], pars->GC_51, pars->GC_58, amp[426]); 
  FFV2_5_0(w[1], w[144], w[29], pars->GC_51, pars->GC_58, amp[427]); 
  FFV2_5_0(w[1], w[145], w[29], pars->GC_51, pars->GC_58, amp[428]); 
  FFV1_0(w[1], w[127], w[140], pars->GC_11, amp[429]); 
  FFV1_0(w[1], w[127], w[141], pars->GC_11, amp[430]); 
  FFV1_0(w[1], w[127], w[142], pars->GC_11, amp[431]); 
  FFV1_0(w[147], w[11], w[6], pars->GC_11, amp[432]); 
  FFV1_0(w[148], w[7], w[125], pars->GC_11, amp[433]); 
  FFV1_0(w[148], w[98], w[6], pars->GC_11, amp[434]); 
  FFV1_0(w[149], w[98], w[9], pars->GC_2, amp[435]); 
  FFV1_0(w[149], w[11], w[93], pars->GC_11, amp[436]); 
  FFV1_0(w[146], w[11], w[125], pars->GC_11, amp[437]); 
  FFV1_0(w[147], w[30], w[6], pars->GC_11, amp[438]); 
  FFV1_0(w[150], w[7], w[125], pars->GC_11, amp[439]); 
  FFV1_0(w[150], w[98], w[6], pars->GC_11, amp[440]); 
  FFV2_5_0(w[149], w[98], w[29], pars->GC_51, pars->GC_58, amp[441]); 
  FFV1_0(w[149], w[30], w[93], pars->GC_11, amp[442]); 
  FFV1_0(w[146], w[30], w[125], pars->GC_11, amp[443]); 
  FFV1_0(w[147], w[48], w[9], pars->GC_2, amp[444]); 
  FFV1_0(w[148], w[48], w[93], pars->GC_11, amp[445]); 
  FFV2_5_0(w[147], w[48], w[29], pars->GC_51, pars->GC_58, amp[446]); 
  FFV1_0(w[150], w[48], w[93], pars->GC_11, amp[447]); 
  FFV1_0(w[151], w[11], w[5], pars->GC_11, amp[448]); 
  FFV1_0(w[148], w[7], w[129], pars->GC_11, amp[449]); 
  FFV1_0(w[148], w[62], w[5], pars->GC_11, amp[450]); 
  FFV1_0(w[152], w[62], w[9], pars->GC_2, amp[451]); 
  FFV1_0(w[152], w[11], w[57], pars->GC_11, amp[452]); 
  FFV1_0(w[146], w[11], w[129], pars->GC_11, amp[453]); 
  FFV1_0(w[151], w[30], w[5], pars->GC_11, amp[454]); 
  FFV1_0(w[150], w[7], w[129], pars->GC_11, amp[455]); 
  FFV1_0(w[150], w[62], w[5], pars->GC_11, amp[456]); 
  FFV2_5_0(w[152], w[62], w[29], pars->GC_51, pars->GC_58, amp[457]); 
  FFV1_0(w[152], w[30], w[57], pars->GC_11, amp[458]); 
  FFV1_0(w[146], w[30], w[129], pars->GC_11, amp[459]); 
  FFV1_0(w[151], w[36], w[9], pars->GC_2, amp[460]); 
  FFV1_0(w[148], w[36], w[57], pars->GC_11, amp[461]); 
  FFV2_5_0(w[151], w[36], w[29], pars->GC_51, pars->GC_58, amp[462]); 
  FFV1_0(w[150], w[36], w[57], pars->GC_11, amp[463]); 
  FFV1_0(w[148], w[101], w[6], pars->GC_11, amp[464]); 
  FFV1_0(w[148], w[67], w[5], pars->GC_11, amp[465]); 
  FFV1_0(w[152], w[66], w[6], pars->GC_11, amp[466]); 
  FFV1_0(w[152], w[67], w[9], pars->GC_2, amp[467]); 
  FFV1_0(w[149], w[66], w[5], pars->GC_11, amp[468]); 
  FFV1_0(w[149], w[101], w[9], pars->GC_2, amp[469]); 
  FFV1_0(w[150], w[101], w[6], pars->GC_11, amp[470]); 
  FFV1_0(w[150], w[67], w[5], pars->GC_11, amp[471]); 
  FFV1_0(w[152], w[68], w[6], pars->GC_11, amp[472]); 
  FFV2_5_0(w[152], w[67], w[29], pars->GC_51, pars->GC_58, amp[473]); 
  FFV1_0(w[149], w[68], w[5], pars->GC_11, amp[474]); 
  FFV2_5_0(w[149], w[101], w[29], pars->GC_51, pars->GC_58, amp[475]); 
  FFV1_0(w[148], w[63], w[32], pars->GC_11, amp[476]); 
  FFV1_0(w[153], w[63], w[9], pars->GC_2, amp[477]); 
  FFV1_0(w[150], w[63], w[32], pars->GC_11, amp[478]); 
  FFV2_5_0(w[153], w[63], w[29], pars->GC_51, pars->GC_58, amp[479]); 
  FFV1_0(w[155], w[11], w[6], pars->GC_11, amp[480]); 
  FFV1_0(w[156], w[11], w[5], pars->GC_11, amp[481]); 
  FFV1_0(w[157], w[11], w[6], pars->GC_11, amp[482]); 
  FFV1_0(w[158], w[11], w[4], pars->GC_11, amp[483]); 
  FFV1_0(w[159], w[11], w[5], pars->GC_11, amp[484]); 
  FFV1_0(w[160], w[11], w[4], pars->GC_11, amp[485]); 
  FFV1_0(w[155], w[30], w[6], pars->GC_11, amp[486]); 
  FFV1_0(w[156], w[30], w[5], pars->GC_11, amp[487]); 
  FFV1_0(w[157], w[30], w[6], pars->GC_11, amp[488]); 
  FFV1_0(w[158], w[30], w[4], pars->GC_11, amp[489]); 
  FFV1_0(w[159], w[30], w[5], pars->GC_11, amp[490]); 
  FFV1_0(w[160], w[30], w[4], pars->GC_11, amp[491]); 
  FFV1_0(w[154], w[11], w[32], pars->GC_11, amp[492]); 
  FFV1_0(w[154], w[35], w[9], pars->GC_2, amp[493]); 
  FFV1_0(w[148], w[7], w[139], pars->GC_11, amp[494]); 
  FFV1_0(w[148], w[35], w[4], pars->GC_11, amp[495]); 
  FFV1_0(w[153], w[11], w[4], pars->GC_11, amp[496]); 
  FFV1_0(w[146], w[11], w[139], pars->GC_11, amp[497]); 
  FFV1_0(w[154], w[30], w[32], pars->GC_11, amp[498]); 
  FFV2_5_0(w[154], w[35], w[29], pars->GC_51, pars->GC_58, amp[499]); 
  FFV1_0(w[150], w[7], w[139], pars->GC_11, amp[500]); 
  FFV1_0(w[150], w[35], w[4], pars->GC_11, amp[501]); 
  FFV1_0(w[153], w[30], w[4], pars->GC_11, amp[502]); 
  FFV1_0(w[146], w[30], w[139], pars->GC_11, amp[503]); 
  FFV1_0(w[154], w[38], w[6], pars->GC_11, amp[504]); 
  FFV1_0(w[154], w[39], w[9], pars->GC_2, amp[505]); 
  FFV1_0(w[148], w[120], w[6], pars->GC_11, amp[506]); 
  FFV1_0(w[148], w[39], w[4], pars->GC_11, amp[507]); 
  FFV1_0(w[149], w[120], w[9], pars->GC_2, amp[508]); 
  FFV1_0(w[149], w[38], w[4], pars->GC_11, amp[509]); 
  FFV1_0(w[154], w[40], w[6], pars->GC_11, amp[510]); 
  FFV2_5_0(w[154], w[39], w[29], pars->GC_51, pars->GC_58, amp[511]); 
  FFV1_0(w[150], w[120], w[6], pars->GC_11, amp[512]); 
  FFV1_0(w[150], w[39], w[4], pars->GC_11, amp[513]); 
  FFV2_5_0(w[149], w[120], w[29], pars->GC_51, pars->GC_58, amp[514]); 
  FFV1_0(w[149], w[40], w[4], pars->GC_11, amp[515]); 
  FFV1_0(w[154], w[50], w[5], pars->GC_11, amp[516]); 
  FFV1_0(w[154], w[51], w[9], pars->GC_2, amp[517]); 
  FFV1_0(w[148], w[90], w[5], pars->GC_11, amp[518]); 
  FFV1_0(w[148], w[51], w[4], pars->GC_11, amp[519]); 
  FFV1_0(w[152], w[90], w[9], pars->GC_2, amp[520]); 
  FFV1_0(w[152], w[50], w[4], pars->GC_11, amp[521]); 
  FFV1_0(w[154], w[52], w[5], pars->GC_11, amp[522]); 
  FFV2_5_0(w[154], w[51], w[29], pars->GC_51, pars->GC_58, amp[523]); 
  FFV1_0(w[150], w[90], w[5], pars->GC_11, amp[524]); 
  FFV1_0(w[150], w[51], w[4], pars->GC_11, amp[525]); 
  FFV2_5_0(w[152], w[90], w[29], pars->GC_51, pars->GC_58, amp[526]); 
  FFV1_0(w[152], w[52], w[4], pars->GC_11, amp[527]); 
  FFV1_0(w[161], w[7], w[9], pars->GC_2, amp[528]); 
  FFV1_0(w[162], w[7], w[9], pars->GC_2, amp[529]); 
  FFV1_0(w[163], w[7], w[9], pars->GC_2, amp[530]); 
  FFV1_0(w[148], w[7], w[140], pars->GC_11, amp[531]); 
  FFV1_0(w[148], w[7], w[141], pars->GC_11, amp[532]); 
  FFV1_0(w[148], w[7], w[142], pars->GC_11, amp[533]); 
  FFV2_5_0(w[161], w[7], w[29], pars->GC_51, pars->GC_58, amp[534]); 
  FFV2_5_0(w[162], w[7], w[29], pars->GC_51, pars->GC_58, amp[535]); 
  FFV2_5_0(w[163], w[7], w[29], pars->GC_51, pars->GC_58, amp[536]); 
  FFV1_0(w[150], w[7], w[140], pars->GC_11, amp[537]); 
  FFV1_0(w[150], w[7], w[141], pars->GC_11, amp[538]); 
  FFV1_0(w[150], w[7], w[142], pars->GC_11, amp[539]); 
  FFV1_0(w[1], w[11], w[165], pars->GC_11, amp[540]); 
  FFV1_0(w[166], w[11], w[6], pars->GC_11, amp[541]); 
  FFV1_0(w[14], w[7], w[165], pars->GC_11, amp[542]); 
  FFV1_0(w[14], w[167], w[6], pars->GC_11, amp[543]); 
  FFV1_0(w[1], w[11], w[168], pars->GC_11, amp[544]); 
  FFV1_0(w[1], w[169], w[125], pars->GC_11, amp[545]); 
  FFV1_0(w[14], w[7], w[168], pars->GC_11, amp[546]); 
  FFV1_0(w[170], w[7], w[125], pars->GC_11, amp[547]); 
  FFV1_0(w[14], w[171], w[6], pars->GC_11, amp[548]); 
  FFV1_0(w[170], w[98], w[6], pars->GC_11, amp[549]); 
  FFV1_0(w[172], w[11], w[6], pars->GC_11, amp[550]); 
  FFV1_0(w[96], w[169], w[6], pars->GC_11, amp[551]); 
  FFV1_0(w[1], w[11], w[173], pars->GC_11, amp[552]); 
  FFV1_0(w[1], w[11], w[174], pars->GC_11, amp[553]); 
  FFV1_0(w[1], w[11], w[175], pars->GC_11, amp[554]); 
  FFV1_0(w[14], w[7], w[173], pars->GC_11, amp[555]); 
  FFV1_0(w[14], w[7], w[174], pars->GC_11, amp[556]); 
  FFV1_0(w[14], w[7], w[175], pars->GC_11, amp[557]); 
  FFV1_0(w[1], w[30], w[165], pars->GC_11, amp[558]); 
  FFV1_0(w[166], w[30], w[6], pars->GC_11, amp[559]); 
  FFV1_0(w[31], w[7], w[165], pars->GC_11, amp[560]); 
  FFV1_0(w[31], w[167], w[6], pars->GC_11, amp[561]); 
  FFV1_0(w[1], w[30], w[168], pars->GC_11, amp[562]); 
  FFV1_0(w[1], w[176], w[125], pars->GC_11, amp[563]); 
  FFV1_0(w[31], w[7], w[168], pars->GC_11, amp[564]); 
  FFV1_0(w[177], w[7], w[125], pars->GC_11, amp[565]); 
  FFV1_0(w[31], w[171], w[6], pars->GC_11, amp[566]); 
  FFV1_0(w[177], w[98], w[6], pars->GC_11, amp[567]); 
  FFV1_0(w[172], w[30], w[6], pars->GC_11, amp[568]); 
  FFV1_0(w[96], w[176], w[6], pars->GC_11, amp[569]); 
  FFV1_0(w[1], w[30], w[173], pars->GC_11, amp[570]); 
  FFV1_0(w[1], w[30], w[174], pars->GC_11, amp[571]); 
  FFV1_0(w[1], w[30], w[175], pars->GC_11, amp[572]); 
  FFV1_0(w[31], w[7], w[173], pars->GC_11, amp[573]); 
  FFV1_0(w[31], w[7], w[174], pars->GC_11, amp[574]); 
  FFV1_0(w[31], w[7], w[175], pars->GC_11, amp[575]); 
  FFV1_0(w[1], w[50], w[164], pars->GC_11, amp[576]); 
  FFV1_0(w[14], w[48], w[164], pars->GC_11, amp[577]); 
  FFV1_0(w[96], w[178], w[9], pars->GC_2, amp[578]); 
  FFV1_0(w[14], w[178], w[93], pars->GC_11, amp[579]); 
  FFV1_0(w[14], w[179], w[0], pars->GC_11, amp[580]); 
  FFV1_0(w[96], w[50], w[0], pars->GC_11, amp[581]); 
  FFV1_0(w[1], w[52], w[164], pars->GC_11, amp[582]); 
  FFV1_0(w[31], w[48], w[164], pars->GC_11, amp[583]); 
  FFV2_5_0(w[96], w[178], w[29], pars->GC_51, pars->GC_58, amp[584]); 
  FFV1_0(w[31], w[178], w[93], pars->GC_11, amp[585]); 
  FFV1_0(w[31], w[179], w[0], pars->GC_11, amp[586]); 
  FFV1_0(w[96], w[52], w[0], pars->GC_11, amp[587]); 
  FFV1_0(w[53], w[7], w[164], pars->GC_11, amp[588]); 
  FFV1_0(w[41], w[11], w[164], pars->GC_11, amp[589]); 
  FFV1_0(w[180], w[98], w[9], pars->GC_2, amp[590]); 
  FFV1_0(w[180], w[11], w[93], pars->GC_11, amp[591]); 
  FFV1_0(w[181], w[11], w[0], pars->GC_11, amp[592]); 
  FFV1_0(w[53], w[98], w[0], pars->GC_11, amp[593]); 
  FFV1_0(w[55], w[7], w[164], pars->GC_11, amp[594]); 
  FFV1_0(w[41], w[30], w[164], pars->GC_11, amp[595]); 
  FFV2_5_0(w[180], w[98], w[29], pars->GC_51, pars->GC_58, amp[596]); 
  FFV1_0(w[180], w[30], w[93], pars->GC_11, amp[597]); 
  FFV1_0(w[181], w[30], w[0], pars->GC_11, amp[598]); 
  FFV1_0(w[55], w[98], w[0], pars->GC_11, amp[599]); 
  FFV1_0(w[1], w[11], w[183], pars->GC_11, amp[600]); 
  FFV1_0(w[184], w[11], w[5], pars->GC_11, amp[601]); 
  FFV1_0(w[14], w[7], w[183], pars->GC_11, amp[602]); 
  FFV1_0(w[14], w[185], w[5], pars->GC_11, amp[603]); 
  FFV1_0(w[1], w[11], w[186], pars->GC_11, amp[604]); 
  FFV1_0(w[1], w[169], w[129], pars->GC_11, amp[605]); 
  FFV1_0(w[14], w[7], w[186], pars->GC_11, amp[606]); 
  FFV1_0(w[170], w[7], w[129], pars->GC_11, amp[607]); 
  FFV1_0(w[14], w[187], w[5], pars->GC_11, amp[608]); 
  FFV1_0(w[170], w[62], w[5], pars->GC_11, amp[609]); 
  FFV1_0(w[188], w[11], w[5], pars->GC_11, amp[610]); 
  FFV1_0(w[60], w[169], w[5], pars->GC_11, amp[611]); 
  FFV1_0(w[1], w[11], w[189], pars->GC_11, amp[612]); 
  FFV1_0(w[1], w[11], w[190], pars->GC_11, amp[613]); 
  FFV1_0(w[1], w[11], w[191], pars->GC_11, amp[614]); 
  FFV1_0(w[14], w[7], w[189], pars->GC_11, amp[615]); 
  FFV1_0(w[14], w[7], w[190], pars->GC_11, amp[616]); 
  FFV1_0(w[14], w[7], w[191], pars->GC_11, amp[617]); 
  FFV1_0(w[1], w[30], w[183], pars->GC_11, amp[618]); 
  FFV1_0(w[184], w[30], w[5], pars->GC_11, amp[619]); 
  FFV1_0(w[31], w[7], w[183], pars->GC_11, amp[620]); 
  FFV1_0(w[31], w[185], w[5], pars->GC_11, amp[621]); 
  FFV1_0(w[1], w[30], w[186], pars->GC_11, amp[622]); 
  FFV1_0(w[1], w[176], w[129], pars->GC_11, amp[623]); 
  FFV1_0(w[31], w[7], w[186], pars->GC_11, amp[624]); 
  FFV1_0(w[177], w[7], w[129], pars->GC_11, amp[625]); 
  FFV1_0(w[31], w[187], w[5], pars->GC_11, amp[626]); 
  FFV1_0(w[177], w[62], w[5], pars->GC_11, amp[627]); 
  FFV1_0(w[188], w[30], w[5], pars->GC_11, amp[628]); 
  FFV1_0(w[60], w[176], w[5], pars->GC_11, amp[629]); 
  FFV1_0(w[1], w[30], w[189], pars->GC_11, amp[630]); 
  FFV1_0(w[1], w[30], w[190], pars->GC_11, amp[631]); 
  FFV1_0(w[1], w[30], w[191], pars->GC_11, amp[632]); 
  FFV1_0(w[31], w[7], w[189], pars->GC_11, amp[633]); 
  FFV1_0(w[31], w[7], w[190], pars->GC_11, amp[634]); 
  FFV1_0(w[31], w[7], w[191], pars->GC_11, amp[635]); 
  FFV1_0(w[1], w[38], w[182], pars->GC_11, amp[636]); 
  FFV1_0(w[14], w[36], w[182], pars->GC_11, amp[637]); 
  FFV1_0(w[60], w[192], w[9], pars->GC_2, amp[638]); 
  FFV1_0(w[14], w[192], w[57], pars->GC_11, amp[639]); 
  FFV1_0(w[14], w[193], w[0], pars->GC_11, amp[640]); 
  FFV1_0(w[60], w[38], w[0], pars->GC_11, amp[641]); 
  FFV1_0(w[1], w[40], w[182], pars->GC_11, amp[642]); 
  FFV1_0(w[31], w[36], w[182], pars->GC_11, amp[643]); 
  FFV2_5_0(w[60], w[192], w[29], pars->GC_51, pars->GC_58, amp[644]); 
  FFV1_0(w[31], w[192], w[57], pars->GC_11, amp[645]); 
  FFV1_0(w[31], w[193], w[0], pars->GC_11, amp[646]); 
  FFV1_0(w[60], w[40], w[0], pars->GC_11, amp[647]); 
  FFV1_0(w[45], w[7], w[182], pars->GC_11, amp[648]); 
  FFV1_0(w[43], w[11], w[182], pars->GC_11, amp[649]); 
  FFV1_0(w[194], w[62], w[9], pars->GC_2, amp[650]); 
  FFV1_0(w[194], w[11], w[57], pars->GC_11, amp[651]); 
  FFV1_0(w[195], w[11], w[0], pars->GC_11, amp[652]); 
  FFV1_0(w[45], w[62], w[0], pars->GC_11, amp[653]); 
  FFV1_0(w[47], w[7], w[182], pars->GC_11, amp[654]); 
  FFV1_0(w[43], w[30], w[182], pars->GC_11, amp[655]); 
  FFV2_5_0(w[194], w[62], w[29], pars->GC_51, pars->GC_58, amp[656]); 
  FFV1_0(w[194], w[30], w[57], pars->GC_11, amp[657]); 
  FFV1_0(w[195], w[30], w[0], pars->GC_11, amp[658]); 
  FFV1_0(w[47], w[62], w[0], pars->GC_11, amp[659]); 
  FFV1_0(w[14], w[197], w[6], pars->GC_11, amp[660]); 
  FFV1_0(w[14], w[198], w[5], pars->GC_11, amp[661]); 
  FFV1_0(w[14], w[199], w[6], pars->GC_11, amp[662]); 
  FFV1_0(w[170], w[101], w[6], pars->GC_11, amp[663]); 
  FFV1_0(w[14], w[200], w[5], pars->GC_11, amp[664]); 
  FFV1_0(w[170], w[67], w[5], pars->GC_11, amp[665]); 
  FFV1_0(w[31], w[197], w[6], pars->GC_11, amp[666]); 
  FFV1_0(w[31], w[198], w[5], pars->GC_11, amp[667]); 
  FFV1_0(w[31], w[199], w[6], pars->GC_11, amp[668]); 
  FFV1_0(w[177], w[101], w[6], pars->GC_11, amp[669]); 
  FFV1_0(w[31], w[200], w[5], pars->GC_11, amp[670]); 
  FFV1_0(w[177], w[67], w[5], pars->GC_11, amp[671]); 
  FFV1_0(w[14], w[196], w[32], pars->GC_11, amp[672]); 
  FFV1_0(w[34], w[196], w[9], pars->GC_2, amp[673]); 
  FFV1_0(w[1], w[66], w[201], pars->GC_11, amp[674]); 
  FFV1_0(w[14], w[63], w[201], pars->GC_11, amp[675]); 
  FFV1_0(w[34], w[66], w[0], pars->GC_11, amp[676]); 
  FFV1_0(w[14], w[202], w[0], pars->GC_11, amp[677]); 
  FFV1_0(w[31], w[196], w[32], pars->GC_11, amp[678]); 
  FFV2_5_0(w[34], w[196], w[29], pars->GC_51, pars->GC_58, amp[679]); 
  FFV1_0(w[1], w[68], w[201], pars->GC_11, amp[680]); 
  FFV1_0(w[31], w[63], w[201], pars->GC_11, amp[681]); 
  FFV1_0(w[34], w[68], w[0], pars->GC_11, amp[682]); 
  FFV1_0(w[31], w[202], w[0], pars->GC_11, amp[683]); 
  FFV1_0(w[45], w[196], w[6], pars->GC_11, amp[684]); 
  FFV1_0(w[46], w[196], w[9], pars->GC_2, amp[685]); 
  FFV1_0(w[194], w[66], w[6], pars->GC_11, amp[686]); 
  FFV1_0(w[194], w[67], w[9], pars->GC_2, amp[687]); 
  FFV1_0(w[46], w[66], w[0], pars->GC_11, amp[688]); 
  FFV1_0(w[45], w[67], w[0], pars->GC_11, amp[689]); 
  FFV1_0(w[47], w[196], w[6], pars->GC_11, amp[690]); 
  FFV2_5_0(w[46], w[196], w[29], pars->GC_51, pars->GC_58, amp[691]); 
  FFV1_0(w[194], w[68], w[6], pars->GC_11, amp[692]); 
  FFV2_5_0(w[194], w[67], w[29], pars->GC_51, pars->GC_58, amp[693]); 
  FFV1_0(w[46], w[68], w[0], pars->GC_11, amp[694]); 
  FFV1_0(w[47], w[67], w[0], pars->GC_11, amp[695]); 
  FFV1_0(w[53], w[196], w[5], pars->GC_11, amp[696]); 
  FFV1_0(w[54], w[196], w[9], pars->GC_2, amp[697]); 
  FFV1_0(w[180], w[66], w[5], pars->GC_11, amp[698]); 
  FFV1_0(w[180], w[101], w[9], pars->GC_2, amp[699]); 
  FFV1_0(w[54], w[66], w[0], pars->GC_11, amp[700]); 
  FFV1_0(w[53], w[101], w[0], pars->GC_11, amp[701]); 
  FFV1_0(w[55], w[196], w[5], pars->GC_11, amp[702]); 
  FFV2_5_0(w[54], w[196], w[29], pars->GC_51, pars->GC_58, amp[703]); 
  FFV1_0(w[180], w[68], w[5], pars->GC_11, amp[704]); 
  FFV2_5_0(w[180], w[101], w[29], pars->GC_51, pars->GC_58, amp[705]); 
  FFV1_0(w[54], w[68], w[0], pars->GC_11, amp[706]); 
  FFV1_0(w[55], w[101], w[0], pars->GC_11, amp[707]); 
  FFV1_0(w[204], w[11], w[6], pars->GC_11, amp[708]); 
  FFV1_0(w[205], w[11], w[5], pars->GC_11, amp[709]); 
  FFV1_0(w[206], w[11], w[6], pars->GC_11, amp[710]); 
  FFV1_0(w[104], w[169], w[6], pars->GC_11, amp[711]); 
  FFV1_0(w[207], w[11], w[5], pars->GC_11, amp[712]); 
  FFV1_0(w[73], w[169], w[5], pars->GC_11, amp[713]); 
  FFV1_0(w[204], w[30], w[6], pars->GC_11, amp[714]); 
  FFV1_0(w[205], w[30], w[5], pars->GC_11, amp[715]); 
  FFV1_0(w[206], w[30], w[6], pars->GC_11, amp[716]); 
  FFV1_0(w[104], w[176], w[6], pars->GC_11, amp[717]); 
  FFV1_0(w[207], w[30], w[5], pars->GC_11, amp[718]); 
  FFV1_0(w[73], w[176], w[5], pars->GC_11, amp[719]); 
  FFV1_0(w[203], w[11], w[32], pars->GC_11, amp[720]); 
  FFV1_0(w[203], w[35], w[9], pars->GC_2, amp[721]); 
  FFV1_0(w[72], w[7], w[201], pars->GC_11, amp[722]); 
  FFV1_0(w[70], w[11], w[201], pars->GC_11, amp[723]); 
  FFV1_0(w[72], w[35], w[0], pars->GC_11, amp[724]); 
  FFV1_0(w[208], w[11], w[0], pars->GC_11, amp[725]); 
  FFV1_0(w[203], w[30], w[32], pars->GC_11, amp[726]); 
  FFV2_5_0(w[203], w[35], w[29], pars->GC_51, pars->GC_58, amp[727]); 
  FFV1_0(w[74], w[7], w[201], pars->GC_11, amp[728]); 
  FFV1_0(w[70], w[30], w[201], pars->GC_11, amp[729]); 
  FFV1_0(w[74], w[35], w[0], pars->GC_11, amp[730]); 
  FFV1_0(w[208], w[30], w[0], pars->GC_11, amp[731]); 
  FFV1_0(w[203], w[38], w[6], pars->GC_11, amp[732]); 
  FFV1_0(w[203], w[39], w[9], pars->GC_2, amp[733]); 
  FFV1_0(w[72], w[192], w[6], pars->GC_11, amp[734]); 
  FFV1_0(w[73], w[192], w[9], pars->GC_2, amp[735]); 
  FFV1_0(w[72], w[39], w[0], pars->GC_11, amp[736]); 
  FFV1_0(w[73], w[38], w[0], pars->GC_11, amp[737]); 
  FFV1_0(w[203], w[40], w[6], pars->GC_11, amp[738]); 
  FFV2_5_0(w[203], w[39], w[29], pars->GC_51, pars->GC_58, amp[739]); 
  FFV1_0(w[74], w[192], w[6], pars->GC_11, amp[740]); 
  FFV2_5_0(w[73], w[192], w[29], pars->GC_51, pars->GC_58, amp[741]); 
  FFV1_0(w[74], w[39], w[0], pars->GC_11, amp[742]); 
  FFV1_0(w[73], w[40], w[0], pars->GC_11, amp[743]); 
  FFV1_0(w[203], w[50], w[5], pars->GC_11, amp[744]); 
  FFV1_0(w[203], w[51], w[9], pars->GC_2, amp[745]); 
  FFV1_0(w[72], w[178], w[5], pars->GC_11, amp[746]); 
  FFV1_0(w[104], w[178], w[9], pars->GC_2, amp[747]); 
  FFV1_0(w[72], w[51], w[0], pars->GC_11, amp[748]); 
  FFV1_0(w[104], w[50], w[0], pars->GC_11, amp[749]); 
  FFV1_0(w[203], w[52], w[5], pars->GC_11, amp[750]); 
  FFV2_5_0(w[203], w[51], w[29], pars->GC_51, pars->GC_58, amp[751]); 
  FFV1_0(w[74], w[178], w[5], pars->GC_11, amp[752]); 
  FFV2_5_0(w[104], w[178], w[29], pars->GC_51, pars->GC_58, amp[753]); 
  FFV1_0(w[74], w[51], w[0], pars->GC_11, amp[754]); 
  FFV1_0(w[104], w[52], w[0], pars->GC_11, amp[755]); 
  FFV1_0(w[1], w[11], w[209], pars->GC_11, amp[756]); 
  FFV1_0(w[210], w[11], w[4], pars->GC_11, amp[757]); 
  FFV1_0(w[14], w[7], w[209], pars->GC_11, amp[758]); 
  FFV1_0(w[14], w[211], w[4], pars->GC_11, amp[759]); 
  FFV1_0(w[1], w[11], w[212], pars->GC_11, amp[760]); 
  FFV1_0(w[1], w[169], w[139], pars->GC_11, amp[761]); 
  FFV1_0(w[14], w[7], w[212], pars->GC_11, amp[762]); 
  FFV1_0(w[170], w[7], w[139], pars->GC_11, amp[763]); 
  FFV1_0(w[34], w[169], w[4], pars->GC_11, amp[764]); 
  FFV1_0(w[213], w[11], w[4], pars->GC_11, amp[765]); 
  FFV1_0(w[170], w[35], w[4], pars->GC_11, amp[766]); 
  FFV1_0(w[14], w[214], w[4], pars->GC_11, amp[767]); 
  FFV1_0(w[1], w[11], w[215], pars->GC_11, amp[768]); 
  FFV1_0(w[1], w[11], w[216], pars->GC_11, amp[769]); 
  FFV1_0(w[1], w[11], w[217], pars->GC_11, amp[770]); 
  FFV1_0(w[14], w[7], w[215], pars->GC_11, amp[771]); 
  FFV1_0(w[14], w[7], w[216], pars->GC_11, amp[772]); 
  FFV1_0(w[14], w[7], w[217], pars->GC_11, amp[773]); 
  FFV1_0(w[1], w[30], w[209], pars->GC_11, amp[774]); 
  FFV1_0(w[210], w[30], w[4], pars->GC_11, amp[775]); 
  FFV1_0(w[31], w[7], w[209], pars->GC_11, amp[776]); 
  FFV1_0(w[31], w[211], w[4], pars->GC_11, amp[777]); 
  FFV1_0(w[1], w[30], w[212], pars->GC_11, amp[778]); 
  FFV1_0(w[1], w[176], w[139], pars->GC_11, amp[779]); 
  FFV1_0(w[31], w[7], w[212], pars->GC_11, amp[780]); 
  FFV1_0(w[177], w[7], w[139], pars->GC_11, amp[781]); 
  FFV1_0(w[34], w[176], w[4], pars->GC_11, amp[782]); 
  FFV1_0(w[213], w[30], w[4], pars->GC_11, amp[783]); 
  FFV1_0(w[177], w[35], w[4], pars->GC_11, amp[784]); 
  FFV1_0(w[31], w[214], w[4], pars->GC_11, amp[785]); 
  FFV1_0(w[1], w[30], w[215], pars->GC_11, amp[786]); 
  FFV1_0(w[1], w[30], w[216], pars->GC_11, amp[787]); 
  FFV1_0(w[1], w[30], w[217], pars->GC_11, amp[788]); 
  FFV1_0(w[31], w[7], w[215], pars->GC_11, amp[789]); 
  FFV1_0(w[31], w[7], w[216], pars->GC_11, amp[790]); 
  FFV1_0(w[31], w[7], w[217], pars->GC_11, amp[791]); 
  FFV1_0(w[14], w[218], w[6], pars->GC_11, amp[792]); 
  FFV1_0(w[14], w[219], w[4], pars->GC_11, amp[793]); 
  FFV1_0(w[14], w[220], w[6], pars->GC_11, amp[794]); 
  FFV1_0(w[170], w[120], w[6], pars->GC_11, amp[795]); 
  FFV1_0(w[170], w[39], w[4], pars->GC_11, amp[796]); 
  FFV1_0(w[14], w[221], w[4], pars->GC_11, amp[797]); 
  FFV1_0(w[31], w[218], w[6], pars->GC_11, amp[798]); 
  FFV1_0(w[31], w[219], w[4], pars->GC_11, amp[799]); 
  FFV1_0(w[31], w[220], w[6], pars->GC_11, amp[800]); 
  FFV1_0(w[177], w[120], w[6], pars->GC_11, amp[801]); 
  FFV1_0(w[177], w[39], w[4], pars->GC_11, amp[802]); 
  FFV1_0(w[31], w[221], w[4], pars->GC_11, amp[803]); 
  FFV1_0(w[91], w[192], w[9], pars->GC_2, amp[804]); 
  FFV1_0(w[53], w[192], w[4], pars->GC_11, amp[805]); 
  FFV1_0(w[180], w[120], w[9], pars->GC_2, amp[806]); 
  FFV1_0(w[180], w[38], w[4], pars->GC_11, amp[807]); 
  FFV1_0(w[53], w[120], w[0], pars->GC_11, amp[808]); 
  FFV1_0(w[91], w[38], w[0], pars->GC_11, amp[809]); 
  FFV2_5_0(w[91], w[192], w[29], pars->GC_51, pars->GC_58, amp[810]); 
  FFV1_0(w[55], w[192], w[4], pars->GC_11, amp[811]); 
  FFV2_5_0(w[180], w[120], w[29], pars->GC_51, pars->GC_58, amp[812]); 
  FFV1_0(w[180], w[40], w[4], pars->GC_11, amp[813]); 
  FFV1_0(w[55], w[120], w[0], pars->GC_11, amp[814]); 
  FFV1_0(w[91], w[40], w[0], pars->GC_11, amp[815]); 
  FFV1_0(w[222], w[11], w[6], pars->GC_11, amp[816]); 
  FFV1_0(w[223], w[11], w[4], pars->GC_11, amp[817]); 
  FFV1_0(w[224], w[11], w[6], pars->GC_11, amp[818]); 
  FFV1_0(w[121], w[169], w[6], pars->GC_11, amp[819]); 
  FFV1_0(w[46], w[169], w[4], pars->GC_11, amp[820]); 
  FFV1_0(w[225], w[11], w[4], pars->GC_11, amp[821]); 
  FFV1_0(w[222], w[30], w[6], pars->GC_11, amp[822]); 
  FFV1_0(w[223], w[30], w[4], pars->GC_11, amp[823]); 
  FFV1_0(w[224], w[30], w[6], pars->GC_11, amp[824]); 
  FFV1_0(w[121], w[176], w[6], pars->GC_11, amp[825]); 
  FFV1_0(w[46], w[176], w[4], pars->GC_11, amp[826]); 
  FFV1_0(w[225], w[30], w[4], pars->GC_11, amp[827]); 
  FFV1_0(w[194], w[90], w[9], pars->GC_2, amp[828]); 
  FFV1_0(w[194], w[50], w[4], pars->GC_11, amp[829]); 
  FFV1_0(w[121], w[178], w[9], pars->GC_2, amp[830]); 
  FFV1_0(w[45], w[178], w[4], pars->GC_11, amp[831]); 
  FFV1_0(w[121], w[50], w[0], pars->GC_11, amp[832]); 
  FFV1_0(w[45], w[90], w[0], pars->GC_11, amp[833]); 
  FFV2_5_0(w[194], w[90], w[29], pars->GC_51, pars->GC_58, amp[834]); 
  FFV1_0(w[194], w[52], w[4], pars->GC_11, amp[835]); 
  FFV2_5_0(w[121], w[178], w[29], pars->GC_51, pars->GC_58, amp[836]); 
  FFV1_0(w[47], w[178], w[4], pars->GC_11, amp[837]); 
  FFV1_0(w[121], w[52], w[0], pars->GC_11, amp[838]); 
  FFV1_0(w[47], w[90], w[0], pars->GC_11, amp[839]); 
  FFV1_0(w[14], w[226], w[5], pars->GC_11, amp[840]); 
  FFV1_0(w[14], w[227], w[4], pars->GC_11, amp[841]); 
  FFV1_0(w[14], w[228], w[5], pars->GC_11, amp[842]); 
  FFV1_0(w[170], w[90], w[5], pars->GC_11, amp[843]); 
  FFV1_0(w[170], w[51], w[4], pars->GC_11, amp[844]); 
  FFV1_0(w[14], w[229], w[4], pars->GC_11, amp[845]); 
  FFV1_0(w[31], w[226], w[5], pars->GC_11, amp[846]); 
  FFV1_0(w[31], w[227], w[4], pars->GC_11, amp[847]); 
  FFV1_0(w[31], w[228], w[5], pars->GC_11, amp[848]); 
  FFV1_0(w[177], w[90], w[5], pars->GC_11, amp[849]); 
  FFV1_0(w[177], w[51], w[4], pars->GC_11, amp[850]); 
  FFV1_0(w[31], w[229], w[4], pars->GC_11, amp[851]); 
  FFV1_0(w[230], w[11], w[5], pars->GC_11, amp[852]); 
  FFV1_0(w[231], w[11], w[4], pars->GC_11, amp[853]); 
  FFV1_0(w[232], w[11], w[5], pars->GC_11, amp[854]); 
  FFV1_0(w[91], w[169], w[5], pars->GC_11, amp[855]); 
  FFV1_0(w[54], w[169], w[4], pars->GC_11, amp[856]); 
  FFV1_0(w[233], w[11], w[4], pars->GC_11, amp[857]); 
  FFV1_0(w[230], w[30], w[5], pars->GC_11, amp[858]); 
  FFV1_0(w[231], w[30], w[4], pars->GC_11, amp[859]); 
  FFV1_0(w[232], w[30], w[5], pars->GC_11, amp[860]); 
  FFV1_0(w[91], w[176], w[5], pars->GC_11, amp[861]); 
  FFV1_0(w[54], w[176], w[4], pars->GC_11, amp[862]); 
  FFV1_0(w[233], w[30], w[4], pars->GC_11, amp[863]); 
  FFV1_0(w[1], w[11], w[237], pars->GC_11, amp[864]); 
  FFV1_0(w[1], w[11], w[238], pars->GC_11, amp[865]); 
  FFV1_0(w[1], w[11], w[239], pars->GC_11, amp[866]); 
  FFV1_0(w[14], w[7], w[237], pars->GC_11, amp[867]); 
  FFV1_0(w[14], w[7], w[238], pars->GC_11, amp[868]); 
  FFV1_0(w[14], w[7], w[239], pars->GC_11, amp[869]); 
  FFV1_0(w[14], w[240], w[6], pars->GC_11, amp[870]); 
  FFV1_0(w[14], w[241], w[6], pars->GC_11, amp[871]); 
  FFV1_0(w[14], w[242], w[6], pars->GC_11, amp[872]); 
  FFV1_0(w[243], w[11], w[6], pars->GC_11, amp[873]); 
  FFV1_0(w[244], w[11], w[6], pars->GC_11, amp[874]); 
  FFV1_0(w[245], w[11], w[6], pars->GC_11, amp[875]); 
  FFV1_0(w[1], w[30], w[237], pars->GC_11, amp[876]); 
  FFV1_0(w[1], w[30], w[238], pars->GC_11, amp[877]); 
  FFV1_0(w[1], w[30], w[239], pars->GC_11, amp[878]); 
  FFV1_0(w[31], w[7], w[237], pars->GC_11, amp[879]); 
  FFV1_0(w[31], w[7], w[238], pars->GC_11, amp[880]); 
  FFV1_0(w[31], w[7], w[239], pars->GC_11, amp[881]); 
  FFV1_0(w[31], w[240], w[6], pars->GC_11, amp[882]); 
  FFV1_0(w[31], w[241], w[6], pars->GC_11, amp[883]); 
  FFV1_0(w[31], w[242], w[6], pars->GC_11, amp[884]); 
  FFV1_0(w[243], w[30], w[6], pars->GC_11, amp[885]); 
  FFV1_0(w[244], w[30], w[6], pars->GC_11, amp[886]); 
  FFV1_0(w[245], w[30], w[6], pars->GC_11, amp[887]); 
  FFV1_0(w[1], w[246], w[9], pars->GC_2, amp[888]); 
  FFV1_0(w[1], w[247], w[9], pars->GC_2, amp[889]); 
  FFV1_0(w[1], w[248], w[9], pars->GC_2, amp[890]); 
  FFV1_0(w[243], w[48], w[9], pars->GC_2, amp[891]); 
  FFV1_0(w[244], w[48], w[9], pars->GC_2, amp[892]); 
  FFV1_0(w[245], w[48], w[9], pars->GC_2, amp[893]); 
  FFV2_5_0(w[1], w[246], w[29], pars->GC_51, pars->GC_58, amp[894]); 
  FFV2_5_0(w[1], w[247], w[29], pars->GC_51, pars->GC_58, amp[895]); 
  FFV2_5_0(w[1], w[248], w[29], pars->GC_51, pars->GC_58, amp[896]); 
  FFV2_5_0(w[243], w[48], w[29], pars->GC_51, pars->GC_58, amp[897]); 
  FFV2_5_0(w[244], w[48], w[29], pars->GC_51, pars->GC_58, amp[898]); 
  FFV2_5_0(w[245], w[48], w[29], pars->GC_51, pars->GC_58, amp[899]); 
  FFV1_0(w[249], w[7], w[9], pars->GC_2, amp[900]); 
  FFV1_0(w[250], w[7], w[9], pars->GC_2, amp[901]); 
  FFV1_0(w[251], w[7], w[9], pars->GC_2, amp[902]); 
  FFV1_0(w[41], w[240], w[9], pars->GC_2, amp[903]); 
  FFV1_0(w[41], w[241], w[9], pars->GC_2, amp[904]); 
  FFV1_0(w[41], w[242], w[9], pars->GC_2, amp[905]); 
  FFV2_5_0(w[249], w[7], w[29], pars->GC_51, pars->GC_58, amp[906]); 
  FFV2_5_0(w[250], w[7], w[29], pars->GC_51, pars->GC_58, amp[907]); 
  FFV2_5_0(w[251], w[7], w[29], pars->GC_51, pars->GC_58, amp[908]); 
  FFV2_5_0(w[41], w[240], w[29], pars->GC_51, pars->GC_58, amp[909]); 
  FFV2_5_0(w[41], w[241], w[29], pars->GC_51, pars->GC_58, amp[910]); 
  FFV2_5_0(w[41], w[242], w[29], pars->GC_51, pars->GC_58, amp[911]); 
  FFV1_0(w[1], w[11], w[255], pars->GC_11, amp[912]); 
  FFV1_0(w[1], w[11], w[256], pars->GC_11, amp[913]); 
  FFV1_0(w[1], w[11], w[257], pars->GC_11, amp[914]); 
  FFV1_0(w[14], w[7], w[255], pars->GC_11, amp[915]); 
  FFV1_0(w[14], w[7], w[256], pars->GC_11, amp[916]); 
  FFV1_0(w[14], w[7], w[257], pars->GC_11, amp[917]); 
  FFV1_0(w[14], w[258], w[5], pars->GC_11, amp[918]); 
  FFV1_0(w[14], w[259], w[5], pars->GC_11, amp[919]); 
  FFV1_0(w[14], w[260], w[5], pars->GC_11, amp[920]); 
  FFV1_0(w[261], w[11], w[5], pars->GC_11, amp[921]); 
  FFV1_0(w[262], w[11], w[5], pars->GC_11, amp[922]); 
  FFV1_0(w[263], w[11], w[5], pars->GC_11, amp[923]); 
  FFV1_0(w[1], w[30], w[255], pars->GC_11, amp[924]); 
  FFV1_0(w[1], w[30], w[256], pars->GC_11, amp[925]); 
  FFV1_0(w[1], w[30], w[257], pars->GC_11, amp[926]); 
  FFV1_0(w[31], w[7], w[255], pars->GC_11, amp[927]); 
  FFV1_0(w[31], w[7], w[256], pars->GC_11, amp[928]); 
  FFV1_0(w[31], w[7], w[257], pars->GC_11, amp[929]); 
  FFV1_0(w[31], w[258], w[5], pars->GC_11, amp[930]); 
  FFV1_0(w[31], w[259], w[5], pars->GC_11, amp[931]); 
  FFV1_0(w[31], w[260], w[5], pars->GC_11, amp[932]); 
  FFV1_0(w[261], w[30], w[5], pars->GC_11, amp[933]); 
  FFV1_0(w[262], w[30], w[5], pars->GC_11, amp[934]); 
  FFV1_0(w[263], w[30], w[5], pars->GC_11, amp[935]); 
  FFV1_0(w[1], w[264], w[9], pars->GC_2, amp[936]); 
  FFV1_0(w[1], w[265], w[9], pars->GC_2, amp[937]); 
  FFV1_0(w[1], w[266], w[9], pars->GC_2, amp[938]); 
  FFV1_0(w[261], w[36], w[9], pars->GC_2, amp[939]); 
  FFV1_0(w[262], w[36], w[9], pars->GC_2, amp[940]); 
  FFV1_0(w[263], w[36], w[9], pars->GC_2, amp[941]); 
  FFV2_5_0(w[1], w[264], w[29], pars->GC_51, pars->GC_58, amp[942]); 
  FFV2_5_0(w[1], w[265], w[29], pars->GC_51, pars->GC_58, amp[943]); 
  FFV2_5_0(w[1], w[266], w[29], pars->GC_51, pars->GC_58, amp[944]); 
  FFV2_5_0(w[261], w[36], w[29], pars->GC_51, pars->GC_58, amp[945]); 
  FFV2_5_0(w[262], w[36], w[29], pars->GC_51, pars->GC_58, amp[946]); 
  FFV2_5_0(w[263], w[36], w[29], pars->GC_51, pars->GC_58, amp[947]); 
  FFV1_0(w[267], w[7], w[9], pars->GC_2, amp[948]); 
  FFV1_0(w[268], w[7], w[9], pars->GC_2, amp[949]); 
  FFV1_0(w[269], w[7], w[9], pars->GC_2, amp[950]); 
  FFV1_0(w[43], w[258], w[9], pars->GC_2, amp[951]); 
  FFV1_0(w[43], w[259], w[9], pars->GC_2, amp[952]); 
  FFV1_0(w[43], w[260], w[9], pars->GC_2, amp[953]); 
  FFV2_5_0(w[267], w[7], w[29], pars->GC_51, pars->GC_58, amp[954]); 
  FFV2_5_0(w[268], w[7], w[29], pars->GC_51, pars->GC_58, amp[955]); 
  FFV2_5_0(w[269], w[7], w[29], pars->GC_51, pars->GC_58, amp[956]); 
  FFV2_5_0(w[43], w[258], w[29], pars->GC_51, pars->GC_58, amp[957]); 
  FFV2_5_0(w[43], w[259], w[29], pars->GC_51, pars->GC_58, amp[958]); 
  FFV2_5_0(w[43], w[260], w[29], pars->GC_51, pars->GC_58, amp[959]); 
  FFV1_0(w[1], w[273], w[9], pars->GC_2, amp[960]); 
  FFV1_0(w[1], w[274], w[9], pars->GC_2, amp[961]); 
  FFV1_0(w[1], w[275], w[9], pars->GC_2, amp[962]); 
  FFV1_0(w[276], w[63], w[9], pars->GC_2, amp[963]); 
  FFV1_0(w[277], w[63], w[9], pars->GC_2, amp[964]); 
  FFV1_0(w[278], w[63], w[9], pars->GC_2, amp[965]); 
  FFV2_5_0(w[1], w[273], w[29], pars->GC_51, pars->GC_58, amp[966]); 
  FFV2_5_0(w[1], w[274], w[29], pars->GC_51, pars->GC_58, amp[967]); 
  FFV2_5_0(w[1], w[275], w[29], pars->GC_51, pars->GC_58, amp[968]); 
  FFV2_5_0(w[276], w[63], w[29], pars->GC_51, pars->GC_58, amp[969]); 
  FFV2_5_0(w[277], w[63], w[29], pars->GC_51, pars->GC_58, amp[970]); 
  FFV2_5_0(w[278], w[63], w[29], pars->GC_51, pars->GC_58, amp[971]); 
  FFV1_0(w[279], w[7], w[9], pars->GC_2, amp[972]); 
  FFV1_0(w[280], w[7], w[9], pars->GC_2, amp[973]); 
  FFV1_0(w[281], w[7], w[9], pars->GC_2, amp[974]); 
  FFV1_0(w[70], w[282], w[9], pars->GC_2, amp[975]); 
  FFV1_0(w[70], w[283], w[9], pars->GC_2, amp[976]); 
  FFV1_0(w[70], w[284], w[9], pars->GC_2, amp[977]); 
  FFV2_5_0(w[279], w[7], w[29], pars->GC_51, pars->GC_58, amp[978]); 
  FFV2_5_0(w[280], w[7], w[29], pars->GC_51, pars->GC_58, amp[979]); 
  FFV2_5_0(w[281], w[7], w[29], pars->GC_51, pars->GC_58, amp[980]); 
  FFV2_5_0(w[70], w[282], w[29], pars->GC_51, pars->GC_58, amp[981]); 
  FFV2_5_0(w[70], w[283], w[29], pars->GC_51, pars->GC_58, amp[982]); 
  FFV2_5_0(w[70], w[284], w[29], pars->GC_51, pars->GC_58, amp[983]); 
  FFV1_0(w[1], w[11], w[285], pars->GC_11, amp[984]); 
  FFV1_0(w[1], w[11], w[286], pars->GC_11, amp[985]); 
  FFV1_0(w[1], w[11], w[287], pars->GC_11, amp[986]); 
  FFV1_0(w[14], w[7], w[285], pars->GC_11, amp[987]); 
  FFV1_0(w[14], w[7], w[286], pars->GC_11, amp[988]); 
  FFV1_0(w[14], w[7], w[287], pars->GC_11, amp[989]); 
  FFV1_0(w[14], w[282], w[4], pars->GC_11, amp[990]); 
  FFV1_0(w[14], w[283], w[4], pars->GC_11, amp[991]); 
  FFV1_0(w[14], w[284], w[4], pars->GC_11, amp[992]); 
  FFV1_0(w[276], w[11], w[4], pars->GC_11, amp[993]); 
  FFV1_0(w[277], w[11], w[4], pars->GC_11, amp[994]); 
  FFV1_0(w[278], w[11], w[4], pars->GC_11, amp[995]); 
  FFV1_0(w[1], w[30], w[285], pars->GC_11, amp[996]); 
  FFV1_0(w[1], w[30], w[286], pars->GC_11, amp[997]); 
  FFV1_0(w[1], w[30], w[287], pars->GC_11, amp[998]); 
  FFV1_0(w[31], w[7], w[285], pars->GC_11, amp[999]); 
  FFV1_0(w[31], w[7], w[286], pars->GC_11, amp[1000]); 
  FFV1_0(w[31], w[7], w[287], pars->GC_11, amp[1001]); 
  FFV1_0(w[31], w[282], w[4], pars->GC_11, amp[1002]); 
  FFV1_0(w[31], w[283], w[4], pars->GC_11, amp[1003]); 
  FFV1_0(w[31], w[284], w[4], pars->GC_11, amp[1004]); 
  FFV1_0(w[276], w[30], w[4], pars->GC_11, amp[1005]); 
  FFV1_0(w[277], w[30], w[4], pars->GC_11, amp[1006]); 
  FFV1_0(w[278], w[30], w[4], pars->GC_11, amp[1007]); 
  FFV1_0(w[1], w[11], w[288], pars->GC_11, amp[1008]); 
  FFV1_0(w[1], w[11], w[289], pars->GC_11, amp[1009]); 
  FFV1_0(w[1], w[11], w[290], pars->GC_11, amp[1010]); 
  FFV1_0(w[14], w[7], w[288], pars->GC_11, amp[1011]); 
  FFV1_0(w[14], w[7], w[289], pars->GC_11, amp[1012]); 
  FFV1_0(w[14], w[7], w[290], pars->GC_11, amp[1013]); 
  FFV1_0(w[14], w[291], w[0], pars->GC_11, amp[1014]); 
  FFV1_0(w[14], w[292], w[0], pars->GC_11, amp[1015]); 
  FFV1_0(w[14], w[293], w[0], pars->GC_11, amp[1016]); 
  FFV1_0(w[294], w[11], w[0], pars->GC_11, amp[1017]); 
  FFV1_0(w[295], w[11], w[0], pars->GC_11, amp[1018]); 
  FFV1_0(w[296], w[11], w[0], pars->GC_11, amp[1019]); 
  FFV1_0(w[1], w[30], w[288], pars->GC_11, amp[1020]); 
  FFV1_0(w[1], w[30], w[289], pars->GC_11, amp[1021]); 
  FFV1_0(w[1], w[30], w[290], pars->GC_11, amp[1022]); 
  FFV1_0(w[31], w[7], w[288], pars->GC_11, amp[1023]); 
  FFV1_0(w[31], w[7], w[289], pars->GC_11, amp[1024]); 
  FFV1_0(w[31], w[7], w[290], pars->GC_11, amp[1025]); 
  FFV1_0(w[31], w[291], w[0], pars->GC_11, amp[1026]); 
  FFV1_0(w[31], w[292], w[0], pars->GC_11, amp[1027]); 
  FFV1_0(w[31], w[293], w[0], pars->GC_11, amp[1028]); 
  FFV1_0(w[294], w[30], w[0], pars->GC_11, amp[1029]); 
  FFV1_0(w[295], w[30], w[0], pars->GC_11, amp[1030]); 
  FFV1_0(w[296], w[30], w[0], pars->GC_11, amp[1031]); 

}
double CPPProcess::matrix_gu_taptamgggu_no_h() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 1032; 
  const int ncolor = 24; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {54, 54, 54, 54, 54, 54, 54, 54, 54, 54,
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54};
  static const double cf[ncolor][ncolor] = {{512, -64, -64, 8, 8, 80, -64, 8,
      8, -1, -1, -10, 8, -1, 80, -10, 71, 62, -1, -10, -10, 62, 62, -28}, {-64,
      512, 8, 80, -64, 8, 8, -64, -1, -10, 8, -1, -1, -10, -10, 62, 62, -28, 8,
      -1, 80, -10, 71, 62}, {-64, 8, 512, -64, 80, 8, 8, -1, 80, -10, 71, 62,
      -64, 8, 8, -1, -1, -10, -10, -1, 62, -28, -10, 62}, {8, 80, -64, 512, 8,
      -64, -1, -10, -10, 62, 62, -28, 8, -64, -1, -10, 8, -1, -1, 8, 71, 62,
      80, -10}, {8, -64, 80, 8, 512, -64, -1, 8, 71, 62, 80, -10, -10, -1, 62,
      -28, -10, 62, -64, 8, 8, -1, -1, -10}, {80, 8, 8, -64, -64, 512, -10, -1,
      62, -28, -10, 62, -1, 8, 71, 62, 80, -10, 8, -64, -1, -10, 8, -1}, {-64,
      8, 8, -1, -1, -10, 512, -64, -64, 8, 8, 80, 80, -10, 8, -1, 62, 71, -10,
      62, -1, -10, -28, 62}, {8, -64, -1, -10, 8, -1, -64, 512, 8, 80, -64, 8,
      -10, 62, -1, -10, -28, 62, 80, -10, 8, -1, 62, 71}, {8, -1, 80, -10, 71,
      62, -64, 8, 512, -64, 80, 8, 8, -1, -64, 8, -10, -1, 62, -28, -10, -1,
      62, -10}, {-1, -10, -10, 62, 62, -28, 8, 80, -64, 512, 8, -64, -1, -10,
      8, -64, -1, 8, 71, 62, -1, 8, -10, 80}, {-1, 8, 71, 62, 80, -10, 8, -64,
      80, 8, 512, -64, 62, -28, -10, -1, 62, -10, 8, -1, -64, 8, -10, -1},
      {-10, -1, 62, -28, -10, 62, 80, 8, 8, -64, -64, 512, 71, 62, -1, 8, -10,
      80, -1, -10, 8, -64, -1, 8}, {8, -1, -64, 8, -10, -1, 80, -10, 8, -1, 62,
      71, 512, -64, -64, 8, 8, 80, 62, -10, -28, 62, -1, -10}, {-1, -10, 8,
      -64, -1, 8, -10, 62, -1, -10, -28, 62, -64, 512, 8, 80, -64, 8, -10, 80,
      62, 71, 8, -1}, {80, -10, 8, -1, 62, 71, 8, -1, -64, 8, -10, -1, -64, 8,
      512, -64, 80, 8, -28, 62, 62, -10, -10, -1}, {-10, 62, -1, -10, -28, 62,
      -1, -10, 8, -64, -1, 8, 8, 80, -64, 512, 8, -64, 62, 71, -10, 80, -1, 8},
      {71, 62, -1, 8, -10, 80, 62, -28, -10, -1, 62, -10, 8, -64, 80, 8, 512,
      -64, -1, 8, -10, -1, -64, 8}, {62, -28, -10, -1, 62, -10, 71, 62, -1, 8,
      -10, 80, 80, 8, 8, -64, -64, 512, -10, -1, -1, 8, 8, -64}, {-1, 8, -10,
      -1, -64, 8, -10, 80, 62, 71, 8, -1, 62, -10, -28, 62, -1, -10, 512, -64,
      -64, 8, 8, 80}, {-10, -1, -1, 8, 8, -64, 62, -10, -28, 62, -1, -10, -10,
      80, 62, 71, 8, -1, -64, 512, 8, 80, -64, 8}, {-10, 80, 62, 71, 8, -1, -1,
      8, -10, -1, -64, 8, -28, 62, 62, -10, -10, -1, -64, 8, 512, -64, 80, 8},
      {62, -10, -28, 62, -1, -10, -10, -1, -1, 8, 8, -64, 62, 71, -10, 80, -1,
      8, 8, 80, -64, 512, 8, -64}, {62, 71, -10, 80, -1, 8, -28, 62, 62, -10,
      -10, -1, -1, 8, -10, -1, -64, 8, 8, -64, 80, 8, 512, -64}, {-28, 62, 62,
      -10, -10, -1, 62, 71, -10, 80, -1, 8, -10, -1, -1, 8, 8, -64, 80, 8, 8,
      -64, -64, 512}};

  // Calculate color flows
  jamp[0] = +std::complex<double> (0, 1) * amp[0] + std::complex<double> (0, 1)
      * amp[2] - amp[3] - std::complex<double> (0, 1) * amp[8] -
      std::complex<double> (0, 1) * amp[14] + std::complex<double> (0, 1) *
      amp[12] - std::complex<double> (0, 1) * amp[17] + std::complex<double>
      (0, 1) * amp[15] + std::complex<double> (0, 1) * amp[18] +
      std::complex<double> (0, 1) * amp[20] - amp[21] - std::complex<double>
      (0, 1) * amp[26] - std::complex<double> (0, 1) * amp[32] +
      std::complex<double> (0, 1) * amp[30] - std::complex<double> (0, 1) *
      amp[35] + std::complex<double> (0, 1) * amp[33] + std::complex<double>
      (0, 1) * amp[36] + std::complex<double> (0, 1) * amp[37] - amp[38] -
      amp[39] - amp[42] + std::complex<double> (0, 1) * amp[44] +
      std::complex<double> (0, 1) * amp[45] - amp[46] - amp[47] - amp[50] -
      amp[96] - amp[97] - std::complex<double> (0, 1) * amp[99] -
      std::complex<double> (0, 1) * amp[100] - std::complex<double> (0, 1) *
      amp[101] - amp[102] - amp[103] - std::complex<double> (0, 1) * amp[105] -
      std::complex<double> (0, 1) * amp[106] - std::complex<double> (0, 1) *
      amp[107] - std::complex<double> (0, 1) * amp[324] - amp[325] - amp[329] -
      std::complex<double> (0, 1) * amp[330] - amp[331] - amp[335] -
      std::complex<double> (0, 1) * amp[336] - std::complex<double> (0, 1) *
      amp[337] - std::complex<double> (0, 1) * amp[338] - std::complex<double>
      (0, 1) * amp[339] + amp[372] + amp[378] - std::complex<double> (0, 1) *
      amp[384] - std::complex<double> (0, 1) * amp[385] - amp[386] -
      std::complex<double> (0, 1) * amp[387] - amp[389] - std::complex<double>
      (0, 1) * amp[390] - std::complex<double> (0, 1) * amp[391] - amp[392] -
      std::complex<double> (0, 1) * amp[393] - amp[395] + amp[408] + amp[409] +
      amp[411] + amp[414] + amp[415] + amp[417] - amp[420] + amp[422] +
      amp[425] - amp[423] - amp[426] + amp[428] + amp[431] - amp[429] +
      std::complex<double> (0, 1) * amp[540] + std::complex<double> (0, 1) *
      amp[542] - amp[543] + std::complex<double> (0, 1) * amp[544] - amp[545] +
      std::complex<double> (0, 1) * amp[546] - std::complex<double> (0, 1) *
      amp[554] + std::complex<double> (0, 1) * amp[552] - std::complex<double>
      (0, 1) * amp[557] + std::complex<double> (0, 1) * amp[555] +
      std::complex<double> (0, 1) * amp[558] + std::complex<double> (0, 1) *
      amp[560] - amp[561] + std::complex<double> (0, 1) * amp[562] - amp[563] +
      std::complex<double> (0, 1) * amp[564] - std::complex<double> (0, 1) *
      amp[572] + std::complex<double> (0, 1) * amp[570] - std::complex<double>
      (0, 1) * amp[575] + std::complex<double> (0, 1) * amp[573] - amp[588] -
      amp[589] - std::complex<double> (0, 1) * amp[592] - amp[594] - amp[595] -
      std::complex<double> (0, 1) * amp[598] + std::complex<double> (0, 1) *
      amp[760] - amp[761] + std::complex<double> (0, 1) * amp[762] -
      std::complex<double> (0, 1) * amp[764] - std::complex<double> (0, 1) *
      amp[770] + std::complex<double> (0, 1) * amp[768] - std::complex<double>
      (0, 1) * amp[773] + std::complex<double> (0, 1) * amp[771] +
      std::complex<double> (0, 1) * amp[778] - amp[779] + std::complex<double>
      (0, 1) * amp[780] - std::complex<double> (0, 1) * amp[782] -
      std::complex<double> (0, 1) * amp[788] + std::complex<double> (0, 1) *
      amp[786] - std::complex<double> (0, 1) * amp[791] + std::complex<double>
      (0, 1) * amp[789] + amp[856] + amp[862] - std::complex<double> (0, 1) *
      amp[866] + std::complex<double> (0, 1) * amp[864] - std::complex<double>
      (0, 1) * amp[869] + std::complex<double> (0, 1) * amp[867] + amp[872] -
      amp[870] - std::complex<double> (0, 1) * amp[878] + std::complex<double>
      (0, 1) * amp[876] - std::complex<double> (0, 1) * amp[881] +
      std::complex<double> (0, 1) * amp[879] + amp[884] - amp[882] - amp[900] +
      amp[902] - amp[903] + amp[905] - amp[906] + amp[908] - amp[909] +
      amp[911] - std::complex<double> (0, 1) * amp[1010] + std::complex<double>
      (0, 1) * amp[1008] - std::complex<double> (0, 1) * amp[1013] +
      std::complex<double> (0, 1) * amp[1011] + amp[1019] - amp[1017] -
      std::complex<double> (0, 1) * amp[1022] + std::complex<double> (0, 1) *
      amp[1020] - std::complex<double> (0, 1) * amp[1025] +
      std::complex<double> (0, 1) * amp[1023] + amp[1031] - amp[1029];
  jamp[1] = +std::complex<double> (0, 1) * amp[4] + std::complex<double> (0, 1)
      * amp[6] - amp[7] - std::complex<double> (0, 1) * amp[9] -
      std::complex<double> (0, 1) * amp[12] - std::complex<double> (0, 1) *
      amp[13] - std::complex<double> (0, 1) * amp[15] - std::complex<double>
      (0, 1) * amp[16] + std::complex<double> (0, 1) * amp[22] +
      std::complex<double> (0, 1) * amp[24] - amp[25] - std::complex<double>
      (0, 1) * amp[27] - std::complex<double> (0, 1) * amp[30] -
      std::complex<double> (0, 1) * amp[31] - std::complex<double> (0, 1) *
      amp[33] - std::complex<double> (0, 1) * amp[34] - std::complex<double>
      (0, 1) * amp[36] - std::complex<double> (0, 1) * amp[37] + amp[38] +
      amp[39] + amp[42] - std::complex<double> (0, 1) * amp[44] -
      std::complex<double> (0, 1) * amp[45] + amp[46] + amp[47] + amp[50] -
      amp[69] - amp[70] - std::complex<double> (0, 1) * amp[71] -
      std::complex<double> (0, 1) * amp[72] - std::complex<double> (0, 1) *
      amp[73] - amp[75] - amp[76] - std::complex<double> (0, 1) * amp[77] -
      std::complex<double> (0, 1) * amp[78] - std::complex<double> (0, 1) *
      amp[79] - std::complex<double> (0, 1) * amp[340] - amp[341] - amp[345] -
      std::complex<double> (0, 1) * amp[346] - amp[347] - amp[351] -
      std::complex<double> (0, 1) * amp[352] - std::complex<double> (0, 1) *
      amp[353] - std::complex<double> (0, 1) * amp[354] - std::complex<double>
      (0, 1) * amp[355] + amp[373] + amp[379] + std::complex<double> (0, 1) *
      amp[384] + std::complex<double> (0, 1) * amp[385] + amp[386] +
      std::complex<double> (0, 1) * amp[387] + amp[389] + std::complex<double>
      (0, 1) * amp[390] + std::complex<double> (0, 1) * amp[391] + amp[392] +
      std::complex<double> (0, 1) * amp[393] + amp[395] + amp[396] + amp[397] +
      amp[399] + amp[402] + amp[403] + amp[405] + amp[420] + amp[421] +
      amp[423] + amp[424] + amp[426] + amp[427] + amp[429] + amp[430] +
      std::complex<double> (0, 1) * amp[600] + std::complex<double> (0, 1) *
      amp[602] - amp[603] + std::complex<double> (0, 1) * amp[604] - amp[605] +
      std::complex<double> (0, 1) * amp[606] - std::complex<double> (0, 1) *
      amp[614] + std::complex<double> (0, 1) * amp[612] - std::complex<double>
      (0, 1) * amp[617] + std::complex<double> (0, 1) * amp[615] +
      std::complex<double> (0, 1) * amp[618] + std::complex<double> (0, 1) *
      amp[620] - amp[621] + std::complex<double> (0, 1) * amp[622] - amp[623] +
      std::complex<double> (0, 1) * amp[624] - std::complex<double> (0, 1) *
      amp[632] + std::complex<double> (0, 1) * amp[630] - std::complex<double>
      (0, 1) * amp[635] + std::complex<double> (0, 1) * amp[633] - amp[648] -
      amp[649] - std::complex<double> (0, 1) * amp[652] - amp[654] - amp[655] -
      std::complex<double> (0, 1) * amp[658] - std::complex<double> (0, 1) *
      amp[760] + amp[761] - std::complex<double> (0, 1) * amp[762] +
      std::complex<double> (0, 1) * amp[764] + std::complex<double> (0, 1) *
      amp[770] - std::complex<double> (0, 1) * amp[768] + std::complex<double>
      (0, 1) * amp[773] - std::complex<double> (0, 1) * amp[771] -
      std::complex<double> (0, 1) * amp[778] + amp[779] - std::complex<double>
      (0, 1) * amp[780] + std::complex<double> (0, 1) * amp[782] +
      std::complex<double> (0, 1) * amp[788] - std::complex<double> (0, 1) *
      amp[786] + std::complex<double> (0, 1) * amp[791] - std::complex<double>
      (0, 1) * amp[789] + amp[820] + amp[826] - std::complex<double> (0, 1) *
      amp[914] + std::complex<double> (0, 1) * amp[912] - std::complex<double>
      (0, 1) * amp[917] + std::complex<double> (0, 1) * amp[915] + amp[920] -
      amp[918] - std::complex<double> (0, 1) * amp[926] + std::complex<double>
      (0, 1) * amp[924] - std::complex<double> (0, 1) * amp[929] +
      std::complex<double> (0, 1) * amp[927] + amp[932] - amp[930] - amp[948] +
      amp[950] - amp[951] + amp[953] - amp[954] + amp[956] - amp[957] +
      amp[959] - std::complex<double> (0, 1) * amp[1008] - std::complex<double>
      (0, 1) * amp[1009] - std::complex<double> (0, 1) * amp[1011] -
      std::complex<double> (0, 1) * amp[1012] + amp[1017] + amp[1018] -
      std::complex<double> (0, 1) * amp[1020] - std::complex<double> (0, 1) *
      amp[1021] - std::complex<double> (0, 1) * amp[1023] -
      std::complex<double> (0, 1) * amp[1024] + amp[1029] + amp[1030];
  jamp[2] = +std::complex<double> (0, 1) * amp[108] + std::complex<double> (0,
      1) * amp[109] - amp[110] - amp[111] - amp[115] + std::complex<double> (0,
      1) * amp[116] + std::complex<double> (0, 1) * amp[117] - amp[118] -
      amp[119] - amp[123] + std::complex<double> (0, 1) * amp[156] +
      std::complex<double> (0, 1) * amp[158] - amp[159] - std::complex<double>
      (0, 1) * amp[164] - std::complex<double> (0, 1) * amp[170] +
      std::complex<double> (0, 1) * amp[168] - std::complex<double> (0, 1) *
      amp[173] + std::complex<double> (0, 1) * amp[171] + std::complex<double>
      (0, 1) * amp[174] + std::complex<double> (0, 1) * amp[176] - amp[177] -
      std::complex<double> (0, 1) * amp[182] - std::complex<double> (0, 1) *
      amp[188] + std::complex<double> (0, 1) * amp[186] - std::complex<double>
      (0, 1) * amp[191] + std::complex<double> (0, 1) * amp[189] - amp[204] -
      amp[205] - std::complex<double> (0, 1) * amp[207] - std::complex<double>
      (0, 1) * amp[208] - std::complex<double> (0, 1) * amp[209] - amp[210] -
      amp[211] - std::complex<double> (0, 1) * amp[213] - std::complex<double>
      (0, 1) * amp[214] - std::complex<double> (0, 1) * amp[215] +
      std::complex<double> (0, 1) * amp[324] + amp[325] + amp[329] +
      std::complex<double> (0, 1) * amp[330] + amp[331] + amp[335] +
      std::complex<double> (0, 1) * amp[336] + std::complex<double> (0, 1) *
      amp[337] + std::complex<double> (0, 1) * amp[338] + std::complex<double>
      (0, 1) * amp[339] + amp[341] - std::complex<double> (0, 1) * amp[342] -
      std::complex<double> (0, 1) * amp[343] - std::complex<double> (0, 1) *
      amp[344] + amp[345] + amp[347] - std::complex<double> (0, 1) * amp[348] -
      std::complex<double> (0, 1) * amp[349] - std::complex<double> (0, 1) *
      amp[350] + amp[351] + amp[374] + amp[380] + amp[410] + amp[412] +
      amp[413] + amp[416] + amp[418] + amp[419] - amp[421] - amp[422] -
      amp[425] - amp[424] - amp[427] - amp[428] - amp[431] - amp[430] -
      std::complex<double> (0, 1) * amp[540] - std::complex<double> (0, 1) *
      amp[542] + amp[543] - std::complex<double> (0, 1) * amp[544] + amp[545] -
      std::complex<double> (0, 1) * amp[546] + std::complex<double> (0, 1) *
      amp[554] - std::complex<double> (0, 1) * amp[552] + std::complex<double>
      (0, 1) * amp[557] - std::complex<double> (0, 1) * amp[555] -
      std::complex<double> (0, 1) * amp[558] - std::complex<double> (0, 1) *
      amp[560] + amp[561] - std::complex<double> (0, 1) * amp[562] + amp[563] -
      std::complex<double> (0, 1) * amp[564] + std::complex<double> (0, 1) *
      amp[572] - std::complex<double> (0, 1) * amp[570] + std::complex<double>
      (0, 1) * amp[575] - std::complex<double> (0, 1) * amp[573] + amp[588] +
      amp[589] + std::complex<double> (0, 1) * amp[592] + amp[594] + amp[595] +
      std::complex<double> (0, 1) * amp[598] - std::complex<double> (0, 1) *
      amp[604] + amp[605] - std::complex<double> (0, 1) * amp[606] -
      std::complex<double> (0, 1) * amp[611] - std::complex<double> (0, 1) *
      amp[612] - std::complex<double> (0, 1) * amp[613] - std::complex<double>
      (0, 1) * amp[615] - std::complex<double> (0, 1) * amp[616] -
      std::complex<double> (0, 1) * amp[622] + amp[623] - std::complex<double>
      (0, 1) * amp[624] - std::complex<double> (0, 1) * amp[629] -
      std::complex<double> (0, 1) * amp[630] - std::complex<double> (0, 1) *
      amp[631] - std::complex<double> (0, 1) * amp[633] - std::complex<double>
      (0, 1) * amp[634] + amp[855] + amp[861] - std::complex<double> (0, 1) *
      amp[864] - std::complex<double> (0, 1) * amp[865] - std::complex<double>
      (0, 1) * amp[867] - std::complex<double> (0, 1) * amp[868] + amp[870] +
      amp[871] - std::complex<double> (0, 1) * amp[876] - std::complex<double>
      (0, 1) * amp[877] - std::complex<double> (0, 1) * amp[879] -
      std::complex<double> (0, 1) * amp[880] + amp[882] + amp[883] + amp[900] +
      amp[901] + amp[903] + amp[904] + amp[906] + amp[907] + amp[909] +
      amp[910] + std::complex<double> (0, 1) * amp[1010] + std::complex<double>
      (0, 1) * amp[1009] + std::complex<double> (0, 1) * amp[1013] +
      std::complex<double> (0, 1) * amp[1012] - amp[1019] - amp[1018] +
      std::complex<double> (0, 1) * amp[1022] + std::complex<double> (0, 1) *
      amp[1021] + std::complex<double> (0, 1) * amp[1025] +
      std::complex<double> (0, 1) * amp[1024] - amp[1031] - amp[1030];
  jamp[3] = -std::complex<double> (0, 1) * amp[108] - std::complex<double> (0,
      1) * amp[109] + amp[110] + amp[111] + amp[115] - std::complex<double> (0,
      1) * amp[116] - std::complex<double> (0, 1) * amp[117] + amp[118] +
      amp[119] + amp[123] - amp[141] - amp[142] - std::complex<double> (0, 1) *
      amp[143] - std::complex<double> (0, 1) * amp[144] - std::complex<double>
      (0, 1) * amp[145] - amp[147] - amp[148] - std::complex<double> (0, 1) *
      amp[149] - std::complex<double> (0, 1) * amp[150] - std::complex<double>
      (0, 1) * amp[151] + std::complex<double> (0, 1) * amp[160] +
      std::complex<double> (0, 1) * amp[162] - amp[163] - std::complex<double>
      (0, 1) * amp[165] - std::complex<double> (0, 1) * amp[168] -
      std::complex<double> (0, 1) * amp[169] - std::complex<double> (0, 1) *
      amp[171] - std::complex<double> (0, 1) * amp[172] + std::complex<double>
      (0, 1) * amp[178] + std::complex<double> (0, 1) * amp[180] - amp[181] -
      std::complex<double> (0, 1) * amp[183] - std::complex<double> (0, 1) *
      amp[186] - std::complex<double> (0, 1) * amp[187] - std::complex<double>
      (0, 1) * amp[189] - std::complex<double> (0, 1) * amp[190] - amp[341] +
      std::complex<double> (0, 1) * amp[342] + std::complex<double> (0, 1) *
      amp[343] + std::complex<double> (0, 1) * amp[344] - amp[345] - amp[347] +
      std::complex<double> (0, 1) * amp[348] + std::complex<double> (0, 1) *
      amp[349] + std::complex<double> (0, 1) * amp[350] - amp[351] + amp[357] +
      amp[358] + amp[359] + amp[363] + amp[364] + amp[365] -
      std::complex<double> (0, 1) * amp[368] - std::complex<double> (0, 1) *
      amp[369] - std::complex<double> (0, 1) * amp[370] - std::complex<double>
      (0, 1) * amp[371] + amp[375] + amp[381] + amp[386] - std::complex<double>
      (0, 1) * amp[388] + amp[389] + amp[392] - std::complex<double> (0, 1) *
      amp[394] + amp[395] + amp[420] + amp[421] + amp[423] + amp[424] +
      amp[426] + amp[427] + amp[429] + amp[430] + std::complex<double> (0, 1) *
      amp[604] - amp[605] + std::complex<double> (0, 1) * amp[606] +
      std::complex<double> (0, 1) * amp[611] + std::complex<double> (0, 1) *
      amp[612] + std::complex<double> (0, 1) * amp[613] + std::complex<double>
      (0, 1) * amp[615] + std::complex<double> (0, 1) * amp[616] +
      std::complex<double> (0, 1) * amp[622] - amp[623] + std::complex<double>
      (0, 1) * amp[624] + std::complex<double> (0, 1) * amp[629] +
      std::complex<double> (0, 1) * amp[630] + std::complex<double> (0, 1) *
      amp[631] + std::complex<double> (0, 1) * amp[633] + std::complex<double>
      (0, 1) * amp[634] + amp[713] + amp[719] - amp[722] - amp[723] -
      std::complex<double> (0, 1) * amp[725] - amp[728] - amp[729] -
      std::complex<double> (0, 1) * amp[731] + std::complex<double> (0, 1) *
      amp[756] + std::complex<double> (0, 1) * amp[758] - amp[759] -
      std::complex<double> (0, 1) * amp[760] + amp[761] - std::complex<double>
      (0, 1) * amp[762] - std::complex<double> (0, 1) * amp[768] -
      std::complex<double> (0, 1) * amp[769] - std::complex<double> (0, 1) *
      amp[771] - std::complex<double> (0, 1) * amp[772] + std::complex<double>
      (0, 1) * amp[774] + std::complex<double> (0, 1) * amp[776] - amp[777] -
      std::complex<double> (0, 1) * amp[778] + amp[779] - std::complex<double>
      (0, 1) * amp[780] - std::complex<double> (0, 1) * amp[786] -
      std::complex<double> (0, 1) * amp[787] - std::complex<double> (0, 1) *
      amp[789] - std::complex<double> (0, 1) * amp[790] - amp[972] + amp[974] -
      amp[975] + amp[977] - amp[978] + amp[980] - amp[981] + amp[983] -
      std::complex<double> (0, 1) * amp[986] + std::complex<double> (0, 1) *
      amp[984] - std::complex<double> (0, 1) * amp[989] + std::complex<double>
      (0, 1) * amp[987] + amp[992] - amp[990] - std::complex<double> (0, 1) *
      amp[998] + std::complex<double> (0, 1) * amp[996] - std::complex<double>
      (0, 1) * amp[1001] + std::complex<double> (0, 1) * amp[999] + amp[1004] -
      amp[1002] - std::complex<double> (0, 1) * amp[1008] -
      std::complex<double> (0, 1) * amp[1009] - std::complex<double> (0, 1) *
      amp[1011] - std::complex<double> (0, 1) * amp[1012] + amp[1017] +
      amp[1018] - std::complex<double> (0, 1) * amp[1020] -
      std::complex<double> (0, 1) * amp[1021] - std::complex<double> (0, 1) *
      amp[1023] - std::complex<double> (0, 1) * amp[1024] + amp[1029] +
      amp[1030];
  jamp[4] = +std::complex<double> (0, 1) * amp[216] + std::complex<double> (0,
      1) * amp[217] - amp[218] - amp[219] - amp[223] + std::complex<double> (0,
      1) * amp[224] + std::complex<double> (0, 1) * amp[225] - amp[226] -
      amp[227] - amp[231] + std::complex<double> (0, 1) * amp[264] +
      std::complex<double> (0, 1) * amp[266] - amp[267] - std::complex<double>
      (0, 1) * amp[272] - std::complex<double> (0, 1) * amp[278] +
      std::complex<double> (0, 1) * amp[276] - std::complex<double> (0, 1) *
      amp[281] + std::complex<double> (0, 1) * amp[279] + std::complex<double>
      (0, 1) * amp[282] + std::complex<double> (0, 1) * amp[284] - amp[285] -
      std::complex<double> (0, 1) * amp[290] - std::complex<double> (0, 1) *
      amp[296] + std::complex<double> (0, 1) * amp[294] - std::complex<double>
      (0, 1) * amp[299] + std::complex<double> (0, 1) * amp[297] - amp[312] -
      amp[313] - std::complex<double> (0, 1) * amp[315] - std::complex<double>
      (0, 1) * amp[316] - std::complex<double> (0, 1) * amp[317] - amp[318] -
      amp[319] - std::complex<double> (0, 1) * amp[321] - std::complex<double>
      (0, 1) * amp[322] - std::complex<double> (0, 1) * amp[323] + amp[325] -
      std::complex<double> (0, 1) * amp[326] - std::complex<double> (0, 1) *
      amp[327] - std::complex<double> (0, 1) * amp[328] + amp[329] + amp[331] -
      std::complex<double> (0, 1) * amp[332] - std::complex<double> (0, 1) *
      amp[333] - std::complex<double> (0, 1) * amp[334] + amp[335] +
      std::complex<double> (0, 1) * amp[340] + amp[341] + amp[345] +
      std::complex<double> (0, 1) * amp[346] + amp[347] + amp[351] +
      std::complex<double> (0, 1) * amp[352] + std::complex<double> (0, 1) *
      amp[353] + std::complex<double> (0, 1) * amp[354] + std::complex<double>
      (0, 1) * amp[355] + amp[376] + amp[382] + amp[398] + amp[400] + amp[401]
      + amp[404] + amp[406] + amp[407] - amp[421] - amp[422] - amp[425] -
      amp[424] - amp[427] - amp[428] - amp[431] - amp[430] -
      std::complex<double> (0, 1) * amp[544] + amp[545] - std::complex<double>
      (0, 1) * amp[546] - std::complex<double> (0, 1) * amp[551] -
      std::complex<double> (0, 1) * amp[552] - std::complex<double> (0, 1) *
      amp[553] - std::complex<double> (0, 1) * amp[555] - std::complex<double>
      (0, 1) * amp[556] - std::complex<double> (0, 1) * amp[562] + amp[563] -
      std::complex<double> (0, 1) * amp[564] - std::complex<double> (0, 1) *
      amp[569] - std::complex<double> (0, 1) * amp[570] - std::complex<double>
      (0, 1) * amp[571] - std::complex<double> (0, 1) * amp[573] -
      std::complex<double> (0, 1) * amp[574] - std::complex<double> (0, 1) *
      amp[600] - std::complex<double> (0, 1) * amp[602] + amp[603] -
      std::complex<double> (0, 1) * amp[604] + amp[605] - std::complex<double>
      (0, 1) * amp[606] + std::complex<double> (0, 1) * amp[614] -
      std::complex<double> (0, 1) * amp[612] + std::complex<double> (0, 1) *
      amp[617] - std::complex<double> (0, 1) * amp[615] - std::complex<double>
      (0, 1) * amp[618] - std::complex<double> (0, 1) * amp[620] + amp[621] -
      std::complex<double> (0, 1) * amp[622] + amp[623] - std::complex<double>
      (0, 1) * amp[624] + std::complex<double> (0, 1) * amp[632] -
      std::complex<double> (0, 1) * amp[630] + std::complex<double> (0, 1) *
      amp[635] - std::complex<double> (0, 1) * amp[633] + amp[648] + amp[649] +
      std::complex<double> (0, 1) * amp[652] + amp[654] + amp[655] +
      std::complex<double> (0, 1) * amp[658] + amp[819] + amp[825] -
      std::complex<double> (0, 1) * amp[912] - std::complex<double> (0, 1) *
      amp[913] - std::complex<double> (0, 1) * amp[915] - std::complex<double>
      (0, 1) * amp[916] + amp[918] + amp[919] - std::complex<double> (0, 1) *
      amp[924] - std::complex<double> (0, 1) * amp[925] - std::complex<double>
      (0, 1) * amp[927] - std::complex<double> (0, 1) * amp[928] + amp[930] +
      amp[931] + amp[948] + amp[949] + amp[951] + amp[952] + amp[954] +
      amp[955] + amp[957] + amp[958] + std::complex<double> (0, 1) * amp[1010]
      + std::complex<double> (0, 1) * amp[1009] + std::complex<double> (0, 1) *
      amp[1013] + std::complex<double> (0, 1) * amp[1012] - amp[1019] -
      amp[1018] + std::complex<double> (0, 1) * amp[1022] +
      std::complex<double> (0, 1) * amp[1021] + std::complex<double> (0, 1) *
      amp[1025] + std::complex<double> (0, 1) * amp[1024] - amp[1031] -
      amp[1030];
  jamp[5] = -std::complex<double> (0, 1) * amp[216] - std::complex<double> (0,
      1) * amp[217] + amp[218] + amp[219] + amp[223] - std::complex<double> (0,
      1) * amp[224] - std::complex<double> (0, 1) * amp[225] + amp[226] +
      amp[227] + amp[231] - amp[249] - amp[250] - std::complex<double> (0, 1) *
      amp[251] - std::complex<double> (0, 1) * amp[252] - std::complex<double>
      (0, 1) * amp[253] - amp[255] - amp[256] - std::complex<double> (0, 1) *
      amp[257] - std::complex<double> (0, 1) * amp[258] - std::complex<double>
      (0, 1) * amp[259] + std::complex<double> (0, 1) * amp[268] +
      std::complex<double> (0, 1) * amp[270] - amp[271] - std::complex<double>
      (0, 1) * amp[273] - std::complex<double> (0, 1) * amp[276] -
      std::complex<double> (0, 1) * amp[277] - std::complex<double> (0, 1) *
      amp[279] - std::complex<double> (0, 1) * amp[280] + std::complex<double>
      (0, 1) * amp[286] + std::complex<double> (0, 1) * amp[288] - amp[289] -
      std::complex<double> (0, 1) * amp[291] - std::complex<double> (0, 1) *
      amp[294] - std::complex<double> (0, 1) * amp[295] - std::complex<double>
      (0, 1) * amp[297] - std::complex<double> (0, 1) * amp[298] - amp[325] +
      std::complex<double> (0, 1) * amp[326] + std::complex<double> (0, 1) *
      amp[327] + std::complex<double> (0, 1) * amp[328] - amp[329] - amp[331] +
      std::complex<double> (0, 1) * amp[332] + std::complex<double> (0, 1) *
      amp[333] + std::complex<double> (0, 1) * amp[334] - amp[335] + amp[356] +
      amp[360] + amp[361] + amp[362] + amp[366] + amp[367] +
      std::complex<double> (0, 1) * amp[368] + std::complex<double> (0, 1) *
      amp[369] + std::complex<double> (0, 1) * amp[370] + std::complex<double>
      (0, 1) * amp[371] + amp[377] + amp[383] - amp[386] + std::complex<double>
      (0, 1) * amp[388] - amp[389] - amp[392] + std::complex<double> (0, 1) *
      amp[394] - amp[395] - amp[420] + amp[422] + amp[425] - amp[423] -
      amp[426] + amp[428] + amp[431] - amp[429] + std::complex<double> (0, 1) *
      amp[544] - amp[545] + std::complex<double> (0, 1) * amp[546] +
      std::complex<double> (0, 1) * amp[551] + std::complex<double> (0, 1) *
      amp[552] + std::complex<double> (0, 1) * amp[553] + std::complex<double>
      (0, 1) * amp[555] + std::complex<double> (0, 1) * amp[556] +
      std::complex<double> (0, 1) * amp[562] - amp[563] + std::complex<double>
      (0, 1) * amp[564] + std::complex<double> (0, 1) * amp[569] +
      std::complex<double> (0, 1) * amp[570] + std::complex<double> (0, 1) *
      amp[571] + std::complex<double> (0, 1) * amp[573] + std::complex<double>
      (0, 1) * amp[574] + amp[711] + amp[717] + amp[722] + amp[723] +
      std::complex<double> (0, 1) * amp[725] + amp[728] + amp[729] +
      std::complex<double> (0, 1) * amp[731] - std::complex<double> (0, 1) *
      amp[756] - std::complex<double> (0, 1) * amp[758] + amp[759] +
      std::complex<double> (0, 1) * amp[760] - amp[761] + std::complex<double>
      (0, 1) * amp[762] + std::complex<double> (0, 1) * amp[768] +
      std::complex<double> (0, 1) * amp[769] + std::complex<double> (0, 1) *
      amp[771] + std::complex<double> (0, 1) * amp[772] - std::complex<double>
      (0, 1) * amp[774] - std::complex<double> (0, 1) * amp[776] + amp[777] +
      std::complex<double> (0, 1) * amp[778] - amp[779] + std::complex<double>
      (0, 1) * amp[780] + std::complex<double> (0, 1) * amp[786] +
      std::complex<double> (0, 1) * amp[787] + std::complex<double> (0, 1) *
      amp[789] + std::complex<double> (0, 1) * amp[790] + amp[972] + amp[973] +
      amp[975] + amp[976] + amp[978] + amp[979] + amp[981] + amp[982] -
      std::complex<double> (0, 1) * amp[984] - std::complex<double> (0, 1) *
      amp[985] - std::complex<double> (0, 1) * amp[987] - std::complex<double>
      (0, 1) * amp[988] + amp[990] + amp[991] - std::complex<double> (0, 1) *
      amp[996] - std::complex<double> (0, 1) * amp[997] - std::complex<double>
      (0, 1) * amp[999] - std::complex<double> (0, 1) * amp[1000] + amp[1002] +
      amp[1003] - std::complex<double> (0, 1) * amp[1010] +
      std::complex<double> (0, 1) * amp[1008] - std::complex<double> (0, 1) *
      amp[1013] + std::complex<double> (0, 1) * amp[1011] + amp[1019] -
      amp[1017] - std::complex<double> (0, 1) * amp[1022] +
      std::complex<double> (0, 1) * amp[1020] - std::complex<double> (0, 1) *
      amp[1025] + std::complex<double> (0, 1) * amp[1023] + amp[1031] -
      amp[1029];
  jamp[6] = -std::complex<double> (0, 1) * amp[0] - std::complex<double> (0, 1)
      * amp[2] + amp[3] + std::complex<double> (0, 1) * amp[8] +
      std::complex<double> (0, 1) * amp[14] - std::complex<double> (0, 1) *
      amp[12] + std::complex<double> (0, 1) * amp[17] - std::complex<double>
      (0, 1) * amp[15] - std::complex<double> (0, 1) * amp[18] -
      std::complex<double> (0, 1) * amp[20] + amp[21] + std::complex<double>
      (0, 1) * amp[26] + std::complex<double> (0, 1) * amp[32] -
      std::complex<double> (0, 1) * amp[30] + std::complex<double> (0, 1) *
      amp[35] - std::complex<double> (0, 1) * amp[33] - std::complex<double>
      (0, 1) * amp[36] - std::complex<double> (0, 1) * amp[37] + amp[38] +
      amp[39] + amp[42] - std::complex<double> (0, 1) * amp[44] -
      std::complex<double> (0, 1) * amp[45] + amp[46] + amp[47] + amp[50] +
      amp[96] + amp[97] + std::complex<double> (0, 1) * amp[99] +
      std::complex<double> (0, 1) * amp[100] + std::complex<double> (0, 1) *
      amp[101] + amp[102] + amp[103] + std::complex<double> (0, 1) * amp[105] +
      std::complex<double> (0, 1) * amp[106] + std::complex<double> (0, 1) *
      amp[107] - std::complex<double> (0, 1) * amp[124] - amp[125] - amp[126] -
      std::complex<double> (0, 1) * amp[130] - amp[131] - amp[132] -
      std::complex<double> (0, 1) * amp[136] - std::complex<double> (0, 1) *
      amp[137] - std::complex<double> (0, 1) * amp[138] - std::complex<double>
      (0, 1) * amp[139] - std::complex<double> (0, 1) * amp[156] -
      std::complex<double> (0, 1) * amp[158] + amp[159] - std::complex<double>
      (0, 1) * amp[160] - amp[161] - std::complex<double> (0, 1) * amp[162] +
      std::complex<double> (0, 1) * amp[170] + std::complex<double> (0, 1) *
      amp[169] + std::complex<double> (0, 1) * amp[173] + std::complex<double>
      (0, 1) * amp[172] - std::complex<double> (0, 1) * amp[174] -
      std::complex<double> (0, 1) * amp[176] + amp[177] - std::complex<double>
      (0, 1) * amp[178] - amp[179] - std::complex<double> (0, 1) * amp[180] +
      std::complex<double> (0, 1) * amp[188] + std::complex<double> (0, 1) *
      amp[187] + std::complex<double> (0, 1) * amp[191] + std::complex<double>
      (0, 1) * amp[190] + amp[204] + amp[205] - std::complex<double> (0, 1) *
      amp[206] + amp[210] + amp[211] - std::complex<double> (0, 1) * amp[212] +
      amp[660] + amp[666] - std::complex<double> (0, 1) * amp[672] -
      std::complex<double> (0, 1) * amp[673] - amp[674] - amp[675] -
      std::complex<double> (0, 1) * amp[676] - std::complex<double> (0, 1) *
      amp[678] - std::complex<double> (0, 1) * amp[679] - amp[680] - amp[681] -
      std::complex<double> (0, 1) * amp[682] + amp[696] + amp[697] + amp[700] +
      amp[702] + amp[703] + amp[706] - std::complex<double> (0, 1) * amp[756] -
      amp[757] - std::complex<double> (0, 1) * amp[758] - std::complex<double>
      (0, 1) * amp[765] + std::complex<double> (0, 1) * amp[770] +
      std::complex<double> (0, 1) * amp[769] + std::complex<double> (0, 1) *
      amp[773] + std::complex<double> (0, 1) * amp[772] - std::complex<double>
      (0, 1) * amp[774] - amp[775] - std::complex<double> (0, 1) * amp[776] -
      std::complex<double> (0, 1) * amp[783] + std::complex<double> (0, 1) *
      amp[788] + std::complex<double> (0, 1) * amp[787] + std::complex<double>
      (0, 1) * amp[791] + std::complex<double> (0, 1) * amp[790] + amp[857] +
      amp[863] + std::complex<double> (0, 1) * amp[866] + std::complex<double>
      (0, 1) * amp[865] + std::complex<double> (0, 1) * amp[869] +
      std::complex<double> (0, 1) * amp[868] - amp[872] - amp[871] +
      std::complex<double> (0, 1) * amp[878] + std::complex<double> (0, 1) *
      amp[877] + std::complex<double> (0, 1) * amp[881] + std::complex<double>
      (0, 1) * amp[880] - amp[884] - amp[883] - amp[901] - amp[902] - amp[904]
      - amp[905] - amp[907] - amp[908] - amp[910] - amp[911] - amp[960] +
      amp[962] + amp[965] - amp[963] - amp[966] + amp[968] + amp[971] -
      amp[969] + std::complex<double> (0, 1) * amp[986] - std::complex<double>
      (0, 1) * amp[984] + std::complex<double> (0, 1) * amp[989] -
      std::complex<double> (0, 1) * amp[987] + amp[995] - amp[993] +
      std::complex<double> (0, 1) * amp[998] - std::complex<double> (0, 1) *
      amp[996] + std::complex<double> (0, 1) * amp[1001] - std::complex<double>
      (0, 1) * amp[999] + amp[1007] - amp[1005];
  jamp[7] = -std::complex<double> (0, 1) * amp[4] - std::complex<double> (0, 1)
      * amp[6] + amp[7] + std::complex<double> (0, 1) * amp[9] +
      std::complex<double> (0, 1) * amp[12] + std::complex<double> (0, 1) *
      amp[13] + std::complex<double> (0, 1) * amp[15] + std::complex<double>
      (0, 1) * amp[16] - std::complex<double> (0, 1) * amp[22] -
      std::complex<double> (0, 1) * amp[24] + amp[25] + std::complex<double>
      (0, 1) * amp[27] + std::complex<double> (0, 1) * amp[30] +
      std::complex<double> (0, 1) * amp[31] + std::complex<double> (0, 1) *
      amp[33] + std::complex<double> (0, 1) * amp[34] + std::complex<double>
      (0, 1) * amp[36] + std::complex<double> (0, 1) * amp[37] - amp[38] -
      amp[39] - amp[42] + std::complex<double> (0, 1) * amp[44] +
      std::complex<double> (0, 1) * amp[45] - amp[46] - amp[47] - amp[50] +
      amp[69] + amp[70] + std::complex<double> (0, 1) * amp[71] +
      std::complex<double> (0, 1) * amp[72] + std::complex<double> (0, 1) *
      amp[73] + amp[75] + amp[76] + std::complex<double> (0, 1) * amp[77] +
      std::complex<double> (0, 1) * amp[78] + std::complex<double> (0, 1) *
      amp[79] - std::complex<double> (0, 1) * amp[232] - amp[233] - amp[234] -
      std::complex<double> (0, 1) * amp[238] - amp[239] - amp[240] -
      std::complex<double> (0, 1) * amp[244] - std::complex<double> (0, 1) *
      amp[245] - std::complex<double> (0, 1) * amp[246] - std::complex<double>
      (0, 1) * amp[247] - std::complex<double> (0, 1) * amp[264] -
      std::complex<double> (0, 1) * amp[266] + amp[267] - std::complex<double>
      (0, 1) * amp[268] - amp[269] - std::complex<double> (0, 1) * amp[270] +
      std::complex<double> (0, 1) * amp[278] + std::complex<double> (0, 1) *
      amp[277] + std::complex<double> (0, 1) * amp[281] + std::complex<double>
      (0, 1) * amp[280] - std::complex<double> (0, 1) * amp[282] -
      std::complex<double> (0, 1) * amp[284] + amp[285] - std::complex<double>
      (0, 1) * amp[286] - amp[287] - std::complex<double> (0, 1) * amp[288] +
      std::complex<double> (0, 1) * amp[296] + std::complex<double> (0, 1) *
      amp[295] + std::complex<double> (0, 1) * amp[299] + std::complex<double>
      (0, 1) * amp[298] + amp[312] + amp[313] - std::complex<double> (0, 1) *
      amp[314] + amp[318] + amp[319] - std::complex<double> (0, 1) * amp[320] +
      amp[661] + amp[667] + std::complex<double> (0, 1) * amp[672] +
      std::complex<double> (0, 1) * amp[673] + amp[674] + amp[675] +
      std::complex<double> (0, 1) * amp[676] + std::complex<double> (0, 1) *
      amp[678] + std::complex<double> (0, 1) * amp[679] + amp[680] + amp[681] +
      std::complex<double> (0, 1) * amp[682] + amp[684] + amp[685] + amp[688] +
      amp[690] + amp[691] + amp[694] + std::complex<double> (0, 1) * amp[756] +
      amp[757] + std::complex<double> (0, 1) * amp[758] + std::complex<double>
      (0, 1) * amp[765] - std::complex<double> (0, 1) * amp[770] -
      std::complex<double> (0, 1) * amp[769] - std::complex<double> (0, 1) *
      amp[773] - std::complex<double> (0, 1) * amp[772] + std::complex<double>
      (0, 1) * amp[774] + amp[775] + std::complex<double> (0, 1) * amp[776] +
      std::complex<double> (0, 1) * amp[783] - std::complex<double> (0, 1) *
      amp[788] - std::complex<double> (0, 1) * amp[787] - std::complex<double>
      (0, 1) * amp[791] - std::complex<double> (0, 1) * amp[790] + amp[821] +
      amp[827] + std::complex<double> (0, 1) * amp[914] + std::complex<double>
      (0, 1) * amp[913] + std::complex<double> (0, 1) * amp[917] +
      std::complex<double> (0, 1) * amp[916] - amp[920] - amp[919] +
      std::complex<double> (0, 1) * amp[926] + std::complex<double> (0, 1) *
      amp[925] + std::complex<double> (0, 1) * amp[929] + std::complex<double>
      (0, 1) * amp[928] - amp[932] - amp[931] - amp[949] - amp[950] - amp[952]
      - amp[953] - amp[955] - amp[956] - amp[958] - amp[959] + amp[960] +
      amp[961] + amp[963] + amp[964] + amp[966] + amp[967] + amp[969] +
      amp[970] + std::complex<double> (0, 1) * amp[984] + std::complex<double>
      (0, 1) * amp[985] + std::complex<double> (0, 1) * amp[987] +
      std::complex<double> (0, 1) * amp[988] + amp[993] + amp[994] +
      std::complex<double> (0, 1) * amp[996] + std::complex<double> (0, 1) *
      amp[997] + std::complex<double> (0, 1) * amp[999] + std::complex<double>
      (0, 1) * amp[1000] + amp[1005] + amp[1006];
  jamp[8] = +std::complex<double> (0, 1) * amp[124] + amp[125] + amp[126] +
      std::complex<double> (0, 1) * amp[130] + amp[131] + amp[132] +
      std::complex<double> (0, 1) * amp[136] + std::complex<double> (0, 1) *
      amp[137] + std::complex<double> (0, 1) * amp[138] + std::complex<double>
      (0, 1) * amp[139] + std::complex<double> (0, 1) * amp[156] +
      std::complex<double> (0, 1) * amp[158] - amp[159] + std::complex<double>
      (0, 1) * amp[160] + amp[161] + std::complex<double> (0, 1) * amp[162] -
      std::complex<double> (0, 1) * amp[170] - std::complex<double> (0, 1) *
      amp[169] - std::complex<double> (0, 1) * amp[173] - std::complex<double>
      (0, 1) * amp[172] + std::complex<double> (0, 1) * amp[174] +
      std::complex<double> (0, 1) * amp[176] - amp[177] + std::complex<double>
      (0, 1) * amp[178] + amp[179] + std::complex<double> (0, 1) * amp[180] -
      std::complex<double> (0, 1) * amp[188] - std::complex<double> (0, 1) *
      amp[187] - std::complex<double> (0, 1) * amp[191] - std::complex<double>
      (0, 1) * amp[190] - amp[204] - amp[205] + std::complex<double> (0, 1) *
      amp[206] - amp[210] - amp[211] + std::complex<double> (0, 1) * amp[212] -
      std::complex<double> (0, 1) * amp[216] - std::complex<double> (0, 1) *
      amp[217] - amp[220] - amp[221] - amp[222] - std::complex<double> (0, 1) *
      amp[224] - std::complex<double> (0, 1) * amp[225] - amp[228] - amp[229] -
      amp[230] + amp[233] + amp[234] - std::complex<double> (0, 1) * amp[235] -
      std::complex<double> (0, 1) * amp[236] - std::complex<double> (0, 1) *
      amp[237] + amp[239] + amp[240] - std::complex<double> (0, 1) * amp[241] -
      std::complex<double> (0, 1) * amp[242] - std::complex<double> (0, 1) *
      amp[243] + std::complex<double> (0, 1) * amp[268] + amp[269] +
      std::complex<double> (0, 1) * amp[270] - std::complex<double> (0, 1) *
      amp[275] - std::complex<double> (0, 1) * amp[276] - std::complex<double>
      (0, 1) * amp[277] - std::complex<double> (0, 1) * amp[279] -
      std::complex<double> (0, 1) * amp[280] + std::complex<double> (0, 1) *
      amp[286] + amp[287] + std::complex<double> (0, 1) * amp[288] -
      std::complex<double> (0, 1) * amp[293] - std::complex<double> (0, 1) *
      amp[294] - std::complex<double> (0, 1) * amp[295] - std::complex<double>
      (0, 1) * amp[297] - std::complex<double> (0, 1) * amp[298] -
      std::complex<double> (0, 1) * amp[540] - std::complex<double> (0, 1) *
      amp[542] + amp[543] - std::complex<double> (0, 1) * amp[548] +
      std::complex<double> (0, 1) * amp[554] + std::complex<double> (0, 1) *
      amp[553] + std::complex<double> (0, 1) * amp[557] + std::complex<double>
      (0, 1) * amp[556] - std::complex<double> (0, 1) * amp[558] -
      std::complex<double> (0, 1) * amp[560] + amp[561] - std::complex<double>
      (0, 1) * amp[566] + std::complex<double> (0, 1) * amp[572] +
      std::complex<double> (0, 1) * amp[571] + std::complex<double> (0, 1) *
      amp[575] + std::complex<double> (0, 1) * amp[574] + amp[588] + amp[589] -
      std::complex<double> (0, 1) * amp[590] - std::complex<double> (0, 1) *
      amp[591] - std::complex<double> (0, 1) * amp[593] + amp[594] + amp[595] -
      std::complex<double> (0, 1) * amp[596] - std::complex<double> (0, 1) *
      amp[597] - std::complex<double> (0, 1) * amp[599] + amp[662] + amp[668] +
      amp[698] + amp[699] + amp[701] + amp[704] + amp[705] + amp[707] +
      amp[853] + amp[859] - std::complex<double> (0, 1) * amp[864] -
      std::complex<double> (0, 1) * amp[865] - std::complex<double> (0, 1) *
      amp[867] - std::complex<double> (0, 1) * amp[868] + amp[870] + amp[871] -
      std::complex<double> (0, 1) * amp[876] - std::complex<double> (0, 1) *
      amp[877] - std::complex<double> (0, 1) * amp[879] - std::complex<double>
      (0, 1) * amp[880] + amp[882] + amp[883] + amp[900] + amp[901] + amp[903]
      + amp[904] + amp[906] + amp[907] + amp[909] + amp[910] - amp[961] -
      amp[962] - amp[965] - amp[964] - amp[967] - amp[968] - amp[971] -
      amp[970] - std::complex<double> (0, 1) * amp[986] - std::complex<double>
      (0, 1) * amp[985] - std::complex<double> (0, 1) * amp[989] -
      std::complex<double> (0, 1) * amp[988] - amp[995] - amp[994] -
      std::complex<double> (0, 1) * amp[998] - std::complex<double> (0, 1) *
      amp[997] - std::complex<double> (0, 1) * amp[1001] - std::complex<double>
      (0, 1) * amp[1000] - amp[1007] - amp[1006];
  jamp[9] = +std::complex<double> (0, 1) * amp[216] + std::complex<double> (0,
      1) * amp[217] + amp[220] + amp[221] + amp[222] + std::complex<double> (0,
      1) * amp[224] + std::complex<double> (0, 1) * amp[225] + amp[228] +
      amp[229] + amp[230] - amp[233] - amp[234] + std::complex<double> (0, 1) *
      amp[235] + std::complex<double> (0, 1) * amp[236] + std::complex<double>
      (0, 1) * amp[237] - amp[239] - amp[240] + std::complex<double> (0, 1) *
      amp[241] + std::complex<double> (0, 1) * amp[242] + std::complex<double>
      (0, 1) * amp[243] - std::complex<double> (0, 1) * amp[268] - amp[269] -
      std::complex<double> (0, 1) * amp[270] + std::complex<double> (0, 1) *
      amp[275] + std::complex<double> (0, 1) * amp[276] + std::complex<double>
      (0, 1) * amp[277] + std::complex<double> (0, 1) * amp[279] +
      std::complex<double> (0, 1) * amp[280] - std::complex<double> (0, 1) *
      amp[286] - amp[287] - std::complex<double> (0, 1) * amp[288] +
      std::complex<double> (0, 1) * amp[293] + std::complex<double> (0, 1) *
      amp[294] + std::complex<double> (0, 1) * amp[295] + std::complex<double>
      (0, 1) * amp[297] + std::complex<double> (0, 1) * amp[298] - amp[433] -
      std::complex<double> (0, 1) * amp[434] - std::complex<double> (0, 1) *
      amp[435] - std::complex<double> (0, 1) * amp[436] - amp[437] - amp[439] -
      std::complex<double> (0, 1) * amp[440] - std::complex<double> (0, 1) *
      amp[441] - std::complex<double> (0, 1) * amp[442] - amp[443] + amp[464] +
      amp[468] + amp[469] + amp[470] + amp[474] + amp[475] -
      std::complex<double> (0, 1) * amp[476] - std::complex<double> (0, 1) *
      amp[477] - std::complex<double> (0, 1) * amp[478] - std::complex<double>
      (0, 1) * amp[479] + amp[485] + amp[491] - amp[494] - std::complex<double>
      (0, 1) * amp[496] - amp[497] - amp[500] - std::complex<double> (0, 1) *
      amp[502] - amp[503] - amp[528] + amp[530] + amp[533] - amp[531] -
      amp[534] + amp[536] + amp[539] - amp[537] - std::complex<double> (0, 1) *
      amp[544] - std::complex<double> (0, 1) * amp[546] - amp[547] -
      std::complex<double> (0, 1) * amp[549] - std::complex<double> (0, 1) *
      amp[552] - std::complex<double> (0, 1) * amp[553] - std::complex<double>
      (0, 1) * amp[555] - std::complex<double> (0, 1) * amp[556] -
      std::complex<double> (0, 1) * amp[562] - std::complex<double> (0, 1) *
      amp[564] - amp[565] - std::complex<double> (0, 1) * amp[567] -
      std::complex<double> (0, 1) * amp[570] - std::complex<double> (0, 1) *
      amp[571] - std::complex<double> (0, 1) * amp[573] - std::complex<double>
      (0, 1) * amp[574] + amp[663] + amp[669] + amp[674] + amp[675] -
      std::complex<double> (0, 1) * amp[677] + amp[680] + amp[681] -
      std::complex<double> (0, 1) * amp[683] + std::complex<double> (0, 1) *
      amp[756] + amp[757] + std::complex<double> (0, 1) * amp[758] -
      std::complex<double> (0, 1) * amp[760] - std::complex<double> (0, 1) *
      amp[762] - amp[763] - std::complex<double> (0, 1) * amp[768] -
      std::complex<double> (0, 1) * amp[769] - std::complex<double> (0, 1) *
      amp[771] - std::complex<double> (0, 1) * amp[772] + std::complex<double>
      (0, 1) * amp[774] + amp[775] + std::complex<double> (0, 1) * amp[776] -
      std::complex<double> (0, 1) * amp[778] - std::complex<double> (0, 1) *
      amp[780] - amp[781] - std::complex<double> (0, 1) * amp[786] -
      std::complex<double> (0, 1) * amp[787] - std::complex<double> (0, 1) *
      amp[789] - std::complex<double> (0, 1) * amp[790] + amp[960] + amp[961] +
      amp[963] + amp[964] + amp[966] + amp[967] + amp[969] + amp[970] +
      std::complex<double> (0, 1) * amp[984] + std::complex<double> (0, 1) *
      amp[985] + std::complex<double> (0, 1) * amp[987] + std::complex<double>
      (0, 1) * amp[988] + amp[993] + amp[994] + std::complex<double> (0, 1) *
      amp[996] + std::complex<double> (0, 1) * amp[997] + std::complex<double>
      (0, 1) * amp[999] + std::complex<double> (0, 1) * amp[1000] + amp[1005] +
      amp[1006] + std::complex<double> (0, 1) * amp[1010] -
      std::complex<double> (0, 1) * amp[1008] + std::complex<double> (0, 1) *
      amp[1013] - std::complex<double> (0, 1) * amp[1011] + amp[1016] -
      amp[1014] + std::complex<double> (0, 1) * amp[1022] -
      std::complex<double> (0, 1) * amp[1020] + std::complex<double> (0, 1) *
      amp[1025] - std::complex<double> (0, 1) * amp[1023] + amp[1028] -
      amp[1026];
  jamp[10] = -std::complex<double> (0, 1) * amp[108] - std::complex<double> (0,
      1) * amp[109] - amp[112] - amp[113] - amp[114] - std::complex<double> (0,
      1) * amp[116] - std::complex<double> (0, 1) * amp[117] - amp[120] -
      amp[121] - amp[122] + amp[125] + amp[126] - std::complex<double> (0, 1) *
      amp[127] - std::complex<double> (0, 1) * amp[128] - std::complex<double>
      (0, 1) * amp[129] + amp[131] + amp[132] - std::complex<double> (0, 1) *
      amp[133] - std::complex<double> (0, 1) * amp[134] - std::complex<double>
      (0, 1) * amp[135] + std::complex<double> (0, 1) * amp[160] + amp[161] +
      std::complex<double> (0, 1) * amp[162] - std::complex<double> (0, 1) *
      amp[167] - std::complex<double> (0, 1) * amp[168] - std::complex<double>
      (0, 1) * amp[169] - std::complex<double> (0, 1) * amp[171] -
      std::complex<double> (0, 1) * amp[172] + std::complex<double> (0, 1) *
      amp[178] + amp[179] + std::complex<double> (0, 1) * amp[180] -
      std::complex<double> (0, 1) * amp[185] - std::complex<double> (0, 1) *
      amp[186] - std::complex<double> (0, 1) * amp[187] - std::complex<double>
      (0, 1) * amp[189] - std::complex<double> (0, 1) * amp[190] +
      std::complex<double> (0, 1) * amp[232] + amp[233] + amp[234] +
      std::complex<double> (0, 1) * amp[238] + amp[239] + amp[240] +
      std::complex<double> (0, 1) * amp[244] + std::complex<double> (0, 1) *
      amp[245] + std::complex<double> (0, 1) * amp[246] + std::complex<double>
      (0, 1) * amp[247] + std::complex<double> (0, 1) * amp[264] +
      std::complex<double> (0, 1) * amp[266] - amp[267] + std::complex<double>
      (0, 1) * amp[268] + amp[269] + std::complex<double> (0, 1) * amp[270] -
      std::complex<double> (0, 1) * amp[278] - std::complex<double> (0, 1) *
      amp[277] - std::complex<double> (0, 1) * amp[281] - std::complex<double>
      (0, 1) * amp[280] + std::complex<double> (0, 1) * amp[282] +
      std::complex<double> (0, 1) * amp[284] - amp[285] + std::complex<double>
      (0, 1) * amp[286] + amp[287] + std::complex<double> (0, 1) * amp[288] -
      std::complex<double> (0, 1) * amp[296] - std::complex<double> (0, 1) *
      amp[295] - std::complex<double> (0, 1) * amp[299] - std::complex<double>
      (0, 1) * amp[298] - amp[312] - amp[313] + std::complex<double> (0, 1) *
      amp[314] - amp[318] - amp[319] + std::complex<double> (0, 1) * amp[320] -
      std::complex<double> (0, 1) * amp[600] - std::complex<double> (0, 1) *
      amp[602] + amp[603] - std::complex<double> (0, 1) * amp[608] +
      std::complex<double> (0, 1) * amp[614] + std::complex<double> (0, 1) *
      amp[613] + std::complex<double> (0, 1) * amp[617] + std::complex<double>
      (0, 1) * amp[616] - std::complex<double> (0, 1) * amp[618] -
      std::complex<double> (0, 1) * amp[620] + amp[621] - std::complex<double>
      (0, 1) * amp[626] + std::complex<double> (0, 1) * amp[632] +
      std::complex<double> (0, 1) * amp[631] + std::complex<double> (0, 1) *
      amp[635] + std::complex<double> (0, 1) * amp[634] + amp[648] + amp[649] -
      std::complex<double> (0, 1) * amp[650] - std::complex<double> (0, 1) *
      amp[651] - std::complex<double> (0, 1) * amp[653] + amp[654] + amp[655] -
      std::complex<double> (0, 1) * amp[656] - std::complex<double> (0, 1) *
      amp[657] - std::complex<double> (0, 1) * amp[659] + amp[664] + amp[670] +
      amp[686] + amp[687] + amp[689] + amp[692] + amp[693] + amp[695] +
      amp[817] + amp[823] - std::complex<double> (0, 1) * amp[912] -
      std::complex<double> (0, 1) * amp[913] - std::complex<double> (0, 1) *
      amp[915] - std::complex<double> (0, 1) * amp[916] + amp[918] + amp[919] -
      std::complex<double> (0, 1) * amp[924] - std::complex<double> (0, 1) *
      amp[925] - std::complex<double> (0, 1) * amp[927] - std::complex<double>
      (0, 1) * amp[928] + amp[930] + amp[931] + amp[948] + amp[949] + amp[951]
      + amp[952] + amp[954] + amp[955] + amp[957] + amp[958] - amp[961] -
      amp[962] - amp[965] - amp[964] - amp[967] - amp[968] - amp[971] -
      amp[970] - std::complex<double> (0, 1) * amp[986] - std::complex<double>
      (0, 1) * amp[985] - std::complex<double> (0, 1) * amp[989] -
      std::complex<double> (0, 1) * amp[988] - amp[995] - amp[994] -
      std::complex<double> (0, 1) * amp[998] - std::complex<double> (0, 1) *
      amp[997] - std::complex<double> (0, 1) * amp[1001] - std::complex<double>
      (0, 1) * amp[1000] - amp[1007] - amp[1006];
  jamp[11] = +std::complex<double> (0, 1) * amp[108] + std::complex<double> (0,
      1) * amp[109] + amp[112] + amp[113] + amp[114] + std::complex<double> (0,
      1) * amp[116] + std::complex<double> (0, 1) * amp[117] + amp[120] +
      amp[121] + amp[122] - amp[125] - amp[126] + std::complex<double> (0, 1) *
      amp[127] + std::complex<double> (0, 1) * amp[128] + std::complex<double>
      (0, 1) * amp[129] - amp[131] - amp[132] + std::complex<double> (0, 1) *
      amp[133] + std::complex<double> (0, 1) * amp[134] + std::complex<double>
      (0, 1) * amp[135] - std::complex<double> (0, 1) * amp[160] - amp[161] -
      std::complex<double> (0, 1) * amp[162] + std::complex<double> (0, 1) *
      amp[167] + std::complex<double> (0, 1) * amp[168] + std::complex<double>
      (0, 1) * amp[169] + std::complex<double> (0, 1) * amp[171] +
      std::complex<double> (0, 1) * amp[172] - std::complex<double> (0, 1) *
      amp[178] - amp[179] - std::complex<double> (0, 1) * amp[180] +
      std::complex<double> (0, 1) * amp[185] + std::complex<double> (0, 1) *
      amp[186] + std::complex<double> (0, 1) * amp[187] + std::complex<double>
      (0, 1) * amp[189] + std::complex<double> (0, 1) * amp[190] - amp[449] -
      std::complex<double> (0, 1) * amp[450] - std::complex<double> (0, 1) *
      amp[451] - std::complex<double> (0, 1) * amp[452] - amp[453] - amp[455] -
      std::complex<double> (0, 1) * amp[456] - std::complex<double> (0, 1) *
      amp[457] - std::complex<double> (0, 1) * amp[458] - amp[459] + amp[465] +
      amp[466] + amp[467] + amp[471] + amp[472] + amp[473] +
      std::complex<double> (0, 1) * amp[476] + std::complex<double> (0, 1) *
      amp[477] + std::complex<double> (0, 1) * amp[478] + std::complex<double>
      (0, 1) * amp[479] + amp[483] + amp[489] + amp[494] + std::complex<double>
      (0, 1) * amp[496] + amp[497] + amp[500] + std::complex<double> (0, 1) *
      amp[502] + amp[503] + amp[528] + amp[529] + amp[531] + amp[532] +
      amp[534] + amp[535] + amp[537] + amp[538] - std::complex<double> (0, 1) *
      amp[604] - std::complex<double> (0, 1) * amp[606] - amp[607] -
      std::complex<double> (0, 1) * amp[609] - std::complex<double> (0, 1) *
      amp[612] - std::complex<double> (0, 1) * amp[613] - std::complex<double>
      (0, 1) * amp[615] - std::complex<double> (0, 1) * amp[616] -
      std::complex<double> (0, 1) * amp[622] - std::complex<double> (0, 1) *
      amp[624] - amp[625] - std::complex<double> (0, 1) * amp[627] -
      std::complex<double> (0, 1) * amp[630] - std::complex<double> (0, 1) *
      amp[631] - std::complex<double> (0, 1) * amp[633] - std::complex<double>
      (0, 1) * amp[634] + amp[665] + amp[671] - amp[674] - amp[675] +
      std::complex<double> (0, 1) * amp[677] - amp[680] - amp[681] +
      std::complex<double> (0, 1) * amp[683] - std::complex<double> (0, 1) *
      amp[756] - amp[757] - std::complex<double> (0, 1) * amp[758] +
      std::complex<double> (0, 1) * amp[760] + std::complex<double> (0, 1) *
      amp[762] + amp[763] + std::complex<double> (0, 1) * amp[768] +
      std::complex<double> (0, 1) * amp[769] + std::complex<double> (0, 1) *
      amp[771] + std::complex<double> (0, 1) * amp[772] - std::complex<double>
      (0, 1) * amp[774] - amp[775] - std::complex<double> (0, 1) * amp[776] +
      std::complex<double> (0, 1) * amp[778] + std::complex<double> (0, 1) *
      amp[780] + amp[781] + std::complex<double> (0, 1) * amp[786] +
      std::complex<double> (0, 1) * amp[787] + std::complex<double> (0, 1) *
      amp[789] + std::complex<double> (0, 1) * amp[790] - amp[960] + amp[962] +
      amp[965] - amp[963] - amp[966] + amp[968] + amp[971] - amp[969] +
      std::complex<double> (0, 1) * amp[986] - std::complex<double> (0, 1) *
      amp[984] + std::complex<double> (0, 1) * amp[989] - std::complex<double>
      (0, 1) * amp[987] + amp[995] - amp[993] + std::complex<double> (0, 1) *
      amp[998] - std::complex<double> (0, 1) * amp[996] + std::complex<double>
      (0, 1) * amp[1001] - std::complex<double> (0, 1) * amp[999] + amp[1007] -
      amp[1005] + std::complex<double> (0, 1) * amp[1008] +
      std::complex<double> (0, 1) * amp[1009] + std::complex<double> (0, 1) *
      amp[1011] + std::complex<double> (0, 1) * amp[1012] + amp[1014] +
      amp[1015] + std::complex<double> (0, 1) * amp[1020] +
      std::complex<double> (0, 1) * amp[1021] + std::complex<double> (0, 1) *
      amp[1023] + std::complex<double> (0, 1) * amp[1024] + amp[1026] +
      amp[1027];
  jamp[12] = -std::complex<double> (0, 1) * amp[0] - std::complex<double> (0,
      1) * amp[2] + amp[3] - std::complex<double> (0, 1) * amp[4] - amp[5] -
      std::complex<double> (0, 1) * amp[6] + std::complex<double> (0, 1) *
      amp[14] + std::complex<double> (0, 1) * amp[13] + std::complex<double>
      (0, 1) * amp[17] + std::complex<double> (0, 1) * amp[16] -
      std::complex<double> (0, 1) * amp[18] - std::complex<double> (0, 1) *
      amp[20] + amp[21] - std::complex<double> (0, 1) * amp[22] - amp[23] -
      std::complex<double> (0, 1) * amp[24] + std::complex<double> (0, 1) *
      amp[32] + std::complex<double> (0, 1) * amp[31] + std::complex<double>
      (0, 1) * amp[35] + std::complex<double> (0, 1) * amp[34] -
      std::complex<double> (0, 1) * amp[52] - amp[53] - amp[54] -
      std::complex<double> (0, 1) * amp[58] - amp[59] - amp[60] -
      std::complex<double> (0, 1) * amp[64] - std::complex<double> (0, 1) *
      amp[65] - std::complex<double> (0, 1) * amp[66] - std::complex<double>
      (0, 1) * amp[67] + amp[96] + amp[97] - std::complex<double> (0, 1) *
      amp[98] + amp[102] + amp[103] - std::complex<double> (0, 1) * amp[104] -
      std::complex<double> (0, 1) * amp[108] - std::complex<double> (0, 1) *
      amp[109] + amp[110] + amp[111] + amp[115] - std::complex<double> (0, 1) *
      amp[116] - std::complex<double> (0, 1) * amp[117] + amp[118] + amp[119] +
      amp[123] - std::complex<double> (0, 1) * amp[156] - std::complex<double>
      (0, 1) * amp[158] + amp[159] + std::complex<double> (0, 1) * amp[164] +
      std::complex<double> (0, 1) * amp[170] - std::complex<double> (0, 1) *
      amp[168] + std::complex<double> (0, 1) * amp[173] - std::complex<double>
      (0, 1) * amp[171] - std::complex<double> (0, 1) * amp[174] -
      std::complex<double> (0, 1) * amp[176] + amp[177] + std::complex<double>
      (0, 1) * amp[182] + std::complex<double> (0, 1) * amp[188] -
      std::complex<double> (0, 1) * amp[186] + std::complex<double> (0, 1) *
      amp[191] - std::complex<double> (0, 1) * amp[189] + amp[204] + amp[205] +
      std::complex<double> (0, 1) * amp[207] + std::complex<double> (0, 1) *
      amp[208] + std::complex<double> (0, 1) * amp[209] + amp[210] + amp[211] +
      std::complex<double> (0, 1) * amp[213] + std::complex<double> (0, 1) *
      amp[214] + std::complex<double> (0, 1) * amp[215] - std::complex<double>
      (0, 1) * amp[600] - amp[601] - std::complex<double> (0, 1) * amp[602] -
      std::complex<double> (0, 1) * amp[610] + std::complex<double> (0, 1) *
      amp[614] + std::complex<double> (0, 1) * amp[613] + std::complex<double>
      (0, 1) * amp[617] + std::complex<double> (0, 1) * amp[616] -
      std::complex<double> (0, 1) * amp[618] - amp[619] - std::complex<double>
      (0, 1) * amp[620] - std::complex<double> (0, 1) * amp[628] +
      std::complex<double> (0, 1) * amp[632] + std::complex<double> (0, 1) *
      amp[631] + std::complex<double> (0, 1) * amp[635] + std::complex<double>
      (0, 1) * amp[634] - amp[636] - amp[637] - std::complex<double> (0, 1) *
      amp[638] - std::complex<double> (0, 1) * amp[639] - std::complex<double>
      (0, 1) * amp[641] - amp[642] - amp[643] - std::complex<double> (0, 1) *
      amp[644] - std::complex<double> (0, 1) * amp[645] - std::complex<double>
      (0, 1) * amp[647] + amp[792] + amp[798] + amp[804] + amp[805] + amp[809]
      + amp[810] + amp[811] + amp[815] + amp[854] + amp[860] +
      std::complex<double> (0, 1) * amp[866] + std::complex<double> (0, 1) *
      amp[865] + std::complex<double> (0, 1) * amp[869] + std::complex<double>
      (0, 1) * amp[868] - amp[872] - amp[871] + std::complex<double> (0, 1) *
      amp[878] + std::complex<double> (0, 1) * amp[877] + std::complex<double>
      (0, 1) * amp[881] + std::complex<double> (0, 1) * amp[880] - amp[884] -
      amp[883] - amp[901] - amp[902] - amp[904] - amp[905] - amp[907] -
      amp[908] - amp[910] - amp[911] + std::complex<double> (0, 1) * amp[914] -
      std::complex<double> (0, 1) * amp[912] + std::complex<double> (0, 1) *
      amp[917] - std::complex<double> (0, 1) * amp[915] + amp[923] - amp[921] +
      std::complex<double> (0, 1) * amp[926] - std::complex<double> (0, 1) *
      amp[924] + std::complex<double> (0, 1) * amp[929] - std::complex<double>
      (0, 1) * amp[927] + amp[935] - amp[933] - amp[936] + amp[938] + amp[941]
      - amp[939] - amp[942] + amp[944] + amp[947] - amp[945];
  jamp[13] = +std::complex<double> (0, 1) * amp[108] + std::complex<double> (0,
      1) * amp[109] - amp[110] - amp[111] - amp[115] + std::complex<double> (0,
      1) * amp[116] + std::complex<double> (0, 1) * amp[117] - amp[118] -
      amp[119] - amp[123] + amp[141] + amp[142] + std::complex<double> (0, 1) *
      amp[143] + std::complex<double> (0, 1) * amp[144] + std::complex<double>
      (0, 1) * amp[145] + amp[147] + amp[148] + std::complex<double> (0, 1) *
      amp[149] + std::complex<double> (0, 1) * amp[150] + std::complex<double>
      (0, 1) * amp[151] - std::complex<double> (0, 1) * amp[160] -
      std::complex<double> (0, 1) * amp[162] + amp[163] + std::complex<double>
      (0, 1) * amp[165] + std::complex<double> (0, 1) * amp[168] +
      std::complex<double> (0, 1) * amp[169] + std::complex<double> (0, 1) *
      amp[171] + std::complex<double> (0, 1) * amp[172] - std::complex<double>
      (0, 1) * amp[178] - std::complex<double> (0, 1) * amp[180] + amp[181] +
      std::complex<double> (0, 1) * amp[183] + std::complex<double> (0, 1) *
      amp[186] + std::complex<double> (0, 1) * amp[187] + std::complex<double>
      (0, 1) * amp[189] + std::complex<double> (0, 1) * amp[190] -
      std::complex<double> (0, 1) * amp[248] + amp[249] + amp[250] -
      std::complex<double> (0, 1) * amp[254] + amp[255] + amp[256] -
      std::complex<double> (0, 1) * amp[260] - std::complex<double> (0, 1) *
      amp[261] - std::complex<double> (0, 1) * amp[262] - std::complex<double>
      (0, 1) * amp[263] - std::complex<double> (0, 1) * amp[264] - amp[265] -
      std::complex<double> (0, 1) * amp[266] - std::complex<double> (0, 1) *
      amp[268] - std::complex<double> (0, 1) * amp[270] + amp[271] +
      std::complex<double> (0, 1) * amp[278] + std::complex<double> (0, 1) *
      amp[277] + std::complex<double> (0, 1) * amp[281] + std::complex<double>
      (0, 1) * amp[280] - std::complex<double> (0, 1) * amp[282] - amp[283] -
      std::complex<double> (0, 1) * amp[284] - std::complex<double> (0, 1) *
      amp[286] - std::complex<double> (0, 1) * amp[288] + amp[289] +
      std::complex<double> (0, 1) * amp[296] + std::complex<double> (0, 1) *
      amp[295] + std::complex<double> (0, 1) * amp[299] + std::complex<double>
      (0, 1) * amp[298] - amp[300] - amp[301] - std::complex<double> (0, 1) *
      amp[302] - amp[306] - amp[307] - std::complex<double> (0, 1) * amp[308] +
      std::complex<double> (0, 1) * amp[600] + amp[601] + std::complex<double>
      (0, 1) * amp[602] + std::complex<double> (0, 1) * amp[610] -
      std::complex<double> (0, 1) * amp[614] - std::complex<double> (0, 1) *
      amp[613] - std::complex<double> (0, 1) * amp[617] - std::complex<double>
      (0, 1) * amp[616] + std::complex<double> (0, 1) * amp[618] + amp[619] +
      std::complex<double> (0, 1) * amp[620] + std::complex<double> (0, 1) *
      amp[628] - std::complex<double> (0, 1) * amp[632] - std::complex<double>
      (0, 1) * amp[631] - std::complex<double> (0, 1) * amp[635] -
      std::complex<double> (0, 1) * amp[634] + amp[636] + amp[637] +
      std::complex<double> (0, 1) * amp[638] + std::complex<double> (0, 1) *
      amp[639] + std::complex<double> (0, 1) * amp[641] + amp[642] + amp[643] +
      std::complex<double> (0, 1) * amp[644] + std::complex<double> (0, 1) *
      amp[645] + std::complex<double> (0, 1) * amp[647] + amp[712] + amp[718] +
      amp[734] + amp[735] + amp[737] + amp[740] + amp[741] + amp[743] +
      amp[793] + amp[799] + std::complex<double> (0, 1) * amp[912] +
      std::complex<double> (0, 1) * amp[913] + std::complex<double> (0, 1) *
      amp[915] + std::complex<double> (0, 1) * amp[916] + amp[921] + amp[922] +
      std::complex<double> (0, 1) * amp[924] + std::complex<double> (0, 1) *
      amp[925] + std::complex<double> (0, 1) * amp[927] + std::complex<double>
      (0, 1) * amp[928] + amp[933] + amp[934] + amp[936] + amp[937] + amp[939]
      + amp[940] + amp[942] + amp[943] + amp[945] + amp[946] - amp[973] -
      amp[974] - amp[976] - amp[977] - amp[979] - amp[980] - amp[982] -
      amp[983] + std::complex<double> (0, 1) * amp[986] + std::complex<double>
      (0, 1) * amp[985] + std::complex<double> (0, 1) * amp[989] +
      std::complex<double> (0, 1) * amp[988] - amp[992] - amp[991] +
      std::complex<double> (0, 1) * amp[998] + std::complex<double> (0, 1) *
      amp[997] + std::complex<double> (0, 1) * amp[1001] + std::complex<double>
      (0, 1) * amp[1000] - amp[1004] - amp[1003];
  jamp[14] = +std::complex<double> (0, 1) * amp[0] + std::complex<double> (0,
      1) * amp[2] - amp[3] + std::complex<double> (0, 1) * amp[4] + amp[5] +
      std::complex<double> (0, 1) * amp[6] - std::complex<double> (0, 1) *
      amp[14] - std::complex<double> (0, 1) * amp[13] - std::complex<double>
      (0, 1) * amp[17] - std::complex<double> (0, 1) * amp[16] +
      std::complex<double> (0, 1) * amp[18] + std::complex<double> (0, 1) *
      amp[20] - amp[21] + std::complex<double> (0, 1) * amp[22] + amp[23] +
      std::complex<double> (0, 1) * amp[24] - std::complex<double> (0, 1) *
      amp[32] - std::complex<double> (0, 1) * amp[31] - std::complex<double>
      (0, 1) * amp[35] - std::complex<double> (0, 1) * amp[34] +
      std::complex<double> (0, 1) * amp[52] + amp[53] + amp[54] +
      std::complex<double> (0, 1) * amp[58] + amp[59] + amp[60] +
      std::complex<double> (0, 1) * amp[64] + std::complex<double> (0, 1) *
      amp[65] + std::complex<double> (0, 1) * amp[66] + std::complex<double>
      (0, 1) * amp[67] - amp[96] - amp[97] + std::complex<double> (0, 1) *
      amp[98] - amp[102] - amp[103] + std::complex<double> (0, 1) * amp[104] +
      std::complex<double> (0, 1) * amp[216] + std::complex<double> (0, 1) *
      amp[217] + amp[220] + amp[221] + amp[222] + std::complex<double> (0, 1) *
      amp[224] + std::complex<double> (0, 1) * amp[225] + amp[228] + amp[229] +
      amp[230] + std::complex<double> (0, 1) * amp[264] + amp[265] +
      std::complex<double> (0, 1) * amp[266] - std::complex<double> (0, 1) *
      amp[274] - std::complex<double> (0, 1) * amp[278] + std::complex<double>
      (0, 1) * amp[276] - std::complex<double> (0, 1) * amp[281] +
      std::complex<double> (0, 1) * amp[279] + std::complex<double> (0, 1) *
      amp[282] + amp[283] + std::complex<double> (0, 1) * amp[284] -
      std::complex<double> (0, 1) * amp[292] - std::complex<double> (0, 1) *
      amp[296] + std::complex<double> (0, 1) * amp[294] - std::complex<double>
      (0, 1) * amp[299] + std::complex<double> (0, 1) * amp[297] + amp[300] +
      amp[301] - std::complex<double> (0, 1) * amp[303] - std::complex<double>
      (0, 1) * amp[304] - std::complex<double> (0, 1) * amp[305] + amp[306] +
      amp[307] - std::complex<double> (0, 1) * amp[309] - std::complex<double>
      (0, 1) * amp[310] - std::complex<double> (0, 1) * amp[311] +
      std::complex<double> (0, 1) * amp[540] + std::complex<double> (0, 1) *
      amp[542] - amp[543] + std::complex<double> (0, 1) * amp[548] -
      std::complex<double> (0, 1) * amp[554] - std::complex<double> (0, 1) *
      amp[553] - std::complex<double> (0, 1) * amp[557] - std::complex<double>
      (0, 1) * amp[556] + std::complex<double> (0, 1) * amp[558] +
      std::complex<double> (0, 1) * amp[560] - amp[561] + std::complex<double>
      (0, 1) * amp[566] - std::complex<double> (0, 1) * amp[572] -
      std::complex<double> (0, 1) * amp[571] - std::complex<double> (0, 1) *
      amp[575] - std::complex<double> (0, 1) * amp[574] - amp[588] - amp[589] +
      std::complex<double> (0, 1) * amp[590] + std::complex<double> (0, 1) *
      amp[591] + std::complex<double> (0, 1) * amp[593] - amp[594] - amp[595] +
      std::complex<double> (0, 1) * amp[596] + std::complex<double> (0, 1) *
      amp[597] + std::complex<double> (0, 1) * amp[599] + amp[794] + amp[800] +
      amp[806] + amp[807] + amp[808] + amp[812] + amp[813] + amp[814] +
      amp[852] + amp[858] - std::complex<double> (0, 1) * amp[866] +
      std::complex<double> (0, 1) * amp[864] - std::complex<double> (0, 1) *
      amp[869] + std::complex<double> (0, 1) * amp[867] + amp[872] - amp[870] -
      std::complex<double> (0, 1) * amp[878] + std::complex<double> (0, 1) *
      amp[876] - std::complex<double> (0, 1) * amp[881] + std::complex<double>
      (0, 1) * amp[879] + amp[884] - amp[882] - amp[900] + amp[902] - amp[903]
      + amp[905] - amp[906] + amp[908] - amp[909] + amp[911] -
      std::complex<double> (0, 1) * amp[914] - std::complex<double> (0, 1) *
      amp[913] - std::complex<double> (0, 1) * amp[917] - std::complex<double>
      (0, 1) * amp[916] - amp[923] - amp[922] - std::complex<double> (0, 1) *
      amp[926] - std::complex<double> (0, 1) * amp[925] - std::complex<double>
      (0, 1) * amp[929] - std::complex<double> (0, 1) * amp[928] - amp[935] -
      amp[934] - amp[937] - amp[938] - amp[941] - amp[940] - amp[943] -
      amp[944] - amp[947] - amp[946];
  jamp[15] = -std::complex<double> (0, 1) * amp[216] - std::complex<double> (0,
      1) * amp[217] - amp[220] - amp[221] - amp[222] - std::complex<double> (0,
      1) * amp[224] - std::complex<double> (0, 1) * amp[225] - amp[228] -
      amp[229] - amp[230] - std::complex<double> (0, 1) * amp[264] - amp[265] -
      std::complex<double> (0, 1) * amp[266] + std::complex<double> (0, 1) *
      amp[274] + std::complex<double> (0, 1) * amp[278] - std::complex<double>
      (0, 1) * amp[276] + std::complex<double> (0, 1) * amp[281] -
      std::complex<double> (0, 1) * amp[279] - std::complex<double> (0, 1) *
      amp[282] - amp[283] - std::complex<double> (0, 1) * amp[284] +
      std::complex<double> (0, 1) * amp[292] + std::complex<double> (0, 1) *
      amp[296] - std::complex<double> (0, 1) * amp[294] + std::complex<double>
      (0, 1) * amp[299] - std::complex<double> (0, 1) * amp[297] - amp[300] -
      amp[301] + std::complex<double> (0, 1) * amp[303] + std::complex<double>
      (0, 1) * amp[304] + std::complex<double> (0, 1) * amp[305] - amp[306] -
      amp[307] + std::complex<double> (0, 1) * amp[309] + std::complex<double>
      (0, 1) * amp[310] + std::complex<double> (0, 1) * amp[311] + amp[433] +
      std::complex<double> (0, 1) * amp[434] + std::complex<double> (0, 1) *
      amp[435] + std::complex<double> (0, 1) * amp[436] + amp[437] + amp[439] +
      std::complex<double> (0, 1) * amp[440] + std::complex<double> (0, 1) *
      amp[441] + std::complex<double> (0, 1) * amp[442] + amp[443] -
      std::complex<double> (0, 1) * amp[448] + amp[449] + amp[453] -
      std::complex<double> (0, 1) * amp[454] + amp[455] + amp[459] -
      std::complex<double> (0, 1) * amp[460] - std::complex<double> (0, 1) *
      amp[461] - std::complex<double> (0, 1) * amp[462] - std::complex<double>
      (0, 1) * amp[463] + amp[484] + amp[490] + amp[506] + amp[508] + amp[509]
      + amp[512] + amp[514] + amp[515] - amp[529] - amp[530] - amp[533] -
      amp[532] - amp[535] - amp[536] - amp[539] - amp[538] +
      std::complex<double> (0, 1) * amp[544] + std::complex<double> (0, 1) *
      amp[546] + amp[547] + std::complex<double> (0, 1) * amp[549] +
      std::complex<double> (0, 1) * amp[552] + std::complex<double> (0, 1) *
      amp[553] + std::complex<double> (0, 1) * amp[555] + std::complex<double>
      (0, 1) * amp[556] + std::complex<double> (0, 1) * amp[562] +
      std::complex<double> (0, 1) * amp[564] + amp[565] + std::complex<double>
      (0, 1) * amp[567] + std::complex<double> (0, 1) * amp[570] +
      std::complex<double> (0, 1) * amp[571] + std::complex<double> (0, 1) *
      amp[573] + std::complex<double> (0, 1) * amp[574] + std::complex<double>
      (0, 1) * amp[600] + amp[601] + std::complex<double> (0, 1) * amp[602] +
      std::complex<double> (0, 1) * amp[604] + std::complex<double> (0, 1) *
      amp[606] + amp[607] - std::complex<double> (0, 1) * amp[614] +
      std::complex<double> (0, 1) * amp[612] - std::complex<double> (0, 1) *
      amp[617] + std::complex<double> (0, 1) * amp[615] + std::complex<double>
      (0, 1) * amp[618] + amp[619] + std::complex<double> (0, 1) * amp[620] +
      std::complex<double> (0, 1) * amp[622] + std::complex<double> (0, 1) *
      amp[624] + amp[625] - std::complex<double> (0, 1) * amp[632] +
      std::complex<double> (0, 1) * amp[630] - std::complex<double> (0, 1) *
      amp[635] + std::complex<double> (0, 1) * amp[633] + amp[636] + amp[637] -
      std::complex<double> (0, 1) * amp[640] + amp[642] + amp[643] -
      std::complex<double> (0, 1) * amp[646] + amp[795] + amp[801] +
      std::complex<double> (0, 1) * amp[912] + std::complex<double> (0, 1) *
      amp[913] + std::complex<double> (0, 1) * amp[915] + std::complex<double>
      (0, 1) * amp[916] + amp[921] + amp[922] + std::complex<double> (0, 1) *
      amp[924] + std::complex<double> (0, 1) * amp[925] + std::complex<double>
      (0, 1) * amp[927] + std::complex<double> (0, 1) * amp[928] + amp[933] +
      amp[934] + amp[936] + amp[937] + amp[939] + amp[940] + amp[942] +
      amp[943] + amp[945] + amp[946] - std::complex<double> (0, 1) * amp[1010]
      - std::complex<double> (0, 1) * amp[1009] - std::complex<double> (0, 1) *
      amp[1013] - std::complex<double> (0, 1) * amp[1012] - amp[1016] -
      amp[1015] - std::complex<double> (0, 1) * amp[1022] -
      std::complex<double> (0, 1) * amp[1021] - std::complex<double> (0, 1) *
      amp[1025] - std::complex<double> (0, 1) * amp[1024] - amp[1028] -
      amp[1027];
  jamp[16] = +std::complex<double> (0, 1) * amp[4] + amp[5] +
      std::complex<double> (0, 1) * amp[6] - std::complex<double> (0, 1) *
      amp[11] - std::complex<double> (0, 1) * amp[12] - std::complex<double>
      (0, 1) * amp[13] - std::complex<double> (0, 1) * amp[15] -
      std::complex<double> (0, 1) * amp[16] + std::complex<double> (0, 1) *
      amp[22] + amp[23] + std::complex<double> (0, 1) * amp[24] -
      std::complex<double> (0, 1) * amp[29] - std::complex<double> (0, 1) *
      amp[30] - std::complex<double> (0, 1) * amp[31] - std::complex<double>
      (0, 1) * amp[33] - std::complex<double> (0, 1) * amp[34] -
      std::complex<double> (0, 1) * amp[36] - std::complex<double> (0, 1) *
      amp[37] - amp[40] - amp[41] - amp[43] - std::complex<double> (0, 1) *
      amp[44] - std::complex<double> (0, 1) * amp[45] - amp[48] - amp[49] -
      amp[51] + amp[53] + amp[54] - std::complex<double> (0, 1) * amp[55] -
      std::complex<double> (0, 1) * amp[56] - std::complex<double> (0, 1) *
      amp[57] + amp[59] + amp[60] - std::complex<double> (0, 1) * amp[61] -
      std::complex<double> (0, 1) * amp[62] - std::complex<double> (0, 1) *
      amp[63] + std::complex<double> (0, 1) * amp[248] - amp[249] - amp[250] +
      std::complex<double> (0, 1) * amp[254] - amp[255] - amp[256] +
      std::complex<double> (0, 1) * amp[260] + std::complex<double> (0, 1) *
      amp[261] + std::complex<double> (0, 1) * amp[262] + std::complex<double>
      (0, 1) * amp[263] + std::complex<double> (0, 1) * amp[264] + amp[265] +
      std::complex<double> (0, 1) * amp[266] + std::complex<double> (0, 1) *
      amp[268] + std::complex<double> (0, 1) * amp[270] - amp[271] -
      std::complex<double> (0, 1) * amp[278] - std::complex<double> (0, 1) *
      amp[277] - std::complex<double> (0, 1) * amp[281] - std::complex<double>
      (0, 1) * amp[280] + std::complex<double> (0, 1) * amp[282] + amp[283] +
      std::complex<double> (0, 1) * amp[284] + std::complex<double> (0, 1) *
      amp[286] + std::complex<double> (0, 1) * amp[288] - amp[289] -
      std::complex<double> (0, 1) * amp[296] - std::complex<double> (0, 1) *
      amp[295] - std::complex<double> (0, 1) * amp[299] - std::complex<double>
      (0, 1) * amp[298] + amp[300] + amp[301] + std::complex<double> (0, 1) *
      amp[302] + amp[306] + amp[307] + std::complex<double> (0, 1) * amp[308] +
      amp[709] + amp[715] - std::complex<double> (0, 1) * amp[720] -
      std::complex<double> (0, 1) * amp[721] + amp[722] + amp[723] -
      std::complex<double> (0, 1) * amp[724] - std::complex<double> (0, 1) *
      amp[726] - std::complex<double> (0, 1) * amp[727] + amp[728] + amp[729] -
      std::complex<double> (0, 1) * amp[730] + amp[732] + amp[733] + amp[736] +
      amp[738] + amp[739] + amp[742] - std::complex<double> (0, 1) * amp[756] -
      std::complex<double> (0, 1) * amp[758] + amp[759] - std::complex<double>
      (0, 1) * amp[767] + std::complex<double> (0, 1) * amp[770] +
      std::complex<double> (0, 1) * amp[769] + std::complex<double> (0, 1) *
      amp[773] + std::complex<double> (0, 1) * amp[772] - std::complex<double>
      (0, 1) * amp[774] - std::complex<double> (0, 1) * amp[776] + amp[777] -
      std::complex<double> (0, 1) * amp[785] + std::complex<double> (0, 1) *
      amp[788] + std::complex<double> (0, 1) * amp[787] + std::complex<double>
      (0, 1) * amp[791] + std::complex<double> (0, 1) * amp[790] + amp[797] +
      amp[803] - std::complex<double> (0, 1) * amp[914] - std::complex<double>
      (0, 1) * amp[913] - std::complex<double> (0, 1) * amp[917] -
      std::complex<double> (0, 1) * amp[916] - amp[923] - amp[922] -
      std::complex<double> (0, 1) * amp[926] - std::complex<double> (0, 1) *
      amp[925] - std::complex<double> (0, 1) * amp[929] - std::complex<double>
      (0, 1) * amp[928] - amp[935] - amp[934] - amp[937] - amp[938] - amp[941]
      - amp[940] - amp[943] - amp[944] - amp[947] - amp[946] + amp[972] +
      amp[973] + amp[975] + amp[976] + amp[978] + amp[979] + amp[981] +
      amp[982] - std::complex<double> (0, 1) * amp[984] - std::complex<double>
      (0, 1) * amp[985] - std::complex<double> (0, 1) * amp[987] -
      std::complex<double> (0, 1) * amp[988] + amp[990] + amp[991] -
      std::complex<double> (0, 1) * amp[996] - std::complex<double> (0, 1) *
      amp[997] - std::complex<double> (0, 1) * amp[999] - std::complex<double>
      (0, 1) * amp[1000] + amp[1002] + amp[1003];
  jamp[17] = -std::complex<double> (0, 1) * amp[4] - amp[5] -
      std::complex<double> (0, 1) * amp[6] + std::complex<double> (0, 1) *
      amp[11] + std::complex<double> (0, 1) * amp[12] + std::complex<double>
      (0, 1) * amp[13] + std::complex<double> (0, 1) * amp[15] +
      std::complex<double> (0, 1) * amp[16] - std::complex<double> (0, 1) *
      amp[22] - amp[23] - std::complex<double> (0, 1) * amp[24] +
      std::complex<double> (0, 1) * amp[29] + std::complex<double> (0, 1) *
      amp[30] + std::complex<double> (0, 1) * amp[31] + std::complex<double>
      (0, 1) * amp[33] + std::complex<double> (0, 1) * amp[34] +
      std::complex<double> (0, 1) * amp[36] + std::complex<double> (0, 1) *
      amp[37] + amp[40] + amp[41] + amp[43] + std::complex<double> (0, 1) *
      amp[44] + std::complex<double> (0, 1) * amp[45] + amp[48] + amp[49] +
      amp[51] - amp[53] - amp[54] + std::complex<double> (0, 1) * amp[55] +
      std::complex<double> (0, 1) * amp[56] + std::complex<double> (0, 1) *
      amp[57] - amp[59] - amp[60] + std::complex<double> (0, 1) * amp[61] +
      std::complex<double> (0, 1) * amp[62] + std::complex<double> (0, 1) *
      amp[63] + std::complex<double> (0, 1) * amp[448] - amp[449] - amp[453] +
      std::complex<double> (0, 1) * amp[454] - amp[455] - amp[459] +
      std::complex<double> (0, 1) * amp[460] + std::complex<double> (0, 1) *
      amp[461] + std::complex<double> (0, 1) * amp[462] + std::complex<double>
      (0, 1) * amp[463] + amp[481] + amp[487] - std::complex<double> (0, 1) *
      amp[492] - std::complex<double> (0, 1) * amp[493] + amp[494] -
      std::complex<double> (0, 1) * amp[495] + amp[497] - std::complex<double>
      (0, 1) * amp[498] - std::complex<double> (0, 1) * amp[499] + amp[500] -
      std::complex<double> (0, 1) * amp[501] + amp[503] + amp[504] + amp[505] +
      amp[507] + amp[510] + amp[511] + amp[513] + amp[528] + amp[529] +
      amp[531] + amp[532] + amp[534] + amp[535] + amp[537] + amp[538] -
      std::complex<double> (0, 1) * amp[600] - amp[601] - std::complex<double>
      (0, 1) * amp[602] - std::complex<double> (0, 1) * amp[604] -
      std::complex<double> (0, 1) * amp[606] - amp[607] + std::complex<double>
      (0, 1) * amp[614] - std::complex<double> (0, 1) * amp[612] +
      std::complex<double> (0, 1) * amp[617] - std::complex<double> (0, 1) *
      amp[615] - std::complex<double> (0, 1) * amp[618] - amp[619] -
      std::complex<double> (0, 1) * amp[620] - std::complex<double> (0, 1) *
      amp[622] - std::complex<double> (0, 1) * amp[624] - amp[625] +
      std::complex<double> (0, 1) * amp[632] - std::complex<double> (0, 1) *
      amp[630] + std::complex<double> (0, 1) * amp[635] - std::complex<double>
      (0, 1) * amp[633] - amp[636] - amp[637] + std::complex<double> (0, 1) *
      amp[640] - amp[642] - amp[643] + std::complex<double> (0, 1) * amp[646] +
      std::complex<double> (0, 1) * amp[760] + std::complex<double> (0, 1) *
      amp[762] + amp[763] - std::complex<double> (0, 1) * amp[766] -
      std::complex<double> (0, 1) * amp[770] + std::complex<double> (0, 1) *
      amp[768] - std::complex<double> (0, 1) * amp[773] + std::complex<double>
      (0, 1) * amp[771] + std::complex<double> (0, 1) * amp[778] +
      std::complex<double> (0, 1) * amp[780] + amp[781] - std::complex<double>
      (0, 1) * amp[784] - std::complex<double> (0, 1) * amp[788] +
      std::complex<double> (0, 1) * amp[786] - std::complex<double> (0, 1) *
      amp[791] + std::complex<double> (0, 1) * amp[789] + amp[796] + amp[802] +
      std::complex<double> (0, 1) * amp[914] - std::complex<double> (0, 1) *
      amp[912] + std::complex<double> (0, 1) * amp[917] - std::complex<double>
      (0, 1) * amp[915] + amp[923] - amp[921] + std::complex<double> (0, 1) *
      amp[926] - std::complex<double> (0, 1) * amp[924] + std::complex<double>
      (0, 1) * amp[929] - std::complex<double> (0, 1) * amp[927] + amp[935] -
      amp[933] - amp[936] + amp[938] + amp[941] - amp[939] - amp[942] +
      amp[944] + amp[947] - amp[945] + std::complex<double> (0, 1) * amp[1008]
      + std::complex<double> (0, 1) * amp[1009] + std::complex<double> (0, 1) *
      amp[1011] + std::complex<double> (0, 1) * amp[1012] + amp[1014] +
      amp[1015] + std::complex<double> (0, 1) * amp[1020] +
      std::complex<double> (0, 1) * amp[1021] + std::complex<double> (0, 1) *
      amp[1023] + std::complex<double> (0, 1) * amp[1024] + amp[1026] +
      amp[1027];
  jamp[18] = -std::complex<double> (0, 1) * amp[0] - amp[1] -
      std::complex<double> (0, 1) * amp[2] - std::complex<double> (0, 1) *
      amp[4] - std::complex<double> (0, 1) * amp[6] + amp[7] +
      std::complex<double> (0, 1) * amp[14] + std::complex<double> (0, 1) *
      amp[13] + std::complex<double> (0, 1) * amp[17] + std::complex<double>
      (0, 1) * amp[16] - std::complex<double> (0, 1) * amp[18] - amp[19] -
      std::complex<double> (0, 1) * amp[20] - std::complex<double> (0, 1) *
      amp[22] - std::complex<double> (0, 1) * amp[24] + amp[25] +
      std::complex<double> (0, 1) * amp[32] + std::complex<double> (0, 1) *
      amp[31] + std::complex<double> (0, 1) * amp[35] + std::complex<double>
      (0, 1) * amp[34] - std::complex<double> (0, 1) * amp[68] + amp[69] +
      amp[70] - std::complex<double> (0, 1) * amp[74] + amp[75] + amp[76] -
      std::complex<double> (0, 1) * amp[80] - std::complex<double> (0, 1) *
      amp[81] - std::complex<double> (0, 1) * amp[82] - std::complex<double>
      (0, 1) * amp[83] - amp[84] - amp[85] - std::complex<double> (0, 1) *
      amp[86] - amp[90] - amp[91] - std::complex<double> (0, 1) * amp[92] -
      std::complex<double> (0, 1) * amp[216] - std::complex<double> (0, 1) *
      amp[217] + amp[218] + amp[219] + amp[223] - std::complex<double> (0, 1) *
      amp[224] - std::complex<double> (0, 1) * amp[225] + amp[226] + amp[227] +
      amp[231] - std::complex<double> (0, 1) * amp[264] - std::complex<double>
      (0, 1) * amp[266] + amp[267] + std::complex<double> (0, 1) * amp[272] +
      std::complex<double> (0, 1) * amp[278] - std::complex<double> (0, 1) *
      amp[276] + std::complex<double> (0, 1) * amp[281] - std::complex<double>
      (0, 1) * amp[279] - std::complex<double> (0, 1) * amp[282] -
      std::complex<double> (0, 1) * amp[284] + amp[285] + std::complex<double>
      (0, 1) * amp[290] + std::complex<double> (0, 1) * amp[296] -
      std::complex<double> (0, 1) * amp[294] + std::complex<double> (0, 1) *
      amp[299] - std::complex<double> (0, 1) * amp[297] + amp[312] + amp[313] +
      std::complex<double> (0, 1) * amp[315] + std::complex<double> (0, 1) *
      amp[316] + std::complex<double> (0, 1) * amp[317] + amp[318] + amp[319] +
      std::complex<double> (0, 1) * amp[321] + std::complex<double> (0, 1) *
      amp[322] + std::complex<double> (0, 1) * amp[323] - std::complex<double>
      (0, 1) * amp[540] - amp[541] - std::complex<double> (0, 1) * amp[542] -
      std::complex<double> (0, 1) * amp[550] + std::complex<double> (0, 1) *
      amp[554] + std::complex<double> (0, 1) * amp[553] + std::complex<double>
      (0, 1) * amp[557] + std::complex<double> (0, 1) * amp[556] -
      std::complex<double> (0, 1) * amp[558] - amp[559] - std::complex<double>
      (0, 1) * amp[560] - std::complex<double> (0, 1) * amp[568] +
      std::complex<double> (0, 1) * amp[572] + std::complex<double> (0, 1) *
      amp[571] + std::complex<double> (0, 1) * amp[575] + std::complex<double>
      (0, 1) * amp[574] - amp[576] - amp[577] - std::complex<double> (0, 1) *
      amp[578] - std::complex<double> (0, 1) * amp[579] - std::complex<double>
      (0, 1) * amp[581] - amp[582] - amp[583] - std::complex<double> (0, 1) *
      amp[584] - std::complex<double> (0, 1) * amp[585] - std::complex<double>
      (0, 1) * amp[587] + amp[818] + amp[824] + amp[830] + amp[831] + amp[832]
      + amp[836] + amp[837] + amp[838] + amp[840] + amp[846] +
      std::complex<double> (0, 1) * amp[866] - std::complex<double> (0, 1) *
      amp[864] + std::complex<double> (0, 1) * amp[869] - std::complex<double>
      (0, 1) * amp[867] + amp[875] - amp[873] + std::complex<double> (0, 1) *
      amp[878] - std::complex<double> (0, 1) * amp[876] + std::complex<double>
      (0, 1) * amp[881] - std::complex<double> (0, 1) * amp[879] + amp[887] -
      amp[885] - amp[888] + amp[890] + amp[893] - amp[891] - amp[894] +
      amp[896] + amp[899] - amp[897] + std::complex<double> (0, 1) * amp[914] +
      std::complex<double> (0, 1) * amp[913] + std::complex<double> (0, 1) *
      amp[917] + std::complex<double> (0, 1) * amp[916] - amp[920] - amp[919] +
      std::complex<double> (0, 1) * amp[926] + std::complex<double> (0, 1) *
      amp[925] + std::complex<double> (0, 1) * amp[929] + std::complex<double>
      (0, 1) * amp[928] - amp[932] - amp[931] - amp[949] - amp[950] - amp[952]
      - amp[953] - amp[955] - amp[956] - amp[958] - amp[959];
  jamp[19] = -std::complex<double> (0, 1) * amp[140] + amp[141] + amp[142] -
      std::complex<double> (0, 1) * amp[146] + amp[147] + amp[148] -
      std::complex<double> (0, 1) * amp[152] - std::complex<double> (0, 1) *
      amp[153] - std::complex<double> (0, 1) * amp[154] - std::complex<double>
      (0, 1) * amp[155] - std::complex<double> (0, 1) * amp[156] - amp[157] -
      std::complex<double> (0, 1) * amp[158] - std::complex<double> (0, 1) *
      amp[160] - std::complex<double> (0, 1) * amp[162] + amp[163] +
      std::complex<double> (0, 1) * amp[170] + std::complex<double> (0, 1) *
      amp[169] + std::complex<double> (0, 1) * amp[173] + std::complex<double>
      (0, 1) * amp[172] - std::complex<double> (0, 1) * amp[174] - amp[175] -
      std::complex<double> (0, 1) * amp[176] - std::complex<double> (0, 1) *
      amp[178] - std::complex<double> (0, 1) * amp[180] + amp[181] +
      std::complex<double> (0, 1) * amp[188] + std::complex<double> (0, 1) *
      amp[187] + std::complex<double> (0, 1) * amp[191] + std::complex<double>
      (0, 1) * amp[190] - amp[192] - amp[193] - std::complex<double> (0, 1) *
      amp[194] - amp[198] - amp[199] - std::complex<double> (0, 1) * amp[200] +
      std::complex<double> (0, 1) * amp[216] + std::complex<double> (0, 1) *
      amp[217] - amp[218] - amp[219] - amp[223] + std::complex<double> (0, 1) *
      amp[224] + std::complex<double> (0, 1) * amp[225] - amp[226] - amp[227] -
      amp[231] + amp[249] + amp[250] + std::complex<double> (0, 1) * amp[251] +
      std::complex<double> (0, 1) * amp[252] + std::complex<double> (0, 1) *
      amp[253] + amp[255] + amp[256] + std::complex<double> (0, 1) * amp[257] +
      std::complex<double> (0, 1) * amp[258] + std::complex<double> (0, 1) *
      amp[259] - std::complex<double> (0, 1) * amp[268] - std::complex<double>
      (0, 1) * amp[270] + amp[271] + std::complex<double> (0, 1) * amp[273] +
      std::complex<double> (0, 1) * amp[276] + std::complex<double> (0, 1) *
      amp[277] + std::complex<double> (0, 1) * amp[279] + std::complex<double>
      (0, 1) * amp[280] - std::complex<double> (0, 1) * amp[286] -
      std::complex<double> (0, 1) * amp[288] + amp[289] + std::complex<double>
      (0, 1) * amp[291] + std::complex<double> (0, 1) * amp[294] +
      std::complex<double> (0, 1) * amp[295] + std::complex<double> (0, 1) *
      amp[297] + std::complex<double> (0, 1) * amp[298] + std::complex<double>
      (0, 1) * amp[540] + amp[541] + std::complex<double> (0, 1) * amp[542] +
      std::complex<double> (0, 1) * amp[550] - std::complex<double> (0, 1) *
      amp[554] - std::complex<double> (0, 1) * amp[553] - std::complex<double>
      (0, 1) * amp[557] - std::complex<double> (0, 1) * amp[556] +
      std::complex<double> (0, 1) * amp[558] + amp[559] + std::complex<double>
      (0, 1) * amp[560] + std::complex<double> (0, 1) * amp[568] -
      std::complex<double> (0, 1) * amp[572] - std::complex<double> (0, 1) *
      amp[571] - std::complex<double> (0, 1) * amp[575] - std::complex<double>
      (0, 1) * amp[574] + amp[576] + amp[577] + std::complex<double> (0, 1) *
      amp[578] + std::complex<double> (0, 1) * amp[579] + std::complex<double>
      (0, 1) * amp[581] + amp[582] + amp[583] + std::complex<double> (0, 1) *
      amp[584] + std::complex<double> (0, 1) * amp[585] + std::complex<double>
      (0, 1) * amp[587] + amp[710] + amp[716] + amp[746] + amp[747] + amp[749]
      + amp[752] + amp[753] + amp[755] + amp[841] + amp[847] +
      std::complex<double> (0, 1) * amp[864] + std::complex<double> (0, 1) *
      amp[865] + std::complex<double> (0, 1) * amp[867] + std::complex<double>
      (0, 1) * amp[868] + amp[873] + amp[874] + std::complex<double> (0, 1) *
      amp[876] + std::complex<double> (0, 1) * amp[877] + std::complex<double>
      (0, 1) * amp[879] + std::complex<double> (0, 1) * amp[880] + amp[885] +
      amp[886] + amp[888] + amp[889] + amp[891] + amp[892] + amp[894] +
      amp[895] + amp[897] + amp[898] - amp[973] - amp[974] - amp[976] -
      amp[977] - amp[979] - amp[980] - amp[982] - amp[983] +
      std::complex<double> (0, 1) * amp[986] + std::complex<double> (0, 1) *
      amp[985] + std::complex<double> (0, 1) * amp[989] + std::complex<double>
      (0, 1) * amp[988] - amp[992] - amp[991] + std::complex<double> (0, 1) *
      amp[998] + std::complex<double> (0, 1) * amp[997] + std::complex<double>
      (0, 1) * amp[1001] + std::complex<double> (0, 1) * amp[1000] - amp[1004]
      - amp[1003];
  jamp[20] = +std::complex<double> (0, 1) * amp[0] + amp[1] +
      std::complex<double> (0, 1) * amp[2] + std::complex<double> (0, 1) *
      amp[4] + std::complex<double> (0, 1) * amp[6] - amp[7] -
      std::complex<double> (0, 1) * amp[14] - std::complex<double> (0, 1) *
      amp[13] - std::complex<double> (0, 1) * amp[17] - std::complex<double>
      (0, 1) * amp[16] + std::complex<double> (0, 1) * amp[18] + amp[19] +
      std::complex<double> (0, 1) * amp[20] + std::complex<double> (0, 1) *
      amp[22] + std::complex<double> (0, 1) * amp[24] - amp[25] -
      std::complex<double> (0, 1) * amp[32] - std::complex<double> (0, 1) *
      amp[31] - std::complex<double> (0, 1) * amp[35] - std::complex<double>
      (0, 1) * amp[34] + std::complex<double> (0, 1) * amp[68] - amp[69] -
      amp[70] + std::complex<double> (0, 1) * amp[74] - amp[75] - amp[76] +
      std::complex<double> (0, 1) * amp[80] + std::complex<double> (0, 1) *
      amp[81] + std::complex<double> (0, 1) * amp[82] + std::complex<double>
      (0, 1) * amp[83] + amp[84] + amp[85] + std::complex<double> (0, 1) *
      amp[86] + amp[90] + amp[91] + std::complex<double> (0, 1) * amp[92] +
      std::complex<double> (0, 1) * amp[108] + std::complex<double> (0, 1) *
      amp[109] + amp[112] + amp[113] + amp[114] + std::complex<double> (0, 1) *
      amp[116] + std::complex<double> (0, 1) * amp[117] + amp[120] + amp[121] +
      amp[122] + std::complex<double> (0, 1) * amp[156] + amp[157] +
      std::complex<double> (0, 1) * amp[158] - std::complex<double> (0, 1) *
      amp[166] - std::complex<double> (0, 1) * amp[170] + std::complex<double>
      (0, 1) * amp[168] - std::complex<double> (0, 1) * amp[173] +
      std::complex<double> (0, 1) * amp[171] + std::complex<double> (0, 1) *
      amp[174] + amp[175] + std::complex<double> (0, 1) * amp[176] -
      std::complex<double> (0, 1) * amp[184] - std::complex<double> (0, 1) *
      amp[188] + std::complex<double> (0, 1) * amp[186] - std::complex<double>
      (0, 1) * amp[191] + std::complex<double> (0, 1) * amp[189] + amp[192] +
      amp[193] - std::complex<double> (0, 1) * amp[195] - std::complex<double>
      (0, 1) * amp[196] - std::complex<double> (0, 1) * amp[197] + amp[198] +
      amp[199] - std::complex<double> (0, 1) * amp[201] - std::complex<double>
      (0, 1) * amp[202] - std::complex<double> (0, 1) * amp[203] +
      std::complex<double> (0, 1) * amp[600] + std::complex<double> (0, 1) *
      amp[602] - amp[603] + std::complex<double> (0, 1) * amp[608] -
      std::complex<double> (0, 1) * amp[614] - std::complex<double> (0, 1) *
      amp[613] - std::complex<double> (0, 1) * amp[617] - std::complex<double>
      (0, 1) * amp[616] + std::complex<double> (0, 1) * amp[618] +
      std::complex<double> (0, 1) * amp[620] - amp[621] + std::complex<double>
      (0, 1) * amp[626] - std::complex<double> (0, 1) * amp[632] -
      std::complex<double> (0, 1) * amp[631] - std::complex<double> (0, 1) *
      amp[635] - std::complex<double> (0, 1) * amp[634] - amp[648] - amp[649] +
      std::complex<double> (0, 1) * amp[650] + std::complex<double> (0, 1) *
      amp[651] + std::complex<double> (0, 1) * amp[653] - amp[654] - amp[655] +
      std::complex<double> (0, 1) * amp[656] + std::complex<double> (0, 1) *
      amp[657] + std::complex<double> (0, 1) * amp[659] + amp[816] + amp[822] +
      amp[828] + amp[829] + amp[833] + amp[834] + amp[835] + amp[839] +
      amp[842] + amp[848] - std::complex<double> (0, 1) * amp[866] -
      std::complex<double> (0, 1) * amp[865] - std::complex<double> (0, 1) *
      amp[869] - std::complex<double> (0, 1) * amp[868] - amp[875] - amp[874] -
      std::complex<double> (0, 1) * amp[878] - std::complex<double> (0, 1) *
      amp[877] - std::complex<double> (0, 1) * amp[881] - std::complex<double>
      (0, 1) * amp[880] - amp[887] - amp[886] - amp[889] - amp[890] - amp[893]
      - amp[892] - amp[895] - amp[896] - amp[899] - amp[898] -
      std::complex<double> (0, 1) * amp[914] + std::complex<double> (0, 1) *
      amp[912] - std::complex<double> (0, 1) * amp[917] + std::complex<double>
      (0, 1) * amp[915] + amp[920] - amp[918] - std::complex<double> (0, 1) *
      amp[926] + std::complex<double> (0, 1) * amp[924] - std::complex<double>
      (0, 1) * amp[929] + std::complex<double> (0, 1) * amp[927] + amp[932] -
      amp[930] - amp[948] + amp[950] - amp[951] + amp[953] - amp[954] +
      amp[956] - amp[957] + amp[959];
  jamp[21] = -std::complex<double> (0, 1) * amp[108] - std::complex<double> (0,
      1) * amp[109] - amp[112] - amp[113] - amp[114] - std::complex<double> (0,
      1) * amp[116] - std::complex<double> (0, 1) * amp[117] - amp[120] -
      amp[121] - amp[122] - std::complex<double> (0, 1) * amp[156] - amp[157] -
      std::complex<double> (0, 1) * amp[158] + std::complex<double> (0, 1) *
      amp[166] + std::complex<double> (0, 1) * amp[170] - std::complex<double>
      (0, 1) * amp[168] + std::complex<double> (0, 1) * amp[173] -
      std::complex<double> (0, 1) * amp[171] - std::complex<double> (0, 1) *
      amp[174] - amp[175] - std::complex<double> (0, 1) * amp[176] +
      std::complex<double> (0, 1) * amp[184] + std::complex<double> (0, 1) *
      amp[188] - std::complex<double> (0, 1) * amp[186] + std::complex<double>
      (0, 1) * amp[191] - std::complex<double> (0, 1) * amp[189] - amp[192] -
      amp[193] + std::complex<double> (0, 1) * amp[195] + std::complex<double>
      (0, 1) * amp[196] + std::complex<double> (0, 1) * amp[197] - amp[198] -
      amp[199] + std::complex<double> (0, 1) * amp[201] + std::complex<double>
      (0, 1) * amp[202] + std::complex<double> (0, 1) * amp[203] -
      std::complex<double> (0, 1) * amp[432] + amp[433] + amp[437] -
      std::complex<double> (0, 1) * amp[438] + amp[439] + amp[443] -
      std::complex<double> (0, 1) * amp[444] - std::complex<double> (0, 1) *
      amp[445] - std::complex<double> (0, 1) * amp[446] - std::complex<double>
      (0, 1) * amp[447] + amp[449] + std::complex<double> (0, 1) * amp[450] +
      std::complex<double> (0, 1) * amp[451] + std::complex<double> (0, 1) *
      amp[452] + amp[453] + amp[455] + std::complex<double> (0, 1) * amp[456] +
      std::complex<double> (0, 1) * amp[457] + std::complex<double> (0, 1) *
      amp[458] + amp[459] + amp[482] + amp[488] + amp[518] + amp[520] +
      amp[521] + amp[524] + amp[526] + amp[527] - amp[529] - amp[530] -
      amp[533] - amp[532] - amp[535] - amp[536] - amp[539] - amp[538] +
      std::complex<double> (0, 1) * amp[540] + amp[541] + std::complex<double>
      (0, 1) * amp[542] + std::complex<double> (0, 1) * amp[544] +
      std::complex<double> (0, 1) * amp[546] + amp[547] - std::complex<double>
      (0, 1) * amp[554] + std::complex<double> (0, 1) * amp[552] -
      std::complex<double> (0, 1) * amp[557] + std::complex<double> (0, 1) *
      amp[555] + std::complex<double> (0, 1) * amp[558] + amp[559] +
      std::complex<double> (0, 1) * amp[560] + std::complex<double> (0, 1) *
      amp[562] + std::complex<double> (0, 1) * amp[564] + amp[565] -
      std::complex<double> (0, 1) * amp[572] + std::complex<double> (0, 1) *
      amp[570] - std::complex<double> (0, 1) * amp[575] + std::complex<double>
      (0, 1) * amp[573] + amp[576] + amp[577] - std::complex<double> (0, 1) *
      amp[580] + amp[582] + amp[583] - std::complex<double> (0, 1) * amp[586] +
      std::complex<double> (0, 1) * amp[604] + std::complex<double> (0, 1) *
      amp[606] + amp[607] + std::complex<double> (0, 1) * amp[609] +
      std::complex<double> (0, 1) * amp[612] + std::complex<double> (0, 1) *
      amp[613] + std::complex<double> (0, 1) * amp[615] + std::complex<double>
      (0, 1) * amp[616] + std::complex<double> (0, 1) * amp[622] +
      std::complex<double> (0, 1) * amp[624] + amp[625] + std::complex<double>
      (0, 1) * amp[627] + std::complex<double> (0, 1) * amp[630] +
      std::complex<double> (0, 1) * amp[631] + std::complex<double> (0, 1) *
      amp[633] + std::complex<double> (0, 1) * amp[634] + amp[843] + amp[849] +
      std::complex<double> (0, 1) * amp[864] + std::complex<double> (0, 1) *
      amp[865] + std::complex<double> (0, 1) * amp[867] + std::complex<double>
      (0, 1) * amp[868] + amp[873] + amp[874] + std::complex<double> (0, 1) *
      amp[876] + std::complex<double> (0, 1) * amp[877] + std::complex<double>
      (0, 1) * amp[879] + std::complex<double> (0, 1) * amp[880] + amp[885] +
      amp[886] + amp[888] + amp[889] + amp[891] + amp[892] + amp[894] +
      amp[895] + amp[897] + amp[898] - std::complex<double> (0, 1) * amp[1010]
      - std::complex<double> (0, 1) * amp[1009] - std::complex<double> (0, 1) *
      amp[1013] - std::complex<double> (0, 1) * amp[1012] - amp[1016] -
      amp[1015] - std::complex<double> (0, 1) * amp[1022] -
      std::complex<double> (0, 1) * amp[1021] - std::complex<double> (0, 1) *
      amp[1025] - std::complex<double> (0, 1) * amp[1024] - amp[1028] -
      amp[1027];
  jamp[22] = +std::complex<double> (0, 1) * amp[0] + amp[1] +
      std::complex<double> (0, 1) * amp[2] - std::complex<double> (0, 1) *
      amp[10] - std::complex<double> (0, 1) * amp[14] + std::complex<double>
      (0, 1) * amp[12] - std::complex<double> (0, 1) * amp[17] +
      std::complex<double> (0, 1) * amp[15] + std::complex<double> (0, 1) *
      amp[18] + amp[19] + std::complex<double> (0, 1) * amp[20] -
      std::complex<double> (0, 1) * amp[28] - std::complex<double> (0, 1) *
      amp[32] + std::complex<double> (0, 1) * amp[30] - std::complex<double>
      (0, 1) * amp[35] + std::complex<double> (0, 1) * amp[33] +
      std::complex<double> (0, 1) * amp[36] + std::complex<double> (0, 1) *
      amp[37] + amp[40] + amp[41] + amp[43] + std::complex<double> (0, 1) *
      amp[44] + std::complex<double> (0, 1) * amp[45] + amp[48] + amp[49] +
      amp[51] + amp[84] + amp[85] - std::complex<double> (0, 1) * amp[87] -
      std::complex<double> (0, 1) * amp[88] - std::complex<double> (0, 1) *
      amp[89] + amp[90] + amp[91] - std::complex<double> (0, 1) * amp[93] -
      std::complex<double> (0, 1) * amp[94] - std::complex<double> (0, 1) *
      amp[95] + std::complex<double> (0, 1) * amp[140] - amp[141] - amp[142] +
      std::complex<double> (0, 1) * amp[146] - amp[147] - amp[148] +
      std::complex<double> (0, 1) * amp[152] + std::complex<double> (0, 1) *
      amp[153] + std::complex<double> (0, 1) * amp[154] + std::complex<double>
      (0, 1) * amp[155] + std::complex<double> (0, 1) * amp[156] + amp[157] +
      std::complex<double> (0, 1) * amp[158] + std::complex<double> (0, 1) *
      amp[160] + std::complex<double> (0, 1) * amp[162] - amp[163] -
      std::complex<double> (0, 1) * amp[170] - std::complex<double> (0, 1) *
      amp[169] - std::complex<double> (0, 1) * amp[173] - std::complex<double>
      (0, 1) * amp[172] + std::complex<double> (0, 1) * amp[174] + amp[175] +
      std::complex<double> (0, 1) * amp[176] + std::complex<double> (0, 1) *
      amp[178] + std::complex<double> (0, 1) * amp[180] - amp[181] -
      std::complex<double> (0, 1) * amp[188] - std::complex<double> (0, 1) *
      amp[187] - std::complex<double> (0, 1) * amp[191] - std::complex<double>
      (0, 1) * amp[190] + amp[192] + amp[193] + std::complex<double> (0, 1) *
      amp[194] + amp[198] + amp[199] + std::complex<double> (0, 1) * amp[200] +
      amp[708] + amp[714] + std::complex<double> (0, 1) * amp[720] +
      std::complex<double> (0, 1) * amp[721] - amp[722] - amp[723] +
      std::complex<double> (0, 1) * amp[724] + std::complex<double> (0, 1) *
      amp[726] + std::complex<double> (0, 1) * amp[727] - amp[728] - amp[729] +
      std::complex<double> (0, 1) * amp[730] + amp[744] + amp[745] + amp[748] +
      amp[750] + amp[751] + amp[754] + std::complex<double> (0, 1) * amp[756] +
      std::complex<double> (0, 1) * amp[758] - amp[759] + std::complex<double>
      (0, 1) * amp[767] - std::complex<double> (0, 1) * amp[770] -
      std::complex<double> (0, 1) * amp[769] - std::complex<double> (0, 1) *
      amp[773] - std::complex<double> (0, 1) * amp[772] + std::complex<double>
      (0, 1) * amp[774] + std::complex<double> (0, 1) * amp[776] - amp[777] +
      std::complex<double> (0, 1) * amp[785] - std::complex<double> (0, 1) *
      amp[788] - std::complex<double> (0, 1) * amp[787] - std::complex<double>
      (0, 1) * amp[791] - std::complex<double> (0, 1) * amp[790] + amp[845] +
      amp[851] - std::complex<double> (0, 1) * amp[866] - std::complex<double>
      (0, 1) * amp[865] - std::complex<double> (0, 1) * amp[869] -
      std::complex<double> (0, 1) * amp[868] - amp[875] - amp[874] -
      std::complex<double> (0, 1) * amp[878] - std::complex<double> (0, 1) *
      amp[877] - std::complex<double> (0, 1) * amp[881] - std::complex<double>
      (0, 1) * amp[880] - amp[887] - amp[886] - amp[889] - amp[890] - amp[893]
      - amp[892] - amp[895] - amp[896] - amp[899] - amp[898] - amp[972] +
      amp[974] - amp[975] + amp[977] - amp[978] + amp[980] - amp[981] +
      amp[983] - std::complex<double> (0, 1) * amp[986] + std::complex<double>
      (0, 1) * amp[984] - std::complex<double> (0, 1) * amp[989] +
      std::complex<double> (0, 1) * amp[987] + amp[992] - amp[990] -
      std::complex<double> (0, 1) * amp[998] + std::complex<double> (0, 1) *
      amp[996] - std::complex<double> (0, 1) * amp[1001] + std::complex<double>
      (0, 1) * amp[999] + amp[1004] - amp[1002];
  jamp[23] = -std::complex<double> (0, 1) * amp[0] - amp[1] -
      std::complex<double> (0, 1) * amp[2] + std::complex<double> (0, 1) *
      amp[10] + std::complex<double> (0, 1) * amp[14] - std::complex<double>
      (0, 1) * amp[12] + std::complex<double> (0, 1) * amp[17] -
      std::complex<double> (0, 1) * amp[15] - std::complex<double> (0, 1) *
      amp[18] - amp[19] - std::complex<double> (0, 1) * amp[20] +
      std::complex<double> (0, 1) * amp[28] + std::complex<double> (0, 1) *
      amp[32] - std::complex<double> (0, 1) * amp[30] + std::complex<double>
      (0, 1) * amp[35] - std::complex<double> (0, 1) * amp[33] -
      std::complex<double> (0, 1) * amp[36] - std::complex<double> (0, 1) *
      amp[37] - amp[40] - amp[41] - amp[43] - std::complex<double> (0, 1) *
      amp[44] - std::complex<double> (0, 1) * amp[45] - amp[48] - amp[49] -
      amp[51] - amp[84] - amp[85] + std::complex<double> (0, 1) * amp[87] +
      std::complex<double> (0, 1) * amp[88] + std::complex<double> (0, 1) *
      amp[89] - amp[90] - amp[91] + std::complex<double> (0, 1) * amp[93] +
      std::complex<double> (0, 1) * amp[94] + std::complex<double> (0, 1) *
      amp[95] + std::complex<double> (0, 1) * amp[432] - amp[433] - amp[437] +
      std::complex<double> (0, 1) * amp[438] - amp[439] - amp[443] +
      std::complex<double> (0, 1) * amp[444] + std::complex<double> (0, 1) *
      amp[445] + std::complex<double> (0, 1) * amp[446] + std::complex<double>
      (0, 1) * amp[447] + amp[480] + amp[486] + std::complex<double> (0, 1) *
      amp[492] + std::complex<double> (0, 1) * amp[493] - amp[494] +
      std::complex<double> (0, 1) * amp[495] - amp[497] + std::complex<double>
      (0, 1) * amp[498] + std::complex<double> (0, 1) * amp[499] - amp[500] +
      std::complex<double> (0, 1) * amp[501] - amp[503] + amp[516] + amp[517] +
      amp[519] + amp[522] + amp[523] + amp[525] - amp[528] + amp[530] +
      amp[533] - amp[531] - amp[534] + amp[536] + amp[539] - amp[537] -
      std::complex<double> (0, 1) * amp[540] - amp[541] - std::complex<double>
      (0, 1) * amp[542] - std::complex<double> (0, 1) * amp[544] -
      std::complex<double> (0, 1) * amp[546] - amp[547] + std::complex<double>
      (0, 1) * amp[554] - std::complex<double> (0, 1) * amp[552] +
      std::complex<double> (0, 1) * amp[557] - std::complex<double> (0, 1) *
      amp[555] - std::complex<double> (0, 1) * amp[558] - amp[559] -
      std::complex<double> (0, 1) * amp[560] - std::complex<double> (0, 1) *
      amp[562] - std::complex<double> (0, 1) * amp[564] - amp[565] +
      std::complex<double> (0, 1) * amp[572] - std::complex<double> (0, 1) *
      amp[570] + std::complex<double> (0, 1) * amp[575] - std::complex<double>
      (0, 1) * amp[573] - amp[576] - amp[577] + std::complex<double> (0, 1) *
      amp[580] - amp[582] - amp[583] + std::complex<double> (0, 1) * amp[586] -
      std::complex<double> (0, 1) * amp[760] - std::complex<double> (0, 1) *
      amp[762] - amp[763] + std::complex<double> (0, 1) * amp[766] +
      std::complex<double> (0, 1) * amp[770] - std::complex<double> (0, 1) *
      amp[768] + std::complex<double> (0, 1) * amp[773] - std::complex<double>
      (0, 1) * amp[771] - std::complex<double> (0, 1) * amp[778] -
      std::complex<double> (0, 1) * amp[780] - amp[781] + std::complex<double>
      (0, 1) * amp[784] + std::complex<double> (0, 1) * amp[788] -
      std::complex<double> (0, 1) * amp[786] + std::complex<double> (0, 1) *
      amp[791] - std::complex<double> (0, 1) * amp[789] + amp[844] + amp[850] +
      std::complex<double> (0, 1) * amp[866] - std::complex<double> (0, 1) *
      amp[864] + std::complex<double> (0, 1) * amp[869] - std::complex<double>
      (0, 1) * amp[867] + amp[875] - amp[873] + std::complex<double> (0, 1) *
      amp[878] - std::complex<double> (0, 1) * amp[876] + std::complex<double>
      (0, 1) * amp[881] - std::complex<double> (0, 1) * amp[879] + amp[887] -
      amp[885] - amp[888] + amp[890] + amp[893] - amp[891] - amp[894] +
      amp[896] + amp[899] - amp[897] + std::complex<double> (0, 1) * amp[1010]
      - std::complex<double> (0, 1) * amp[1008] + std::complex<double> (0, 1) *
      amp[1013] - std::complex<double> (0, 1) * amp[1011] + amp[1016] -
      amp[1014] + std::complex<double> (0, 1) * amp[1022] -
      std::complex<double> (0, 1) * amp[1020] + std::complex<double> (0, 1) *
      amp[1025] - std::complex<double> (0, 1) * amp[1023] + amp[1028] -
      amp[1026];

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



