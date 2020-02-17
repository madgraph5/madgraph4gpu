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
// Process: u u > ta+ ta- u u u u~ WEIGHTED<=8 / h
// Process: c c > ta+ ta- c c c c~ WEIGHTED<=8 / h

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
  const int denominators[nprocesses] = {216}; 

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
        t[0] = matrix_uu_taptamuuuux_no_h(); 

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
      t[0] = matrix_uu_taptamuuuux_no_h(); 

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
  if(id1 == 4 && id2 == 4)
  {
    // Add matrix elements for processes with beams (4, 4)
    return matrix_element[0]; 
  }
  else if(id1 == 2 && id2 == 2)
  {
    // Add matrix elements for processes with beams (2, 2)
    return matrix_element[0]; 
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
  oxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  oxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]); 
  oxxxxx(p[perm[6]], mME[6], hel[6], +1, w[6]); 
  ixxxxx(p[perm[7]], mME[7], hel[7], -1, w[7]); 
  FFV1P0_3(w[0], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[8]); 
  FFV1P0_3(w[1], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[9]); 
  FFV1P0_3(w[2], w[3], pars->GC_3, pars->ZERO, pars->ZERO, w[10]); 
  VVV1P0_1(w[8], w[9], pars->GC_10, pars->ZERO, pars->ZERO, w[11]); 
  FFV1_1(w[6], w[10], pars->GC_2, pars->ZERO, pars->ZERO, w[12]); 
  FFV1_2(w[7], w[10], pars->GC_2, pars->ZERO, pars->ZERO, w[13]); 
  FFV1_1(w[6], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[14]); 
  FFV1_2(w[7], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[15]); 
  FFV1_2(w[7], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[16]); 
  FFV1_1(w[6], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[17]); 
  FFV2_4_3(w[2], w[3], -pars->GC_51, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ,
      w[18]);
  FFV2_5_1(w[6], w[18], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[19]);
  FFV2_5_2(w[7], w[18], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[20]);
  FFV1P0_3(w[1], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[21]); 
  VVV1P0_1(w[8], w[21], pars->GC_10, pars->ZERO, pars->ZERO, w[22]); 
  FFV1_1(w[5], w[10], pars->GC_2, pars->ZERO, pars->ZERO, w[23]); 
  FFV1_1(w[5], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[24]); 
  FFV1_2(w[7], w[21], pars->GC_11, pars->ZERO, pars->ZERO, w[25]); 
  FFV1_1(w[5], w[21], pars->GC_11, pars->ZERO, pars->ZERO, w[26]); 
  FFV2_5_1(w[5], w[18], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[27]);
  FFV1_2(w[1], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[28]); 
  FFV1P0_3(w[28], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[29]); 
  FFV1P0_3(w[28], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[30]); 
  FFV1_2(w[1], w[10], pars->GC_2, pars->ZERO, pars->ZERO, w[31]); 
  FFV1P0_3(w[7], w[24], pars->GC_11, pars->ZERO, pars->ZERO, w[32]); 
  FFV1P0_3(w[1], w[24], pars->GC_11, pars->ZERO, pars->ZERO, w[33]); 
  FFV1P0_3(w[7], w[14], pars->GC_11, pars->ZERO, pars->ZERO, w[34]); 
  FFV1P0_3(w[1], w[14], pars->GC_11, pars->ZERO, pars->ZERO, w[35]); 
  FFV1P0_3(w[16], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[36]); 
  FFV1P0_3(w[16], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[37]); 
  FFV2_5_2(w[1], w[18], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[38]);
  FFV1P0_3(w[7], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[39]); 
  FFV1_1(w[6], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[40]); 
  VVV1P0_1(w[8], w[39], pars->GC_10, pars->ZERO, pars->ZERO, w[41]); 
  FFV1_2(w[1], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[42]); 
  FFV1P0_3(w[7], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[43]); 
  FFV1_1(w[5], w[43], pars->GC_11, pars->ZERO, pars->ZERO, w[44]); 
  FFV1_2(w[1], w[43], pars->GC_11, pars->ZERO, pars->ZERO, w[45]); 
  VVV1P0_1(w[8], w[43], pars->GC_10, pars->ZERO, pars->ZERO, w[46]); 
  FFV1P0_3(w[0], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[47]); 
  FFV1P0_3(w[1], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[48]); 
  VVV1P0_1(w[47], w[48], pars->GC_10, pars->ZERO, pars->ZERO, w[49]); 
  FFV1_1(w[6], w[47], pars->GC_11, pars->ZERO, pars->ZERO, w[50]); 
  FFV1_2(w[7], w[48], pars->GC_11, pars->ZERO, pars->ZERO, w[51]); 
  FFV1_2(w[7], w[47], pars->GC_11, pars->ZERO, pars->ZERO, w[52]); 
  FFV1_1(w[6], w[48], pars->GC_11, pars->ZERO, pars->ZERO, w[53]); 
  VVV1P0_1(w[47], w[21], pars->GC_10, pars->ZERO, pars->ZERO, w[54]); 
  FFV1_1(w[4], w[10], pars->GC_2, pars->ZERO, pars->ZERO, w[55]); 
  FFV1_1(w[4], w[47], pars->GC_11, pars->ZERO, pars->ZERO, w[56]); 
  FFV1_1(w[4], w[21], pars->GC_11, pars->ZERO, pars->ZERO, w[57]); 
  FFV2_5_1(w[4], w[18], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[58]);
  FFV1_2(w[1], w[47], pars->GC_11, pars->ZERO, pars->ZERO, w[59]); 
  FFV1P0_3(w[59], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[60]); 
  FFV1P0_3(w[59], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[61]); 
  FFV1P0_3(w[7], w[56], pars->GC_11, pars->ZERO, pars->ZERO, w[62]); 
  FFV1P0_3(w[1], w[56], pars->GC_11, pars->ZERO, pars->ZERO, w[63]); 
  FFV1P0_3(w[7], w[50], pars->GC_11, pars->ZERO, pars->ZERO, w[64]); 
  FFV1P0_3(w[1], w[50], pars->GC_11, pars->ZERO, pars->ZERO, w[65]); 
  FFV1P0_3(w[52], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[66]); 
  FFV1P0_3(w[52], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[67]); 
  FFV1P0_3(w[7], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[68]); 
  FFV1_1(w[6], w[68], pars->GC_11, pars->ZERO, pars->ZERO, w[69]); 
  VVV1P0_1(w[47], w[68], pars->GC_10, pars->ZERO, pars->ZERO, w[70]); 
  FFV1_2(w[1], w[68], pars->GC_11, pars->ZERO, pars->ZERO, w[71]); 
  FFV1_1(w[4], w[43], pars->GC_11, pars->ZERO, pars->ZERO, w[72]); 
  VVV1P0_1(w[47], w[43], pars->GC_10, pars->ZERO, pars->ZERO, w[73]); 
  FFV1P0_3(w[0], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[74]); 
  VVV1P0_1(w[74], w[48], pars->GC_10, pars->ZERO, pars->ZERO, w[75]); 
  FFV1_1(w[5], w[74], pars->GC_11, pars->ZERO, pars->ZERO, w[76]); 
  FFV1_2(w[7], w[74], pars->GC_11, pars->ZERO, pars->ZERO, w[77]); 
  FFV1_1(w[5], w[48], pars->GC_11, pars->ZERO, pars->ZERO, w[78]); 
  VVV1P0_1(w[74], w[9], pars->GC_10, pars->ZERO, pars->ZERO, w[79]); 
  FFV1_1(w[4], w[74], pars->GC_11, pars->ZERO, pars->ZERO, w[80]); 
  FFV1_1(w[4], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[81]); 
  FFV1_2(w[1], w[74], pars->GC_11, pars->ZERO, pars->ZERO, w[82]); 
  FFV1P0_3(w[82], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[83]); 
  FFV1P0_3(w[82], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[84]); 
  FFV1P0_3(w[7], w[80], pars->GC_11, pars->ZERO, pars->ZERO, w[85]); 
  FFV1P0_3(w[1], w[80], pars->GC_11, pars->ZERO, pars->ZERO, w[86]); 
  FFV1P0_3(w[7], w[76], pars->GC_11, pars->ZERO, pars->ZERO, w[87]); 
  FFV1P0_3(w[1], w[76], pars->GC_11, pars->ZERO, pars->ZERO, w[88]); 
  FFV1P0_3(w[77], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[89]); 
  FFV1P0_3(w[77], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[90]); 
  FFV1_1(w[5], w[68], pars->GC_11, pars->ZERO, pars->ZERO, w[91]); 
  VVV1P0_1(w[74], w[68], pars->GC_10, pars->ZERO, pars->ZERO, w[92]); 
  FFV1_1(w[4], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[93]); 
  VVV1P0_1(w[74], w[39], pars->GC_10, pars->ZERO, pars->ZERO, w[94]); 
  FFV1_2(w[0], w[48], pars->GC_11, pars->ZERO, pars->ZERO, w[95]); 
  FFV1P0_3(w[95], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[96]); 
  FFV1P0_3(w[95], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[97]); 
  FFV1_2(w[0], w[10], pars->GC_2, pars->ZERO, pars->ZERO, w[98]); 
  FFV1P0_3(w[98], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[99]); 
  FFV1P0_3(w[98], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[100]); 
  FFV1P0_3(w[0], w[78], pars->GC_11, pars->ZERO, pars->ZERO, w[101]); 
  FFV1P0_3(w[0], w[12], pars->GC_11, pars->ZERO, pars->ZERO, w[102]); 
  FFV1P0_3(w[0], w[53], pars->GC_11, pars->ZERO, pars->ZERO, w[103]); 
  FFV1P0_3(w[0], w[23], pars->GC_11, pars->ZERO, pars->ZERO, w[104]); 
  FFV2_5_2(w[0], w[18], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[105]);
  FFV1P0_3(w[105], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[106]); 
  FFV1P0_3(w[105], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[107]); 
  FFV1P0_3(w[0], w[19], pars->GC_11, pars->ZERO, pars->ZERO, w[108]); 
  FFV1P0_3(w[0], w[27], pars->GC_11, pars->ZERO, pars->ZERO, w[109]); 
  VVV1P0_1(w[48], w[39], pars->GC_10, pars->ZERO, pars->ZERO, w[110]); 
  FFV1_2(w[0], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[111]); 
  VVV1P0_1(w[48], w[43], pars->GC_10, pars->ZERO, pars->ZERO, w[112]); 
  FFV1_2(w[0], w[43], pars->GC_11, pars->ZERO, pars->ZERO, w[113]); 
  FFV1_2(w[0], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[114]); 
  FFV1P0_3(w[114], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[115]); 
  FFV1P0_3(w[114], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[116]); 
  FFV1P0_3(w[98], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[117]); 
  FFV1P0_3(w[0], w[81], pars->GC_11, pars->ZERO, pars->ZERO, w[118]); 
  FFV1P0_3(w[0], w[17], pars->GC_11, pars->ZERO, pars->ZERO, w[119]); 
  FFV1P0_3(w[0], w[55], pars->GC_11, pars->ZERO, pars->ZERO, w[120]); 
  FFV1P0_3(w[105], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[121]); 
  FFV1P0_3(w[0], w[58], pars->GC_11, pars->ZERO, pars->ZERO, w[122]); 
  VVV1P0_1(w[9], w[68], pars->GC_10, pars->ZERO, pars->ZERO, w[123]); 
  FFV1_2(w[0], w[68], pars->GC_11, pars->ZERO, pars->ZERO, w[124]); 
  VVV1P0_1(w[9], w[43], pars->GC_10, pars->ZERO, pars->ZERO, w[125]); 
  FFV1_2(w[0], w[21], pars->GC_11, pars->ZERO, pars->ZERO, w[126]); 
  FFV1P0_3(w[126], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[127]); 
  FFV1P0_3(w[126], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[128]); 
  FFV1P0_3(w[0], w[57], pars->GC_11, pars->ZERO, pars->ZERO, w[129]); 
  FFV1P0_3(w[0], w[26], pars->GC_11, pars->ZERO, pars->ZERO, w[130]); 
  VVV1P0_1(w[21], w[68], pars->GC_10, pars->ZERO, pars->ZERO, w[131]); 
  VVV1P0_1(w[21], w[39], pars->GC_10, pars->ZERO, pars->ZERO, w[132]); 
  FFV1P0_3(w[124], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[133]); 
  FFV1P0_3(w[124], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[134]); 
  FFV1P0_3(w[0], w[91], pars->GC_11, pars->ZERO, pars->ZERO, w[135]); 
  FFV1P0_3(w[0], w[69], pars->GC_11, pars->ZERO, pars->ZERO, w[136]); 
  FFV1P0_3(w[111], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[137]); 
  FFV1P0_3(w[111], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[138]); 
  FFV1P0_3(w[0], w[93], pars->GC_11, pars->ZERO, pars->ZERO, w[139]); 
  FFV1P0_3(w[0], w[40], pars->GC_11, pars->ZERO, pars->ZERO, w[140]); 
  FFV1P0_3(w[113], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[141]); 
  FFV1P0_3(w[113], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[142]); 
  FFV1P0_3(w[0], w[72], pars->GC_11, pars->ZERO, pars->ZERO, w[143]); 
  FFV1P0_3(w[0], w[44], pars->GC_11, pars->ZERO, pars->ZERO, w[144]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[7], w[12], w[11], pars->GC_11, amp[0]); 
  FFV1_0(w[13], w[6], w[11], pars->GC_11, amp[1]); 
  FFV1_0(w[15], w[14], w[10], pars->GC_2, amp[2]); 
  FFV1_0(w[13], w[14], w[9], pars->GC_11, amp[3]); 
  FFV1_0(w[16], w[17], w[10], pars->GC_2, amp[4]); 
  FFV1_0(w[16], w[12], w[9], pars->GC_11, amp[5]); 
  FFV1_0(w[13], w[17], w[8], pars->GC_11, amp[6]); 
  FFV1_0(w[15], w[12], w[8], pars->GC_11, amp[7]); 
  FFV1_0(w[7], w[19], w[11], pars->GC_11, amp[8]); 
  FFV1_0(w[20], w[6], w[11], pars->GC_11, amp[9]); 
  FFV2_5_0(w[15], w[14], w[18], pars->GC_51, pars->GC_58, amp[10]); 
  FFV1_0(w[20], w[14], w[9], pars->GC_11, amp[11]); 
  FFV2_5_0(w[16], w[17], w[18], pars->GC_51, pars->GC_58, amp[12]); 
  FFV1_0(w[16], w[19], w[9], pars->GC_11, amp[13]); 
  FFV1_0(w[20], w[17], w[8], pars->GC_11, amp[14]); 
  FFV1_0(w[15], w[19], w[8], pars->GC_11, amp[15]); 
  FFV1_0(w[7], w[23], w[22], pars->GC_11, amp[16]); 
  FFV1_0(w[13], w[5], w[22], pars->GC_11, amp[17]); 
  FFV1_0(w[25], w[24], w[10], pars->GC_2, amp[18]); 
  FFV1_0(w[13], w[24], w[21], pars->GC_11, amp[19]); 
  FFV1_0(w[16], w[26], w[10], pars->GC_2, amp[20]); 
  FFV1_0(w[16], w[23], w[21], pars->GC_11, amp[21]); 
  FFV1_0(w[13], w[26], w[8], pars->GC_11, amp[22]); 
  FFV1_0(w[25], w[23], w[8], pars->GC_11, amp[23]); 
  FFV1_0(w[7], w[27], w[22], pars->GC_11, amp[24]); 
  FFV1_0(w[20], w[5], w[22], pars->GC_11, amp[25]); 
  FFV2_5_0(w[25], w[24], w[18], pars->GC_51, pars->GC_58, amp[26]); 
  FFV1_0(w[20], w[24], w[21], pars->GC_11, amp[27]); 
  FFV2_5_0(w[16], w[26], w[18], pars->GC_51, pars->GC_58, amp[28]); 
  FFV1_0(w[16], w[27], w[21], pars->GC_11, amp[29]); 
  FFV1_0(w[20], w[26], w[8], pars->GC_11, amp[30]); 
  FFV1_0(w[25], w[27], w[8], pars->GC_11, amp[31]); 
  FFV1_0(w[7], w[23], w[29], pars->GC_11, amp[32]); 
  FFV1_0(w[7], w[12], w[30], pars->GC_11, amp[33]); 
  FFV1_0(w[13], w[6], w[30], pars->GC_11, amp[34]); 
  FFV1_0(w[13], w[5], w[29], pars->GC_11, amp[35]); 
  FFV1_0(w[31], w[6], w[32], pars->GC_11, amp[36]); 
  FFV1_0(w[7], w[12], w[33], pars->GC_11, amp[37]); 
  FFV1_0(w[1], w[12], w[32], pars->GC_11, amp[38]); 
  FFV1_0(w[13], w[6], w[33], pars->GC_11, amp[39]); 
  FFV1_0(w[31], w[5], w[34], pars->GC_11, amp[40]); 
  FFV1_0(w[7], w[23], w[35], pars->GC_11, amp[41]); 
  FFV1_0(w[1], w[23], w[34], pars->GC_11, amp[42]); 
  FFV1_0(w[13], w[5], w[35], pars->GC_11, amp[43]); 
  FFV1_0(w[31], w[6], w[36], pars->GC_11, amp[44]); 
  FFV1_0(w[31], w[5], w[37], pars->GC_11, amp[45]); 
  FFV1_0(w[1], w[23], w[37], pars->GC_11, amp[46]); 
  FFV1_0(w[1], w[12], w[36], pars->GC_11, amp[47]); 
  FFV1_0(w[7], w[27], w[29], pars->GC_11, amp[48]); 
  FFV1_0(w[7], w[19], w[30], pars->GC_11, amp[49]); 
  FFV1_0(w[20], w[6], w[30], pars->GC_11, amp[50]); 
  FFV1_0(w[20], w[5], w[29], pars->GC_11, amp[51]); 
  FFV1_0(w[38], w[6], w[32], pars->GC_11, amp[52]); 
  FFV1_0(w[7], w[19], w[33], pars->GC_11, amp[53]); 
  FFV1_0(w[1], w[19], w[32], pars->GC_11, amp[54]); 
  FFV1_0(w[20], w[6], w[33], pars->GC_11, amp[55]); 
  FFV1_0(w[38], w[5], w[34], pars->GC_11, amp[56]); 
  FFV1_0(w[7], w[27], w[35], pars->GC_11, amp[57]); 
  FFV1_0(w[1], w[27], w[34], pars->GC_11, amp[58]); 
  FFV1_0(w[20], w[5], w[35], pars->GC_11, amp[59]); 
  FFV1_0(w[38], w[6], w[36], pars->GC_11, amp[60]); 
  FFV1_0(w[38], w[5], w[37], pars->GC_11, amp[61]); 
  FFV1_0(w[1], w[27], w[37], pars->GC_11, amp[62]); 
  FFV1_0(w[1], w[19], w[36], pars->GC_11, amp[63]); 
  FFV1_0(w[28], w[12], w[39], pars->GC_11, amp[64]); 
  FFV1_0(w[28], w[40], w[10], pars->GC_2, amp[65]); 
  FFV1_0(w[31], w[6], w[41], pars->GC_11, amp[66]); 
  FFV1_0(w[1], w[12], w[41], pars->GC_11, amp[67]); 
  FFV1_0(w[31], w[14], w[39], pars->GC_11, amp[68]); 
  FFV1_0(w[42], w[14], w[10], pars->GC_2, amp[69]); 
  FFV1_0(w[31], w[40], w[8], pars->GC_11, amp[70]); 
  FFV1_0(w[42], w[12], w[8], pars->GC_11, amp[71]); 
  FFV1_0(w[28], w[19], w[39], pars->GC_11, amp[72]); 
  FFV2_5_0(w[28], w[40], w[18], pars->GC_51, pars->GC_58, amp[73]); 
  FFV1_0(w[38], w[6], w[41], pars->GC_11, amp[74]); 
  FFV1_0(w[1], w[19], w[41], pars->GC_11, amp[75]); 
  FFV1_0(w[38], w[14], w[39], pars->GC_11, amp[76]); 
  FFV2_5_0(w[42], w[14], w[18], pars->GC_51, pars->GC_58, amp[77]); 
  FFV1_0(w[38], w[40], w[8], pars->GC_11, amp[78]); 
  FFV1_0(w[42], w[19], w[8], pars->GC_11, amp[79]); 
  FFV1_0(w[28], w[23], w[43], pars->GC_11, amp[80]); 
  FFV1_0(w[28], w[44], w[10], pars->GC_2, amp[81]); 
  FFV1_0(w[31], w[24], w[43], pars->GC_11, amp[82]); 
  FFV1_0(w[45], w[24], w[10], pars->GC_2, amp[83]); 
  FFV1_0(w[31], w[5], w[46], pars->GC_11, amp[84]); 
  FFV1_0(w[1], w[23], w[46], pars->GC_11, amp[85]); 
  FFV1_0(w[31], w[44], w[8], pars->GC_11, amp[86]); 
  FFV1_0(w[45], w[23], w[8], pars->GC_11, amp[87]); 
  FFV1_0(w[28], w[27], w[43], pars->GC_11, amp[88]); 
  FFV2_5_0(w[28], w[44], w[18], pars->GC_51, pars->GC_58, amp[89]); 
  FFV1_0(w[38], w[24], w[43], pars->GC_11, amp[90]); 
  FFV2_5_0(w[45], w[24], w[18], pars->GC_51, pars->GC_58, amp[91]); 
  FFV1_0(w[38], w[5], w[46], pars->GC_11, amp[92]); 
  FFV1_0(w[1], w[27], w[46], pars->GC_11, amp[93]); 
  FFV1_0(w[38], w[44], w[8], pars->GC_11, amp[94]); 
  FFV1_0(w[45], w[27], w[8], pars->GC_11, amp[95]); 
  FFV1_0(w[7], w[12], w[49], pars->GC_11, amp[96]); 
  FFV1_0(w[13], w[6], w[49], pars->GC_11, amp[97]); 
  FFV1_0(w[51], w[50], w[10], pars->GC_2, amp[98]); 
  FFV1_0(w[13], w[50], w[48], pars->GC_11, amp[99]); 
  FFV1_0(w[52], w[53], w[10], pars->GC_2, amp[100]); 
  FFV1_0(w[52], w[12], w[48], pars->GC_11, amp[101]); 
  FFV1_0(w[13], w[53], w[47], pars->GC_11, amp[102]); 
  FFV1_0(w[51], w[12], w[47], pars->GC_11, amp[103]); 
  FFV1_0(w[7], w[19], w[49], pars->GC_11, amp[104]); 
  FFV1_0(w[20], w[6], w[49], pars->GC_11, amp[105]); 
  FFV2_5_0(w[51], w[50], w[18], pars->GC_51, pars->GC_58, amp[106]); 
  FFV1_0(w[20], w[50], w[48], pars->GC_11, amp[107]); 
  FFV2_5_0(w[52], w[53], w[18], pars->GC_51, pars->GC_58, amp[108]); 
  FFV1_0(w[52], w[19], w[48], pars->GC_11, amp[109]); 
  FFV1_0(w[20], w[53], w[47], pars->GC_11, amp[110]); 
  FFV1_0(w[51], w[19], w[47], pars->GC_11, amp[111]); 
  FFV1_0(w[7], w[55], w[54], pars->GC_11, amp[112]); 
  FFV1_0(w[13], w[4], w[54], pars->GC_11, amp[113]); 
  FFV1_0(w[25], w[56], w[10], pars->GC_2, amp[114]); 
  FFV1_0(w[13], w[56], w[21], pars->GC_11, amp[115]); 
  FFV1_0(w[52], w[57], w[10], pars->GC_2, amp[116]); 
  FFV1_0(w[52], w[55], w[21], pars->GC_11, amp[117]); 
  FFV1_0(w[13], w[57], w[47], pars->GC_11, amp[118]); 
  FFV1_0(w[25], w[55], w[47], pars->GC_11, amp[119]); 
  FFV1_0(w[7], w[58], w[54], pars->GC_11, amp[120]); 
  FFV1_0(w[20], w[4], w[54], pars->GC_11, amp[121]); 
  FFV2_5_0(w[25], w[56], w[18], pars->GC_51, pars->GC_58, amp[122]); 
  FFV1_0(w[20], w[56], w[21], pars->GC_11, amp[123]); 
  FFV2_5_0(w[52], w[57], w[18], pars->GC_51, pars->GC_58, amp[124]); 
  FFV1_0(w[52], w[58], w[21], pars->GC_11, amp[125]); 
  FFV1_0(w[20], w[57], w[47], pars->GC_11, amp[126]); 
  FFV1_0(w[25], w[58], w[47], pars->GC_11, amp[127]); 
  FFV1_0(w[7], w[55], w[60], pars->GC_11, amp[128]); 
  FFV1_0(w[7], w[12], w[61], pars->GC_11, amp[129]); 
  FFV1_0(w[13], w[6], w[61], pars->GC_11, amp[130]); 
  FFV1_0(w[13], w[4], w[60], pars->GC_11, amp[131]); 
  FFV1_0(w[31], w[6], w[62], pars->GC_11, amp[132]); 
  FFV1_0(w[7], w[12], w[63], pars->GC_11, amp[133]); 
  FFV1_0(w[1], w[12], w[62], pars->GC_11, amp[134]); 
  FFV1_0(w[13], w[6], w[63], pars->GC_11, amp[135]); 
  FFV1_0(w[31], w[4], w[64], pars->GC_11, amp[136]); 
  FFV1_0(w[7], w[55], w[65], pars->GC_11, amp[137]); 
  FFV1_0(w[1], w[55], w[64], pars->GC_11, amp[138]); 
  FFV1_0(w[13], w[4], w[65], pars->GC_11, amp[139]); 
  FFV1_0(w[31], w[6], w[66], pars->GC_11, amp[140]); 
  FFV1_0(w[31], w[4], w[67], pars->GC_11, amp[141]); 
  FFV1_0(w[1], w[55], w[67], pars->GC_11, amp[142]); 
  FFV1_0(w[1], w[12], w[66], pars->GC_11, amp[143]); 
  FFV1_0(w[7], w[58], w[60], pars->GC_11, amp[144]); 
  FFV1_0(w[7], w[19], w[61], pars->GC_11, amp[145]); 
  FFV1_0(w[20], w[6], w[61], pars->GC_11, amp[146]); 
  FFV1_0(w[20], w[4], w[60], pars->GC_11, amp[147]); 
  FFV1_0(w[38], w[6], w[62], pars->GC_11, amp[148]); 
  FFV1_0(w[7], w[19], w[63], pars->GC_11, amp[149]); 
  FFV1_0(w[1], w[19], w[62], pars->GC_11, amp[150]); 
  FFV1_0(w[20], w[6], w[63], pars->GC_11, amp[151]); 
  FFV1_0(w[38], w[4], w[64], pars->GC_11, amp[152]); 
  FFV1_0(w[7], w[58], w[65], pars->GC_11, amp[153]); 
  FFV1_0(w[1], w[58], w[64], pars->GC_11, amp[154]); 
  FFV1_0(w[20], w[4], w[65], pars->GC_11, amp[155]); 
  FFV1_0(w[38], w[6], w[66], pars->GC_11, amp[156]); 
  FFV1_0(w[38], w[4], w[67], pars->GC_11, amp[157]); 
  FFV1_0(w[1], w[58], w[67], pars->GC_11, amp[158]); 
  FFV1_0(w[1], w[19], w[66], pars->GC_11, amp[159]); 
  FFV1_0(w[59], w[12], w[68], pars->GC_11, amp[160]); 
  FFV1_0(w[59], w[69], w[10], pars->GC_2, amp[161]); 
  FFV1_0(w[31], w[6], w[70], pars->GC_11, amp[162]); 
  FFV1_0(w[1], w[12], w[70], pars->GC_11, amp[163]); 
  FFV1_0(w[31], w[50], w[68], pars->GC_11, amp[164]); 
  FFV1_0(w[71], w[50], w[10], pars->GC_2, amp[165]); 
  FFV1_0(w[31], w[69], w[47], pars->GC_11, amp[166]); 
  FFV1_0(w[71], w[12], w[47], pars->GC_11, amp[167]); 
  FFV1_0(w[59], w[19], w[68], pars->GC_11, amp[168]); 
  FFV2_5_0(w[59], w[69], w[18], pars->GC_51, pars->GC_58, amp[169]); 
  FFV1_0(w[38], w[6], w[70], pars->GC_11, amp[170]); 
  FFV1_0(w[1], w[19], w[70], pars->GC_11, amp[171]); 
  FFV1_0(w[38], w[50], w[68], pars->GC_11, amp[172]); 
  FFV2_5_0(w[71], w[50], w[18], pars->GC_51, pars->GC_58, amp[173]); 
  FFV1_0(w[38], w[69], w[47], pars->GC_11, amp[174]); 
  FFV1_0(w[71], w[19], w[47], pars->GC_11, amp[175]); 
  FFV1_0(w[59], w[55], w[43], pars->GC_11, amp[176]); 
  FFV1_0(w[59], w[72], w[10], pars->GC_2, amp[177]); 
  FFV1_0(w[31], w[56], w[43], pars->GC_11, amp[178]); 
  FFV1_0(w[45], w[56], w[10], pars->GC_2, amp[179]); 
  FFV1_0(w[31], w[4], w[73], pars->GC_11, amp[180]); 
  FFV1_0(w[1], w[55], w[73], pars->GC_11, amp[181]); 
  FFV1_0(w[31], w[72], w[47], pars->GC_11, amp[182]); 
  FFV1_0(w[45], w[55], w[47], pars->GC_11, amp[183]); 
  FFV1_0(w[59], w[58], w[43], pars->GC_11, amp[184]); 
  FFV2_5_0(w[59], w[72], w[18], pars->GC_51, pars->GC_58, amp[185]); 
  FFV1_0(w[38], w[56], w[43], pars->GC_11, amp[186]); 
  FFV2_5_0(w[45], w[56], w[18], pars->GC_51, pars->GC_58, amp[187]); 
  FFV1_0(w[38], w[4], w[73], pars->GC_11, amp[188]); 
  FFV1_0(w[1], w[58], w[73], pars->GC_11, amp[189]); 
  FFV1_0(w[38], w[72], w[47], pars->GC_11, amp[190]); 
  FFV1_0(w[45], w[58], w[47], pars->GC_11, amp[191]); 
  FFV1_0(w[7], w[23], w[75], pars->GC_11, amp[192]); 
  FFV1_0(w[13], w[5], w[75], pars->GC_11, amp[193]); 
  FFV1_0(w[51], w[76], w[10], pars->GC_2, amp[194]); 
  FFV1_0(w[13], w[76], w[48], pars->GC_11, amp[195]); 
  FFV1_0(w[77], w[78], w[10], pars->GC_2, amp[196]); 
  FFV1_0(w[77], w[23], w[48], pars->GC_11, amp[197]); 
  FFV1_0(w[13], w[78], w[74], pars->GC_11, amp[198]); 
  FFV1_0(w[51], w[23], w[74], pars->GC_11, amp[199]); 
  FFV1_0(w[7], w[27], w[75], pars->GC_11, amp[200]); 
  FFV1_0(w[20], w[5], w[75], pars->GC_11, amp[201]); 
  FFV2_5_0(w[51], w[76], w[18], pars->GC_51, pars->GC_58, amp[202]); 
  FFV1_0(w[20], w[76], w[48], pars->GC_11, amp[203]); 
  FFV2_5_0(w[77], w[78], w[18], pars->GC_51, pars->GC_58, amp[204]); 
  FFV1_0(w[77], w[27], w[48], pars->GC_11, amp[205]); 
  FFV1_0(w[20], w[78], w[74], pars->GC_11, amp[206]); 
  FFV1_0(w[51], w[27], w[74], pars->GC_11, amp[207]); 
  FFV1_0(w[7], w[55], w[79], pars->GC_11, amp[208]); 
  FFV1_0(w[13], w[4], w[79], pars->GC_11, amp[209]); 
  FFV1_0(w[15], w[80], w[10], pars->GC_2, amp[210]); 
  FFV1_0(w[13], w[80], w[9], pars->GC_11, amp[211]); 
  FFV1_0(w[77], w[81], w[10], pars->GC_2, amp[212]); 
  FFV1_0(w[77], w[55], w[9], pars->GC_11, amp[213]); 
  FFV1_0(w[13], w[81], w[74], pars->GC_11, amp[214]); 
  FFV1_0(w[15], w[55], w[74], pars->GC_11, amp[215]); 
  FFV1_0(w[7], w[58], w[79], pars->GC_11, amp[216]); 
  FFV1_0(w[20], w[4], w[79], pars->GC_11, amp[217]); 
  FFV2_5_0(w[15], w[80], w[18], pars->GC_51, pars->GC_58, amp[218]); 
  FFV1_0(w[20], w[80], w[9], pars->GC_11, amp[219]); 
  FFV2_5_0(w[77], w[81], w[18], pars->GC_51, pars->GC_58, amp[220]); 
  FFV1_0(w[77], w[58], w[9], pars->GC_11, amp[221]); 
  FFV1_0(w[20], w[81], w[74], pars->GC_11, amp[222]); 
  FFV1_0(w[15], w[58], w[74], pars->GC_11, amp[223]); 
  FFV1_0(w[7], w[55], w[83], pars->GC_11, amp[224]); 
  FFV1_0(w[7], w[23], w[84], pars->GC_11, amp[225]); 
  FFV1_0(w[13], w[5], w[84], pars->GC_11, amp[226]); 
  FFV1_0(w[13], w[4], w[83], pars->GC_11, amp[227]); 
  FFV1_0(w[31], w[5], w[85], pars->GC_11, amp[228]); 
  FFV1_0(w[7], w[23], w[86], pars->GC_11, amp[229]); 
  FFV1_0(w[1], w[23], w[85], pars->GC_11, amp[230]); 
  FFV1_0(w[13], w[5], w[86], pars->GC_11, amp[231]); 
  FFV1_0(w[31], w[4], w[87], pars->GC_11, amp[232]); 
  FFV1_0(w[7], w[55], w[88], pars->GC_11, amp[233]); 
  FFV1_0(w[1], w[55], w[87], pars->GC_11, amp[234]); 
  FFV1_0(w[13], w[4], w[88], pars->GC_11, amp[235]); 
  FFV1_0(w[31], w[5], w[89], pars->GC_11, amp[236]); 
  FFV1_0(w[31], w[4], w[90], pars->GC_11, amp[237]); 
  FFV1_0(w[1], w[55], w[90], pars->GC_11, amp[238]); 
  FFV1_0(w[1], w[23], w[89], pars->GC_11, amp[239]); 
  FFV1_0(w[7], w[58], w[83], pars->GC_11, amp[240]); 
  FFV1_0(w[7], w[27], w[84], pars->GC_11, amp[241]); 
  FFV1_0(w[20], w[5], w[84], pars->GC_11, amp[242]); 
  FFV1_0(w[20], w[4], w[83], pars->GC_11, amp[243]); 
  FFV1_0(w[38], w[5], w[85], pars->GC_11, amp[244]); 
  FFV1_0(w[7], w[27], w[86], pars->GC_11, amp[245]); 
  FFV1_0(w[1], w[27], w[85], pars->GC_11, amp[246]); 
  FFV1_0(w[20], w[5], w[86], pars->GC_11, amp[247]); 
  FFV1_0(w[38], w[4], w[87], pars->GC_11, amp[248]); 
  FFV1_0(w[7], w[58], w[88], pars->GC_11, amp[249]); 
  FFV1_0(w[1], w[58], w[87], pars->GC_11, amp[250]); 
  FFV1_0(w[20], w[4], w[88], pars->GC_11, amp[251]); 
  FFV1_0(w[38], w[5], w[89], pars->GC_11, amp[252]); 
  FFV1_0(w[38], w[4], w[90], pars->GC_11, amp[253]); 
  FFV1_0(w[1], w[58], w[90], pars->GC_11, amp[254]); 
  FFV1_0(w[1], w[27], w[89], pars->GC_11, amp[255]); 
  FFV1_0(w[82], w[23], w[68], pars->GC_11, amp[256]); 
  FFV1_0(w[82], w[91], w[10], pars->GC_2, amp[257]); 
  FFV1_0(w[31], w[5], w[92], pars->GC_11, amp[258]); 
  FFV1_0(w[1], w[23], w[92], pars->GC_11, amp[259]); 
  FFV1_0(w[31], w[76], w[68], pars->GC_11, amp[260]); 
  FFV1_0(w[71], w[76], w[10], pars->GC_2, amp[261]); 
  FFV1_0(w[31], w[91], w[74], pars->GC_11, amp[262]); 
  FFV1_0(w[71], w[23], w[74], pars->GC_11, amp[263]); 
  FFV1_0(w[82], w[27], w[68], pars->GC_11, amp[264]); 
  FFV2_5_0(w[82], w[91], w[18], pars->GC_51, pars->GC_58, amp[265]); 
  FFV1_0(w[38], w[5], w[92], pars->GC_11, amp[266]); 
  FFV1_0(w[1], w[27], w[92], pars->GC_11, amp[267]); 
  FFV1_0(w[38], w[76], w[68], pars->GC_11, amp[268]); 
  FFV2_5_0(w[71], w[76], w[18], pars->GC_51, pars->GC_58, amp[269]); 
  FFV1_0(w[38], w[91], w[74], pars->GC_11, amp[270]); 
  FFV1_0(w[71], w[27], w[74], pars->GC_11, amp[271]); 
  FFV1_0(w[82], w[55], w[39], pars->GC_11, amp[272]); 
  FFV1_0(w[82], w[93], w[10], pars->GC_2, amp[273]); 
  FFV1_0(w[31], w[80], w[39], pars->GC_11, amp[274]); 
  FFV1_0(w[42], w[80], w[10], pars->GC_2, amp[275]); 
  FFV1_0(w[31], w[4], w[94], pars->GC_11, amp[276]); 
  FFV1_0(w[1], w[55], w[94], pars->GC_11, amp[277]); 
  FFV1_0(w[31], w[93], w[74], pars->GC_11, amp[278]); 
  FFV1_0(w[42], w[55], w[74], pars->GC_11, amp[279]); 
  FFV1_0(w[82], w[58], w[39], pars->GC_11, amp[280]); 
  FFV2_5_0(w[82], w[93], w[18], pars->GC_51, pars->GC_58, amp[281]); 
  FFV1_0(w[38], w[80], w[39], pars->GC_11, amp[282]); 
  FFV2_5_0(w[42], w[80], w[18], pars->GC_51, pars->GC_58, amp[283]); 
  FFV1_0(w[38], w[4], w[94], pars->GC_11, amp[284]); 
  FFV1_0(w[1], w[58], w[94], pars->GC_11, amp[285]); 
  FFV1_0(w[38], w[93], w[74], pars->GC_11, amp[286]); 
  FFV1_0(w[42], w[58], w[74], pars->GC_11, amp[287]); 
  FFV1_0(w[7], w[23], w[96], pars->GC_11, amp[288]); 
  FFV1_0(w[7], w[12], w[97], pars->GC_11, amp[289]); 
  FFV1_0(w[13], w[6], w[97], pars->GC_11, amp[290]); 
  FFV1_0(w[13], w[5], w[96], pars->GC_11, amp[291]); 
  FFV1_0(w[7], w[78], w[99], pars->GC_11, amp[292]); 
  FFV1_0(w[7], w[53], w[100], pars->GC_11, amp[293]); 
  FFV1_0(w[51], w[6], w[100], pars->GC_11, amp[294]); 
  FFV1_0(w[51], w[5], w[99], pars->GC_11, amp[295]); 
  FFV1_0(w[7], w[12], w[101], pars->GC_11, amp[296]); 
  FFV1_0(w[7], w[78], w[102], pars->GC_11, amp[297]); 
  FFV1_0(w[13], w[6], w[101], pars->GC_11, amp[298]); 
  FFV1_0(w[7], w[23], w[103], pars->GC_11, amp[299]); 
  FFV1_0(w[7], w[53], w[104], pars->GC_11, amp[300]); 
  FFV1_0(w[13], w[5], w[103], pars->GC_11, amp[301]); 
  FFV1_0(w[51], w[6], w[104], pars->GC_11, amp[302]); 
  FFV1_0(w[51], w[5], w[102], pars->GC_11, amp[303]); 
  FFV1_0(w[7], w[27], w[96], pars->GC_11, amp[304]); 
  FFV1_0(w[7], w[19], w[97], pars->GC_11, amp[305]); 
  FFV1_0(w[20], w[6], w[97], pars->GC_11, amp[306]); 
  FFV1_0(w[20], w[5], w[96], pars->GC_11, amp[307]); 
  FFV1_0(w[7], w[78], w[106], pars->GC_11, amp[308]); 
  FFV1_0(w[7], w[53], w[107], pars->GC_11, amp[309]); 
  FFV1_0(w[51], w[6], w[107], pars->GC_11, amp[310]); 
  FFV1_0(w[51], w[5], w[106], pars->GC_11, amp[311]); 
  FFV1_0(w[7], w[19], w[101], pars->GC_11, amp[312]); 
  FFV1_0(w[7], w[78], w[108], pars->GC_11, amp[313]); 
  FFV1_0(w[20], w[6], w[101], pars->GC_11, amp[314]); 
  FFV1_0(w[7], w[27], w[103], pars->GC_11, amp[315]); 
  FFV1_0(w[7], w[53], w[109], pars->GC_11, amp[316]); 
  FFV1_0(w[20], w[5], w[103], pars->GC_11, amp[317]); 
  FFV1_0(w[51], w[6], w[109], pars->GC_11, amp[318]); 
  FFV1_0(w[51], w[5], w[108], pars->GC_11, amp[319]); 
  FFV1_0(w[95], w[12], w[39], pars->GC_11, amp[320]); 
  FFV1_0(w[95], w[40], w[10], pars->GC_2, amp[321]); 
  FFV1_0(w[98], w[6], w[110], pars->GC_11, amp[322]); 
  FFV1_0(w[98], w[53], w[39], pars->GC_11, amp[323]); 
  FFV1_0(w[98], w[40], w[48], pars->GC_11, amp[324]); 
  FFV1_0(w[111], w[53], w[10], pars->GC_2, amp[325]); 
  FFV1_0(w[111], w[12], w[48], pars->GC_11, amp[326]); 
  FFV1_0(w[0], w[12], w[110], pars->GC_11, amp[327]); 
  FFV1_0(w[95], w[19], w[39], pars->GC_11, amp[328]); 
  FFV2_5_0(w[95], w[40], w[18], pars->GC_51, pars->GC_58, amp[329]); 
  FFV1_0(w[105], w[6], w[110], pars->GC_11, amp[330]); 
  FFV1_0(w[105], w[53], w[39], pars->GC_11, amp[331]); 
  FFV1_0(w[105], w[40], w[48], pars->GC_11, amp[332]); 
  FFV2_5_0(w[111], w[53], w[18], pars->GC_51, pars->GC_58, amp[333]); 
  FFV1_0(w[111], w[19], w[48], pars->GC_11, amp[334]); 
  FFV1_0(w[0], w[19], w[110], pars->GC_11, amp[335]); 
  FFV1_0(w[95], w[23], w[43], pars->GC_11, amp[336]); 
  FFV1_0(w[95], w[44], w[10], pars->GC_2, amp[337]); 
  FFV1_0(w[98], w[78], w[43], pars->GC_11, amp[338]); 
  FFV1_0(w[98], w[5], w[112], pars->GC_11, amp[339]); 
  FFV1_0(w[98], w[44], w[48], pars->GC_11, amp[340]); 
  FFV1_0(w[113], w[78], w[10], pars->GC_2, amp[341]); 
  FFV1_0(w[113], w[23], w[48], pars->GC_11, amp[342]); 
  FFV1_0(w[0], w[23], w[112], pars->GC_11, amp[343]); 
  FFV1_0(w[95], w[27], w[43], pars->GC_11, amp[344]); 
  FFV2_5_0(w[95], w[44], w[18], pars->GC_51, pars->GC_58, amp[345]); 
  FFV1_0(w[105], w[78], w[43], pars->GC_11, amp[346]); 
  FFV1_0(w[105], w[5], w[112], pars->GC_11, amp[347]); 
  FFV1_0(w[105], w[44], w[48], pars->GC_11, amp[348]); 
  FFV2_5_0(w[113], w[78], w[18], pars->GC_51, pars->GC_58, amp[349]); 
  FFV1_0(w[113], w[27], w[48], pars->GC_11, amp[350]); 
  FFV1_0(w[0], w[27], w[112], pars->GC_11, amp[351]); 
  FFV1_0(w[7], w[55], w[115], pars->GC_11, amp[352]); 
  FFV1_0(w[7], w[12], w[116], pars->GC_11, amp[353]); 
  FFV1_0(w[13], w[6], w[116], pars->GC_11, amp[354]); 
  FFV1_0(w[13], w[4], w[115], pars->GC_11, amp[355]); 
  FFV1_0(w[7], w[81], w[99], pars->GC_11, amp[356]); 
  FFV1_0(w[7], w[17], w[117], pars->GC_11, amp[357]); 
  FFV1_0(w[15], w[6], w[117], pars->GC_11, amp[358]); 
  FFV1_0(w[15], w[4], w[99], pars->GC_11, amp[359]); 
  FFV1_0(w[7], w[12], w[118], pars->GC_11, amp[360]); 
  FFV1_0(w[7], w[81], w[102], pars->GC_11, amp[361]); 
  FFV1_0(w[13], w[6], w[118], pars->GC_11, amp[362]); 
  FFV1_0(w[7], w[55], w[119], pars->GC_11, amp[363]); 
  FFV1_0(w[7], w[17], w[120], pars->GC_11, amp[364]); 
  FFV1_0(w[13], w[4], w[119], pars->GC_11, amp[365]); 
  FFV1_0(w[15], w[6], w[120], pars->GC_11, amp[366]); 
  FFV1_0(w[15], w[4], w[102], pars->GC_11, amp[367]); 
  FFV1_0(w[7], w[58], w[115], pars->GC_11, amp[368]); 
  FFV1_0(w[7], w[19], w[116], pars->GC_11, amp[369]); 
  FFV1_0(w[20], w[6], w[116], pars->GC_11, amp[370]); 
  FFV1_0(w[20], w[4], w[115], pars->GC_11, amp[371]); 
  FFV1_0(w[7], w[81], w[106], pars->GC_11, amp[372]); 
  FFV1_0(w[7], w[17], w[121], pars->GC_11, amp[373]); 
  FFV1_0(w[15], w[6], w[121], pars->GC_11, amp[374]); 
  FFV1_0(w[15], w[4], w[106], pars->GC_11, amp[375]); 
  FFV1_0(w[7], w[19], w[118], pars->GC_11, amp[376]); 
  FFV1_0(w[7], w[81], w[108], pars->GC_11, amp[377]); 
  FFV1_0(w[20], w[6], w[118], pars->GC_11, amp[378]); 
  FFV1_0(w[7], w[58], w[119], pars->GC_11, amp[379]); 
  FFV1_0(w[7], w[17], w[122], pars->GC_11, amp[380]); 
  FFV1_0(w[20], w[4], w[119], pars->GC_11, amp[381]); 
  FFV1_0(w[15], w[6], w[122], pars->GC_11, amp[382]); 
  FFV1_0(w[15], w[4], w[108], pars->GC_11, amp[383]); 
  FFV1_0(w[114], w[12], w[68], pars->GC_11, amp[384]); 
  FFV1_0(w[114], w[69], w[10], pars->GC_2, amp[385]); 
  FFV1_0(w[98], w[6], w[123], pars->GC_11, amp[386]); 
  FFV1_0(w[98], w[17], w[68], pars->GC_11, amp[387]); 
  FFV1_0(w[98], w[69], w[9], pars->GC_11, amp[388]); 
  FFV1_0(w[124], w[17], w[10], pars->GC_2, amp[389]); 
  FFV1_0(w[124], w[12], w[9], pars->GC_11, amp[390]); 
  FFV1_0(w[0], w[12], w[123], pars->GC_11, amp[391]); 
  FFV1_0(w[114], w[19], w[68], pars->GC_11, amp[392]); 
  FFV2_5_0(w[114], w[69], w[18], pars->GC_51, pars->GC_58, amp[393]); 
  FFV1_0(w[105], w[6], w[123], pars->GC_11, amp[394]); 
  FFV1_0(w[105], w[17], w[68], pars->GC_11, amp[395]); 
  FFV1_0(w[105], w[69], w[9], pars->GC_11, amp[396]); 
  FFV2_5_0(w[124], w[17], w[18], pars->GC_51, pars->GC_58, amp[397]); 
  FFV1_0(w[124], w[19], w[9], pars->GC_11, amp[398]); 
  FFV1_0(w[0], w[19], w[123], pars->GC_11, amp[399]); 
  FFV1_0(w[114], w[55], w[43], pars->GC_11, amp[400]); 
  FFV1_0(w[114], w[72], w[10], pars->GC_2, amp[401]); 
  FFV1_0(w[98], w[81], w[43], pars->GC_11, amp[402]); 
  FFV1_0(w[98], w[4], w[125], pars->GC_11, amp[403]); 
  FFV1_0(w[98], w[72], w[9], pars->GC_11, amp[404]); 
  FFV1_0(w[113], w[81], w[10], pars->GC_2, amp[405]); 
  FFV1_0(w[113], w[55], w[9], pars->GC_11, amp[406]); 
  FFV1_0(w[0], w[55], w[125], pars->GC_11, amp[407]); 
  FFV1_0(w[114], w[58], w[43], pars->GC_11, amp[408]); 
  FFV2_5_0(w[114], w[72], w[18], pars->GC_51, pars->GC_58, amp[409]); 
  FFV1_0(w[105], w[81], w[43], pars->GC_11, amp[410]); 
  FFV1_0(w[105], w[4], w[125], pars->GC_11, amp[411]); 
  FFV1_0(w[105], w[72], w[9], pars->GC_11, amp[412]); 
  FFV2_5_0(w[113], w[81], w[18], pars->GC_51, pars->GC_58, amp[413]); 
  FFV1_0(w[113], w[58], w[9], pars->GC_11, amp[414]); 
  FFV1_0(w[0], w[58], w[125], pars->GC_11, amp[415]); 
  FFV1_0(w[7], w[55], w[127], pars->GC_11, amp[416]); 
  FFV1_0(w[7], w[23], w[128], pars->GC_11, amp[417]); 
  FFV1_0(w[13], w[5], w[128], pars->GC_11, amp[418]); 
  FFV1_0(w[13], w[4], w[127], pars->GC_11, amp[419]); 
  FFV1_0(w[7], w[57], w[100], pars->GC_11, amp[420]); 
  FFV1_0(w[7], w[26], w[117], pars->GC_11, amp[421]); 
  FFV1_0(w[25], w[5], w[117], pars->GC_11, amp[422]); 
  FFV1_0(w[25], w[4], w[100], pars->GC_11, amp[423]); 
  FFV1_0(w[7], w[23], w[129], pars->GC_11, amp[424]); 
  FFV1_0(w[7], w[57], w[104], pars->GC_11, amp[425]); 
  FFV1_0(w[13], w[5], w[129], pars->GC_11, amp[426]); 
  FFV1_0(w[7], w[55], w[130], pars->GC_11, amp[427]); 
  FFV1_0(w[7], w[26], w[120], pars->GC_11, amp[428]); 
  FFV1_0(w[13], w[4], w[130], pars->GC_11, amp[429]); 
  FFV1_0(w[25], w[5], w[120], pars->GC_11, amp[430]); 
  FFV1_0(w[25], w[4], w[104], pars->GC_11, amp[431]); 
  FFV1_0(w[7], w[58], w[127], pars->GC_11, amp[432]); 
  FFV1_0(w[7], w[27], w[128], pars->GC_11, amp[433]); 
  FFV1_0(w[20], w[5], w[128], pars->GC_11, amp[434]); 
  FFV1_0(w[20], w[4], w[127], pars->GC_11, amp[435]); 
  FFV1_0(w[7], w[57], w[107], pars->GC_11, amp[436]); 
  FFV1_0(w[7], w[26], w[121], pars->GC_11, amp[437]); 
  FFV1_0(w[25], w[5], w[121], pars->GC_11, amp[438]); 
  FFV1_0(w[25], w[4], w[107], pars->GC_11, amp[439]); 
  FFV1_0(w[7], w[27], w[129], pars->GC_11, amp[440]); 
  FFV1_0(w[7], w[57], w[109], pars->GC_11, amp[441]); 
  FFV1_0(w[20], w[5], w[129], pars->GC_11, amp[442]); 
  FFV1_0(w[7], w[58], w[130], pars->GC_11, amp[443]); 
  FFV1_0(w[7], w[26], w[122], pars->GC_11, amp[444]); 
  FFV1_0(w[20], w[4], w[130], pars->GC_11, amp[445]); 
  FFV1_0(w[25], w[5], w[122], pars->GC_11, amp[446]); 
  FFV1_0(w[25], w[4], w[109], pars->GC_11, amp[447]); 
  FFV1_0(w[126], w[23], w[68], pars->GC_11, amp[448]); 
  FFV1_0(w[126], w[91], w[10], pars->GC_2, amp[449]); 
  FFV1_0(w[98], w[5], w[131], pars->GC_11, amp[450]); 
  FFV1_0(w[98], w[26], w[68], pars->GC_11, amp[451]); 
  FFV1_0(w[98], w[91], w[21], pars->GC_11, amp[452]); 
  FFV1_0(w[124], w[26], w[10], pars->GC_2, amp[453]); 
  FFV1_0(w[124], w[23], w[21], pars->GC_11, amp[454]); 
  FFV1_0(w[0], w[23], w[131], pars->GC_11, amp[455]); 
  FFV1_0(w[126], w[27], w[68], pars->GC_11, amp[456]); 
  FFV2_5_0(w[126], w[91], w[18], pars->GC_51, pars->GC_58, amp[457]); 
  FFV1_0(w[105], w[5], w[131], pars->GC_11, amp[458]); 
  FFV1_0(w[105], w[26], w[68], pars->GC_11, amp[459]); 
  FFV1_0(w[105], w[91], w[21], pars->GC_11, amp[460]); 
  FFV2_5_0(w[124], w[26], w[18], pars->GC_51, pars->GC_58, amp[461]); 
  FFV1_0(w[124], w[27], w[21], pars->GC_11, amp[462]); 
  FFV1_0(w[0], w[27], w[131], pars->GC_11, amp[463]); 
  FFV1_0(w[126], w[55], w[39], pars->GC_11, amp[464]); 
  FFV1_0(w[126], w[93], w[10], pars->GC_2, amp[465]); 
  FFV1_0(w[98], w[57], w[39], pars->GC_11, amp[466]); 
  FFV1_0(w[98], w[4], w[132], pars->GC_11, amp[467]); 
  FFV1_0(w[98], w[93], w[21], pars->GC_11, amp[468]); 
  FFV1_0(w[111], w[57], w[10], pars->GC_2, amp[469]); 
  FFV1_0(w[111], w[55], w[21], pars->GC_11, amp[470]); 
  FFV1_0(w[0], w[55], w[132], pars->GC_11, amp[471]); 
  FFV1_0(w[126], w[58], w[39], pars->GC_11, amp[472]); 
  FFV2_5_0(w[126], w[93], w[18], pars->GC_51, pars->GC_58, amp[473]); 
  FFV1_0(w[105], w[57], w[39], pars->GC_11, amp[474]); 
  FFV1_0(w[105], w[4], w[132], pars->GC_11, amp[475]); 
  FFV1_0(w[105], w[93], w[21], pars->GC_11, amp[476]); 
  FFV2_5_0(w[111], w[57], w[18], pars->GC_51, pars->GC_58, amp[477]); 
  FFV1_0(w[111], w[58], w[21], pars->GC_11, amp[478]); 
  FFV1_0(w[0], w[58], w[132], pars->GC_11, amp[479]); 
  FFV1_0(w[71], w[6], w[100], pars->GC_11, amp[480]); 
  FFV1_0(w[71], w[5], w[99], pars->GC_11, amp[481]); 
  FFV1_0(w[1], w[91], w[99], pars->GC_11, amp[482]); 
  FFV1_0(w[1], w[69], w[100], pars->GC_11, amp[483]); 
  FFV1_0(w[31], w[6], w[133], pars->GC_11, amp[484]); 
  FFV1_0(w[31], w[5], w[134], pars->GC_11, amp[485]); 
  FFV1_0(w[1], w[23], w[134], pars->GC_11, amp[486]); 
  FFV1_0(w[1], w[12], w[133], pars->GC_11, amp[487]); 
  FFV1_0(w[31], w[6], w[135], pars->GC_11, amp[488]); 
  FFV1_0(w[31], w[5], w[136], pars->GC_11, amp[489]); 
  FFV1_0(w[71], w[6], w[104], pars->GC_11, amp[490]); 
  FFV1_0(w[71], w[5], w[102], pars->GC_11, amp[491]); 
  FFV1_0(w[1], w[69], w[104], pars->GC_11, amp[492]); 
  FFV1_0(w[1], w[23], w[136], pars->GC_11, amp[493]); 
  FFV1_0(w[1], w[91], w[102], pars->GC_11, amp[494]); 
  FFV1_0(w[1], w[12], w[135], pars->GC_11, amp[495]); 
  FFV1_0(w[71], w[6], w[107], pars->GC_11, amp[496]); 
  FFV1_0(w[71], w[5], w[106], pars->GC_11, amp[497]); 
  FFV1_0(w[1], w[91], w[106], pars->GC_11, amp[498]); 
  FFV1_0(w[1], w[69], w[107], pars->GC_11, amp[499]); 
  FFV1_0(w[38], w[6], w[133], pars->GC_11, amp[500]); 
  FFV1_0(w[38], w[5], w[134], pars->GC_11, amp[501]); 
  FFV1_0(w[1], w[27], w[134], pars->GC_11, amp[502]); 
  FFV1_0(w[1], w[19], w[133], pars->GC_11, amp[503]); 
  FFV1_0(w[38], w[6], w[135], pars->GC_11, amp[504]); 
  FFV1_0(w[38], w[5], w[136], pars->GC_11, amp[505]); 
  FFV1_0(w[71], w[6], w[109], pars->GC_11, amp[506]); 
  FFV1_0(w[71], w[5], w[108], pars->GC_11, amp[507]); 
  FFV1_0(w[1], w[69], w[109], pars->GC_11, amp[508]); 
  FFV1_0(w[1], w[27], w[136], pars->GC_11, amp[509]); 
  FFV1_0(w[1], w[91], w[108], pars->GC_11, amp[510]); 
  FFV1_0(w[1], w[19], w[135], pars->GC_11, amp[511]); 
  FFV1_0(w[42], w[6], w[117], pars->GC_11, amp[512]); 
  FFV1_0(w[42], w[4], w[99], pars->GC_11, amp[513]); 
  FFV1_0(w[1], w[93], w[99], pars->GC_11, amp[514]); 
  FFV1_0(w[1], w[40], w[117], pars->GC_11, amp[515]); 
  FFV1_0(w[31], w[6], w[137], pars->GC_11, amp[516]); 
  FFV1_0(w[31], w[4], w[138], pars->GC_11, amp[517]); 
  FFV1_0(w[1], w[55], w[138], pars->GC_11, amp[518]); 
  FFV1_0(w[1], w[12], w[137], pars->GC_11, amp[519]); 
  FFV1_0(w[31], w[6], w[139], pars->GC_11, amp[520]); 
  FFV1_0(w[31], w[4], w[140], pars->GC_11, amp[521]); 
  FFV1_0(w[42], w[6], w[120], pars->GC_11, amp[522]); 
  FFV1_0(w[42], w[4], w[102], pars->GC_11, amp[523]); 
  FFV1_0(w[1], w[40], w[120], pars->GC_11, amp[524]); 
  FFV1_0(w[1], w[55], w[140], pars->GC_11, amp[525]); 
  FFV1_0(w[1], w[93], w[102], pars->GC_11, amp[526]); 
  FFV1_0(w[1], w[12], w[139], pars->GC_11, amp[527]); 
  FFV1_0(w[42], w[6], w[121], pars->GC_11, amp[528]); 
  FFV1_0(w[42], w[4], w[106], pars->GC_11, amp[529]); 
  FFV1_0(w[1], w[93], w[106], pars->GC_11, amp[530]); 
  FFV1_0(w[1], w[40], w[121], pars->GC_11, amp[531]); 
  FFV1_0(w[38], w[6], w[137], pars->GC_11, amp[532]); 
  FFV1_0(w[38], w[4], w[138], pars->GC_11, amp[533]); 
  FFV1_0(w[1], w[58], w[138], pars->GC_11, amp[534]); 
  FFV1_0(w[1], w[19], w[137], pars->GC_11, amp[535]); 
  FFV1_0(w[38], w[6], w[139], pars->GC_11, amp[536]); 
  FFV1_0(w[38], w[4], w[140], pars->GC_11, amp[537]); 
  FFV1_0(w[42], w[6], w[122], pars->GC_11, amp[538]); 
  FFV1_0(w[42], w[4], w[108], pars->GC_11, amp[539]); 
  FFV1_0(w[1], w[40], w[122], pars->GC_11, amp[540]); 
  FFV1_0(w[1], w[58], w[140], pars->GC_11, amp[541]); 
  FFV1_0(w[1], w[93], w[108], pars->GC_11, amp[542]); 
  FFV1_0(w[1], w[19], w[139], pars->GC_11, amp[543]); 
  FFV1_0(w[45], w[5], w[117], pars->GC_11, amp[544]); 
  FFV1_0(w[45], w[4], w[100], pars->GC_11, amp[545]); 
  FFV1_0(w[1], w[72], w[100], pars->GC_11, amp[546]); 
  FFV1_0(w[1], w[44], w[117], pars->GC_11, amp[547]); 
  FFV1_0(w[31], w[5], w[141], pars->GC_11, amp[548]); 
  FFV1_0(w[31], w[4], w[142], pars->GC_11, amp[549]); 
  FFV1_0(w[1], w[55], w[142], pars->GC_11, amp[550]); 
  FFV1_0(w[1], w[23], w[141], pars->GC_11, amp[551]); 
  FFV1_0(w[31], w[5], w[143], pars->GC_11, amp[552]); 
  FFV1_0(w[31], w[4], w[144], pars->GC_11, amp[553]); 
  FFV1_0(w[45], w[5], w[120], pars->GC_11, amp[554]); 
  FFV1_0(w[45], w[4], w[104], pars->GC_11, amp[555]); 
  FFV1_0(w[1], w[44], w[120], pars->GC_11, amp[556]); 
  FFV1_0(w[1], w[55], w[144], pars->GC_11, amp[557]); 
  FFV1_0(w[1], w[72], w[104], pars->GC_11, amp[558]); 
  FFV1_0(w[1], w[23], w[143], pars->GC_11, amp[559]); 
  FFV1_0(w[45], w[5], w[121], pars->GC_11, amp[560]); 
  FFV1_0(w[45], w[4], w[107], pars->GC_11, amp[561]); 
  FFV1_0(w[1], w[72], w[107], pars->GC_11, amp[562]); 
  FFV1_0(w[1], w[44], w[121], pars->GC_11, amp[563]); 
  FFV1_0(w[38], w[5], w[141], pars->GC_11, amp[564]); 
  FFV1_0(w[38], w[4], w[142], pars->GC_11, amp[565]); 
  FFV1_0(w[1], w[58], w[142], pars->GC_11, amp[566]); 
  FFV1_0(w[1], w[27], w[141], pars->GC_11, amp[567]); 
  FFV1_0(w[38], w[5], w[143], pars->GC_11, amp[568]); 
  FFV1_0(w[38], w[4], w[144], pars->GC_11, amp[569]); 
  FFV1_0(w[45], w[5], w[122], pars->GC_11, amp[570]); 
  FFV1_0(w[45], w[4], w[109], pars->GC_11, amp[571]); 
  FFV1_0(w[1], w[44], w[122], pars->GC_11, amp[572]); 
  FFV1_0(w[1], w[58], w[144], pars->GC_11, amp[573]); 
  FFV1_0(w[1], w[72], w[109], pars->GC_11, amp[574]); 
  FFV1_0(w[1], w[27], w[143], pars->GC_11, amp[575]); 

}
double CPPProcess::matrix_uu_taptamuuuux_no_h() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 576; 
  const int ncolor = 6; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1, 1, 1, 1, 1, 1}; 
  static const double cf[ncolor][ncolor] = {{27, 9, 9, 3, 3, 9}, {9, 27, 3, 9,
      9, 3}, {9, 3, 27, 9, 9, 3}, {3, 9, 9, 27, 3, 9}, {3, 9, 9, 3, 27, 9}, {9,
      3, 3, 9, 9, 27}};

  // Calculate color flows
  jamp[0] = +1./4. * (+1./9. * amp[2] + 1./9. * amp[3] + 1./9. * amp[4] + 1./9.
      * amp[5] + 1./9. * amp[6] + 1./9. * amp[7] + 1./9. * amp[10] + 1./9. *
      amp[11] + 1./9. * amp[12] + 1./9. * amp[13] + 1./9. * amp[14] + 1./9. *
      amp[15] + 1./3. * amp[18] + 1./3. * amp[19] + 1./3. * amp[20] + 1./3. *
      amp[21] + 1./3. * amp[22] + 1./3. * amp[23] + 1./3. * amp[26] + 1./3. *
      amp[27] + 1./3. * amp[28] + 1./3. * amp[29] + 1./3. * amp[30] + 1./3. *
      amp[31] + 1./3. * amp[32] + 1./9. * amp[33] + 1./9. * amp[34] + 1./3. *
      amp[35] + 1./3. * amp[36] + 1./9. * amp[37] + 1./3. * amp[38] + 1./9. *
      amp[39] + 1./9. * amp[40] + 1./3. * amp[41] + 1./9. * amp[42] + 1./3. *
      amp[43] + 1./3. * amp[44] + 1./9. * amp[45] + 1./9. * amp[46] + 1./3. *
      amp[47] + 1./3. * amp[48] + 1./9. * amp[49] + 1./9. * amp[50] + 1./3. *
      amp[51] + 1./3. * amp[52] + 1./9. * amp[53] + 1./3. * amp[54] + 1./9. *
      amp[55] + 1./9. * amp[56] + 1./3. * amp[57] + 1./9. * amp[58] + 1./3. *
      amp[59] + 1./3. * amp[60] + 1./9. * amp[61] + 1./9. * amp[62] + 1./3. *
      amp[63] + 1./3. * amp[64] + 1./3. * amp[65] + 1./3. * amp[68] + 1./3. *
      amp[69] + 1./3. * amp[70] + 1./3. * amp[71] + 1./3. * amp[72] + 1./3. *
      amp[73] + 1./3. * amp[76] + 1./3. * amp[77] + 1./3. * amp[78] + 1./3. *
      amp[79] + 1./9. * amp[80] + 1./9. * amp[81] + 1./9. * amp[82] + 1./9. *
      amp[83] + 1./9. * amp[86] + 1./9. * amp[87] + 1./9. * amp[88] + 1./9. *
      amp[89] + 1./9. * amp[90] + 1./9. * amp[91] + 1./9. * amp[94] + 1./9. *
      amp[95] - std::complex<double> (0, 1) * amp[112] - std::complex<double>
      (0, 1) * amp[113] + amp[114] + amp[115] + amp[119] - std::complex<double>
      (0, 1) * amp[120] - std::complex<double> (0, 1) * amp[121] + amp[122] +
      amp[123] + amp[127] + amp[128] + 1./3. * amp[129] + 1./3. * amp[130] +
      amp[131] + amp[132] + 1./3. * amp[133] + amp[134] + 1./3. * amp[135] +
      amp[144] + 1./3. * amp[145] + 1./3. * amp[146] + amp[147] + amp[148] +
      1./3. * amp[149] + amp[150] + 1./3. * amp[151] + amp[160] + amp[161] +
      std::complex<double> (0, 1) * amp[162] + std::complex<double> (0, 1) *
      amp[163] + amp[166] + amp[168] + amp[169] + std::complex<double> (0, 1) *
      amp[170] + std::complex<double> (0, 1) * amp[171] + amp[174] + 1./3. *
      amp[176] + 1./3. * amp[177] + 1./3. * amp[178] + 1./3. * amp[179] + 1./3.
      * amp[182] + 1./3. * amp[183] + 1./3. * amp[184] + 1./3. * amp[185] +
      1./3. * amp[186] + 1./3. * amp[187] + 1./3. * amp[190] + 1./3. * amp[191]
      + std::complex<double> (0, 1) * amp[192] + std::complex<double> (0, 1) *
      amp[193] + amp[196] + amp[197] + amp[198] + std::complex<double> (0, 1) *
      amp[200] + std::complex<double> (0, 1) * amp[201] + amp[204] + amp[205] +
      amp[206] + 1./3. * amp[210] + 1./3. * amp[211] + 1./3. * amp[212] + 1./3.
      * amp[213] + 1./3. * amp[214] + 1./3. * amp[215] + 1./3. * amp[218] +
      1./3. * amp[219] + 1./3. * amp[220] + 1./3. * amp[221] + 1./3. * amp[222]
      + 1./3. * amp[223] + 1./3. * amp[228] + amp[229] + 1./3. * amp[230] +
      amp[231] + 1./3. * amp[236] + amp[237] + amp[238] + 1./3. * amp[239] +
      1./3. * amp[244] + amp[245] + 1./3. * amp[246] + amp[247] + 1./3. *
      amp[252] + amp[253] + amp[254] + 1./3. * amp[255] + amp[274] + amp[275] -
      std::complex<double> (0, 1) * amp[276] - std::complex<double> (0, 1) *
      amp[277] + amp[279] + amp[282] + amp[283] - std::complex<double> (0, 1) *
      amp[284] - std::complex<double> (0, 1) * amp[285] + amp[287] + amp[288] +
      1./3. * amp[289] + 1./3. * amp[290] + amp[291] + amp[292] + 1./3. *
      amp[296] + amp[297] + 1./3. * amp[298] + amp[304] + 1./3. * amp[305] +
      1./3. * amp[306] + amp[307] + amp[308] + 1./3. * amp[312] + amp[313] +
      1./3. * amp[314] + amp[320] + amp[321] + std::complex<double> (0, 1) *
      amp[322] + amp[324] + std::complex<double> (0, 1) * amp[327] + amp[328] +
      amp[329] + std::complex<double> (0, 1) * amp[330] + amp[332] +
      std::complex<double> (0, 1) * amp[335] + 1./3. * amp[336] + 1./3. *
      amp[337] + 1./3. * amp[338] + 1./3. * amp[340] + 1./3. * amp[341] + 1./3.
      * amp[342] + 1./3. * amp[344] + 1./3. * amp[345] + 1./3. * amp[346] +
      1./3. * amp[348] + 1./3. * amp[349] + 1./3. * amp[350] + 1./3. * amp[352]
      + 1./9. * amp[353] + 1./9. * amp[354] + 1./3. * amp[355] + 1./3. *
      amp[356] + 1./9. * amp[357] + 1./9. * amp[358] + 1./3. * amp[359] + 1./9.
      * amp[360] + 1./3. * amp[361] + 1./9. * amp[362] + 1./3. * amp[363] +
      1./9. * amp[364] + 1./3. * amp[365] + 1./9. * amp[366] + 1./3. * amp[367]
      + 1./3. * amp[368] + 1./9. * amp[369] + 1./9. * amp[370] + 1./3. *
      amp[371] + 1./3. * amp[372] + 1./9. * amp[373] + 1./9. * amp[374] + 1./3.
      * amp[375] + 1./9. * amp[376] + 1./3. * amp[377] + 1./9. * amp[378] +
      1./3. * amp[379] + 1./9. * amp[380] + 1./3. * amp[381] + 1./9. * amp[382]
      + 1./3. * amp[383] + 1./3. * amp[384] + 1./3. * amp[385] + 1./3. *
      amp[387] + 1./3. * amp[388] + 1./3. * amp[389] + 1./3. * amp[390] + 1./3.
      * amp[392] + 1./3. * amp[393] + 1./3. * amp[395] + 1./3. * amp[396] +
      1./3. * amp[397] + 1./3. * amp[398] + 1./9. * amp[400] + 1./9. * amp[401]
      + 1./9. * amp[402] + 1./9. * amp[404] + 1./9. * amp[405] + 1./9. *
      amp[406] + 1./9. * amp[408] + 1./9. * amp[409] + 1./9. * amp[410] + 1./9.
      * amp[412] + 1./9. * amp[413] + 1./9. * amp[414] + 1./3. * amp[421] +
      1./3. * amp[422] + amp[423] + amp[427] + 1./3. * amp[428] + amp[429] +
      1./3. * amp[430] + amp[431] + 1./3. * amp[437] + 1./3. * amp[438] +
      amp[439] + amp[443] + 1./3. * amp[444] + amp[445] + 1./3. * amp[446] +
      amp[447] - std::complex<double> (0, 1) * amp[450] + amp[451] + amp[453] +
      amp[454] - std::complex<double> (0, 1) * amp[455] - std::complex<double>
      (0, 1) * amp[458] + amp[459] + amp[461] + amp[462] - std::complex<double>
      (0, 1) * amp[463] + amp[483] + amp[484] + 1./3. * amp[485] + 1./3. *
      amp[486] + amp[487] + 1./3. * amp[489] + amp[492] + 1./3. * amp[493] +
      amp[499] + amp[500] + 1./3. * amp[501] + 1./3. * amp[502] + amp[503] +
      1./3. * amp[505] + amp[508] + 1./3. * amp[509] + 1./3. * amp[512] +
      amp[513] + 1./3. * amp[515] + amp[521] + 1./3. * amp[522] + amp[523] +
      1./3. * amp[524] + amp[525] + 1./3. * amp[528] + amp[529] + 1./3. *
      amp[531] + amp[537] + 1./3. * amp[538] + amp[539] + 1./3. * amp[540] +
      amp[541] + 1./9. * amp[544] + 1./3. * amp[545] + 1./3. * amp[546] + 1./9.
      * amp[547] + 1./9. * amp[548] + 1./3. * amp[549] + 1./3. * amp[550] +
      1./9. * amp[551] + 1./9. * amp[552] + 1./3. * amp[553] + 1./9. * amp[554]
      + 1./3. * amp[555] + 1./9. * amp[556] + 1./3. * amp[557] + 1./3. *
      amp[558] + 1./9. * amp[559] + 1./9. * amp[560] + 1./3. * amp[561] + 1./3.
      * amp[562] + 1./9. * amp[563] + 1./9. * amp[564] + 1./3. * amp[565] +
      1./3. * amp[566] + 1./9. * amp[567] + 1./9. * amp[568] + 1./3. * amp[569]
      + 1./9. * amp[570] + 1./3. * amp[571] + 1./9. * amp[572] + 1./3. *
      amp[573] + 1./3. * amp[574] + 1./9. * amp[575]);
  jamp[1] = +1./4. * (-1./3. * amp[2] - 1./3. * amp[3] - 1./3. * amp[4] - 1./3.
      * amp[5] - 1./3. * amp[6] - 1./3. * amp[7] - 1./3. * amp[10] - 1./3. *
      amp[11] - 1./3. * amp[12] - 1./3. * amp[13] - 1./3. * amp[14] - 1./3. *
      amp[15] - 1./9. * amp[18] - 1./9. * amp[19] - 1./9. * amp[20] - 1./9. *
      amp[21] - 1./9. * amp[22] - 1./9. * amp[23] - 1./9. * amp[26] - 1./9. *
      amp[27] - 1./9. * amp[28] - 1./9. * amp[29] - 1./9. * amp[30] - 1./9. *
      amp[31] - 1./9. * amp[32] - 1./3. * amp[33] - 1./3. * amp[34] - 1./9. *
      amp[35] - 1./9. * amp[36] - 1./3. * amp[37] - 1./9. * amp[38] - 1./3. *
      amp[39] - 1./3. * amp[40] - 1./9. * amp[41] - 1./3. * amp[42] - 1./9. *
      amp[43] - 1./9. * amp[44] - 1./3. * amp[45] - 1./3. * amp[46] - 1./9. *
      amp[47] - 1./9. * amp[48] - 1./3. * amp[49] - 1./3. * amp[50] - 1./9. *
      amp[51] - 1./9. * amp[52] - 1./3. * amp[53] - 1./9. * amp[54] - 1./3. *
      amp[55] - 1./3. * amp[56] - 1./9. * amp[57] - 1./3. * amp[58] - 1./9. *
      amp[59] - 1./9. * amp[60] - 1./3. * amp[61] - 1./3. * amp[62] - 1./9. *
      amp[63] - 1./9. * amp[64] - 1./9. * amp[65] - 1./9. * amp[68] - 1./9. *
      amp[69] - 1./9. * amp[70] - 1./9. * amp[71] - 1./9. * amp[72] - 1./9. *
      amp[73] - 1./9. * amp[76] - 1./9. * amp[77] - 1./9. * amp[78] - 1./9. *
      amp[79] - 1./3. * amp[80] - 1./3. * amp[81] - 1./3. * amp[82] - 1./3. *
      amp[83] - 1./3. * amp[86] - 1./3. * amp[87] - 1./3. * amp[88] - 1./3. *
      amp[89] - 1./3. * amp[90] - 1./3. * amp[91] - 1./3. * amp[94] - 1./3. *
      amp[95] - std::complex<double> (0, 1) * amp[96] - std::complex<double>
      (0, 1) * amp[97] - amp[100] - amp[101] - amp[102] - std::complex<double>
      (0, 1) * amp[104] - std::complex<double> (0, 1) * amp[105] - amp[108] -
      amp[109] - amp[110] - 1./3. * amp[114] - 1./3. * amp[115] - 1./3. *
      amp[116] - 1./3. * amp[117] - 1./3. * amp[118] - 1./3. * amp[119] - 1./3.
      * amp[122] - 1./3. * amp[123] - 1./3. * amp[124] - 1./3. * amp[125] -
      1./3. * amp[126] - 1./3. * amp[127] - 1./3. * amp[132] - amp[133] - 1./3.
      * amp[134] - amp[135] - 1./3. * amp[140] - amp[141] - amp[142] - 1./3. *
      amp[143] - 1./3. * amp[148] - amp[149] - 1./3. * amp[150] - amp[151] -
      1./3. * amp[156] - amp[157] - amp[158] - 1./3. * amp[159] - amp[178] -
      amp[179] + std::complex<double> (0, 1) * amp[180] + std::complex<double>
      (0, 1) * amp[181] - amp[183] - amp[186] - amp[187] + std::complex<double>
      (0, 1) * amp[188] + std::complex<double> (0, 1) * amp[189] - amp[191] +
      std::complex<double> (0, 1) * amp[208] + std::complex<double> (0, 1) *
      amp[209] - amp[210] - amp[211] - amp[215] + std::complex<double> (0, 1) *
      amp[216] + std::complex<double> (0, 1) * amp[217] - amp[218] - amp[219] -
      amp[223] - amp[224] - 1./3. * amp[225] - 1./3. * amp[226] - amp[227] -
      amp[228] - 1./3. * amp[229] - amp[230] - 1./3. * amp[231] - amp[240] -
      1./3. * amp[241] - 1./3. * amp[242] - amp[243] - amp[244] - 1./3. *
      amp[245] - amp[246] - 1./3. * amp[247] - amp[256] - amp[257] -
      std::complex<double> (0, 1) * amp[258] - std::complex<double> (0, 1) *
      amp[259] - amp[262] - amp[264] - amp[265] - std::complex<double> (0, 1) *
      amp[266] - std::complex<double> (0, 1) * amp[267] - amp[270] - 1./3. *
      amp[272] - 1./3. * amp[273] - 1./3. * amp[274] - 1./3. * amp[275] - 1./3.
      * amp[278] - 1./3. * amp[279] - 1./3. * amp[280] - 1./3. * amp[281] -
      1./3. * amp[282] - 1./3. * amp[283] - 1./3. * amp[286] - 1./3. * amp[287]
      - 1./3. * amp[288] - amp[289] - amp[290] - 1./3. * amp[291] - amp[293] -
      1./3. * amp[299] - amp[300] - 1./3. * amp[301] - 1./3. * amp[304] -
      amp[305] - amp[306] - 1./3. * amp[307] - amp[309] - 1./3. * amp[315] -
      amp[316] - 1./3. * amp[317] - 1./3. * amp[320] - 1./3. * amp[321] - 1./3.
      * amp[323] - 1./3. * amp[324] - 1./3. * amp[325] - 1./3. * amp[326] -
      1./3. * amp[328] - 1./3. * amp[329] - 1./3. * amp[331] - 1./3. * amp[332]
      - 1./3. * amp[333] - 1./3. * amp[334] - amp[336] - amp[337] -
      std::complex<double> (0, 1) * amp[339] - amp[340] - std::complex<double>
      (0, 1) * amp[343] - amp[344] - amp[345] - std::complex<double> (0, 1) *
      amp[347] - amp[348] - std::complex<double> (0, 1) * amp[351] - 1./3. *
      amp[357] - 1./3. * amp[358] - amp[359] - amp[363] - 1./3. * amp[364] -
      amp[365] - 1./3. * amp[366] - amp[367] - 1./3. * amp[373] - 1./3. *
      amp[374] - amp[375] - amp[379] - 1./3. * amp[380] - amp[381] - 1./3. *
      amp[382] - amp[383] + std::complex<double> (0, 1) * amp[386] - amp[387] -
      amp[389] - amp[390] + std::complex<double> (0, 1) * amp[391] +
      std::complex<double> (0, 1) * amp[394] - amp[395] - amp[397] - amp[398] +
      std::complex<double> (0, 1) * amp[399] - 1./3. * amp[416] - 1./9. *
      amp[417] - 1./9. * amp[418] - 1./3. * amp[419] - 1./3. * amp[420] - 1./9.
      * amp[421] - 1./9. * amp[422] - 1./3. * amp[423] - 1./9. * amp[424] -
      1./3. * amp[425] - 1./9. * amp[426] - 1./3. * amp[427] - 1./9. * amp[428]
      - 1./3. * amp[429] - 1./9. * amp[430] - 1./3. * amp[431] - 1./3. *
      amp[432] - 1./9. * amp[433] - 1./9. * amp[434] - 1./3. * amp[435] - 1./3.
      * amp[436] - 1./9. * amp[437] - 1./9. * amp[438] - 1./3. * amp[439] -
      1./9. * amp[440] - 1./3. * amp[441] - 1./9. * amp[442] - 1./3. * amp[443]
      - 1./9. * amp[444] - 1./3. * amp[445] - 1./9. * amp[446] - 1./3. *
      amp[447] - 1./3. * amp[448] - 1./3. * amp[449] - 1./3. * amp[451] - 1./3.
      * amp[452] - 1./3. * amp[453] - 1./3. * amp[454] - 1./3. * amp[456] -
      1./3. * amp[457] - 1./3. * amp[459] - 1./3. * amp[460] - 1./3. * amp[461]
      - 1./3. * amp[462] - 1./9. * amp[464] - 1./9. * amp[465] - 1./9. *
      amp[466] - 1./9. * amp[468] - 1./9. * amp[469] - 1./9. * amp[470] - 1./9.
      * amp[472] - 1./9. * amp[473] - 1./9. * amp[474] - 1./9. * amp[476] -
      1./9. * amp[477] - 1./9. * amp[478] - amp[482] - 1./3. * amp[484] -
      amp[485] - amp[486] - 1./3. * amp[487] - 1./3. * amp[488] - amp[494] -
      1./3. * amp[495] - amp[498] - 1./3. * amp[500] - amp[501] - amp[502] -
      1./3. * amp[503] - 1./3. * amp[504] - amp[510] - 1./3. * amp[511] - 1./9.
      * amp[512] - 1./3. * amp[513] - 1./3. * amp[514] - 1./9. * amp[515] -
      1./9. * amp[516] - 1./3. * amp[517] - 1./3. * amp[518] - 1./9. * amp[519]
      - 1./9. * amp[520] - 1./3. * amp[521] - 1./9. * amp[522] - 1./3. *
      amp[523] - 1./9. * amp[524] - 1./3. * amp[525] - 1./3. * amp[526] - 1./9.
      * amp[527] - 1./9. * amp[528] - 1./3. * amp[529] - 1./3. * amp[530] -
      1./9. * amp[531] - 1./9. * amp[532] - 1./3. * amp[533] - 1./3. * amp[534]
      - 1./9. * amp[535] - 1./9. * amp[536] - 1./3. * amp[537] - 1./9. *
      amp[538] - 1./3. * amp[539] - 1./9. * amp[540] - 1./3. * amp[541] - 1./3.
      * amp[542] - 1./9. * amp[543] - 1./3. * amp[544] - amp[545] - 1./3. *
      amp[547] - amp[553] - 1./3. * amp[554] - amp[555] - 1./3. * amp[556] -
      amp[557] - 1./3. * amp[560] - amp[561] - 1./3. * amp[563] - amp[569] -
      1./3. * amp[570] - amp[571] - 1./3. * amp[572] - amp[573]);
  jamp[2] = +1./4. * (+std::complex<double> (0, 1) * amp[16] +
      std::complex<double> (0, 1) * amp[17] - amp[18] - amp[19] - amp[23] +
      std::complex<double> (0, 1) * amp[24] + std::complex<double> (0, 1) *
      amp[25] - amp[26] - amp[27] - amp[31] - amp[32] - 1./3. * amp[33] - 1./3.
      * amp[34] - amp[35] - amp[36] - 1./3. * amp[37] - amp[38] - 1./3. *
      amp[39] - amp[48] - 1./3. * amp[49] - 1./3. * amp[50] - amp[51] - amp[52]
      - 1./3. * amp[53] - amp[54] - 1./3. * amp[55] - amp[64] - amp[65] -
      std::complex<double> (0, 1) * amp[66] - std::complex<double> (0, 1) *
      amp[67] - amp[70] - amp[72] - amp[73] - std::complex<double> (0, 1) *
      amp[74] - std::complex<double> (0, 1) * amp[75] - amp[78] - 1./3. *
      amp[80] - 1./3. * amp[81] - 1./3. * amp[82] - 1./3. * amp[83] - 1./3. *
      amp[86] - 1./3. * amp[87] - 1./3. * amp[88] - 1./3. * amp[89] - 1./3. *
      amp[90] - 1./3. * amp[91] - 1./3. * amp[94] - 1./3. * amp[95] - 1./9. *
      amp[98] - 1./9. * amp[99] - 1./9. * amp[100] - 1./9. * amp[101] - 1./9. *
      amp[102] - 1./9. * amp[103] - 1./9. * amp[106] - 1./9. * amp[107] - 1./9.
      * amp[108] - 1./9. * amp[109] - 1./9. * amp[110] - 1./9. * amp[111] -
      1./3. * amp[114] - 1./3. * amp[115] - 1./3. * amp[116] - 1./3. * amp[117]
      - 1./3. * amp[118] - 1./3. * amp[119] - 1./3. * amp[122] - 1./3. *
      amp[123] - 1./3. * amp[124] - 1./3. * amp[125] - 1./3. * amp[126] - 1./3.
      * amp[127] - 1./3. * amp[128] - 1./9. * amp[129] - 1./9. * amp[130] -
      1./3. * amp[131] - 1./3. * amp[132] - 1./9. * amp[133] - 1./3. * amp[134]
      - 1./9. * amp[135] - 1./9. * amp[136] - 1./3. * amp[137] - 1./9. *
      amp[138] - 1./3. * amp[139] - 1./3. * amp[140] - 1./9. * amp[141] - 1./9.
      * amp[142] - 1./3. * amp[143] - 1./3. * amp[144] - 1./9. * amp[145] -
      1./9. * amp[146] - 1./3. * amp[147] - 1./3. * amp[148] - 1./9. * amp[149]
      - 1./3. * amp[150] - 1./9. * amp[151] - 1./9. * amp[152] - 1./3. *
      amp[153] - 1./9. * amp[154] - 1./3. * amp[155] - 1./3. * amp[156] - 1./9.
      * amp[157] - 1./9. * amp[158] - 1./3. * amp[159] - 1./3. * amp[160] -
      1./3. * amp[161] - 1./3. * amp[164] - 1./3. * amp[165] - 1./3. * amp[166]
      - 1./3. * amp[167] - 1./3. * amp[168] - 1./3. * amp[169] - 1./3. *
      amp[172] - 1./3. * amp[173] - 1./3. * amp[174] - 1./3. * amp[175] - 1./9.
      * amp[176] - 1./9. * amp[177] - 1./9. * amp[178] - 1./9. * amp[179] -
      1./9. * amp[182] - 1./9. * amp[183] - 1./9. * amp[184] - 1./9. * amp[185]
      - 1./9. * amp[186] - 1./9. * amp[187] - 1./9. * amp[190] - 1./9. *
      amp[191] - 1./3. * amp[194] - 1./3. * amp[195] - 1./3. * amp[196] - 1./3.
      * amp[197] - 1./3. * amp[198] - 1./3. * amp[199] - 1./3. * amp[202] -
      1./3. * amp[203] - 1./3. * amp[204] - 1./3. * amp[205] - 1./3. * amp[206]
      - 1./3. * amp[207] - std::complex<double> (0, 1) * amp[208] -
      std::complex<double> (0, 1) * amp[209] - amp[212] - amp[213] - amp[214] -
      std::complex<double> (0, 1) * amp[216] - std::complex<double> (0, 1) *
      amp[217] - amp[220] - amp[221] - amp[222] - 1./3. * amp[232] - amp[233] -
      1./3. * amp[234] - amp[235] - amp[236] - 1./3. * amp[237] - 1./3. *
      amp[238] - amp[239] - 1./3. * amp[248] - amp[249] - 1./3. * amp[250] -
      amp[251] - amp[252] - 1./3. * amp[253] - 1./3. * amp[254] - amp[255] +
      std::complex<double> (0, 1) * amp[258] + std::complex<double> (0, 1) *
      amp[259] - amp[260] - amp[261] - amp[263] + std::complex<double> (0, 1) *
      amp[266] + std::complex<double> (0, 1) * amp[267] - amp[268] - amp[269] -
      amp[271] - 1./3. * amp[288] - 1./9. * amp[289] - 1./9. * amp[290] - 1./3.
      * amp[291] - 1./3. * amp[292] - 1./9. * amp[293] - 1./9. * amp[294] -
      1./3. * amp[295] - 1./9. * amp[296] - 1./3. * amp[297] - 1./9. * amp[298]
      - 1./3. * amp[299] - 1./9. * amp[300] - 1./3. * amp[301] - 1./9. *
      amp[302] - 1./3. * amp[303] - 1./3. * amp[304] - 1./9. * amp[305] - 1./9.
      * amp[306] - 1./3. * amp[307] - 1./3. * amp[308] - 1./9. * amp[309] -
      1./9. * amp[310] - 1./3. * amp[311] - 1./9. * amp[312] - 1./3. * amp[313]
      - 1./9. * amp[314] - 1./3. * amp[315] - 1./9. * amp[316] - 1./3. *
      amp[317] - 1./9. * amp[318] - 1./3. * amp[319] - 1./3. * amp[320] - 1./3.
      * amp[321] - 1./3. * amp[323] - 1./3. * amp[324] - 1./3. * amp[325] -
      1./3. * amp[326] - 1./3. * amp[328] - 1./3. * amp[329] - 1./3. * amp[331]
      - 1./3. * amp[332] - 1./3. * amp[333] - 1./3. * amp[334] - 1./9. *
      amp[336] - 1./9. * amp[337] - 1./9. * amp[338] - 1./9. * amp[340] - 1./9.
      * amp[341] - 1./9. * amp[342] - 1./9. * amp[344] - 1./9. * amp[345] -
      1./9. * amp[346] - 1./9. * amp[348] - 1./9. * amp[349] - 1./9. * amp[350]
      - amp[352] - 1./3. * amp[353] - 1./3. * amp[354] - amp[355] - amp[356] -
      1./3. * amp[360] - amp[361] - 1./3. * amp[362] - amp[368] - 1./3. *
      amp[369] - 1./3. * amp[370] - amp[371] - amp[372] - 1./3. * amp[376] -
      amp[377] - 1./3. * amp[378] - amp[384] - amp[385] - std::complex<double>
      (0, 1) * amp[386] - amp[388] - std::complex<double> (0, 1) * amp[391] -
      amp[392] - amp[393] - std::complex<double> (0, 1) * amp[394] - amp[396] -
      std::complex<double> (0, 1) * amp[399] - 1./3. * amp[400] - 1./3. *
      amp[401] - 1./3. * amp[402] - 1./3. * amp[404] - 1./3. * amp[405] - 1./3.
      * amp[406] - 1./3. * amp[408] - 1./3. * amp[409] - 1./3. * amp[410] -
      1./3. * amp[412] - 1./3. * amp[413] - 1./3. * amp[414] - 1./3. * amp[420]
      - amp[422] - 1./3. * amp[423] - amp[424] - 1./3. * amp[425] - amp[426] -
      amp[430] - 1./3. * amp[431] - 1./3. * amp[436] - amp[438] - 1./3. *
      amp[439] - amp[440] - 1./3. * amp[441] - amp[442] - amp[446] - 1./3. *
      amp[447] - amp[466] + std::complex<double> (0, 1) * amp[467] - amp[469] -
      amp[470] + std::complex<double> (0, 1) * amp[471] - amp[474] +
      std::complex<double> (0, 1) * amp[475] - amp[477] - amp[478] +
      std::complex<double> (0, 1) * amp[479] - 1./3. * amp[480] - amp[481] -
      1./3. * amp[483] - amp[489] - 1./3. * amp[490] - amp[491] - 1./3. *
      amp[492] - amp[493] - 1./3. * amp[496] - amp[497] - 1./3. * amp[499] -
      amp[505] - 1./3. * amp[506] - amp[507] - 1./3. * amp[508] - amp[509] -
      amp[515] - amp[516] - 1./3. * amp[517] - 1./3. * amp[518] - amp[519] -
      1./3. * amp[521] - amp[524] - 1./3. * amp[525] - amp[531] - amp[532] -
      1./3. * amp[533] - 1./3. * amp[534] - amp[535] - 1./3. * amp[537] -
      amp[540] - 1./3. * amp[541] - 1./3. * amp[544] - 1./9. * amp[545] - 1./9.
      * amp[546] - 1./3. * amp[547] - 1./3. * amp[548] - 1./9. * amp[549] -
      1./9. * amp[550] - 1./3. * amp[551] - 1./3. * amp[552] - 1./9. * amp[553]
      - 1./3. * amp[554] - 1./9. * amp[555] - 1./3. * amp[556] - 1./9. *
      amp[557] - 1./9. * amp[558] - 1./3. * amp[559] - 1./3. * amp[560] - 1./9.
      * amp[561] - 1./9. * amp[562] - 1./3. * amp[563] - 1./3. * amp[564] -
      1./9. * amp[565] - 1./9. * amp[566] - 1./3. * amp[567] - 1./3. * amp[568]
      - 1./9. * amp[569] - 1./3. * amp[570] - 1./9. * amp[571] - 1./3. *
      amp[572] - 1./9. * amp[573] - 1./9. * amp[574] - 1./3. * amp[575]);
  jamp[3] = +1./4. * (-std::complex<double> (0, 1) * amp[0] -
      std::complex<double> (0, 1) * amp[1] + amp[2] + amp[3] + amp[7] -
      std::complex<double> (0, 1) * amp[8] - std::complex<double> (0, 1) *
      amp[9] + amp[10] + amp[11] + amp[15] + 1./3. * amp[32] + amp[33] +
      amp[34] + 1./3. * amp[35] + amp[40] + 1./3. * amp[41] + amp[42] + 1./3. *
      amp[43] + 1./3. * amp[48] + amp[49] + amp[50] + 1./3. * amp[51] + amp[56]
      + 1./3. * amp[57] + amp[58] + 1./3. * amp[59] + 1./3. * amp[64] + 1./3. *
      amp[65] + 1./3. * amp[68] + 1./3. * amp[69] + 1./3. * amp[70] + 1./3. *
      amp[71] + 1./3. * amp[72] + 1./3. * amp[73] + 1./3. * amp[76] + 1./3. *
      amp[77] + 1./3. * amp[78] + 1./3. * amp[79] + amp[80] + amp[81] +
      std::complex<double> (0, 1) * amp[84] + std::complex<double> (0, 1) *
      amp[85] + amp[86] + amp[88] + amp[89] + std::complex<double> (0, 1) *
      amp[92] + std::complex<double> (0, 1) * amp[93] + amp[94] + 1./3. *
      amp[98] + 1./3. * amp[99] + 1./3. * amp[100] + 1./3. * amp[101] + 1./3. *
      amp[102] + 1./3. * amp[103] + 1./3. * amp[106] + 1./3. * amp[107] + 1./3.
      * amp[108] + 1./3. * amp[109] + 1./3. * amp[110] + 1./3. * amp[111] +
      std::complex<double> (0, 1) * amp[112] + std::complex<double> (0, 1) *
      amp[113] + amp[116] + amp[117] + amp[118] + std::complex<double> (0, 1) *
      amp[120] + std::complex<double> (0, 1) * amp[121] + amp[124] + amp[125] +
      amp[126] + 1./3. * amp[136] + amp[137] + 1./3. * amp[138] + amp[139] +
      amp[140] + 1./3. * amp[141] + 1./3. * amp[142] + amp[143] + 1./3. *
      amp[152] + amp[153] + 1./3. * amp[154] + amp[155] + amp[156] + 1./3. *
      amp[157] + 1./3. * amp[158] + amp[159] - std::complex<double> (0, 1) *
      amp[162] - std::complex<double> (0, 1) * amp[163] + amp[164] + amp[165] +
      amp[167] - std::complex<double> (0, 1) * amp[170] - std::complex<double>
      (0, 1) * amp[171] + amp[172] + amp[173] + amp[175] + 1./9. * amp[194] +
      1./9. * amp[195] + 1./9. * amp[196] + 1./9. * amp[197] + 1./9. * amp[198]
      + 1./9. * amp[199] + 1./9. * amp[202] + 1./9. * amp[203] + 1./9. *
      amp[204] + 1./9. * amp[205] + 1./9. * amp[206] + 1./9. * amp[207] + 1./3.
      * amp[210] + 1./3. * amp[211] + 1./3. * amp[212] + 1./3. * amp[213] +
      1./3. * amp[214] + 1./3. * amp[215] + 1./3. * amp[218] + 1./3. * amp[219]
      + 1./3. * amp[220] + 1./3. * amp[221] + 1./3. * amp[222] + 1./3. *
      amp[223] + 1./3. * amp[224] + 1./9. * amp[225] + 1./9. * amp[226] + 1./3.
      * amp[227] + 1./3. * amp[228] + 1./9. * amp[229] + 1./3. * amp[230] +
      1./9. * amp[231] + 1./9. * amp[232] + 1./3. * amp[233] + 1./9. * amp[234]
      + 1./3. * amp[235] + 1./3. * amp[236] + 1./9. * amp[237] + 1./9. *
      amp[238] + 1./3. * amp[239] + 1./3. * amp[240] + 1./9. * amp[241] + 1./9.
      * amp[242] + 1./3. * amp[243] + 1./3. * amp[244] + 1./9. * amp[245] +
      1./3. * amp[246] + 1./9. * amp[247] + 1./9. * amp[248] + 1./3. * amp[249]
      + 1./9. * amp[250] + 1./3. * amp[251] + 1./3. * amp[252] + 1./9. *
      amp[253] + 1./9. * amp[254] + 1./3. * amp[255] + 1./3. * amp[256] + 1./3.
      * amp[257] + 1./3. * amp[260] + 1./3. * amp[261] + 1./3. * amp[262] +
      1./3. * amp[263] + 1./3. * amp[264] + 1./3. * amp[265] + 1./3. * amp[268]
      + 1./3. * amp[269] + 1./3. * amp[270] + 1./3. * amp[271] + 1./9. *
      amp[272] + 1./9. * amp[273] + 1./9. * amp[274] + 1./9. * amp[275] + 1./9.
      * amp[278] + 1./9. * amp[279] + 1./9. * amp[280] + 1./9. * amp[281] +
      1./9. * amp[282] + 1./9. * amp[283] + 1./9. * amp[286] + 1./9. * amp[287]
      + 1./9. * amp[288] + 1./3. * amp[289] + 1./3. * amp[290] + 1./9. *
      amp[291] + 1./9. * amp[292] + 1./3. * amp[293] + 1./3. * amp[294] + 1./9.
      * amp[295] + 1./3. * amp[296] + 1./9. * amp[297] + 1./3. * amp[298] +
      1./9. * amp[299] + 1./3. * amp[300] + 1./9. * amp[301] + 1./3. * amp[302]
      + 1./9. * amp[303] + 1./9. * amp[304] + 1./3. * amp[305] + 1./3. *
      amp[306] + 1./9. * amp[307] + 1./9. * amp[308] + 1./3. * amp[309] + 1./3.
      * amp[310] + 1./9. * amp[311] + 1./3. * amp[312] + 1./9. * amp[313] +
      1./3. * amp[314] + 1./9. * amp[315] + 1./3. * amp[316] + 1./9. * amp[317]
      + 1./3. * amp[318] + 1./9. * amp[319] + 1./9. * amp[320] + 1./9. *
      amp[321] + 1./9. * amp[323] + 1./9. * amp[324] + 1./9. * amp[325] + 1./9.
      * amp[326] + 1./9. * amp[328] + 1./9. * amp[329] + 1./9. * amp[331] +
      1./9. * amp[332] + 1./9. * amp[333] + 1./9. * amp[334] + 1./3. * amp[336]
      + 1./3. * amp[337] + 1./3. * amp[338] + 1./3. * amp[340] + 1./3. *
      amp[341] + 1./3. * amp[342] + 1./3. * amp[344] + 1./3. * amp[345] + 1./3.
      * amp[346] + 1./3. * amp[348] + 1./3. * amp[349] + 1./3. * amp[350] +
      1./3. * amp[356] + amp[358] + 1./3. * amp[359] + amp[360] + 1./3. *
      amp[361] + amp[362] + amp[366] + 1./3. * amp[367] + 1./3. * amp[372] +
      amp[374] + 1./3. * amp[375] + amp[376] + 1./3. * amp[377] + amp[378] +
      amp[382] + 1./3. * amp[383] + amp[402] - std::complex<double> (0, 1) *
      amp[403] + amp[405] + amp[406] - std::complex<double> (0, 1) * amp[407] +
      amp[410] - std::complex<double> (0, 1) * amp[411] + amp[413] + amp[414] -
      std::complex<double> (0, 1) * amp[415] + amp[416] + 1./3. * amp[417] +
      1./3. * amp[418] + amp[419] + amp[420] + 1./3. * amp[424] + amp[425] +
      1./3. * amp[426] + amp[432] + 1./3. * amp[433] + 1./3. * amp[434] +
      amp[435] + amp[436] + 1./3. * amp[440] + amp[441] + 1./3. * amp[442] +
      amp[448] + amp[449] + std::complex<double> (0, 1) * amp[450] + amp[452] +
      std::complex<double> (0, 1) * amp[455] + amp[456] + amp[457] +
      std::complex<double> (0, 1) * amp[458] + amp[460] + std::complex<double>
      (0, 1) * amp[463] + 1./3. * amp[464] + 1./3. * amp[465] + 1./3. *
      amp[466] + 1./3. * amp[468] + 1./3. * amp[469] + 1./3. * amp[470] + 1./3.
      * amp[472] + 1./3. * amp[473] + 1./3. * amp[474] + 1./3. * amp[476] +
      1./3. * amp[477] + 1./3. * amp[478] + amp[480] + 1./3. * amp[481] + 1./3.
      * amp[482] + amp[488] + amp[490] + 1./3. * amp[491] + 1./3. * amp[494] +
      amp[495] + amp[496] + 1./3. * amp[497] + 1./3. * amp[498] + amp[504] +
      amp[506] + 1./3. * amp[507] + 1./3. * amp[510] + amp[511] + 1./3. *
      amp[512] + 1./9. * amp[513] + 1./9. * amp[514] + 1./3. * amp[515] + 1./3.
      * amp[516] + 1./9. * amp[517] + 1./9. * amp[518] + 1./3. * amp[519] +
      1./3. * amp[520] + 1./9. * amp[521] + 1./3. * amp[522] + 1./9. * amp[523]
      + 1./3. * amp[524] + 1./9. * amp[525] + 1./9. * amp[526] + 1./3. *
      amp[527] + 1./3. * amp[528] + 1./9. * amp[529] + 1./9. * amp[530] + 1./3.
      * amp[531] + 1./3. * amp[532] + 1./9. * amp[533] + 1./9. * amp[534] +
      1./3. * amp[535] + 1./3. * amp[536] + 1./9. * amp[537] + 1./3. * amp[538]
      + 1./9. * amp[539] + 1./3. * amp[540] + 1./9. * amp[541] + 1./9. *
      amp[542] + 1./3. * amp[543] + amp[547] + amp[548] + 1./3. * amp[549] +
      1./3. * amp[550] + amp[551] + 1./3. * amp[553] + amp[556] + 1./3. *
      amp[557] + amp[563] + amp[564] + 1./3. * amp[565] + 1./3. * amp[566] +
      amp[567] + 1./3. * amp[569] + amp[572] + 1./3. * amp[573]);
  jamp[4] = +1./4. * (+std::complex<double> (0, 1) * amp[0] +
      std::complex<double> (0, 1) * amp[1] + amp[4] + amp[5] + amp[6] +
      std::complex<double> (0, 1) * amp[8] + std::complex<double> (0, 1) *
      amp[9] + amp[12] + amp[13] + amp[14] + 1./3. * amp[18] + 1./3. * amp[19]
      + 1./3. * amp[20] + 1./3. * amp[21] + 1./3. * amp[22] + 1./3. * amp[23] +
      1./3. * amp[26] + 1./3. * amp[27] + 1./3. * amp[28] + 1./3. * amp[29] +
      1./3. * amp[30] + 1./3. * amp[31] + 1./3. * amp[36] + amp[37] + 1./3. *
      amp[38] + amp[39] + 1./3. * amp[44] + amp[45] + amp[46] + 1./3. * amp[47]
      + 1./3. * amp[52] + amp[53] + 1./3. * amp[54] + amp[55] + 1./3. * amp[60]
      + amp[61] + amp[62] + 1./3. * amp[63] + amp[82] + amp[83] -
      std::complex<double> (0, 1) * amp[84] - std::complex<double> (0, 1) *
      amp[85] + amp[87] + amp[90] + amp[91] - std::complex<double> (0, 1) *
      amp[92] - std::complex<double> (0, 1) * amp[93] + amp[95] + 1./3. *
      amp[98] + 1./3. * amp[99] + 1./3. * amp[100] + 1./3. * amp[101] + 1./3. *
      amp[102] + 1./3. * amp[103] + 1./3. * amp[106] + 1./3. * amp[107] + 1./3.
      * amp[108] + 1./3. * amp[109] + 1./3. * amp[110] + 1./3. * amp[111] +
      1./9. * amp[114] + 1./9. * amp[115] + 1./9. * amp[116] + 1./9. * amp[117]
      + 1./9. * amp[118] + 1./9. * amp[119] + 1./9. * amp[122] + 1./9. *
      amp[123] + 1./9. * amp[124] + 1./9. * amp[125] + 1./9. * amp[126] + 1./9.
      * amp[127] + 1./9. * amp[128] + 1./3. * amp[129] + 1./3. * amp[130] +
      1./9. * amp[131] + 1./9. * amp[132] + 1./3. * amp[133] + 1./9. * amp[134]
      + 1./3. * amp[135] + 1./3. * amp[136] + 1./9. * amp[137] + 1./3. *
      amp[138] + 1./9. * amp[139] + 1./9. * amp[140] + 1./3. * amp[141] + 1./3.
      * amp[142] + 1./9. * amp[143] + 1./9. * amp[144] + 1./3. * amp[145] +
      1./3. * amp[146] + 1./9. * amp[147] + 1./9. * amp[148] + 1./3. * amp[149]
      + 1./9. * amp[150] + 1./3. * amp[151] + 1./3. * amp[152] + 1./9. *
      amp[153] + 1./3. * amp[154] + 1./9. * amp[155] + 1./9. * amp[156] + 1./3.
      * amp[157] + 1./3. * amp[158] + 1./9. * amp[159] + 1./9. * amp[160] +
      1./9. * amp[161] + 1./9. * amp[164] + 1./9. * amp[165] + 1./9. * amp[166]
      + 1./9. * amp[167] + 1./9. * amp[168] + 1./9. * amp[169] + 1./9. *
      amp[172] + 1./9. * amp[173] + 1./9. * amp[174] + 1./9. * amp[175] + 1./3.
      * amp[176] + 1./3. * amp[177] + 1./3. * amp[178] + 1./3. * amp[179] +
      1./3. * amp[182] + 1./3. * amp[183] + 1./3. * amp[184] + 1./3. * amp[185]
      + 1./3. * amp[186] + 1./3. * amp[187] + 1./3. * amp[190] + 1./3. *
      amp[191] - std::complex<double> (0, 1) * amp[192] - std::complex<double>
      (0, 1) * amp[193] + amp[194] + amp[195] + amp[199] - std::complex<double>
      (0, 1) * amp[200] - std::complex<double> (0, 1) * amp[201] + amp[202] +
      amp[203] + amp[207] + 1./3. * amp[224] + amp[225] + amp[226] + 1./3. *
      amp[227] + amp[232] + 1./3. * amp[233] + amp[234] + 1./3. * amp[235] +
      1./3. * amp[240] + amp[241] + amp[242] + 1./3. * amp[243] + amp[248] +
      1./3. * amp[249] + amp[250] + 1./3. * amp[251] + 1./3. * amp[256] + 1./3.
      * amp[257] + 1./3. * amp[260] + 1./3. * amp[261] + 1./3. * amp[262] +
      1./3. * amp[263] + 1./3. * amp[264] + 1./3. * amp[265] + 1./3. * amp[268]
      + 1./3. * amp[269] + 1./3. * amp[270] + 1./3. * amp[271] + amp[272] +
      amp[273] + std::complex<double> (0, 1) * amp[276] + std::complex<double>
      (0, 1) * amp[277] + amp[278] + amp[280] + amp[281] + std::complex<double>
      (0, 1) * amp[284] + std::complex<double> (0, 1) * amp[285] + amp[286] +
      1./3. * amp[293] + 1./3. * amp[294] + amp[295] + amp[299] + 1./3. *
      amp[300] + amp[301] + 1./3. * amp[302] + amp[303] + 1./3. * amp[309] +
      1./3. * amp[310] + amp[311] + amp[315] + 1./3. * amp[316] + amp[317] +
      1./3. * amp[318] + amp[319] - std::complex<double> (0, 1) * amp[322] +
      amp[323] + amp[325] + amp[326] - std::complex<double> (0, 1) * amp[327] -
      std::complex<double> (0, 1) * amp[330] + amp[331] + amp[333] + amp[334] -
      std::complex<double> (0, 1) * amp[335] + 1./3. * amp[352] + amp[353] +
      amp[354] + 1./3. * amp[355] + amp[357] + 1./3. * amp[363] + amp[364] +
      1./3. * amp[365] + 1./3. * amp[368] + amp[369] + amp[370] + 1./3. *
      amp[371] + amp[373] + 1./3. * amp[379] + amp[380] + 1./3. * amp[381] +
      1./3. * amp[384] + 1./3. * amp[385] + 1./3. * amp[387] + 1./3. * amp[388]
      + 1./3. * amp[389] + 1./3. * amp[390] + 1./3. * amp[392] + 1./3. *
      amp[393] + 1./3. * amp[395] + 1./3. * amp[396] + 1./3. * amp[397] + 1./3.
      * amp[398] + amp[400] + amp[401] + std::complex<double> (0, 1) * amp[403]
      + amp[404] + std::complex<double> (0, 1) * amp[407] + amp[408] + amp[409]
      + std::complex<double> (0, 1) * amp[411] + amp[412] +
      std::complex<double> (0, 1) * amp[415] + 1./9. * amp[416] + 1./3. *
      amp[417] + 1./3. * amp[418] + 1./9. * amp[419] + 1./9. * amp[420] + 1./3.
      * amp[421] + 1./3. * amp[422] + 1./9. * amp[423] + 1./3. * amp[424] +
      1./9. * amp[425] + 1./3. * amp[426] + 1./9. * amp[427] + 1./3. * amp[428]
      + 1./9. * amp[429] + 1./3. * amp[430] + 1./9. * amp[431] + 1./9. *
      amp[432] + 1./3. * amp[433] + 1./3. * amp[434] + 1./9. * amp[435] + 1./9.
      * amp[436] + 1./3. * amp[437] + 1./3. * amp[438] + 1./9. * amp[439] +
      1./3. * amp[440] + 1./9. * amp[441] + 1./3. * amp[442] + 1./9. * amp[443]
      + 1./3. * amp[444] + 1./9. * amp[445] + 1./3. * amp[446] + 1./9. *
      amp[447] + 1./9. * amp[448] + 1./9. * amp[449] + 1./9. * amp[451] + 1./9.
      * amp[452] + 1./9. * amp[453] + 1./9. * amp[454] + 1./9. * amp[456] +
      1./9. * amp[457] + 1./9. * amp[459] + 1./9. * amp[460] + 1./9. * amp[461]
      + 1./9. * amp[462] + 1./3. * amp[464] + 1./3. * amp[465] + 1./3. *
      amp[466] + 1./3. * amp[468] + 1./3. * amp[469] + 1./3. * amp[470] + 1./3.
      * amp[472] + 1./3. * amp[473] + 1./3. * amp[474] + 1./3. * amp[476] +
      1./3. * amp[477] + 1./3. * amp[478] + 1./9. * amp[480] + 1./3. * amp[481]
      + 1./3. * amp[482] + 1./9. * amp[483] + 1./9. * amp[484] + 1./3. *
      amp[485] + 1./3. * amp[486] + 1./9. * amp[487] + 1./9. * amp[488] + 1./3.
      * amp[489] + 1./9. * amp[490] + 1./3. * amp[491] + 1./9. * amp[492] +
      1./3. * amp[493] + 1./3. * amp[494] + 1./9. * amp[495] + 1./9. * amp[496]
      + 1./3. * amp[497] + 1./3. * amp[498] + 1./9. * amp[499] + 1./9. *
      amp[500] + 1./3. * amp[501] + 1./3. * amp[502] + 1./9. * amp[503] + 1./9.
      * amp[504] + 1./3. * amp[505] + 1./9. * amp[506] + 1./3. * amp[507] +
      1./9. * amp[508] + 1./3. * amp[509] + 1./3. * amp[510] + 1./9. * amp[511]
      + amp[514] + 1./3. * amp[516] + amp[517] + amp[518] + 1./3. * amp[519] +
      1./3. * amp[520] + amp[526] + 1./3. * amp[527] + amp[530] + 1./3. *
      amp[532] + amp[533] + amp[534] + 1./3. * amp[535] + 1./3. * amp[536] +
      amp[542] + 1./3. * amp[543] + amp[544] + 1./3. * amp[545] + 1./3. *
      amp[546] + amp[552] + amp[554] + 1./3. * amp[555] + 1./3. * amp[558] +
      amp[559] + amp[560] + 1./3. * amp[561] + 1./3. * amp[562] + amp[568] +
      amp[570] + 1./3. * amp[571] + 1./3. * amp[574] + amp[575]);
  jamp[5] = +1./4. * (-1./3. * amp[2] - 1./3. * amp[3] - 1./3. * amp[4] - 1./3.
      * amp[5] - 1./3. * amp[6] - 1./3. * amp[7] - 1./3. * amp[10] - 1./3. *
      amp[11] - 1./3. * amp[12] - 1./3. * amp[13] - 1./3. * amp[14] - 1./3. *
      amp[15] - std::complex<double> (0, 1) * amp[16] - std::complex<double>
      (0, 1) * amp[17] - amp[20] - amp[21] - amp[22] - std::complex<double> (0,
      1) * amp[24] - std::complex<double> (0, 1) * amp[25] - amp[28] - amp[29]
      - amp[30] - 1./3. * amp[40] - amp[41] - 1./3. * amp[42] - amp[43] -
      amp[44] - 1./3. * amp[45] - 1./3. * amp[46] - amp[47] - 1./3. * amp[56] -
      amp[57] - 1./3. * amp[58] - amp[59] - amp[60] - 1./3. * amp[61] - 1./3. *
      amp[62] - amp[63] + std::complex<double> (0, 1) * amp[66] +
      std::complex<double> (0, 1) * amp[67] - amp[68] - amp[69] - amp[71] +
      std::complex<double> (0, 1) * amp[74] + std::complex<double> (0, 1) *
      amp[75] - amp[76] - amp[77] - amp[79] + std::complex<double> (0, 1) *
      amp[96] + std::complex<double> (0, 1) * amp[97] - amp[98] - amp[99] -
      amp[103] + std::complex<double> (0, 1) * amp[104] + std::complex<double>
      (0, 1) * amp[105] - amp[106] - amp[107] - amp[111] - 1./3. * amp[128] -
      amp[129] - amp[130] - 1./3. * amp[131] - amp[136] - 1./3. * amp[137] -
      amp[138] - 1./3. * amp[139] - 1./3. * amp[144] - amp[145] - amp[146] -
      1./3. * amp[147] - amp[152] - 1./3. * amp[153] - amp[154] - 1./3. *
      amp[155] - 1./3. * amp[160] - 1./3. * amp[161] - 1./3. * amp[164] - 1./3.
      * amp[165] - 1./3. * amp[166] - 1./3. * amp[167] - 1./3. * amp[168] -
      1./3. * amp[169] - 1./3. * amp[172] - 1./3. * amp[173] - 1./3. * amp[174]
      - 1./3. * amp[175] - amp[176] - amp[177] - std::complex<double> (0, 1) *
      amp[180] - std::complex<double> (0, 1) * amp[181] - amp[182] - amp[184] -
      amp[185] - std::complex<double> (0, 1) * amp[188] - std::complex<double>
      (0, 1) * amp[189] - amp[190] - 1./3. * amp[194] - 1./3. * amp[195] -
      1./3. * amp[196] - 1./3. * amp[197] - 1./3. * amp[198] - 1./3. * amp[199]
      - 1./3. * amp[202] - 1./3. * amp[203] - 1./3. * amp[204] - 1./3. *
      amp[205] - 1./3. * amp[206] - 1./3. * amp[207] - 1./9. * amp[210] - 1./9.
      * amp[211] - 1./9. * amp[212] - 1./9. * amp[213] - 1./9. * amp[214] -
      1./9. * amp[215] - 1./9. * amp[218] - 1./9. * amp[219] - 1./9. * amp[220]
      - 1./9. * amp[221] - 1./9. * amp[222] - 1./9. * amp[223] - 1./9. *
      amp[224] - 1./3. * amp[225] - 1./3. * amp[226] - 1./9. * amp[227] - 1./9.
      * amp[228] - 1./3. * amp[229] - 1./9. * amp[230] - 1./3. * amp[231] -
      1./3. * amp[232] - 1./9. * amp[233] - 1./3. * amp[234] - 1./9. * amp[235]
      - 1./9. * amp[236] - 1./3. * amp[237] - 1./3. * amp[238] - 1./9. *
      amp[239] - 1./9. * amp[240] - 1./3. * amp[241] - 1./3. * amp[242] - 1./9.
      * amp[243] - 1./9. * amp[244] - 1./3. * amp[245] - 1./9. * amp[246] -
      1./3. * amp[247] - 1./3. * amp[248] - 1./9. * amp[249] - 1./3. * amp[250]
      - 1./9. * amp[251] - 1./9. * amp[252] - 1./3. * amp[253] - 1./3. *
      amp[254] - 1./9. * amp[255] - 1./9. * amp[256] - 1./9. * amp[257] - 1./9.
      * amp[260] - 1./9. * amp[261] - 1./9. * amp[262] - 1./9. * amp[263] -
      1./9. * amp[264] - 1./9. * amp[265] - 1./9. * amp[268] - 1./9. * amp[269]
      - 1./9. * amp[270] - 1./9. * amp[271] - 1./3. * amp[272] - 1./3. *
      amp[273] - 1./3. * amp[274] - 1./3. * amp[275] - 1./3. * amp[278] - 1./3.
      * amp[279] - 1./3. * amp[280] - 1./3. * amp[281] - 1./3. * amp[282] -
      1./3. * amp[283] - 1./3. * amp[286] - 1./3. * amp[287] - 1./3. * amp[292]
      - amp[294] - 1./3. * amp[295] - amp[296] - 1./3. * amp[297] - amp[298] -
      amp[302] - 1./3. * amp[303] - 1./3. * amp[308] - amp[310] - 1./3. *
      amp[311] - amp[312] - 1./3. * amp[313] - amp[314] - amp[318] - 1./3. *
      amp[319] - amp[338] + std::complex<double> (0, 1) * amp[339] - amp[341] -
      amp[342] + std::complex<double> (0, 1) * amp[343] - amp[346] +
      std::complex<double> (0, 1) * amp[347] - amp[349] - amp[350] +
      std::complex<double> (0, 1) * amp[351] - 1./9. * amp[352] - 1./3. *
      amp[353] - 1./3. * amp[354] - 1./9. * amp[355] - 1./9. * amp[356] - 1./3.
      * amp[357] - 1./3. * amp[358] - 1./9. * amp[359] - 1./3. * amp[360] -
      1./9. * amp[361] - 1./3. * amp[362] - 1./9. * amp[363] - 1./3. * amp[364]
      - 1./9. * amp[365] - 1./3. * amp[366] - 1./9. * amp[367] - 1./9. *
      amp[368] - 1./3. * amp[369] - 1./3. * amp[370] - 1./9. * amp[371] - 1./9.
      * amp[372] - 1./3. * amp[373] - 1./3. * amp[374] - 1./9. * amp[375] -
      1./3. * amp[376] - 1./9. * amp[377] - 1./3. * amp[378] - 1./9. * amp[379]
      - 1./3. * amp[380] - 1./9. * amp[381] - 1./3. * amp[382] - 1./9. *
      amp[383] - 1./9. * amp[384] - 1./9. * amp[385] - 1./9. * amp[387] - 1./9.
      * amp[388] - 1./9. * amp[389] - 1./9. * amp[390] - 1./9. * amp[392] -
      1./9. * amp[393] - 1./9. * amp[395] - 1./9. * amp[396] - 1./9. * amp[397]
      - 1./9. * amp[398] - 1./3. * amp[400] - 1./3. * amp[401] - 1./3. *
      amp[402] - 1./3. * amp[404] - 1./3. * amp[405] - 1./3. * amp[406] - 1./3.
      * amp[408] - 1./3. * amp[409] - 1./3. * amp[410] - 1./3. * amp[412] -
      1./3. * amp[413] - 1./3. * amp[414] - 1./3. * amp[416] - amp[417] -
      amp[418] - 1./3. * amp[419] - amp[421] - 1./3. * amp[427] - amp[428] -
      1./3. * amp[429] - 1./3. * amp[432] - amp[433] - amp[434] - 1./3. *
      amp[435] - amp[437] - 1./3. * amp[443] - amp[444] - 1./3. * amp[445] -
      1./3. * amp[448] - 1./3. * amp[449] - 1./3. * amp[451] - 1./3. * amp[452]
      - 1./3. * amp[453] - 1./3. * amp[454] - 1./3. * amp[456] - 1./3. *
      amp[457] - 1./3. * amp[459] - 1./3. * amp[460] - 1./3. * amp[461] - 1./3.
      * amp[462] - amp[464] - amp[465] - std::complex<double> (0, 1) * amp[467]
      - amp[468] - std::complex<double> (0, 1) * amp[471] - amp[472] - amp[473]
      - std::complex<double> (0, 1) * amp[475] - amp[476] -
      std::complex<double> (0, 1) * amp[479] - 1./3. * amp[480] - 1./9. *
      amp[481] - 1./9. * amp[482] - 1./3. * amp[483] - 1./3. * amp[484] - 1./9.
      * amp[485] - 1./9. * amp[486] - 1./3. * amp[487] - 1./3. * amp[488] -
      1./9. * amp[489] - 1./3. * amp[490] - 1./9. * amp[491] - 1./3. * amp[492]
      - 1./9. * amp[493] - 1./9. * amp[494] - 1./3. * amp[495] - 1./3. *
      amp[496] - 1./9. * amp[497] - 1./9. * amp[498] - 1./3. * amp[499] - 1./3.
      * amp[500] - 1./9. * amp[501] - 1./9. * amp[502] - 1./3. * amp[503] -
      1./3. * amp[504] - 1./9. * amp[505] - 1./3. * amp[506] - 1./9. * amp[507]
      - 1./3. * amp[508] - 1./9. * amp[509] - 1./9. * amp[510] - 1./3. *
      amp[511] - amp[512] - 1./3. * amp[513] - 1./3. * amp[514] - amp[520] -
      amp[522] - 1./3. * amp[523] - 1./3. * amp[526] - amp[527] - amp[528] -
      1./3. * amp[529] - 1./3. * amp[530] - amp[536] - amp[538] - 1./3. *
      amp[539] - 1./3. * amp[542] - amp[543] - amp[546] - 1./3. * amp[548] -
      amp[549] - amp[550] - 1./3. * amp[551] - 1./3. * amp[552] - amp[558] -
      1./3. * amp[559] - amp[562] - 1./3. * amp[564] - amp[565] - amp[566] -
      1./3. * amp[567] - 1./3. * amp[568] - amp[574] - 1./3. * amp[575]);

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



