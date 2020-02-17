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
// Process: u u > e+ e- g g u u WEIGHTED<=8 / h
// Process: u u > mu+ mu- g g u u WEIGHTED<=8 / h
// Process: c c > e+ e- g g c c WEIGHTED<=8 / h
// Process: c c > mu+ mu- g g c c WEIGHTED<=8 / h

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
  const int denominators[nprocesses] = {144}; 

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
        t[0] = matrix_uu_epemgguu_no_h(); 

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
      t[0] = matrix_uu_epemgguu_no_h(); 

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
    return matrix_element[0] * 2; 
  }
  else if(id1 == 2 && id2 == 2)
  {
    // Add matrix elements for processes with beams (2, 2)
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
  FFV1_1(w[7], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[11]); 
  FFV1P0_3(w[0], w[10], pars->GC_11, pars->ZERO, pars->ZERO, w[12]); 
  FFV1P0_3(w[1], w[10], pars->GC_11, pars->ZERO, pars->ZERO, w[13]); 
  FFV1_2(w[0], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[14]); 
  FFV1_2(w[1], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[15]); 
  FFV1_1(w[7], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[16]); 
  FFV1_1(w[6], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[17]); 
  FFV1P0_3(w[0], w[16], pars->GC_11, pars->ZERO, pars->ZERO, w[18]); 
  FFV1P0_3(w[1], w[16], pars->GC_11, pars->ZERO, pars->ZERO, w[19]); 
  FFV1_2(w[0], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[20]); 
  FFV1P0_3(w[20], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[21]); 
  FFV1P0_3(w[20], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[22]); 
  FFV1_2(w[1], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[23]); 
  FFV1P0_3(w[23], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[24]); 
  FFV1P0_3(w[23], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[25]); 
  FFV2_4_3(w[2], w[3], -pars->GC_51, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ,
      w[26]);
  FFV2_5_1(w[7], w[26], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[27]);
  FFV2_5_2(w[0], w[26], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[28]);
  FFV2_5_2(w[1], w[26], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[29]);
  FFV2_5_1(w[6], w[26], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[30]);
  FFV1P0_3(w[0], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[31]); 
  VVV1P0_1(w[8], w[31], pars->GC_10, pars->ZERO, pars->ZERO, w[32]); 
  FFV1_2(w[1], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[33]); 
  FFV1_1(w[7], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[34]); 
  FFV1P0_3(w[1], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[35]); 
  VVV1P0_1(w[8], w[35], pars->GC_10, pars->ZERO, pars->ZERO, w[36]); 
  FFV1_2(w[0], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[37]); 
  FFV1_1(w[7], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[38]); 
  FFV1P0_3(w[0], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[39]); 
  FFV1_2(w[1], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[40]); 
  VVV1P0_1(w[8], w[39], pars->GC_10, pars->ZERO, pars->ZERO, w[41]); 
  FFV1_1(w[6], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[42]); 
  FFV1P0_3(w[1], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[43]); 
  FFV1_2(w[0], w[43], pars->GC_11, pars->ZERO, pars->ZERO, w[44]); 
  VVV1P0_1(w[8], w[43], pars->GC_10, pars->ZERO, pars->ZERO, w[45]); 
  FFV1_1(w[6], w[43], pars->GC_11, pars->ZERO, pars->ZERO, w[46]); 
  FFV1_1(w[6], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[47]); 
  FFV1_1(w[7], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[48]); 
  FFV1_1(w[47], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[49]); 
  FFV1P0_3(w[0], w[48], pars->GC_11, pars->ZERO, pars->ZERO, w[50]); 
  FFV1P0_3(w[1], w[48], pars->GC_11, pars->ZERO, pars->ZERO, w[51]); 
  FFV1P0_3(w[0], w[47], pars->GC_11, pars->ZERO, pars->ZERO, w[52]); 
  FFV1_1(w[48], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[53]); 
  FFV1P0_3(w[1], w[47], pars->GC_11, pars->ZERO, pars->ZERO, w[54]); 
  FFV2_5_1(w[47], w[26], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[55]);
  FFV2_5_1(w[48], w[26], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[56]);
  FFV1_2(w[0], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[57]); 
  FFV1P0_3(w[57], w[47], pars->GC_11, pars->ZERO, pars->ZERO, w[58]); 
  FFV1P0_3(w[57], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[59]); 
  FFV1_2(w[57], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[60]); 
  FFV2_5_2(w[57], w[26], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[61]);
  FFV1_1(w[47], w[43], pars->GC_11, pars->ZERO, pars->ZERO, w[62]); 
  FFV1_2(w[1], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[63]); 
  FFV1P0_3(w[63], w[47], pars->GC_11, pars->ZERO, pars->ZERO, w[64]); 
  FFV1P0_3(w[63], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[65]); 
  FFV1_2(w[63], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[66]); 
  FFV2_5_2(w[63], w[26], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[67]);
  FFV1_1(w[47], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[68]); 
  FFV1_1(w[47], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[69]); 
  FFV1P0_3(w[0], w[69], pars->GC_11, pars->ZERO, pars->ZERO, w[70]); 
  FFV1P0_3(w[1], w[69], pars->GC_11, pars->ZERO, pars->ZERO, w[71]); 
  VVV1P0_1(w[52], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[72]); 
  FFV1_2(w[1], w[52], pars->GC_11, pars->ZERO, pars->ZERO, w[73]); 
  FFV1_1(w[7], w[52], pars->GC_11, pars->ZERO, pars->ZERO, w[74]); 
  VVV1P0_1(w[54], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[75]); 
  FFV1_2(w[0], w[54], pars->GC_11, pars->ZERO, pars->ZERO, w[76]); 
  FFV1_1(w[7], w[54], pars->GC_11, pars->ZERO, pars->ZERO, w[77]); 
  VVV1P0_1(w[5], w[39], pars->GC_10, pars->ZERO, pars->ZERO, w[78]); 
  VVV1P0_1(w[5], w[43], pars->GC_10, pars->ZERO, pars->ZERO, w[79]); 
  FFV1_1(w[7], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[80]); 
  FFV1_1(w[6], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[81]); 
  FFV1_1(w[80], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[82]); 
  FFV1P0_3(w[0], w[81], pars->GC_11, pars->ZERO, pars->ZERO, w[83]); 
  FFV1P0_3(w[1], w[81], pars->GC_11, pars->ZERO, pars->ZERO, w[84]); 
  FFV1P0_3(w[0], w[80], pars->GC_11, pars->ZERO, pars->ZERO, w[85]); 
  FFV1_1(w[81], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[86]); 
  FFV1P0_3(w[1], w[80], pars->GC_11, pars->ZERO, pars->ZERO, w[87]); 
  FFV2_5_1(w[80], w[26], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[88]);
  FFV2_5_1(w[81], w[26], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[89]);
  FFV1P0_3(w[57], w[80], pars->GC_11, pars->ZERO, pars->ZERO, w[90]); 
  FFV1P0_3(w[57], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[91]); 
  FFV1_1(w[80], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[92]); 
  FFV1P0_3(w[63], w[80], pars->GC_11, pars->ZERO, pars->ZERO, w[93]); 
  FFV1P0_3(w[63], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[94]); 
  FFV1_1(w[80], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[95]); 
  FFV1_1(w[80], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[96]); 
  FFV1P0_3(w[0], w[96], pars->GC_11, pars->ZERO, pars->ZERO, w[97]); 
  FFV1P0_3(w[1], w[96], pars->GC_11, pars->ZERO, pars->ZERO, w[98]); 
  VVV1P0_1(w[85], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[99]); 
  FFV1_2(w[1], w[85], pars->GC_11, pars->ZERO, pars->ZERO, w[100]); 
  FFV1_1(w[6], w[85], pars->GC_11, pars->ZERO, pars->ZERO, w[101]); 
  VVV1P0_1(w[87], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[102]); 
  FFV1_2(w[0], w[87], pars->GC_11, pars->ZERO, pars->ZERO, w[103]); 
  FFV1_1(w[6], w[87], pars->GC_11, pars->ZERO, pars->ZERO, w[104]); 
  VVV1P0_1(w[5], w[31], pars->GC_10, pars->ZERO, pars->ZERO, w[105]); 
  VVV1P0_1(w[5], w[35], pars->GC_10, pars->ZERO, pars->ZERO, w[106]); 
  FFV1_2(w[0], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[107]); 
  FFV1P0_3(w[107], w[81], pars->GC_11, pars->ZERO, pars->ZERO, w[108]); 
  FFV1_2(w[107], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[109]); 
  FFV1P0_3(w[107], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[110]); 
  FFV2_5_2(w[107], w[26], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[111]);
  FFV1_2(w[107], w[43], pars->GC_11, pars->ZERO, pars->ZERO, w[112]); 
  FFV1P0_3(w[107], w[48], pars->GC_11, pars->ZERO, pars->ZERO, w[113]); 
  FFV1P0_3(w[107], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[114]); 
  FFV1_2(w[107], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[115]); 
  FFV1_2(w[107], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[116]); 
  FFV1P0_3(w[116], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[117]); 
  FFV1P0_3(w[116], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[118]); 
  VVV1P0_1(w[114], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[119]); 
  FFV1_2(w[1], w[114], pars->GC_11, pars->ZERO, pars->ZERO, w[120]); 
  FFV1_1(w[7], w[114], pars->GC_11, pars->ZERO, pars->ZERO, w[121]); 
  VVV1P0_1(w[110], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[122]); 
  FFV1_2(w[1], w[110], pars->GC_11, pars->ZERO, pars->ZERO, w[123]); 
  FFV1_1(w[6], w[110], pars->GC_11, pars->ZERO, pars->ZERO, w[124]); 
  FFV1_2(w[1], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[125]); 
  FFV1P0_3(w[125], w[81], pars->GC_11, pars->ZERO, pars->ZERO, w[126]); 
  FFV1_2(w[125], w[9], pars->GC_2, pars->ZERO, pars->ZERO, w[127]); 
  FFV1P0_3(w[125], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[128]); 
  FFV2_5_2(w[125], w[26], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[129]);
  FFV1_2(w[125], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[130]); 
  FFV1P0_3(w[125], w[48], pars->GC_11, pars->ZERO, pars->ZERO, w[131]); 
  FFV1P0_3(w[125], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[132]); 
  FFV1_2(w[125], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[133]); 
  FFV1_2(w[125], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[134]); 
  FFV1P0_3(w[134], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[135]); 
  FFV1P0_3(w[134], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[136]); 
  VVV1P0_1(w[132], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[137]); 
  FFV1_2(w[0], w[132], pars->GC_11, pars->ZERO, pars->ZERO, w[138]); 
  FFV1_1(w[7], w[132], pars->GC_11, pars->ZERO, pars->ZERO, w[139]); 
  VVV1P0_1(w[128], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[140]); 
  FFV1_2(w[0], w[128], pars->GC_11, pars->ZERO, pars->ZERO, w[141]); 
  FFV1_1(w[6], w[128], pars->GC_11, pars->ZERO, pars->ZERO, w[142]); 
  FFV1_1(w[81], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[143]); 
  FFV1P0_3(w[0], w[143], pars->GC_11, pars->ZERO, pars->ZERO, w[144]); 
  FFV1P0_3(w[1], w[143], pars->GC_11, pars->ZERO, pars->ZERO, w[145]); 
  VVV1P0_1(w[4], w[83], pars->GC_10, pars->ZERO, pars->ZERO, w[146]); 
  FFV1_1(w[11], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[147]); 
  FFV1_2(w[15], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[148]); 
  VVV1P0_1(w[4], w[84], pars->GC_10, pars->ZERO, pars->ZERO, w[149]); 
  FFV1_2(w[14], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[150]); 
  FFV1_1(w[27], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[151]); 
  FFV1_2(w[29], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[152]); 
  FFV1_2(w[28], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[153]); 
  VVV1P0_1(w[4], w[39], pars->GC_10, pars->ZERO, pars->ZERO, w[154]); 
  FFV1_1(w[81], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[155]); 
  VVV1P0_1(w[4], w[43], pars->GC_10, pars->ZERO, pars->ZERO, w[156]); 
  FFV1_1(w[81], w[43], pars->GC_11, pars->ZERO, pars->ZERO, w[157]); 
  FFV1_1(w[48], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[158]); 
  FFV1P0_3(w[0], w[158], pars->GC_11, pars->ZERO, pars->ZERO, w[159]); 
  FFV1P0_3(w[1], w[158], pars->GC_11, pars->ZERO, pars->ZERO, w[160]); 
  VVV1P0_1(w[4], w[50], pars->GC_10, pars->ZERO, pars->ZERO, w[161]); 
  FFV1_1(w[17], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[162]); 
  VVV1P0_1(w[4], w[51], pars->GC_10, pars->ZERO, pars->ZERO, w[163]); 
  FFV1_1(w[30], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[164]); 
  VVV1P0_1(w[4], w[31], pars->GC_10, pars->ZERO, pars->ZERO, w[165]); 
  FFV1_1(w[48], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[166]); 
  VVV1P0_1(w[4], w[35], pars->GC_10, pars->ZERO, pars->ZERO, w[167]); 
  FFV1_1(w[48], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[168]); 
  FFV1_2(w[57], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[169]); 
  FFV1P0_3(w[169], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[170]); 
  FFV1P0_3(w[169], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[171]); 
  VVV1P0_1(w[4], w[91], pars->GC_10, pars->ZERO, pars->ZERO, w[172]); 
  VVV1P0_1(w[4], w[59], pars->GC_10, pars->ZERO, pars->ZERO, w[173]); 
  FFV1_2(w[57], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[174]); 
  FFV1_2(w[57], w[43], pars->GC_11, pars->ZERO, pars->ZERO, w[175]); 
  FFV1_2(w[63], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[176]); 
  FFV1P0_3(w[176], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[177]); 
  FFV1P0_3(w[176], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[178]); 
  VVV1P0_1(w[4], w[94], pars->GC_10, pars->ZERO, pars->ZERO, w[179]); 
  VVV1P0_1(w[4], w[65], pars->GC_10, pars->ZERO, pars->ZERO, w[180]); 
  FFV1_2(w[63], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[181]); 
  FFV1_2(w[63], w[39], pars->GC_11, pars->ZERO, pars->ZERO, w[182]); 
  VVV1P0_1(w[165], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[183]); 
  FFV1_2(w[1], w[165], pars->GC_11, pars->ZERO, pars->ZERO, w[184]); 
  FFV1_1(w[7], w[165], pars->GC_11, pars->ZERO, pars->ZERO, w[185]); 
  VVV1P0_1(w[4], w[105], pars->GC_10, pars->ZERO, pars->ZERO, w[186]); 
  FFV1_2(w[33], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[187]); 
  FFV1_1(w[34], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[188]); 
  VVVV1P0_1(w[4], w[5], w[31], pars->GC_12, pars->ZERO, pars->ZERO, w[189]); 
  VVVV3P0_1(w[4], w[5], w[31], pars->GC_12, pars->ZERO, pars->ZERO, w[190]); 
  VVVV4P0_1(w[4], w[5], w[31], pars->GC_12, pars->ZERO, pars->ZERO, w[191]); 
  VVV1P0_1(w[167], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[192]); 
  FFV1_2(w[0], w[167], pars->GC_11, pars->ZERO, pars->ZERO, w[193]); 
  FFV1_1(w[7], w[167], pars->GC_11, pars->ZERO, pars->ZERO, w[194]); 
  VVV1P0_1(w[4], w[106], pars->GC_10, pars->ZERO, pars->ZERO, w[195]); 
  FFV1_2(w[37], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[196]); 
  FFV1_1(w[38], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[197]); 
  VVVV1P0_1(w[4], w[5], w[35], pars->GC_12, pars->ZERO, pars->ZERO, w[198]); 
  VVVV3P0_1(w[4], w[5], w[35], pars->GC_12, pars->ZERO, pars->ZERO, w[199]); 
  VVVV4P0_1(w[4], w[5], w[35], pars->GC_12, pars->ZERO, pars->ZERO, w[200]); 
  VVV1P0_1(w[154], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[201]); 
  FFV1_2(w[1], w[154], pars->GC_11, pars->ZERO, pars->ZERO, w[202]); 
  FFV1_1(w[6], w[154], pars->GC_11, pars->ZERO, pars->ZERO, w[203]); 
  VVV1P0_1(w[4], w[78], pars->GC_10, pars->ZERO, pars->ZERO, w[204]); 
  FFV1_2(w[40], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[205]); 
  FFV1_1(w[42], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[206]); 
  VVVV1P0_1(w[4], w[5], w[39], pars->GC_12, pars->ZERO, pars->ZERO, w[207]); 
  VVVV3P0_1(w[4], w[5], w[39], pars->GC_12, pars->ZERO, pars->ZERO, w[208]); 
  VVVV4P0_1(w[4], w[5], w[39], pars->GC_12, pars->ZERO, pars->ZERO, w[209]); 
  VVV1P0_1(w[156], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[210]); 
  FFV1_2(w[0], w[156], pars->GC_11, pars->ZERO, pars->ZERO, w[211]); 
  FFV1_1(w[6], w[156], pars->GC_11, pars->ZERO, pars->ZERO, w[212]); 
  VVV1P0_1(w[4], w[79], pars->GC_10, pars->ZERO, pars->ZERO, w[213]); 
  FFV1_2(w[44], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[214]); 
  FFV1_1(w[46], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[215]); 
  VVVV1P0_1(w[4], w[5], w[43], pars->GC_12, pars->ZERO, pars->ZERO, w[216]); 
  VVVV3P0_1(w[4], w[5], w[43], pars->GC_12, pars->ZERO, pars->ZERO, w[217]); 
  VVVV4P0_1(w[4], w[5], w[43], pars->GC_12, pars->ZERO, pars->ZERO, w[218]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[1], w[11], w[12], pars->GC_11, amp[0]); 
  FFV1_0(w[0], w[11], w[13], pars->GC_11, amp[1]); 
  FFV1_0(w[14], w[7], w[13], pars->GC_11, amp[2]); 
  FFV1_0(w[15], w[7], w[12], pars->GC_11, amp[3]); 
  FFV1_0(w[1], w[17], w[18], pars->GC_11, amp[4]); 
  FFV1_0(w[0], w[17], w[19], pars->GC_11, amp[5]); 
  FFV1_0(w[14], w[6], w[19], pars->GC_11, amp[6]); 
  FFV1_0(w[15], w[6], w[18], pars->GC_11, amp[7]); 
  FFV1_0(w[1], w[17], w[21], pars->GC_11, amp[8]); 
  FFV1_0(w[1], w[11], w[22], pars->GC_11, amp[9]); 
  FFV1_0(w[15], w[7], w[22], pars->GC_11, amp[10]); 
  FFV1_0(w[15], w[6], w[21], pars->GC_11, amp[11]); 
  FFV1_0(w[0], w[17], w[24], pars->GC_11, amp[12]); 
  FFV1_0(w[0], w[11], w[25], pars->GC_11, amp[13]); 
  FFV1_0(w[14], w[7], w[25], pars->GC_11, amp[14]); 
  FFV1_0(w[14], w[6], w[24], pars->GC_11, amp[15]); 
  FFV1_0(w[1], w[27], w[12], pars->GC_11, amp[16]); 
  FFV1_0(w[0], w[27], w[13], pars->GC_11, amp[17]); 
  FFV1_0(w[28], w[7], w[13], pars->GC_11, amp[18]); 
  FFV1_0(w[29], w[7], w[12], pars->GC_11, amp[19]); 
  FFV1_0(w[1], w[30], w[18], pars->GC_11, amp[20]); 
  FFV1_0(w[0], w[30], w[19], pars->GC_11, amp[21]); 
  FFV1_0(w[28], w[6], w[19], pars->GC_11, amp[22]); 
  FFV1_0(w[29], w[6], w[18], pars->GC_11, amp[23]); 
  FFV1_0(w[1], w[30], w[21], pars->GC_11, amp[24]); 
  FFV1_0(w[1], w[27], w[22], pars->GC_11, amp[25]); 
  FFV1_0(w[29], w[7], w[22], pars->GC_11, amp[26]); 
  FFV1_0(w[29], w[6], w[21], pars->GC_11, amp[27]); 
  FFV1_0(w[0], w[30], w[24], pars->GC_11, amp[28]); 
  FFV1_0(w[0], w[27], w[25], pars->GC_11, amp[29]); 
  FFV1_0(w[28], w[7], w[25], pars->GC_11, amp[30]); 
  FFV1_0(w[28], w[6], w[24], pars->GC_11, amp[31]); 
  FFV1_0(w[1], w[11], w[32], pars->GC_11, amp[32]); 
  FFV1_0(w[15], w[7], w[32], pars->GC_11, amp[33]); 
  FFV1_0(w[15], w[16], w[31], pars->GC_11, amp[34]); 
  FFV1_0(w[33], w[16], w[9], pars->GC_2, amp[35]); 
  FFV1_0(w[23], w[11], w[31], pars->GC_11, amp[36]); 
  FFV1_0(w[23], w[34], w[9], pars->GC_2, amp[37]); 
  FFV1_0(w[33], w[11], w[8], pars->GC_11, amp[38]); 
  FFV1_0(w[15], w[34], w[8], pars->GC_11, amp[39]); 
  FFV1_0(w[1], w[27], w[32], pars->GC_11, amp[40]); 
  FFV1_0(w[29], w[7], w[32], pars->GC_11, amp[41]); 
  FFV1_0(w[29], w[16], w[31], pars->GC_11, amp[42]); 
  FFV2_5_0(w[33], w[16], w[26], pars->GC_51, pars->GC_58, amp[43]); 
  FFV1_0(w[23], w[27], w[31], pars->GC_11, amp[44]); 
  FFV2_5_0(w[23], w[34], w[26], pars->GC_51, pars->GC_58, amp[45]); 
  FFV1_0(w[33], w[27], w[8], pars->GC_11, amp[46]); 
  FFV1_0(w[29], w[34], w[8], pars->GC_11, amp[47]); 
  FFV1_0(w[0], w[11], w[36], pars->GC_11, amp[48]); 
  FFV1_0(w[14], w[7], w[36], pars->GC_11, amp[49]); 
  FFV1_0(w[14], w[16], w[35], pars->GC_11, amp[50]); 
  FFV1_0(w[37], w[16], w[9], pars->GC_2, amp[51]); 
  FFV1_0(w[20], w[11], w[35], pars->GC_11, amp[52]); 
  FFV1_0(w[20], w[38], w[9], pars->GC_2, amp[53]); 
  FFV1_0(w[37], w[11], w[8], pars->GC_11, amp[54]); 
  FFV1_0(w[14], w[38], w[8], pars->GC_11, amp[55]); 
  FFV1_0(w[0], w[27], w[36], pars->GC_11, amp[56]); 
  FFV1_0(w[28], w[7], w[36], pars->GC_11, amp[57]); 
  FFV1_0(w[28], w[16], w[35], pars->GC_11, amp[58]); 
  FFV2_5_0(w[37], w[16], w[26], pars->GC_51, pars->GC_58, amp[59]); 
  FFV1_0(w[20], w[27], w[35], pars->GC_11, amp[60]); 
  FFV2_5_0(w[20], w[38], w[26], pars->GC_51, pars->GC_58, amp[61]); 
  FFV1_0(w[37], w[27], w[8], pars->GC_11, amp[62]); 
  FFV1_0(w[28], w[38], w[8], pars->GC_11, amp[63]); 
  FFV1_0(w[15], w[10], w[39], pars->GC_11, amp[64]); 
  FFV1_0(w[40], w[10], w[9], pars->GC_2, amp[65]); 
  FFV1_0(w[1], w[17], w[41], pars->GC_11, amp[66]); 
  FFV1_0(w[15], w[6], w[41], pars->GC_11, amp[67]); 
  FFV1_0(w[23], w[17], w[39], pars->GC_11, amp[68]); 
  FFV1_0(w[23], w[42], w[9], pars->GC_2, amp[69]); 
  FFV1_0(w[40], w[17], w[8], pars->GC_11, amp[70]); 
  FFV1_0(w[15], w[42], w[8], pars->GC_11, amp[71]); 
  FFV1_0(w[29], w[10], w[39], pars->GC_11, amp[72]); 
  FFV2_5_0(w[40], w[10], w[26], pars->GC_51, pars->GC_58, amp[73]); 
  FFV1_0(w[1], w[30], w[41], pars->GC_11, amp[74]); 
  FFV1_0(w[29], w[6], w[41], pars->GC_11, amp[75]); 
  FFV1_0(w[23], w[30], w[39], pars->GC_11, amp[76]); 
  FFV2_5_0(w[23], w[42], w[26], pars->GC_51, pars->GC_58, amp[77]); 
  FFV1_0(w[40], w[30], w[8], pars->GC_11, amp[78]); 
  FFV1_0(w[29], w[42], w[8], pars->GC_11, amp[79]); 
  FFV1_0(w[14], w[10], w[43], pars->GC_11, amp[80]); 
  FFV1_0(w[44], w[10], w[9], pars->GC_2, amp[81]); 
  FFV1_0(w[0], w[17], w[45], pars->GC_11, amp[82]); 
  FFV1_0(w[14], w[6], w[45], pars->GC_11, amp[83]); 
  FFV1_0(w[20], w[17], w[43], pars->GC_11, amp[84]); 
  FFV1_0(w[20], w[46], w[9], pars->GC_2, amp[85]); 
  FFV1_0(w[44], w[17], w[8], pars->GC_11, amp[86]); 
  FFV1_0(w[14], w[46], w[8], pars->GC_11, amp[87]); 
  FFV1_0(w[28], w[10], w[43], pars->GC_11, amp[88]); 
  FFV2_5_0(w[44], w[10], w[26], pars->GC_51, pars->GC_58, amp[89]); 
  FFV1_0(w[0], w[30], w[45], pars->GC_11, amp[90]); 
  FFV1_0(w[28], w[6], w[45], pars->GC_11, amp[91]); 
  FFV1_0(w[20], w[30], w[43], pars->GC_11, amp[92]); 
  FFV2_5_0(w[20], w[46], w[26], pars->GC_51, pars->GC_58, amp[93]); 
  FFV1_0(w[44], w[30], w[8], pars->GC_11, amp[94]); 
  FFV1_0(w[28], w[46], w[8], pars->GC_11, amp[95]); 
  FFV1_0(w[1], w[49], w[50], pars->GC_11, amp[96]); 
  FFV1_0(w[0], w[49], w[51], pars->GC_11, amp[97]); 
  FFV1_0(w[1], w[53], w[52], pars->GC_11, amp[98]); 
  FFV1_0(w[15], w[48], w[52], pars->GC_11, amp[99]); 
  FFV1_0(w[0], w[53], w[54], pars->GC_11, amp[100]); 
  FFV1_0(w[14], w[48], w[54], pars->GC_11, amp[101]); 
  FFV1_0(w[15], w[47], w[50], pars->GC_11, amp[102]); 
  FFV1_0(w[14], w[47], w[51], pars->GC_11, amp[103]); 
  FFV1_0(w[1], w[55], w[50], pars->GC_11, amp[104]); 
  FFV1_0(w[0], w[55], w[51], pars->GC_11, amp[105]); 
  FFV1_0(w[1], w[56], w[52], pars->GC_11, amp[106]); 
  FFV1_0(w[29], w[48], w[52], pars->GC_11, amp[107]); 
  FFV1_0(w[0], w[56], w[54], pars->GC_11, amp[108]); 
  FFV1_0(w[28], w[48], w[54], pars->GC_11, amp[109]); 
  FFV1_0(w[29], w[47], w[50], pars->GC_11, amp[110]); 
  FFV1_0(w[28], w[47], w[51], pars->GC_11, amp[111]); 
  FFV1_0(w[1], w[11], w[58], pars->GC_11, amp[112]); 
  FFV1_0(w[15], w[7], w[58], pars->GC_11, amp[113]); 
  FFV1_0(w[1], w[49], w[59], pars->GC_11, amp[114]); 
  FFV1_0(w[60], w[7], w[54], pars->GC_11, amp[115]); 
  FFV1_0(w[57], w[11], w[54], pars->GC_11, amp[116]); 
  FFV1_0(w[15], w[47], w[59], pars->GC_11, amp[117]); 
  FFV1_0(w[1], w[27], w[58], pars->GC_11, amp[118]); 
  FFV1_0(w[29], w[7], w[58], pars->GC_11, amp[119]); 
  FFV1_0(w[1], w[55], w[59], pars->GC_11, amp[120]); 
  FFV1_0(w[61], w[7], w[54], pars->GC_11, amp[121]); 
  FFV1_0(w[57], w[27], w[54], pars->GC_11, amp[122]); 
  FFV1_0(w[29], w[47], w[59], pars->GC_11, amp[123]); 
  FFV1_0(w[57], w[49], w[43], pars->GC_11, amp[124]); 
  FFV1_0(w[57], w[62], w[9], pars->GC_2, amp[125]); 
  FFV1_0(w[57], w[55], w[43], pars->GC_11, amp[126]); 
  FFV2_5_0(w[57], w[62], w[26], pars->GC_51, pars->GC_58, amp[127]); 
  FFV1_0(w[0], w[11], w[64], pars->GC_11, amp[128]); 
  FFV1_0(w[14], w[7], w[64], pars->GC_11, amp[129]); 
  FFV1_0(w[0], w[49], w[65], pars->GC_11, amp[130]); 
  FFV1_0(w[66], w[7], w[52], pars->GC_11, amp[131]); 
  FFV1_0(w[63], w[11], w[52], pars->GC_11, amp[132]); 
  FFV1_0(w[14], w[47], w[65], pars->GC_11, amp[133]); 
  FFV1_0(w[0], w[27], w[64], pars->GC_11, amp[134]); 
  FFV1_0(w[28], w[7], w[64], pars->GC_11, amp[135]); 
  FFV1_0(w[0], w[55], w[65], pars->GC_11, amp[136]); 
  FFV1_0(w[67], w[7], w[52], pars->GC_11, amp[137]); 
  FFV1_0(w[63], w[27], w[52], pars->GC_11, amp[138]); 
  FFV1_0(w[28], w[47], w[65], pars->GC_11, amp[139]); 
  FFV1_0(w[63], w[49], w[39], pars->GC_11, amp[140]); 
  FFV1_0(w[63], w[68], w[9], pars->GC_2, amp[141]); 
  FFV1_0(w[63], w[55], w[39], pars->GC_11, amp[142]); 
  FFV2_5_0(w[63], w[68], w[26], pars->GC_51, pars->GC_58, amp[143]); 
  FFV1_0(w[1], w[11], w[70], pars->GC_11, amp[144]); 
  FFV1_0(w[0], w[11], w[71], pars->GC_11, amp[145]); 
  FFV1_0(w[14], w[7], w[71], pars->GC_11, amp[146]); 
  FFV1_0(w[15], w[7], w[70], pars->GC_11, amp[147]); 
  FFV1_0(w[1], w[11], w[72], pars->GC_11, amp[148]); 
  FFV1_0(w[73], w[11], w[5], pars->GC_11, amp[149]); 
  FFV1_0(w[15], w[7], w[72], pars->GC_11, amp[150]); 
  FFV1_0(w[15], w[74], w[5], pars->GC_11, amp[151]); 
  FFV1_0(w[0], w[11], w[75], pars->GC_11, amp[152]); 
  FFV1_0(w[76], w[11], w[5], pars->GC_11, amp[153]); 
  FFV1_0(w[14], w[7], w[75], pars->GC_11, amp[154]); 
  FFV1_0(w[14], w[77], w[5], pars->GC_11, amp[155]); 
  FFV1_0(w[1], w[27], w[70], pars->GC_11, amp[156]); 
  FFV1_0(w[0], w[27], w[71], pars->GC_11, amp[157]); 
  FFV1_0(w[28], w[7], w[71], pars->GC_11, amp[158]); 
  FFV1_0(w[29], w[7], w[70], pars->GC_11, amp[159]); 
  FFV1_0(w[1], w[27], w[72], pars->GC_11, amp[160]); 
  FFV1_0(w[73], w[27], w[5], pars->GC_11, amp[161]); 
  FFV1_0(w[29], w[7], w[72], pars->GC_11, amp[162]); 
  FFV1_0(w[29], w[74], w[5], pars->GC_11, amp[163]); 
  FFV1_0(w[0], w[27], w[75], pars->GC_11, amp[164]); 
  FFV1_0(w[76], w[27], w[5], pars->GC_11, amp[165]); 
  FFV1_0(w[28], w[7], w[75], pars->GC_11, amp[166]); 
  FFV1_0(w[28], w[77], w[5], pars->GC_11, amp[167]); 
  FFV1_0(w[15], w[69], w[39], pars->GC_11, amp[168]); 
  FFV1_0(w[40], w[69], w[9], pars->GC_2, amp[169]); 
  FFV1_0(w[1], w[49], w[78], pars->GC_11, amp[170]); 
  FFV1_0(w[40], w[49], w[5], pars->GC_11, amp[171]); 
  FFV1_0(w[15], w[68], w[5], pars->GC_11, amp[172]); 
  FFV1_0(w[15], w[47], w[78], pars->GC_11, amp[173]); 
  FFV1_0(w[29], w[69], w[39], pars->GC_11, amp[174]); 
  FFV2_5_0(w[40], w[69], w[26], pars->GC_51, pars->GC_58, amp[175]); 
  FFV1_0(w[1], w[55], w[78], pars->GC_11, amp[176]); 
  FFV1_0(w[40], w[55], w[5], pars->GC_11, amp[177]); 
  FFV1_0(w[29], w[68], w[5], pars->GC_11, amp[178]); 
  FFV1_0(w[29], w[47], w[78], pars->GC_11, amp[179]); 
  FFV1_0(w[14], w[69], w[43], pars->GC_11, amp[180]); 
  FFV1_0(w[44], w[69], w[9], pars->GC_2, amp[181]); 
  FFV1_0(w[0], w[49], w[79], pars->GC_11, amp[182]); 
  FFV1_0(w[44], w[49], w[5], pars->GC_11, amp[183]); 
  FFV1_0(w[14], w[62], w[5], pars->GC_11, amp[184]); 
  FFV1_0(w[14], w[47], w[79], pars->GC_11, amp[185]); 
  FFV1_0(w[28], w[69], w[43], pars->GC_11, amp[186]); 
  FFV2_5_0(w[44], w[69], w[26], pars->GC_51, pars->GC_58, amp[187]); 
  FFV1_0(w[0], w[55], w[79], pars->GC_11, amp[188]); 
  FFV1_0(w[44], w[55], w[5], pars->GC_11, amp[189]); 
  FFV1_0(w[28], w[62], w[5], pars->GC_11, amp[190]); 
  FFV1_0(w[28], w[47], w[79], pars->GC_11, amp[191]); 
  FFV1_0(w[1], w[82], w[83], pars->GC_11, amp[192]); 
  FFV1_0(w[0], w[82], w[84], pars->GC_11, amp[193]); 
  FFV1_0(w[1], w[86], w[85], pars->GC_11, amp[194]); 
  FFV1_0(w[15], w[81], w[85], pars->GC_11, amp[195]); 
  FFV1_0(w[0], w[86], w[87], pars->GC_11, amp[196]); 
  FFV1_0(w[14], w[81], w[87], pars->GC_11, amp[197]); 
  FFV1_0(w[15], w[80], w[83], pars->GC_11, amp[198]); 
  FFV1_0(w[14], w[80], w[84], pars->GC_11, amp[199]); 
  FFV1_0(w[1], w[88], w[83], pars->GC_11, amp[200]); 
  FFV1_0(w[0], w[88], w[84], pars->GC_11, amp[201]); 
  FFV1_0(w[1], w[89], w[85], pars->GC_11, amp[202]); 
  FFV1_0(w[29], w[81], w[85], pars->GC_11, amp[203]); 
  FFV1_0(w[0], w[89], w[87], pars->GC_11, amp[204]); 
  FFV1_0(w[28], w[81], w[87], pars->GC_11, amp[205]); 
  FFV1_0(w[29], w[80], w[83], pars->GC_11, amp[206]); 
  FFV1_0(w[28], w[80], w[84], pars->GC_11, amp[207]); 
  FFV1_0(w[1], w[17], w[90], pars->GC_11, amp[208]); 
  FFV1_0(w[15], w[6], w[90], pars->GC_11, amp[209]); 
  FFV1_0(w[1], w[82], w[91], pars->GC_11, amp[210]); 
  FFV1_0(w[60], w[6], w[87], pars->GC_11, amp[211]); 
  FFV1_0(w[57], w[17], w[87], pars->GC_11, amp[212]); 
  FFV1_0(w[15], w[80], w[91], pars->GC_11, amp[213]); 
  FFV1_0(w[1], w[30], w[90], pars->GC_11, amp[214]); 
  FFV1_0(w[29], w[6], w[90], pars->GC_11, amp[215]); 
  FFV1_0(w[1], w[88], w[91], pars->GC_11, amp[216]); 
  FFV1_0(w[61], w[6], w[87], pars->GC_11, amp[217]); 
  FFV1_0(w[57], w[30], w[87], pars->GC_11, amp[218]); 
  FFV1_0(w[29], w[80], w[91], pars->GC_11, amp[219]); 
  FFV1_0(w[57], w[82], w[35], pars->GC_11, amp[220]); 
  FFV1_0(w[57], w[92], w[9], pars->GC_2, amp[221]); 
  FFV1_0(w[57], w[88], w[35], pars->GC_11, amp[222]); 
  FFV2_5_0(w[57], w[92], w[26], pars->GC_51, pars->GC_58, amp[223]); 
  FFV1_0(w[0], w[17], w[93], pars->GC_11, amp[224]); 
  FFV1_0(w[14], w[6], w[93], pars->GC_11, amp[225]); 
  FFV1_0(w[0], w[82], w[94], pars->GC_11, amp[226]); 
  FFV1_0(w[66], w[6], w[85], pars->GC_11, amp[227]); 
  FFV1_0(w[63], w[17], w[85], pars->GC_11, amp[228]); 
  FFV1_0(w[14], w[80], w[94], pars->GC_11, amp[229]); 
  FFV1_0(w[0], w[30], w[93], pars->GC_11, amp[230]); 
  FFV1_0(w[28], w[6], w[93], pars->GC_11, amp[231]); 
  FFV1_0(w[0], w[88], w[94], pars->GC_11, amp[232]); 
  FFV1_0(w[67], w[6], w[85], pars->GC_11, amp[233]); 
  FFV1_0(w[63], w[30], w[85], pars->GC_11, amp[234]); 
  FFV1_0(w[28], w[80], w[94], pars->GC_11, amp[235]); 
  FFV1_0(w[63], w[82], w[31], pars->GC_11, amp[236]); 
  FFV1_0(w[63], w[95], w[9], pars->GC_2, amp[237]); 
  FFV1_0(w[63], w[88], w[31], pars->GC_11, amp[238]); 
  FFV2_5_0(w[63], w[95], w[26], pars->GC_51, pars->GC_58, amp[239]); 
  FFV1_0(w[1], w[17], w[97], pars->GC_11, amp[240]); 
  FFV1_0(w[0], w[17], w[98], pars->GC_11, amp[241]); 
  FFV1_0(w[14], w[6], w[98], pars->GC_11, amp[242]); 
  FFV1_0(w[15], w[6], w[97], pars->GC_11, amp[243]); 
  FFV1_0(w[1], w[17], w[99], pars->GC_11, amp[244]); 
  FFV1_0(w[100], w[17], w[5], pars->GC_11, amp[245]); 
  FFV1_0(w[15], w[6], w[99], pars->GC_11, amp[246]); 
  FFV1_0(w[15], w[101], w[5], pars->GC_11, amp[247]); 
  FFV1_0(w[0], w[17], w[102], pars->GC_11, amp[248]); 
  FFV1_0(w[103], w[17], w[5], pars->GC_11, amp[249]); 
  FFV1_0(w[14], w[6], w[102], pars->GC_11, amp[250]); 
  FFV1_0(w[14], w[104], w[5], pars->GC_11, amp[251]); 
  FFV1_0(w[1], w[30], w[97], pars->GC_11, amp[252]); 
  FFV1_0(w[0], w[30], w[98], pars->GC_11, amp[253]); 
  FFV1_0(w[28], w[6], w[98], pars->GC_11, amp[254]); 
  FFV1_0(w[29], w[6], w[97], pars->GC_11, amp[255]); 
  FFV1_0(w[1], w[30], w[99], pars->GC_11, amp[256]); 
  FFV1_0(w[100], w[30], w[5], pars->GC_11, amp[257]); 
  FFV1_0(w[29], w[6], w[99], pars->GC_11, amp[258]); 
  FFV1_0(w[29], w[101], w[5], pars->GC_11, amp[259]); 
  FFV1_0(w[0], w[30], w[102], pars->GC_11, amp[260]); 
  FFV1_0(w[103], w[30], w[5], pars->GC_11, amp[261]); 
  FFV1_0(w[28], w[6], w[102], pars->GC_11, amp[262]); 
  FFV1_0(w[28], w[104], w[5], pars->GC_11, amp[263]); 
  FFV1_0(w[15], w[96], w[31], pars->GC_11, amp[264]); 
  FFV1_0(w[33], w[96], w[9], pars->GC_2, amp[265]); 
  FFV1_0(w[1], w[82], w[105], pars->GC_11, amp[266]); 
  FFV1_0(w[33], w[82], w[5], pars->GC_11, amp[267]); 
  FFV1_0(w[15], w[95], w[5], pars->GC_11, amp[268]); 
  FFV1_0(w[15], w[80], w[105], pars->GC_11, amp[269]); 
  FFV1_0(w[29], w[96], w[31], pars->GC_11, amp[270]); 
  FFV2_5_0(w[33], w[96], w[26], pars->GC_51, pars->GC_58, amp[271]); 
  FFV1_0(w[1], w[88], w[105], pars->GC_11, amp[272]); 
  FFV1_0(w[33], w[88], w[5], pars->GC_11, amp[273]); 
  FFV1_0(w[29], w[95], w[5], pars->GC_11, amp[274]); 
  FFV1_0(w[29], w[80], w[105], pars->GC_11, amp[275]); 
  FFV1_0(w[14], w[96], w[35], pars->GC_11, amp[276]); 
  FFV1_0(w[37], w[96], w[9], pars->GC_2, amp[277]); 
  FFV1_0(w[0], w[82], w[106], pars->GC_11, amp[278]); 
  FFV1_0(w[37], w[82], w[5], pars->GC_11, amp[279]); 
  FFV1_0(w[14], w[92], w[5], pars->GC_11, amp[280]); 
  FFV1_0(w[14], w[80], w[106], pars->GC_11, amp[281]); 
  FFV1_0(w[28], w[96], w[35], pars->GC_11, amp[282]); 
  FFV2_5_0(w[37], w[96], w[26], pars->GC_51, pars->GC_58, amp[283]); 
  FFV1_0(w[0], w[88], w[106], pars->GC_11, amp[284]); 
  FFV1_0(w[37], w[88], w[5], pars->GC_11, amp[285]); 
  FFV1_0(w[28], w[92], w[5], pars->GC_11, amp[286]); 
  FFV1_0(w[28], w[80], w[106], pars->GC_11, amp[287]); 
  FFV1_0(w[1], w[11], w[108], pars->GC_11, amp[288]); 
  FFV1_0(w[15], w[7], w[108], pars->GC_11, amp[289]); 
  FFV1_0(w[109], w[7], w[84], pars->GC_11, amp[290]); 
  FFV1_0(w[1], w[86], w[110], pars->GC_11, amp[291]); 
  FFV1_0(w[15], w[81], w[110], pars->GC_11, amp[292]); 
  FFV1_0(w[107], w[11], w[84], pars->GC_11, amp[293]); 
  FFV1_0(w[1], w[27], w[108], pars->GC_11, amp[294]); 
  FFV1_0(w[29], w[7], w[108], pars->GC_11, amp[295]); 
  FFV1_0(w[111], w[7], w[84], pars->GC_11, amp[296]); 
  FFV1_0(w[1], w[89], w[110], pars->GC_11, amp[297]); 
  FFV1_0(w[29], w[81], w[110], pars->GC_11, amp[298]); 
  FFV1_0(w[107], w[27], w[84], pars->GC_11, amp[299]); 
  FFV1_0(w[109], w[81], w[43], pars->GC_11, amp[300]); 
  FFV1_0(w[112], w[81], w[9], pars->GC_2, amp[301]); 
  FFV1_0(w[111], w[81], w[43], pars->GC_11, amp[302]); 
  FFV2_5_0(w[112], w[81], w[26], pars->GC_51, pars->GC_58, amp[303]); 
  FFV1_0(w[1], w[17], w[113], pars->GC_11, amp[304]); 
  FFV1_0(w[15], w[6], w[113], pars->GC_11, amp[305]); 
  FFV1_0(w[109], w[6], w[51], pars->GC_11, amp[306]); 
  FFV1_0(w[1], w[53], w[114], pars->GC_11, amp[307]); 
  FFV1_0(w[15], w[48], w[114], pars->GC_11, amp[308]); 
  FFV1_0(w[107], w[17], w[51], pars->GC_11, amp[309]); 
  FFV1_0(w[1], w[30], w[113], pars->GC_11, amp[310]); 
  FFV1_0(w[29], w[6], w[113], pars->GC_11, amp[311]); 
  FFV1_0(w[111], w[6], w[51], pars->GC_11, amp[312]); 
  FFV1_0(w[1], w[56], w[114], pars->GC_11, amp[313]); 
  FFV1_0(w[29], w[48], w[114], pars->GC_11, amp[314]); 
  FFV1_0(w[107], w[30], w[51], pars->GC_11, amp[315]); 
  FFV1_0(w[109], w[48], w[35], pars->GC_11, amp[316]); 
  FFV1_0(w[115], w[48], w[9], pars->GC_2, amp[317]); 
  FFV1_0(w[111], w[48], w[35], pars->GC_11, amp[318]); 
  FFV2_5_0(w[115], w[48], w[26], pars->GC_51, pars->GC_58, amp[319]); 
  FFV1_0(w[109], w[7], w[94], pars->GC_11, amp[320]); 
  FFV1_0(w[109], w[6], w[65], pars->GC_11, amp[321]); 
  FFV1_0(w[66], w[7], w[114], pars->GC_11, amp[322]); 
  FFV1_0(w[63], w[11], w[114], pars->GC_11, amp[323]); 
  FFV1_0(w[66], w[6], w[110], pars->GC_11, amp[324]); 
  FFV1_0(w[63], w[17], w[110], pars->GC_11, amp[325]); 
  FFV1_0(w[107], w[11], w[94], pars->GC_11, amp[326]); 
  FFV1_0(w[107], w[17], w[65], pars->GC_11, amp[327]); 
  FFV1_0(w[111], w[7], w[94], pars->GC_11, amp[328]); 
  FFV1_0(w[111], w[6], w[65], pars->GC_11, amp[329]); 
  FFV1_0(w[67], w[7], w[114], pars->GC_11, amp[330]); 
  FFV1_0(w[63], w[27], w[114], pars->GC_11, amp[331]); 
  FFV1_0(w[67], w[6], w[110], pars->GC_11, amp[332]); 
  FFV1_0(w[63], w[30], w[110], pars->GC_11, amp[333]); 
  FFV1_0(w[107], w[27], w[94], pars->GC_11, amp[334]); 
  FFV1_0(w[107], w[30], w[65], pars->GC_11, amp[335]); 
  FFV1_0(w[1], w[17], w[117], pars->GC_11, amp[336]); 
  FFV1_0(w[1], w[11], w[118], pars->GC_11, amp[337]); 
  FFV1_0(w[15], w[7], w[118], pars->GC_11, amp[338]); 
  FFV1_0(w[15], w[6], w[117], pars->GC_11, amp[339]); 
  FFV1_0(w[1], w[11], w[119], pars->GC_11, amp[340]); 
  FFV1_0(w[120], w[11], w[5], pars->GC_11, amp[341]); 
  FFV1_0(w[15], w[7], w[119], pars->GC_11, amp[342]); 
  FFV1_0(w[15], w[121], w[5], pars->GC_11, amp[343]); 
  FFV1_0(w[1], w[17], w[122], pars->GC_11, amp[344]); 
  FFV1_0(w[123], w[17], w[5], pars->GC_11, amp[345]); 
  FFV1_0(w[15], w[6], w[122], pars->GC_11, amp[346]); 
  FFV1_0(w[15], w[124], w[5], pars->GC_11, amp[347]); 
  FFV1_0(w[1], w[30], w[117], pars->GC_11, amp[348]); 
  FFV1_0(w[1], w[27], w[118], pars->GC_11, amp[349]); 
  FFV1_0(w[29], w[7], w[118], pars->GC_11, amp[350]); 
  FFV1_0(w[29], w[6], w[117], pars->GC_11, amp[351]); 
  FFV1_0(w[1], w[27], w[119], pars->GC_11, amp[352]); 
  FFV1_0(w[120], w[27], w[5], pars->GC_11, amp[353]); 
  FFV1_0(w[29], w[7], w[119], pars->GC_11, amp[354]); 
  FFV1_0(w[29], w[121], w[5], pars->GC_11, amp[355]); 
  FFV1_0(w[1], w[30], w[122], pars->GC_11, amp[356]); 
  FFV1_0(w[123], w[30], w[5], pars->GC_11, amp[357]); 
  FFV1_0(w[29], w[6], w[122], pars->GC_11, amp[358]); 
  FFV1_0(w[29], w[124], w[5], pars->GC_11, amp[359]); 
  FFV1_0(w[116], w[11], w[35], pars->GC_11, amp[360]); 
  FFV1_0(w[116], w[38], w[9], pars->GC_2, amp[361]); 
  FFV1_0(w[109], w[7], w[106], pars->GC_11, amp[362]); 
  FFV1_0(w[109], w[38], w[5], pars->GC_11, amp[363]); 
  FFV1_0(w[115], w[11], w[5], pars->GC_11, amp[364]); 
  FFV1_0(w[107], w[11], w[106], pars->GC_11, amp[365]); 
  FFV1_0(w[116], w[27], w[35], pars->GC_11, amp[366]); 
  FFV2_5_0(w[116], w[38], w[26], pars->GC_51, pars->GC_58, amp[367]); 
  FFV1_0(w[111], w[7], w[106], pars->GC_11, amp[368]); 
  FFV1_0(w[111], w[38], w[5], pars->GC_11, amp[369]); 
  FFV1_0(w[115], w[27], w[5], pars->GC_11, amp[370]); 
  FFV1_0(w[107], w[27], w[106], pars->GC_11, amp[371]); 
  FFV1_0(w[116], w[17], w[43], pars->GC_11, amp[372]); 
  FFV1_0(w[116], w[46], w[9], pars->GC_2, amp[373]); 
  FFV1_0(w[109], w[6], w[79], pars->GC_11, amp[374]); 
  FFV1_0(w[109], w[46], w[5], pars->GC_11, amp[375]); 
  FFV1_0(w[112], w[17], w[5], pars->GC_11, amp[376]); 
  FFV1_0(w[107], w[17], w[79], pars->GC_11, amp[377]); 
  FFV1_0(w[116], w[30], w[43], pars->GC_11, amp[378]); 
  FFV2_5_0(w[116], w[46], w[26], pars->GC_51, pars->GC_58, amp[379]); 
  FFV1_0(w[111], w[6], w[79], pars->GC_11, amp[380]); 
  FFV1_0(w[111], w[46], w[5], pars->GC_11, amp[381]); 
  FFV1_0(w[112], w[30], w[5], pars->GC_11, amp[382]); 
  FFV1_0(w[107], w[30], w[79], pars->GC_11, amp[383]); 
  FFV1_0(w[0], w[11], w[126], pars->GC_11, amp[384]); 
  FFV1_0(w[14], w[7], w[126], pars->GC_11, amp[385]); 
  FFV1_0(w[127], w[7], w[83], pars->GC_11, amp[386]); 
  FFV1_0(w[0], w[86], w[128], pars->GC_11, amp[387]); 
  FFV1_0(w[14], w[81], w[128], pars->GC_11, amp[388]); 
  FFV1_0(w[125], w[11], w[83], pars->GC_11, amp[389]); 
  FFV1_0(w[0], w[27], w[126], pars->GC_11, amp[390]); 
  FFV1_0(w[28], w[7], w[126], pars->GC_11, amp[391]); 
  FFV1_0(w[129], w[7], w[83], pars->GC_11, amp[392]); 
  FFV1_0(w[0], w[89], w[128], pars->GC_11, amp[393]); 
  FFV1_0(w[28], w[81], w[128], pars->GC_11, amp[394]); 
  FFV1_0(w[125], w[27], w[83], pars->GC_11, amp[395]); 
  FFV1_0(w[127], w[81], w[39], pars->GC_11, amp[396]); 
  FFV1_0(w[130], w[81], w[9], pars->GC_2, amp[397]); 
  FFV1_0(w[129], w[81], w[39], pars->GC_11, amp[398]); 
  FFV2_5_0(w[130], w[81], w[26], pars->GC_51, pars->GC_58, amp[399]); 
  FFV1_0(w[0], w[17], w[131], pars->GC_11, amp[400]); 
  FFV1_0(w[14], w[6], w[131], pars->GC_11, amp[401]); 
  FFV1_0(w[127], w[6], w[50], pars->GC_11, amp[402]); 
  FFV1_0(w[0], w[53], w[132], pars->GC_11, amp[403]); 
  FFV1_0(w[14], w[48], w[132], pars->GC_11, amp[404]); 
  FFV1_0(w[125], w[17], w[50], pars->GC_11, amp[405]); 
  FFV1_0(w[0], w[30], w[131], pars->GC_11, amp[406]); 
  FFV1_0(w[28], w[6], w[131], pars->GC_11, amp[407]); 
  FFV1_0(w[129], w[6], w[50], pars->GC_11, amp[408]); 
  FFV1_0(w[0], w[56], w[132], pars->GC_11, amp[409]); 
  FFV1_0(w[28], w[48], w[132], pars->GC_11, amp[410]); 
  FFV1_0(w[125], w[30], w[50], pars->GC_11, amp[411]); 
  FFV1_0(w[127], w[48], w[31], pars->GC_11, amp[412]); 
  FFV1_0(w[133], w[48], w[9], pars->GC_2, amp[413]); 
  FFV1_0(w[129], w[48], w[31], pars->GC_11, amp[414]); 
  FFV2_5_0(w[133], w[48], w[26], pars->GC_51, pars->GC_58, amp[415]); 
  FFV1_0(w[127], w[7], w[91], pars->GC_11, amp[416]); 
  FFV1_0(w[127], w[6], w[59], pars->GC_11, amp[417]); 
  FFV1_0(w[60], w[7], w[132], pars->GC_11, amp[418]); 
  FFV1_0(w[57], w[11], w[132], pars->GC_11, amp[419]); 
  FFV1_0(w[60], w[6], w[128], pars->GC_11, amp[420]); 
  FFV1_0(w[57], w[17], w[128], pars->GC_11, amp[421]); 
  FFV1_0(w[125], w[11], w[91], pars->GC_11, amp[422]); 
  FFV1_0(w[125], w[17], w[59], pars->GC_11, amp[423]); 
  FFV1_0(w[129], w[7], w[91], pars->GC_11, amp[424]); 
  FFV1_0(w[129], w[6], w[59], pars->GC_11, amp[425]); 
  FFV1_0(w[61], w[7], w[132], pars->GC_11, amp[426]); 
  FFV1_0(w[57], w[27], w[132], pars->GC_11, amp[427]); 
  FFV1_0(w[61], w[6], w[128], pars->GC_11, amp[428]); 
  FFV1_0(w[57], w[30], w[128], pars->GC_11, amp[429]); 
  FFV1_0(w[125], w[27], w[91], pars->GC_11, amp[430]); 
  FFV1_0(w[125], w[30], w[59], pars->GC_11, amp[431]); 
  FFV1_0(w[0], w[17], w[135], pars->GC_11, amp[432]); 
  FFV1_0(w[0], w[11], w[136], pars->GC_11, amp[433]); 
  FFV1_0(w[14], w[7], w[136], pars->GC_11, amp[434]); 
  FFV1_0(w[14], w[6], w[135], pars->GC_11, amp[435]); 
  FFV1_0(w[0], w[11], w[137], pars->GC_11, amp[436]); 
  FFV1_0(w[138], w[11], w[5], pars->GC_11, amp[437]); 
  FFV1_0(w[14], w[7], w[137], pars->GC_11, amp[438]); 
  FFV1_0(w[14], w[139], w[5], pars->GC_11, amp[439]); 
  FFV1_0(w[0], w[17], w[140], pars->GC_11, amp[440]); 
  FFV1_0(w[141], w[17], w[5], pars->GC_11, amp[441]); 
  FFV1_0(w[14], w[6], w[140], pars->GC_11, amp[442]); 
  FFV1_0(w[14], w[142], w[5], pars->GC_11, amp[443]); 
  FFV1_0(w[0], w[30], w[135], pars->GC_11, amp[444]); 
  FFV1_0(w[0], w[27], w[136], pars->GC_11, amp[445]); 
  FFV1_0(w[28], w[7], w[136], pars->GC_11, amp[446]); 
  FFV1_0(w[28], w[6], w[135], pars->GC_11, amp[447]); 
  FFV1_0(w[0], w[27], w[137], pars->GC_11, amp[448]); 
  FFV1_0(w[138], w[27], w[5], pars->GC_11, amp[449]); 
  FFV1_0(w[28], w[7], w[137], pars->GC_11, amp[450]); 
  FFV1_0(w[28], w[139], w[5], pars->GC_11, amp[451]); 
  FFV1_0(w[0], w[30], w[140], pars->GC_11, amp[452]); 
  FFV1_0(w[141], w[30], w[5], pars->GC_11, amp[453]); 
  FFV1_0(w[28], w[6], w[140], pars->GC_11, amp[454]); 
  FFV1_0(w[28], w[142], w[5], pars->GC_11, amp[455]); 
  FFV1_0(w[134], w[11], w[31], pars->GC_11, amp[456]); 
  FFV1_0(w[134], w[34], w[9], pars->GC_2, amp[457]); 
  FFV1_0(w[127], w[7], w[105], pars->GC_11, amp[458]); 
  FFV1_0(w[127], w[34], w[5], pars->GC_11, amp[459]); 
  FFV1_0(w[133], w[11], w[5], pars->GC_11, amp[460]); 
  FFV1_0(w[125], w[11], w[105], pars->GC_11, amp[461]); 
  FFV1_0(w[134], w[27], w[31], pars->GC_11, amp[462]); 
  FFV2_5_0(w[134], w[34], w[26], pars->GC_51, pars->GC_58, amp[463]); 
  FFV1_0(w[129], w[7], w[105], pars->GC_11, amp[464]); 
  FFV1_0(w[129], w[34], w[5], pars->GC_11, amp[465]); 
  FFV1_0(w[133], w[27], w[5], pars->GC_11, amp[466]); 
  FFV1_0(w[125], w[27], w[105], pars->GC_11, amp[467]); 
  FFV1_0(w[134], w[17], w[39], pars->GC_11, amp[468]); 
  FFV1_0(w[134], w[42], w[9], pars->GC_2, amp[469]); 
  FFV1_0(w[127], w[6], w[78], pars->GC_11, amp[470]); 
  FFV1_0(w[127], w[42], w[5], pars->GC_11, amp[471]); 
  FFV1_0(w[130], w[17], w[5], pars->GC_11, amp[472]); 
  FFV1_0(w[125], w[17], w[78], pars->GC_11, amp[473]); 
  FFV1_0(w[134], w[30], w[39], pars->GC_11, amp[474]); 
  FFV2_5_0(w[134], w[42], w[26], pars->GC_51, pars->GC_58, amp[475]); 
  FFV1_0(w[129], w[6], w[78], pars->GC_11, amp[476]); 
  FFV1_0(w[129], w[42], w[5], pars->GC_11, amp[477]); 
  FFV1_0(w[130], w[30], w[5], pars->GC_11, amp[478]); 
  FFV1_0(w[125], w[30], w[78], pars->GC_11, amp[479]); 
  FFV1_0(w[1], w[11], w[144], pars->GC_11, amp[480]); 
  FFV1_0(w[0], w[11], w[145], pars->GC_11, amp[481]); 
  FFV1_0(w[14], w[7], w[145], pars->GC_11, amp[482]); 
  FFV1_0(w[15], w[7], w[144], pars->GC_11, amp[483]); 
  FFV1_0(w[1], w[11], w[146], pars->GC_11, amp[484]); 
  FFV1_0(w[1], w[147], w[83], pars->GC_11, amp[485]); 
  FFV1_0(w[15], w[7], w[146], pars->GC_11, amp[486]); 
  FFV1_0(w[148], w[7], w[83], pars->GC_11, amp[487]); 
  FFV1_0(w[0], w[11], w[149], pars->GC_11, amp[488]); 
  FFV1_0(w[0], w[147], w[84], pars->GC_11, amp[489]); 
  FFV1_0(w[14], w[7], w[149], pars->GC_11, amp[490]); 
  FFV1_0(w[150], w[7], w[84], pars->GC_11, amp[491]); 
  FFV1_0(w[1], w[27], w[144], pars->GC_11, amp[492]); 
  FFV1_0(w[0], w[27], w[145], pars->GC_11, amp[493]); 
  FFV1_0(w[28], w[7], w[145], pars->GC_11, amp[494]); 
  FFV1_0(w[29], w[7], w[144], pars->GC_11, amp[495]); 
  FFV1_0(w[1], w[27], w[146], pars->GC_11, amp[496]); 
  FFV1_0(w[1], w[151], w[83], pars->GC_11, amp[497]); 
  FFV1_0(w[29], w[7], w[146], pars->GC_11, amp[498]); 
  FFV1_0(w[152], w[7], w[83], pars->GC_11, amp[499]); 
  FFV1_0(w[0], w[27], w[149], pars->GC_11, amp[500]); 
  FFV1_0(w[0], w[151], w[84], pars->GC_11, amp[501]); 
  FFV1_0(w[28], w[7], w[149], pars->GC_11, amp[502]); 
  FFV1_0(w[153], w[7], w[84], pars->GC_11, amp[503]); 
  FFV1_0(w[15], w[143], w[39], pars->GC_11, amp[504]); 
  FFV1_0(w[40], w[143], w[9], pars->GC_2, amp[505]); 
  FFV1_0(w[1], w[86], w[154], pars->GC_11, amp[506]); 
  FFV1_0(w[15], w[81], w[154], pars->GC_11, amp[507]); 
  FFV1_0(w[40], w[86], w[4], pars->GC_11, amp[508]); 
  FFV1_0(w[15], w[155], w[4], pars->GC_11, amp[509]); 
  FFV1_0(w[29], w[143], w[39], pars->GC_11, amp[510]); 
  FFV2_5_0(w[40], w[143], w[26], pars->GC_51, pars->GC_58, amp[511]); 
  FFV1_0(w[1], w[89], w[154], pars->GC_11, amp[512]); 
  FFV1_0(w[29], w[81], w[154], pars->GC_11, amp[513]); 
  FFV1_0(w[40], w[89], w[4], pars->GC_11, amp[514]); 
  FFV1_0(w[29], w[155], w[4], pars->GC_11, amp[515]); 
  FFV1_0(w[14], w[143], w[43], pars->GC_11, amp[516]); 
  FFV1_0(w[44], w[143], w[9], pars->GC_2, amp[517]); 
  FFV1_0(w[0], w[86], w[156], pars->GC_11, amp[518]); 
  FFV1_0(w[14], w[81], w[156], pars->GC_11, amp[519]); 
  FFV1_0(w[44], w[86], w[4], pars->GC_11, amp[520]); 
  FFV1_0(w[14], w[157], w[4], pars->GC_11, amp[521]); 
  FFV1_0(w[28], w[143], w[43], pars->GC_11, amp[522]); 
  FFV2_5_0(w[44], w[143], w[26], pars->GC_51, pars->GC_58, amp[523]); 
  FFV1_0(w[0], w[89], w[156], pars->GC_11, amp[524]); 
  FFV1_0(w[28], w[81], w[156], pars->GC_11, amp[525]); 
  FFV1_0(w[44], w[89], w[4], pars->GC_11, amp[526]); 
  FFV1_0(w[28], w[157], w[4], pars->GC_11, amp[527]); 
  FFV1_0(w[1], w[17], w[159], pars->GC_11, amp[528]); 
  FFV1_0(w[0], w[17], w[160], pars->GC_11, amp[529]); 
  FFV1_0(w[14], w[6], w[160], pars->GC_11, amp[530]); 
  FFV1_0(w[15], w[6], w[159], pars->GC_11, amp[531]); 
  FFV1_0(w[1], w[17], w[161], pars->GC_11, amp[532]); 
  FFV1_0(w[1], w[162], w[50], pars->GC_11, amp[533]); 
  FFV1_0(w[15], w[6], w[161], pars->GC_11, amp[534]); 
  FFV1_0(w[148], w[6], w[50], pars->GC_11, amp[535]); 
  FFV1_0(w[0], w[17], w[163], pars->GC_11, amp[536]); 
  FFV1_0(w[0], w[162], w[51], pars->GC_11, amp[537]); 
  FFV1_0(w[14], w[6], w[163], pars->GC_11, amp[538]); 
  FFV1_0(w[150], w[6], w[51], pars->GC_11, amp[539]); 
  FFV1_0(w[1], w[30], w[159], pars->GC_11, amp[540]); 
  FFV1_0(w[0], w[30], w[160], pars->GC_11, amp[541]); 
  FFV1_0(w[28], w[6], w[160], pars->GC_11, amp[542]); 
  FFV1_0(w[29], w[6], w[159], pars->GC_11, amp[543]); 
  FFV1_0(w[1], w[30], w[161], pars->GC_11, amp[544]); 
  FFV1_0(w[1], w[164], w[50], pars->GC_11, amp[545]); 
  FFV1_0(w[29], w[6], w[161], pars->GC_11, amp[546]); 
  FFV1_0(w[152], w[6], w[50], pars->GC_11, amp[547]); 
  FFV1_0(w[0], w[30], w[163], pars->GC_11, amp[548]); 
  FFV1_0(w[0], w[164], w[51], pars->GC_11, amp[549]); 
  FFV1_0(w[28], w[6], w[163], pars->GC_11, amp[550]); 
  FFV1_0(w[153], w[6], w[51], pars->GC_11, amp[551]); 
  FFV1_0(w[15], w[158], w[31], pars->GC_11, amp[552]); 
  FFV1_0(w[33], w[158], w[9], pars->GC_2, amp[553]); 
  FFV1_0(w[1], w[53], w[165], pars->GC_11, amp[554]); 
  FFV1_0(w[15], w[48], w[165], pars->GC_11, amp[555]); 
  FFV1_0(w[33], w[53], w[4], pars->GC_11, amp[556]); 
  FFV1_0(w[15], w[166], w[4], pars->GC_11, amp[557]); 
  FFV1_0(w[29], w[158], w[31], pars->GC_11, amp[558]); 
  FFV2_5_0(w[33], w[158], w[26], pars->GC_51, pars->GC_58, amp[559]); 
  FFV1_0(w[1], w[56], w[165], pars->GC_11, amp[560]); 
  FFV1_0(w[29], w[48], w[165], pars->GC_11, amp[561]); 
  FFV1_0(w[33], w[56], w[4], pars->GC_11, amp[562]); 
  FFV1_0(w[29], w[166], w[4], pars->GC_11, amp[563]); 
  FFV1_0(w[14], w[158], w[35], pars->GC_11, amp[564]); 
  FFV1_0(w[37], w[158], w[9], pars->GC_2, amp[565]); 
  FFV1_0(w[0], w[53], w[167], pars->GC_11, amp[566]); 
  FFV1_0(w[14], w[48], w[167], pars->GC_11, amp[567]); 
  FFV1_0(w[37], w[53], w[4], pars->GC_11, amp[568]); 
  FFV1_0(w[14], w[168], w[4], pars->GC_11, amp[569]); 
  FFV1_0(w[28], w[158], w[35], pars->GC_11, amp[570]); 
  FFV2_5_0(w[37], w[158], w[26], pars->GC_51, pars->GC_58, amp[571]); 
  FFV1_0(w[0], w[56], w[167], pars->GC_11, amp[572]); 
  FFV1_0(w[28], w[48], w[167], pars->GC_11, amp[573]); 
  FFV1_0(w[37], w[56], w[4], pars->GC_11, amp[574]); 
  FFV1_0(w[28], w[168], w[4], pars->GC_11, amp[575]); 
  FFV1_0(w[1], w[17], w[170], pars->GC_11, amp[576]); 
  FFV1_0(w[1], w[11], w[171], pars->GC_11, amp[577]); 
  FFV1_0(w[15], w[7], w[171], pars->GC_11, amp[578]); 
  FFV1_0(w[15], w[6], w[170], pars->GC_11, amp[579]); 
  FFV1_0(w[1], w[11], w[172], pars->GC_11, amp[580]); 
  FFV1_0(w[1], w[147], w[91], pars->GC_11, amp[581]); 
  FFV1_0(w[15], w[7], w[172], pars->GC_11, amp[582]); 
  FFV1_0(w[148], w[7], w[91], pars->GC_11, amp[583]); 
  FFV1_0(w[1], w[17], w[173], pars->GC_11, amp[584]); 
  FFV1_0(w[1], w[162], w[59], pars->GC_11, amp[585]); 
  FFV1_0(w[15], w[6], w[173], pars->GC_11, amp[586]); 
  FFV1_0(w[148], w[6], w[59], pars->GC_11, amp[587]); 
  FFV1_0(w[1], w[30], w[170], pars->GC_11, amp[588]); 
  FFV1_0(w[1], w[27], w[171], pars->GC_11, amp[589]); 
  FFV1_0(w[29], w[7], w[171], pars->GC_11, amp[590]); 
  FFV1_0(w[29], w[6], w[170], pars->GC_11, amp[591]); 
  FFV1_0(w[1], w[27], w[172], pars->GC_11, amp[592]); 
  FFV1_0(w[1], w[151], w[91], pars->GC_11, amp[593]); 
  FFV1_0(w[29], w[7], w[172], pars->GC_11, amp[594]); 
  FFV1_0(w[152], w[7], w[91], pars->GC_11, amp[595]); 
  FFV1_0(w[1], w[30], w[173], pars->GC_11, amp[596]); 
  FFV1_0(w[1], w[164], w[59], pars->GC_11, amp[597]); 
  FFV1_0(w[29], w[6], w[173], pars->GC_11, amp[598]); 
  FFV1_0(w[152], w[6], w[59], pars->GC_11, amp[599]); 
  FFV1_0(w[169], w[11], w[35], pars->GC_11, amp[600]); 
  FFV1_0(w[169], w[38], w[9], pars->GC_2, amp[601]); 
  FFV1_0(w[60], w[7], w[167], pars->GC_11, amp[602]); 
  FFV1_0(w[57], w[11], w[167], pars->GC_11, amp[603]); 
  FFV1_0(w[60], w[38], w[4], pars->GC_11, amp[604]); 
  FFV1_0(w[174], w[11], w[4], pars->GC_11, amp[605]); 
  FFV1_0(w[169], w[27], w[35], pars->GC_11, amp[606]); 
  FFV2_5_0(w[169], w[38], w[26], pars->GC_51, pars->GC_58, amp[607]); 
  FFV1_0(w[61], w[7], w[167], pars->GC_11, amp[608]); 
  FFV1_0(w[57], w[27], w[167], pars->GC_11, amp[609]); 
  FFV1_0(w[61], w[38], w[4], pars->GC_11, amp[610]); 
  FFV1_0(w[174], w[27], w[4], pars->GC_11, amp[611]); 
  FFV1_0(w[169], w[17], w[43], pars->GC_11, amp[612]); 
  FFV1_0(w[169], w[46], w[9], pars->GC_2, amp[613]); 
  FFV1_0(w[60], w[6], w[156], pars->GC_11, amp[614]); 
  FFV1_0(w[57], w[17], w[156], pars->GC_11, amp[615]); 
  FFV1_0(w[60], w[46], w[4], pars->GC_11, amp[616]); 
  FFV1_0(w[175], w[17], w[4], pars->GC_11, amp[617]); 
  FFV1_0(w[169], w[30], w[43], pars->GC_11, amp[618]); 
  FFV2_5_0(w[169], w[46], w[26], pars->GC_51, pars->GC_58, amp[619]); 
  FFV1_0(w[61], w[6], w[156], pars->GC_11, amp[620]); 
  FFV1_0(w[57], w[30], w[156], pars->GC_11, amp[621]); 
  FFV1_0(w[61], w[46], w[4], pars->GC_11, amp[622]); 
  FFV1_0(w[175], w[30], w[4], pars->GC_11, amp[623]); 
  FFV1_0(w[0], w[17], w[177], pars->GC_11, amp[624]); 
  FFV1_0(w[0], w[11], w[178], pars->GC_11, amp[625]); 
  FFV1_0(w[14], w[7], w[178], pars->GC_11, amp[626]); 
  FFV1_0(w[14], w[6], w[177], pars->GC_11, amp[627]); 
  FFV1_0(w[0], w[11], w[179], pars->GC_11, amp[628]); 
  FFV1_0(w[0], w[147], w[94], pars->GC_11, amp[629]); 
  FFV1_0(w[14], w[7], w[179], pars->GC_11, amp[630]); 
  FFV1_0(w[150], w[7], w[94], pars->GC_11, amp[631]); 
  FFV1_0(w[0], w[17], w[180], pars->GC_11, amp[632]); 
  FFV1_0(w[0], w[162], w[65], pars->GC_11, amp[633]); 
  FFV1_0(w[14], w[6], w[180], pars->GC_11, amp[634]); 
  FFV1_0(w[150], w[6], w[65], pars->GC_11, amp[635]); 
  FFV1_0(w[0], w[30], w[177], pars->GC_11, amp[636]); 
  FFV1_0(w[0], w[27], w[178], pars->GC_11, amp[637]); 
  FFV1_0(w[28], w[7], w[178], pars->GC_11, amp[638]); 
  FFV1_0(w[28], w[6], w[177], pars->GC_11, amp[639]); 
  FFV1_0(w[0], w[27], w[179], pars->GC_11, amp[640]); 
  FFV1_0(w[0], w[151], w[94], pars->GC_11, amp[641]); 
  FFV1_0(w[28], w[7], w[179], pars->GC_11, amp[642]); 
  FFV1_0(w[153], w[7], w[94], pars->GC_11, amp[643]); 
  FFV1_0(w[0], w[30], w[180], pars->GC_11, amp[644]); 
  FFV1_0(w[0], w[164], w[65], pars->GC_11, amp[645]); 
  FFV1_0(w[28], w[6], w[180], pars->GC_11, amp[646]); 
  FFV1_0(w[153], w[6], w[65], pars->GC_11, amp[647]); 
  FFV1_0(w[176], w[11], w[31], pars->GC_11, amp[648]); 
  FFV1_0(w[176], w[34], w[9], pars->GC_2, amp[649]); 
  FFV1_0(w[66], w[7], w[165], pars->GC_11, amp[650]); 
  FFV1_0(w[63], w[11], w[165], pars->GC_11, amp[651]); 
  FFV1_0(w[66], w[34], w[4], pars->GC_11, amp[652]); 
  FFV1_0(w[181], w[11], w[4], pars->GC_11, amp[653]); 
  FFV1_0(w[176], w[27], w[31], pars->GC_11, amp[654]); 
  FFV2_5_0(w[176], w[34], w[26], pars->GC_51, pars->GC_58, amp[655]); 
  FFV1_0(w[67], w[7], w[165], pars->GC_11, amp[656]); 
  FFV1_0(w[63], w[27], w[165], pars->GC_11, amp[657]); 
  FFV1_0(w[67], w[34], w[4], pars->GC_11, amp[658]); 
  FFV1_0(w[181], w[27], w[4], pars->GC_11, amp[659]); 
  FFV1_0(w[176], w[17], w[39], pars->GC_11, amp[660]); 
  FFV1_0(w[176], w[42], w[9], pars->GC_2, amp[661]); 
  FFV1_0(w[66], w[6], w[154], pars->GC_11, amp[662]); 
  FFV1_0(w[63], w[17], w[154], pars->GC_11, amp[663]); 
  FFV1_0(w[66], w[42], w[4], pars->GC_11, amp[664]); 
  FFV1_0(w[182], w[17], w[4], pars->GC_11, amp[665]); 
  FFV1_0(w[176], w[30], w[39], pars->GC_11, amp[666]); 
  FFV2_5_0(w[176], w[42], w[26], pars->GC_51, pars->GC_58, amp[667]); 
  FFV1_0(w[67], w[6], w[154], pars->GC_11, amp[668]); 
  FFV1_0(w[63], w[30], w[154], pars->GC_11, amp[669]); 
  FFV1_0(w[67], w[42], w[4], pars->GC_11, amp[670]); 
  FFV1_0(w[182], w[30], w[4], pars->GC_11, amp[671]); 
  FFV1_0(w[1], w[11], w[183], pars->GC_11, amp[672]); 
  FFV1_0(w[184], w[11], w[5], pars->GC_11, amp[673]); 
  FFV1_0(w[15], w[7], w[183], pars->GC_11, amp[674]); 
  FFV1_0(w[15], w[185], w[5], pars->GC_11, amp[675]); 
  FFV1_0(w[1], w[11], w[186], pars->GC_11, amp[676]); 
  FFV1_0(w[1], w[147], w[105], pars->GC_11, amp[677]); 
  FFV1_0(w[15], w[7], w[186], pars->GC_11, amp[678]); 
  FFV1_0(w[148], w[7], w[105], pars->GC_11, amp[679]); 
  FFV1_0(w[33], w[147], w[5], pars->GC_11, amp[680]); 
  FFV1_0(w[187], w[11], w[5], pars->GC_11, amp[681]); 
  FFV1_0(w[148], w[34], w[5], pars->GC_11, amp[682]); 
  FFV1_0(w[15], w[188], w[5], pars->GC_11, amp[683]); 
  FFV1_0(w[1], w[11], w[189], pars->GC_11, amp[684]); 
  FFV1_0(w[1], w[11], w[190], pars->GC_11, amp[685]); 
  FFV1_0(w[1], w[11], w[191], pars->GC_11, amp[686]); 
  FFV1_0(w[15], w[7], w[189], pars->GC_11, amp[687]); 
  FFV1_0(w[15], w[7], w[190], pars->GC_11, amp[688]); 
  FFV1_0(w[15], w[7], w[191], pars->GC_11, amp[689]); 
  FFV1_0(w[1], w[27], w[183], pars->GC_11, amp[690]); 
  FFV1_0(w[184], w[27], w[5], pars->GC_11, amp[691]); 
  FFV1_0(w[29], w[7], w[183], pars->GC_11, amp[692]); 
  FFV1_0(w[29], w[185], w[5], pars->GC_11, amp[693]); 
  FFV1_0(w[1], w[27], w[186], pars->GC_11, amp[694]); 
  FFV1_0(w[1], w[151], w[105], pars->GC_11, amp[695]); 
  FFV1_0(w[29], w[7], w[186], pars->GC_11, amp[696]); 
  FFV1_0(w[152], w[7], w[105], pars->GC_11, amp[697]); 
  FFV1_0(w[33], w[151], w[5], pars->GC_11, amp[698]); 
  FFV1_0(w[187], w[27], w[5], pars->GC_11, amp[699]); 
  FFV1_0(w[152], w[34], w[5], pars->GC_11, amp[700]); 
  FFV1_0(w[29], w[188], w[5], pars->GC_11, amp[701]); 
  FFV1_0(w[1], w[27], w[189], pars->GC_11, amp[702]); 
  FFV1_0(w[1], w[27], w[190], pars->GC_11, amp[703]); 
  FFV1_0(w[1], w[27], w[191], pars->GC_11, amp[704]); 
  FFV1_0(w[29], w[7], w[189], pars->GC_11, amp[705]); 
  FFV1_0(w[29], w[7], w[190], pars->GC_11, amp[706]); 
  FFV1_0(w[29], w[7], w[191], pars->GC_11, amp[707]); 
  FFV1_0(w[0], w[11], w[192], pars->GC_11, amp[708]); 
  FFV1_0(w[193], w[11], w[5], pars->GC_11, amp[709]); 
  FFV1_0(w[14], w[7], w[192], pars->GC_11, amp[710]); 
  FFV1_0(w[14], w[194], w[5], pars->GC_11, amp[711]); 
  FFV1_0(w[0], w[11], w[195], pars->GC_11, amp[712]); 
  FFV1_0(w[0], w[147], w[106], pars->GC_11, amp[713]); 
  FFV1_0(w[14], w[7], w[195], pars->GC_11, amp[714]); 
  FFV1_0(w[150], w[7], w[106], pars->GC_11, amp[715]); 
  FFV1_0(w[37], w[147], w[5], pars->GC_11, amp[716]); 
  FFV1_0(w[196], w[11], w[5], pars->GC_11, amp[717]); 
  FFV1_0(w[150], w[38], w[5], pars->GC_11, amp[718]); 
  FFV1_0(w[14], w[197], w[5], pars->GC_11, amp[719]); 
  FFV1_0(w[0], w[11], w[198], pars->GC_11, amp[720]); 
  FFV1_0(w[0], w[11], w[199], pars->GC_11, amp[721]); 
  FFV1_0(w[0], w[11], w[200], pars->GC_11, amp[722]); 
  FFV1_0(w[14], w[7], w[198], pars->GC_11, amp[723]); 
  FFV1_0(w[14], w[7], w[199], pars->GC_11, amp[724]); 
  FFV1_0(w[14], w[7], w[200], pars->GC_11, amp[725]); 
  FFV1_0(w[0], w[27], w[192], pars->GC_11, amp[726]); 
  FFV1_0(w[193], w[27], w[5], pars->GC_11, amp[727]); 
  FFV1_0(w[28], w[7], w[192], pars->GC_11, amp[728]); 
  FFV1_0(w[28], w[194], w[5], pars->GC_11, amp[729]); 
  FFV1_0(w[0], w[27], w[195], pars->GC_11, amp[730]); 
  FFV1_0(w[0], w[151], w[106], pars->GC_11, amp[731]); 
  FFV1_0(w[28], w[7], w[195], pars->GC_11, amp[732]); 
  FFV1_0(w[153], w[7], w[106], pars->GC_11, amp[733]); 
  FFV1_0(w[37], w[151], w[5], pars->GC_11, amp[734]); 
  FFV1_0(w[196], w[27], w[5], pars->GC_11, amp[735]); 
  FFV1_0(w[153], w[38], w[5], pars->GC_11, amp[736]); 
  FFV1_0(w[28], w[197], w[5], pars->GC_11, amp[737]); 
  FFV1_0(w[0], w[27], w[198], pars->GC_11, amp[738]); 
  FFV1_0(w[0], w[27], w[199], pars->GC_11, amp[739]); 
  FFV1_0(w[0], w[27], w[200], pars->GC_11, amp[740]); 
  FFV1_0(w[28], w[7], w[198], pars->GC_11, amp[741]); 
  FFV1_0(w[28], w[7], w[199], pars->GC_11, amp[742]); 
  FFV1_0(w[28], w[7], w[200], pars->GC_11, amp[743]); 
  FFV1_0(w[1], w[17], w[201], pars->GC_11, amp[744]); 
  FFV1_0(w[202], w[17], w[5], pars->GC_11, amp[745]); 
  FFV1_0(w[15], w[6], w[201], pars->GC_11, amp[746]); 
  FFV1_0(w[15], w[203], w[5], pars->GC_11, amp[747]); 
  FFV1_0(w[1], w[17], w[204], pars->GC_11, amp[748]); 
  FFV1_0(w[1], w[162], w[78], pars->GC_11, amp[749]); 
  FFV1_0(w[15], w[6], w[204], pars->GC_11, amp[750]); 
  FFV1_0(w[148], w[6], w[78], pars->GC_11, amp[751]); 
  FFV1_0(w[40], w[162], w[5], pars->GC_11, amp[752]); 
  FFV1_0(w[205], w[17], w[5], pars->GC_11, amp[753]); 
  FFV1_0(w[148], w[42], w[5], pars->GC_11, amp[754]); 
  FFV1_0(w[15], w[206], w[5], pars->GC_11, amp[755]); 
  FFV1_0(w[1], w[17], w[207], pars->GC_11, amp[756]); 
  FFV1_0(w[1], w[17], w[208], pars->GC_11, amp[757]); 
  FFV1_0(w[1], w[17], w[209], pars->GC_11, amp[758]); 
  FFV1_0(w[15], w[6], w[207], pars->GC_11, amp[759]); 
  FFV1_0(w[15], w[6], w[208], pars->GC_11, amp[760]); 
  FFV1_0(w[15], w[6], w[209], pars->GC_11, amp[761]); 
  FFV1_0(w[1], w[30], w[201], pars->GC_11, amp[762]); 
  FFV1_0(w[202], w[30], w[5], pars->GC_11, amp[763]); 
  FFV1_0(w[29], w[6], w[201], pars->GC_11, amp[764]); 
  FFV1_0(w[29], w[203], w[5], pars->GC_11, amp[765]); 
  FFV1_0(w[1], w[30], w[204], pars->GC_11, amp[766]); 
  FFV1_0(w[1], w[164], w[78], pars->GC_11, amp[767]); 
  FFV1_0(w[29], w[6], w[204], pars->GC_11, amp[768]); 
  FFV1_0(w[152], w[6], w[78], pars->GC_11, amp[769]); 
  FFV1_0(w[40], w[164], w[5], pars->GC_11, amp[770]); 
  FFV1_0(w[205], w[30], w[5], pars->GC_11, amp[771]); 
  FFV1_0(w[152], w[42], w[5], pars->GC_11, amp[772]); 
  FFV1_0(w[29], w[206], w[5], pars->GC_11, amp[773]); 
  FFV1_0(w[1], w[30], w[207], pars->GC_11, amp[774]); 
  FFV1_0(w[1], w[30], w[208], pars->GC_11, amp[775]); 
  FFV1_0(w[1], w[30], w[209], pars->GC_11, amp[776]); 
  FFV1_0(w[29], w[6], w[207], pars->GC_11, amp[777]); 
  FFV1_0(w[29], w[6], w[208], pars->GC_11, amp[778]); 
  FFV1_0(w[29], w[6], w[209], pars->GC_11, amp[779]); 
  FFV1_0(w[0], w[17], w[210], pars->GC_11, amp[780]); 
  FFV1_0(w[211], w[17], w[5], pars->GC_11, amp[781]); 
  FFV1_0(w[14], w[6], w[210], pars->GC_11, amp[782]); 
  FFV1_0(w[14], w[212], w[5], pars->GC_11, amp[783]); 
  FFV1_0(w[0], w[17], w[213], pars->GC_11, amp[784]); 
  FFV1_0(w[0], w[162], w[79], pars->GC_11, amp[785]); 
  FFV1_0(w[14], w[6], w[213], pars->GC_11, amp[786]); 
  FFV1_0(w[150], w[6], w[79], pars->GC_11, amp[787]); 
  FFV1_0(w[44], w[162], w[5], pars->GC_11, amp[788]); 
  FFV1_0(w[214], w[17], w[5], pars->GC_11, amp[789]); 
  FFV1_0(w[150], w[46], w[5], pars->GC_11, amp[790]); 
  FFV1_0(w[14], w[215], w[5], pars->GC_11, amp[791]); 
  FFV1_0(w[0], w[17], w[216], pars->GC_11, amp[792]); 
  FFV1_0(w[0], w[17], w[217], pars->GC_11, amp[793]); 
  FFV1_0(w[0], w[17], w[218], pars->GC_11, amp[794]); 
  FFV1_0(w[14], w[6], w[216], pars->GC_11, amp[795]); 
  FFV1_0(w[14], w[6], w[217], pars->GC_11, amp[796]); 
  FFV1_0(w[14], w[6], w[218], pars->GC_11, amp[797]); 
  FFV1_0(w[0], w[30], w[210], pars->GC_11, amp[798]); 
  FFV1_0(w[211], w[30], w[5], pars->GC_11, amp[799]); 
  FFV1_0(w[28], w[6], w[210], pars->GC_11, amp[800]); 
  FFV1_0(w[28], w[212], w[5], pars->GC_11, amp[801]); 
  FFV1_0(w[0], w[30], w[213], pars->GC_11, amp[802]); 
  FFV1_0(w[0], w[164], w[79], pars->GC_11, amp[803]); 
  FFV1_0(w[28], w[6], w[213], pars->GC_11, amp[804]); 
  FFV1_0(w[153], w[6], w[79], pars->GC_11, amp[805]); 
  FFV1_0(w[44], w[164], w[5], pars->GC_11, amp[806]); 
  FFV1_0(w[214], w[30], w[5], pars->GC_11, amp[807]); 
  FFV1_0(w[153], w[46], w[5], pars->GC_11, amp[808]); 
  FFV1_0(w[28], w[215], w[5], pars->GC_11, amp[809]); 
  FFV1_0(w[0], w[30], w[216], pars->GC_11, amp[810]); 
  FFV1_0(w[0], w[30], w[217], pars->GC_11, amp[811]); 
  FFV1_0(w[0], w[30], w[218], pars->GC_11, amp[812]); 
  FFV1_0(w[28], w[6], w[216], pars->GC_11, amp[813]); 
  FFV1_0(w[28], w[6], w[217], pars->GC_11, amp[814]); 
  FFV1_0(w[28], w[6], w[218], pars->GC_11, amp[815]); 

}
double CPPProcess::matrix_uu_epemgguu_no_h() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 816; 
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
  jamp[0] = +1./2. * (-1./3. * std::complex<double> (0, 1) * amp[0] -
      std::complex<double> (0, 1) * amp[1] - std::complex<double> (0, 1) *
      amp[2] - 1./3. * std::complex<double> (0, 1) * amp[3] -
      std::complex<double> (0, 1) * amp[8] - 1./3. * std::complex<double> (0,
      1) * amp[9] - 1./3. * std::complex<double> (0, 1) * amp[10] -
      std::complex<double> (0, 1) * amp[11] - 1./3. * std::complex<double> (0,
      1) * amp[16] - std::complex<double> (0, 1) * amp[17] -
      std::complex<double> (0, 1) * amp[18] - 1./3. * std::complex<double> (0,
      1) * amp[19] - std::complex<double> (0, 1) * amp[24] - 1./3. *
      std::complex<double> (0, 1) * amp[25] - 1./3. * std::complex<double> (0,
      1) * amp[26] - std::complex<double> (0, 1) * amp[27] + amp[48] + amp[49]
      - std::complex<double> (0, 1) * amp[52] - std::complex<double> (0, 1) *
      amp[53] - std::complex<double> (0, 1) * amp[55] + amp[56] + amp[57] -
      std::complex<double> (0, 1) * amp[60] - std::complex<double> (0, 1) *
      amp[61] - std::complex<double> (0, 1) * amp[63] - std::complex<double>
      (0, 1) * amp[64] - std::complex<double> (0, 1) * amp[65] - amp[66] -
      amp[67] - std::complex<double> (0, 1) * amp[70] - std::complex<double>
      (0, 1) * amp[72] - std::complex<double> (0, 1) * amp[73] - amp[74] -
      amp[75] - std::complex<double> (0, 1) * amp[78] - 1./3. *
      std::complex<double> (0, 1) * amp[80] - 1./3. * std::complex<double> (0,
      1) * amp[81] - 1./3. * std::complex<double> (0, 1) * amp[84] - 1./3. *
      std::complex<double> (0, 1) * amp[85] - 1./3. * std::complex<double> (0,
      1) * amp[86] - 1./3. * std::complex<double> (0, 1) * amp[87] - 1./3. *
      std::complex<double> (0, 1) * amp[88] - 1./3. * std::complex<double> (0,
      1) * amp[89] - 1./3. * std::complex<double> (0, 1) * amp[92] - 1./3. *
      std::complex<double> (0, 1) * amp[93] - 1./3. * std::complex<double> (0,
      1) * amp[94] - 1./3. * std::complex<double> (0, 1) * amp[95] + 1./3. *
      amp[112] + 1./3. * amp[113] + amp[114] + amp[115] + amp[116] + amp[117] +
      1./3. * amp[118] + 1./3. * amp[119] + amp[120] + amp[121] + amp[122] +
      amp[123] + 1./3. * amp[124] + 1./3. * amp[125] + 1./3. * amp[126] + 1./3.
      * amp[127] + 1./3. * amp[144] + amp[145] + amp[146] + 1./3. * amp[147] -
      std::complex<double> (0, 1) * amp[152] - std::complex<double> (0, 1) *
      amp[154] + amp[155] + 1./3. * amp[156] + amp[157] + amp[158] + 1./3. *
      amp[159] - std::complex<double> (0, 1) * amp[164] - std::complex<double>
      (0, 1) * amp[166] + amp[167] + amp[168] + amp[169] - std::complex<double>
      (0, 1) * amp[170] + amp[171] - std::complex<double> (0, 1) * amp[173] +
      amp[174] + amp[175] - std::complex<double> (0, 1) * amp[176] + amp[177] -
      std::complex<double> (0, 1) * amp[179] + 1./3. * amp[180] + 1./3. *
      amp[181] + 1./3. * amp[183] + 1./3. * amp[184] + 1./3. * amp[186] + 1./3.
      * amp[187] + 1./3. * amp[189] + 1./3. * amp[190] + amp[576] + 1./3. *
      amp[577] + 1./3. * amp[578] + amp[579] - std::complex<double> (0, 1) *
      amp[584] + amp[585] - std::complex<double> (0, 1) * amp[586] + amp[588] +
      1./3. * amp[589] + 1./3. * amp[590] + amp[591] - std::complex<double> (0,
      1) * amp[596] + amp[597] - std::complex<double> (0, 1) * amp[598] +
      amp[600] + amp[601] + std::complex<double> (0, 1) * amp[602] +
      std::complex<double> (0, 1) * amp[603] + amp[604] + amp[606] + amp[607] +
      std::complex<double> (0, 1) * amp[608] + std::complex<double> (0, 1) *
      amp[609] + amp[610] + 1./3. * amp[612] + 1./3. * amp[613] + 1./3. *
      amp[616] + 1./3. * amp[617] + 1./3. * amp[618] + 1./3. * amp[619] + 1./3.
      * amp[622] + 1./3. * amp[623] + amp[708] + amp[710] +
      std::complex<double> (0, 1) * amp[711] + amp[719] - amp[722] - amp[721] -
      amp[725] - amp[724] + amp[726] + amp[728] + std::complex<double> (0, 1) *
      amp[729] + amp[737] - amp[740] - amp[739] - amp[743] - amp[742] -
      amp[748] - std::complex<double> (0, 1) * amp[749] - amp[750] + amp[752] +
      amp[758] - amp[756] + amp[761] - amp[759] - amp[766] -
      std::complex<double> (0, 1) * amp[767] - amp[768] + amp[770] + amp[776] -
      amp[774] + amp[779] - amp[777] + 1./3. * amp[788] + 1./3. * amp[791] +
      1./3. * amp[806] + 1./3. * amp[809]);
  jamp[1] = +1./2. * (+std::complex<double> (0, 1) * amp[0] + 1./3. *
      std::complex<double> (0, 1) * amp[1] + 1./3. * std::complex<double> (0,
      1) * amp[2] + std::complex<double> (0, 1) * amp[3] + std::complex<double>
      (0, 1) * amp[12] + 1./3. * std::complex<double> (0, 1) * amp[13] + 1./3.
      * std::complex<double> (0, 1) * amp[14] + std::complex<double> (0, 1) *
      amp[15] + std::complex<double> (0, 1) * amp[16] + 1./3. *
      std::complex<double> (0, 1) * amp[17] + 1./3. * std::complex<double> (0,
      1) * amp[18] + std::complex<double> (0, 1) * amp[19] +
      std::complex<double> (0, 1) * amp[28] + 1./3. * std::complex<double> (0,
      1) * amp[29] + 1./3. * std::complex<double> (0, 1) * amp[30] +
      std::complex<double> (0, 1) * amp[31] - amp[32] - amp[33] +
      std::complex<double> (0, 1) * amp[36] + std::complex<double> (0, 1) *
      amp[37] + std::complex<double> (0, 1) * amp[39] - amp[40] - amp[41] +
      std::complex<double> (0, 1) * amp[44] + std::complex<double> (0, 1) *
      amp[45] + std::complex<double> (0, 1) * amp[47] + 1./3. *
      std::complex<double> (0, 1) * amp[64] + 1./3. * std::complex<double> (0,
      1) * amp[65] + 1./3. * std::complex<double> (0, 1) * amp[68] + 1./3. *
      std::complex<double> (0, 1) * amp[69] + 1./3. * std::complex<double> (0,
      1) * amp[70] + 1./3. * std::complex<double> (0, 1) * amp[71] + 1./3. *
      std::complex<double> (0, 1) * amp[72] + 1./3. * std::complex<double> (0,
      1) * amp[73] + 1./3. * std::complex<double> (0, 1) * amp[76] + 1./3. *
      std::complex<double> (0, 1) * amp[77] + 1./3. * std::complex<double> (0,
      1) * amp[78] + 1./3. * std::complex<double> (0, 1) * amp[79] +
      std::complex<double> (0, 1) * amp[80] + std::complex<double> (0, 1) *
      amp[81] + amp[82] + amp[83] + std::complex<double> (0, 1) * amp[86] +
      std::complex<double> (0, 1) * amp[88] + std::complex<double> (0, 1) *
      amp[89] + amp[90] + amp[91] + std::complex<double> (0, 1) * amp[94] -
      1./3. * amp[128] - 1./3. * amp[129] - amp[130] - amp[131] - amp[132] -
      amp[133] - 1./3. * amp[134] - 1./3. * amp[135] - amp[136] - amp[137] -
      amp[138] - amp[139] - 1./3. * amp[140] - 1./3. * amp[141] - 1./3. *
      amp[142] - 1./3. * amp[143] - amp[144] - 1./3. * amp[145] - 1./3. *
      amp[146] - amp[147] + std::complex<double> (0, 1) * amp[148] +
      std::complex<double> (0, 1) * amp[150] - amp[151] - amp[156] - 1./3. *
      amp[157] - 1./3. * amp[158] - amp[159] + std::complex<double> (0, 1) *
      amp[160] + std::complex<double> (0, 1) * amp[162] - amp[163] - 1./3. *
      amp[168] - 1./3. * amp[169] - 1./3. * amp[171] - 1./3. * amp[172] - 1./3.
      * amp[174] - 1./3. * amp[175] - 1./3. * amp[177] - 1./3. * amp[178] -
      amp[180] - amp[181] + std::complex<double> (0, 1) * amp[182] - amp[183] +
      std::complex<double> (0, 1) * amp[185] - amp[186] - amp[187] +
      std::complex<double> (0, 1) * amp[188] - amp[189] + std::complex<double>
      (0, 1) * amp[191] - amp[624] - 1./3. * amp[625] - 1./3. * amp[626] -
      amp[627] + std::complex<double> (0, 1) * amp[632] - amp[633] +
      std::complex<double> (0, 1) * amp[634] - amp[636] - 1./3. * amp[637] -
      1./3. * amp[638] - amp[639] + std::complex<double> (0, 1) * amp[644] -
      amp[645] + std::complex<double> (0, 1) * amp[646] - amp[648] - amp[649] -
      std::complex<double> (0, 1) * amp[650] - std::complex<double> (0, 1) *
      amp[651] - amp[652] - amp[654] - amp[655] - std::complex<double> (0, 1) *
      amp[656] - std::complex<double> (0, 1) * amp[657] - amp[658] - 1./3. *
      amp[660] - 1./3. * amp[661] - 1./3. * amp[664] - 1./3. * amp[665] - 1./3.
      * amp[666] - 1./3. * amp[667] - 1./3. * amp[670] - 1./3. * amp[671] -
      amp[672] - amp[674] - std::complex<double> (0, 1) * amp[675] - amp[683] +
      amp[686] + amp[685] + amp[689] + amp[688] - amp[690] - amp[692] -
      std::complex<double> (0, 1) * amp[693] - amp[701] + amp[704] + amp[703] +
      amp[707] + amp[706] - 1./3. * amp[752] - 1./3. * amp[755] - 1./3. *
      amp[770] - 1./3. * amp[773] + amp[784] + std::complex<double> (0, 1) *
      amp[785] + amp[786] - amp[788] - amp[794] + amp[792] - amp[797] +
      amp[795] + amp[802] + std::complex<double> (0, 1) * amp[803] + amp[804] -
      amp[806] - amp[812] + amp[810] - amp[815] + amp[813]);
  jamp[2] = +1./2. * (+1./3. * std::complex<double> (0, 1) * amp[4] +
      std::complex<double> (0, 1) * amp[5] + std::complex<double> (0, 1) *
      amp[6] + 1./3. * std::complex<double> (0, 1) * amp[7] + 1./3. *
      std::complex<double> (0, 1) * amp[8] + std::complex<double> (0, 1) *
      amp[9] + std::complex<double> (0, 1) * amp[10] + 1./3. *
      std::complex<double> (0, 1) * amp[11] + 1./3. * std::complex<double> (0,
      1) * amp[20] + std::complex<double> (0, 1) * amp[21] +
      std::complex<double> (0, 1) * amp[22] + 1./3. * std::complex<double> (0,
      1) * amp[23] + 1./3. * std::complex<double> (0, 1) * amp[24] +
      std::complex<double> (0, 1) * amp[25] + std::complex<double> (0, 1) *
      amp[26] + 1./3. * std::complex<double> (0, 1) * amp[27] + amp[32] +
      amp[33] + std::complex<double> (0, 1) * amp[34] + std::complex<double>
      (0, 1) * amp[35] + std::complex<double> (0, 1) * amp[38] + amp[40] +
      amp[41] + std::complex<double> (0, 1) * amp[42] + std::complex<double>
      (0, 1) * amp[43] + std::complex<double> (0, 1) * amp[46] + 1./3. *
      std::complex<double> (0, 1) * amp[50] + 1./3. * std::complex<double> (0,
      1) * amp[51] + 1./3. * std::complex<double> (0, 1) * amp[52] + 1./3. *
      std::complex<double> (0, 1) * amp[53] + 1./3. * std::complex<double> (0,
      1) * amp[54] + 1./3. * std::complex<double> (0, 1) * amp[55] + 1./3. *
      std::complex<double> (0, 1) * amp[58] + 1./3. * std::complex<double> (0,
      1) * amp[59] + 1./3. * std::complex<double> (0, 1) * amp[60] + 1./3. *
      std::complex<double> (0, 1) * amp[61] + 1./3. * std::complex<double> (0,
      1) * amp[62] + 1./3. * std::complex<double> (0, 1) * amp[63] - amp[82] -
      amp[83] + std::complex<double> (0, 1) * amp[84] + std::complex<double>
      (0, 1) * amp[85] + std::complex<double> (0, 1) * amp[87] - amp[90] -
      amp[91] + std::complex<double> (0, 1) * amp[92] + std::complex<double>
      (0, 1) * amp[93] + std::complex<double> (0, 1) * amp[95] - 1./3. *
      amp[208] - 1./3. * amp[209] - amp[210] - amp[211] - amp[212] - amp[213] -
      1./3. * amp[214] - 1./3. * amp[215] - amp[216] - amp[217] - amp[218] -
      amp[219] - 1./3. * amp[220] - 1./3. * amp[221] - 1./3. * amp[222] - 1./3.
      * amp[223] - 1./3. * amp[240] - amp[241] - amp[242] - 1./3. * amp[243] +
      std::complex<double> (0, 1) * amp[248] + std::complex<double> (0, 1) *
      amp[250] - amp[251] - 1./3. * amp[252] - amp[253] - amp[254] - 1./3. *
      amp[255] + std::complex<double> (0, 1) * amp[260] + std::complex<double>
      (0, 1) * amp[262] - amp[263] - amp[264] - amp[265] + std::complex<double>
      (0, 1) * amp[266] - amp[267] + std::complex<double> (0, 1) * amp[269] -
      amp[270] - amp[271] + std::complex<double> (0, 1) * amp[272] - amp[273] +
      std::complex<double> (0, 1) * amp[275] - 1./3. * amp[276] - 1./3. *
      amp[277] - 1./3. * amp[279] - 1./3. * amp[280] - 1./3. * amp[282] - 1./3.
      * amp[283] - 1./3. * amp[285] - 1./3. * amp[286] - 1./3. * amp[576] -
      amp[577] - amp[578] - 1./3. * amp[579] + std::complex<double> (0, 1) *
      amp[580] - amp[581] + std::complex<double> (0, 1) * amp[582] - 1./3. *
      amp[588] - amp[589] - amp[590] - 1./3. * amp[591] + std::complex<double>
      (0, 1) * amp[592] - amp[593] + std::complex<double> (0, 1) * amp[594] -
      1./3. * amp[600] - 1./3. * amp[601] - 1./3. * amp[604] - 1./3. * amp[605]
      - 1./3. * amp[606] - 1./3. * amp[607] - 1./3. * amp[610] - 1./3. *
      amp[611] - amp[612] - amp[613] - std::complex<double> (0, 1) * amp[614] -
      std::complex<double> (0, 1) * amp[615] - amp[616] - amp[618] - amp[619] -
      std::complex<double> (0, 1) * amp[620] - std::complex<double> (0, 1) *
      amp[621] - amp[622] + amp[676] + std::complex<double> (0, 1) * amp[677] +
      amp[678] - amp[680] - amp[686] + amp[684] - amp[689] + amp[687] +
      amp[694] + std::complex<double> (0, 1) * amp[695] + amp[696] - amp[698] -
      amp[704] + amp[702] - amp[707] + amp[705] - 1./3. * amp[716] - 1./3. *
      amp[719] - 1./3. * amp[734] - 1./3. * amp[737] - amp[780] - amp[782] -
      std::complex<double> (0, 1) * amp[783] - amp[791] + amp[794] + amp[793] +
      amp[797] + amp[796] - amp[798] - amp[800] - std::complex<double> (0, 1) *
      amp[801] - amp[809] + amp[812] + amp[811] + amp[815] + amp[814]);
  jamp[3] = +1./2. * (-std::complex<double> (0, 1) * amp[4] - 1./3. *
      std::complex<double> (0, 1) * amp[5] - 1./3. * std::complex<double> (0,
      1) * amp[6] - std::complex<double> (0, 1) * amp[7] - 1./3. *
      std::complex<double> (0, 1) * amp[12] - std::complex<double> (0, 1) *
      amp[13] - std::complex<double> (0, 1) * amp[14] - 1./3. *
      std::complex<double> (0, 1) * amp[15] - std::complex<double> (0, 1) *
      amp[20] - 1./3. * std::complex<double> (0, 1) * amp[21] - 1./3. *
      std::complex<double> (0, 1) * amp[22] - std::complex<double> (0, 1) *
      amp[23] - 1./3. * std::complex<double> (0, 1) * amp[28] -
      std::complex<double> (0, 1) * amp[29] - std::complex<double> (0, 1) *
      amp[30] - 1./3. * std::complex<double> (0, 1) * amp[31] - 1./3. *
      std::complex<double> (0, 1) * amp[34] - 1./3. * std::complex<double> (0,
      1) * amp[35] - 1./3. * std::complex<double> (0, 1) * amp[36] - 1./3. *
      std::complex<double> (0, 1) * amp[37] - 1./3. * std::complex<double> (0,
      1) * amp[38] - 1./3. * std::complex<double> (0, 1) * amp[39] - 1./3. *
      std::complex<double> (0, 1) * amp[42] - 1./3. * std::complex<double> (0,
      1) * amp[43] - 1./3. * std::complex<double> (0, 1) * amp[44] - 1./3. *
      std::complex<double> (0, 1) * amp[45] - 1./3. * std::complex<double> (0,
      1) * amp[46] - 1./3. * std::complex<double> (0, 1) * amp[47] - amp[48] -
      amp[49] - std::complex<double> (0, 1) * amp[50] - std::complex<double>
      (0, 1) * amp[51] - std::complex<double> (0, 1) * amp[54] - amp[56] -
      amp[57] - std::complex<double> (0, 1) * amp[58] - std::complex<double>
      (0, 1) * amp[59] - std::complex<double> (0, 1) * amp[62] + amp[66] +
      amp[67] - std::complex<double> (0, 1) * amp[68] - std::complex<double>
      (0, 1) * amp[69] - std::complex<double> (0, 1) * amp[71] + amp[74] +
      amp[75] - std::complex<double> (0, 1) * amp[76] - std::complex<double>
      (0, 1) * amp[77] - std::complex<double> (0, 1) * amp[79] + 1./3. *
      amp[224] + 1./3. * amp[225] + amp[226] + amp[227] + amp[228] + amp[229] +
      1./3. * amp[230] + 1./3. * amp[231] + amp[232] + amp[233] + amp[234] +
      amp[235] + 1./3. * amp[236] + 1./3. * amp[237] + 1./3. * amp[238] + 1./3.
      * amp[239] + amp[240] + 1./3. * amp[241] + 1./3. * amp[242] + amp[243] -
      std::complex<double> (0, 1) * amp[244] - std::complex<double> (0, 1) *
      amp[246] + amp[247] + amp[252] + 1./3. * amp[253] + 1./3. * amp[254] +
      amp[255] - std::complex<double> (0, 1) * amp[256] - std::complex<double>
      (0, 1) * amp[258] + amp[259] + 1./3. * amp[264] + 1./3. * amp[265] +
      1./3. * amp[267] + 1./3. * amp[268] + 1./3. * amp[270] + 1./3. * amp[271]
      + 1./3. * amp[273] + 1./3. * amp[274] + amp[276] + amp[277] -
      std::complex<double> (0, 1) * amp[278] + amp[279] - std::complex<double>
      (0, 1) * amp[281] + amp[282] + amp[283] - std::complex<double> (0, 1) *
      amp[284] + amp[285] - std::complex<double> (0, 1) * amp[287] + 1./3. *
      amp[624] + amp[625] + amp[626] + 1./3. * amp[627] - std::complex<double>
      (0, 1) * amp[628] + amp[629] - std::complex<double> (0, 1) * amp[630] +
      1./3. * amp[636] + amp[637] + amp[638] + 1./3. * amp[639] -
      std::complex<double> (0, 1) * amp[640] + amp[641] - std::complex<double>
      (0, 1) * amp[642] + 1./3. * amp[648] + 1./3. * amp[649] + 1./3. *
      amp[652] + 1./3. * amp[653] + 1./3. * amp[654] + 1./3. * amp[655] + 1./3.
      * amp[658] + 1./3. * amp[659] + amp[660] + amp[661] +
      std::complex<double> (0, 1) * amp[662] + std::complex<double> (0, 1) *
      amp[663] + amp[664] + amp[666] + amp[667] + std::complex<double> (0, 1) *
      amp[668] + std::complex<double> (0, 1) * amp[669] + amp[670] + 1./3. *
      amp[680] + 1./3. * amp[683] + 1./3. * amp[698] + 1./3. * amp[701] -
      amp[712] - std::complex<double> (0, 1) * amp[713] - amp[714] + amp[716] +
      amp[722] - amp[720] + amp[725] - amp[723] - amp[730] -
      std::complex<double> (0, 1) * amp[731] - amp[732] + amp[734] + amp[740] -
      amp[738] + amp[743] - amp[741] + amp[744] + amp[746] +
      std::complex<double> (0, 1) * amp[747] + amp[755] - amp[758] - amp[757] -
      amp[761] - amp[760] + amp[762] + amp[764] + std::complex<double> (0, 1) *
      amp[765] + amp[773] - amp[776] - amp[775] - amp[779] - amp[778]);
  jamp[4] = +1./2. * (+amp[96] + 1./3. * amp[97] + 1./3. * amp[98] + 1./3. *
      amp[99] + amp[100] + amp[101] + amp[102] + 1./3. * amp[103] + amp[104] +
      1./3. * amp[105] + 1./3. * amp[106] + 1./3. * amp[107] + amp[108] +
      amp[109] + amp[110] + 1./3. * amp[111] + amp[128] + amp[129] + 1./3. *
      amp[130] + 1./3. * amp[131] + 1./3. * amp[132] + 1./3. * amp[133] +
      amp[134] + amp[135] + 1./3. * amp[136] + 1./3. * amp[137] + 1./3. *
      amp[138] + 1./3. * amp[139] + amp[140] + amp[141] + amp[142] + amp[143] +
      1./3. * amp[149] + 1./3. * amp[151] + std::complex<double> (0, 1) *
      amp[152] + amp[153] + std::complex<double> (0, 1) * amp[154] + 1./3. *
      amp[161] + 1./3. * amp[163] + std::complex<double> (0, 1) * amp[164] +
      amp[165] + std::complex<double> (0, 1) * amp[166] + std::complex<double>
      (0, 1) * amp[170] + amp[172] + std::complex<double> (0, 1) * amp[173] +
      std::complex<double> (0, 1) * amp[176] + amp[178] + std::complex<double>
      (0, 1) * amp[179] + amp[304] + amp[305] + 1./3. * amp[306] + 1./3. *
      amp[307] + 1./3. * amp[308] + 1./3. * amp[309] + amp[310] + amp[311] +
      1./3. * amp[312] + 1./3. * amp[313] + 1./3. * amp[314] + 1./3. * amp[315]
      + amp[316] + amp[317] + amp[318] + amp[319] + amp[320] + 1./3. * amp[321]
      + 1./3. * amp[322] + 1./3. * amp[323] + amp[324] + amp[325] + amp[326] +
      1./3. * amp[327] + amp[328] + 1./3. * amp[329] + 1./3. * amp[330] + 1./3.
      * amp[331] + amp[332] + amp[333] + amp[334] + 1./3. * amp[335] + 1./3. *
      amp[341] + 1./3. * amp[343] - std::complex<double> (0, 1) * amp[344] -
      std::complex<double> (0, 1) * amp[346] + amp[347] + 1./3. * amp[353] +
      1./3. * amp[355] - std::complex<double> (0, 1) * amp[356] -
      std::complex<double> (0, 1) * amp[358] + amp[359] - std::complex<double>
      (0, 1) * amp[362] + amp[364] - std::complex<double> (0, 1) * amp[365] -
      std::complex<double> (0, 1) * amp[368] + amp[370] - std::complex<double>
      (0, 1) * amp[371] - std::complex<double> (0, 1) * amp[532] + amp[533] -
      std::complex<double> (0, 1) * amp[534] + 1./3. * amp[537] + 1./3. *
      amp[539] - std::complex<double> (0, 1) * amp[544] + amp[545] -
      std::complex<double> (0, 1) * amp[546] + 1./3. * amp[549] + 1./3. *
      amp[551] + std::complex<double> (0, 1) * amp[566] + std::complex<double>
      (0, 1) * amp[567] + amp[569] + std::complex<double> (0, 1) * amp[572] +
      std::complex<double> (0, 1) * amp[573] + amp[575] + std::complex<double>
      (0, 1) * amp[628] + std::complex<double> (0, 1) * amp[630] + amp[631] +
      1./3. * amp[633] + 1./3. * amp[635] + std::complex<double> (0, 1) *
      amp[640] + std::complex<double> (0, 1) * amp[642] + amp[643] + 1./3. *
      amp[645] + 1./3. * amp[647] - std::complex<double> (0, 1) * amp[662] -
      std::complex<double> (0, 1) * amp[663] + amp[665] - std::complex<double>
      (0, 1) * amp[668] - std::complex<double> (0, 1) * amp[669] + amp[671] -
      amp[708] + std::complex<double> (0, 1) * amp[709] - amp[710] + amp[712] +
      amp[714] - std::complex<double> (0, 1) * amp[715] + amp[720] + amp[721] +
      amp[723] + amp[724] - amp[726] + std::complex<double> (0, 1) * amp[727] -
      amp[728] + amp[730] + amp[732] - std::complex<double> (0, 1) * amp[733] +
      amp[738] + amp[739] + amp[741] + amp[742] - amp[744] - amp[746] -
      std::complex<double> (0, 1) * amp[747] + amp[748] + std::complex<double>
      (0, 1) * amp[749] + amp[750] + amp[756] + amp[757] + amp[759] + amp[760]
      - amp[762] - amp[764] - std::complex<double> (0, 1) * amp[765] + amp[766]
      + std::complex<double> (0, 1) * amp[767] + amp[768] + amp[774] + amp[775]
      + amp[777] + amp[778]);
  jamp[5] = +1./2. * (-1./3. * amp[96] - amp[97] - amp[98] - amp[99] - 1./3. *
      amp[100] - 1./3. * amp[101] - 1./3. * amp[102] - amp[103] - 1./3. *
      amp[104] - amp[105] - amp[106] - amp[107] - 1./3. * amp[108] - 1./3. *
      amp[109] - 1./3. * amp[110] - amp[111] - amp[112] - amp[113] - 1./3. *
      amp[114] - 1./3. * amp[115] - 1./3. * amp[116] - 1./3. * amp[117] -
      amp[118] - amp[119] - 1./3. * amp[120] - 1./3. * amp[121] - 1./3. *
      amp[122] - 1./3. * amp[123] - amp[124] - amp[125] - amp[126] - amp[127] -
      std::complex<double> (0, 1) * amp[148] - amp[149] - std::complex<double>
      (0, 1) * amp[150] - 1./3. * amp[153] - 1./3. * amp[155] -
      std::complex<double> (0, 1) * amp[160] - amp[161] - std::complex<double>
      (0, 1) * amp[162] - 1./3. * amp[165] - 1./3. * amp[167] -
      std::complex<double> (0, 1) * amp[182] - amp[184] - std::complex<double>
      (0, 1) * amp[185] - std::complex<double> (0, 1) * amp[188] - amp[190] -
      std::complex<double> (0, 1) * amp[191] - amp[400] - amp[401] - 1./3. *
      amp[402] - 1./3. * amp[403] - 1./3. * amp[404] - 1./3. * amp[405] -
      amp[406] - amp[407] - 1./3. * amp[408] - 1./3. * amp[409] - 1./3. *
      amp[410] - 1./3. * amp[411] - amp[412] - amp[413] - amp[414] - amp[415] -
      amp[416] - 1./3. * amp[417] - 1./3. * amp[418] - 1./3. * amp[419] -
      amp[420] - amp[421] - amp[422] - 1./3. * amp[423] - amp[424] - 1./3. *
      amp[425] - 1./3. * amp[426] - 1./3. * amp[427] - amp[428] - amp[429] -
      amp[430] - 1./3. * amp[431] - 1./3. * amp[437] - 1./3. * amp[439] +
      std::complex<double> (0, 1) * amp[440] + std::complex<double> (0, 1) *
      amp[442] - amp[443] - 1./3. * amp[449] - 1./3. * amp[451] +
      std::complex<double> (0, 1) * amp[452] + std::complex<double> (0, 1) *
      amp[454] - amp[455] + std::complex<double> (0, 1) * amp[458] - amp[460] +
      std::complex<double> (0, 1) * amp[461] + std::complex<double> (0, 1) *
      amp[464] - amp[466] + std::complex<double> (0, 1) * amp[467] - 1./3. *
      amp[533] - 1./3. * amp[535] + std::complex<double> (0, 1) * amp[536] -
      amp[537] + std::complex<double> (0, 1) * amp[538] - 1./3. * amp[545] -
      1./3. * amp[547] + std::complex<double> (0, 1) * amp[548] - amp[549] +
      std::complex<double> (0, 1) * amp[550] - std::complex<double> (0, 1) *
      amp[554] - std::complex<double> (0, 1) * amp[555] - amp[557] -
      std::complex<double> (0, 1) * amp[560] - std::complex<double> (0, 1) *
      amp[561] - amp[563] - std::complex<double> (0, 1) * amp[580] -
      std::complex<double> (0, 1) * amp[582] - amp[583] - 1./3. * amp[585] -
      1./3. * amp[587] - std::complex<double> (0, 1) * amp[592] -
      std::complex<double> (0, 1) * amp[594] - amp[595] - 1./3. * amp[597] -
      1./3. * amp[599] + std::complex<double> (0, 1) * amp[614] +
      std::complex<double> (0, 1) * amp[615] - amp[617] + std::complex<double>
      (0, 1) * amp[620] + std::complex<double> (0, 1) * amp[621] - amp[623] +
      amp[672] - std::complex<double> (0, 1) * amp[673] + amp[674] - amp[676] -
      amp[678] + std::complex<double> (0, 1) * amp[679] - amp[684] - amp[685] -
      amp[687] - amp[688] + amp[690] - std::complex<double> (0, 1) * amp[691] +
      amp[692] - amp[694] - amp[696] + std::complex<double> (0, 1) * amp[697] -
      amp[702] - amp[703] - amp[705] - amp[706] + amp[780] + amp[782] +
      std::complex<double> (0, 1) * amp[783] - amp[784] - std::complex<double>
      (0, 1) * amp[785] - amp[786] - amp[792] - amp[793] - amp[795] - amp[796]
      + amp[798] + amp[800] + std::complex<double> (0, 1) * amp[801] - amp[802]
      - std::complex<double> (0, 1) * amp[803] - amp[804] - amp[810] - amp[811]
      - amp[813] - amp[814]);
  jamp[6] = +1./2. * (-amp[192] - 1./3. * amp[193] - 1./3. * amp[194] - 1./3. *
      amp[195] - amp[196] - amp[197] - amp[198] - 1./3. * amp[199] - amp[200] -
      1./3. * amp[201] - 1./3. * amp[202] - 1./3. * amp[203] - amp[204] -
      amp[205] - amp[206] - 1./3. * amp[207] - amp[224] - amp[225] - 1./3. *
      amp[226] - 1./3. * amp[227] - 1./3. * amp[228] - 1./3. * amp[229] -
      amp[230] - amp[231] - 1./3. * amp[232] - 1./3. * amp[233] - 1./3. *
      amp[234] - 1./3. * amp[235] - amp[236] - amp[237] - amp[238] - amp[239] -
      1./3. * amp[245] - 1./3. * amp[247] - std::complex<double> (0, 1) *
      amp[248] - amp[249] - std::complex<double> (0, 1) * amp[250] - 1./3. *
      amp[257] - 1./3. * amp[259] - std::complex<double> (0, 1) * amp[260] -
      amp[261] - std::complex<double> (0, 1) * amp[262] - std::complex<double>
      (0, 1) * amp[266] - amp[268] - std::complex<double> (0, 1) * amp[269] -
      std::complex<double> (0, 1) * amp[272] - amp[274] - std::complex<double>
      (0, 1) * amp[275] - amp[288] - amp[289] - 1./3. * amp[290] - 1./3. *
      amp[291] - 1./3. * amp[292] - 1./3. * amp[293] - amp[294] - amp[295] -
      1./3. * amp[296] - 1./3. * amp[297] - 1./3. * amp[298] - 1./3. * amp[299]
      - amp[300] - amp[301] - amp[302] - amp[303] - 1./3. * amp[320] - amp[321]
      - amp[322] - amp[323] - 1./3. * amp[324] - 1./3. * amp[325] - 1./3. *
      amp[326] - amp[327] - 1./3. * amp[328] - amp[329] - amp[330] - amp[331] -
      1./3. * amp[332] - 1./3. * amp[333] - 1./3. * amp[334] - amp[335] +
      std::complex<double> (0, 1) * amp[340] + std::complex<double> (0, 1) *
      amp[342] - amp[343] - 1./3. * amp[345] - 1./3. * amp[347] +
      std::complex<double> (0, 1) * amp[352] + std::complex<double> (0, 1) *
      amp[354] - amp[355] - 1./3. * amp[357] - 1./3. * amp[359] +
      std::complex<double> (0, 1) * amp[374] - amp[376] + std::complex<double>
      (0, 1) * amp[377] + std::complex<double> (0, 1) * amp[380] - amp[382] +
      std::complex<double> (0, 1) * amp[383] + std::complex<double> (0, 1) *
      amp[484] - amp[485] + std::complex<double> (0, 1) * amp[486] - 1./3. *
      amp[489] - 1./3. * amp[491] + std::complex<double> (0, 1) * amp[496] -
      amp[497] + std::complex<double> (0, 1) * amp[498] - 1./3. * amp[501] -
      1./3. * amp[503] - std::complex<double> (0, 1) * amp[518] -
      std::complex<double> (0, 1) * amp[519] - amp[521] - std::complex<double>
      (0, 1) * amp[524] - std::complex<double> (0, 1) * amp[525] - amp[527] -
      1./3. * amp[629] - 1./3. * amp[631] - std::complex<double> (0, 1) *
      amp[632] - std::complex<double> (0, 1) * amp[634] - amp[635] - 1./3. *
      amp[641] - 1./3. * amp[643] - std::complex<double> (0, 1) * amp[644] -
      std::complex<double> (0, 1) * amp[646] - amp[647] + std::complex<double>
      (0, 1) * amp[650] + std::complex<double> (0, 1) * amp[651] - amp[653] +
      std::complex<double> (0, 1) * amp[656] + std::complex<double> (0, 1) *
      amp[657] - amp[659] + amp[672] + amp[674] + std::complex<double> (0, 1) *
      amp[675] - amp[676] - std::complex<double> (0, 1) * amp[677] - amp[678] -
      amp[684] - amp[685] - amp[687] - amp[688] + amp[690] + amp[692] +
      std::complex<double> (0, 1) * amp[693] - amp[694] - std::complex<double>
      (0, 1) * amp[695] - amp[696] - amp[702] - amp[703] - amp[705] - amp[706]
      + amp[780] - std::complex<double> (0, 1) * amp[781] + amp[782] - amp[784]
      - amp[786] + std::complex<double> (0, 1) * amp[787] - amp[792] - amp[793]
      - amp[795] - amp[796] + amp[798] - std::complex<double> (0, 1) * amp[799]
      + amp[800] - amp[802] - amp[804] + std::complex<double> (0, 1) * amp[805]
      - amp[810] - amp[811] - amp[813] - amp[814]);
  jamp[7] = +1./2. * (+1./3. * amp[192] + amp[193] + amp[194] + amp[195] +
      1./3. * amp[196] + 1./3. * amp[197] + 1./3. * amp[198] + amp[199] + 1./3.
      * amp[200] + amp[201] + amp[202] + amp[203] + 1./3. * amp[204] + 1./3. *
      amp[205] + 1./3. * amp[206] + amp[207] + amp[208] + amp[209] + 1./3. *
      amp[210] + 1./3. * amp[211] + 1./3. * amp[212] + 1./3. * amp[213] +
      amp[214] + amp[215] + 1./3. * amp[216] + 1./3. * amp[217] + 1./3. *
      amp[218] + 1./3. * amp[219] + amp[220] + amp[221] + amp[222] + amp[223] +
      std::complex<double> (0, 1) * amp[244] + amp[245] + std::complex<double>
      (0, 1) * amp[246] + 1./3. * amp[249] + 1./3. * amp[251] +
      std::complex<double> (0, 1) * amp[256] + amp[257] + std::complex<double>
      (0, 1) * amp[258] + 1./3. * amp[261] + 1./3. * amp[263] +
      std::complex<double> (0, 1) * amp[278] + amp[280] + std::complex<double>
      (0, 1) * amp[281] + std::complex<double> (0, 1) * amp[284] + amp[286] +
      std::complex<double> (0, 1) * amp[287] + amp[384] + amp[385] + 1./3. *
      amp[386] + 1./3. * amp[387] + 1./3. * amp[388] + 1./3. * amp[389] +
      amp[390] + amp[391] + 1./3. * amp[392] + 1./3. * amp[393] + 1./3. *
      amp[394] + 1./3. * amp[395] + amp[396] + amp[397] + amp[398] + amp[399] +
      1./3. * amp[416] + amp[417] + amp[418] + amp[419] + 1./3. * amp[420] +
      1./3. * amp[421] + 1./3. * amp[422] + amp[423] + 1./3. * amp[424] +
      amp[425] + amp[426] + amp[427] + 1./3. * amp[428] + 1./3. * amp[429] +
      1./3. * amp[430] + amp[431] - std::complex<double> (0, 1) * amp[436] -
      std::complex<double> (0, 1) * amp[438] + amp[439] + 1./3. * amp[441] +
      1./3. * amp[443] - std::complex<double> (0, 1) * amp[448] -
      std::complex<double> (0, 1) * amp[450] + amp[451] + 1./3. * amp[453] +
      1./3. * amp[455] - std::complex<double> (0, 1) * amp[470] + amp[472] -
      std::complex<double> (0, 1) * amp[473] - std::complex<double> (0, 1) *
      amp[476] + amp[478] - std::complex<double> (0, 1) * amp[479] + 1./3. *
      amp[485] + 1./3. * amp[487] - std::complex<double> (0, 1) * amp[488] +
      amp[489] - std::complex<double> (0, 1) * amp[490] + 1./3. * amp[497] +
      1./3. * amp[499] - std::complex<double> (0, 1) * amp[500] + amp[501] -
      std::complex<double> (0, 1) * amp[502] + std::complex<double> (0, 1) *
      amp[506] + std::complex<double> (0, 1) * amp[507] + amp[509] +
      std::complex<double> (0, 1) * amp[512] + std::complex<double> (0, 1) *
      amp[513] + amp[515] + 1./3. * amp[581] + 1./3. * amp[583] +
      std::complex<double> (0, 1) * amp[584] + std::complex<double> (0, 1) *
      amp[586] + amp[587] + 1./3. * amp[593] + 1./3. * amp[595] +
      std::complex<double> (0, 1) * amp[596] + std::complex<double> (0, 1) *
      amp[598] + amp[599] - std::complex<double> (0, 1) * amp[602] -
      std::complex<double> (0, 1) * amp[603] + amp[605] - std::complex<double>
      (0, 1) * amp[608] - std::complex<double> (0, 1) * amp[609] + amp[611] -
      amp[708] - amp[710] - std::complex<double> (0, 1) * amp[711] + amp[712] +
      std::complex<double> (0, 1) * amp[713] + amp[714] + amp[720] + amp[721] +
      amp[723] + amp[724] - amp[726] - amp[728] - std::complex<double> (0, 1) *
      amp[729] + amp[730] + std::complex<double> (0, 1) * amp[731] + amp[732] +
      amp[738] + amp[739] + amp[741] + amp[742] - amp[744] +
      std::complex<double> (0, 1) * amp[745] - amp[746] + amp[748] + amp[750] -
      std::complex<double> (0, 1) * amp[751] + amp[756] + amp[757] + amp[759] +
      amp[760] - amp[762] + std::complex<double> (0, 1) * amp[763] - amp[764] +
      amp[766] + amp[768] - std::complex<double> (0, 1) * amp[769] + amp[774] +
      amp[775] + amp[777] + amp[778]);
  jamp[8] = +1./2. * (+1./3. * std::complex<double> (0, 1) * amp[0] +
      std::complex<double> (0, 1) * amp[1] + std::complex<double> (0, 1) *
      amp[2] + 1./3. * std::complex<double> (0, 1) * amp[3] +
      std::complex<double> (0, 1) * amp[8] + 1./3. * std::complex<double> (0,
      1) * amp[9] + 1./3. * std::complex<double> (0, 1) * amp[10] +
      std::complex<double> (0, 1) * amp[11] + 1./3. * std::complex<double> (0,
      1) * amp[16] + std::complex<double> (0, 1) * amp[17] +
      std::complex<double> (0, 1) * amp[18] + 1./3. * std::complex<double> (0,
      1) * amp[19] + std::complex<double> (0, 1) * amp[24] + 1./3. *
      std::complex<double> (0, 1) * amp[25] + 1./3. * std::complex<double> (0,
      1) * amp[26] + std::complex<double> (0, 1) * amp[27] - amp[48] - amp[49]
      + std::complex<double> (0, 1) * amp[52] + std::complex<double> (0, 1) *
      amp[53] + std::complex<double> (0, 1) * amp[55] - amp[56] - amp[57] +
      std::complex<double> (0, 1) * amp[60] + std::complex<double> (0, 1) *
      amp[61] + std::complex<double> (0, 1) * amp[63] + std::complex<double>
      (0, 1) * amp[64] + std::complex<double> (0, 1) * amp[65] + amp[66] +
      amp[67] + std::complex<double> (0, 1) * amp[70] + std::complex<double>
      (0, 1) * amp[72] + std::complex<double> (0, 1) * amp[73] + amp[74] +
      amp[75] + std::complex<double> (0, 1) * amp[78] + 1./3. *
      std::complex<double> (0, 1) * amp[80] + 1./3. * std::complex<double> (0,
      1) * amp[81] + 1./3. * std::complex<double> (0, 1) * amp[84] + 1./3. *
      std::complex<double> (0, 1) * amp[85] + 1./3. * std::complex<double> (0,
      1) * amp[86] + 1./3. * std::complex<double> (0, 1) * amp[87] + 1./3. *
      std::complex<double> (0, 1) * amp[88] + 1./3. * std::complex<double> (0,
      1) * amp[89] + 1./3. * std::complex<double> (0, 1) * amp[92] + 1./3. *
      std::complex<double> (0, 1) * amp[93] + 1./3. * std::complex<double> (0,
      1) * amp[94] + 1./3. * std::complex<double> (0, 1) * amp[95] + 1./3. *
      amp[288] + 1./3. * amp[289] + amp[290] + amp[291] + amp[292] + amp[293] +
      1./3. * amp[294] + 1./3. * amp[295] + amp[296] + amp[297] + amp[298] +
      amp[299] + 1./3. * amp[300] + 1./3. * amp[301] + 1./3. * amp[302] + 1./3.
      * amp[303] + amp[336] + 1./3. * amp[337] + 1./3. * amp[338] + amp[339] +
      std::complex<double> (0, 1) * amp[344] + amp[345] + std::complex<double>
      (0, 1) * amp[346] + amp[348] + 1./3. * amp[349] + 1./3. * amp[350] +
      amp[351] + std::complex<double> (0, 1) * amp[356] + amp[357] +
      std::complex<double> (0, 1) * amp[358] + amp[360] + amp[361] +
      std::complex<double> (0, 1) * amp[362] + amp[363] + std::complex<double>
      (0, 1) * amp[365] + amp[366] + amp[367] + std::complex<double> (0, 1) *
      amp[368] + amp[369] + std::complex<double> (0, 1) * amp[371] + 1./3. *
      amp[372] + 1./3. * amp[373] + 1./3. * amp[375] + 1./3. * amp[376] + 1./3.
      * amp[378] + 1./3. * amp[379] + 1./3. * amp[381] + 1./3. * amp[382] +
      1./3. * amp[480] + amp[481] + amp[482] + 1./3. * amp[483] +
      std::complex<double> (0, 1) * amp[488] + std::complex<double> (0, 1) *
      amp[490] + amp[491] + 1./3. * amp[492] + amp[493] + amp[494] + 1./3. *
      amp[495] + std::complex<double> (0, 1) * amp[500] + std::complex<double>
      (0, 1) * amp[502] + amp[503] + amp[504] + amp[505] - std::complex<double>
      (0, 1) * amp[506] - std::complex<double> (0, 1) * amp[507] + amp[508] +
      amp[510] + amp[511] - std::complex<double> (0, 1) * amp[512] -
      std::complex<double> (0, 1) * amp[513] + amp[514] + 1./3. * amp[516] +
      1./3. * amp[517] + 1./3. * amp[520] + 1./3. * amp[521] + 1./3. * amp[522]
      + 1./3. * amp[523] + 1./3. * amp[526] + 1./3. * amp[527] - amp[712] -
      amp[714] + std::complex<double> (0, 1) * amp[715] + amp[718] + amp[722] -
      amp[720] + amp[725] - amp[723] - amp[730] - amp[732] +
      std::complex<double> (0, 1) * amp[733] + amp[736] + amp[740] - amp[738] +
      amp[743] - amp[741] + amp[744] - std::complex<double> (0, 1) * amp[745] +
      amp[746] + amp[753] - amp[758] - amp[757] - amp[761] - amp[760] +
      amp[762] - std::complex<double> (0, 1) * amp[763] + amp[764] + amp[771] -
      amp[776] - amp[775] - amp[779] - amp[778] + 1./3. * amp[789] + 1./3. *
      amp[790] + 1./3. * amp[807] + 1./3. * amp[808]);
  jamp[9] = +1./2. * (-std::complex<double> (0, 1) * amp[0] - 1./3. *
      std::complex<double> (0, 1) * amp[1] - 1./3. * std::complex<double> (0,
      1) * amp[2] - std::complex<double> (0, 1) * amp[3] - std::complex<double>
      (0, 1) * amp[12] - 1./3. * std::complex<double> (0, 1) * amp[13] - 1./3.
      * std::complex<double> (0, 1) * amp[14] - std::complex<double> (0, 1) *
      amp[15] - std::complex<double> (0, 1) * amp[16] - 1./3. *
      std::complex<double> (0, 1) * amp[17] - 1./3. * std::complex<double> (0,
      1) * amp[18] - std::complex<double> (0, 1) * amp[19] -
      std::complex<double> (0, 1) * amp[28] - 1./3. * std::complex<double> (0,
      1) * amp[29] - 1./3. * std::complex<double> (0, 1) * amp[30] -
      std::complex<double> (0, 1) * amp[31] + amp[32] + amp[33] -
      std::complex<double> (0, 1) * amp[36] - std::complex<double> (0, 1) *
      amp[37] - std::complex<double> (0, 1) * amp[39] + amp[40] + amp[41] -
      std::complex<double> (0, 1) * amp[44] - std::complex<double> (0, 1) *
      amp[45] - std::complex<double> (0, 1) * amp[47] - 1./3. *
      std::complex<double> (0, 1) * amp[64] - 1./3. * std::complex<double> (0,
      1) * amp[65] - 1./3. * std::complex<double> (0, 1) * amp[68] - 1./3. *
      std::complex<double> (0, 1) * amp[69] - 1./3. * std::complex<double> (0,
      1) * amp[70] - 1./3. * std::complex<double> (0, 1) * amp[71] - 1./3. *
      std::complex<double> (0, 1) * amp[72] - 1./3. * std::complex<double> (0,
      1) * amp[73] - 1./3. * std::complex<double> (0, 1) * amp[76] - 1./3. *
      std::complex<double> (0, 1) * amp[77] - 1./3. * std::complex<double> (0,
      1) * amp[78] - 1./3. * std::complex<double> (0, 1) * amp[79] -
      std::complex<double> (0, 1) * amp[80] - std::complex<double> (0, 1) *
      amp[81] - amp[82] - amp[83] - std::complex<double> (0, 1) * amp[86] -
      std::complex<double> (0, 1) * amp[88] - std::complex<double> (0, 1) *
      amp[89] - amp[90] - amp[91] - std::complex<double> (0, 1) * amp[94] -
      1./3. * amp[384] - 1./3. * amp[385] - amp[386] - amp[387] - amp[388] -
      amp[389] - 1./3. * amp[390] - 1./3. * amp[391] - amp[392] - amp[393] -
      amp[394] - amp[395] - 1./3. * amp[396] - 1./3. * amp[397] - 1./3. *
      amp[398] - 1./3. * amp[399] - amp[432] - 1./3. * amp[433] - 1./3. *
      amp[434] - amp[435] - std::complex<double> (0, 1) * amp[440] - amp[441] -
      std::complex<double> (0, 1) * amp[442] - amp[444] - 1./3. * amp[445] -
      1./3. * amp[446] - amp[447] - std::complex<double> (0, 1) * amp[452] -
      amp[453] - std::complex<double> (0, 1) * amp[454] - amp[456] - amp[457] -
      std::complex<double> (0, 1) * amp[458] - amp[459] - std::complex<double>
      (0, 1) * amp[461] - amp[462] - amp[463] - std::complex<double> (0, 1) *
      amp[464] - amp[465] - std::complex<double> (0, 1) * amp[467] - 1./3. *
      amp[468] - 1./3. * amp[469] - 1./3. * amp[471] - 1./3. * amp[472] - 1./3.
      * amp[474] - 1./3. * amp[475] - 1./3. * amp[477] - 1./3. * amp[478] -
      amp[480] - 1./3. * amp[481] - 1./3. * amp[482] - amp[483] -
      std::complex<double> (0, 1) * amp[484] - std::complex<double> (0, 1) *
      amp[486] - amp[487] - amp[492] - 1./3. * amp[493] - 1./3. * amp[494] -
      amp[495] - std::complex<double> (0, 1) * amp[496] - std::complex<double>
      (0, 1) * amp[498] - amp[499] - 1./3. * amp[504] - 1./3. * amp[505] -
      1./3. * amp[508] - 1./3. * amp[509] - 1./3. * amp[510] - 1./3. * amp[511]
      - 1./3. * amp[514] - 1./3. * amp[515] - amp[516] - amp[517] +
      std::complex<double> (0, 1) * amp[518] + std::complex<double> (0, 1) *
      amp[519] - amp[520] - amp[522] - amp[523] + std::complex<double> (0, 1) *
      amp[524] + std::complex<double> (0, 1) * amp[525] - amp[526] + amp[676] +
      amp[678] - std::complex<double> (0, 1) * amp[679] - amp[682] - amp[686] +
      amp[684] - amp[689] + amp[687] + amp[694] + amp[696] -
      std::complex<double> (0, 1) * amp[697] - amp[700] - amp[704] + amp[702] -
      amp[707] + amp[705] - 1./3. * amp[753] - 1./3. * amp[754] - 1./3. *
      amp[771] - 1./3. * amp[772] - amp[780] + std::complex<double> (0, 1) *
      amp[781] - amp[782] - amp[789] + amp[794] + amp[793] + amp[797] +
      amp[796] - amp[798] + std::complex<double> (0, 1) * amp[799] - amp[800] -
      amp[807] + amp[812] + amp[811] + amp[815] + amp[814]);
  jamp[10] = +1./2. * (-1./3. * std::complex<double> (0, 1) * amp[4] -
      std::complex<double> (0, 1) * amp[5] - std::complex<double> (0, 1) *
      amp[6] - 1./3. * std::complex<double> (0, 1) * amp[7] - 1./3. *
      std::complex<double> (0, 1) * amp[8] - std::complex<double> (0, 1) *
      amp[9] - std::complex<double> (0, 1) * amp[10] - 1./3. *
      std::complex<double> (0, 1) * amp[11] - 1./3. * std::complex<double> (0,
      1) * amp[20] - std::complex<double> (0, 1) * amp[21] -
      std::complex<double> (0, 1) * amp[22] - 1./3. * std::complex<double> (0,
      1) * amp[23] - 1./3. * std::complex<double> (0, 1) * amp[24] -
      std::complex<double> (0, 1) * amp[25] - std::complex<double> (0, 1) *
      amp[26] - 1./3. * std::complex<double> (0, 1) * amp[27] - amp[32] -
      amp[33] - std::complex<double> (0, 1) * amp[34] - std::complex<double>
      (0, 1) * amp[35] - std::complex<double> (0, 1) * amp[38] - amp[40] -
      amp[41] - std::complex<double> (0, 1) * amp[42] - std::complex<double>
      (0, 1) * amp[43] - std::complex<double> (0, 1) * amp[46] - 1./3. *
      std::complex<double> (0, 1) * amp[50] - 1./3. * std::complex<double> (0,
      1) * amp[51] - 1./3. * std::complex<double> (0, 1) * amp[52] - 1./3. *
      std::complex<double> (0, 1) * amp[53] - 1./3. * std::complex<double> (0,
      1) * amp[54] - 1./3. * std::complex<double> (0, 1) * amp[55] - 1./3. *
      std::complex<double> (0, 1) * amp[58] - 1./3. * std::complex<double> (0,
      1) * amp[59] - 1./3. * std::complex<double> (0, 1) * amp[60] - 1./3. *
      std::complex<double> (0, 1) * amp[61] - 1./3. * std::complex<double> (0,
      1) * amp[62] - 1./3. * std::complex<double> (0, 1) * amp[63] + amp[82] +
      amp[83] - std::complex<double> (0, 1) * amp[84] - std::complex<double>
      (0, 1) * amp[85] - std::complex<double> (0, 1) * amp[87] + amp[90] +
      amp[91] - std::complex<double> (0, 1) * amp[92] - std::complex<double>
      (0, 1) * amp[93] - std::complex<double> (0, 1) * amp[95] - 1./3. *
      amp[304] - 1./3. * amp[305] - amp[306] - amp[307] - amp[308] - amp[309] -
      1./3. * amp[310] - 1./3. * amp[311] - amp[312] - amp[313] - amp[314] -
      amp[315] - 1./3. * amp[316] - 1./3. * amp[317] - 1./3. * amp[318] - 1./3.
      * amp[319] - 1./3. * amp[336] - amp[337] - amp[338] - 1./3. * amp[339] -
      std::complex<double> (0, 1) * amp[340] - amp[341] - std::complex<double>
      (0, 1) * amp[342] - 1./3. * amp[348] - amp[349] - amp[350] - 1./3. *
      amp[351] - std::complex<double> (0, 1) * amp[352] - amp[353] -
      std::complex<double> (0, 1) * amp[354] - 1./3. * amp[360] - 1./3. *
      amp[361] - 1./3. * amp[363] - 1./3. * amp[364] - 1./3. * amp[366] - 1./3.
      * amp[367] - 1./3. * amp[369] - 1./3. * amp[370] - amp[372] - amp[373] -
      std::complex<double> (0, 1) * amp[374] - amp[375] - std::complex<double>
      (0, 1) * amp[377] - amp[378] - amp[379] - std::complex<double> (0, 1) *
      amp[380] - amp[381] - std::complex<double> (0, 1) * amp[383] - 1./3. *
      amp[528] - amp[529] - amp[530] - 1./3. * amp[531] - std::complex<double>
      (0, 1) * amp[536] - std::complex<double> (0, 1) * amp[538] - amp[539] -
      1./3. * amp[540] - amp[541] - amp[542] - 1./3. * amp[543] -
      std::complex<double> (0, 1) * amp[548] - std::complex<double> (0, 1) *
      amp[550] - amp[551] - amp[552] - amp[553] + std::complex<double> (0, 1) *
      amp[554] + std::complex<double> (0, 1) * amp[555] - amp[556] - amp[558] -
      amp[559] + std::complex<double> (0, 1) * amp[560] + std::complex<double>
      (0, 1) * amp[561] - amp[562] - 1./3. * amp[564] - 1./3. * amp[565] -
      1./3. * amp[568] - 1./3. * amp[569] - 1./3. * amp[570] - 1./3. * amp[571]
      - 1./3. * amp[574] - 1./3. * amp[575] - amp[672] + std::complex<double>
      (0, 1) * amp[673] - amp[674] - amp[681] + amp[686] + amp[685] + amp[689]
      + amp[688] - amp[690] + std::complex<double> (0, 1) * amp[691] - amp[692]
      - amp[699] + amp[704] + amp[703] + amp[707] + amp[706] - 1./3. * amp[717]
      - 1./3. * amp[718] - 1./3. * amp[735] - 1./3. * amp[736] + amp[784] +
      amp[786] - std::complex<double> (0, 1) * amp[787] - amp[790] - amp[794] +
      amp[792] - amp[797] + amp[795] + amp[802] + amp[804] -
      std::complex<double> (0, 1) * amp[805] - amp[808] - amp[812] + amp[810] -
      amp[815] + amp[813]);
  jamp[11] = +1./2. * (+std::complex<double> (0, 1) * amp[4] + 1./3. *
      std::complex<double> (0, 1) * amp[5] + 1./3. * std::complex<double> (0,
      1) * amp[6] + std::complex<double> (0, 1) * amp[7] + 1./3. *
      std::complex<double> (0, 1) * amp[12] + std::complex<double> (0, 1) *
      amp[13] + std::complex<double> (0, 1) * amp[14] + 1./3. *
      std::complex<double> (0, 1) * amp[15] + std::complex<double> (0, 1) *
      amp[20] + 1./3. * std::complex<double> (0, 1) * amp[21] + 1./3. *
      std::complex<double> (0, 1) * amp[22] + std::complex<double> (0, 1) *
      amp[23] + 1./3. * std::complex<double> (0, 1) * amp[28] +
      std::complex<double> (0, 1) * amp[29] + std::complex<double> (0, 1) *
      amp[30] + 1./3. * std::complex<double> (0, 1) * amp[31] + 1./3. *
      std::complex<double> (0, 1) * amp[34] + 1./3. * std::complex<double> (0,
      1) * amp[35] + 1./3. * std::complex<double> (0, 1) * amp[36] + 1./3. *
      std::complex<double> (0, 1) * amp[37] + 1./3. * std::complex<double> (0,
      1) * amp[38] + 1./3. * std::complex<double> (0, 1) * amp[39] + 1./3. *
      std::complex<double> (0, 1) * amp[42] + 1./3. * std::complex<double> (0,
      1) * amp[43] + 1./3. * std::complex<double> (0, 1) * amp[44] + 1./3. *
      std::complex<double> (0, 1) * amp[45] + 1./3. * std::complex<double> (0,
      1) * amp[46] + 1./3. * std::complex<double> (0, 1) * amp[47] + amp[48] +
      amp[49] + std::complex<double> (0, 1) * amp[50] + std::complex<double>
      (0, 1) * amp[51] + std::complex<double> (0, 1) * amp[54] + amp[56] +
      amp[57] + std::complex<double> (0, 1) * amp[58] + std::complex<double>
      (0, 1) * amp[59] + std::complex<double> (0, 1) * amp[62] - amp[66] -
      amp[67] + std::complex<double> (0, 1) * amp[68] + std::complex<double>
      (0, 1) * amp[69] + std::complex<double> (0, 1) * amp[71] - amp[74] -
      amp[75] + std::complex<double> (0, 1) * amp[76] + std::complex<double>
      (0, 1) * amp[77] + std::complex<double> (0, 1) * amp[79] + 1./3. *
      amp[400] + 1./3. * amp[401] + amp[402] + amp[403] + amp[404] + amp[405] +
      1./3. * amp[406] + 1./3. * amp[407] + amp[408] + amp[409] + amp[410] +
      amp[411] + 1./3. * amp[412] + 1./3. * amp[413] + 1./3. * amp[414] + 1./3.
      * amp[415] + 1./3. * amp[432] + amp[433] + amp[434] + 1./3. * amp[435] +
      std::complex<double> (0, 1) * amp[436] + amp[437] + std::complex<double>
      (0, 1) * amp[438] + 1./3. * amp[444] + amp[445] + amp[446] + 1./3. *
      amp[447] + std::complex<double> (0, 1) * amp[448] + amp[449] +
      std::complex<double> (0, 1) * amp[450] + 1./3. * amp[456] + 1./3. *
      amp[457] + 1./3. * amp[459] + 1./3. * amp[460] + 1./3. * amp[462] + 1./3.
      * amp[463] + 1./3. * amp[465] + 1./3. * amp[466] + amp[468] + amp[469] +
      std::complex<double> (0, 1) * amp[470] + amp[471] + std::complex<double>
      (0, 1) * amp[473] + amp[474] + amp[475] + std::complex<double> (0, 1) *
      amp[476] + amp[477] + std::complex<double> (0, 1) * amp[479] + amp[528] +
      1./3. * amp[529] + 1./3. * amp[530] + amp[531] + std::complex<double> (0,
      1) * amp[532] + std::complex<double> (0, 1) * amp[534] + amp[535] +
      amp[540] + 1./3. * amp[541] + 1./3. * amp[542] + amp[543] +
      std::complex<double> (0, 1) * amp[544] + std::complex<double> (0, 1) *
      amp[546] + amp[547] + 1./3. * amp[552] + 1./3. * amp[553] + 1./3. *
      amp[556] + 1./3. * amp[557] + 1./3. * amp[558] + 1./3. * amp[559] + 1./3.
      * amp[562] + 1./3. * amp[563] + amp[564] + amp[565] -
      std::complex<double> (0, 1) * amp[566] - std::complex<double> (0, 1) *
      amp[567] + amp[568] + amp[570] + amp[571] - std::complex<double> (0, 1) *
      amp[572] - std::complex<double> (0, 1) * amp[573] + amp[574] + 1./3. *
      amp[681] + 1./3. * amp[682] + 1./3. * amp[699] + 1./3. * amp[700] +
      amp[708] - std::complex<double> (0, 1) * amp[709] + amp[710] + amp[717] -
      amp[722] - amp[721] - amp[725] - amp[724] + amp[726] -
      std::complex<double> (0, 1) * amp[727] + amp[728] + amp[735] - amp[740] -
      amp[739] - amp[743] - amp[742] - amp[748] - amp[750] +
      std::complex<double> (0, 1) * amp[751] + amp[754] + amp[758] - amp[756] +
      amp[761] - amp[759] - amp[766] - amp[768] + std::complex<double> (0, 1) *
      amp[769] + amp[772] + amp[776] - amp[774] + amp[779] - amp[777]);

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



