//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.0, 2020-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"
#include "HelAmps_sm.h"
#include <algorithm> // perf stats
#include <iostream>
#include <numeric> // perf stats

using namespace MG5_sm;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1

//--------------------------------------------------------------------------
// Initialize process.

void CPPProcess::initProc(string param_card_name, bool verb) {
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm::getInstance();
  SLHAReader slha(param_card_name, verb);
  pars->setIndependentParameters(slha);
  pars->setIndependentCouplings();
  if (verb) {
    pars->printIndependentParameters();
    pars->printIndependentCouplings();
  }
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO);
  mME.push_back(pars->ZERO);
  mME.push_back(pars->ZERO);
  mME.push_back(pars->ZERO);
  jamp2[0] = new double[1];
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void CPPProcess::sigmaKin(bool ppar) {
  // Set the parameters which change event by event
  pars->setDependentParameters();
  pars->setDependentCouplings();
  static bool firsttime = true;
  if (firsttime && ppar) {
    pars->printDependentParameters();
    pars->printDependentCouplings();
    firsttime = false;
  }

  m_timer.Start();

  // Reset color flows
  for (int i = 0; i < 1; i++)
    jamp2[0][i] = 0.;

  // Local variables and constants
  const int ncomb = 16;
  static bool goodhel[ncomb] = {ncomb * false};
  static int ntry = 0, sum_hel = 0, ngood = 0;
  static int igood[ncomb];
  static int jhel;
  std::complex<double> **wfs;
  double t[nprocesses];
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {
      {-1, -1, -1, -1}, {-1, -1, -1, 1}, {-1, -1, 1, -1}, {-1, -1, 1, 1},
      {-1, 1, -1, -1},  {-1, 1, -1, 1},  {-1, 1, 1, -1},  {-1, 1, 1, 1},
      {1, -1, -1, -1},  {1, -1, -1, 1},  {1, -1, 1, -1},  {1, -1, 1, 1},
      {1, 1, -1, -1},   {1, 1, -1, 1},   {1, 1, 1, -1},   {1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {4};

  ntry = ntry + 1;

  // Reset the matrix elements
  for (int i = 0; i < nprocesses; i++) {
    matrix_element[i] = 0.;
  }
  // Define permutation
  int perm[nexternal];
  for (int i = 0; i < nexternal; i++) {
    perm[i] = i;
  }

  if (sum_hel == 0 || ntry < 10) {
    // Calculate the matrix element for all helicities
    for (int ihel = 0; ihel < ncomb; ihel++) {
      if (goodhel[ihel] || ntry < 2) {

        /*
        std::cout << std::endl
                  << std::endl
                  << "<<<<< " << ihel << " " << ihel << " " << ihel << " "
                  << ihel << " " << ihel << " " << ihel << " " << ihel << " "
                  << ihel << " " << ihel << " " << ihel << " " << ihel << " "
                  << ihel << " " << ihel << " " << ihel << " " << ihel << " "
                  << ihel << " "
                  << " >>>>>>>>" << std::endl;
                  */

        calculate_wavefunctions(perm, helicities[ihel]);
        t[0] = matrix_1_epem_mupmum();
        // std::cout << "t[0] oben: " << t[0] << std::endl;

        double tsum = 0;
        for (int iproc = 0; iproc < nprocesses; iproc++) {
          matrix_element[iproc] += t[iproc];
          tsum += t[iproc];
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel]) {
          goodhel[ihel] = true;
          ngood++;
          igood[ngood] = ihel;
        }
      }
    }
    jhel = 0;
    sum_hel = min(sum_hel, ngood);
  } else {
    // Only use the "good" helicities
    for (int j = 0; j < sum_hel; j++) {
      jhel++;
      if (jhel >= ngood)
        jhel = 0;
      double hwgt = double(ngood) / double(sum_hel);
      int ihel = igood[jhel];
      calculate_wavefunctions(perm, helicities[ihel]);
      t[0] = matrix_1_epem_mupmum();
      // std::cout << "t[0] unten: " << t[0] << std::endl;

      for (int iproc = 0; iproc < nprocesses; iproc++) {
        matrix_element[iproc] += t[iproc] * hwgt;
      }
    }
  }

  //if ( ntry == 10 ) for (int ihel = 0; ihel < ncomb; ihel++ ) printf( "sigmakin: ihel %2d %d\n", ihel, goodhel[ihel] );

  for (int i = 0; i < nprocesses; i++)
    matrix_element[i] /= denominators[i];

  float gputime = m_timer.GetDuration();
  m_wavetimes.push_back(gputime);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double CPPProcess::sigmaHat() {
  // Select between the different processes
  if (id1 == -11 && id2 == 11) {
    // Add matrix elements for processes with beams (-11, 11)
    return matrix_element[0];
  } else {
    // Return 0 if not correct initial state assignment
    return 0.;
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void CPPProcess::calculate_wavefunctions(const int perm[], const int hel[]) {
  // Calculate wavefunctions for all processes
  int i, j;

  /*
  std::cout << "<<< w: " << std::endl;
  for (int i = 0; i < 6; ++i) {
    std::cout << "w" << i << ": ";
    for (int j = 0; j < 18; ++j) {
      if (w[i][j].real() || w[i][j].imag())
        std::cout << w[i][j] << " ";
      else
        std::cout << "0 ";
    }
    std::cout << std::endl;
  }
  */

  // Calculate all wavefunctions
  oxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
  ixxxxx(p[perm[1]], mME[1], hel[1], +1, w[1]);
  ixxxxx(p[perm[2]], mME[2], hel[2], -1, w[2]);
  oxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]);
  FFV1P0_3(w[1], w[0], pars->GC_3, pars->ZERO, pars->ZERO, w[4]);
  FFV2_4_3(w[1], w[0], -pars->GC_51, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ,
           w[5]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[2], w[3], w[4], pars->GC_3, amp[0]);
  FFV2_4_0(w[2], w[3], w[5], -pars->GC_51, pars->GC_59, amp[1]);

  // std::cout << "Wave function time: " << gputime << std::endl;

  /*
  std::cout << ">>> w: " << std::endl;
  for (int i = 0; i < 6; ++i) {
    std::cout << "w" << i << ": ";
    for (int j = 0; j < 18; ++j) {
      if (w[i][j].real() || w[i][j].imag())
        std::cout << w[i][j] << " ";
      else
        std::cout << "0 ";
    }
    std::cout << std::endl;
  }

  std::cout << ">>>>>>>> tamp: ";
  for (int x = 0; x < namplitudes; ++x) {
    std::cout << amp[x] << " ";
  }
  std::cout << std::endl;
  */
}

double CPPProcess::matrix_1_epem_mupmum() {
  int i, j;
  // Local variables
  const int ngraphs = 2;
  const int ncolor = 1;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  // The color matrix;
  static const double denom[ncolor] = {1};
  static const double cf[ncolor][ncolor] = {{1}};

  // Calculate color flows
  jamp[0] = -amp[0] - amp[1];

  // Sum and square the color flows to get the matrix element
  double matrix = 0;
  for (i = 0; i < ncolor; i++) {
    ztemp = 0.;
    for (j = 0; j < ncolor; j++)
      ztemp = ztemp + cf[i][j] * jamp[j];
    matrix = matrix + real(ztemp * conj(jamp[i])) / denom[i];
  }

  // Store the leading color flows for choice of color
  for (i = 0; i < ncolor; i++)
    jamp2[0][i] += real(jamp[i] * conj(jamp[i]));

  return matrix;
}

void CPPProcess::printPerformanceStats() {
  float sum = std::accumulate(m_wavetimes.begin(), m_wavetimes.end(), 0.0);
  int numelems = m_wavetimes.size();
  float mean = sum / numelems;
  float sq_sum = std::inner_product(m_wavetimes.begin(), m_wavetimes.end(),
                                    m_wavetimes.begin(), 0.0);
  float stdev = std::sqrt(sq_sum / numelems - mean * mean);
  std::vector<float>::iterator mintime =
      std::min_element(m_wavetimes.begin(), m_wavetimes.end());
  std::vector<float>::iterator maxtime =
      std::max_element(m_wavetimes.begin(), m_wavetimes.end());

  std::cout << "***********************************" << std::endl
            << "NumberOfEntries       = " << numelems << std::endl
            << std::scientific 
            << "TotalTimeInWaveFuncs  = " << sum << " sec" << std::endl
            << "MeanTimeInWaveFuncs   = " << mean << " sec" << std::endl
            << "StdDevTimeInWaveFuncs = " << stdev << " sec" << std::endl
            << "MinTimeInWaveFuncs    = " << *mintime << " sec" << std::endl
            << "MaxTimeInWaveFuncs    = " << *maxtime << " sec" << std::endl
            << "-----------------------------------" << std::endl
            << "NumMatrixElements     = " << numelems << std::endl
            << "MatrixElementsPerSec  = " << 1/mean << " sec^-1" << std::endl;  
}
