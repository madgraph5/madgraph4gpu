#include "CPPProcess.h"
#include "HelAmps_sm.h"

#include <algorithm>
#include <iostream>
#include <complex>

CPPProcess::CPPProcess(int numiterations, int groups, int items,
                       bool verbose, bool debug)
    : m_numiterations(numiterations), work_groups(groups),
      work_items(items), m_verbose(verbose), m_debug(debug),
      dim(groups * items), mME(4, 0.00) {

  static double tmpmME[4] = {0.0, 0.0, 0.0, 0.0};

  // Helicities for the process - nodim
  static int tHel[ncomb][nexternal] = {
      {-1, -1, -1, -1}, {-1, -1, -1, 1}, {-1, -1, 1, -1}, {-1, -1, 1, 1},
      {-1, 1, -1, -1},  {-1, 1, -1, 1},  {-1, 1, 1, -1},  {-1, 1, 1, 1},
      {1, -1, -1, -1},  {1, -1, -1, 1},  {1, -1, 1, -1},  {1, -1, 1, 1},
      {1, 1, -1, -1},   {1, 1, -1, 1},   {1, 1, 1, -1},   {1, 1, 1, 1}};

  // perm - nodim
  static int perm[nexternal] = {0, 1, 2, 3};
}

CPPProcess::~CPPProcess() {}

const std::vector<double> &CPPProcess::getMasses() const { return mME; }

void CPPProcess::initProc(std::string param_card_name) {

  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm::getInstance();
  SLHAReader slha(param_card_name, m_verbose);
  pars->setIndependentParameters(slha);
  pars->setIndependentCouplings();
  IPC[0] = pars->GC_3;
  IPC[1] = pars->GC_51;
  IPC[2] = pars->GC_59;
  IPD[0] = pars->mdl_MZ;
  IPD[1] = pars->mdl_WZ;
  if (m_verbose) {
    pars->printIndependentParameters();
    pars->printIndependentCouplings();
  }
  // jamp2[0] = new double[1];
}

void CPPProcess::preSigmaKin() {
  // Set the parameters which change event by event
  pars->setDependentParameters();
  pars->setDependentCouplings();
  static bool firsttime = true;
  if (firsttime && m_verbose) {
    pars->printDependentParameters();
    pars->printDependentCouplings();
    firsttime = false;
  }
}
