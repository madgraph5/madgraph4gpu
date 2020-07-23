//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.0, 2020-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_epem_mupmum_H
#define MG5_Sigma_sm_epem_mupmum_H

#include <complex>
#include <vector>

#include "Parameters_sm.h"

#include <complex>

//==========================================================================
// A class for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1
//--------------------------------------------------------------------------

class CPPProcess {
public:
  // Constructor.

  CPPProcess(int numiterations, int groups, int items,
             bool verbose = false, bool debug = false);

  ~CPPProcess();

  // Initialize process.
  virtual void initProc(std::string param_card_name);

  // everything from sigmaKin which has to be left outside the kernel
  void preSigmaKin();

  // Info on the subprocess.
  virtual std::string name() const { return "e+ e- > mu+ mu- (sm)"; }

  virtual int code() const { return 1; }

  const std::vector<double> &getMasses() const;

  void setInitial(int inid1, int inid2) {
    id1 = inid1;
    id2 = inid2;
  }

  int getDim() const { return dim; }

  int getNIOParticles() const { return nioparticles; }

  // Constants for array limits
  static const int ninitial = 2;
  static const int nexternal = 4;
  static const int nprocesses = 1;

  std::complex<double> IPC[3];
  double IPD[2];

private:
  int m_numiterations;
  // gpu variables
  int work_groups;
  int work_items;
  int dim; // work_groups * work_items;

  // print verbose info
  bool m_verbose;

  // print debug info
  bool m_debug;

  static const int nwavefuncs = 6;

  static const int namplitudes = 2;

  static const int ncomb = 16;

  static const int wrows = 6; // was 18;

  static const int nioparticles = 4;

  std::complex<double> **amp; // [dim][namplitudes];

  // Pointer to the model parameters
  Parameters_sm *pars;

  // vector with external particle masses
  std::vector<double> mME;

  // Initial particle ids
  int id1, id2;
};

#endif // MG5_Sigma_sm_epem_mupmum_H
