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
#include "timer.h"

#include <thrust/complex.h>

#define TIMERTYPE std::chrono::high_resolution_clock

//==========================================================================
// A class for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1
//--------------------------------------------------------------------------

struct processMem {
  static const int tnwavefuncs = 6;
  static const int twrows = 18;
  static const int tnamplitudes = 2;
  thrust::complex<double> *tamp; //[tnamplitudes];
  thrust::complex<double> **tw;  //[twrows][tnwavefuncs];
};

class CPPProcess {
public:
  // Constructor.

  CPPProcess();
  CPPProcess(processMem *pm, bool verbose=false, bool debug=false);

  ~CPPProcess();

  void resetGPUMemory();

  // Initialize process.
  virtual void initProc(std::string param_card_name);

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Info on the subprocess.
  virtual std::string name() const { return "e+ e- > mu+ mu- (sm)"; }

  virtual int code() const { return 1; }

  // thrust::host_vector<double> getMasses();
  const std::vector<double> &getMasses() const;

  // Get and set momenta for matrix element evaluation
  std::vector<double *> getMomenta() { return p; }                                                                      

  void setMomenta(std::vector<double *> &momenta); //{ p = momenta; }                                                     
  //void setMomenta(std::vector<double *> &momenta);

  void setInitial(int inid1, int inid2) {
    id1 = inid1;
    id2 = inid2;
  }

  // Get matrix element vector
  const double *getMatrixElements() const { return matrix_element; }

  // Constants for array limits
  static const int ninitial = 2;
  static const int nexternal = 4;
  static const int nprocesses = 1;

private:

  Timer<TIMERTYPE> m_timer;

  // print verbose info
  bool m_verbose;
  
  // print debug info
  bool m_debug;
  
  // Private functions to calculate the matrix element for all subprocesses
  // Calculate wavefunctions
  void calculate_wavefunctions(const int perm[], const int hel[]);

  processMem *m;

  static const int nwavefuncs = 6;

  thrust::complex<double> w[nwavefuncs][18];

  static const int namplitudes = 2;

  thrust::complex<double> amp[namplitudes];

  double matrix_1_epem_mupmum();

  // Store the matrix element value from sigmaKin
  double matrix_element[nprocesses];

  // Color flows, used when selecting color
  double *jamp2[nprocesses];

  // Pointer to the model parameters
  Parameters_sm *pars;

  // vector with external particle masses
  std::vector<double> mME;

  // vector with momenta (to be changed each event)
   std::vector<double *> p;

  // Initial particle ids
  int id1, id2;
};

#endif // MG5_Sigma_sm_epem_mupmum_H
