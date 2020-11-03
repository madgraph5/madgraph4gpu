//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.0, 2020-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_ckm_no_b_mass_uxdx_epemsuxdxsx_H
#define MG5_Sigma_sm_ckm_no_b_mass_uxdx_epemsuxdxsx_H

#include <complex> 
#include <vector> 

#include "Parameters_sm_ckm_no_b_mass.h"

using namespace std; 

//==========================================================================
// A class for calculating the matrix elements for
// Process: u~ d~ > e+ e- s u~ d~ s~ WEIGHTED<=8 / h
// Process: u~ d~ > e+ e- d~ b u~ b~ WEIGHTED<=8 / h
// Process: u~ d~ > mu+ mu- s u~ d~ s~ WEIGHTED<=8 / h
// Process: u~ d~ > mu+ mu- d~ b u~ b~ WEIGHTED<=8 / h
// Process: u~ s~ > e+ e- d u~ s~ d~ WEIGHTED<=8 / h
// Process: u~ s~ > e+ e- s~ b u~ b~ WEIGHTED<=8 / h
// Process: u~ s~ > mu+ mu- d u~ s~ d~ WEIGHTED<=8 / h
// Process: u~ s~ > mu+ mu- s~ b u~ b~ WEIGHTED<=8 / h
// Process: u~ b~ > e+ e- d u~ b~ d~ WEIGHTED<=8 / h
// Process: u~ b~ > e+ e- s u~ b~ s~ WEIGHTED<=8 / h
// Process: u~ b~ > mu+ mu- d u~ b~ d~ WEIGHTED<=8 / h
// Process: u~ b~ > mu+ mu- s u~ b~ s~ WEIGHTED<=8 / h
// Process: c~ d~ > e+ e- s c~ d~ s~ WEIGHTED<=8 / h
// Process: c~ d~ > e+ e- d~ b c~ b~ WEIGHTED<=8 / h
// Process: c~ d~ > mu+ mu- s c~ d~ s~ WEIGHTED<=8 / h
// Process: c~ d~ > mu+ mu- d~ b c~ b~ WEIGHTED<=8 / h
// Process: c~ s~ > e+ e- d c~ s~ d~ WEIGHTED<=8 / h
// Process: c~ s~ > e+ e- s~ b c~ b~ WEIGHTED<=8 / h
// Process: c~ s~ > mu+ mu- d c~ s~ d~ WEIGHTED<=8 / h
// Process: c~ s~ > mu+ mu- s~ b c~ b~ WEIGHTED<=8 / h
// Process: c~ b~ > e+ e- d c~ b~ d~ WEIGHTED<=8 / h
// Process: c~ b~ > e+ e- s c~ b~ s~ WEIGHTED<=8 / h
// Process: c~ b~ > mu+ mu- d c~ b~ d~ WEIGHTED<=8 / h
// Process: c~ b~ > mu+ mu- s c~ b~ s~ WEIGHTED<=8 / h
//--------------------------------------------------------------------------

class CPPProcess
{
  public:

    // Constructor.
    CPPProcess() {}

    // Initialize process.
    virtual void initProc(string param_card_name); 

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin(); 

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat(); 

    // Info on the subprocess.
    virtual string name() const {return "u~ d~ > e+ e- s u~ d~ s~ (sm_ckm_no_b_mass)";}

    virtual int code() const {return 0;}

    const vector<double> & getMasses() const {return mME;}

    // Get and set momenta for matrix element evaluation
    vector < double * > getMomenta(){return p;}
    void setMomenta(vector < double * > & momenta){p = momenta;}
    void setInitial(int inid1, int inid2){id1 = inid1; id2 = inid2;}

    // Get matrix element vector
    const double * getMatrixElements() const {return matrix_element;}

    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 8; 
    static const int nprocesses = 2; 

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    static const int nwavefuncs = 52; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 96; 
    std::complex<double> amp[namplitudes]; 
    double matrix_uxdx_epemsuxdxsx_no_h(); 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // Pointer to the model parameters
    Parameters_sm_ckm_no_b_mass * pars; 

    // vector with external particle masses
    vector<double> mME; 

    // vector with momenta (to be changed each event)
    vector < double * > p; 
    // Initial particle ids
    int id1, id2; 

}; 


#endif  // MG5_Sigma_sm_ckm_no_b_mass_uxdx_epemsuxdxsx_H
