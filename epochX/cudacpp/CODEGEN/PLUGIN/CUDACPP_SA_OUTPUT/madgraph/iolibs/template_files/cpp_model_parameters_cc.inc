// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: J. Teig, A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
%(info_lines)s
//==========================================================================

#include "Parameters_%(model_name)s.h"

#include <iomanip>
#include <iostream>

#ifndef MGONGPU_HARDCODE_PARAM

// Initialize static instance
Parameters_%(model_name)s* Parameters_%(model_name)s::instance = 0;

// Function to get static instance - only one instance per program
Parameters_%(model_name)s*
Parameters_%(model_name)s::getInstance()
{
  if( instance == 0 )
    instance = new Parameters_%(model_name)s();
  return instance;
}

void
Parameters_%(model_name)s::setIndependentParameters( SLHAReader& slha )
{
  zero = 0; // define "zero"
  ZERO = 0; // define "zero"
  //std::vector<int> indices(2, 0); // prepare a vector for indices
  %(set_independent_parameters)s
}

void
Parameters_%(model_name)s::setIndependentCouplings()
{
  %(set_independent_couplings)s
}

/*
void
Parameters_%(model_name)s::setDependentParameters() // now computed event-by-event (running alphas #373)
{
  %(set_dependent_parameters)s
}

void
Parameters_%(model_name)s::setDependentCouplings() // now computed event-by-event (running alphas #373)
{
  %(set_dependent_couplings)s
}
*/

#endif

// Routines for printing out parameters
void
Parameters_%(model_name)s::printIndependentParameters()
{
  std::cout << "%(model_name)s model parameters independent of event kinematics:" << std::endl;
  std::cout << "(Warning: aS in the runcard is ignored because event-by-event Gs are hardcoded or retrieved from Fortran)" << std::endl;
  %(print_independent_parameters)s
}

void
Parameters_%(model_name)s::printIndependentCouplings()
{
  std::cout << "%(model_name)s model couplings independent of event kinematics:" << std::endl;
  %(print_independent_couplings)s
}

/*
void
Parameters_%(model_name)s::printDependentParameters() // now computed event-by-event (running alphas #373)
{
  std::cout << "%(model_name)s model parameters dependent on event kinematics:" << std::endl;
  %(print_dependent_parameters)s
}

void
Parameters_%(model_name)s::printDependentCouplings() // now computed event-by-event (running alphas #373)
{
  std::cout << "%(model_name)s model couplings dependent on event kinematics:" << std::endl;
  %(print_dependent_couplings)s
}
*/
