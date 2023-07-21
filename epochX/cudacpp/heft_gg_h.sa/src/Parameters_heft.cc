// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.5.0_lo_vect, 2023-06-09
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "Parameters_heft.h"

#include <iomanip>
#include <iostream>

#ifdef __CUDACC__
using namespace mg5amcGpu;
#else
using namespace mg5amcCpu;
#endif

#ifndef MGONGPU_HARDCODE_PARAM

// Initialize static instance
Parameters_heft* Parameters_heft::instance = 0;

// Function to get static instance - only one instance per program
Parameters_heft*
Parameters_heft::getInstance()
{
  if( instance == 0 )
    instance = new Parameters_heft();
  return instance;
}

void
Parameters_heft::setIndependentParameters( SLHAReader& slha )
{
  zero = 0; // define "zero"
  ZERO = 0; // define "zero"
  //std::vector<int> indices(2, 0); // prepare a vector for indices
  mdl_WH1 = slha.get_block_entry( "decay", 9000006, 6.382339e-03 );
  mdl_WH = slha.get_block_entry( "decay", 25, 6.382339e-03 );
  mdl_WW = slha.get_block_entry( "decay", 24, 2.047600e+00 );
  mdl_WZ = slha.get_block_entry( "decay", 23, 2.441404e+00 );
  mdl_WT = slha.get_block_entry( "decay", 6, 1.491500e+00 );
  mdl_ymtau = slha.get_block_entry( "yukawa", 15, 1.777000e+00 );
  mdl_ymt = slha.get_block_entry( "yukawa", 6, 1.645000e+02 );
  mdl_ymb = slha.get_block_entry( "yukawa", 5, 4.200000e+00 );
  //aS = slha.get_block_entry( "sminputs", 3, 1.180000e-01 ); // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  mdl_Gf = slha.get_block_entry( "sminputs", 2, 1.166390e-05 );
  aEWM1 = slha.get_block_entry( "sminputs", 1, 1.325070e+02 );
  mdl_MP = slha.get_block_entry( "mass", 9000006, 1.250001e+02 );
  mdl_MH = slha.get_block_entry( "mass", 25, 1.250000e+02 );
  mdl_MZ = slha.get_block_entry( "mass", 23, 9.118800e+01 );
  mdl_MTA = slha.get_block_entry( "mass", 15, 1.777000e+00 );
  mdl_MT = slha.get_block_entry( "mass", 6, 1.730000e+02 );
  mdl_MB = slha.get_block_entry( "mass", 5, 4.700000e+00 );
  mdl_conjg__CKM3x3 = 1.;
  mdl_complexi = cxsmpl<double>( 0., 1. );
  mdl_MZ__exp__2 = ( ( mdl_MZ ) * ( mdl_MZ ) );
  mdl_MZ__exp__4 = ( ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) );
  mdl_sqrt__2 = sqrt( 2. );
  mdl_MH__exp__4 = ( ( mdl_MH ) * ( mdl_MH ) * ( mdl_MH ) * ( mdl_MH ) );
  mdl_MT__exp__4 = ( ( mdl_MT ) * ( mdl_MT ) * ( mdl_MT ) * ( mdl_MT ) );
  mdl_MH__exp__2 = ( ( mdl_MH ) * ( mdl_MH ) );
  mdl_MT__exp__2 = ( ( mdl_MT ) * ( mdl_MT ) );
  mdl_MH__exp__12 = pow( mdl_MH, 12. );
  mdl_MH__exp__10 = pow( mdl_MH, 10. );
  mdl_MH__exp__8 = pow( mdl_MH, 8. );
  mdl_MH__exp__6 = pow( mdl_MH, 6. );
  mdl_MT__exp__6 = pow( mdl_MT, 6. );
  mdl_aEW = 1. / aEWM1;
  mdl_MW = sqrt( mdl_MZ__exp__2 / 2. + sqrt( mdl_MZ__exp__4 / 4. - ( mdl_aEW * M_PI * mdl_MZ__exp__2 ) / ( mdl_Gf * mdl_sqrt__2 ) ) );
  mdl_sqrt__aEW = sqrt( mdl_aEW );
  mdl_ee = 2. * mdl_sqrt__aEW * sqrt( M_PI );
  mdl_MW__exp__2 = ( ( mdl_MW ) * ( mdl_MW ) );
  mdl_sw2 = 1. - mdl_MW__exp__2 / mdl_MZ__exp__2;
  mdl_cw = sqrt( 1. - mdl_sw2 );
  mdl_sqrt__sw2 = sqrt( mdl_sw2 );
  mdl_sw = mdl_sqrt__sw2;
  mdl_g1 = mdl_ee / mdl_cw;
  mdl_gw = mdl_ee / mdl_sw;
  mdl_v = ( 2. * mdl_MW * mdl_sw ) / mdl_ee;
  mdl_ee__exp__2 = ( ( mdl_ee ) * ( mdl_ee ) );
  mdl_MW__exp__12 = pow( mdl_MW, 12. );
  mdl_MW__exp__10 = pow( mdl_MW, 10. );
  mdl_MW__exp__8 = pow( mdl_MW, 8. );
  mdl_MW__exp__6 = pow( mdl_MW, 6. );
  mdl_MW__exp__4 = ( ( mdl_MW ) * ( mdl_MW ) * ( mdl_MW ) * ( mdl_MW ) );
  mdl_AH = ( 47. * mdl_ee__exp__2 * ( 1. - ( 2. * mdl_MH__exp__4 ) / ( 987. * mdl_MT__exp__4 ) - ( 14. * mdl_MH__exp__2 ) / ( 705. * mdl_MT__exp__2 ) + ( 213. * mdl_MH__exp__12 ) / ( 2.634632e7 * mdl_MW__exp__12 ) + ( 5. * mdl_MH__exp__10 ) / ( 119756. * mdl_MW__exp__10 ) + ( 41. * mdl_MH__exp__8 ) / ( 180950. * mdl_MW__exp__8 ) + ( 87. * mdl_MH__exp__6 ) / ( 65800. * mdl_MW__exp__6 ) + ( 57. * mdl_MH__exp__4 ) / ( 6580. * mdl_MW__exp__4 ) + ( 33. * mdl_MH__exp__2 ) / ( 470. * mdl_MW__exp__2 ) ) ) / ( 72. * ( ( M_PI ) * ( M_PI ) ) * mdl_v );
  mdl_v__exp__2 = ( ( mdl_v ) * ( mdl_v ) );
  mdl_lam = mdl_MH__exp__2 / ( 2. * mdl_v__exp__2 );
  mdl_yb = ( mdl_ymb * mdl_sqrt__2 ) / mdl_v;
  mdl_yt = ( mdl_ymt * mdl_sqrt__2 ) / mdl_v;
  mdl_ytau = ( mdl_ymtau * mdl_sqrt__2 ) / mdl_v;
  mdl_muH = sqrt( mdl_lam * mdl_v__exp__2 );
  mdl_gw__exp__2 = ( ( mdl_gw ) * ( mdl_gw ) );
  mdl_cw__exp__2 = ( ( mdl_cw ) * ( mdl_cw ) );
  mdl_sw__exp__2 = ( ( mdl_sw ) * ( mdl_sw ) );
}

void
Parameters_heft::setIndependentCouplings()
{
  // (none)
}

/*
void
Parameters_heft::setDependentParameters() // now computed event-by-event (running alphas #373)
{
  mdl_sqrt__aS = sqrt( aS );
  G = 2. * mdl_sqrt__aS * sqrt( M_PI );
  mdl_G__exp__2 = ( ( G ) * ( G ) );
  mdl_GH = -( mdl_G__exp__2 * ( 1. + ( 13. * mdl_MH__exp__6 ) / ( 16800. * mdl_MT__exp__6 ) + mdl_MH__exp__4 / ( 168. * mdl_MT__exp__4 ) + ( 7. * mdl_MH__exp__2 ) / ( 120. * mdl_MT__exp__2 ) ) ) / ( 12. * ( ( M_PI ) * ( M_PI ) ) * mdl_v );
  mdl_Gphi = -( mdl_G__exp__2 * ( 1. + mdl_MH__exp__6 / ( 560. * mdl_MT__exp__6 ) + mdl_MH__exp__4 / ( 90. * mdl_MT__exp__4 ) + mdl_MH__exp__2 / ( 12. * mdl_MT__exp__2 ) ) ) / ( 8. * ( ( M_PI ) * ( M_PI ) ) * mdl_v );
}

void
Parameters_heft::setDependentCouplings() // now computed event-by-event (running alphas #373)
{
  GC_13 = -( mdl_complexi * mdl_GH );
}
*/

#endif

// Routines for printing out parameters
void
Parameters_heft::printIndependentParameters()
{
  std::cout << "heft model parameters independent of event kinematics:" << std::endl;
  std::cout << "(Warning: aS in the runcard is ignored because event-by-event Gs are hardcoded or retrieved from Fortran)" << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WH1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WH1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WZ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WZ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WT = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WT << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ymtau = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ymtau << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ymt = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ymt << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ymb = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ymb << std::endl;
  //std::cout << std::setw( 20 ) << "aS = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << aS << std::endl; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  std::cout << std::setw( 20 ) << "mdl_Gf = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_Gf << std::endl;
  std::cout << std::setw( 20 ) << "aEWM1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << aEWM1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MP = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MP << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MTA = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MTA << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MT = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MT << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MB = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MB << std::endl;
  std::cout << std::setw( 20 ) << "mdl_conjg__CKM3x3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_conjg__CKM3x3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_complexi = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_complexi << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MT__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MT__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MT__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MT__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__12 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__12 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__10 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__10 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__6 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MT__exp__6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MT__exp__6 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_aEW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_aEW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__aEW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__aEW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ee = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ee << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sw2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sw2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cw = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cw << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__sw2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__sw2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sw = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sw << std::endl;
  std::cout << std::setw( 20 ) << "mdl_g1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_g1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gw = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gw << std::endl;
  std::cout << std::setw( 20 ) << "mdl_v = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_v << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ee__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ee__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__12 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__12 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__10 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__10 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__6 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_AH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_AH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_v__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_v__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_lam = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_lam << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yb = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yb << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yt = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yt << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ytau = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ytau << std::endl;
  std::cout << std::setw( 20 ) << "mdl_muH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_muH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gw__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gw__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cw__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cw__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sw__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sw__exp__2 << std::endl;
}

void
Parameters_heft::printIndependentCouplings()
{
  std::cout << "heft model couplings independent of event kinematics:" << std::endl;
  // (none)
}

/*
void
Parameters_heft::printDependentParameters() // now computed event-by-event (running alphas #373)
{
  std::cout << "heft model parameters dependent on event kinematics:" << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__aS = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__aS << std::endl;
  std::cout << std::setw( 20 ) << "G = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << G << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_GH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_GH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_Gphi = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_Gphi << std::endl;
}

void
Parameters_heft::printDependentCouplings() // now computed event-by-event (running alphas #373)
{
  std::cout << "heft model couplings dependent on event kinematics:" << std::endl;
  std::cout << std::setw( 20 ) << "GC_13 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_13 << std::endl;
}
*/
