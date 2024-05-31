// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: J. Teig, A. Valassi (2021-2024) for the MG5aMC CUDACPP plugin.
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.5.3_lo_vect, 2023-12-23
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "Parameters_EWdim6.h"

#include <iomanip>
#include <iostream>

#ifdef MGONGPUCPP_GPUIMPL
using namespace mg5amcGpu;
#else
using namespace mg5amcCpu;
#endif

#ifndef MGONGPU_HARDCODE_PARAM

// Initialize static instance
Parameters_EWdim6* Parameters_EWdim6::instance = 0;

// Function to get static instance - only one instance per program
Parameters_EWdim6*
Parameters_EWdim6::getInstance()
{
  if( instance == 0 )
    instance = new Parameters_EWdim6();
  return instance;
}

void
Parameters_EWdim6::setIndependentParameters( SLHAReader& slha )
{
  zero = 0;                         // define "zero"
  ZERO = 0;                         // define "zero"
  std::vector<int> indices( 2, 0 ); // prepare a vector for indices
  mdl_WH1 = slha.get_block_entry( "decay", 9000006, 6.382339e-03 );
  mdl_WH = slha.get_block_entry( "decay", 25, 6.382339e-03 );
  mdl_WW = slha.get_block_entry( "decay", 24, 2.085000e+00 );
  mdl_WZ = slha.get_block_entry( "decay", 23, 2.495200e+00 );
  mdl_WTau = slha.get_block_entry( "decay", 15, 2.270000e-12 );
  mdl_WT = slha.get_block_entry( "decay", 6, 1.508336e+00 );
  mdl_ymtau = slha.get_block_entry( "yukawa", 15, 1.777000e+00 );
  mdl_ymt = slha.get_block_entry( "yukawa", 6, 1.720000e+02 );
  mdl_ymb = slha.get_block_entry( "yukawa", 5, 4.700000e+00 );
  //aS = slha.get_block_entry( "sminputs", 3, 1.184000e-01 ); // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  mdl_Gf = slha.get_block_entry( "sminputs", 2, 1.166370e-05 );
  aEWM1 = slha.get_block_entry( "sminputs", 1, 1.279000e+02 );
  mdl_MP = slha.get_block_entry( "mass", 9000006, 1.250120e+02 );
  mdl_MH = slha.get_block_entry( "mass", 25, 1.250000e+02 );
  mdl_MZ = slha.get_block_entry( "mass", 23, 9.118760e+01 );
  mdl_MTA = slha.get_block_entry( "mass", 15, 1.777000e+00 );
  mdl_MT = slha.get_block_entry( "mass", 6, 1.720000e+02 );
  mdl_MB = slha.get_block_entry( "mass", 5, 4.700000e+00 );
  mdl_CPWL2 = slha.get_block_entry( "dim6", 5, 9.999999e+03 );
  mdl_CPWWWL2 = slha.get_block_entry( "dim6", 4, 9.999999e+02 );
  mdl_CBL2 = slha.get_block_entry( "dim6", 3, 9.999999e+01 );
  mdl_CWL2 = slha.get_block_entry( "dim6", 2, 9.999999e+00 );
  mdl_CWWWL2 = slha.get_block_entry( "dim6", 1, 1.000000e+00 );
  mdl_conjg__CKM1x1 = 1.;
  mdl_CphidL2__exp__2 = 0.;
  mdl_CphidL2 = 0.;
  mdl_CphiWL2 = 0.;
  mdl_complexi = cxsmpl<double>( 0., 1. );
  mdl_nb__2__exp__0_25 = pow( 2., 0.25 );
  mdl_vev = sqrt( 1. / mdl_Gf ) / mdl_nb__2__exp__0_25;
  mdl_MZ__exp__2 = ( ( mdl_MZ ) * ( mdl_MZ ) );
  mdl_MZ__exp__4 = ( ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) );
  mdl_sqrt__2 = sqrt( 2. );
  mdl_MH__exp__6 = pow( mdl_MH, 6. );
  mdl_MT__exp__6 = pow( mdl_MT, 6. );
  mdl_MH__exp__4 = ( ( mdl_MH ) * ( mdl_MH ) * ( mdl_MH ) * ( mdl_MH ) );
  mdl_MT__exp__4 = ( ( mdl_MT ) * ( mdl_MT ) * ( mdl_MT ) * ( mdl_MT ) );
  mdl_MH__exp__2 = ( ( mdl_MH ) * ( mdl_MH ) );
  mdl_MT__exp__2 = ( ( mdl_MT ) * ( mdl_MT ) );
  mdl_vev__exp__2 = ( ( mdl_vev ) * ( mdl_vev ) );
  mdl_vev__exp__4 = ( ( mdl_vev ) * ( mdl_vev ) * ( mdl_vev ) * ( mdl_vev ) );
  mdl_lam = ( 4. * mdl_MH__exp__2 - ( mdl_CphidL2 * mdl_MH__exp__2 * mdl_vev__exp__2 ) / 250000. ) / ( 8. * mdl_vev__exp__2 * ( 1. - ( 3. * mdl_CphidL2 * mdl_vev__exp__2 ) / 1.e6 + ( mdl_CphidL2__exp__2 * mdl_vev__exp__4 ) / 5.e11 ) );
  mdl_yb = ( mdl_ymb * mdl_sqrt__2 ) / mdl_vev;
  mdl_yt = ( mdl_ymt * mdl_sqrt__2 ) / mdl_vev;
  mdl_ytau = ( mdl_ymtau * mdl_sqrt__2 ) / mdl_vev;
  mdl_muH = sqrt( ( -4. * mdl_lam * mdl_vev__exp__2 + ( mdl_CphidL2 * mdl_lam * mdl_vev__exp__4 ) / 250000. ) / ( -1. + ( mdl_CphidL2 * mdl_vev__exp__2 ) / 1.e6 ) ) / 2.;
  mdl_MH__exp__12 = pow( mdl_MH, 12. );
  mdl_MH__exp__10 = pow( mdl_MH, 10. );
  mdl_MH__exp__8 = pow( mdl_MH, 8. );
  mdl_vev__exp__3 = ( ( mdl_vev ) * ( mdl_vev ) * ( mdl_vev ) );
  mdl_aEW = 1. / aEWM1;
  mdl_MW = sqrt( mdl_MZ__exp__2 / 2. + sqrt( mdl_MZ__exp__4 / 4. - ( mdl_aEW * M_PI * mdl_MZ__exp__2 ) / ( mdl_Gf * mdl_sqrt__2 ) ) );
  mdl_sqrt__aEW = sqrt( mdl_aEW );
  mdl_ee = 2. * mdl_sqrt__aEW * sqrt( M_PI );
  mdl_MW__exp__2 = ( ( mdl_MW ) * ( mdl_MW ) );
  mdl_sw2 = 1. - mdl_MW__exp__2 / mdl_MZ__exp__2;
  mdl_ee__exp__2 = ( ( mdl_ee ) * ( mdl_ee ) );
  mdl_MW__exp__12 = pow( mdl_MW, 12. );
  mdl_MW__exp__10 = pow( mdl_MW, 10. );
  mdl_MW__exp__8 = pow( mdl_MW, 8. );
  mdl_MW__exp__6 = pow( mdl_MW, 6. );
  mdl_MW__exp__4 = ( ( mdl_MW ) * ( mdl_MW ) * ( mdl_MW ) * ( mdl_MW ) );
  mdl_AH = ( 47. * mdl_ee__exp__2 * ( 1. - ( 2. * mdl_MH__exp__4 ) / ( 987. * mdl_MT__exp__4 ) - ( 14. * mdl_MH__exp__2 ) / ( 705. * mdl_MT__exp__2 ) + ( 213. * mdl_MH__exp__12 ) / ( 2.634632e7 * mdl_MW__exp__12 ) + ( 5. * mdl_MH__exp__10 ) / ( 119756. * mdl_MW__exp__10 ) + ( 41. * mdl_MH__exp__8 ) / ( 180950. * mdl_MW__exp__8 ) + ( 87. * mdl_MH__exp__6 ) / ( 65800. * mdl_MW__exp__6 ) + ( 57. * mdl_MH__exp__4 ) / ( 6580. * mdl_MW__exp__4 ) + ( 33. * mdl_MH__exp__2 ) / ( 470. * mdl_MW__exp__2 ) ) ) / ( 72. * ( ( M_PI ) * ( M_PI ) ) * mdl_vev );
  mdl_cw = sqrt( 1. - mdl_sw2 );
  mdl_sqrt__sw2 = sqrt( mdl_sw2 );
  mdl_sw = mdl_sqrt__sw2;
  mdl_g1 = mdl_ee / mdl_cw;
  mdl_gw = mdl_ee / mdl_sw;
  mdl_ee__exp__4 = ( ( mdl_ee ) * ( mdl_ee ) * ( mdl_ee ) * ( mdl_ee ) );
  mdl_sw__exp__4 = ( ( mdl_sw ) * ( mdl_sw ) * ( mdl_sw ) * ( mdl_sw ) );
  mdl_ee__exp__3 = ( ( mdl_ee ) * ( mdl_ee ) * ( mdl_ee ) );
  mdl_sw__exp__3 = ( ( mdl_sw ) * ( mdl_sw ) * ( mdl_sw ) );
  mdl_cw__exp__2 = ( ( mdl_cw ) * ( mdl_cw ) );
  mdl_sw__exp__2 = ( ( mdl_sw ) * ( mdl_sw ) );
  mdl_ee__exp__6 = pow( mdl_ee, 6. );
  mdl_sw__exp__6 = pow( mdl_sw, 6. );
  mdl_ee__exp__5 = pow( mdl_ee, 5. );
  mdl_sw__exp__5 = pow( mdl_sw, 5. );
  mdl_cw__exp__3 = ( ( mdl_cw ) * ( mdl_cw ) * ( mdl_cw ) );
  // BSM parameters that do not depend on alphaS but are needed in the computation of alphaS-dependent couplings;
  mdl_bsmIndepParam[0] = mdl_MH;
  mdl_bsmIndepParam[1] = mdl_MT;
  mdl_bsmIndepParam[2] = mdl_vev;
  mdl_bsmIndepParam[3] = mdl_MH__exp__6;
  mdl_bsmIndepParam[4] = mdl_MT__exp__6;
  mdl_bsmIndepParam[5] = mdl_MH__exp__4;
  mdl_bsmIndepParam[6] = mdl_MT__exp__4;
  mdl_bsmIndepParam[7] = mdl_MH__exp__2;
  mdl_bsmIndepParam[8] = mdl_MT__exp__2;
}

void
Parameters_EWdim6::setIndependentCouplings()
{
  GC_38 = -( mdl_CPWWWL2 * mdl_cw * mdl_ee__exp__3 * mdl_complexi ) / ( 2.e6 * mdl_sw__exp__3 );
  GC_39 = ( 3. * mdl_cw * mdl_CWWWL2 * mdl_ee__exp__3 * mdl_complexi ) / ( 2.e6 * mdl_sw__exp__3 );
  GC_109 = ( mdl_CPWL2 * mdl_cw * mdl_ee__exp__3 * mdl_complexi * mdl_vev__exp__2 ) / ( 4.e6 * mdl_sw__exp__3 );
  GC_110 = -( mdl_cw * mdl_CWL2 * mdl_ee__exp__3 * mdl_complexi * mdl_vev__exp__2 ) / ( 8.e6 * mdl_sw__exp__3 );
  GC_113 = ( mdl_CBL2 * mdl_ee__exp__3 * mdl_complexi * mdl_vev__exp__2 ) / ( 8.e6 * mdl_cw * mdl_sw );
  GC_114 = ( mdl_CPWL2 * mdl_ee__exp__3 * mdl_complexi * mdl_vev__exp__2 ) / ( 4.e6 * mdl_cw * mdl_sw );
  GC_115 = -( mdl_CWL2 * mdl_ee__exp__3 * mdl_complexi * mdl_vev__exp__2 ) / ( 8.e6 * mdl_cw * mdl_sw );
  GC_145 = ( mdl_ee * mdl_complexi * mdl_conjg__CKM1x1 ) / ( mdl_sw * mdl_sqrt__2 );
}

/*
void
Parameters_EWdim6::setDependentParameters() // now computed event-by-event (running alphas #373)
{
  mdl_sqrt__aS = sqrt( aS );
  G = 2. * mdl_sqrt__aS * sqrt( M_PI );
  mdl_G__exp__2 = ( ( G ) * ( G ) );
  mdl_GH = -( mdl_G__exp__2 * ( 1. + ( 13. * mdl_MH__exp__6 ) / ( 16800. * mdl_MT__exp__6 ) + mdl_MH__exp__4 / ( 168. * mdl_MT__exp__4 ) + ( 7. * mdl_MH__exp__2 ) / ( 120. * mdl_MT__exp__2 ) ) ) / ( 12. * ( ( M_PI ) * ( M_PI ) ) * mdl_vev );
  mdl_Gphi = -( mdl_G__exp__2 * ( 1. + mdl_MH__exp__6 / ( 560. * mdl_MT__exp__6 ) + mdl_MH__exp__4 / ( 90. * mdl_MT__exp__4 ) + mdl_MH__exp__2 / ( 12. * mdl_MT__exp__2 ) ) ) / ( 8. * ( ( M_PI ) * ( M_PI ) ) * mdl_vev );
}

void
Parameters_EWdim6::setDependentCouplings() // now computed event-by-event (running alphas #373)
{
  // (none)
}
*/

#endif

// Routines for printing out parameters
void
Parameters_EWdim6::printIndependentParameters()
{
  std::cout << "EWdim6 model parameters independent of event kinematics:" << std::endl;
  std::cout << "(Warning: aS in the runcard is ignored because event-by-event Gs are hardcoded or retrieved from Fortran)" << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WH1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WH1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WZ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WZ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WTau = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WTau << std::endl;
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
  std::cout << std::setw( 20 ) << "mdl_CPWL2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CPWL2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_CPWWWL2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CPWWWL2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_CBL2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CBL2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_CWL2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CWL2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_CWWWL2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CWWWL2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_conjg__CKM1x1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_conjg__CKM1x1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_CphidL2__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CphidL2__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_CphidL2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CphidL2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_CphiWL2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CphiWL2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_complexi = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_complexi << std::endl;
  std::cout << std::setw( 20 ) << "mdl_nb__2__exp__0_25 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_nb__2__exp__0_25 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_vev = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_vev << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__6 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MT__exp__6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MT__exp__6 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MT__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MT__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MT__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MT__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_vev__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_vev__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_vev__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_vev__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_lam = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_lam << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yb = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yb << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yt = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yt << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ytau = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ytau << std::endl;
  std::cout << std::setw( 20 ) << "mdl_muH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_muH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__12 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__12 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__10 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__10 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_vev__exp__3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_vev__exp__3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_aEW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_aEW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__aEW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__aEW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ee = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ee << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sw2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sw2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ee__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ee__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__12 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__12 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__10 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__10 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__6 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_AH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_AH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cw = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cw << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__sw2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__sw2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sw = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sw << std::endl;
  std::cout << std::setw( 20 ) << "mdl_g1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_g1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gw = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gw << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ee__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ee__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sw__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sw__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ee__exp__3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ee__exp__3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sw__exp__3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sw__exp__3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cw__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cw__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sw__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sw__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ee__exp__6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ee__exp__6 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sw__exp__6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sw__exp__6 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ee__exp__5 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ee__exp__5 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sw__exp__5 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sw__exp__5 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cw__exp__3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cw__exp__3 << std::endl;
}

void
Parameters_EWdim6::printIndependentCouplings()
{
  std::cout << "EWdim6 model couplings independent of event kinematics:" << std::endl;
  std::cout << std::setw( 20 ) << "GC_38 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_38 << std::endl;
  std::cout << std::setw( 20 ) << "GC_39 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_39 << std::endl;
  std::cout << std::setw( 20 ) << "GC_109 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_109 << std::endl;
  std::cout << std::setw( 20 ) << "GC_110 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_110 << std::endl;
  std::cout << std::setw( 20 ) << "GC_113 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_113 << std::endl;
  std::cout << std::setw( 20 ) << "GC_114 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_114 << std::endl;
  std::cout << std::setw( 20 ) << "GC_115 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_115 << std::endl;
  std::cout << std::setw( 20 ) << "GC_145 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_145 << std::endl;
}

/*
void
Parameters_EWdim6::printDependentParameters() // now computed event-by-event (running alphas #373)
{
  std::cout << "EWdim6 model parameters dependent on event kinematics:" << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__aS = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__aS << std::endl;
  std::cout << std::setw( 20 ) << "G = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << G << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_GH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_GH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_Gphi = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_Gphi << std::endl;
}

void
Parameters_EWdim6::printDependentCouplings() // now computed event-by-event (running alphas #373)
{
  std::cout << "EWdim6 model couplings dependent on event kinematics:" << std::endl;
  // (none)
}
*/
