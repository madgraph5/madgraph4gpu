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

#include "Parameters_loop_sm_no_b_mass.h"

#include <iomanip>
#include <iostream>

#ifndef MGONGPU_HARDCODE_PARAM

// Initialize static instance
Parameters_loop_sm_no_b_mass* Parameters_loop_sm_no_b_mass::instance = 0;

// Function to get static instance - only one instance per program
Parameters_loop_sm_no_b_mass*
Parameters_loop_sm_no_b_mass::getInstance()
{
  if( instance == 0 )
    instance = new Parameters_loop_sm_no_b_mass();
  return instance;
}

void
Parameters_loop_sm_no_b_mass::setIndependentParameters( SLHAReader& slha )
{
  zero = 0; // define "zero"
  ZERO = 0; // define "zero"
  //std::vector<int> indices(2, 0); // prepare a vector for indices
  mdl_WH = slha.get_block_entry( "decay", 25, 6.382339e-03 );
  mdl_WW = slha.get_block_entry( "decay", 24, 2.047600e+00 );
  mdl_WZ = slha.get_block_entry( "decay", 23, 2.441404e+00 );
  mdl_WT = slha.get_block_entry( "decay", 6, 1.491500e+00 );
  mdl_ymtau = slha.get_block_entry( "yukawa", 15, 1.777000e+00 );
  mdl_ymt = slha.get_block_entry( "yukawa", 6, 1.730000e+02 );
  //aS = slha.get_block_entry( "sminputs", 3, 1.180000e-01 ); // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  mdl_Gf = slha.get_block_entry( "sminputs", 2, 1.166390e-05 );
  aEWM1 = slha.get_block_entry( "sminputs", 1, 1.325070e+02 );
  mdl_MH = slha.get_block_entry( "mass", 25, 1.250000e+02 );
  mdl_MZ = slha.get_block_entry( "mass", 23, 9.118800e+01 );
  mdl_MTA = slha.get_block_entry( "mass", 15, 1.777000e+00 );
  mdl_MT = slha.get_block_entry( "mass", 6, 1.730000e+02 );
  MU_R = slha.get_block_entry( "loop", 1, 9.118800e+01 );
  mdl_lhv = 1.;
  mdl_conjg__CKM3x3 = 1.;
  mdl_conjg__CKM22 = 1.;
  mdl_I4x33 = 0.;
  mdl_I1x33 = 0.;
  mdl_CKM3x3 = 1.;
  mdl_CKM33 = 1.;
  mdl_CKM22 = 1.;
  mdl_Ncol = 3.;
  mdl_CA = 3.;
  mdl_TF = 0.5;
  mdl_CF = ( 4. / 3. );
  mdl_complexi = cxsmpl<double>( 0., 1. );
  mdl_MZ__exp__2 = ( ( mdl_MZ ) * ( mdl_MZ ) );
  mdl_MZ__exp__4 = ( ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) );
  mdl_sqrt__2 = sqrt( 2. );
  mdl_MH__exp__2 = ( ( mdl_MH ) * ( mdl_MH ) );
  mdl_Ncol__exp__2 = ( ( mdl_Ncol ) * ( mdl_Ncol ) );
  mdl_MT__exp__2 = ( ( mdl_MT ) * ( mdl_MT ) );
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
  mdl_v__exp__2 = ( ( mdl_v ) * ( mdl_v ) );
  mdl_lam = mdl_MH__exp__2 / ( 2. * mdl_v__exp__2 );
  mdl_yt = ( mdl_ymt * mdl_sqrt__2 ) / mdl_v;
  mdl_ytau = ( mdl_ymtau * mdl_sqrt__2 ) / mdl_v;
  mdl_muH = sqrt( mdl_lam * mdl_v__exp__2 );
  mdl_AxialZUp = ( 3. / 2. ) * ( -( mdl_ee * mdl_sw ) / ( 6. * mdl_cw ) ) - ( 1. / 2. ) * ( ( mdl_cw * mdl_ee ) / ( 2. * mdl_sw ) );
  mdl_AxialZDown = ( -1. / 2. ) * ( -( mdl_cw * mdl_ee ) / ( 2. * mdl_sw ) ) + ( -3. / 2. ) * ( -( mdl_ee * mdl_sw ) / ( 6. * mdl_cw ) );
  mdl_VectorZUp = ( 1. / 2. ) * ( ( mdl_cw * mdl_ee ) / ( 2. * mdl_sw ) ) + ( 5. / 2. ) * ( -( mdl_ee * mdl_sw ) / ( 6. * mdl_cw ) );
  mdl_VectorZDown = ( 1. / 2. ) * ( -( mdl_cw * mdl_ee ) / ( 2. * mdl_sw ) ) + ( -1. / 2. ) * ( -( mdl_ee * mdl_sw ) / ( 6. * mdl_cw ) );
  mdl_VectorAUp = ( 2. * mdl_ee ) / 3.;
  mdl_VectorADown = -( mdl_ee ) / 3.;
  mdl_VectorWmDxU = ( 1. / 2. ) * ( ( mdl_ee ) / ( mdl_sw * mdl_sqrt__2 ) );
  mdl_AxialWmDxU = ( -1. / 2. ) * ( ( mdl_ee ) / ( mdl_sw * mdl_sqrt__2 ) );
  mdl_VectorWpUxD = ( 1. / 2. ) * ( ( mdl_ee ) / ( mdl_sw * mdl_sqrt__2 ) );
  mdl_AxialWpUxD = -( 1. / 2. ) * ( ( mdl_ee ) / ( mdl_sw * mdl_sqrt__2 ) );
  mdl_I2x33 = mdl_yt * mdl_conjg__CKM3x3;
  mdl_I3x33 = mdl_CKM3x3 * mdl_yt;
  mdl_Vector_tbGp = mdl_I1x33 - mdl_I2x33;
  mdl_Axial_tbGp = -mdl_I2x33 - mdl_I1x33;
  mdl_Vector_tbGm = mdl_I3x33 - mdl_I4x33;
  mdl_Axial_tbGm = -mdl_I4x33 - mdl_I3x33;
  mdl_gw__exp__2 = ( ( mdl_gw ) * ( mdl_gw ) );
  mdl_cw__exp__2 = ( ( mdl_cw ) * ( mdl_cw ) );
  mdl_ee__exp__2 = ( ( mdl_ee ) * ( mdl_ee ) );
  mdl_sw__exp__2 = ( ( mdl_sw ) * ( mdl_sw ) );
  mdl_yt__exp__2 = ( ( mdl_yt ) * ( mdl_yt ) );
  mdl_MU_R__exp__2 = ( ( MU_R ) * ( MU_R ) );
}

void
Parameters_loop_sm_no_b_mass::setIndependentCouplings()
{
  GC_11 = ( mdl_ee * mdl_complexi ) / ( mdl_sw * mdl_sqrt__2 );
}

/*
void
Parameters_loop_sm_no_b_mass::setDependentParameters() // now computed event-by-event (running alphas #373)
{
  mdl_sqrt__aS = sqrt( aS );
  G = 2. * mdl_sqrt__aS * sqrt( M_PI );
  mdl_G__exp__4 = ( ( G ) * ( G ) * ( G ) * ( G ) );
  mdl_RGR2_FIN_ = -( 3. / 2. ) * mdl_G__exp__4 / ( 96. * ( ( M_PI ) * ( M_PI ) ) );
  mdl_G__exp__2 = ( ( G ) * ( G ) );
  mdl_R2MixedFactor_FIN_ = -( mdl_G__exp__2 * ( 1. + mdl_lhv ) * ( mdl_Ncol__exp__2 - 1. ) ) / ( 2. * mdl_Ncol * 16. * ( ( M_PI ) * ( M_PI ) ) );
  mdl_G_UVg_1EPS_ = -( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 11. * mdl_CA;
  mdl_G_UVq_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF;
  mdl_G_UVc_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF;
  mdl_G_UVb_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF;
  mdl_G_UVt_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF;
  mdl_GWcft_UV_t_1EPS_ = COND( mdl_MT, 0., - ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF );
  mdl_tWcft_UV_1EPS_ = COND( mdl_MT, 0., - ( ( mdl_G__exp__2 ) / ( 2. * 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF );
  mdl_tMass_UV_1EPS_ = COND( mdl_MT, 0., mdl_complexi * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * mdl_MT );
  mdl_UV_yuk_c_1EPS_ = -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * 2.;
  mdl_UV_yuk_b_1EPS_ = -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * 2.;
  mdl_UV_yuk_t_1EPS_ = -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * 2.;
  mdl_R2_GGGpGm_factor_FIN_ = mdl_G__exp__2 / ( 8. * ( ( M_PI ) * ( M_PI ) ) );
  mdl_R2_GGG0G0_factor_FIN_ = mdl_G__exp__2 / ( 8. * ( ( M_PI ) * ( M_PI ) ) );
  mdl_G__exp__3 = ( ( G ) * ( G ) * ( G ) );
  mdl_G_UVt_FIN_ = COND( mdl_MT, 0., - ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) );
  mdl_GWcft_UV_t_FIN_ = COND( mdl_MT, 0., ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) );
  mdl_tWcft_UV_FIN_ = COND( mdl_MT, 0., - ( ( mdl_G__exp__2 ) / ( 2. * 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * mdl_CF * ( 4. - 3. * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) ) );
  mdl_tMass_UV_FIN_ = COND( mdl_MT, 0., mdl_complexi * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * mdl_CF * ( 4. - 3. * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) ) * mdl_MT );
  mdl_UV_yuk_t_FIN_ = COND( mdl_MT, 0., - ( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * mdl_CF * ( -3. * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) + 4. ) * 2. );
}

void
Parameters_loop_sm_no_b_mass::setDependentCouplings() // now computed event-by-event (running alphas #373)
{
  GC_4 = -G;
  GC_5 = mdl_complexi * G;
}
*/

#endif

// Routines for printing out parameters
void
Parameters_loop_sm_no_b_mass::printIndependentParameters()
{
  std::cout << "loop_sm_no_b_mass model parameters independent of event kinematics:" << std::endl;
  std::cout << "(Warning: aS in the runcard is ignored because event-by-event Gs are hardcoded or retrieved from Fortran)" << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WZ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WZ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WT = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WT << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ymtau = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ymtau << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ymt = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ymt << std::endl;
  //std::cout << std::setw( 20 ) << "aS = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << aS << std::endl; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  std::cout << std::setw( 20 ) << "mdl_Gf = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_Gf << std::endl;
  std::cout << std::setw( 20 ) << "aEWM1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << aEWM1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MTA = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MTA << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MT = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MT << std::endl;
  std::cout << std::setw( 20 ) << "MU_R = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << MU_R << std::endl;
  std::cout << std::setw( 20 ) << "mdl_lhv = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_lhv << std::endl;
  std::cout << std::setw( 20 ) << "mdl_conjg__CKM3x3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_conjg__CKM3x3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_conjg__CKM22 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_conjg__CKM22 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_I4x33 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_I4x33 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_I1x33 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_I1x33 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_CKM3x3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CKM3x3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_CKM33 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CKM33 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_CKM22 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CKM22 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_Ncol = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_Ncol << std::endl;
  std::cout << std::setw( 20 ) << "mdl_CA = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CA << std::endl;
  std::cout << std::setw( 20 ) << "mdl_TF = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_TF << std::endl;
  std::cout << std::setw( 20 ) << "mdl_CF = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_CF << std::endl;
  std::cout << std::setw( 20 ) << "mdl_complexi = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_complexi << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_Ncol__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_Ncol__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MT__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MT__exp__2 << std::endl;
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
  std::cout << std::setw( 20 ) << "mdl_v__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_v__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_lam = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_lam << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yt = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yt << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ytau = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ytau << std::endl;
  std::cout << std::setw( 20 ) << "mdl_muH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_muH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_AxialZUp = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_AxialZUp << std::endl;
  std::cout << std::setw( 20 ) << "mdl_AxialZDown = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_AxialZDown << std::endl;
  std::cout << std::setw( 20 ) << "mdl_VectorZUp = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_VectorZUp << std::endl;
  std::cout << std::setw( 20 ) << "mdl_VectorZDown = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_VectorZDown << std::endl;
  std::cout << std::setw( 20 ) << "mdl_VectorAUp = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_VectorAUp << std::endl;
  std::cout << std::setw( 20 ) << "mdl_VectorADown = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_VectorADown << std::endl;
  std::cout << std::setw( 20 ) << "mdl_VectorWmDxU = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_VectorWmDxU << std::endl;
  std::cout << std::setw( 20 ) << "mdl_AxialWmDxU = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_AxialWmDxU << std::endl;
  std::cout << std::setw( 20 ) << "mdl_VectorWpUxD = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_VectorWpUxD << std::endl;
  std::cout << std::setw( 20 ) << "mdl_AxialWpUxD = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_AxialWpUxD << std::endl;
  std::cout << std::setw( 20 ) << "mdl_I2x33 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_I2x33 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_I3x33 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_I3x33 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_Vector_tbGp = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_Vector_tbGp << std::endl;
  std::cout << std::setw( 20 ) << "mdl_Axial_tbGp = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_Axial_tbGp << std::endl;
  std::cout << std::setw( 20 ) << "mdl_Vector_tbGm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_Vector_tbGm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_Axial_tbGm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_Axial_tbGm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gw__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gw__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cw__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cw__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ee__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ee__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sw__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sw__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yt__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yt__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MU_R__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MU_R__exp__2 << std::endl;
}

void
Parameters_loop_sm_no_b_mass::printIndependentCouplings()
{
  std::cout << "loop_sm_no_b_mass model couplings independent of event kinematics:" << std::endl;
  std::cout << std::setw( 20 ) << "GC_11 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_11 << std::endl;
}

/*
void
Parameters_loop_sm_no_b_mass::printDependentParameters() // now computed event-by-event (running alphas #373)
{
  std::cout << "loop_sm_no_b_mass model parameters dependent on event kinematics:" << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__aS = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__aS << std::endl;
  std::cout << std::setw( 20 ) << "G = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << G << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_RGR2_FIN_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_RGR2_FIN_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_R2MixedFactor_FIN_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_R2MixedFactor_FIN_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G_UVg_1EPS_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G_UVg_1EPS_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G_UVq_1EPS_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G_UVq_1EPS_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G_UVc_1EPS_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G_UVc_1EPS_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G_UVb_1EPS_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G_UVb_1EPS_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G_UVt_1EPS_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G_UVt_1EPS_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_GWcft_UV_t_1EPS_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_GWcft_UV_t_1EPS_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_tWcft_UV_1EPS_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_tWcft_UV_1EPS_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_tMass_UV_1EPS_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_tMass_UV_1EPS_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_UV_yuk_c_1EPS_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_UV_yuk_c_1EPS_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_UV_yuk_b_1EPS_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_UV_yuk_b_1EPS_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_UV_yuk_t_1EPS_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_UV_yuk_t_1EPS_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_R2_GGGpGm_factor_FIN_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_R2_GGGpGm_factor_FIN_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_R2_GGG0G0_factor_FIN_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_R2_GGG0G0_factor_FIN_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G__exp__3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G__exp__3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G_UVt_FIN_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G_UVt_FIN_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_GWcft_UV_t_FIN_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_GWcft_UV_t_FIN_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_tWcft_UV_FIN_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_tWcft_UV_FIN_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_tMass_UV_FIN_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_tMass_UV_FIN_ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_UV_yuk_t_FIN_ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_UV_yuk_t_FIN_ << std::endl;
}

void
Parameters_loop_sm_no_b_mass::printDependentCouplings() // now computed event-by-event (running alphas #373)
{
  std::cout << "loop_sm_no_b_mass model couplings dependent on event kinematics:" << std::endl;
  std::cout << std::setw( 20 ) << "GC_4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_4 << std::endl;
  std::cout << std::setw( 20 ) << "GC_5 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_5 << std::endl;
}
*/
