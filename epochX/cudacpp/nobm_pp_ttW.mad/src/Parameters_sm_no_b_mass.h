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

#ifndef Parameters_loop_sm_no_b_mass_H
#define Parameters_loop_sm_no_b_mass_H

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"
#include "mgOnGpuVectors.h"

//==========================================================================

#ifndef MGONGPU_HARDCODE_PARAM // this is only supported in SM processes (e.g. not in EFT models) for the moment (#439)
#error This non-SM physics process only supports MGONGPU_HARDCODE_PARAM builds (#439): please run "make HRDCOD=1"

#include "read_slha.h"

class Parameters_loop_sm_no_b_mass
{
public:

  static Parameters_loop_sm_no_b_mass* getInstance();

  // Define "zero"
  double zero, ZERO;

  // Model parameters independent of aS
  //double aS; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  double mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau, mdl_ymt, mdl_Gf, aEWM1, mdl_MH, mdl_MZ, mdl_MTA, mdl_MT, MU_R, mdl_lhv, mdl_conjg__CKM3x3, mdl_conjg__CKM22, mdl_I4x33, mdl_I1x33, mdl_CKM3x3, mdl_CKM33, mdl_CKM22, mdl_Ncol, mdl_CA, mdl_TF, mdl_CF, mdl_MZ__exp__2, mdl_MZ__exp__4, mdl_sqrt__2, mdl_MH__exp__2, mdl_Ncol__exp__2, mdl_MT__exp__2, mdl_aEW, mdl_MW, mdl_sqrt__aEW, mdl_ee, mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw, mdl_v, mdl_v__exp__2, mdl_lam, mdl_yt, mdl_ytau, mdl_muH, mdl_AxialZUp, mdl_AxialZDown, mdl_VectorZUp, mdl_VectorZDown, mdl_VectorAUp, mdl_VectorADown, mdl_VectorWmDxU, mdl_AxialWmDxU, mdl_VectorWpUxD, mdl_AxialWpUxD, mdl_gw__exp__2, mdl_cw__exp__2, mdl_ee__exp__2, mdl_sw__exp__2, mdl_yt__exp__2, mdl_MU_R__exp__2;
  cxsmpl<double> mdl_complexi, mdl_I2x33, mdl_I3x33, mdl_Vector_tbGp, mdl_Axial_tbGp, mdl_Vector_tbGm, mdl_Axial_tbGm;

  // Model couplings independent of aS
  cxsmpl<double> GC_11;

  // Model parameters dependent on aS
  //double mdl_sqrt__aS, G, mdl_G__exp__4, mdl_RGR2_FIN_, mdl_G__exp__2, mdl_R2MixedFactor_FIN_, mdl_G_UVg_1EPS_, mdl_G_UVq_1EPS_, mdl_G_UVc_1EPS_, mdl_G_UVb_1EPS_, mdl_G_UVt_1EPS_, mdl_GWcft_UV_t_1EPS_, mdl_tWcft_UV_1EPS_, mdl_UV_yuk_c_1EPS_, mdl_UV_yuk_b_1EPS_, mdl_UV_yuk_t_1EPS_, mdl_R2_GGGpGm_factor_FIN_, mdl_R2_GGG0G0_factor_FIN_, mdl_G__exp__3, mdl_G_UVt_FIN_, mdl_GWcft_UV_t_FIN_, mdl_tWcft_UV_FIN_, mdl_UV_yuk_t_FIN_; // now computed event-by-event (running alphas #373)
  //cxsmpl<double> mdl_tMass_UV_1EPS_, mdl_tMass_UV_FIN_; // now computed event-by-event (running alphas #373)

  // Model couplings dependent on aS
  //cxsmpl<double> GC_4, GC_5; // now computed event-by-event (running alphas #373)

  // Set parameters that are unchanged during the run
  void setIndependentParameters( SLHAReader& slha );

  // Set couplings that are unchanged during the run
  void setIndependentCouplings();

  // Set parameters that are changed event by event
  //void setDependentParameters(); // now computed event-by-event (running alphas #373)

  // Set couplings that are changed event by event
  //void setDependentCouplings(); // now computed event-by-event (running alphas #373)

  // Print parameters that are unchanged during the run
  void printIndependentParameters();

  // Print couplings that are unchanged during the run
  void printIndependentCouplings();

  // Print parameters that are changed event by event
  //void printDependentParameters(); // now computed event-by-event (running alphas #373)

  // Print couplings that are changed event by event
  //void printDependentCouplings(); // now computed event-by-event (running alphas #373)

private:

  static Parameters_loop_sm_no_b_mass* instance;
};

#else

#include <cassert>
#include <limits>

// Hardcoded constexpr physics parameters
namespace Parameters_loop_sm_no_b_mass // keep the same name rather than HardcodedParameters_loop_sm_no_b_mass for simplicity
{
  // Constexpr implementation of sqrt (see https://stackoverflow.com/a/34134071)
  double constexpr sqrtNewtonRaphson( double x, double curr, double prev )
  {
    return curr == prev ? curr : sqrtNewtonRaphson( x, 0.5 * ( curr + x / curr ), curr );
  }
  double constexpr constexpr_sqrt( double x )
  {
    return x >= 0 // && x < std::numeric_limits<double>::infinity() // avoid -Wtautological-constant-compare warning in fast math
      ? sqrtNewtonRaphson( x, x, 0 )
      : std::numeric_limits<double>::quiet_NaN();
  }

  // Constexpr implementation of floor (see https://stackoverflow.com/a/66146159)
  constexpr int constexpr_floor( double d )
  {
    const int i = static_cast<int>( d );
    return d < i ? i - 1 : i;
  }

  // Constexpr implementation of pow
  constexpr double constexpr_pow( double base, double exp )
  {
    // NB(1): this implementation of constexpr_pow requires exponent >= 0
    assert( exp >= 0 ); // NB would fail at compile time with "error: call to non-‘constexpr’ function ‘void __assert_fail'"
    // NB(2): this implementation of constexpr_pow requires an integer exponent
    const int iexp = constexpr_floor( exp );
    assert( static_cast<double>( iexp ) == exp ); // NB would fail at compile time with "error: call to non-‘constexpr’ function ‘void __assert_fail'"
    // Iterative implementation of pow if exp is a non negative integer
    return iexp == 0 ? 1 : base * constexpr_pow( base, iexp - 1 );
  }

  // Model parameters independent of aS
  constexpr double zero = 0;
  constexpr double ZERO = 0;
  constexpr double mdl_WH = 6.382339e-03;
  constexpr double mdl_WW = 2.047600e+00;
  constexpr double mdl_WZ = 2.441404e+00;
  constexpr double mdl_WT = 1.491500e+00;
  constexpr double mdl_ymtau = 1.777000e+00;
  constexpr double mdl_ymt = 1.730000e+02;
  //constexpr double aS = 1.180000e-01; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  constexpr double mdl_Gf = 1.166390e-05;
  constexpr double aEWM1 = 1.325070e+02;
  constexpr double mdl_MH = 1.250000e+02;
  constexpr double mdl_MZ = 9.118800e+01;
  constexpr double mdl_MTA = 1.777000e+00;
  constexpr double mdl_MT = 1.730000e+02;
  constexpr double MU_R = 9.118800e+01;
  constexpr double mdl_lhv = 1.;
  constexpr double mdl_conjg__CKM3x3 = 1.;
  constexpr double mdl_conjg__CKM22 = 1.;
  constexpr double mdl_I4x33 = 0.;
  constexpr double mdl_I1x33 = 0.;
  constexpr double mdl_CKM3x3 = 1.;
  constexpr double mdl_CKM33 = 1.;
  constexpr double mdl_CKM22 = 1.;
  constexpr double mdl_Ncol = 3.;
  constexpr double mdl_CA = 3.;
  constexpr double mdl_TF = 0.5;
  constexpr double mdl_CF = ( 4. / 3. );
  constexpr cxsmpl<double> mdl_complexi = cxsmpl<double>( 0., 1. );
  constexpr double mdl_MZ__exp__2 = ( ( mdl_MZ ) * ( mdl_MZ ) );
  constexpr double mdl_MZ__exp__4 = ( ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) );
  constexpr double mdl_sqrt__2 = constexpr_sqrt( 2. );
  constexpr double mdl_MH__exp__2 = ( ( mdl_MH ) * ( mdl_MH ) );
  constexpr double mdl_Ncol__exp__2 = ( ( mdl_Ncol ) * ( mdl_Ncol ) );
  constexpr double mdl_MT__exp__2 = ( ( mdl_MT ) * ( mdl_MT ) );
  constexpr double mdl_aEW = 1. / aEWM1;
  constexpr double mdl_MW = constexpr_sqrt( mdl_MZ__exp__2 / 2. + constexpr_sqrt( mdl_MZ__exp__4 / 4. - ( mdl_aEW * M_PI * mdl_MZ__exp__2 ) / ( mdl_Gf * mdl_sqrt__2 ) ) );
  constexpr double mdl_sqrt__aEW = constexpr_sqrt( mdl_aEW );
  constexpr double mdl_ee = 2. * mdl_sqrt__aEW * constexpr_sqrt( M_PI );
  constexpr double mdl_MW__exp__2 = ( ( mdl_MW ) * ( mdl_MW ) );
  constexpr double mdl_sw2 = 1. - mdl_MW__exp__2 / mdl_MZ__exp__2;
  constexpr double mdl_cw = constexpr_sqrt( 1. - mdl_sw2 );
  constexpr double mdl_sqrt__sw2 = constexpr_sqrt( mdl_sw2 );
  constexpr double mdl_sw = mdl_sqrt__sw2;
  constexpr double mdl_g1 = mdl_ee / mdl_cw;
  constexpr double mdl_gw = mdl_ee / mdl_sw;
  constexpr double mdl_v = ( 2. * mdl_MW * mdl_sw ) / mdl_ee;
  constexpr double mdl_v__exp__2 = ( ( mdl_v ) * ( mdl_v ) );
  constexpr double mdl_lam = mdl_MH__exp__2 / ( 2. * mdl_v__exp__2 );
  constexpr double mdl_yt = ( mdl_ymt * mdl_sqrt__2 ) / mdl_v;
  constexpr double mdl_ytau = ( mdl_ymtau * mdl_sqrt__2 ) / mdl_v;
  constexpr double mdl_muH = constexpr_sqrt( mdl_lam * mdl_v__exp__2 );
  constexpr double mdl_AxialZUp = ( 3. / 2. ) * ( -( mdl_ee * mdl_sw ) / ( 6. * mdl_cw ) ) - ( 1. / 2. ) * ( ( mdl_cw * mdl_ee ) / ( 2. * mdl_sw ) );
  constexpr double mdl_AxialZDown = ( -1. / 2. ) * ( -( mdl_cw * mdl_ee ) / ( 2. * mdl_sw ) ) + ( -3. / 2. ) * ( -( mdl_ee * mdl_sw ) / ( 6. * mdl_cw ) );
  constexpr double mdl_VectorZUp = ( 1. / 2. ) * ( ( mdl_cw * mdl_ee ) / ( 2. * mdl_sw ) ) + ( 5. / 2. ) * ( -( mdl_ee * mdl_sw ) / ( 6. * mdl_cw ) );
  constexpr double mdl_VectorZDown = ( 1. / 2. ) * ( -( mdl_cw * mdl_ee ) / ( 2. * mdl_sw ) ) + ( -1. / 2. ) * ( -( mdl_ee * mdl_sw ) / ( 6. * mdl_cw ) );
  constexpr double mdl_VectorAUp = ( 2. * mdl_ee ) / 3.;
  constexpr double mdl_VectorADown = -( mdl_ee ) / 3.;
  constexpr double mdl_VectorWmDxU = ( 1. / 2. ) * ( ( mdl_ee ) / ( mdl_sw * mdl_sqrt__2 ) );
  constexpr double mdl_AxialWmDxU = ( -1. / 2. ) * ( ( mdl_ee ) / ( mdl_sw * mdl_sqrt__2 ) );
  constexpr double mdl_VectorWpUxD = ( 1. / 2. ) * ( ( mdl_ee ) / ( mdl_sw * mdl_sqrt__2 ) );
  constexpr double mdl_AxialWpUxD = -( 1. / 2. ) * ( ( mdl_ee ) / ( mdl_sw * mdl_sqrt__2 ) );
  constexpr cxsmpl<double> mdl_I2x33 = mdl_yt * mdl_conjg__CKM3x3;
  constexpr cxsmpl<double> mdl_I3x33 = mdl_CKM3x3 * mdl_yt;
  constexpr cxsmpl<double> mdl_Vector_tbGp = mdl_I1x33 - mdl_I2x33;
  constexpr cxsmpl<double> mdl_Axial_tbGp = -mdl_I2x33 - mdl_I1x33;
  constexpr cxsmpl<double> mdl_Vector_tbGm = mdl_I3x33 - mdl_I4x33;
  constexpr cxsmpl<double> mdl_Axial_tbGm = -mdl_I4x33 - mdl_I3x33;
  constexpr double mdl_gw__exp__2 = ( ( mdl_gw ) * ( mdl_gw ) );
  constexpr double mdl_cw__exp__2 = ( ( mdl_cw ) * ( mdl_cw ) );
  constexpr double mdl_ee__exp__2 = ( ( mdl_ee ) * ( mdl_ee ) );
  constexpr double mdl_sw__exp__2 = ( ( mdl_sw ) * ( mdl_sw ) );
  constexpr double mdl_yt__exp__2 = ( ( mdl_yt ) * ( mdl_yt ) );
  constexpr double mdl_MU_R__exp__2 = ( ( MU_R ) * ( MU_R ) );

  // Model couplings independent of aS
  constexpr cxsmpl<double> GC_11 = ( mdl_ee * mdl_complexi ) / ( mdl_sw * mdl_sqrt__2 );

  // Model parameters dependent on aS
  //constexpr double mdl_sqrt__aS = //constexpr_sqrt( aS ); // now computed event-by-event (running alphas #373)
  //constexpr double G = 2. * mdl_sqrt__aS * //constexpr_sqrt( M_PI ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_G__exp__4 = ( ( G ) * ( G ) * ( G ) * ( G ) ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_RGR2_FIN_ = -( 3. / 2. ) * mdl_G__exp__4 / ( 96. * ( ( M_PI ) * ( M_PI ) ) ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_G__exp__2 = ( ( G ) * ( G ) ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_R2MixedFactor_FIN_ = -( mdl_G__exp__2 * ( 1. + mdl_lhv ) * ( mdl_Ncol__exp__2 - 1. ) ) / ( 2. * mdl_Ncol * 16. * ( ( M_PI ) * ( M_PI ) ) ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_G_UVg_1EPS_ = -( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 11. * mdl_CA; // now computed event-by-event (running alphas #373)
  //constexpr double mdl_G_UVq_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF; // now computed event-by-event (running alphas #373)
  //constexpr double mdl_G_UVc_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF; // now computed event-by-event (running alphas #373)
  //constexpr double mdl_G_UVb_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF; // now computed event-by-event (running alphas #373)
  //constexpr double mdl_G_UVt_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF; // now computed event-by-event (running alphas #373)
  //constexpr double mdl_GWcft_UV_t_1EPS_ = COND( mdl_MT, 0., -( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_tWcft_UV_1EPS_ = COND( mdl_MT, 0., -( ( mdl_G__exp__2 ) / ( 2. * 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF ); // now computed event-by-event (running alphas #373)
  //constexpr cxsmpl<double> mdl_tMass_UV_1EPS_ = COND( mdl_MT, 0., mdl_complexi * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * mdl_MT ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_UV_yuk_c_1EPS_ = -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * 2.; // now computed event-by-event (running alphas #373)
  //constexpr double mdl_UV_yuk_b_1EPS_ = -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * 2.; // now computed event-by-event (running alphas #373)
  //constexpr double mdl_UV_yuk_t_1EPS_ = -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * 2.; // now computed event-by-event (running alphas #373)
  //constexpr double mdl_R2_GGGpGm_factor_FIN_ = mdl_G__exp__2 / ( 8. * ( ( M_PI ) * ( M_PI ) ) ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_R2_GGG0G0_factor_FIN_ = mdl_G__exp__2 / ( 8. * ( ( M_PI ) * ( M_PI ) ) ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_G__exp__3 = ( ( G ) * ( G ) * ( G ) ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_G_UVt_FIN_ = COND( mdl_MT, 0., -( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_GWcft_UV_t_FIN_ = COND( mdl_MT, 0., ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_tWcft_UV_FIN_ = COND( mdl_MT, 0., -( ( mdl_G__exp__2 ) / ( 2. * 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * mdl_CF * ( 4. - 3. * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) ) ); // now computed event-by-event (running alphas #373)
  //constexpr cxsmpl<double> mdl_tMass_UV_FIN_ = COND( mdl_MT, 0., mdl_complexi * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * mdl_CF * ( 4. - 3. * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) ) * mdl_MT ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_UV_yuk_t_FIN_ = COND( mdl_MT, 0., -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * mdl_CF * ( -3. * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) + 4. ) * 2. ); // now computed event-by-event (running alphas #373)

  // Model couplings dependent on aS
  //constexpr cxsmpl<double> GC_4 = -G; // now computed event-by-event (running alphas #373)
  //constexpr cxsmpl<double> GC_5 = mdl_complexi * G; // now computed event-by-event (running alphas #373)

  // Print parameters that are unchanged during the run
  void printIndependentParameters();

  // Print couplings that are unchanged during the run
  void printIndependentCouplings();

  // Print parameters that are changed event by event
  //void printDependentParameters(); // now computed event-by-event (running alphas #373)

  // Print couplings that are changed event by event
  //void printDependentCouplings(); // now computed event-by-event (running alphas #373)
}

#endif

//==========================================================================

namespace Parameters_loop_sm_no_b_mass_dependentCouplings
{
  constexpr size_t ndcoup = 2; // #couplings that vary event by event because they depend on the running alphas QCD
  constexpr size_t idcoup_GC_4 = 0;
  constexpr size_t idcoup_GC_5 = 1;
  struct DependentCouplings_sv
  {
    cxtype_sv GC_4;
    cxtype_sv GC_5;
  };
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"  // e.g. <<warning: unused variable ‘mdl_G__exp__2’ [-Wunused-variable]>>
#pragma GCC diagnostic ignored "-Wunused-parameter" // e.g. <<warning: unused parameter ‘G’ [-Wunused-parameter]>>
#ifdef __CUDACC__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 177 // e.g. <<warning #177-D: variable "mdl_G__exp__2" was declared but never referenced>>
#endif
  __host__ __device__ inline const DependentCouplings_sv computeDependentCouplings_fromG( const fptype_sv& G_sv )
  {
#ifdef MGONGPU_HARDCODE_PARAM
    using namespace Parameters_loop_sm_no_b_mass;
#endif
    // NB: hardcode cxtype cI(0,1) instead of cxtype (or hardcoded cxsmpl) mdl_complexi (which exists in Parameters_loop_sm_no_b_mass) because:
    // (1) mdl_complexi is always (0,1); (2) mdl_complexi is undefined in device code; (3) need cxsmpl conversion to cxtype in code below
    const cxtype cI( 0., 1. );
    DependentCouplings_sv out;
    // Begin non-SM (e.g. EFT) implementation - special handling of vectors of floats (#439)
#if not( defined MGONGPU_CPPSIMD && defined MGONGPU_FPTYPE_FLOAT )
    {
      const fptype_sv& G = G_sv;
      // Model parameters dependent on aS
      //const fptype_sv mdl_sqrt__aS = constexpr_sqrt( aS );
      //const fptype_sv G = 2. * mdl_sqrt__aS * constexpr_sqrt( M_PI );
      const fptype_sv mdl_G__exp__4 = ( ( G ) * ( G ) * ( G ) * ( G ) );
      const fptype_sv mdl_RGR2_FIN_ = -( 3. / 2. ) * mdl_G__exp__4 / ( 96. * ( ( M_PI ) * ( M_PI ) ) );
      const fptype_sv mdl_G__exp__2 = ( ( G ) * ( G ) );
      const fptype_sv mdl_R2MixedFactor_FIN_ = -( mdl_G__exp__2 * ( 1. + mdl_lhv ) * ( mdl_Ncol__exp__2 - 1. ) ) / ( 2. * mdl_Ncol * 16. * ( ( M_PI ) * ( M_PI ) ) );
      const fptype_sv mdl_G_UVg_1EPS_ = -( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 11. * mdl_CA;
      const fptype_sv mdl_G_UVq_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF;
      const fptype_sv mdl_G_UVc_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF;
      const fptype_sv mdl_G_UVb_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF;
      const fptype_sv mdl_G_UVt_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF;
      const fptype_sv mdl_GWcft_UV_t_1EPS_ = COND( mdl_MT, 0., -( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF );
      const fptype_sv mdl_tWcft_UV_1EPS_ = COND( mdl_MT, 0., -( ( mdl_G__exp__2 ) / ( 2. * 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF );
      constexpr cxsmpl<double> mdl_tMass_UV_1EPS_ = COND( mdl_MT, 0., mdl_complexi * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * mdl_MT );
      const fptype_sv mdl_UV_yuk_c_1EPS_ = -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * 2.;
      const fptype_sv mdl_UV_yuk_b_1EPS_ = -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * 2.;
      const fptype_sv mdl_UV_yuk_t_1EPS_ = -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * 2.;
      const fptype_sv mdl_R2_GGGpGm_factor_FIN_ = mdl_G__exp__2 / ( 8. * ( ( M_PI ) * ( M_PI ) ) );
      const fptype_sv mdl_R2_GGG0G0_factor_FIN_ = mdl_G__exp__2 / ( 8. * ( ( M_PI ) * ( M_PI ) ) );
      const fptype_sv mdl_G__exp__3 = ( ( G ) * ( G ) * ( G ) );
      const fptype_sv mdl_G_UVt_FIN_ = COND( mdl_MT, 0., -( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) );
      const fptype_sv mdl_GWcft_UV_t_FIN_ = COND( mdl_MT, 0., ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) );
      const fptype_sv mdl_tWcft_UV_FIN_ = COND( mdl_MT, 0., -( ( mdl_G__exp__2 ) / ( 2. * 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * mdl_CF * ( 4. - 3. * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) ) );
      constexpr cxsmpl<double> mdl_tMass_UV_FIN_ = COND( mdl_MT, 0., mdl_complexi * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * mdl_CF * ( 4. - 3. * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) ) * mdl_MT );
      const fptype_sv mdl_UV_yuk_t_FIN_ = COND( mdl_MT, 0., -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * mdl_CF * ( -3. * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) + 4. ) * 2. );
      // Model couplings dependent on aS
      out.GC_4 = -G;
      out.GC_5 = cI * G;
    }
#else
    // ** NB #439: special handling is necessary ONLY FOR VECTORS OF FLOATS (variable Gs are vector floats, fixed parameters are scalar doubles)
    // Use an explicit loop to avoid <<error: conversion of scalar ‘double’ to vector ‘fptype_sv’ {aka ‘__vector(8) float’} involves truncation>>
    // Problems may come e.g. in EFTs from multiplying a vector float (related to aS-dependent G) by a scalar double (aS-independent parameters)
    fptype_v GC_4r_v;
    fptype_v GC_4i_v;
    fptype_v GC_5r_v;
    fptype_v GC_5i_v;
    for( int i = 0; i < neppV; i++ )
    {
      const fptype& G = G_sv[i];
      // Model parameters dependent on aS
      //const fptype mdl_sqrt__aS = constexpr_sqrt( aS );
      //const fptype G = 2. * mdl_sqrt__aS * constexpr_sqrt( M_PI );
      const fptype mdl_G__exp__4 = ( ( G ) * ( G ) * ( G ) * ( G ) );
      const fptype mdl_RGR2_FIN_ = -( 3. / 2. ) * mdl_G__exp__4 / ( 96. * ( ( M_PI ) * ( M_PI ) ) );
      const fptype mdl_G__exp__2 = ( ( G ) * ( G ) );
      const fptype mdl_R2MixedFactor_FIN_ = -( mdl_G__exp__2 * ( 1. + mdl_lhv ) * ( mdl_Ncol__exp__2 - 1. ) ) / ( 2. * mdl_Ncol * 16. * ( ( M_PI ) * ( M_PI ) ) );
      const fptype mdl_G_UVg_1EPS_ = -( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 11. * mdl_CA;
      const fptype mdl_G_UVq_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF;
      const fptype mdl_G_UVc_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF;
      const fptype mdl_G_UVb_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF;
      const fptype mdl_G_UVt_1EPS_ = ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF;
      const fptype mdl_GWcft_UV_t_1EPS_ = COND( mdl_MT, 0., -( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF );
      const fptype mdl_tWcft_UV_1EPS_ = COND( mdl_MT, 0., -( ( mdl_G__exp__2 ) / ( 2. * 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF );
      constexpr cxsmpl<double> mdl_tMass_UV_1EPS_ = COND( mdl_MT, 0., mdl_complexi * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * mdl_MT );
      const fptype mdl_UV_yuk_c_1EPS_ = -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * 2.;
      const fptype mdl_UV_yuk_b_1EPS_ = -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * 2.;
      const fptype mdl_UV_yuk_t_1EPS_ = -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * 3. * mdl_CF * 2.;
      const fptype mdl_R2_GGGpGm_factor_FIN_ = mdl_G__exp__2 / ( 8. * ( ( M_PI ) * ( M_PI ) ) );
      const fptype mdl_R2_GGG0G0_factor_FIN_ = mdl_G__exp__2 / ( 8. * ( ( M_PI ) * ( M_PI ) ) );
      const fptype mdl_G__exp__3 = ( ( G ) * ( G ) * ( G ) );
      const fptype mdl_G_UVt_FIN_ = COND( mdl_MT, 0., -( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) );
      const fptype mdl_GWcft_UV_t_FIN_ = COND( mdl_MT, 0., ( ( mdl_G__exp__2 ) / ( 2. * 48. * ( ( M_PI ) * ( M_PI ) ) ) ) * 4. * mdl_TF * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) );
      const fptype mdl_tWcft_UV_FIN_ = COND( mdl_MT, 0., -( ( mdl_G__exp__2 ) / ( 2. * 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * mdl_CF * ( 4. - 3. * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) ) );
      constexpr cxsmpl<double> mdl_tMass_UV_FIN_ = COND( mdl_MT, 0., mdl_complexi * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * mdl_CF * ( 4. - 3. * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) ) * mdl_MT );
      const fptype mdl_UV_yuk_t_FIN_ = COND( mdl_MT, 0., -( 1. / 2. ) * ( ( mdl_G__exp__2 ) / ( 16. * ( ( M_PI ) * ( M_PI ) ) ) ) * mdl_CF * ( -3. * reglog( mdl_MT__exp__2 / mdl_MU_R__exp__2 ) + 4. ) * 2. );
      // Model couplings dependent on aS
      const cxtype GC_4 = -G;
      const cxtype GC_5 = cI * G;
      GC_4r_v[i] = cxreal( GC_4 );
      GC_4i_v[i] = cximag( GC_4 );
      GC_5r_v[i] = cxreal( GC_5 );
      GC_5i_v[i] = cximag( GC_5 );
    }
    out.GC_4 = cxtype_v( GC_4r_v, GC_4i_v );
    out.GC_5 = cxtype_v( GC_5r_v, GC_5i_v );
#endif
    // End non-SM (e.g. EFT) implementation - special handling of vectors of floats (#439)
    return out;
  }
#ifdef __CUDACC__
#pragma GCC diagnostic pop
#pragma nv_diagnostic pop
#endif
}

//==========================================================================

namespace Parameters_loop_sm_no_b_mass_independentCouplings
{
  constexpr size_t nicoup = 1; // #couplings that are fixed for all events because they do not depend on the running alphas QCD
  //constexpr size_t ixcoup_GC_11 = 0 + Parameters_loop_sm_no_b_mass_dependentCouplings::ndcoup; // out of ndcoup+nicoup
}

//==========================================================================

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
#pragma GCC diagnostic push
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wunused-but-set-variable" // e.g. <<warning: variable ‘couplings_sv’ set but not used [-Wunused-but-set-variable]>>
#endif
  // Compute the output couplings (e.g. gc10 and gc11) from the input gs
  template<class G_ACCESS, class C_ACCESS>
  __device__ inline void
  G2COUP( const fptype gs[],
          fptype couplings[] )
  {
    mgDebug( 0, __FUNCTION__ );
    using namespace Parameters_loop_sm_no_b_mass_dependentCouplings;
    const fptype_sv& gs_sv = G_ACCESS::kernelAccessConst( gs );
    DependentCouplings_sv couplings_sv = computeDependentCouplings_fromG( gs_sv );
    fptype* GC_4s = C_ACCESS::idcoupAccessBuffer( couplings, idcoup_GC_4 );
    fptype* GC_5s = C_ACCESS::idcoupAccessBuffer( couplings, idcoup_GC_5 );
    cxtype_sv_ref GC_4s_sv = C_ACCESS::kernelAccess( GC_4s );
    cxtype_sv_ref GC_5s_sv = C_ACCESS::kernelAccess( GC_5s );
    GC_4s_sv = couplings_sv.GC_4;
    GC_5s_sv = couplings_sv.GC_5;
    mgDebug( 1, __FUNCTION__ );
    return;
  }
#pragma GCC diagnostic pop
}

//==========================================================================

#endif // Parameters_loop_sm_no_b_mass_H
