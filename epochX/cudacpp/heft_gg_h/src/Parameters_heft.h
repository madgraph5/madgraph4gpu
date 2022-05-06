//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.3.1_lo_vect, 2022-01-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_heft_H
#define Parameters_heft_H

#include "mgOnGpuCxtypes.h"

#ifndef MGONGPU_HARDCODE_CIPC

#include "read_slha.h"

class Parameters_heft
{
public:

  static Parameters_heft* getInstance();

  // Define "zero"
  double zero, ZERO;

  // Model parameters independent of aS
  double mdl_WH1, mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau, mdl_ymt, mdl_ymb, aS, mdl_Gf, aEWM1, mdl_MP, mdl_MH, mdl_MZ, mdl_MTA, mdl_MT, mdl_MB, mdl_conjg__CKM3x3, mdl_MZ__exp__2, mdl_MZ__exp__4, mdl_sqrt__2, mdl_MH__exp__4, mdl_MT__exp__4, mdl_MH__exp__2, mdl_MT__exp__2, mdl_MH__exp__12, mdl_MH__exp__10, mdl_MH__exp__8, mdl_MH__exp__6, mdl_MT__exp__6, mdl_aEW, mdl_MW, mdl_sqrt__aEW, mdl_ee, mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw, mdl_v, mdl_ee__exp__2, mdl_MW__exp__12, mdl_MW__exp__10, mdl_MW__exp__8, mdl_MW__exp__6, mdl_MW__exp__4, mdl_AH, mdl_v__exp__2, mdl_lam, mdl_yb, mdl_yt, mdl_ytau, mdl_muH, mdl_gw__exp__2, mdl_cw__exp__2, mdl_sw__exp__2;
  cxsmpl<double> mdl_complexi;

  // Model couplings independent of aS
  // (none)

  // Model parameters dependent on aS
  double mdl_sqrt__aS, G, mdl_G__exp__2, mdl_GH, mdl_Gphi;

  // Model couplings dependent on aS
  cxsmpl<double> GC_13;

  // Set parameters that are unchanged during the run
  void setIndependentParameters( SLHAReader& slha );

  // Set couplings that are unchanged during the run
  void setIndependentCouplings();

  // Set parameters that are changed event by event
  void setDependentParameters();

  // Set couplings that are changed event by event
  void setDependentCouplings();

  // Print parameters that are unchanged during the run
  void printIndependentParameters();

  // Print couplings that are unchanged during the run
  void printIndependentCouplings();

  // Print parameters that are changed event by event
  void printDependentParameters();

  // Print couplings that are changed event by event
  void printDependentCouplings();

private:

  static Parameters_heft* instance;
};

#else

#include <limits>

// Hardcoded constexpr physics parameters
namespace Parameters_heft // keep the same name rather than HardcodedParameters_heft for simplicity
{
  // Constexpr implementation of sqrt (see https://stackoverflow.com/a/34134071)
  double constexpr detailSqrtNewtonRaphson( double x, double curr, double prev )
  {
    return curr == prev ? curr : detailSqrtNewtonRaphson( x, 0.5 * ( curr + x / curr ), curr );
  }
  double constexpr sqrtNR( double x )
  {
    return x >= 0 // && x < std::numeric_limits<double>::infinity() // avoid -Wtautological-constant-compare warning in fast math
      ? detailSqrtNewtonRaphson( x, x, 0 )
      : std::numeric_limits<double>::quiet_NaN();
  }

  // Constexpr implementation of pow
  template <typename T>
  constexpr T constexpr_pow( T base, unsigned int iexp )
  {
    return iexp == 0 ? 1 : base * constexpr_pow( base, iexp - 1 );
  }

  // Model parameters independent of aS
  constexpr double zero = 0;
  constexpr double ZERO = 0;
  constexpr double mdl_WH1 = 6.382339e-03;
  constexpr double mdl_WH = 6.382339e-03;
  constexpr double mdl_WW = 2.047600e+00;
  constexpr double mdl_WZ = 2.441404e+00;
  constexpr double mdl_WT = 1.491500e+00;
  constexpr double mdl_ymtau = 1.777000e+00;
  constexpr double mdl_ymt = 1.645000e+02;
  constexpr double mdl_ymb = 4.200000e+00;
  constexpr double aS = 1.180000e-01;
  constexpr double mdl_Gf = 1.166390e-05;
  constexpr double aEWM1 = 1.325070e+02;
  constexpr double mdl_MP = 1.250001e+02;
  constexpr double mdl_MH = 1.250000e+02;
  constexpr double mdl_MZ = 9.118800e+01;
  constexpr double mdl_MTA = 1.777000e+00;
  constexpr double mdl_MT = 1.730000e+02;
  constexpr double mdl_MB = 4.700000e+00;
  constexpr double mdl_conjg__CKM3x3 = 1.;
  constexpr cxsmpl<double> mdl_complexi = cxsmpl<double>( 0., 1. );
  constexpr double mdl_MZ__exp__2 = ( ( mdl_MZ ) * ( mdl_MZ ) );
  constexpr double mdl_MZ__exp__4 = ( ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) );
  constexpr double mdl_sqrt__2 = sqrtNR( 2. );
  constexpr double mdl_MH__exp__4 = ( ( mdl_MH ) * ( mdl_MH ) * ( mdl_MH ) * ( mdl_MH ) );
  constexpr double mdl_MT__exp__4 = ( ( mdl_MT ) * ( mdl_MT ) * ( mdl_MT ) * ( mdl_MT ) );
  constexpr double mdl_MH__exp__2 = ( ( mdl_MH ) * ( mdl_MH ) );
  constexpr double mdl_MT__exp__2 = ( ( mdl_MT ) * ( mdl_MT ) );
  constexpr double mdl_MH__exp__12 = constexpr_pow( mdl_MH, 12 );
  constexpr double mdl_MH__exp__10 = constexpr_pow( mdl_MH, 10 );
  constexpr double mdl_MH__exp__8 = constexpr_pow( mdl_MH, 8 );
  constexpr double mdl_MH__exp__6 = constexpr_pow( mdl_MH, 6 );
  constexpr double mdl_MT__exp__6 = constexpr_pow( mdl_MT, 6 );
  constexpr double mdl_aEW = 1. / aEWM1;
  constexpr double mdl_MW = sqrtNR( mdl_MZ__exp__2 / 2. + sqrtNR( mdl_MZ__exp__4 / 4. - ( mdl_aEW * M_PI * mdl_MZ__exp__2 ) / ( mdl_Gf * mdl_sqrt__2 ) ) );
  constexpr double mdl_sqrt__aEW = sqrtNR( mdl_aEW );
  constexpr double mdl_ee = 2. * mdl_sqrt__aEW * sqrtNR( M_PI );
  constexpr double mdl_MW__exp__2 = ( ( mdl_MW ) * ( mdl_MW ) );
  constexpr double mdl_sw2 = 1. - mdl_MW__exp__2 / mdl_MZ__exp__2;
  constexpr double mdl_cw = sqrtNR( 1. - mdl_sw2 );
  constexpr double mdl_sqrt__sw2 = sqrtNR( mdl_sw2 );
  constexpr double mdl_sw = mdl_sqrt__sw2;
  constexpr double mdl_g1 = mdl_ee / mdl_cw;
  constexpr double mdl_gw = mdl_ee / mdl_sw;
  constexpr double mdl_v = ( 2. * mdl_MW * mdl_sw ) / mdl_ee;
  constexpr double mdl_ee__exp__2 = ( ( mdl_ee ) * ( mdl_ee ) );
  constexpr double mdl_MW__exp__12 = constexpr_pow( mdl_MW, 12 );
  constexpr double mdl_MW__exp__10 = constexpr_pow( mdl_MW, 10 );
  constexpr double mdl_MW__exp__8 = constexpr_pow( mdl_MW, 8 );
  constexpr double mdl_MW__exp__6 = constexpr_pow( mdl_MW, 6 );
  constexpr double mdl_MW__exp__4 = ( ( mdl_MW ) * ( mdl_MW ) * ( mdl_MW ) * ( mdl_MW ) );
  constexpr double mdl_AH = ( 47. * mdl_ee__exp__2 * ( 1. - ( 2. * mdl_MH__exp__4 ) / ( 987. * mdl_MT__exp__4 ) - ( 14. * mdl_MH__exp__2 ) / ( 705. * mdl_MT__exp__2 ) + ( 213. * mdl_MH__exp__12 ) / ( 2.634632e7 * mdl_MW__exp__12 ) + ( 5. * mdl_MH__exp__10 ) / ( 119756. * mdl_MW__exp__10 ) + ( 41. * mdl_MH__exp__8 ) / ( 180950. * mdl_MW__exp__8 ) + ( 87. * mdl_MH__exp__6 ) / ( 65800. * mdl_MW__exp__6 ) + ( 57. * mdl_MH__exp__4 ) / ( 6580. * mdl_MW__exp__4 ) + ( 33. * mdl_MH__exp__2 ) / ( 470. * mdl_MW__exp__2 ) ) ) / ( 72. * ( ( M_PI ) * ( M_PI ) ) * mdl_v );
  constexpr double mdl_v__exp__2 = ( ( mdl_v ) * ( mdl_v ) );
  constexpr double mdl_lam = mdl_MH__exp__2 / ( 2. * mdl_v__exp__2 );
  constexpr double mdl_yb = ( mdl_ymb * mdl_sqrt__2 ) / mdl_v;
  constexpr double mdl_yt = ( mdl_ymt * mdl_sqrt__2 ) / mdl_v;
  constexpr double mdl_ytau = ( mdl_ymtau * mdl_sqrt__2 ) / mdl_v;
  constexpr double mdl_muH = sqrtNR( mdl_lam * mdl_v__exp__2 );
  constexpr double mdl_gw__exp__2 = ( ( mdl_gw ) * ( mdl_gw ) );
  constexpr double mdl_cw__exp__2 = ( ( mdl_cw ) * ( mdl_cw ) );
  constexpr double mdl_sw__exp__2 = ( ( mdl_sw ) * ( mdl_sw ) );

  // Model couplings independent of aS
  // (none)

  // Model parameters dependent on aS
  constexpr double mdl_sqrt__aS = sqrtNR( aS );
  constexpr double G = 2. * mdl_sqrt__aS * sqrtNR( M_PI );
  constexpr double mdl_G__exp__2 = ( ( G ) * ( G ) );
  constexpr double mdl_GH = -( mdl_G__exp__2 * ( 1. + ( 13. * mdl_MH__exp__6 ) / ( 16800. * mdl_MT__exp__6 ) + mdl_MH__exp__4 / ( 168. * mdl_MT__exp__4 ) + ( 7. * mdl_MH__exp__2 ) / ( 120. * mdl_MT__exp__2 ) ) ) / ( 12. * ( ( M_PI ) * ( M_PI ) ) * mdl_v );
  constexpr double mdl_Gphi = -( mdl_G__exp__2 * ( 1. + mdl_MH__exp__6 / ( 560. * mdl_MT__exp__6 ) + mdl_MH__exp__4 / ( 90. * mdl_MT__exp__4 ) + mdl_MH__exp__2 / ( 12. * mdl_MT__exp__2 ) ) ) / ( 8. * ( ( M_PI ) * ( M_PI ) ) * mdl_v );

  // Model couplings dependent on aS
  constexpr cxsmpl<double> GC_13 = -( mdl_complexi * mdl_GH );

  // Print parameters that are unchanged during the run
  void printIndependentParameters();

  // Print couplings that are unchanged during the run
  void printIndependentCouplings();

  // Print parameters that are changed event by event
  void printDependentParameters();

  // Print couplings that are changed event by event
  void printDependentCouplings();
}

#endif

#endif // Parameters_heft_H
