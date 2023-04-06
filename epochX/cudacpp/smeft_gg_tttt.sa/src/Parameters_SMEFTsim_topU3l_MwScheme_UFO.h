//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.5.0_lo_vect, 2023-01-26
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_SMEFTsim_topU3l_MwScheme_UFO_H
#define Parameters_SMEFTsim_topU3l_MwScheme_UFO_H

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"
#include "mgOnGpuVectors.h"

//==========================================================================

#ifndef MGONGPU_HARDCODE_PARAM // this is only supported in SM processes (e.g. not in EFT models) for the moment (#439)
#error This non-SM physics process only supports MGONGPU_HARDCODE_PARAM builds (#439): please run "make HRDCOD=1"

#include "read_slha.h"

class Parameters_SMEFTsim_topU3l_MwScheme_UFO
{
public:

  static Parameters_SMEFTsim_topU3l_MwScheme_UFO* getInstance();

  // Define "zero"
  double zero, ZERO;

  // Model parameters independent of aS
  //double aS; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  double mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau, mdl_ymm, mdl_yme, mdl_ymt, mdl_ymb, mdl_ymc, mdl_yms, mdl_ymup, mdl_ymdo, mdl_linearPropCorrections, mdl_Gf, mdl_MW, mdl_LambdaSMEFT, mdl_cleQt3Im, mdl_cleQt1Im, mdl_cleju3Im, mdl_cleju1Im, mdl_clebQIm, mdl_cledjIm, mdl_ceBIm, mdl_ceWIm, mdl_ceHIm, mdl_cQtQb8Im, mdl_cQtQb1Im, mdl_cjtQd8Im, mdl_cjtQd1Im, mdl_cQujb8Im, mdl_cQujb1Im, mdl_cjuQb8Im, mdl_cjuQb1Im, mdl_cQtjd8Im, mdl_cQtjd1Im, mdl_cjujd81Im, mdl_cjujd11Im, mdl_cjujd8Im, mdl_cjujd1Im, mdl_cjQbd8Im, mdl_cjQbd1Im, mdl_cjQtu8Im, mdl_cjQtu1Im, mdl_cutbd8Im, mdl_cutbd1Im, mdl_cHtbIm, mdl_cHudIm, mdl_cbHIm, mdl_cdHIm, mdl_ctHIm, mdl_cuHIm, mdl_cbBIm, mdl_cdBIm, mdl_cbWIm, mdl_cdWIm, mdl_cbGIm, mdl_cdGIm, mdl_ctBIm, mdl_cuBIm, mdl_ctWIm, mdl_cuWIm, mdl_ctGIm, mdl_cuGIm, mdl_cHWBtil, mdl_cHBtil, mdl_cHWtil, mdl_cHGtil, mdl_cWtil, mdl_cGtil, mdl_cleQt3Re, mdl_cleju3Re, mdl_cleQt1Re, mdl_cleju1Re, mdl_clebQRe, mdl_cledjRe, mdl_cle, mdl_cbl, mdl_cld, mdl_ctl, mdl_clu, mdl_cQe, mdl_cje, mdl_cbe, mdl_ced, mdl_cte, mdl_ceu, mdl_cee, mdl_cQl3, mdl_cQl1, mdl_clj3, mdl_clj1, mdl_cll1, mdl_cll, mdl_cHe, mdl_cHl3, mdl_cHl1, mdl_ceBRe, mdl_ceWRe, mdl_ceHRe, mdl_cQtQb8Re, mdl_cQtQb1Re, mdl_cjtQd8Re, mdl_cjtQd1Re, mdl_cQujb8Re, mdl_cQujb1Re, mdl_cjuQb8Re, mdl_cjuQb1Re, mdl_cQtjd8Re, mdl_cQtjd1Re, mdl_cjujd81Re, mdl_cjujd11Re, mdl_cjujd8Re, mdl_cjujd1Re, mdl_cjQbd8Re, mdl_cjQbd1Re, mdl_cjQtu8Re, mdl_cjQtu1Re, mdl_cQb8, mdl_cQb1, mdl_cbj8, mdl_cbj1, mdl_cQd8, mdl_cQd1, mdl_cjd8, mdl_cjd1, mdl_cQt8, mdl_cQt1, mdl_ctj8, mdl_ctj1, mdl_cQu8, mdl_cju8, mdl_cQu1, mdl_cju1, mdl_cutbd8Re, mdl_cutbd1Re, mdl_cbu8, mdl_ctd8, mdl_ctb8, mdl_cud8, mdl_cbu1, mdl_ctd1, mdl_ctb1, mdl_cud1, mdl_cbd8, mdl_cbd1, mdl_cbb, mdl_cdd8, mdl_cdd1, mdl_ctu8, mdl_ctu1, mdl_ctt, mdl_cuu8, mdl_cuu1, mdl_cQQ8, mdl_cQQ1, mdl_cQj38, mdl_cQj31, mdl_cQj18, mdl_cQj11, mdl_cjj38, mdl_cjj31, mdl_cjj18, mdl_cjj11, mdl_cHtbRe, mdl_cHudRe, mdl_cHbq, mdl_cHd, mdl_cHt, mdl_cHu, mdl_cHQ3, mdl_cHj3, mdl_cHQ1, mdl_cHj1, mdl_cbBRe, mdl_cdBRe, mdl_cbWRe, mdl_cdWRe, mdl_cbGRe, mdl_cdGRe, mdl_ctBRe, mdl_cuBRe, mdl_ctWRe, mdl_cuWRe, mdl_ctGRe, mdl_cuGRe, mdl_cbHRe, mdl_cdHRe, mdl_ctHRe, mdl_cuHRe, mdl_cHWB, mdl_cHB, mdl_cHW, mdl_cHG, mdl_cHDD, mdl_cHbox, mdl_cH, mdl_cW, mdl_cG, mdl_MH, mdl_MZ, mdl_MTA, mdl_MMU, mdl_Me, mdl_MT, mdl_MB, mdl_MC, mdl_MS, mdl_MU, mdl_MD, mdl_MWsm, mdl_MW__exp__2, mdl_MZ__exp__2, mdl_sqrt__2, mdl_nb__2__exp__0_25, mdl_MH__exp__2, mdl_sth2, mdl_nb__10__exp___m_40, mdl_propCorr, mdl_MZ1, mdl_MH1, mdl_MT1, mdl_WZ1, mdl_WW1, mdl_WH1, mdl_WT1, mdl_cth, mdl_MW1, mdl_sqrt__sth2, mdl_sth, mdl_LambdaSMEFT__exp__2, mdl_MT__exp__2, mdl_MH__exp__6, mdl_MWsm__exp__6, mdl_MH__exp__4, mdl_MWsm__exp__4, mdl_MWsm__exp__2, mdl_MZ__exp__4, mdl_MZ__exp__6, mdl_cth__exp__2, mdl_sth__exp__2, mdl_MB__exp__2, mdl_MZ__exp__3, mdl_sth__exp__4, mdl_sth__exp__6, mdl_sth__exp__3, mdl_sth__exp__5, mdl_propCorr__exp__2, mdl_propCorr__exp__3, mdl_propCorr__exp__4, mdl_cth__exp__3, mdl_aEW, mdl_sqrt__Gf, mdl_vevhat, mdl_lam, mdl_sqrt__aEW, mdl_ee, mdl_yb, mdl_yc, mdl_ydo, mdl_ye, mdl_ym, mdl_ys, mdl_yt, mdl_ytau, mdl_yup, mdl_vevhat__exp__2, mdl_dGf, mdl_dkH, mdl_vevT, mdl_g1, mdl_gw, mdl_ee__exp__2, mdl_gHaa, mdl_gHza, mdl_dMZ2, mdl_dMH2, mdl_dgw, mdl_barlam, mdl_dWT, mdl_dWW, mdl_gwsh, mdl_vev, mdl_dg1, mdl_dWHc, mdl_dWHb, mdl_dWHta, mdl_dWZ, mdl_g1sh, mdl_ee__exp__3, mdl_vevhat__exp__3;
  cxsmpl<double> mdl_complexi, mdl_cuH, mdl_ctHH, mdl_cdH, mdl_cbH, mdl_cuG, mdl_ctG, mdl_cuW, mdl_ctW, mdl_cuB, mdl_ctB, mdl_cdG, mdl_cbG, mdl_cdW, mdl_cbW, mdl_cdB, mdl_cbBB, mdl_cHud, mdl_cHtb, mdl_cutbd1, mdl_cutbd8, mdl_cjQtu1, mdl_cjQtu8, mdl_cjQbd1, mdl_cjQbd8, mdl_cjujd1, mdl_cjujd8, mdl_cjujd11, mdl_cjujd81, mdl_cQtjd1, mdl_cQtjd8, mdl_cjuQb1, mdl_cjuQb8, mdl_cQujb1, mdl_cQujb8, mdl_cjtQd1, mdl_cjtQd8, mdl_cQtQb1, mdl_cQtQb8, mdl_ceH, mdl_ceW, mdl_ceB, mdl_cledj, mdl_clebQ, mdl_cleju1, mdl_cleju3, mdl_cleQt1, mdl_cleQt3, mdl_conjg__cbH, mdl_conjg__ctHH, mdl_yb0, mdl_yt0;

  // Model couplings independent of aS
  // (none)

  // Model parameters dependent on aS
  //double mdl_sqrt__aS, G, mdl_gHgg2, mdl_gHgg4, mdl_gHgg5, mdl_G__exp__2, mdl_gHgg1, mdl_gHgg3, mdl_dWH; // now computed event-by-event (running alphas #373)
  //cxsmpl<double> mdl_G__exp__3; // now computed event-by-event (running alphas #373)

  // Model couplings dependent on aS
  //cxsmpl<double> GC_6, GC_7, GC_8; // now computed event-by-event (running alphas #373)

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

  static Parameters_SMEFTsim_topU3l_MwScheme_UFO* instance;
};

#else

#include <cassert>
#include <limits>

// Hardcoded constexpr physics parameters
namespace Parameters_SMEFTsim_topU3l_MwScheme_UFO // keep the same name rather than HardcodedParameters_SMEFTsim_topU3l_MwScheme_UFO for simplicity
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
  constexpr double mdl_WH = 4.070000e - 03;
  constexpr double mdl_WW = 2.085000e + 00;
  constexpr double mdl_WZ = 2.495200e + 00;
  constexpr double mdl_WT = 1.330000e + 00;
  constexpr double mdl_ymtau = 1.777000e + 00;
  constexpr double mdl_ymm = 1.056600e - 01;
  constexpr double mdl_yme = 5.110000e - 04;
  constexpr double mdl_ymt = 1.727600e + 02;
  constexpr double mdl_ymb = 4.180000e + 00;
  constexpr double mdl_ymc = 1.270000e + 00;
  constexpr double mdl_yms = 9.300000e - 02;
  constexpr double mdl_ymup = 2.160000e - 03;
  constexpr double mdl_ymdo = 4.670000e - 03;
  constexpr double mdl_linearPropCorrections = 0.000000e + 00;
  //constexpr double aS = 1.179000e - 01; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  constexpr double mdl_Gf = 1.166379e - 05;
  constexpr double mdl_MW = 8.038700e + 01;
  constexpr double mdl_LambdaSMEFT = 1.000000e + 03;
  constexpr double mdl_cleQt3Im = 0.000000e + 00;
  constexpr double mdl_cleQt1Im = 0.000000e + 00;
  constexpr double mdl_cleju3Im = 0.000000e + 00;
  constexpr double mdl_cleju1Im = 0.000000e + 00;
  constexpr double mdl_clebQIm = 0.000000e + 00;
  constexpr double mdl_cledjIm = 0.000000e + 00;
  constexpr double mdl_ceBIm = 0.000000e + 00;
  constexpr double mdl_ceWIm = 0.000000e + 00;
  constexpr double mdl_ceHIm = 0.000000e + 00;
  constexpr double mdl_cQtQb8Im = 0.000000e + 00;
  constexpr double mdl_cQtQb1Im = 0.000000e + 00;
  constexpr double mdl_cjtQd8Im = 0.000000e + 00;
  constexpr double mdl_cjtQd1Im = 0.000000e + 00;
  constexpr double mdl_cQujb8Im = 0.000000e + 00;
  constexpr double mdl_cQujb1Im = 0.000000e + 00;
  constexpr double mdl_cjuQb8Im = 0.000000e + 00;
  constexpr double mdl_cjuQb1Im = 0.000000e + 00;
  constexpr double mdl_cQtjd8Im = 0.000000e + 00;
  constexpr double mdl_cQtjd1Im = 0.000000e + 00;
  constexpr double mdl_cjujd81Im = 0.000000e + 00;
  constexpr double mdl_cjujd11Im = 0.000000e + 00;
  constexpr double mdl_cjujd8Im = 0.000000e + 00;
  constexpr double mdl_cjujd1Im = 0.000000e + 00;
  constexpr double mdl_cjQbd8Im = 0.000000e + 00;
  constexpr double mdl_cjQbd1Im = 0.000000e + 00;
  constexpr double mdl_cjQtu8Im = 0.000000e + 00;
  constexpr double mdl_cjQtu1Im = 0.000000e + 00;
  constexpr double mdl_cutbd8Im = 0.000000e + 00;
  constexpr double mdl_cutbd1Im = 0.000000e + 00;
  constexpr double mdl_cHtbIm = 0.000000e + 00;
  constexpr double mdl_cHudIm = 0.000000e + 00;
  constexpr double mdl_cbHIm = 0.000000e + 00;
  constexpr double mdl_cdHIm = 0.000000e + 00;
  constexpr double mdl_ctHIm = 0.000000e + 00;
  constexpr double mdl_cuHIm = 0.000000e + 00;
  constexpr double mdl_cbBIm = 0.000000e + 00;
  constexpr double mdl_cdBIm = 0.000000e + 00;
  constexpr double mdl_cbWIm = 0.000000e + 00;
  constexpr double mdl_cdWIm = 0.000000e + 00;
  constexpr double mdl_cbGIm = 0.000000e + 00;
  constexpr double mdl_cdGIm = 0.000000e + 00;
  constexpr double mdl_ctBIm = 0.000000e + 00;
  constexpr double mdl_cuBIm = 0.000000e + 00;
  constexpr double mdl_ctWIm = 0.000000e + 00;
  constexpr double mdl_cuWIm = 0.000000e + 00;
  constexpr double mdl_ctGIm = 0.000000e + 00;
  constexpr double mdl_cuGIm = 0.000000e + 00;
  constexpr double mdl_cHWBtil = 0.000000e + 00;
  constexpr double mdl_cHBtil = 0.000000e + 00;
  constexpr double mdl_cHWtil = 0.000000e + 00;
  constexpr double mdl_cHGtil = 0.000000e + 00;
  constexpr double mdl_cWtil = 0.000000e + 00;
  constexpr double mdl_cGtil = 0.000000e + 00;
  constexpr double mdl_cleQt3Re = 0.000000e + 00;
  constexpr double mdl_cleju3Re = 0.000000e + 00;
  constexpr double mdl_cleQt1Re = 0.000000e + 00;
  constexpr double mdl_cleju1Re = 0.000000e + 00;
  constexpr double mdl_clebQRe = 0.000000e + 00;
  constexpr double mdl_cledjRe = 0.000000e + 00;
  constexpr double mdl_cle = 0.000000e + 00;
  constexpr double mdl_cbl = 0.000000e + 00;
  constexpr double mdl_cld = 0.000000e + 00;
  constexpr double mdl_ctl = 0.000000e + 00;
  constexpr double mdl_clu = 0.000000e + 00;
  constexpr double mdl_cQe = 0.000000e + 00;
  constexpr double mdl_cje = 0.000000e + 00;
  constexpr double mdl_cbe = 0.000000e + 00;
  constexpr double mdl_ced = 0.000000e + 00;
  constexpr double mdl_cte = 0.000000e + 00;
  constexpr double mdl_ceu = 0.000000e + 00;
  constexpr double mdl_cee = 0.000000e + 00;
  constexpr double mdl_cQl3 = 0.000000e + 00;
  constexpr double mdl_cQl1 = 0.000000e + 00;
  constexpr double mdl_clj3 = 0.000000e + 00;
  constexpr double mdl_clj1 = 0.000000e + 00;
  constexpr double mdl_cll1 = 0.000000e + 00;
  constexpr double mdl_cll = 0.000000e + 00;
  constexpr double mdl_cHe = 0.000000e + 00;
  constexpr double mdl_cHl3 = 0.000000e + 00;
  constexpr double mdl_cHl1 = 0.000000e + 00;
  constexpr double mdl_ceBRe = 0.000000e + 00;
  constexpr double mdl_ceWRe = 0.000000e + 00;
  constexpr double mdl_ceHRe = 0.000000e + 00;
  constexpr double mdl_cQtQb8Re = 0.000000e + 00;
  constexpr double mdl_cQtQb1Re = 0.000000e + 00;
  constexpr double mdl_cjtQd8Re = 0.000000e + 00;
  constexpr double mdl_cjtQd1Re = 0.000000e + 00;
  constexpr double mdl_cQujb8Re = 0.000000e + 00;
  constexpr double mdl_cQujb1Re = 0.000000e + 00;
  constexpr double mdl_cjuQb8Re = 0.000000e + 00;
  constexpr double mdl_cjuQb1Re = 0.000000e + 00;
  constexpr double mdl_cQtjd8Re = 0.000000e + 00;
  constexpr double mdl_cQtjd1Re = 0.000000e + 00;
  constexpr double mdl_cjujd81Re = 0.000000e + 00;
  constexpr double mdl_cjujd11Re = 0.000000e + 00;
  constexpr double mdl_cjujd8Re = 0.000000e + 00;
  constexpr double mdl_cjujd1Re = 0.000000e + 00;
  constexpr double mdl_cjQbd8Re = 0.000000e + 00;
  constexpr double mdl_cjQbd1Re = 0.000000e + 00;
  constexpr double mdl_cjQtu8Re = 0.000000e + 00;
  constexpr double mdl_cjQtu1Re = 0.000000e + 00;
  constexpr double mdl_cQb8 = 0.000000e + 00;
  constexpr double mdl_cQb1 = 0.000000e + 00;
  constexpr double mdl_cbj8 = 0.000000e + 00;
  constexpr double mdl_cbj1 = 0.000000e + 00;
  constexpr double mdl_cQd8 = 0.000000e + 00;
  constexpr double mdl_cQd1 = 0.000000e + 00;
  constexpr double mdl_cjd8 = 0.000000e + 00;
  constexpr double mdl_cjd1 = 0.000000e + 00;
  constexpr double mdl_cQt8 = 0.000000e + 00;
  constexpr double mdl_cQt1 = 0.000000e + 00;
  constexpr double mdl_ctj8 = 0.000000e + 00;
  constexpr double mdl_ctj1 = 0.000000e + 00;
  constexpr double mdl_cQu8 = 0.000000e + 00;
  constexpr double mdl_cju8 = 0.000000e + 00;
  constexpr double mdl_cQu1 = 0.000000e + 00;
  constexpr double mdl_cju1 = 0.000000e + 00;
  constexpr double mdl_cutbd8Re = 0.000000e + 00;
  constexpr double mdl_cutbd1Re = 0.000000e + 00;
  constexpr double mdl_cbu8 = 0.000000e + 00;
  constexpr double mdl_ctd8 = 0.000000e + 00;
  constexpr double mdl_ctb8 = 0.000000e + 00;
  constexpr double mdl_cud8 = 0.000000e + 00;
  constexpr double mdl_cbu1 = 0.000000e + 00;
  constexpr double mdl_ctd1 = 0.000000e + 00;
  constexpr double mdl_ctb1 = 0.000000e + 00;
  constexpr double mdl_cud1 = 0.000000e + 00;
  constexpr double mdl_cbd8 = 0.000000e + 00;
  constexpr double mdl_cbd1 = 0.000000e + 00;
  constexpr double mdl_cbb = 0.000000e + 00;
  constexpr double mdl_cdd8 = 0.000000e + 00;
  constexpr double mdl_cdd1 = 0.000000e + 00;
  constexpr double mdl_ctu8 = 0.000000e + 00;
  constexpr double mdl_ctu1 = 0.000000e + 00;
  constexpr double mdl_ctt = 0.000000e + 00;
  constexpr double mdl_cuu8 = 0.000000e + 00;
  constexpr double mdl_cuu1 = 0.000000e + 00;
  constexpr double mdl_cQQ8 = 0.000000e + 00;
  constexpr double mdl_cQQ1 = 0.000000e + 00;
  constexpr double mdl_cQj38 = 0.000000e + 00;
  constexpr double mdl_cQj31 = 0.000000e + 00;
  constexpr double mdl_cQj18 = 0.000000e + 00;
  constexpr double mdl_cQj11 = 0.000000e + 00;
  constexpr double mdl_cjj38 = 0.000000e + 00;
  constexpr double mdl_cjj31 = 0.000000e + 00;
  constexpr double mdl_cjj18 = 0.000000e + 00;
  constexpr double mdl_cjj11 = 0.000000e + 00;
  constexpr double mdl_cHtbRe = 0.000000e + 00;
  constexpr double mdl_cHudRe = 0.000000e + 00;
  constexpr double mdl_cHbq = 0.000000e + 00;
  constexpr double mdl_cHd = 0.000000e + 00;
  constexpr double mdl_cHt = 0.000000e + 00;
  constexpr double mdl_cHu = 0.000000e + 00;
  constexpr double mdl_cHQ3 = 0.000000e + 00;
  constexpr double mdl_cHj3 = 0.000000e + 00;
  constexpr double mdl_cHQ1 = 0.000000e + 00;
  constexpr double mdl_cHj1 = 0.000000e + 00;
  constexpr double mdl_cbBRe = 0.000000e + 00;
  constexpr double mdl_cdBRe = 0.000000e + 00;
  constexpr double mdl_cbWRe = 0.000000e + 00;
  constexpr double mdl_cdWRe = 0.000000e + 00;
  constexpr double mdl_cbGRe = 0.000000e + 00;
  constexpr double mdl_cdGRe = 0.000000e + 00;
  constexpr double mdl_ctBRe = 0.000000e + 00;
  constexpr double mdl_cuBRe = 0.000000e + 00;
  constexpr double mdl_ctWRe = 0.000000e + 00;
  constexpr double mdl_cuWRe = 0.000000e + 00;
  constexpr double mdl_ctGRe = 0.000000e + 00;
  constexpr double mdl_cuGRe = 0.000000e + 00;
  constexpr double mdl_cbHRe = 0.000000e + 00;
  constexpr double mdl_cdHRe = 0.000000e + 00;
  constexpr double mdl_ctHRe = 0.000000e + 00;
  constexpr double mdl_cuHRe = 0.000000e + 00;
  constexpr double mdl_cHWB = 0.000000e + 00;
  constexpr double mdl_cHB = 0.000000e + 00;
  constexpr double mdl_cHW = 0.000000e + 00;
  constexpr double mdl_cHG = 0.000000e + 00;
  constexpr double mdl_cHDD = 0.000000e + 00;
  constexpr double mdl_cHbox = 0.000000e + 00;
  constexpr double mdl_cH = 0.000000e + 00;
  constexpr double mdl_cW = 0.000000e + 00;
  constexpr double mdl_cG = 0.000000e + 00;
  constexpr double mdl_MH = 1.250900e + 02;
  constexpr double mdl_MZ = 9.118760e + 01;
  constexpr double mdl_MTA = 1.777000e + 00;
  constexpr double mdl_MMU = 1.056600e - 01;
  constexpr double mdl_Me = 5.110000e - 04;
  constexpr double mdl_MT = 1.727600e + 02;
  constexpr double mdl_MB = 4.180000e + 00;
  constexpr double mdl_MC = 1.270000e + 00;
  constexpr double mdl_MS = 9.300000e - 02;
  constexpr double mdl_MU = 2.160000e - 03;
  constexpr double mdl_MD = 4.670000e - 03;
  constexpr cxsmpl<double> mdl_complexi = cxsmpl<double>( 0., 1. );
  constexpr cxsmpl<double> mdl_cuH = mdl_cuHRe + mdl_cuHIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_ctHH = mdl_ctHRe + mdl_ctHIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cdH = mdl_cdHRe + mdl_cdHIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cbH = mdl_cbHRe + mdl_cbHIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cuG = mdl_cuGRe + mdl_cuGIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_ctG = mdl_ctGRe + mdl_ctGIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cuW = mdl_cuWRe + mdl_cuWIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_ctW = mdl_ctWRe + mdl_ctWIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cuB = mdl_cuBRe + mdl_cuBIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_ctB = mdl_ctBRe + mdl_ctBIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cdG = mdl_cdGRe + mdl_cdGIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cbG = mdl_cbGRe + mdl_cbGIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cdW = mdl_cdWRe + mdl_cdWIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cbW = mdl_cbWRe + mdl_cbWIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cdB = mdl_cdBRe + mdl_cdBIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cbBB = mdl_cbBRe + mdl_cbBIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cHud = mdl_cHudRe + mdl_cHudIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cHtb = mdl_cHtbRe + mdl_cHtbIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cutbd1 = mdl_cutbd1Re + mdl_cutbd1Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cutbd8 = mdl_cutbd8Re + mdl_cutbd8Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cjQtu1 = mdl_cjQtu1Re + mdl_cjQtu1Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cjQtu8 = mdl_cjQtu8Re + mdl_cjQtu8Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cjQbd1 = mdl_cjQbd1Re + mdl_cjQbd1Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cjQbd8 = mdl_cjQbd8Re + mdl_cjQbd8Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cjujd1 = mdl_cjujd1Re + mdl_cjujd1Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cjujd8 = mdl_cjujd8Re + mdl_cjujd8Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cjujd11 = mdl_cjujd11Re + mdl_cjujd11Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cjujd81 = mdl_cjujd81Re + mdl_cjujd81Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cQtjd1 = mdl_cQtjd1Re + mdl_cQtjd1Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cQtjd8 = mdl_cQtjd8Re + mdl_cQtjd8Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cjuQb1 = mdl_cjuQb1Re + mdl_cjuQb1Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cjuQb8 = mdl_cjuQb8Re + mdl_cjuQb8Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cQujb1 = mdl_cQujb1Re + mdl_cQujb1Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cQujb8 = mdl_cQujb8Re + mdl_cQujb8Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cjtQd1 = mdl_cjtQd1Re + mdl_cjtQd1Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cjtQd8 = mdl_cjtQd8Re + mdl_cjtQd8Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cQtQb1 = mdl_cQtQb1Re + mdl_cQtQb1Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cQtQb8 = mdl_cQtQb8Re + mdl_cQtQb8Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_ceH = mdl_ceHRe + mdl_ceHIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_ceW = mdl_ceWRe + mdl_ceWIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_ceB = mdl_ceBRe + mdl_ceBIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cledj = mdl_cledjRe + mdl_cledjIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_clebQ = mdl_clebQRe + mdl_clebQIm * mdl_complexi;
  constexpr cxsmpl<double> mdl_cleju1 = mdl_cleju1Re + mdl_cleju1Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cleju3 = mdl_cleju3Re + mdl_cleju3Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cleQt1 = mdl_cleQt1Re + mdl_cleQt1Im * mdl_complexi;
  constexpr cxsmpl<double> mdl_cleQt3 = mdl_cleQt3Re + mdl_cleQt3Im * mdl_complexi;
  constexpr double mdl_MWsm = mdl_MW;
  constexpr double mdl_MW__exp__2 = ( ( mdl_MW ) * ( mdl_MW ) );
  constexpr double mdl_MZ__exp__2 = ( ( mdl_MZ ) * ( mdl_MZ ) );
  constexpr double mdl_sqrt__2 = constexpr_sqrt( 2. );
  constexpr double mdl_nb__2__exp__0_25 = constexpr_pow( 2., 0.25 );
  constexpr double mdl_MH__exp__2 = ( ( mdl_MH ) * ( mdl_MH ) );
  constexpr double mdl_sth2 = 1. - mdl_MW__exp__2 / mdl_MZ__exp__2;
  constexpr double mdl_nb__10__exp___m_40 = constexpr_pow( 10., -40. );
  constexpr double mdl_propCorr = ABS( mdl_linearPropCorrections ) / ( ABS( mdl_linearPropCorrections ) + mdl_nb__10__exp___m_40 );
  constexpr double mdl_MZ1 = mdl_MZ;
  constexpr double mdl_MH1 = mdl_MH;
  constexpr double mdl_MT1 = mdl_MT;
  constexpr double mdl_WZ1 = mdl_WZ;
  constexpr double mdl_WW1 = mdl_WW;
  constexpr double mdl_WH1 = mdl_WH;
  constexpr double mdl_WT1 = mdl_WT;
  constexpr double mdl_cth = constexpr_sqrt( 1. - mdl_sth2 );
  constexpr double mdl_MW1 = mdl_MWsm;
  constexpr double mdl_sqrt__sth2 = constexpr_sqrt( mdl_sth2 );
  constexpr double mdl_sth = mdl_sqrt__sth2;
  constexpr double mdl_LambdaSMEFT__exp__2 = ( ( mdl_LambdaSMEFT ) * ( mdl_LambdaSMEFT ) );
  constexpr cxsmpl<double> mdl_conjg__cbH = conj( mdl_cbH );
  constexpr cxsmpl<double> mdl_conjg__ctHH = conj( mdl_ctHH );
  constexpr double mdl_MT__exp__2 = ( ( mdl_MT ) * ( mdl_MT ) );
  constexpr double mdl_MH__exp__6 = constexpr_pow( mdl_MH, 6. );
  constexpr double mdl_MWsm__exp__6 = constexpr_pow( mdl_MWsm, 6. );
  constexpr double mdl_MH__exp__4 = ( ( mdl_MH ) * ( mdl_MH ) * ( mdl_MH ) * ( mdl_MH ) );
  constexpr double mdl_MWsm__exp__4 = ( ( mdl_MWsm ) * ( mdl_MWsm ) * ( mdl_MWsm ) * ( mdl_MWsm ) );
  constexpr double mdl_MWsm__exp__2 = ( ( mdl_MWsm ) * ( mdl_MWsm ) );
  constexpr double mdl_MZ__exp__4 = ( ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) );
  constexpr double mdl_MZ__exp__6 = constexpr_pow( mdl_MZ, 6. );
  constexpr double mdl_cth__exp__2 = ( ( mdl_cth ) * ( mdl_cth ) );
  constexpr double mdl_sth__exp__2 = ( ( mdl_sth ) * ( mdl_sth ) );
  constexpr double mdl_MB__exp__2 = ( ( mdl_MB ) * ( mdl_MB ) );
  constexpr double mdl_MZ__exp__3 = ( ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) );
  constexpr double mdl_sth__exp__4 = ( ( mdl_sth ) * ( mdl_sth ) * ( mdl_sth ) * ( mdl_sth ) );
  constexpr double mdl_sth__exp__6 = constexpr_pow( mdl_sth, 6. );
  constexpr double mdl_sth__exp__3 = ( ( mdl_sth ) * ( mdl_sth ) * ( mdl_sth ) );
  constexpr double mdl_sth__exp__5 = constexpr_pow( mdl_sth, 5. );
  constexpr double mdl_propCorr__exp__2 = ( ( mdl_propCorr ) * ( mdl_propCorr ) );
  constexpr double mdl_propCorr__exp__3 = ( ( mdl_propCorr ) * ( mdl_propCorr ) * ( mdl_propCorr ) );
  constexpr double mdl_propCorr__exp__4 = ( ( mdl_propCorr ) * ( mdl_propCorr ) * ( mdl_propCorr ) * ( mdl_propCorr ) );
  constexpr double mdl_cth__exp__3 = ( ( mdl_cth ) * ( mdl_cth ) * ( mdl_cth ) );
  constexpr double mdl_aEW = ( mdl_Gf * mdl_MW__exp__2 * ( 1. - mdl_MW__exp__2 / mdl_MZ__exp__2 ) * mdl_sqrt__2 ) / M_PI;
  constexpr double mdl_sqrt__Gf = constexpr_sqrt( mdl_Gf );
  constexpr double mdl_vevhat = 1. / ( mdl_nb__2__exp__0_25 * mdl_sqrt__Gf );
  constexpr double mdl_lam = ( mdl_Gf * mdl_MH__exp__2 ) / mdl_sqrt__2;
  constexpr double mdl_sqrt__aEW = constexpr_sqrt( mdl_aEW );
  constexpr double mdl_ee = 2. * mdl_sqrt__aEW * constexpr_sqrt( M_PI );
  constexpr double mdl_yb = ( mdl_ymb * mdl_sqrt__2 ) / mdl_vevhat;
  constexpr double mdl_yc = ( mdl_ymc * mdl_sqrt__2 ) / mdl_vevhat;
  constexpr double mdl_ydo = ( mdl_ymdo * mdl_sqrt__2 ) / mdl_vevhat;
  constexpr double mdl_ye = ( mdl_yme * mdl_sqrt__2 ) / mdl_vevhat;
  constexpr double mdl_ym = ( mdl_ymm * mdl_sqrt__2 ) / mdl_vevhat;
  constexpr double mdl_ys = ( mdl_yms * mdl_sqrt__2 ) / mdl_vevhat;
  constexpr double mdl_yt = ( mdl_ymt * mdl_sqrt__2 ) / mdl_vevhat;
  constexpr double mdl_ytau = ( mdl_ymtau * mdl_sqrt__2 ) / mdl_vevhat;
  constexpr double mdl_yup = ( mdl_ymup * mdl_sqrt__2 ) / mdl_vevhat;
  constexpr double mdl_vevhat__exp__2 = ( ( mdl_vevhat ) * ( mdl_vevhat ) );
  constexpr double mdl_dGf = ( ( 2. * mdl_cHl3 - mdl_cll1 ) * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2;
  constexpr double mdl_dkH = ( ( mdl_cHbox - mdl_cHDD / 4. ) * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2;
  constexpr double mdl_vevT = ( 1. + mdl_dGf / 2. ) * mdl_vevhat;
  constexpr double mdl_g1 = mdl_ee / mdl_cth;
  constexpr double mdl_gw = mdl_ee / mdl_sth;
  constexpr cxsmpl<double> mdl_yb0 = ( 1. - mdl_dGf / 2. ) * mdl_yb + ( mdl_vevhat__exp__2 * mdl_conjg__cbH ) / ( 2. * mdl_LambdaSMEFT__exp__2 );
  constexpr cxsmpl<double> mdl_yt0 = ( 1. - mdl_dGf / 2. ) * mdl_yt + ( mdl_vevhat__exp__2 * mdl_conjg__ctHH ) / ( 2. * mdl_LambdaSMEFT__exp__2 );
  constexpr double mdl_ee__exp__2 = ( ( mdl_ee ) * ( mdl_ee ) );
  constexpr double mdl_gHaa = ( mdl_ee__exp__2 * ( -1.75 + ( 4. * ( 0.3333333333333333 + ( 7. * mdl_MH__exp__2 ) / ( 360. * mdl_MT__exp__2 ) ) ) / 3. - ( 29. * mdl_MH__exp__6 ) / ( 16800. * mdl_MWsm__exp__6 ) - ( 19. * mdl_MH__exp__4 ) / ( 1680. * mdl_MWsm__exp__4 ) - ( 11. * mdl_MH__exp__2 ) / ( 120. * mdl_MWsm__exp__2 ) ) ) / ( 8. * ( ( M_PI ) * ( M_PI ) ) );
  constexpr double mdl_gHza = ( mdl_ee__exp__2 * ( ( ( 0.4583333333333333 + ( 29. * mdl_MH__exp__6 ) / ( 100800. * mdl_MWsm__exp__6 ) + ( 19. * mdl_MH__exp__4 ) / ( 10080. * mdl_MWsm__exp__4 ) + ( 11. * mdl_MH__exp__2 ) / ( 720. * mdl_MWsm__exp__2 ) + ( mdl_MH__exp__4 * mdl_MZ__exp__2 ) / ( 2100. * mdl_MWsm__exp__6 ) + ( mdl_MH__exp__2 * mdl_MZ__exp__2 ) / ( 280. * mdl_MWsm__exp__4 ) + ( 7. * mdl_MZ__exp__2 ) / ( 180. * mdl_MWsm__exp__2 ) + ( 67. * mdl_MH__exp__2 * mdl_MZ__exp__4 ) / ( 100800. * mdl_MWsm__exp__6 ) + ( 53. * mdl_MZ__exp__4 ) / ( 10080. * mdl_MWsm__exp__4 ) + ( 43. * mdl_MZ__exp__6 ) / ( 50400. * mdl_MWsm__exp__6 ) - ( 31. * mdl_cth__exp__2 ) / ( 24. * mdl_sth__exp__2 ) - ( 29. * mdl_cth__exp__2 * mdl_MH__exp__6 ) / ( 20160. * mdl_MWsm__exp__6 * mdl_sth__exp__2 ) - ( 19. * mdl_cth__exp__2 * mdl_MH__exp__4 ) / ( 2016. * mdl_MWsm__exp__4 * mdl_sth__exp__2 ) - ( 11. * mdl_cth__exp__2 * mdl_MH__exp__2 ) / ( 144. * mdl_MWsm__exp__2 * mdl_sth__exp__2 ) - ( mdl_cth__exp__2 * mdl_MH__exp__4 * mdl_MZ__exp__2 ) / ( 560. * mdl_MWsm__exp__6 * mdl_sth__exp__2 ) - ( 31. * mdl_cth__exp__2 * mdl_MH__exp__2 * mdl_MZ__exp__2 ) / ( 2520. * mdl_MWsm__exp__4 * mdl_sth__exp__2 ) - ( mdl_cth__exp__2 * mdl_MZ__exp__2 ) / ( 9. * mdl_MWsm__exp__2 * mdl_sth__exp__2 ) - ( 43. * mdl_cth__exp__2 * mdl_MH__exp__2 * mdl_MZ__exp__4 ) / ( 20160. * mdl_MWsm__exp__6 * mdl_sth__exp__2 ) - ( 17. * mdl_cth__exp__2 * mdl_MZ__exp__4 ) / ( 1120. * mdl_MWsm__exp__4 * mdl_sth__exp__2 ) - ( 5. * mdl_cth__exp__2 * mdl_MZ__exp__6 ) / ( 2016. * mdl_MWsm__exp__6 * mdl_sth__exp__2 ) ) * mdl_sth ) / mdl_cth + ( ( 0.3333333333333333 + ( 7. * mdl_MH__exp__2 ) / ( 360. * mdl_MT__exp__2 ) + ( 11. * mdl_MZ__exp__2 ) / ( 360. * mdl_MT__exp__2 ) ) * ( 0.5 - ( 4. * mdl_sth__exp__2 ) / 3. ) ) / ( mdl_cth * mdl_sth ) ) ) / ( 4. * ( ( M_PI ) * ( M_PI ) ) );
  constexpr double mdl_dMZ2 = ( ( mdl_cHDD / 2. + 2. * mdl_cHWB * mdl_cth * mdl_sth ) * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2;
  constexpr double mdl_dMH2 = 2. * mdl_dkH - ( 3. * mdl_cH * mdl_vevhat__exp__2 ) / ( 2. * mdl_lam * mdl_LambdaSMEFT__exp__2 );
  constexpr double mdl_dgw = -mdl_dGf / 2.;
  constexpr double mdl_barlam = ( 1. - mdl_dGf - mdl_dMH2 ) * mdl_lam;
  constexpr double mdl_dWT = 2. * mdl_WT * ( mdl_dgw + ( mdl_vevhat * ( mdl_ee * ( 3. * mdl_cHtbRe * mdl_MB * mdl_MT * mdl_MWsm__exp__2 + mdl_cHQ3 * ( ( ( mdl_MB__exp__2 - mdl_MT__exp__2 ) * ( mdl_MB__exp__2 - mdl_MT__exp__2 ) ) + ( mdl_MB__exp__2 + mdl_MT__exp__2 ) * mdl_MWsm__exp__2 - 2. * mdl_MWsm__exp__4 ) ) * mdl_vevhat + 6. * mdl_MWsm__exp__2 * ( mdl_ctWRe * mdl_MT * ( mdl_MB__exp__2 - mdl_MT__exp__2 + mdl_MWsm__exp__2 ) + mdl_cbWRe * mdl_MB * ( -mdl_MB__exp__2 + mdl_MT__exp__2 + mdl_MWsm__exp__2 ) ) * mdl_sth * mdl_sqrt__2 ) ) / ( mdl_ee * mdl_LambdaSMEFT__exp__2 * ( ( ( mdl_MB__exp__2 - mdl_MT__exp__2 ) * ( mdl_MB__exp__2 - mdl_MT__exp__2 ) ) + ( mdl_MB__exp__2 + mdl_MT__exp__2 ) * mdl_MWsm__exp__2 - 2. * mdl_MWsm__exp__4 ) ) );
  constexpr double mdl_dWW = ( 2. * mdl_dgw + ( 2. * ( 2. * mdl_cHj3 + mdl_cHl3 ) * mdl_vevhat__exp__2 ) / ( 3. * mdl_LambdaSMEFT__exp__2 ) ) * mdl_WW;
  constexpr double mdl_gwsh = ( mdl_ee * ( 1. + mdl_dgw - ( mdl_cHW * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 ) ) / mdl_sth;
  constexpr double mdl_vev = ( 1. - ( 3. * mdl_cH * mdl_vevhat__exp__2 ) / ( 8. * mdl_lam * mdl_LambdaSMEFT__exp__2 ) ) * mdl_vevT;
  constexpr double mdl_dg1 = ( -mdl_dGf - mdl_dMZ2 / mdl_sth__exp__2 ) / 2.;
  constexpr double mdl_dWHc = mdl_yc / ( mdl_yc + mdl_nb__10__exp___m_40 ) * ( -0.02884 * mdl_dGf + ( ( 0.05768 * mdl_cHbox - 0.01442 * mdl_cHDD - 0.05768 * mdl_cuHRe ) * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 );
  constexpr double mdl_dWHb = mdl_yb / ( mdl_yb + mdl_nb__10__exp___m_40 ) * ( mdl_vevhat__exp__2 * ( -1.1618 * mdl_cbHRe ) / ( mdl_LambdaSMEFT__exp__2 * ( mdl_yb + mdl_nb__10__exp___m_40 ) ) - 0.5809 * mdl_dGf + ( mdl_vevhat__exp__2 * ( 1.1618 * mdl_cHbox - 0.29045 * mdl_cHDD ) ) / ( mdl_LambdaSMEFT__exp__2 ) );
  constexpr double mdl_dWHta = mdl_ytau / ( mdl_ytau + mdl_nb__10__exp___m_40 ) * ( -0.06256 * mdl_dGf + mdl_vevhat__exp__2 * ( -0.12512 * mdl_ceHRe + 0.12512 * mdl_cHbox - 0.03128 * mdl_cHDD ) / ( mdl_LambdaSMEFT__exp__2 ) );
  constexpr double mdl_dWZ = mdl_WZ * ( -1. + ( 36. * mdl_cth * mdl_MB * mdl_MZ__exp__2 * mdl_sth * ( mdl_cbWRe * mdl_cth + mdl_cbBRe * mdl_sth ) * ( -3. + 4. * mdl_sth__exp__2 ) * mdl_vevhat * mdl_sqrt__2 * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + mdl_ee * mdl_LambdaSMEFT__exp__2 * ( 2. * mdl_MZ__exp__3 * ( 27. + 54. * mdl_dgw - 54. * ( 1. + mdl_dg1 + mdl_dgw ) * mdl_sth__exp__2 + 76. * ( 1. + 4. * mdl_dg1 - 2. * mdl_dgw ) * mdl_sth__exp__4 + 152. * ( -mdl_dg1 + mdl_dgw ) * mdl_sth__exp__6 ) + mdl_MZ__exp__2 * ( 9. + 18. * mdl_dgw - 6. * ( 2. + mdl_dg1 + 3. * mdl_dgw ) * mdl_sth__exp__2 + 8. * ( 1. + 4. * mdl_dg1 - 2. * mdl_dgw ) * mdl_sth__exp__4 + 16. * ( -mdl_dg1 + mdl_dgw ) * mdl_sth__exp__6 ) * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + mdl_MB__exp__2 * ( -9. - 18. * mdl_dgw - 6. * ( 4. + 11. * mdl_dg1 - 3. * mdl_dgw ) * mdl_sth__exp__2 + 16. * ( 1. + 4. * mdl_dg1 - 2. * mdl_dgw ) * mdl_sth__exp__4 + 32. * ( -mdl_dg1 + mdl_dgw ) * mdl_sth__exp__6 ) * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) ) + 2. * mdl_ee * mdl_vevhat__exp__2 * ( 36. * mdl_cHj3 * mdl_MZ__exp__3 + 18. * mdl_cHl3 * mdl_MZ__exp__3 + 9. * ( 3. * mdl_cHbq - mdl_cHQ1 - mdl_cHQ3 ) * mdl_MB__exp__2 * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + 9. * mdl_cHQ1 * mdl_MZ__exp__2 * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + 9. * mdl_cHQ3 * mdl_MZ__exp__2 * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + 3. * mdl_cHWB * mdl_cth * ( -7. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) * mdl_sth * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + 8. * mdl_cHWB * mdl_cth * mdl_sth__exp__3 * ( 2. * mdl_MB__exp__2 * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + mdl_MZ__exp__2 * ( 19. * mdl_MZ + constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) ) ) - 8. * mdl_cHWB * mdl_cth * mdl_sth__exp__5 * ( 2. * mdl_MB__exp__2 * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + mdl_MZ__exp__2 * ( 19. * mdl_MZ + constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) ) ) - 6. * mdl_sth__exp__2 * ( 2. * ( mdl_cHbq + mdl_cHQ1 + mdl_cHQ3 ) * mdl_MB__exp__2 * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + mdl_MZ__exp__2 * ( ( 2. * mdl_cHd + 3. * mdl_cHe - 2. * mdl_cHj1 + 3. * ( 2. * mdl_cHj3 + mdl_cHl1 + mdl_cHl3 ) - 4. * mdl_cHu ) * mdl_MZ + ( mdl_cHbq + mdl_cHQ1 + mdl_cHQ3 ) * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) ) ) ) ) / ( mdl_ee * mdl_LambdaSMEFT__exp__2 * ( 2. * mdl_MZ__exp__3 * ( 27. - 54. * mdl_sth__exp__2 + 76. * mdl_sth__exp__4 ) + mdl_MZ__exp__2 * ( 9. - 12. * mdl_sth__exp__2 + 8. * mdl_sth__exp__4 ) * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + mdl_MB__exp__2 * ( -9. - 24. * mdl_sth__exp__2 + 16. * mdl_sth__exp__4 ) * constexpr_sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) ) ) );
  constexpr double mdl_g1sh = ( mdl_ee * ( 1. + mdl_dg1 - ( mdl_cHB * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 ) ) / mdl_cth;
  constexpr double mdl_ee__exp__3 = ( ( mdl_ee ) * ( mdl_ee ) * ( mdl_ee ) );
  constexpr double mdl_vevhat__exp__3 = ( ( mdl_vevhat ) * ( mdl_vevhat ) * ( mdl_vevhat ) );

  // Model couplings independent of aS
  // (none)

  // Model parameters dependent on aS
  //constexpr double mdl_sqrt__aS = //constexpr_sqrt( aS ); // now computed event-by-event (running alphas #373)
  //constexpr double G = 2. * mdl_sqrt__aS * //constexpr_sqrt( M_PI ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_gHgg2 = ( -7. * aS ) / ( 720. * M_PI ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_gHgg4 = aS / ( 360. * M_PI ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_gHgg5 = aS / ( 20. * M_PI ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_G__exp__2 = ( ( G ) * ( G ) ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_gHgg1 = mdl_G__exp__2 / ( 48. * ( ( M_PI ) * ( M_PI ) ) ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_gHgg3 = ( aS * G ) / ( 60. * M_PI ); // now computed event-by-event (running alphas #373)
  //constexpr cxsmpl<double> mdl_G__exp__3 = ( ( G ) * ( G ) * ( G ) ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_dWH = mdl_WH * ( -0.24161 * mdl_dGf + 0.96644 * mdl_dgw + 0.4832199999999999 * mdl_dkH - 0.11186509426655467 * mdl_dWW + ( 0.36410378449238195 * mdl_cHj3 * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 + ( 0.17608307708657747 * mdl_cHl3 * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 + ( 0.1636 * mdl_cHG * mdl_MT__exp__2 * mdl_vevhat__exp__2 ) / ( mdl_LambdaSMEFT__exp__2 * ( -0.5 * mdl_gHgg2 * mdl_MH__exp__2 + mdl_gHgg1 * mdl_MT__exp__2 ) ) + ( mdl_cHW * ( -0.35937785117066967 * mdl_gHaa * mdl_gHza + 0.006164 * mdl_cth * mdl_gHaa * mdl_sth + 0.00454 * mdl_gHza * mdl_sth__exp__2 ) * mdl_vevhat__exp__2 ) / ( mdl_gHaa * mdl_gHza * mdl_LambdaSMEFT__exp__2 ) + ( mdl_cHWB * ( -0.00454 * mdl_cth * mdl_gHza * mdl_sth + mdl_gHaa * ( -0.0030819999999999997 + 0.006163999999999999 * mdl_sth__exp__2 ) ) * mdl_vevhat__exp__2 ) / ( mdl_gHaa * mdl_gHza * mdl_LambdaSMEFT__exp__2 ) + ( mdl_cHB * ( -0.006163999999999999 * mdl_cth * mdl_gHaa * mdl_sth - 0.00454 * mdl_gHza * ( -1. + mdl_sth__exp__2 ) ) * mdl_vevhat__exp__2 ) / ( mdl_gHaa * mdl_gHza * mdl_LambdaSMEFT__exp__2 ) + mdl_dWHc + mdl_dWHb + mdl_dWHta ); // now computed event-by-event (running alphas #373)

  // Model couplings dependent on aS
  //constexpr cxsmpl<double> GC_6 = -( mdl_complexi * G ); // now computed event-by-event (running alphas #373)
  //constexpr cxsmpl<double> GC_7 = G; // now computed event-by-event (running alphas #373)
  //constexpr cxsmpl<double> GC_8 = mdl_complexi * mdl_G__exp__2; // now computed event-by-event (running alphas #373)

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

namespace Parameters_SMEFTsim_topU3l_MwScheme_UFO_dependentCouplings
{
  constexpr size_t ndcoup = 3; // #couplings that vary event by event because they depend on the running alphas QCD
  constexpr size_t idcoup_GC_6 = 0;
  constexpr size_t idcoup_GC_7 = 1;
  constexpr size_t idcoup_GC_8 = 2;
  struct DependentCouplings_sv
  {
    cxtype_sv GC_6;
    cxtype_sv GC_7;
    cxtype_sv GC_8;
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
    using namespace Parameters_SMEFTsim_topU3l_MwScheme_UFO;
#endif
    // NB: hardcode cxtype cI(0,1) instead of cxtype (or hardcoded cxsmpl) mdl_complexi (which exists in Parameters_SMEFTsim_topU3l_MwScheme_UFO) because:
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
      const fptype_sv mdl_gHgg2 = ( -7. * aS ) / ( 720. * M_PI );
      const fptype_sv mdl_gHgg4 = aS / ( 360. * M_PI );
      const fptype_sv mdl_gHgg5 = aS / ( 20. * M_PI );
      const fptype_sv mdl_G__exp__2 = ( ( G ) * ( G ) );
      const fptype_sv mdl_gHgg1 = mdl_G__exp__2 / ( 48. * ( ( M_PI ) * ( M_PI ) ) );
      const fptype_sv mdl_gHgg3 = ( aS * G ) / ( 60. * M_PI );
      constexpr cxsmpl<double> mdl_G__exp__3 = ( ( G ) * ( G ) * ( G ) );
      const fptype_sv mdl_dWH = mdl_WH * ( -0.24161 * mdl_dGf + 0.96644 * mdl_dgw + 0.4832199999999999 * mdl_dkH - 0.11186509426655467 * mdl_dWW + ( 0.36410378449238195 * mdl_cHj3 * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 + ( 0.17608307708657747 * mdl_cHl3 * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 + ( 0.1636 * mdl_cHG * mdl_MT__exp__2 * mdl_vevhat__exp__2 ) / ( mdl_LambdaSMEFT__exp__2 * ( -0.5 * mdl_gHgg2 * mdl_MH__exp__2 + mdl_gHgg1 * mdl_MT__exp__2 ) ) + ( mdl_cHW * ( -0.35937785117066967 * mdl_gHaa * mdl_gHza + 0.006164 * mdl_cth * mdl_gHaa * mdl_sth + 0.00454 * mdl_gHza * mdl_sth__exp__2 ) * mdl_vevhat__exp__2 ) / ( mdl_gHaa * mdl_gHza * mdl_LambdaSMEFT__exp__2 ) + ( mdl_cHWB * ( -0.00454 * mdl_cth * mdl_gHza * mdl_sth + mdl_gHaa * ( -0.0030819999999999997 + 0.006163999999999999 * mdl_sth__exp__2 ) ) * mdl_vevhat__exp__2 ) / ( mdl_gHaa * mdl_gHza * mdl_LambdaSMEFT__exp__2 ) + ( mdl_cHB * ( -0.006163999999999999 * mdl_cth * mdl_gHaa * mdl_sth - 0.00454 * mdl_gHza * ( -1. + mdl_sth__exp__2 ) ) * mdl_vevhat__exp__2 ) / ( mdl_gHaa * mdl_gHza * mdl_LambdaSMEFT__exp__2 ) + mdl_dWHc + mdl_dWHb + mdl_dWHta );
      // Model couplings dependent on aS
      out.GC_6 = -( cI * G );
      out.GC_7 = G;
      out.GC_8 = cI * mdl_G__exp__2;
    }
#else
    // ** NB #439: special handling is necessary ONLY FOR VECTORS OF FLOATS (variable Gs are vector floats, fixed parameters are scalar doubles)
    // Use an explicit loop to avoid <<error: conversion of scalar ‘double’ to vector ‘fptype_sv’ {aka ‘__vector(8) float’} involves truncation>>
    // Problems may come e.g. in EFTs from multiplying a vector float (related to aS-dependent G) by a scalar double (aS-independent parameters)
    fptype_v GC_6r_v;
    fptype_v GC_6i_v;
    fptype_v GC_7r_v;
    fptype_v GC_7i_v;
    fptype_v GC_8r_v;
    fptype_v GC_8i_v;
    for( int i = 0; i < neppV; i++ )
    {
      const fptype& G = G_sv[i];
      // Model parameters dependent on aS
      //const fptype mdl_sqrt__aS = constexpr_sqrt( aS );
      //const fptype G = 2. * mdl_sqrt__aS * constexpr_sqrt( M_PI );
      const fptype mdl_gHgg2 = ( -7. * aS ) / ( 720. * M_PI );
      const fptype mdl_gHgg4 = aS / ( 360. * M_PI );
      const fptype mdl_gHgg5 = aS / ( 20. * M_PI );
      const fptype mdl_G__exp__2 = ( ( G ) * ( G ) );
      const fptype mdl_gHgg1 = mdl_G__exp__2 / ( 48. * ( ( M_PI ) * ( M_PI ) ) );
      const fptype mdl_gHgg3 = ( aS * G ) / ( 60. * M_PI );
      constexpr cxsmpl<double> mdl_G__exp__3 = ( ( G ) * ( G ) * ( G ) );
      const fptype mdl_dWH = mdl_WH * ( -0.24161 * mdl_dGf + 0.96644 * mdl_dgw + 0.4832199999999999 * mdl_dkH - 0.11186509426655467 * mdl_dWW + ( 0.36410378449238195 * mdl_cHj3 * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 + ( 0.17608307708657747 * mdl_cHl3 * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 + ( 0.1636 * mdl_cHG * mdl_MT__exp__2 * mdl_vevhat__exp__2 ) / ( mdl_LambdaSMEFT__exp__2 * ( -0.5 * mdl_gHgg2 * mdl_MH__exp__2 + mdl_gHgg1 * mdl_MT__exp__2 ) ) + ( mdl_cHW * ( -0.35937785117066967 * mdl_gHaa * mdl_gHza + 0.006164 * mdl_cth * mdl_gHaa * mdl_sth + 0.00454 * mdl_gHza * mdl_sth__exp__2 ) * mdl_vevhat__exp__2 ) / ( mdl_gHaa * mdl_gHza * mdl_LambdaSMEFT__exp__2 ) + ( mdl_cHWB * ( -0.00454 * mdl_cth * mdl_gHza * mdl_sth + mdl_gHaa * ( -0.0030819999999999997 + 0.006163999999999999 * mdl_sth__exp__2 ) ) * mdl_vevhat__exp__2 ) / ( mdl_gHaa * mdl_gHza * mdl_LambdaSMEFT__exp__2 ) + ( mdl_cHB * ( -0.006163999999999999 * mdl_cth * mdl_gHaa * mdl_sth - 0.00454 * mdl_gHza * ( -1. + mdl_sth__exp__2 ) ) * mdl_vevhat__exp__2 ) / ( mdl_gHaa * mdl_gHza * mdl_LambdaSMEFT__exp__2 ) + mdl_dWHc + mdl_dWHb + mdl_dWHta );
      // Model couplings dependent on aS
      const cxtype GC_6 = -( cI * G );
      const cxtype GC_7 = G;
      const cxtype GC_8 = cI * mdl_G__exp__2;
      GC_6r_v[i] = cxreal( GC_6 );
      GC_6i_v[i] = cximag( GC_6 );
      GC_7r_v[i] = cxreal( GC_7 );
      GC_7i_v[i] = cximag( GC_7 );
      GC_8r_v[i] = cxreal( GC_8 );
      GC_8i_v[i] = cximag( GC_8 );
    }
    out.GC_6 = cxtype_v( GC_6r_v, GC_6i_v );
    out.GC_7 = cxtype_v( GC_7r_v, GC_7i_v );
    out.GC_8 = cxtype_v( GC_8r_v, GC_8i_v );
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

namespace Parameters_SMEFTsim_topU3l_MwScheme_UFO_independentCouplings
{
  constexpr size_t nicoup = 0; // #couplings that are fixed for all events because they do not depend on the running alphas QCD
  // NB: there are no aS-independent couplings in this physics process
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
    using namespace Parameters_SMEFTsim_topU3l_MwScheme_UFO_dependentCouplings;
    const fptype_sv& gs_sv = G_ACCESS::kernelAccessConst( gs );
    DependentCouplings_sv couplings_sv = computeDependentCouplings_fromG( gs_sv );
    fptype* GC_6s = C_ACCESS::idcoupAccessBuffer( couplings, idcoup_GC_6 );
    fptype* GC_7s = C_ACCESS::idcoupAccessBuffer( couplings, idcoup_GC_7 );
    fptype* GC_8s = C_ACCESS::idcoupAccessBuffer( couplings, idcoup_GC_8 );
    cxtype_sv_ref GC_6s_sv = C_ACCESS::kernelAccess( GC_6s );
    cxtype_sv_ref GC_7s_sv = C_ACCESS::kernelAccess( GC_7s );
    cxtype_sv_ref GC_8s_sv = C_ACCESS::kernelAccess( GC_8s );
    GC_6s_sv = couplings_sv.GC_6;
    GC_7s_sv = couplings_sv.GC_7;
    GC_8s_sv = couplings_sv.GC_8;
    mgDebug( 1, __FUNCTION__ );
    return;
  }
#pragma GCC diagnostic pop
}

//==========================================================================

#endif // Parameters_SMEFTsim_topU3l_MwScheme_UFO_H
