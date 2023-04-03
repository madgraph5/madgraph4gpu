//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.5.0_lo_vect, 2023-01-26
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_MSSM_SLHA2_H
#define Parameters_MSSM_SLHA2_H

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"
#include "mgOnGpuVectors.h"

//==========================================================================

#ifndef MGONGPU_HARDCODE_PARAM // this is only supported in SM processes (e.g. not in EFT models) for the moment (#439)
#error This non-SM physics process only supports MGONGPU_HARDCODE_PARAM builds (#439): please run "make HRDCOD=1"

#include "read_slha.h"

class Parameters_MSSM_SLHA2
{
public:

  static Parameters_MSSM_SLHA2* getInstance();

  // Define "zero"
  double zero, ZERO;

  // Model parameters independent of aS
  //double aS; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  double mdl_Wsl6, mdl_Wsl5, mdl_Wsl4, mdl_Wsu6, mdl_Wsd6, mdl_Wsu5, mdl_Wsd5, mdl_Wsu4, mdl_Wsd4, mdl_Wch2, mdl_Wneu4, mdl_Wneu3, mdl_Wch1, mdl_Wneu2, mdl_Wgo, mdl_Wsn3, mdl_Wsl3, mdl_Wsn2, mdl_Wsl2, mdl_Wsn1, mdl_Wsl1, mdl_Wsu3, mdl_Wsd3, mdl_Wsu2, mdl_Wsd2, mdl_Wsu1, mdl_Wsd1, mdl_WH, mdl_WA0, mdl_WH02, mdl_WH01, mdl_WW, mdl_WZ, mdl_WT, mdl_Ryu3x3, mdl_Rye3x3, mdl_Ryd3x3, mdl_RVV2x2, mdl_RVV2x1, mdl_RVV1x2, mdl_RVV1x1, mdl_RCKM3x3, mdl_RCKM2x2, mdl_RCKM1x1, mdl_RRu6x6, mdl_RRu6x3, mdl_RRu5x5, mdl_RRu4x4, mdl_RRu3x6, mdl_RRu3x3, mdl_RRu2x2, mdl_RRu1x1, mdl_RMNS3x3, mdl_RMNS2x2, mdl_RMNS1x1, mdl_RUU2x2, mdl_RUU2x1, mdl_RUU1x2, mdl_RUU1x1, mdl_Rtu3x3, mdl_Rte3x3, mdl_Rtd3x3, mdl_RRn3x3, mdl_RRn2x2, mdl_RRn1x1, aEWM1, mdl_RRl6x6, mdl_RRl6x3, mdl_RRl5x5, mdl_RRl4x4, mdl_RRl3x6, mdl_RRl3x3, mdl_RRl2x2, mdl_RRl1x1, mdl_RNN4x4, mdl_RNN4x3, mdl_RNN4x2, mdl_RNN4x1, mdl_RNN3x4, mdl_RNN3x3, mdl_RNN3x2, mdl_RNN3x1, mdl_RNN2x4, mdl_RNN2x3, mdl_RNN2x2, mdl_RNN2x1, mdl_RNN1x4, mdl_RNN1x3, mdl_RNN1x2, mdl_RNN1x1, mdl_RmU23x3, mdl_RmU21x1, mdl_RmQ23x3, mdl_RmQ21x1, mdl_mHu2, mdl_mHd2, mdl_RMx3, mdl_RMx2, mdl_RMx1, mdl_RmL23x3, mdl_RmL21x1, mdl_RmE23x3, mdl_RmE21x1, mdl_RmD23x3, mdl_RmD21x1, mdl_Msl6, mdl_Msl4, mdl_Msu6, mdl_Msd6, mdl_Msu4, mdl_Msd4, mdl_Mch2, mdl_Mneu4, mdl_Mneu3, mdl_Mch1, mdl_Mneu2, mdl_Mneu1, mdl_Mgo, mdl_Msn3, mdl_Msl3, mdl_Msn1, mdl_Msl1, mdl_Msu3, mdl_Msd3, mdl_Msu1, mdl_Msd1, mdl_MH, mdl_MA0, mdl_MH02, mdl_MH01, mdl_MW, mdl_MZ, mdl_Mta, mdl_MT, mdl_MB, mdl_MA2, mdl_tb, mdl_RMUH, mdl_alp, mdl_RRd6x6, mdl_RRd6x3, mdl_RRd5x5, mdl_RRd4x4, mdl_RRd3x6, mdl_RRd3x3, mdl_RRd2x2, mdl_RRd1x1, mdl_Msd5, mdl_Msd2, mdl_Msu5, mdl_Msu2, mdl_Msl5, mdl_Msl2, mdl_Msn2, mdl_RmU22x2, mdl_RmQ22x2, mdl_RmL22x2, mdl_RmE22x2, mdl_RmD22x2, mdl_conjg__Rn3x3, mdl_conjg__CKM3x3, mdl_Ru4x4, mdl_Ru1x1, mdl_Rn3x3, mdl_Rn1x1, mdl_Rl4x4, mdl_Rl1x1, mdl_Rd4x4, mdl_Rd1x1, mdl_I98x11, mdl_I97x11, mdl_I96x11, mdl_I93x11, mdl_I92x11, mdl_I87x11, mdl_I82x11, mdl_I74x11, mdl_I6x44, mdl_I5x11, mdl_I53x11, mdl_I52x44, mdl_I51x11, mdl_I39x11, mdl_I31x11, mdl_I26x44, mdl_I25x11, mdl_I12x11, mdl_I102x44, mdl_I101x44, mdl_I100x44, mdl_CKM3x3, mdl_atan__tb, mdl_beta, mdl_cw, mdl_MZ__exp__2, mdl_cw__exp__2, mdl_sw, mdl_cos__beta, mdl_sin__beta, mdl_sqrt__2, mdl_sw__exp__2, mdl_cos__alp, mdl_sin__alp, mdl_ee, mdl_gp, mdl_gw, mdl_vev, mdl_vd, mdl_vu, mdl_ee__exp__2;
  cxsmpl<double> mdl_mD21x1, mdl_mD22x2, mdl_mD23x3, mdl_mE21x1, mdl_mE22x2, mdl_mE23x3, mdl_mL21x1, mdl_mL22x2, mdl_mL23x3, mdl_mQ21x1, mdl_mQ22x2, mdl_mQ23x3, mdl_mU21x1, mdl_mU22x2, mdl_mU23x3, mdl_MUH, mdl_Mx1, mdl_Mx2, mdl_Mx3, mdl_NN1x1, mdl_NN1x2, mdl_NN1x3, mdl_NN1x4, mdl_NN2x1, mdl_NN2x2, mdl_NN2x3, mdl_NN2x4, mdl_NN3x1, mdl_NN3x2, mdl_NN3x3, mdl_NN3x4, mdl_NN4x1, mdl_NN4x2, mdl_NN4x3, mdl_NN4x4, mdl_Rd3x3, mdl_Rd3x6, mdl_Rd6x3, mdl_Rd6x6, mdl_Rl3x3, mdl_Rl3x6, mdl_Rl6x3, mdl_Rl6x6, mdl_Ru3x3, mdl_Ru3x6, mdl_Ru6x3, mdl_Ru6x6, mdl_UU1x1, mdl_UU1x2, mdl_UU2x1, mdl_UU2x2, mdl_VV1x1, mdl_VV1x2, mdl_VV2x1, mdl_VV2x2, mdl_td3x3, mdl_te3x3, mdl_tu3x3, mdl_yd3x3, mdl_ye3x3, mdl_yu3x3, mdl_bb, mdl_conjg__yu3x3, mdl_I1x33, mdl_conjg__yd3x3, mdl_I10x33, mdl_I10x36, mdl_conjg__Rd3x6, mdl_I100x33, mdl_I100x36, mdl_conjg__Rd6x6, mdl_I100x63, mdl_I100x66, mdl_conjg__Rl3x6, mdl_I101x33, mdl_I101x36, mdl_conjg__Rl6x6, mdl_I101x63, mdl_I101x66, mdl_conjg__Ru3x6, mdl_I102x33, mdl_I102x36, mdl_conjg__Ru6x6, mdl_I102x63, mdl_I102x66, mdl_I11x33, mdl_I11x36, mdl_conjg__Rd3x3, mdl_I12x33, mdl_I12x36, mdl_conjg__Rd6x3, mdl_I12x63, mdl_I12x66, mdl_I13x33, mdl_I13x36, mdl_I13x63, mdl_I13x66, mdl_conjg__td3x3, mdl_I14x33, mdl_I14x36, mdl_I14x63, mdl_I14x66, mdl_I15x33, mdl_I15x36, mdl_I15x63, mdl_I15x66, mdl_I16x33, mdl_I16x36, mdl_I16x63, mdl_I16x66, mdl_I17x33, mdl_I17x36, mdl_I17x63, mdl_I17x66, mdl_I18x33, mdl_I18x36, mdl_I18x63, mdl_I18x66, mdl_I19x33, mdl_I19x36, mdl_I19x63, mdl_I19x66, mdl_I2x33, mdl_I20x33, mdl_I21x33, mdl_conjg__ye3x3, mdl_I22x33, mdl_I23x33, mdl_I23x36, mdl_conjg__Rl3x3, mdl_I24x33, mdl_conjg__Rl6x3, mdl_I24x36, mdl_I25x33, mdl_I25x36, mdl_I25x63, mdl_I25x66, mdl_I26x33, mdl_I26x36, mdl_I26x63, mdl_I26x66, mdl_I27x33, mdl_I27x36, mdl_I28x33, mdl_I28x36, mdl_I29x33, mdl_I29x36, mdl_I3x33, mdl_I3x36, mdl_I30x33, mdl_I30x36, mdl_I31x33, mdl_I31x36, mdl_I31x63, mdl_I31x66, mdl_I32x33, mdl_I32x36, mdl_I32x63, mdl_I32x66, mdl_conjg__te3x3, mdl_I33x33, mdl_I33x36, mdl_I33x63, mdl_I33x66, mdl_I34x33, mdl_I34x36, mdl_I34x63, mdl_I34x66, mdl_I35x33, mdl_I35x36, mdl_I35x63, mdl_I35x66, mdl_I36x33, mdl_I36x36, mdl_I36x63, mdl_I36x66, mdl_I37x33, mdl_I37x36, mdl_I37x63, mdl_I37x66, mdl_I38x33, mdl_I38x36, mdl_I38x63, mdl_I38x66, mdl_I39x33, mdl_I39x36, mdl_I4x33, mdl_I4x36, mdl_I40x33, mdl_I40x36, mdl_I41x33, mdl_I41x36, mdl_I42x33, mdl_I42x36, mdl_I44x33, mdl_I45x33, mdl_I45x36, mdl_I46x33, mdl_I46x36, mdl_I47x33, mdl_I47x36, mdl_I48x33, mdl_I48x36, mdl_I49x33, mdl_I49x36, mdl_I5x33, mdl_I5x36, mdl_I5x63, mdl_I5x66, mdl_conjg__Ru3x3, mdl_I50x33, mdl_conjg__Ru6x3, mdl_I50x36, mdl_I51x33, mdl_I51x36, mdl_I51x63, mdl_I51x66, mdl_I52x33, mdl_I52x36, mdl_I52x63, mdl_I52x66, mdl_I53x33, mdl_I53x36, mdl_I53x63, mdl_I53x66, mdl_conjg__tu3x3, mdl_I54x33, mdl_I54x36, mdl_I54x63, mdl_I54x66, mdl_I55x33, mdl_I55x36, mdl_I55x63, mdl_I55x66, mdl_I56x33, mdl_I56x36, mdl_I56x63, mdl_I56x66, mdl_I57x33, mdl_I57x36, mdl_I57x63, mdl_I57x66, mdl_I58x33, mdl_I58x36, mdl_I58x63, mdl_I58x66, mdl_I59x33, mdl_I59x36, mdl_I59x63, mdl_I59x66, mdl_I6x33, mdl_I6x36, mdl_I6x63, mdl_I6x66, mdl_I60x33, mdl_I60x36, mdl_I60x63, mdl_I60x66, mdl_I61x33, mdl_I61x36, mdl_I62x33, mdl_I62x36, mdl_I63x33, mdl_I63x36, mdl_I64x33, mdl_I64x36, mdl_I65x33, mdl_I65x36, mdl_I66x33, mdl_I66x36, mdl_I66x63, mdl_I66x66, mdl_I67x33, mdl_I67x36, mdl_I67x63, mdl_I67x66, mdl_I68x33, mdl_I68x36, mdl_I68x63, mdl_I68x66, mdl_I69x33, mdl_I69x36, mdl_I69x63, mdl_I69x66, mdl_I7x33, mdl_I7x36, mdl_I70x33, mdl_I70x36, mdl_I70x63, mdl_I70x66, mdl_I71x33, mdl_I71x36, mdl_I71x63, mdl_I71x66, mdl_I72x33, mdl_I72x36, mdl_I72x63, mdl_I72x66, mdl_I73x33, mdl_I73x36, mdl_I73x63, mdl_I73x66, mdl_I74x33, mdl_I74x36, mdl_I74x63, mdl_I74x66, mdl_I75x33, mdl_I75x36, mdl_I75x63, mdl_I75x66, mdl_I76x33, mdl_I76x36, mdl_I76x63, mdl_I76x66, mdl_I77x33, mdl_I77x36, mdl_I77x63, mdl_I77x66, mdl_I78x33, mdl_I78x36, mdl_I78x63, mdl_I78x66, mdl_I79x33, mdl_I79x36, mdl_I79x63, mdl_I79x66, mdl_I8x33, mdl_I8x36, mdl_I80x33, mdl_I80x36, mdl_I80x63, mdl_I80x66, mdl_I81x33, mdl_I81x36, mdl_I81x63, mdl_I81x66, mdl_I82x33, mdl_I82x36, mdl_I83x33, mdl_I83x36, mdl_I84x33, mdl_I84x36, mdl_I85x33, mdl_I85x36, mdl_I86x33, mdl_I86x36, mdl_I88x33, mdl_I89x33, mdl_I89x36, mdl_I9x33, mdl_I9x36, mdl_I90x33, mdl_I90x36, mdl_I91x33, mdl_I91x36, mdl_I92x33, mdl_I92x36, mdl_I92x63, mdl_I92x66, mdl_I93x33, mdl_I93x36, mdl_I94x33, mdl_I94x36, mdl_I94x63, mdl_I94x66, mdl_I95x33, mdl_I95x36, mdl_I96x33, mdl_I96x36, mdl_I96x63, mdl_I96x66, mdl_I97x33, mdl_I97x36, mdl_I97x63, mdl_I97x66, mdl_I98x33, mdl_I98x36, mdl_I98x63, mdl_I98x66, mdl_I99x33, mdl_complexi, mdl_conjg__NN1x1, mdl_conjg__NN1x2, mdl_conjg__NN1x3, mdl_conjg__NN1x4, mdl_conjg__NN2x1, mdl_conjg__NN2x2, mdl_conjg__NN2x3, mdl_conjg__NN2x4, mdl_conjg__NN3x1, mdl_conjg__NN3x2, mdl_conjg__NN3x3, mdl_conjg__NN3x4, mdl_conjg__NN4x1, mdl_conjg__NN4x2, mdl_conjg__NN4x3, mdl_conjg__NN4x4, mdl_conjg__UU1x1, mdl_conjg__UU1x2, mdl_conjg__UU2x1, mdl_conjg__UU2x2, mdl_conjg__VV1x1, mdl_conjg__VV1x2, mdl_conjg__VV2x1, mdl_conjg__VV2x2, mdl_conjg__MUH;

  // Model couplings independent of aS
  // (none)

  // Model parameters dependent on aS
  //double mdl_sqrt__aS, G; // now computed event-by-event (running alphas #373)
  //cxsmpl<double> mdl_G__exp__2; // now computed event-by-event (running alphas #373)

  // Model couplings dependent on aS
  //cxsmpl<double> GC_6, GC_51; // now computed event-by-event (running alphas #373)

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

  static Parameters_MSSM_SLHA2* instance;
};

#else

#include <cassert>
#include <limits>

// Hardcoded constexpr physics parameters
namespace Parameters_MSSM_SLHA2 // keep the same name rather than HardcodedParameters_MSSM_SLHA2 for simplicity
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
  constexpr double mdl_Wsl6 = 2.699061e-01;
  constexpr double mdl_Wsl5 = 2.161216e-01;
  constexpr double mdl_Wsl4 = 2.161216e-01;
  constexpr double mdl_Wsu6 = 7.373133e+00;
  constexpr double mdl_Wsd6 = 8.015663e-01;
  constexpr double mdl_Wsu5 = 1.152973e+00;
  constexpr double mdl_Wsd5 = 2.858123e-01;
  constexpr double mdl_Wsu4 = 1.152973e+00;
  constexpr double mdl_Wsd4 = 2.858123e-01;
  constexpr double mdl_Wch2 = 2.486895e+00;
  constexpr double mdl_Wneu4 = 2.585851e+00;
  constexpr double mdl_Wneu3 = 1.915985e+00;
  constexpr double mdl_Wch1 = 1.704145e-02;
  constexpr double mdl_Wneu2 = 2.077700e-02;
  constexpr double mdl_Wgo = 5.506754e+00;
  constexpr double mdl_Wsn3 = 1.475190e-01;
  constexpr double mdl_Wsl3 = 1.483273e-01;
  constexpr double mdl_Wsn2 = 1.498816e-01;
  constexpr double mdl_Wsl2 = 2.136822e-01;
  constexpr double mdl_Wsn1 = 1.498816e-01;
  constexpr double mdl_Wsl1 = 2.136822e-01;
  constexpr double mdl_Wsu3 = 2.021596e+00;
  constexpr double mdl_Wsd3 = 3.736276e+00;
  constexpr double mdl_Wsu2 = 5.477195e+00;
  constexpr double mdl_Wsd2 = 5.312788e+00;
  constexpr double mdl_Wsu1 = 5.477195e+00;
  constexpr double mdl_Wsd1 = 5.312788e+00;
  constexpr double mdl_WH = 5.469628e-01;
  constexpr double mdl_WA0 = 6.321785e-01;
  constexpr double mdl_WH02 = 5.748014e-01;
  constexpr double mdl_WH01 = 1.986108e-03;
  constexpr double mdl_WW = 2.002822e+00;
  constexpr double mdl_WZ = 2.411433e+00;
  constexpr double mdl_WT = 1.561950e+00;
  constexpr double mdl_Ryu3x3 = 8.928445e-01;
  constexpr double mdl_Rye3x3 = 1.008908e-01;
  constexpr double mdl_Ryd3x3 = 1.388402e-01;
  constexpr double mdl_RVV2x2 = 9.725578e-01;
  constexpr double mdl_RVV2x1 = 2.326612e-01;
  constexpr double mdl_RVV1x2 = -2.326612e-01;
  constexpr double mdl_RVV1x1 = 9.725578e-01;
  constexpr double mdl_RCKM3x3 = 1.000000e+00;
  constexpr double mdl_RCKM2x2 = 1.000000e+00;
  constexpr double mdl_RCKM1x1 = 1.000000e+00;
  constexpr double mdl_RRu6x6 = -5.536450e-01;
  constexpr double mdl_RRu6x3 = 8.327528e-01;
  constexpr double mdl_RRu5x5 = 1.000000e+00;
  constexpr double mdl_RRu4x4 = 1.000000e+00;
  constexpr double mdl_RRu3x6 = 8.327528e-01;
  constexpr double mdl_RRu3x3 = 5.536450e-01;
  constexpr double mdl_RRu2x2 = 1.000000e+00;
  constexpr double mdl_RRu1x1 = 1.000000e+00;
  constexpr double mdl_RMNS3x3 = 1.000000e+00;
  constexpr double mdl_RMNS2x2 = 1.000000e+00;
  constexpr double mdl_RMNS1x1 = 1.000000e+00;
  constexpr double mdl_RUU2x2 = 9.168349e-01;
  constexpr double mdl_RUU2x1 = 3.992666e-01;
  constexpr double mdl_RUU1x2 = -3.992666e-01;
  constexpr double mdl_RUU1x1 = 9.168349e-01;
  constexpr double mdl_Rtu3x3 = -4.447525e+02;
  constexpr double mdl_Rte3x3 = -2.540197e+01;
  constexpr double mdl_Rtd3x3 = -1.106937e+02;
  constexpr double mdl_RRn3x3 = 1.000000e+00;
  constexpr double mdl_RRn2x2 = 1.000000e+00;
  constexpr double mdl_RRn1x1 = 1.000000e+00;
  //constexpr double aS = 1.180000e-01; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  constexpr double aEWM1 = 1.279340e+02;
  constexpr double mdl_RRl6x6 = -2.824872e-01;
  constexpr double mdl_RRl6x3 = 9.592711e-01;
  constexpr double mdl_RRl5x5 = 1.000000e+00;
  constexpr double mdl_RRl4x4 = 1.000000e+00;
  constexpr double mdl_RRl3x6 = 9.592711e-01;
  constexpr double mdl_RRl3x3 = 2.824872e-01;
  constexpr double mdl_RRl2x2 = 1.000000e+00;
  constexpr double mdl_RRl1x1 = 1.000000e+00;
  constexpr double mdl_RNN4x4 = -6.843778e-01;
  constexpr double mdl_RNN4x3 = 6.492260e-01;
  constexpr double mdl_RNN4x2 = 3.107390e-01;
  constexpr double mdl_RNN4x1 = -1.165071e-01;
  constexpr double mdl_RNN3x4 = 7.102270e-01;
  constexpr double mdl_RNN3x3 = 6.958775e-01;
  constexpr double mdl_RNN3x2 = 8.770049e-02;
  constexpr double mdl_RNN3x1 = -6.033880e-02;
  constexpr double mdl_RNN2x4 = 1.561507e-01;
  constexpr double mdl_RNN2x3 = -2.698467e-01;
  constexpr double mdl_RNN2x2 = 9.449493e-01;
  constexpr double mdl_RNN2x1 = 9.935054e-02;
  constexpr double mdl_RNN1x4 = -5.311861e-02;
  constexpr double mdl_RNN1x3 = 1.464340e-01;
  constexpr double mdl_RNN1x2 = -5.311036e-02;
  constexpr double mdl_RNN1x1 = 9.863644e-01;
  constexpr double mdl_RmU23x3 = 1.791371e+05;
  constexpr double mdl_RmU21x1 = 2.803821e+05;
  constexpr double mdl_RmQ23x3 = 2.487654e+05;
  constexpr double mdl_RmQ21x1 = 2.998367e+05;
  constexpr double mdl_mHu2 = -1.288001e+05;
  constexpr double mdl_mHd2 = 3.233749e+04;
  constexpr double mdl_RMx3 = 5.882630e+02;
  constexpr double mdl_RMx2 = 1.915042e+02;
  constexpr double mdl_RMx1 = 1.013965e+02;
  constexpr double mdl_RmL23x3 = 3.782868e+04;
  constexpr double mdl_RmL21x1 = 3.815567e+04;
  constexpr double mdl_RmE23x3 = 1.796764e+04;
  constexpr double mdl_RmE21x1 = 1.863063e+04;
  constexpr double mdl_RmD23x3 = 2.702620e+05;
  constexpr double mdl_RmD21x1 = 2.736847e+05;
  constexpr double mdl_Msl6 = 2.068678e+02;
  constexpr double mdl_Msl4 = 1.441028e+02;
  constexpr double mdl_Msu6 = 5.857858e+02;
  constexpr double mdl_Msd6 = 5.437267e+02;
  constexpr double mdl_Msu4 = 5.492593e+02;
  constexpr double mdl_Msd4 = 5.452285e+02;
  constexpr double mdl_Mch2 = 3.799393e+02;
  constexpr double mdl_Mneu4 = 3.817294e+02;
  constexpr double mdl_Mneu3 = -3.637560e+02;
  constexpr double mdl_Mch1 = 1.816965e+02;
  constexpr double mdl_Mneu2 = 1.810882e+02;
  constexpr double mdl_Mneu1 = 9.668807e+01;
  constexpr double mdl_Mgo = 6.077137e+02;
  constexpr double mdl_Msn3 = 1.847085e+02;
  constexpr double mdl_Msl3 = 1.344909e+02;
  constexpr double mdl_Msn1 = 1.852583e+02;
  constexpr double mdl_Msl1 = 2.029157e+02;
  constexpr double mdl_Msu3 = 3.996685e+02;
  constexpr double mdl_Msd3 = 5.130652e+02;
  constexpr double mdl_Msu1 = 5.611190e+02;
  constexpr double mdl_Msd1 = 5.684411e+02;
  constexpr double mdl_MH = 4.078790e+02;
  constexpr double mdl_MA0 = 3.995839e+02;
  constexpr double mdl_MH02 = 3.999601e+02;
  constexpr double mdl_MH01 = 1.108991e+02;
  constexpr double mdl_MW = 7.982901e+01;
  constexpr double mdl_MZ = 9.118760e+01;
  constexpr double mdl_Mta = 1.777000e+00;
  constexpr double mdl_MT = 1.750000e+02;
  constexpr double mdl_MB = 4.889917e+00;
  constexpr double mdl_MA2 = 1.664391e+05;
  constexpr double mdl_tb = 9.748624e+00;
  constexpr double mdl_RMUH = 3.576810e+02;
  constexpr double mdl_alp = -1.138252e-01;
  constexpr double mdl_RRd6x6 = 9.387379e-01;
  constexpr double mdl_RRd6x3 = -3.446319e-01;
  constexpr double mdl_RRd5x5 = 1.000000e+00;
  constexpr double mdl_RRd4x4 = 1.000000e+00;
  constexpr double mdl_RRd3x6 = 3.446319e-01;
  constexpr double mdl_RRd3x3 = 9.387379e-01;
  constexpr double mdl_RRd2x2 = 1.000000e+00;
  constexpr double mdl_RRd1x1 = 1.000000e+00;
  constexpr double mdl_Msd5 = 1. * mdl_Msd4;
  constexpr double mdl_Msd2 = 1. * mdl_Msd1;
  constexpr double mdl_Msu5 = 1. * mdl_Msu4;
  constexpr double mdl_Msu2 = 1. * mdl_Msu1;
  constexpr double mdl_Msl5 = 1. * mdl_Msl4;
  constexpr double mdl_Msl2 = 1. * mdl_Msl1;
  constexpr double mdl_Msn2 = 1. * mdl_Msn1;
  constexpr double mdl_RmU22x2 = 1. * mdl_RmU21x1;
  constexpr double mdl_RmQ22x2 = 1. * mdl_RmQ21x1;
  constexpr double mdl_RmL22x2 = 1. * mdl_RmL21x1;
  constexpr double mdl_RmE22x2 = 1. * mdl_RmE21x1;
  constexpr double mdl_RmD22x2 = 1. * mdl_RmD21x1;
  constexpr double mdl_conjg__Rn3x3 = 1.;
  constexpr double mdl_conjg__CKM3x3 = 1.;
  constexpr double mdl_Ru4x4 = 1.;
  constexpr double mdl_Ru1x1 = 1.;
  constexpr double mdl_Rn3x3 = 1.;
  constexpr double mdl_Rn1x1 = 1.;
  constexpr double mdl_Rl4x4 = 1.;
  constexpr double mdl_Rl1x1 = 1.;
  constexpr double mdl_Rd4x4 = 1.;
  constexpr double mdl_Rd1x1 = 1.;
  constexpr double mdl_I98x11 = 1.;
  constexpr double mdl_I97x11 = 1.;
  constexpr double mdl_I96x11 = 1.;
  constexpr double mdl_I93x11 = 1.;
  constexpr double mdl_I92x11 = 1.;
  constexpr double mdl_I87x11 = 1.;
  constexpr double mdl_I82x11 = 1.;
  constexpr double mdl_I74x11 = 1.;
  constexpr double mdl_I6x44 = 1.;
  constexpr double mdl_I5x11 = 1.;
  constexpr double mdl_I53x11 = 1.;
  constexpr double mdl_I52x44 = 1.;
  constexpr double mdl_I51x11 = 1.;
  constexpr double mdl_I39x11 = 1.;
  constexpr double mdl_I31x11 = 1.;
  constexpr double mdl_I26x44 = 1.;
  constexpr double mdl_I25x11 = 1.;
  constexpr double mdl_I12x11 = 1.;
  constexpr double mdl_I102x44 = 1.;
  constexpr double mdl_I101x44 = 1.;
  constexpr double mdl_I100x44 = 1.;
  constexpr double mdl_CKM3x3 = 1.;
  constexpr double mdl_atan__tb = atan( mdl_tb );
  constexpr double mdl_beta = mdl_atan__tb;
  constexpr double mdl_cw = mdl_MW / mdl_MZ;
  constexpr cxsmpl<double> mdl_mD21x1 = mdl_RmD21x1;
  constexpr cxsmpl<double> mdl_mD22x2 = mdl_RmD22x2;
  constexpr cxsmpl<double> mdl_mD23x3 = mdl_RmD23x3;
  constexpr cxsmpl<double> mdl_mE21x1 = mdl_RmE21x1;
  constexpr cxsmpl<double> mdl_mE22x2 = mdl_RmE22x2;
  constexpr cxsmpl<double> mdl_mE23x3 = mdl_RmE23x3;
  constexpr cxsmpl<double> mdl_mL21x1 = mdl_RmL21x1;
  constexpr cxsmpl<double> mdl_mL22x2 = mdl_RmL22x2;
  constexpr cxsmpl<double> mdl_mL23x3 = mdl_RmL23x3;
  constexpr cxsmpl<double> mdl_mQ21x1 = mdl_RmQ21x1;
  constexpr cxsmpl<double> mdl_mQ22x2 = mdl_RmQ22x2;
  constexpr cxsmpl<double> mdl_mQ23x3 = mdl_RmQ23x3;
  constexpr cxsmpl<double> mdl_mU21x1 = mdl_RmU21x1;
  constexpr cxsmpl<double> mdl_mU22x2 = mdl_RmU22x2;
  constexpr cxsmpl<double> mdl_mU23x3 = mdl_RmU23x3;
  constexpr cxsmpl<double> mdl_MUH = mdl_RMUH;
  constexpr cxsmpl<double> mdl_Mx1 = mdl_RMx1;
  constexpr cxsmpl<double> mdl_Mx2 = mdl_RMx2;
  constexpr cxsmpl<double> mdl_Mx3 = mdl_RMx3;
  constexpr cxsmpl<double> mdl_NN1x1 = mdl_RNN1x1;
  constexpr cxsmpl<double> mdl_NN1x2 = mdl_RNN1x2;
  constexpr cxsmpl<double> mdl_NN1x3 = mdl_RNN1x3;
  constexpr cxsmpl<double> mdl_NN1x4 = mdl_RNN1x4;
  constexpr cxsmpl<double> mdl_NN2x1 = mdl_RNN2x1;
  constexpr cxsmpl<double> mdl_NN2x2 = mdl_RNN2x2;
  constexpr cxsmpl<double> mdl_NN2x3 = mdl_RNN2x3;
  constexpr cxsmpl<double> mdl_NN2x4 = mdl_RNN2x4;
  constexpr cxsmpl<double> mdl_NN3x1 = mdl_RNN3x1;
  constexpr cxsmpl<double> mdl_NN3x2 = mdl_RNN3x2;
  constexpr cxsmpl<double> mdl_NN3x3 = mdl_RNN3x3;
  constexpr cxsmpl<double> mdl_NN3x4 = mdl_RNN3x4;
  constexpr cxsmpl<double> mdl_NN4x1 = mdl_RNN4x1;
  constexpr cxsmpl<double> mdl_NN4x2 = mdl_RNN4x2;
  constexpr cxsmpl<double> mdl_NN4x3 = mdl_RNN4x3;
  constexpr cxsmpl<double> mdl_NN4x4 = mdl_RNN4x4;
  constexpr cxsmpl<double> mdl_Rd3x3 = mdl_RRd3x3;
  constexpr cxsmpl<double> mdl_Rd3x6 = mdl_RRd3x6;
  constexpr cxsmpl<double> mdl_Rd6x3 = mdl_RRd6x3;
  constexpr cxsmpl<double> mdl_Rd6x6 = mdl_RRd6x6;
  constexpr cxsmpl<double> mdl_Rl3x3 = mdl_RRl3x3;
  constexpr cxsmpl<double> mdl_Rl3x6 = mdl_RRl3x6;
  constexpr cxsmpl<double> mdl_Rl6x3 = mdl_RRl6x3;
  constexpr cxsmpl<double> mdl_Rl6x6 = mdl_RRl6x6;
  constexpr cxsmpl<double> mdl_Ru3x3 = mdl_RRu3x3;
  constexpr cxsmpl<double> mdl_Ru3x6 = mdl_RRu3x6;
  constexpr cxsmpl<double> mdl_Ru6x3 = mdl_RRu6x3;
  constexpr cxsmpl<double> mdl_Ru6x6 = mdl_RRu6x6;
  constexpr cxsmpl<double> mdl_UU1x1 = mdl_RUU1x1;
  constexpr cxsmpl<double> mdl_UU1x2 = mdl_RUU1x2;
  constexpr cxsmpl<double> mdl_UU2x1 = mdl_RUU2x1;
  constexpr cxsmpl<double> mdl_UU2x2 = mdl_RUU2x2;
  constexpr cxsmpl<double> mdl_VV1x1 = mdl_RVV1x1;
  constexpr cxsmpl<double> mdl_VV1x2 = mdl_RVV1x2;
  constexpr cxsmpl<double> mdl_VV2x1 = mdl_RVV2x1;
  constexpr cxsmpl<double> mdl_VV2x2 = mdl_RVV2x2;
  constexpr cxsmpl<double> mdl_td3x3 = mdl_Rtd3x3;
  constexpr cxsmpl<double> mdl_te3x3 = mdl_Rte3x3;
  constexpr cxsmpl<double> mdl_tu3x3 = mdl_Rtu3x3;
  constexpr cxsmpl<double> mdl_yd3x3 = mdl_Ryd3x3;
  constexpr cxsmpl<double> mdl_ye3x3 = mdl_Rye3x3;
  constexpr cxsmpl<double> mdl_yu3x3 = mdl_Ryu3x3;
  constexpr double mdl_MZ__exp__2 = ( ( mdl_MZ ) * ( mdl_MZ ) );
  constexpr cxsmpl<double> mdl_bb = ( ( -mdl_mHd2 + mdl_mHu2 - mdl_MZ__exp__2 * cos( 2. * mdl_beta ) ) * tan( 2. * mdl_beta ) ) / 2.;
  constexpr double mdl_cw__exp__2 = ( ( mdl_cw ) * ( mdl_cw ) );
  constexpr double mdl_sw = constexpr_sqrt( 1. - mdl_cw__exp__2 );
  constexpr double mdl_cos__beta = cos( mdl_beta );
  constexpr double mdl_sin__beta = sin( mdl_beta );
  constexpr cxsmpl<double> mdl_conjg__yu3x3 = conj( mdl_yu3x3 );
  constexpr cxsmpl<double> mdl_I1x33 = mdl_conjg__CKM3x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_conjg__yd3x3 = conj( mdl_yd3x3 );
  constexpr cxsmpl<double> mdl_I10x33 = mdl_Rd3x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I10x36 = mdl_Rd6x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_conjg__Rd3x6 = conj( mdl_Rd3x6 );
  constexpr cxsmpl<double> mdl_I100x33 = mdl_Rd3x6 * mdl_conjg__Rd3x6;
  constexpr cxsmpl<double> mdl_I100x36 = mdl_Rd6x6 * mdl_conjg__Rd3x6;
  constexpr cxsmpl<double> mdl_conjg__Rd6x6 = conj( mdl_Rd6x6 );
  constexpr cxsmpl<double> mdl_I100x63 = mdl_Rd3x6 * mdl_conjg__Rd6x6;
  constexpr cxsmpl<double> mdl_I100x66 = mdl_Rd6x6 * mdl_conjg__Rd6x6;
  constexpr cxsmpl<double> mdl_conjg__Rl3x6 = conj( mdl_Rl3x6 );
  constexpr cxsmpl<double> mdl_I101x33 = mdl_Rl3x6 * mdl_conjg__Rl3x6;
  constexpr cxsmpl<double> mdl_I101x36 = mdl_Rl6x6 * mdl_conjg__Rl3x6;
  constexpr cxsmpl<double> mdl_conjg__Rl6x6 = conj( mdl_Rl6x6 );
  constexpr cxsmpl<double> mdl_I101x63 = mdl_Rl3x6 * mdl_conjg__Rl6x6;
  constexpr cxsmpl<double> mdl_I101x66 = mdl_Rl6x6 * mdl_conjg__Rl6x6;
  constexpr cxsmpl<double> mdl_conjg__Ru3x6 = conj( mdl_Ru3x6 );
  constexpr cxsmpl<double> mdl_I102x33 = mdl_Ru3x6 * mdl_conjg__Ru3x6;
  constexpr cxsmpl<double> mdl_I102x36 = mdl_Ru6x6 * mdl_conjg__Ru3x6;
  constexpr cxsmpl<double> mdl_conjg__Ru6x6 = conj( mdl_Ru6x6 );
  constexpr cxsmpl<double> mdl_I102x63 = mdl_Ru3x6 * mdl_conjg__Ru6x6;
  constexpr cxsmpl<double> mdl_I102x66 = mdl_Ru6x6 * mdl_conjg__Ru6x6;
  constexpr cxsmpl<double> mdl_I11x33 = mdl_Rd3x6 * mdl_yd3x3;
  constexpr cxsmpl<double> mdl_I11x36 = mdl_Rd6x6 * mdl_yd3x3;
  constexpr cxsmpl<double> mdl_conjg__Rd3x3 = conj( mdl_Rd3x3 );
  constexpr cxsmpl<double> mdl_I12x33 = mdl_Rd3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I12x36 = mdl_Rd6x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_conjg__Rd6x3 = conj( mdl_Rd6x3 );
  constexpr cxsmpl<double> mdl_I12x63 = mdl_Rd3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I12x66 = mdl_Rd6x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I13x33 = mdl_Rd3x6 * mdl_conjg__Rd3x6;
  constexpr cxsmpl<double> mdl_I13x36 = mdl_Rd6x6 * mdl_conjg__Rd3x6;
  constexpr cxsmpl<double> mdl_I13x63 = mdl_Rd3x6 * mdl_conjg__Rd6x6;
  constexpr cxsmpl<double> mdl_I13x66 = mdl_Rd6x6 * mdl_conjg__Rd6x6;
  constexpr cxsmpl<double> mdl_conjg__td3x3 = conj( mdl_td3x3 );
  constexpr cxsmpl<double> mdl_I14x33 = mdl_Rd3x3 * mdl_conjg__Rd3x6 * mdl_conjg__td3x3;
  constexpr cxsmpl<double> mdl_I14x36 = mdl_Rd6x3 * mdl_conjg__Rd3x6 * mdl_conjg__td3x3;
  constexpr cxsmpl<double> mdl_I14x63 = mdl_Rd3x3 * mdl_conjg__Rd6x6 * mdl_conjg__td3x3;
  constexpr cxsmpl<double> mdl_I14x66 = mdl_Rd6x3 * mdl_conjg__Rd6x6 * mdl_conjg__td3x3;
  constexpr cxsmpl<double> mdl_I15x33 = mdl_Rd3x3 * mdl_conjg__Rd3x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I15x36 = mdl_Rd6x3 * mdl_conjg__Rd3x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I15x63 = mdl_Rd3x3 * mdl_conjg__Rd6x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I15x66 = mdl_Rd6x3 * mdl_conjg__Rd6x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I16x33 = mdl_Rd3x6 * mdl_td3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I16x36 = mdl_Rd6x6 * mdl_td3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I16x63 = mdl_Rd3x6 * mdl_td3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I16x66 = mdl_Rd6x6 * mdl_td3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I17x33 = mdl_Rd3x3 * mdl_yd3x3 * mdl_conjg__Rd3x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I17x36 = mdl_Rd6x3 * mdl_yd3x3 * mdl_conjg__Rd3x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I17x63 = mdl_Rd3x3 * mdl_yd3x3 * mdl_conjg__Rd6x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I17x66 = mdl_Rd6x3 * mdl_yd3x3 * mdl_conjg__Rd6x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I18x33 = mdl_Rd3x6 * mdl_yd3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I18x36 = mdl_Rd6x6 * mdl_yd3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I18x63 = mdl_Rd3x6 * mdl_yd3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I18x66 = mdl_Rd6x6 * mdl_yd3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I19x33 = mdl_Rd3x6 * mdl_yd3x3 * mdl_conjg__Rd3x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I19x36 = mdl_Rd6x6 * mdl_yd3x3 * mdl_conjg__Rd3x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I19x63 = mdl_Rd3x6 * mdl_yd3x3 * mdl_conjg__Rd6x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I19x66 = mdl_Rd6x6 * mdl_yd3x3 * mdl_conjg__Rd6x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I2x33 = mdl_yd3x3 * mdl_conjg__CKM3x3;
  constexpr cxsmpl<double> mdl_I20x33 = mdl_CKM3x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I21x33 = mdl_CKM3x3 * mdl_yu3x3;
  constexpr cxsmpl<double> mdl_conjg__ye3x3 = conj( mdl_ye3x3 );
  constexpr cxsmpl<double> mdl_I22x33 = mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I23x33 = mdl_conjg__Rl3x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I23x36 = mdl_conjg__Rl6x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_conjg__Rl3x3 = conj( mdl_Rl3x3 );
  constexpr cxsmpl<double> mdl_I24x33 = mdl_ye3x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_conjg__Rl6x3 = conj( mdl_Rl6x3 );
  constexpr cxsmpl<double> mdl_I24x36 = mdl_ye3x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I25x33 = mdl_Rl3x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I25x36 = mdl_Rl6x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I25x63 = mdl_Rl3x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I25x66 = mdl_Rl6x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I26x33 = mdl_Rl3x6 * mdl_conjg__Rl3x6;
  constexpr cxsmpl<double> mdl_I26x36 = mdl_Rl6x6 * mdl_conjg__Rl3x6;
  constexpr cxsmpl<double> mdl_I26x63 = mdl_Rl3x6 * mdl_conjg__Rl6x6;
  constexpr cxsmpl<double> mdl_I26x66 = mdl_Rl6x6 * mdl_conjg__Rl6x6;
  constexpr cxsmpl<double> mdl_I27x33 = mdl_Rl3x3 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I27x36 = mdl_Rl6x3 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I28x33 = mdl_Rl3x6 * mdl_ye3x3;
  constexpr cxsmpl<double> mdl_I28x36 = mdl_Rl6x6 * mdl_ye3x3;
  constexpr cxsmpl<double> mdl_I29x33 = mdl_Rl3x3;
  constexpr cxsmpl<double> mdl_I29x36 = mdl_Rl6x3;
  constexpr cxsmpl<double> mdl_I3x33 = mdl_conjg__Rd3x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I3x36 = mdl_conjg__Rd6x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I30x33 = mdl_Rl3x6 * mdl_ye3x3;
  constexpr cxsmpl<double> mdl_I30x36 = mdl_Rl6x6 * mdl_ye3x3;
  constexpr cxsmpl<double> mdl_I31x33 = mdl_Rl3x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I31x36 = mdl_Rl6x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I31x63 = mdl_Rl3x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I31x66 = mdl_Rl6x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I32x33 = mdl_Rl3x6 * mdl_conjg__Rl3x6;
  constexpr cxsmpl<double> mdl_I32x36 = mdl_Rl6x6 * mdl_conjg__Rl3x6;
  constexpr cxsmpl<double> mdl_I32x63 = mdl_Rl3x6 * mdl_conjg__Rl6x6;
  constexpr cxsmpl<double> mdl_I32x66 = mdl_Rl6x6 * mdl_conjg__Rl6x6;
  constexpr cxsmpl<double> mdl_conjg__te3x3 = conj( mdl_te3x3 );
  constexpr cxsmpl<double> mdl_I33x33 = mdl_Rl3x3 * mdl_conjg__Rl3x6 * mdl_conjg__te3x3;
  constexpr cxsmpl<double> mdl_I33x36 = mdl_Rl6x3 * mdl_conjg__Rl3x6 * mdl_conjg__te3x3;
  constexpr cxsmpl<double> mdl_I33x63 = mdl_Rl3x3 * mdl_conjg__Rl6x6 * mdl_conjg__te3x3;
  constexpr cxsmpl<double> mdl_I33x66 = mdl_Rl6x3 * mdl_conjg__Rl6x6 * mdl_conjg__te3x3;
  constexpr cxsmpl<double> mdl_I34x33 = mdl_Rl3x3 * mdl_conjg__Rl3x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I34x36 = mdl_Rl6x3 * mdl_conjg__Rl3x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I34x63 = mdl_Rl3x3 * mdl_conjg__Rl6x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I34x66 = mdl_Rl6x3 * mdl_conjg__Rl6x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I35x33 = mdl_Rl3x6 * mdl_te3x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I35x36 = mdl_Rl6x6 * mdl_te3x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I35x63 = mdl_Rl3x6 * mdl_te3x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I35x66 = mdl_Rl6x6 * mdl_te3x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I36x33 = mdl_Rl3x3 * mdl_ye3x3 * mdl_conjg__Rl3x3 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I36x36 = mdl_Rl6x3 * mdl_ye3x3 * mdl_conjg__Rl3x3 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I36x63 = mdl_Rl3x3 * mdl_ye3x3 * mdl_conjg__Rl6x3 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I36x66 = mdl_Rl6x3 * mdl_ye3x3 * mdl_conjg__Rl6x3 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I37x33 = mdl_Rl3x6 * mdl_ye3x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I37x36 = mdl_Rl6x6 * mdl_ye3x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I37x63 = mdl_Rl3x6 * mdl_ye3x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I37x66 = mdl_Rl6x6 * mdl_ye3x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I38x33 = mdl_Rl3x6 * mdl_ye3x3 * mdl_conjg__Rl3x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I38x36 = mdl_Rl6x6 * mdl_ye3x3 * mdl_conjg__Rl3x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I38x63 = mdl_Rl3x6 * mdl_ye3x3 * mdl_conjg__Rl6x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I38x66 = mdl_Rl6x6 * mdl_ye3x3 * mdl_conjg__Rl6x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I39x33 = mdl_Rl3x3 * mdl_conjg__Rn3x3;
  constexpr cxsmpl<double> mdl_I39x36 = mdl_Rl6x3 * mdl_conjg__Rn3x3;
  constexpr cxsmpl<double> mdl_I4x33 = mdl_yd3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I4x36 = mdl_yd3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I40x33 = mdl_Rl3x6 * mdl_te3x3 * mdl_conjg__Rn3x3;
  constexpr cxsmpl<double> mdl_I40x36 = mdl_Rl6x6 * mdl_te3x3 * mdl_conjg__Rn3x3;
  constexpr cxsmpl<double> mdl_I41x33 = mdl_Rl3x3 * mdl_ye3x3 * mdl_conjg__Rn3x3 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I41x36 = mdl_Rl6x3 * mdl_ye3x3 * mdl_conjg__Rn3x3 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I42x33 = mdl_Rl3x6 * mdl_ye3x3 * mdl_conjg__Rn3x3;
  constexpr cxsmpl<double> mdl_I42x36 = mdl_Rl6x6 * mdl_ye3x3 * mdl_conjg__Rn3x3;
  constexpr cxsmpl<double> mdl_I44x33 = mdl_Rn3x3 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I45x33 = mdl_Rn3x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I45x36 = mdl_Rn3x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I46x33 = mdl_Rn3x3 * mdl_conjg__Rl3x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I46x36 = mdl_Rn3x3 * mdl_conjg__Rl6x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I47x33 = mdl_Rn3x3 * mdl_conjg__Rl3x6 * mdl_conjg__te3x3;
  constexpr cxsmpl<double> mdl_I47x36 = mdl_Rn3x3 * mdl_conjg__Rl6x6 * mdl_conjg__te3x3;
  constexpr cxsmpl<double> mdl_I48x33 = mdl_Rn3x3 * mdl_ye3x3 * mdl_conjg__Rl3x3 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I48x36 = mdl_Rn3x3 * mdl_ye3x3 * mdl_conjg__Rl6x3 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I49x33 = mdl_conjg__Ru3x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I49x36 = mdl_conjg__Ru6x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I5x33 = mdl_Rd3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I5x36 = mdl_Rd6x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I5x63 = mdl_Rd3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I5x66 = mdl_Rd6x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_conjg__Ru3x3 = conj( mdl_Ru3x3 );
  constexpr cxsmpl<double> mdl_I50x33 = mdl_yu3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_conjg__Ru6x3 = conj( mdl_Ru6x3 );
  constexpr cxsmpl<double> mdl_I50x36 = mdl_yu3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I51x33 = mdl_Ru3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I51x36 = mdl_Ru6x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I51x63 = mdl_Ru3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I51x66 = mdl_Ru6x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I52x33 = mdl_Ru3x6 * mdl_conjg__Ru3x6;
  constexpr cxsmpl<double> mdl_I52x36 = mdl_Ru6x6 * mdl_conjg__Ru3x6;
  constexpr cxsmpl<double> mdl_I52x63 = mdl_Ru3x6 * mdl_conjg__Ru6x6;
  constexpr cxsmpl<double> mdl_I52x66 = mdl_Ru6x6 * mdl_conjg__Ru6x6;
  constexpr cxsmpl<double> mdl_I53x33 = mdl_Rd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I53x36 = mdl_Rd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I53x63 = mdl_Rd6x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I53x66 = mdl_Rd6x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_conjg__tu3x3 = conj( mdl_tu3x3 );
  constexpr cxsmpl<double> mdl_I54x33 = mdl_Rd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x6 * mdl_conjg__tu3x3;
  constexpr cxsmpl<double> mdl_I54x36 = mdl_Rd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x6 * mdl_conjg__tu3x3;
  constexpr cxsmpl<double> mdl_I54x63 = mdl_Rd6x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x6 * mdl_conjg__tu3x3;
  constexpr cxsmpl<double> mdl_I54x66 = mdl_Rd6x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x6 * mdl_conjg__tu3x3;
  constexpr cxsmpl<double> mdl_I55x33 = mdl_Rd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I55x36 = mdl_Rd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I55x63 = mdl_Rd6x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I55x66 = mdl_Rd6x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I56x33 = mdl_Rd3x6 * mdl_td3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I56x36 = mdl_Rd3x6 * mdl_td3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I56x63 = mdl_Rd6x6 * mdl_td3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I56x66 = mdl_Rd6x6 * mdl_td3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I57x33 = mdl_Rd3x6 * mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I57x36 = mdl_Rd3x6 * mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I57x63 = mdl_Rd6x6 * mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I57x66 = mdl_Rd6x6 * mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I58x33 = mdl_Rd3x3 * mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I58x36 = mdl_Rd3x3 * mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I58x63 = mdl_Rd6x3 * mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I58x66 = mdl_Rd6x3 * mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I59x33 = mdl_Rd3x6 * mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I59x36 = mdl_Rd3x6 * mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I59x63 = mdl_Rd6x6 * mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I59x66 = mdl_Rd6x6 * mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I6x33 = mdl_Rd3x6 * mdl_conjg__Rd3x6;
  constexpr cxsmpl<double> mdl_I6x36 = mdl_Rd6x6 * mdl_conjg__Rd3x6;
  constexpr cxsmpl<double> mdl_I6x63 = mdl_Rd3x6 * mdl_conjg__Rd6x6;
  constexpr cxsmpl<double> mdl_I6x66 = mdl_Rd6x6 * mdl_conjg__Rd6x6;
  constexpr cxsmpl<double> mdl_I60x33 = mdl_Rd3x3 * mdl_yu3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I60x36 = mdl_Rd3x3 * mdl_yu3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I60x63 = mdl_Rd6x3 * mdl_yu3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I60x66 = mdl_Rd6x3 * mdl_yu3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I61x33 = mdl_Ru3x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I61x36 = mdl_Ru6x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I62x33 = mdl_Ru3x6 * mdl_yu3x3;
  constexpr cxsmpl<double> mdl_I62x36 = mdl_Ru6x6 * mdl_yu3x3;
  constexpr cxsmpl<double> mdl_I63x33 = mdl_CKM3x3 * mdl_Ru3x3;
  constexpr cxsmpl<double> mdl_I63x36 = mdl_CKM3x3 * mdl_Ru6x3;
  constexpr cxsmpl<double> mdl_I64x33 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I64x36 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I65x33 = mdl_CKM3x3 * mdl_Ru3x6 * mdl_yu3x3;
  constexpr cxsmpl<double> mdl_I65x36 = mdl_CKM3x3 * mdl_Ru6x6 * mdl_yu3x3;
  constexpr cxsmpl<double> mdl_I66x33 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I66x36 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I66x63 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I66x66 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I67x33 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_conjg__Rd3x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I67x36 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_conjg__Rd3x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I67x63 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_conjg__Rd6x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I67x66 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_conjg__Rd6x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I68x33 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_conjg__Rd3x6 * mdl_conjg__td3x3;
  constexpr cxsmpl<double> mdl_I68x36 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_conjg__Rd3x6 * mdl_conjg__td3x3;
  constexpr cxsmpl<double> mdl_I68x63 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_conjg__Rd6x6 * mdl_conjg__td3x3;
  constexpr cxsmpl<double> mdl_I68x66 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_conjg__Rd6x6 * mdl_conjg__td3x3;
  constexpr cxsmpl<double> mdl_I69x33 = mdl_CKM3x3 * mdl_Ru3x6 * mdl_tu3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I69x36 = mdl_CKM3x3 * mdl_Ru6x6 * mdl_tu3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I69x63 = mdl_CKM3x3 * mdl_Ru3x6 * mdl_tu3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I69x66 = mdl_CKM3x3 * mdl_Ru6x6 * mdl_tu3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I7x33 = mdl_Rd3x3 * mdl_conjg__CKM3x3;
  constexpr cxsmpl<double> mdl_I7x36 = mdl_Rd6x3 * mdl_conjg__CKM3x3;
  constexpr cxsmpl<double> mdl_I70x33 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_yd3x3 * mdl_conjg__Rd3x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I70x36 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_yd3x3 * mdl_conjg__Rd3x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I70x63 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_yd3x3 * mdl_conjg__Rd6x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I70x66 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_yd3x3 * mdl_conjg__Rd6x3 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I71x33 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_yu3x3 * mdl_conjg__Rd3x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I71x36 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_yu3x3 * mdl_conjg__Rd3x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I71x63 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_yu3x3 * mdl_conjg__Rd6x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I71x66 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_yu3x3 * mdl_conjg__Rd6x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I72x33 = mdl_CKM3x3 * mdl_Ru3x6 * mdl_yu3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I72x36 = mdl_CKM3x3 * mdl_Ru6x6 * mdl_yu3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I72x63 = mdl_CKM3x3 * mdl_Ru3x6 * mdl_yu3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I72x66 = mdl_CKM3x3 * mdl_Ru6x6 * mdl_yu3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I73x33 = mdl_CKM3x3 * mdl_Ru3x6 * mdl_yu3x3 * mdl_conjg__Rd3x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I73x36 = mdl_CKM3x3 * mdl_Ru6x6 * mdl_yu3x3 * mdl_conjg__Rd3x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I73x63 = mdl_CKM3x3 * mdl_Ru3x6 * mdl_yu3x3 * mdl_conjg__Rd6x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I73x66 = mdl_CKM3x3 * mdl_Ru6x6 * mdl_yu3x3 * mdl_conjg__Rd6x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I74x33 = mdl_Ru3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I74x36 = mdl_Ru6x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I74x63 = mdl_Ru3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I74x66 = mdl_Ru6x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I75x33 = mdl_Ru3x6 * mdl_conjg__Ru3x6;
  constexpr cxsmpl<double> mdl_I75x36 = mdl_Ru6x6 * mdl_conjg__Ru3x6;
  constexpr cxsmpl<double> mdl_I75x63 = mdl_Ru3x6 * mdl_conjg__Ru6x6;
  constexpr cxsmpl<double> mdl_I75x66 = mdl_Ru6x6 * mdl_conjg__Ru6x6;
  constexpr cxsmpl<double> mdl_I76x33 = mdl_Ru3x3 * mdl_conjg__Ru3x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I76x36 = mdl_Ru6x3 * mdl_conjg__Ru3x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I76x63 = mdl_Ru3x3 * mdl_conjg__Ru6x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I76x66 = mdl_Ru6x3 * mdl_conjg__Ru6x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I77x33 = mdl_Ru3x3 * mdl_conjg__Ru3x6 * mdl_conjg__tu3x3;
  constexpr cxsmpl<double> mdl_I77x36 = mdl_Ru6x3 * mdl_conjg__Ru3x6 * mdl_conjg__tu3x3;
  constexpr cxsmpl<double> mdl_I77x63 = mdl_Ru3x3 * mdl_conjg__Ru6x6 * mdl_conjg__tu3x3;
  constexpr cxsmpl<double> mdl_I77x66 = mdl_Ru6x3 * mdl_conjg__Ru6x6 * mdl_conjg__tu3x3;
  constexpr cxsmpl<double> mdl_I78x33 = mdl_Ru3x6 * mdl_tu3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I78x36 = mdl_Ru6x6 * mdl_tu3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I78x63 = mdl_Ru3x6 * mdl_tu3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I78x66 = mdl_Ru6x6 * mdl_tu3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I79x33 = mdl_Ru3x6 * mdl_yu3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I79x36 = mdl_Ru6x6 * mdl_yu3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I79x63 = mdl_Ru3x6 * mdl_yu3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I79x66 = mdl_Ru6x6 * mdl_yu3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I8x33 = mdl_Rd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I8x36 = mdl_Rd6x3 * mdl_conjg__CKM3x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I80x33 = mdl_Ru3x3 * mdl_yu3x3 * mdl_conjg__Ru3x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I80x36 = mdl_Ru6x3 * mdl_yu3x3 * mdl_conjg__Ru3x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I80x63 = mdl_Ru3x3 * mdl_yu3x3 * mdl_conjg__Ru6x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I80x66 = mdl_Ru6x3 * mdl_yu3x3 * mdl_conjg__Ru6x3 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I81x33 = mdl_Ru3x6 * mdl_yu3x3 * mdl_conjg__Ru3x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I81x36 = mdl_Ru6x6 * mdl_yu3x3 * mdl_conjg__Ru3x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I81x63 = mdl_Ru3x6 * mdl_yu3x3 * mdl_conjg__Ru6x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I81x66 = mdl_Ru6x6 * mdl_yu3x3 * mdl_conjg__Ru6x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I82x33 = mdl_CKM3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I82x36 = mdl_CKM3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I83x33 = mdl_CKM3x3 * mdl_conjg__Rd3x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I83x36 = mdl_CKM3x3 * mdl_conjg__Rd6x6 * mdl_conjg__yd3x3;
  constexpr cxsmpl<double> mdl_I84x33 = mdl_CKM3x3 * mdl_yu3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I84x36 = mdl_CKM3x3 * mdl_yu3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I85x33 = mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I85x36 = mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I86x33 = mdl_conjg__Rl3x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I86x36 = mdl_conjg__Rl6x6 * mdl_conjg__ye3x3;
  constexpr cxsmpl<double> mdl_I88x33 = mdl_ye3x3 * mdl_conjg__Rn3x3;
  constexpr cxsmpl<double> mdl_I89x33 = mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I89x36 = mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I9x33 = mdl_Rd3x6 * mdl_yd3x3 * mdl_conjg__CKM3x3;
  constexpr cxsmpl<double> mdl_I9x36 = mdl_Rd6x6 * mdl_yd3x3 * mdl_conjg__CKM3x3;
  constexpr cxsmpl<double> mdl_I90x33 = mdl_conjg__CKM3x3 * mdl_conjg__Ru3x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I90x36 = mdl_conjg__CKM3x3 * mdl_conjg__Ru6x6 * mdl_conjg__yu3x3;
  constexpr cxsmpl<double> mdl_I91x33 = mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I91x36 = mdl_yd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I92x33 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I92x36 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I92x63 = mdl_CKM3x3 * mdl_Ru3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I92x66 = mdl_CKM3x3 * mdl_Ru6x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I93x33 = mdl_Rn3x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I93x36 = mdl_Rn3x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I94x33 = mdl_Rd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I94x36 = mdl_Rd3x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I94x63 = mdl_Rd6x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I94x66 = mdl_Rd6x3 * mdl_conjg__CKM3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I95x33 = mdl_Rl3x3 * mdl_conjg__Rn3x3;
  constexpr cxsmpl<double> mdl_I95x36 = mdl_Rl6x3 * mdl_conjg__Rn3x3;
  constexpr cxsmpl<double> mdl_I96x33 = mdl_Rd3x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I96x36 = mdl_Rd6x3 * mdl_conjg__Rd3x3;
  constexpr cxsmpl<double> mdl_I96x63 = mdl_Rd3x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I96x66 = mdl_Rd6x3 * mdl_conjg__Rd6x3;
  constexpr cxsmpl<double> mdl_I97x33 = mdl_Rl3x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I97x36 = mdl_Rl6x3 * mdl_conjg__Rl3x3;
  constexpr cxsmpl<double> mdl_I97x63 = mdl_Rl3x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I97x66 = mdl_Rl6x3 * mdl_conjg__Rl6x3;
  constexpr cxsmpl<double> mdl_I98x33 = mdl_Ru3x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I98x36 = mdl_Ru6x3 * mdl_conjg__Ru3x3;
  constexpr cxsmpl<double> mdl_I98x63 = mdl_Ru3x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I98x66 = mdl_Ru6x3 * mdl_conjg__Ru6x3;
  constexpr cxsmpl<double> mdl_I99x33 = mdl_ye3x3;
  constexpr cxsmpl<double> mdl_complexi = cxsmpl<double>( 0., 1. );
  constexpr double mdl_sqrt__2 = constexpr_sqrt( 2. );
  constexpr double mdl_sw__exp__2 = ( ( mdl_sw ) * ( mdl_sw ) );
  constexpr cxsmpl<double> mdl_conjg__NN1x1 = conj( mdl_NN1x1 );
  constexpr cxsmpl<double> mdl_conjg__NN1x2 = conj( mdl_NN1x2 );
  constexpr cxsmpl<double> mdl_conjg__NN1x3 = conj( mdl_NN1x3 );
  constexpr cxsmpl<double> mdl_conjg__NN1x4 = conj( mdl_NN1x4 );
  constexpr cxsmpl<double> mdl_conjg__NN2x1 = conj( mdl_NN2x1 );
  constexpr cxsmpl<double> mdl_conjg__NN2x2 = conj( mdl_NN2x2 );
  constexpr cxsmpl<double> mdl_conjg__NN2x3 = conj( mdl_NN2x3 );
  constexpr cxsmpl<double> mdl_conjg__NN2x4 = conj( mdl_NN2x4 );
  constexpr cxsmpl<double> mdl_conjg__NN3x1 = conj( mdl_NN3x1 );
  constexpr cxsmpl<double> mdl_conjg__NN3x2 = conj( mdl_NN3x2 );
  constexpr cxsmpl<double> mdl_conjg__NN3x3 = conj( mdl_NN3x3 );
  constexpr cxsmpl<double> mdl_conjg__NN3x4 = conj( mdl_NN3x4 );
  constexpr cxsmpl<double> mdl_conjg__NN4x1 = conj( mdl_NN4x1 );
  constexpr cxsmpl<double> mdl_conjg__NN4x2 = conj( mdl_NN4x2 );
  constexpr cxsmpl<double> mdl_conjg__NN4x3 = conj( mdl_NN4x3 );
  constexpr cxsmpl<double> mdl_conjg__NN4x4 = conj( mdl_NN4x4 );
  constexpr cxsmpl<double> mdl_conjg__UU1x1 = conj( mdl_UU1x1 );
  constexpr cxsmpl<double> mdl_conjg__UU1x2 = conj( mdl_UU1x2 );
  constexpr cxsmpl<double> mdl_conjg__UU2x1 = conj( mdl_UU2x1 );
  constexpr cxsmpl<double> mdl_conjg__UU2x2 = conj( mdl_UU2x2 );
  constexpr cxsmpl<double> mdl_conjg__VV1x1 = conj( mdl_VV1x1 );
  constexpr cxsmpl<double> mdl_conjg__VV1x2 = conj( mdl_VV1x2 );
  constexpr cxsmpl<double> mdl_conjg__VV2x1 = conj( mdl_VV2x1 );
  constexpr cxsmpl<double> mdl_conjg__VV2x2 = conj( mdl_VV2x2 );
  constexpr double mdl_cos__alp = cos( mdl_alp );
  constexpr double mdl_sin__alp = sin( mdl_alp );
  constexpr cxsmpl<double> mdl_conjg__MUH = conj( mdl_MUH );
  constexpr double mdl_ee = 2. * constexpr_sqrt( 1. / aEWM1 ) * constexpr_sqrt( M_PI );
  constexpr double mdl_gp = mdl_ee / mdl_cw;
  constexpr double mdl_gw = mdl_ee / mdl_sw;
  constexpr double mdl_vev = ( 2. * mdl_cw * mdl_MZ * mdl_sw ) / mdl_ee;
  constexpr double mdl_vd = mdl_vev * mdl_cos__beta;
  constexpr double mdl_vu = mdl_vev * mdl_sin__beta;
  constexpr double mdl_ee__exp__2 = ( ( mdl_ee ) * ( mdl_ee ) );

  if( mdl_Mneu2 < 0 )
    mdl_Wneu2 = -abs( mdl_Wneu2 );
  if( mdl_Mneu3 < 0 )
    mdl_Wneu3 = -abs( mdl_Wneu3 );
  if( mdl_Mneu4 < 0 )
    mdl_Wneu4 = -abs( mdl_Wneu4 );
  if( mdl_Mgo < 0 )
    mdl_Wgo = -abs( mdl_Wgo );
  // Model couplings independent of aS
  // (none)

  // Model parameters dependent on aS
  //constexpr double mdl_sqrt__aS = //constexpr_sqrt( aS ); // now computed event-by-event (running alphas #373)
  //constexpr double G = 2. * mdl_sqrt__aS * //constexpr_sqrt( M_PI ); // now computed event-by-event (running alphas #373)
  //constexpr cxsmpl<double> mdl_G__exp__2 = ( ( G ) * ( G ) ); // now computed event-by-event (running alphas #373)

  // Model couplings dependent on aS
  //constexpr cxsmpl<double> GC_6 = -G; // now computed event-by-event (running alphas #373)
  //constexpr cxsmpl<double> GC_51 = -( mdl_complexi * G * mdl_I51x11 ); // now computed event-by-event (running alphas #373)

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

namespace Parameters_MSSM_SLHA2_dependentCouplings
{
  constexpr size_t ndcoup = 2; // #couplings that vary event by event because they depend on the running alphas QCD
  constexpr size_t idcoup_GC_6 = 0;
  constexpr size_t idcoup_GC_51 = 1;
  struct DependentCouplings_sv
  {
    cxtype_sv GC_6;
    cxtype_sv GC_51;
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
    using namespace Parameters_MSSM_SLHA2;
#endif
    // NB: hardcode cxtype cI(0,1) instead of cxtype (or hardcoded cxsmpl) mdl_complexi (which exists in Parameters_MSSM_SLHA2) because:
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
      constexpr cxsmpl<double> mdl_G__exp__2 = ( ( G ) * ( G ) );
      // Model couplings dependent on aS
      out.GC_6 = -G;
      out.GC_51 = -( cI * G * mdl_I51x11 );
    }
#else
    // ** NB #439: special handling is necessary ONLY FOR VECTORS OF FLOATS (variable Gs are vector floats, fixed parameters are scalar doubles)
    // Use an explicit loop to avoid <<error: conversion of scalar ‘double’ to vector ‘fptype_sv’ {aka ‘__vector(8) float’} involves truncation>>
    // Problems may come e.g. in EFTs from multiplying a vector float (related to aS-dependent G) by a scalar double (aS-independent parameters)
    fptype_v GC_6r_v;
    fptype_v GC_6i_v;
    fptype_v GC_51r_v;
    fptype_v GC_51i_v;
    for( int i = 0; i < neppV; i++ )
    {
      const fptype& G = G_sv[i];
      // Model parameters dependent on aS
      //const fptype mdl_sqrt__aS = constexpr_sqrt( aS );
      //const fptype G = 2. * mdl_sqrt__aS * constexpr_sqrt( M_PI );
      constexpr cxsmpl<double> mdl_G__exp__2 = ( ( G ) * ( G ) );
      // Model couplings dependent on aS
      const cxtype GC_6 = -G;
      const cxtype GC_51 = -( cI * G * mdl_I51x11 );
      GC_6r_v[i] = cxreal( GC_6 );
      GC_6i_v[i] = cximag( GC_6 );
      GC_51r_v[i] = cxreal( GC_51 );
      GC_51i_v[i] = cximag( GC_51 );
    }
    out.GC_6 = cxtype_v( GC_6r_v, GC_6i_v );
    out.GC_51 = cxtype_v( GC_51r_v, GC_51i_v );
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

namespace Parameters_MSSM_SLHA2_independentCouplings
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
    using namespace Parameters_MSSM_SLHA2_dependentCouplings;
    const fptype_sv& gs_sv = G_ACCESS::kernelAccessConst( gs );
    DependentCouplings_sv couplings_sv = computeDependentCouplings_fromG( gs_sv );
    fptype* GC_6s = C_ACCESS::idcoupAccessBuffer( couplings, idcoup_GC_6 );
    fptype* GC_51s = C_ACCESS::idcoupAccessBuffer( couplings, idcoup_GC_51 );
    cxtype_sv_ref GC_6s_sv = C_ACCESS::kernelAccess( GC_6s );
    cxtype_sv_ref GC_51s_sv = C_ACCESS::kernelAccess( GC_51s );
    GC_6s_sv = couplings_sv.GC_6;
    GC_51s_sv = couplings_sv.GC_51;
    mgDebug( 1, __FUNCTION__ );
    return;
  }
#pragma GCC diagnostic pop
}

//==========================================================================

#endif // Parameters_MSSM_SLHA2_H
