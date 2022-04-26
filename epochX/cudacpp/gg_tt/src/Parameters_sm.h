//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.3.1_lo_vect, 2022-01-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_sm_H
#define Parameters_sm_H

#include "mgOnGpuCxtypes.h"
#include "mgOnGpuVectors.h"

#ifndef MGONGPU_HARDCODE_CIPD

#include "read_slha.h"

//constexpr double G = 2. * mdl_sqrt__aS * sqrtNR( M_PI ); // now computed event-by-event (running alphas #373)

namespace Parameters_sm_dependentCouplings
{
  constexpr size_t idcoup_GC_10 = 0;
  constexpr size_t idcoup_GC_11 = 1;  
  // FIXME? should this use a model-dependent mdl_complexi instead of a hardcoded cxmake(0,1)?
  __host__ __device__ inline const cxtype_sv GC_10_fromG( const fptype_sv& G ){ return -G; }
  __host__ __device__ inline const cxtype_sv GC_11_fromG( const fptype_sv& G ){ return cxmake( 0., 1. ) * G; }  
}

class Parameters_sm
{
public:

  static Parameters_sm* getInstance();

  // Define "zero"
  double zero, ZERO;

  // Model parameters independent of aS
  //double aS; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  double mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau, mdl_ymt, mdl_ymb, mdl_Gf, aEWM1, mdl_MH, mdl_MZ, mdl_MTA, mdl_MT, mdl_MB, mdl_conjg__CKM3x3, mdl_conjg__CKM1x1, mdl_CKM3x3, mdl_MZ__exp__2, mdl_MZ__exp__4, mdl_sqrt__2, mdl_MH__exp__2, mdl_aEW, mdl_MW, mdl_sqrt__aEW, mdl_ee, mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw, mdl_vev, mdl_vev__exp__2, mdl_lam, mdl_yb, mdl_yt, mdl_ytau, mdl_muH, mdl_ee__exp__2, mdl_sw__exp__2, mdl_cw__exp__2;
  cxsmpl<double> mdl_complexi, mdl_I1x33, mdl_I2x33, mdl_I3x33, mdl_I4x33;

  // Model couplings independent of aS
  // (none)

  // Model parameters dependent on aS
  //double mdl_sqrt__aS, G, mdl_G__exp__2; // now computed event-by-event (running alphas #373)

  // Model couplings dependent on aS
  //cxsmpl<double> GC_10, GC_11; // now computed event-by-event (running alphas #373)

  // Set parameters that are unchanged during the run
  void setIndependentParameters( SLHAReader& slha );

  // Set couplings that are unchanged during the run
  void setIndependentCouplings();

  // Set parameters that are changed event by event
  //void setDependentParameters(); // now computed event-by-event (running alphas #373)

  // Set couplings that are changed event by event
  //void setDependentCouplings(); // now done event-by-event (running alphas #373)

  // Print parameters that are unchanged during the run
  void printIndependentParameters();

  // Print couplings that are unchanged during the run
  void printIndependentCouplings();

  // Print parameters that are changed event by event
  //void printDependentParameters(); // now computed event-by-event (running alphas #373)

  // Print couplings that are changed event by event
  //void printDependentCouplings(); // now computed event-by-event (running alphas #373)

private:

  static Parameters_sm* instance;
};

#else

#include <limits>

// Hardcoded constexpr physics parameters
namespace Parameters_sm // keep the same name rather than HardcodedParameters_sm for simplicity
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

  // Model parameters independent of aS
  constexpr double zero = 0;
  constexpr double ZERO = 0;
  constexpr double mdl_WH = 6.382339e-03;
  constexpr double mdl_WW = 2.047600e+00;
  constexpr double mdl_WZ = 2.441404e+00;
  constexpr double mdl_WT = 1.491500e+00;
  constexpr double mdl_ymtau = 1.777000e+00;
  constexpr double mdl_ymt = 1.730000e+02;
  constexpr double mdl_ymb = 4.700000e+00;
  //constexpr double aS = 1.180000e-01; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  constexpr double mdl_Gf = 1.166390e-05;
  constexpr double aEWM1 = 1.325070e+02;
  constexpr double mdl_MH = 1.250000e+02;
  constexpr double mdl_MZ = 9.118800e+01;
  constexpr double mdl_MTA = 1.777000e+00;
  constexpr double mdl_MT = 1.730000e+02;
  constexpr double mdl_MB = 4.700000e+00;
  constexpr double mdl_conjg__CKM3x3 = 1.;
  constexpr double mdl_conjg__CKM1x1 = 1.;
  constexpr double mdl_CKM3x3 = 1.;
  constexpr cxsmpl<double> mdl_complexi = cxsmpl<double>( 0., 1. );
  constexpr double mdl_MZ__exp__2 = ( ( mdl_MZ ) * ( mdl_MZ ) );
  constexpr double mdl_MZ__exp__4 = ( ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) );
  constexpr double mdl_sqrt__2 = sqrtNR( 2. );
  constexpr double mdl_MH__exp__2 = ( ( mdl_MH ) * ( mdl_MH ) );
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
  constexpr double mdl_vev = ( 2. * mdl_MW * mdl_sw ) / mdl_ee;
  constexpr double mdl_vev__exp__2 = ( ( mdl_vev ) * ( mdl_vev ) );
  constexpr double mdl_lam = mdl_MH__exp__2 / ( 2. * mdl_vev__exp__2 );
  constexpr double mdl_yb = ( mdl_ymb * mdl_sqrt__2 ) / mdl_vev;
  constexpr double mdl_yt = ( mdl_ymt * mdl_sqrt__2 ) / mdl_vev;
  constexpr double mdl_ytau = ( mdl_ymtau * mdl_sqrt__2 ) / mdl_vev;
  constexpr double mdl_muH = sqrtNR( mdl_lam * mdl_vev__exp__2 );
  constexpr cxsmpl<double> mdl_I1x33 = mdl_yb * mdl_conjg__CKM3x3;
  constexpr cxsmpl<double> mdl_I2x33 = mdl_yt * mdl_conjg__CKM3x3;
  constexpr cxsmpl<double> mdl_I3x33 = mdl_CKM3x3 * mdl_yt;
  constexpr cxsmpl<double> mdl_I4x33 = mdl_CKM3x3 * mdl_yb;
  constexpr double mdl_ee__exp__2 = ( ( mdl_ee ) * ( mdl_ee ) );
  constexpr double mdl_sw__exp__2 = ( ( mdl_sw ) * ( mdl_sw ) );
  constexpr double mdl_cw__exp__2 = ( ( mdl_cw ) * ( mdl_cw ) );

  // Model couplings independent of aS
  // (none)

  // Model parameters dependent on aS
  //constexpr double mdl_sqrt__aS = sqrtNR( aS ); // now computed event-by-event (running alphas #373)
  //constexpr double G = 2. * mdl_sqrt__aS * sqrtNR( M_PI ); // now computed event-by-event (running alphas #373)
  //constexpr double mdl_G__exp__2 = ( ( G ) * ( G ) ); // now computed event-by-event (running alphas #373)

  // Model couplings dependent on aS
  //constexpr cxsmpl<double> GC_10 = -G; // now computed event-by-event (running alphas #373)
  //constexpr cxsmpl<double> GC_11 = mdl_complexi * G; // now computed event-by-event (running alphas #373)

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

#endif // Parameters_sm_H
