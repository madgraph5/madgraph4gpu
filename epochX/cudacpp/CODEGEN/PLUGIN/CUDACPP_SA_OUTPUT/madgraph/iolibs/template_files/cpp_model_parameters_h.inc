//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
%(info_lines)s
//==========================================================================

#ifndef Parameters_%(model_name)s_H
#define Parameters_%(model_name)s_H

#include <complex>

#ifndef MGONGPU_HARDCODE_CIPC

#include "read_slha.h"

class Parameters_%(model_name)s
{
public:

  static Parameters_%(model_name)s* getInstance();

  // Define "zero"
  double zero, ZERO;

  %(independent_parameters)s

  %(independent_couplings)s

  %(dependent_parameters)s

  %(dependent_couplings)s

  // Set parameters that are unchanged during the run
  void setIndependentParameters(SLHAReader& slha);

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

  static Parameters_%(model_name)s* instance;

};

#else

#include <limits>

// Hardcoded constexpr physics parameters
namespace Parameters_%(model_name)s // keep the same name rather than HardcodedParameters_%(model_name)s for simplicity
{
  // Constexpr implementation of sqrt (see https://stackoverflow.com/a/34134071)
  double constexpr detailSqrtNewtonRaphson( double x, double curr, double prev )
  {
    return curr == prev ? curr : detailSqrtNewtonRaphson( x, 0.5 * ( curr + x / curr ), curr );
  }
  double constexpr sqrtNR(double x)
  {
    return x >= 0 && x < std::numeric_limits<double>::infinity()
                         ? detailSqrtNewtonRaphson(x, x, 0)
                         : std::numeric_limits<double>::quiet_NaN();
  }

  // EVENTUALLY: %%(hardcoded_independent_parameters)s
  // Model parameters independent of aS (for eemumu and ggttgg)
  constexpr double zero = 0;
  constexpr double ZERO = 0;
  constexpr double mdl_WH = 6.382339e-03;
  constexpr double mdl_WW = 2.047600e+00;
  constexpr double mdl_WZ = 2.441404e+00;
  constexpr double mdl_WT = 1.491500e+00;
  constexpr double mdl_ymtau = 1.777000e+00;
  constexpr double mdl_ymt = 1.730000e+02;
  constexpr double mdl_ymb = 4.700000e+00;
  constexpr double aS = 1.180000e-01;
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
  constexpr std::complex<double> mdl_complexi = std::complex<double>(0.,1.);
  constexpr double mdl_MZ__exp__2 = ((mdl_MZ)*(mdl_MZ));
  constexpr double mdl_MZ__exp__4 = ((mdl_MZ)*(mdl_MZ)*(mdl_MZ)*(mdl_MZ));
  constexpr double mdl_sqrt__2 = sqrtNR(2.);
  constexpr double mdl_MH__exp__2 = ((mdl_MH)*(mdl_MH));
  constexpr double mdl_aEW = 1./aEWM1;
  constexpr double mdl_MW = sqrtNR(mdl_MZ__exp__2/2.+sqrtNR(mdl_MZ__exp__4/4.-(mdl_aEW*M_PI*mdl_MZ__exp__2)/(mdl_Gf*mdl_sqrt__2)));
  constexpr double mdl_sqrt__aEW = sqrtNR(mdl_aEW);
  constexpr double mdl_ee = 2.*mdl_sqrt__aEW*sqrtNR(M_PI);
  constexpr double mdl_MW__exp__2 = ((mdl_MW)*(mdl_MW));
  constexpr double mdl_sw2 = 1.-mdl_MW__exp__2/mdl_MZ__exp__2;
  constexpr double mdl_cw = sqrtNR(1.-mdl_sw2);
  constexpr double mdl_sqrt__sw2 = sqrtNR(mdl_sw2);
  constexpr double mdl_sw = mdl_sqrt__sw2;
  constexpr double mdl_g1 = mdl_ee/mdl_cw;
  constexpr double mdl_gw = mdl_ee/mdl_sw;
  constexpr double mdl_vev = (2.*mdl_MW*mdl_sw)/mdl_ee;
  constexpr double mdl_vev__exp__2 = ((mdl_vev)*(mdl_vev));
  constexpr double mdl_lam = mdl_MH__exp__2/(2.*mdl_vev__exp__2);
  constexpr double mdl_yb = (mdl_ymb*mdl_sqrt__2)/mdl_vev;
  constexpr double mdl_yt = (mdl_ymt*mdl_sqrt__2)/mdl_vev;
  constexpr double mdl_ytau = (mdl_ymtau*mdl_sqrt__2)/mdl_vev;
  constexpr double mdl_muH = sqrtNR(mdl_lam*mdl_vev__exp__2);
  constexpr std::complex<double> mdl_I1x33 = mdl_yb*mdl_conjg__CKM3x3;
  constexpr std::complex<double> mdl_I2x33 = mdl_yt*mdl_conjg__CKM3x3;
  constexpr std::complex<double> mdl_I3x33 = mdl_CKM3x3*mdl_yt;
  constexpr std::complex<double> mdl_I4x33 = mdl_CKM3x3*mdl_yb;
  constexpr double mdl_ee__exp__2 = ((mdl_ee)*(mdl_ee));
  constexpr double mdl_sw__exp__2 = ((mdl_sw)*(mdl_sw));
  constexpr double mdl_cw__exp__2 = ((mdl_cw)*(mdl_cw));

  // EVENTUALLY: %%(hardcoded_dependent_parameters)s
  // Model couplings independent of aS (for eemumu)
  //constexpr std::complex<double> GC_3 = -(mdl_ee*mdl_complexi);
  //constexpr std::complex<double> GC_50 = -(mdl_cw*mdl_ee*mdl_complexi)/(2.*mdl_sw);
  //constexpr std::complex<double> GC_59 = (mdl_ee*mdl_complexi*mdl_sw)/(2.*mdl_cw);
  constexpr std::complex<double> GC_3 = std::complex<double>(0.,-mdl_ee); // does not build
  constexpr std::complex<double> GC_50 = std::complex<double>(0.,-(mdl_cw*mdl_ee)/(2.*mdl_sw)); // does not build
  constexpr std::complex<double> GC_59 = std::complex<double>(0.,(mdl_ee*mdl_sw)/(2.*mdl_cw)); // does not build
  // Model couplings independent of aS (for ggttgg)
  // (none)

  // EVENTUALLY: %%(hardcoded_independent_couplings)s
  // Model parameters dependent on aS (for eemumu and ggttgg)
  constexpr double mdl_sqrt__aS = sqrtNR(aS);
  constexpr double G = 2.*mdl_sqrt__aS*sqrtNR(M_PI);
  constexpr double mdl_G__exp__2 = ((G)*(G));
  // Model couplings dependent on aS (for eemumu)
  // (none)

  // EVENTUALLY: %%(hardcoded_dependent_couplings)s
  // Model couplings dependent on aS (for ggttgg)
  constexpr std::complex<double> GC_10 = -G;
  //constexpr std::complex<double> GC_11 = mdl_complexi*G; // does not build
  //constexpr std::complex<double> GC_12 = mdl_complexi*mdl_G__exp__2; // does not build
  constexpr std::complex<double> GC_11 = std::complex<double>(0.,G);
  constexpr std::complex<double> GC_12 = std::complex<double>(0.,mdl_G__exp__2);

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

#endif // Parameters_%(model_name)s_H
