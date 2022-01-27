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

  // Model parameters independent of aS
  constexpr double zero = 0;
  constexpr double ZERO = 0;
  %(hardcoded_independent_parameters)s

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
