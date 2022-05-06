//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
%(info_lines)s
//==========================================================================

#ifndef Parameters_%(model_name)s_H
#define Parameters_%(model_name)s_H

#include "mgOnGpuCxtypes.h"

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

  static Parameters_%(model_name)s* instance;
};

#else

#include <cassert>
#include <limits>

// Hardcoded constexpr physics parameters
namespace Parameters_%(model_name)s // keep the same name rather than HardcodedParameters_%(model_name)s for simplicity
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
  %(hardcoded_independent_parameters)s
  // Model couplings independent of aS
  %(hardcoded_independent_couplings)s
  // Model parameters dependent on aS
  %(hardcoded_dependent_parameters)s
  // Model couplings dependent on aS
  %(hardcoded_dependent_couplings)s
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
