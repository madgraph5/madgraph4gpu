//==========================================================================
// This file has been automatically generated for SYCL standalone by
%(info_lines)s
//==========================================================================

#ifndef Parameters_%(model_name)s_H
#define Parameters_%(model_name)s_H

#include "mgOnGpuCxtypes.h"
#include <cmath>

#ifndef MGONGPU_HARDCODE_PARAM // this is only supported in SM processes (e.g. not in EFT models) for the moment (#439)%(efterror)s

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
  //void printDependentParameters(); // now computed event-by-event (running alphas #373)

  // Print couplings that are changed event by event
  //void printDependentCouplings(); // now computed event-by-event (running alphas #373)
}

#endif

//==========================================================================

namespace Parameters_%(model_name)s_dependentCouplings
{
  constexpr size_t ndcoup = %(ndcoup)d; // #couplings that vary event by event because they depend on the running alphas QCD
%(idcoup)s

  template <typename CXType, typename FPType>
  inline void set_couplings_from_G( CXType* couplings, const FPType G ) {
#ifdef MGONGPU_HARDCODE_PARAM
    using namespace Parameters_%(model_name)s;
#endif
    // NB: hardcode cxtype cI(0,1) instead of cxtype (or hardcoded cxsmpl) mdl_complexi (which exists in Parameters_%(model_name)s) because:
    // (1) mdl_complexi is always (0,1); (2) mdl_complexi is undefined in device code; (3) need cxsmpl conversion to cxtype in code below
    #if defined MGONGPU_COMPLEX_CUTHRUST
        const CXType cI( FPType(0.0), FPType(1.0));
    #elif defined MGONGPU_COMPLEX_CUCOMPLEX
        const CXType cI = CXMAKE_SV_2ARG(FPType(0.0), FPType(1.0)); //FIXME CXMAKE_SV_2ARG not included in this header file, perhaps overload assignment operator = in mgOnGpuTypes.h
    #else
        static constexpr CXType cI(FPType(0.0), FPType(1.0));
    #endif

    // Model parameters dependent on aS
%(dcoupsetdpar)s
    // Model couplings dependent on aS
%(dcoupsetdcoup)s
  }
}

//==========================================================================

namespace Parameters_%(model_name)s_independentCouplings
{
  constexpr size_t nicoup = %(nicoup)d; // #couplings that are fixed for all events because they do not depend on the running alphas QCD
%(iicoup)s

#ifdef MGONGPU_HARDCODE_PARAM
  template <typename CXType, typename FPType>
  constexpr CXType independent_couplings[] {
%(icoupseticoup_hrdcod)s
  };
#endif
}

#endif // Parameters_%(model_name)s_H
