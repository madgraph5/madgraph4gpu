//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
%(info_lines)s
//==========================================================================

#ifndef Parameters_%(model_name)s_H
#define Parameters_%(model_name)s_H

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"
#include "mgOnGpuVectors.h"

//==========================================================================

namespace Parameters_%(model_name)s_dependentCouplings
{
  constexpr size_t ndcoup = %(ndcoup)d; // #couplings that depend on the running alphas QCD and vary event by event
%(idcoup)s
  // FIXME? should this use a model-dependent mdl_complexi instead of a hardcoded cxmake(0,1)?
%(inlinedcoup)s
}

//==========================================================================

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  // Compute the output couplings (e.g. gc10 and gc11) from the input gs
  template<class G_ACCESS, class C_ACCESS>
  __device__ inline void
  G2COUP( const fptype gs[],
          fptype couplings[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& gs_sv = G_ACCESS::kernelAccessConst( gs );
%(dcoupaccessbuffer)s
%(dcoupkernelaccess)s
%(dcoupcompute)s
    mgDebug( 1, __FUNCTION__ );
    return;
  }
}

//==========================================================================

#ifndef MGONGPU_HARDCODE_CIPD

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

#include <limits>

// Hardcoded constexpr physics parameters
namespace Parameters_%(model_name)s // keep the same name rather than HardcodedParameters_%(model_name)s for simplicity
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

#endif // Parameters_%(model_name)s_H
