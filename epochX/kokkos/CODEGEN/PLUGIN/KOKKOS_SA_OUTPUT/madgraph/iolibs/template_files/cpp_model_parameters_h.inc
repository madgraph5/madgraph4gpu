//==========================================================================
// This file has been automatically generated for Kokkos standalone by
%(info_lines)s
//==========================================================================

#ifndef Parameters_%(model_name)s_H
#define Parameters_%(model_name)s_H

#include "mgOnGpuTypes.h"
#include "Kokkos_Complex.hpp"


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


//==========================================================================

namespace Parameters_%(model_name)s_dependentCouplings
{
  constexpr size_t ndcoup = %(ndcoup)d; // #couplings that vary event by event because they depend on the running alphas QCD
%(idcoup)s

  template <typename CXType, typename FPType>
  KOKKOS_INLINE_FUNCTION void set_couplings_from_G( CXType* couplings, const FPType G ) {

    // NB: hardcode cxtype cI(0,1) instead of cxtype (or hardcoded cxsmpl) mdl_complexi (which exists in Parameters_%(model_name)s) because:
    // (1) mdl_complexi is always (0,1); (2) mdl_complexi is undefined in device code; (3) need cxsmpl conversion to cxtype in code below
    static CXType cI( 0., 1. );

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

}

#endif // Parameters_%(model_name)s_H
