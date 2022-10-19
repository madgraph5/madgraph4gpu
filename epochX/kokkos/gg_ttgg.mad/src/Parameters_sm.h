//==========================================================================
// This file has been automatically generated for Kokkos standalone by
// MadGraph5_aMC@NLO v. 3.4.0_lo_vect, 2022-05-06
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_sm_H
#define Parameters_sm_H

#include "mgOnGpuTypes.h"
#include "Kokkos_Complex.hpp"


#include "read_slha.h"

class Parameters_sm
{
public:

static Parameters_sm* getInstance();

// Define "zero"
double zero, ZERO;

// Model parameters independent of aS
  //double aS; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  double mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau, mdl_ymt, mdl_ymb, mdl_Gf, aEWM1, mdl_MH, mdl_MZ, mdl_MTA, mdl_MT, mdl_MB, mdl_conjg__CKM3x3, mdl_conjg__CKM1x1, mdl_CKM3x3, mdl_MZ__exp__2, mdl_MZ__exp__4, mdl_sqrt__2, mdl_MH__exp__2, mdl_aEW, mdl_MW, mdl_sqrt__aEW, mdl_ee, mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw, mdl_vev, mdl_vev__exp__2, mdl_lam, mdl_yb, mdl_yt, mdl_ytau, mdl_muH, mdl_ee__exp__2, mdl_sw__exp__2, mdl_cw__exp__2;
  Kokkos::complex<double> mdl_complexi, mdl_I1x33, mdl_I2x33, mdl_I3x33, mdl_I4x33;

// Model couplings independent of aS
  // (none)

// Model parameters dependent on aS
  //double mdl_sqrt__aS, G, mdl_G__exp__2; // now computed event-by-event (running alphas #373)

// Model couplings dependent on aS
  //Kokkos::complex<double> GC_10, GC_11, GC_12; // now computed event-by-event (running alphas #373)

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

  static Parameters_sm* instance;
};


//==========================================================================

namespace Parameters_sm_dependentCouplings
{
  constexpr size_t ndcoup = 3; // #couplings that vary event by event because they depend on the running alphas QCD
  constexpr size_t idcoup_GC_10 = 0;
  constexpr size_t idcoup_GC_11 = 1;
  constexpr size_t idcoup_GC_12 = 2;

  template <typename CXType, typename FPType>
  KOKKOS_INLINE_FUNCTION void set_couplings_from_G( CXType* couplings, const FPType G ) {

    // NB: hardcode cxtype cI(0,1) instead of cxtype (or hardcoded cxsmpl) mdl_complexi (which exists in Parameters_sm) because:
    // (1) mdl_complexi is always (0,1); (2) mdl_complexi is undefined in device code; (3) need cxsmpl conversion to cxtype in code below
    // static constexpr CXType cI( 0., 1. ); FIXME Kokkos::complex does not have constexpr initializer
    CXType cI( 0., 1. );

      // Model parameters dependent on aS
      //const FPType mdl_sqrt__aS = constexpr_sqrt( aS );
      //const FPType G = 2.*mdl_sqrt__aS*constexpr_sqrt( M_PI );
      const FPType mdl_G__exp__2 = ( ( G )*( G ) );
      // Model couplings dependent on aS
      couplings[idcoup_GC_10] = -G;
      couplings[idcoup_GC_11] = cI*G;
      couplings[idcoup_GC_12] = cI*mdl_G__exp__2;
    }
}

//==========================================================================

namespace Parameters_sm_independentCouplings
{
  constexpr size_t nicoup = 0; // #couplings that are fixed for all events because they do not depend on the running alphas QCD
  // NB: there are no aS-independent couplings in this physics process

}

#endif // Parameters_sm_H
