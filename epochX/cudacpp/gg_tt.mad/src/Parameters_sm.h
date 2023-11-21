// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.5.2_lo_vect, 2023-11-08
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_sm_H
#define Parameters_sm_H

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"
#include "mgOnGpuVectors.h"

//==========================================================================

#ifndef MGONGPU_HARDCODE_PARAM // this is only supported in SM processes (e.g. not in EFT models) for the moment (#439)

#include "read_slha.h"

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
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

} // end namespace mg5amcGpu/mg5amcCpu

#else

#include <cassert>
#include <limits>

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  // Hardcoded constexpr physics parameters
  namespace Parameters_sm // keep the same name rather than HardcodedParameters_sm for simplicity
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
    constexpr double mdl_sqrt__2 = constexpr_sqrt( 2. );
    constexpr double mdl_MH__exp__2 = ( ( mdl_MH ) * ( mdl_MH ) );
    constexpr double mdl_aEW = 1. / aEWM1;
    constexpr double mdl_MW = constexpr_sqrt( mdl_MZ__exp__2 / 2. + constexpr_sqrt( mdl_MZ__exp__4 / 4. - ( mdl_aEW * M_PI * mdl_MZ__exp__2 ) / ( mdl_Gf * mdl_sqrt__2 ) ) );
    constexpr double mdl_sqrt__aEW = constexpr_sqrt( mdl_aEW );
    constexpr double mdl_ee = 2. * mdl_sqrt__aEW * constexpr_sqrt( M_PI );
    constexpr double mdl_MW__exp__2 = ( ( mdl_MW ) * ( mdl_MW ) );
    constexpr double mdl_sw2 = 1. - mdl_MW__exp__2 / mdl_MZ__exp__2;
    constexpr double mdl_cw = constexpr_sqrt( 1. - mdl_sw2 );
    constexpr double mdl_sqrt__sw2 = constexpr_sqrt( mdl_sw2 );
    constexpr double mdl_sw = mdl_sqrt__sw2;
    constexpr double mdl_g1 = mdl_ee / mdl_cw;
    constexpr double mdl_gw = mdl_ee / mdl_sw;
    constexpr double mdl_vev = ( 2. * mdl_MW * mdl_sw ) / mdl_ee;
    constexpr double mdl_vev__exp__2 = ( ( mdl_vev ) * ( mdl_vev ) );
    constexpr double mdl_lam = mdl_MH__exp__2 / ( 2. * mdl_vev__exp__2 );
    constexpr double mdl_yb = ( mdl_ymb * mdl_sqrt__2 ) / mdl_vev;
    constexpr double mdl_yt = ( mdl_ymt * mdl_sqrt__2 ) / mdl_vev;
    constexpr double mdl_ytau = ( mdl_ymtau * mdl_sqrt__2 ) / mdl_vev;
    constexpr double mdl_muH = constexpr_sqrt( mdl_lam * mdl_vev__exp__2 );
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
    //constexpr double mdl_sqrt__aS = //constexpr_sqrt( aS ); // now computed event-by-event (running alphas #373)
    //constexpr double G = 2. * mdl_sqrt__aS * //constexpr_sqrt( M_PI ); // now computed event-by-event (running alphas #373)
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

} // end namespace mg5amcGpu/mg5amcCpu

#endif

//==========================================================================

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  namespace Parameters_sm_dependentCouplings
  {
    constexpr size_t ndcoup = 2; // #couplings that vary event by event because they depend on the running alphas QCD
    constexpr size_t idcoup_GC_10 = 0;
    constexpr size_t idcoup_GC_11 = 1;
    struct DependentCouplings_sv
    {
      cxtype_sv GC_10;
      cxtype_sv GC_11;
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
      using namespace Parameters_sm;
#endif
      // NB: hardcode cxtype cI(0,1) instead of cxtype (or hardcoded cxsmpl) mdl_complexi (which exists in Parameters_sm) because:
      // (1) mdl_complexi is always (0,1); (2) mdl_complexi is undefined in device code; (3) need cxsmpl conversion to cxtype in code below
      const cxtype cI( 0., 1. );
      DependentCouplings_sv out;
      // Begin SM implementation - no special handling of vectors of floats as in EFT (#439)
      {
        const fptype_sv& G = G_sv;
        // Model parameters dependent on aS
        //const fptype_sv mdl_sqrt__aS = constexpr_sqrt( aS );
        //const fptype_sv G = 2. * mdl_sqrt__aS * constexpr_sqrt( M_PI );
        const fptype_sv mdl_G__exp__2 = ( ( G ) * ( G ) );
        // Model couplings dependent on aS
        out.GC_10 = -G;
        out.GC_11 = cI * G;
      }
      // End SM implementation - no special handling of vectors of floats as in EFT (#439)
      return out;
    }
#ifdef __CUDACC__
#pragma GCC diagnostic pop
#pragma nv_diagnostic pop
#endif
  }

  //==========================================================================

  namespace Parameters_sm_independentCouplings
  {
    constexpr size_t nicoup = 0; // #couplings that are fixed for all events because they do not depend on the running alphas QCD
    // NB: there are no aS-independent couplings in this physics process
  }

  //==========================================================================

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
    using namespace Parameters_sm_dependentCouplings;
    const fptype_sv& gs_sv = G_ACCESS::kernelAccessConst( gs );
    DependentCouplings_sv couplings_sv = computeDependentCouplings_fromG( gs_sv );
    fptype* GC_10s = C_ACCESS::idcoupAccessBuffer( couplings, idcoup_GC_10 );
    fptype* GC_11s = C_ACCESS::idcoupAccessBuffer( couplings, idcoup_GC_11 );
    cxtype_sv_ref GC_10s_sv = C_ACCESS::kernelAccess( GC_10s );
    cxtype_sv_ref GC_11s_sv = C_ACCESS::kernelAccess( GC_11s );
    GC_10s_sv = couplings_sv.GC_10;
    GC_11s_sv = couplings_sv.GC_11;
    mgDebug( 1, __FUNCTION__ );
    return;
  }
#pragma GCC diagnostic pop

} // end namespace mg5amcGpu/mg5amcCpu

//==========================================================================

#endif // Parameters_sm_H
