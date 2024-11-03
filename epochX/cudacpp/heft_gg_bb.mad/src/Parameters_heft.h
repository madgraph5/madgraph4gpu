// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2021-2024) for the MG5aMC CUDACPP plugin.
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.6.0, 2024-09-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_heft_H
#define Parameters_heft_H

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"
#include "mgOnGpuVectors.h"

#include "constexpr_math.h"

//==========================================================================

// AV Jan 2024 (PR #625): this ugly #define was the only way I found to avoid creating arrays[nBsm] in CPPProcess.cc if nBsm is 0
// The problem is that nBsm is determined when generating Parameters.h, which happens after CPPProcess.cc has already been generated
// For simplicity, keep this code hardcoded also for SM processes (a nullptr is needed as in the case nBsm == 0)
#define MGONGPUCPP_NBSMINDEPPARAM_GT_0 1

#ifndef MGONGPU_HARDCODE_PARAM

#include "read_slha.h"

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  class Parameters_heft
  {
  public:

    static Parameters_heft* getInstance();

    // Define "zero"
    double zero, ZERO;

    // Model parameters independent of aS
    //double aS; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
    double mdl_WH1, mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau, mdl_ymt, mdl_ymb, mdl_Gf, aEWM1, mdl_MP, mdl_MH, mdl_MZ, mdl_MTA, mdl_MT, mdl_MB, mdl_conjg__CKM3x3, mdl_MZ__exp__2, mdl_MZ__exp__4, mdl_sqrt__2, mdl_MH__exp__4, mdl_MT__exp__4, mdl_MH__exp__2, mdl_MT__exp__2, mdl_MH__exp__12, mdl_MH__exp__10, mdl_MH__exp__8, mdl_MH__exp__6, mdl_MT__exp__6, mdl_aEW, mdl_MW, mdl_sqrt__aEW, mdl_ee, mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw, mdl_v, mdl_ee__exp__2, mdl_MW__exp__12, mdl_MW__exp__10, mdl_MW__exp__8, mdl_MW__exp__6, mdl_MW__exp__4, mdl_AH, mdl_v__exp__2, mdl_lam, mdl_yb, mdl_yt, mdl_ytau, mdl_muH, mdl_gw__exp__2, mdl_cw__exp__2, mdl_sw__exp__2;
    cxsmpl<double> mdl_complexi;

    // Model couplings independent of aS
    cxsmpl<double> GC_74;

    // Model parameters dependent on aS
    //double mdl_sqrt__aS, G, mdl_G__exp__2, mdl_GH, mdl_Gphi; // now computed event-by-event (running alphas #373)

    // Model couplings dependent on aS
    //cxsmpl<double> GC_13, GC_10, GC_11; // now computed event-by-event (running alphas #373)

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

    // BSM parameters that do not depend on alphaS but are needed in the computation of alphaS-dependent couplings;
    static constexpr int nBsmIndepParam = 9;
    double mdl_bsmIndepParam[nBsmIndepParam];

  private:

    static Parameters_heft* instance;
  };

} // end namespace mg5amcGpu/mg5amcCpu

#else

#include <cassert>
#include <limits>

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  // Hardcoded constexpr physics parameters
  namespace Parameters_heft // keep the same name rather than HardcodedParameters_heft for simplicity
  {
    // Model parameters independent of aS
    constexpr double zero = 0;
    constexpr double ZERO = 0;
    constexpr double mdl_WH1 = 6.382339e-03;
    constexpr double mdl_WH = 6.382339e-03;
    constexpr double mdl_WW = 2.047600e+00;
    constexpr double mdl_WZ = 2.441404e+00;
    constexpr double mdl_WT = 1.491500e+00;
    constexpr double mdl_ymtau = 1.777000e+00;
    constexpr double mdl_ymt = 1.645000e+02;
    constexpr double mdl_ymb = 4.200000e+00;
    //constexpr double aS = 1.180000e-01; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
    constexpr double mdl_Gf = 1.166390e-05;
    constexpr double aEWM1 = 1.325070e+02;
    constexpr double mdl_MP = 1.250001e+02;
    __device__ constexpr double mdl_MH = 1.250000e+02;
    constexpr double mdl_MZ = 9.118800e+01;
    constexpr double mdl_MTA = 1.777000e+00;
    __device__ constexpr double mdl_MT = 1.730000e+02;
    constexpr double mdl_MB = 4.700000e+00;
    constexpr double mdl_conjg__CKM3x3 = 1.;
    constexpr cxsmpl<double> mdl_complexi = cxsmpl<double>( 0., 1. );
    constexpr double mdl_MZ__exp__2 = ( ( mdl_MZ ) * ( mdl_MZ ) );
    constexpr double mdl_MZ__exp__4 = ( ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) );
    constexpr double mdl_sqrt__2 = constexpr_sqrt( 2. );
    __device__ constexpr double mdl_MH__exp__4 = ( ( mdl_MH ) * ( mdl_MH ) * ( mdl_MH ) * ( mdl_MH ) );
    __device__ constexpr double mdl_MT__exp__4 = ( ( mdl_MT ) * ( mdl_MT ) * ( mdl_MT ) * ( mdl_MT ) );
    __device__ constexpr double mdl_MH__exp__2 = ( ( mdl_MH ) * ( mdl_MH ) );
    __device__ constexpr double mdl_MT__exp__2 = ( ( mdl_MT ) * ( mdl_MT ) );
    constexpr double mdl_MH__exp__12 = constexpr_pow( mdl_MH, 12. );
    constexpr double mdl_MH__exp__10 = constexpr_pow( mdl_MH, 10. );
    constexpr double mdl_MH__exp__8 = constexpr_pow( mdl_MH, 8. );
    __device__ constexpr double mdl_MH__exp__6 = constexpr_pow( mdl_MH, 6. );
    __device__ constexpr double mdl_MT__exp__6 = constexpr_pow( mdl_MT, 6. );
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
    __device__ constexpr double mdl_v = ( 2. * mdl_MW * mdl_sw ) / mdl_ee;
    constexpr double mdl_ee__exp__2 = ( ( mdl_ee ) * ( mdl_ee ) );
    constexpr double mdl_MW__exp__12 = constexpr_pow( mdl_MW, 12. );
    constexpr double mdl_MW__exp__10 = constexpr_pow( mdl_MW, 10. );
    constexpr double mdl_MW__exp__8 = constexpr_pow( mdl_MW, 8. );
    constexpr double mdl_MW__exp__6 = constexpr_pow( mdl_MW, 6. );
    constexpr double mdl_MW__exp__4 = ( ( mdl_MW ) * ( mdl_MW ) * ( mdl_MW ) * ( mdl_MW ) );
    constexpr double mdl_AH = ( 47. * mdl_ee__exp__2 * ( 1. - ( 2. * mdl_MH__exp__4 ) / ( 987. * mdl_MT__exp__4 ) - ( 14. * mdl_MH__exp__2 ) / ( 705. * mdl_MT__exp__2 ) + ( 213. * mdl_MH__exp__12 ) / ( 2.634632e7 * mdl_MW__exp__12 ) + ( 5. * mdl_MH__exp__10 ) / ( 119756. * mdl_MW__exp__10 ) + ( 41. * mdl_MH__exp__8 ) / ( 180950. * mdl_MW__exp__8 ) + ( 87. * mdl_MH__exp__6 ) / ( 65800. * mdl_MW__exp__6 ) + ( 57. * mdl_MH__exp__4 ) / ( 6580. * mdl_MW__exp__4 ) + ( 33. * mdl_MH__exp__2 ) / ( 470. * mdl_MW__exp__2 ) ) ) / ( 72. * ( ( M_PI ) * ( M_PI ) ) * mdl_v );
    constexpr double mdl_v__exp__2 = ( ( mdl_v ) * ( mdl_v ) );
    constexpr double mdl_lam = mdl_MH__exp__2 / ( 2. * mdl_v__exp__2 );
    constexpr double mdl_yb = ( mdl_ymb * mdl_sqrt__2 ) / mdl_v;
    constexpr double mdl_yt = ( mdl_ymt * mdl_sqrt__2 ) / mdl_v;
    constexpr double mdl_ytau = ( mdl_ymtau * mdl_sqrt__2 ) / mdl_v;
    constexpr double mdl_muH = constexpr_sqrt( mdl_lam * mdl_v__exp__2 );
    constexpr double mdl_gw__exp__2 = ( ( mdl_gw ) * ( mdl_gw ) );
    constexpr double mdl_cw__exp__2 = ( ( mdl_cw ) * ( mdl_cw ) );
    constexpr double mdl_sw__exp__2 = ( ( mdl_sw ) * ( mdl_sw ) );

    // Model couplings independent of aS
    constexpr cxsmpl<double> GC_74 = -( ( mdl_complexi * mdl_yb ) / mdl_sqrt__2 );

    // Model parameters dependent on aS
    //constexpr double mdl_sqrt__aS = constexpr_sqrt( aS ); // now computed event-by-event (running alphas #373)
    //constexpr double G = 2. * mdl_sqrt__aS * constexpr_sqrt( M_PI ); // now computed event-by-event (running alphas #373)
    //constexpr double mdl_G__exp__2 = ( ( G ) * ( G ) ); // now computed event-by-event (running alphas #373)
    //constexpr double mdl_GH = -( mdl_G__exp__2 * ( 1. + ( 13. * mdl_MH__exp__6 ) / ( 16800. * mdl_MT__exp__6 ) + mdl_MH__exp__4 / ( 168. * mdl_MT__exp__4 ) + ( 7. * mdl_MH__exp__2 ) / ( 120. * mdl_MT__exp__2 ) ) ) / ( 12. * ( ( M_PI ) * ( M_PI ) ) * mdl_v ); // now computed event-by-event (running alphas #373)
    //constexpr double mdl_Gphi = -( mdl_G__exp__2 * ( 1. + mdl_MH__exp__6 / ( 560. * mdl_MT__exp__6 ) + mdl_MH__exp__4 / ( 90. * mdl_MT__exp__4 ) + mdl_MH__exp__2 / ( 12. * mdl_MT__exp__2 ) ) ) / ( 8. * ( ( M_PI ) * ( M_PI ) ) * mdl_v ); // now computed event-by-event (running alphas #373)

    // Model couplings dependent on aS
    //constexpr cxsmpl<double> GC_13 = -( mdl_complexi * mdl_GH ); // now computed event-by-event (running alphas #373)
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

    // BSM parameters that do not depend on alphaS but are needed in the computation of alphaS-dependent couplings;
    constexpr int nBsmIndepParam = 9;
    __device__ constexpr double mdl_bsmIndepParam[nBsmIndepParam] = { mdl_MH, mdl_MT, mdl_MH__exp__4, mdl_MT__exp__4, mdl_MH__exp__2, mdl_MT__exp__2, mdl_MH__exp__6, mdl_MT__exp__6, mdl_v };
  }

} // end namespace mg5amcGpu/mg5amcCpu

#endif

//==========================================================================

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  namespace Parameters_heft_dependentCouplings
  {
    constexpr size_t ndcoup = 3; // #couplings that vary event by event because they depend on the running alphas QCD
    constexpr size_t idcoup_GC_13 = 0;
    constexpr size_t idcoup_GC_10 = 1;
    constexpr size_t idcoup_GC_11 = 2;
    struct DependentCouplings_sv
    {
      cxtype_sv GC_13;
      cxtype_sv GC_10;
      cxtype_sv GC_11;
    };
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"        // e.g. <<warning: unused parameter ‘G’ [-Wunused-parameter]>>
#pragma GCC diagnostic ignored "-Wunused-variable"         // e.g. <<warning: unused variable ‘mdl_G__exp__2’ [-Wunused-variable]>>
#pragma GCC diagnostic ignored "-Wunused-but-set-variable" // e.g. <<warning: variable ‘mdl_G__exp__2’ set but not used [-Wunused-but-set-variable]>>
#ifdef MGONGPUCPP_GPUIMPL
#pragma nv_diagnostic push
#pragma nv_diag_suppress 177 // e.g. <<warning #177-D: variable "mdl_G__exp__2" was declared but never referenced>>
#endif
    __host__ __device__ inline const DependentCouplings_sv computeDependentCouplings_fromG( const fptype_sv& G_sv, const double* bsmIndepParamPtr )
    {
#ifdef MGONGPU_HARDCODE_PARAM
      using namespace Parameters_heft;
#else
      const double mdl_MH = bsmIndepParamPtr[0];
      const double mdl_MT = bsmIndepParamPtr[1];
      const double mdl_MH__exp__4 = bsmIndepParamPtr[2];
      const double mdl_MT__exp__4 = bsmIndepParamPtr[3];
      const double mdl_MH__exp__2 = bsmIndepParamPtr[4];
      const double mdl_MT__exp__2 = bsmIndepParamPtr[5];
      const double mdl_MH__exp__6 = bsmIndepParamPtr[6];
      const double mdl_MT__exp__6 = bsmIndepParamPtr[7];
      const double mdl_v = bsmIndepParamPtr[8];
#endif
      // NB: hardcode cxtype cI(0,1) instead of cxtype (or hardcoded cxsmpl) mdl_complexi (which exists in Parameters_heft) because:
      // (1) mdl_complexi is always (0,1); (2) mdl_complexi is undefined in device code; (3) need cxsmpl conversion to cxtype in code below
      const cxtype cI( 0., 1. );
      DependentCouplings_sv out;
#if not( defined MGONGPU_CPPSIMD && defined MGONGPU_FPTYPE_FLOAT )
      // Couplings are (scalar, or vector of) doubles, or scalar floats - default implementation
      {
        const fptype_sv& G = G_sv;
        // Model parameters dependent on aS
        //const fptype_sv mdl_sqrt__aS = constexpr_sqrt( aS );
        //const fptype_sv G = 2. * mdl_sqrt__aS * constexpr_sqrt( M_PI );
        // *** NB Compute all dependent parameters, including aS, in terms of G rather than in terms of aS ***
        const fptype_sv mdl_G__exp__2 = ( ( G ) * ( G ) );
        const fptype_sv mdl_GH = -( mdl_G__exp__2 * ( 1. + ( 13. * mdl_MH__exp__6 ) / ( 16800. * mdl_MT__exp__6 ) + mdl_MH__exp__4 / ( 168. * mdl_MT__exp__4 ) + ( 7. * mdl_MH__exp__2 ) / ( 120. * mdl_MT__exp__2 ) ) ) / ( 12. * ( ( M_PI ) * ( M_PI ) ) * mdl_v );
        const fptype_sv mdl_Gphi = -( mdl_G__exp__2 * ( 1. + mdl_MH__exp__6 / ( 560. * mdl_MT__exp__6 ) + mdl_MH__exp__4 / ( 90. * mdl_MT__exp__4 ) + mdl_MH__exp__2 / ( 12. * mdl_MT__exp__2 ) ) ) / ( 8. * ( ( M_PI ) * ( M_PI ) ) * mdl_v );
        // Model couplings dependent on aS
        out.GC_13 = -( cI * mdl_GH );
        out.GC_10 = -G;
        out.GC_11 = cI * G;
      }
#else
      // Couplings are VECTORS OF FLOATS: #439 special handling is needed (variable Gs are vector floats, fixed parameters are scalar doubles)
      // Use an explicit loop to avoid <<error: conversion of scalar ‘double’ to vector ‘fptype_sv’ {aka ‘__vector(8) float’} involves truncation>>
      // Problems may come e.g. in EFTs from multiplying a vector float (related to aS-dependent G) by a scalar double (aS-independent parameters)
      // (NB in pure SM processes this special handling is not needed, but we keep it here for simplicity, see PR #824)
      fptype_v GC_13r_v;
      fptype_v GC_13i_v;
      fptype_v GC_10r_v;
      fptype_v GC_10i_v;
      fptype_v GC_11r_v;
      fptype_v GC_11i_v;
      for( int i = 0; i < neppV; i++ )
      {
        const fptype& G = G_sv[i];
        // Model parameters dependent on aS
        //const fptype mdl_sqrt__aS = constexpr_sqrt( aS );
        //const fptype G = 2. * mdl_sqrt__aS * constexpr_sqrt( M_PI );
        // *** NB Compute all dependent parameters, including aS, in terms of G rather than in terms of aS ***
        const fptype mdl_G__exp__2 = ( ( G ) * ( G ) );
        const fptype mdl_GH = -( mdl_G__exp__2 * ( 1. + ( 13. * mdl_MH__exp__6 ) / ( 16800. * mdl_MT__exp__6 ) + mdl_MH__exp__4 / ( 168. * mdl_MT__exp__4 ) + ( 7. * mdl_MH__exp__2 ) / ( 120. * mdl_MT__exp__2 ) ) ) / ( 12. * ( ( M_PI ) * ( M_PI ) ) * mdl_v );
        const fptype mdl_Gphi = -( mdl_G__exp__2 * ( 1. + mdl_MH__exp__6 / ( 560. * mdl_MT__exp__6 ) + mdl_MH__exp__4 / ( 90. * mdl_MT__exp__4 ) + mdl_MH__exp__2 / ( 12. * mdl_MT__exp__2 ) ) ) / ( 8. * ( ( M_PI ) * ( M_PI ) ) * mdl_v );
        // Model couplings dependent on aS
        const cxtype GC_13 = -( cI * mdl_GH );
        const cxtype GC_10 = -G;
        const cxtype GC_11 = cI * G;
        GC_13r_v[i] = cxreal( GC_13 );
        GC_13i_v[i] = cximag( GC_13 );
        GC_10r_v[i] = cxreal( GC_10 );
        GC_10i_v[i] = cximag( GC_10 );
        GC_11r_v[i] = cxreal( GC_11 );
        GC_11i_v[i] = cximag( GC_11 );
      }
      out.GC_13 = cxtype_v( GC_13r_v, GC_13i_v );
      out.GC_10 = cxtype_v( GC_10r_v, GC_10i_v );
      out.GC_11 = cxtype_v( GC_11r_v, GC_11i_v );
#endif
      return out;
    }
#ifdef MGONGPUCPP_GPUIMPL
#pragma GCC diagnostic pop
#pragma nv_diagnostic pop
#endif
  }

  //==========================================================================

  namespace Parameters_heft_independentCouplings
  {
    constexpr size_t nicoup = 1; // #couplings that are fixed for all events because they do not depend on the running alphas QCD
    //constexpr size_t ixcoup_GC_74 = 0 + Parameters_heft_dependentCouplings::ndcoup; // out of ndcoup+nicoup
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
          fptype couplings[],
          const double* bsmIndepParamPtr )
  {
    mgDebug( 0, __FUNCTION__ );
    using namespace Parameters_heft_dependentCouplings;
    const fptype_sv& gs_sv = G_ACCESS::kernelAccessConst( gs );
    DependentCouplings_sv couplings_sv = computeDependentCouplings_fromG( gs_sv, bsmIndepParamPtr );
    fptype* GC_13s = C_ACCESS::idcoupAccessBuffer( couplings, idcoup_GC_13 );
    fptype* GC_10s = C_ACCESS::idcoupAccessBuffer( couplings, idcoup_GC_10 );
    fptype* GC_11s = C_ACCESS::idcoupAccessBuffer( couplings, idcoup_GC_11 );
    cxtype_sv_ref GC_13s_sv = C_ACCESS::kernelAccess( GC_13s );
    cxtype_sv_ref GC_10s_sv = C_ACCESS::kernelAccess( GC_10s );
    cxtype_sv_ref GC_11s_sv = C_ACCESS::kernelAccess( GC_11s );
    GC_13s_sv = couplings_sv.GC_13;
    GC_10s_sv = couplings_sv.GC_10;
    GC_11s_sv = couplings_sv.GC_11;
    mgDebug( 1, __FUNCTION__ );
    return;
  }
#pragma GCC diagnostic pop

} // end namespace mg5amcGpu/mg5amcCpu

//==========================================================================

#endif // Parameters_heft_H
