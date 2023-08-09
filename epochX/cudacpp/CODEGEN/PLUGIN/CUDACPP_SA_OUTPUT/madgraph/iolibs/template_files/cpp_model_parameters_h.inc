// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
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

#ifndef MGONGPU_HARDCODE_PARAM // this is only supported in SM processes (e.g. not in EFT models) for the moment (#439)%(efterror)s

#include "read_slha.h"

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
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
  namespace Parameters_%(model_name)s_dependentCouplings
  {
    constexpr size_t ndcoup = %(ndcoup)d; // #couplings that vary event by event because they depend on the running alphas QCD
%(idcoup)s
    struct DependentCouplings_sv
    {
%(dcoupdecl)s
    };
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"  // e.g. <<warning: unused variable ‘mdl_G__exp__2’ [-Wunused-variable]>>
#pragma GCC diagnostic ignored "-Wunused-parameter" // e.g. <<warning: unused parameter ‘G’ [-Wunused-parameter]>>
#ifdef MGONGPUCPP_GPUIMPL
#pragma nv_diagnostic push
#pragma nv_diag_suppress 177 // e.g. <<warning #177-D: variable "mdl_G__exp__2" was declared but never referenced>>
#endif
    __host__ __device__ inline const DependentCouplings_sv computeDependentCouplings_fromG( const fptype_sv& G_sv )
    {
#ifdef MGONGPU_HARDCODE_PARAM
      using namespace Parameters_%(model_name)s;
#endif
      // NB: hardcode cxtype cI(0,1) instead of cxtype (or hardcoded cxsmpl) mdl_complexi (which exists in Parameters_%(model_name)s) because:
      // (1) mdl_complexi is always (0,1); (2) mdl_complexi is undefined in device code; (3) need cxsmpl conversion to cxtype in code below
      const cxtype cI( 0., 1. );
      DependentCouplings_sv out;
%(eftspecial1)s
      {
        const fptype_sv& G = G_sv;
        // Model parameters dependent on aS
%(dcoupsetdpar)s
        // Model couplings dependent on aS
%(dcoupsetdcoup)s
      }
%(eftspecial2)s
      return out;
    }
#ifdef MGONGPUCPP_GPUIMPL
#pragma GCC diagnostic pop
#pragma nv_diagnostic pop
#endif
  }

  //==========================================================================

  namespace Parameters_%(model_name)s_independentCouplings
  {
    constexpr size_t nicoup = %(nicoup)d; // #couplings that are fixed for all events because they do not depend on the running alphas QCD
%(iicoup)s
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
    using namespace Parameters_%(model_name)s_dependentCouplings;
    const fptype_sv& gs_sv = G_ACCESS::kernelAccessConst( gs );
    DependentCouplings_sv couplings_sv = computeDependentCouplings_fromG( gs_sv );
%(dcoupaccessbuffer)s%(dcoupkernelaccess)s%(dcoupcompute)s
    mgDebug( 1, __FUNCTION__ );
    return;
  }
#pragma GCC diagnostic pop

} // end namespace mg5amcGpu/mg5amcCpu

//==========================================================================

#endif // Parameters_%(model_name)s_H
