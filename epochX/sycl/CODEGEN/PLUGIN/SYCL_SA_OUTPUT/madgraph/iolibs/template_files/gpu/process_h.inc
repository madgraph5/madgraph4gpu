// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: O. Mattelaer, S. Roiser, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.
//==========================================================================
// Copyright (C) 2021-2023 Argonne National Laboratory.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.
//==========================================================================
// This file has been automatically generated for SYCL standalone by
%(info_lines)s
//==========================================================================

#ifndef MG5_%(process_file_name)s_H
#define MG5_%(process_file_name)s_H

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include "Parameters_%(model_name)s.h"

#include <vector>

//--------------------------------------------------------------------------

namespace Proc
{
  namespace dependentCouplings = Parameters_%(model_name)s_dependentCouplings;
  namespace independentCouplings = Parameters_%(model_name)s_independentCouplings;

  template <typename T>
%(all_helicities)s
  };

#ifdef MGONGPU_HARDCODE_PARAM
  template <typename FPType>
%(independent_parameters_hrdcod)s
#endif
}

%(process_class_definitions)s
  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin_getGoodHel( const vector4* __restrict__ allmomenta, // input: momenta[nevt*npar*4]
                            bool* isGoodHel,                        // output: isGoodHel[ncomb] - device array
                            const signed char* __restrict__ cHel,
                            const cxtype_sv* __restrict__ COUPs,
                            const fptype* __restrict__ cIPD
                            );

  //--------------------------------------------------------------------------

  size_t sigmaKin_setGoodHel( const bool* isGoodHel, size_t* goodHel ); // input: isGoodHel[ncomb] - host array

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  fptype_sv sigmaKin( const vector4* __restrict__ allmomenta, // input: momenta[]
                      const fptype_sv* __restrict__ rndhel,   // input: random numbers[] for helicity selection
                      const fptype_sv* __restrict__ rndcol,   // input: random numbers[] for color selection
                      int_sv* __restrict__ selhel,            // output: helicity selection[]
                      int_sv* __restrict__ selcol,            // output: color selection[]
                      #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                          const size_t channelId,             // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
                      #endif
                      const signed char* __restrict__ cHel,
                      const cxtype_sv* __restrict__ COUPs,
                      const fptype* __restrict__ cIPD,
                      const size_t* __restrict__ cNGoodHel,
                      const size_t* __restrict__ cGoodHel
                    );

  //--------------------------------------------------------------------------
}

#endif // MG5_%(process_file_name)s_H
