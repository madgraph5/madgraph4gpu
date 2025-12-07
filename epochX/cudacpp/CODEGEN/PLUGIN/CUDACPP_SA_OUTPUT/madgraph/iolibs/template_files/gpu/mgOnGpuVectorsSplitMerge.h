// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Nov 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#ifndef MGONGPUVECTORSSPLITMERGE_H
#define MGONGPUVECTORSSPLITMERGE_H 1

#include "mgOnGpuVectors.h"

//==========================================================================

#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
namespace mg5amcCpu
{
  //--------------------------------------------------------------------------

  inline fptype2_v
  fpvmerge( const fptype_v& v1, const fptype_v& v2 )
  {
    // This code is not very efficient! It makes mixed precision FFV/color not faster than double on C++ (#537).
    // I considered various alternatives, including
    // - in gcc12 and clang, __builtin_shufflevector (works with different vector lengths, BUT the same fptype...)
    // - casting vector(4)double to vector(4)float and then assigning via reinterpret_cast... but how to do the cast?
    // Probably the best solution is intrinsics?
    // - see https://stackoverflow.com/questions/5139363
    // - see https://stackoverflow.com/questions/54518744
    /*
    fptype2_v out;
    for( int ieppV = 0; ieppV < neppV; ieppV++ )
    {
      out[ieppV] = v1[ieppV];
      out[ieppV+neppV] = v2[ieppV];
    }
    return out;
    */
#if MGONGPU_CPPSIMD == 2
    fptype2_v out =
      { (fptype2)v1[0], (fptype2)v1[1], (fptype2)v2[0], (fptype2)v2[1] };
#elif MGONGPU_CPPSIMD == 4
    fptype2_v out =
      { (fptype2)v1[0], (fptype2)v1[1], (fptype2)v1[2], (fptype2)v1[3], (fptype2)v2[0], (fptype2)v2[1], (fptype2)v2[2], (fptype2)v2[3] };
#elif MGONGPU_CPPSIMD == 8
    fptype2_v out =
      { (fptype2)v1[0], (fptype2)v1[1], (fptype2)v1[2], (fptype2)v1[3], (fptype2)v1[4], (fptype2)v1[5], (fptype2)v1[6], (fptype2)v1[7], (fptype2)v2[0], (fptype2)v2[1], (fptype2)v2[2], (fptype2)v2[3], (fptype2)v2[4], (fptype2)v2[5], (fptype2)v2[6], (fptype2)v2[7] };
#endif
    return out;
  }

  //--------------------------------------------------------------------------

  inline fptype_v
  fpvsplit0( const fptype2_v& v )
  {
    /*
    fptype_v out = {}; // see #594
    for( int ieppV = 0; ieppV < neppV; ieppV++ )
    {
      out[ieppV] = v[ieppV];
    }
    */
#if MGONGPU_CPPSIMD == 2
    fptype_v out =
      { (fptype)v[0], (fptype)v[1] };
#elif MGONGPU_CPPSIMD == 4
    fptype_v out =
      { (fptype)v[0], (fptype)v[1], (fptype)v[2], (fptype)v[3] };
#elif MGONGPU_CPPSIMD == 8
    fptype_v out =
      { (fptype)v[0], (fptype)v[1], (fptype)v[2], (fptype)v[3], (fptype)v[4], (fptype)v[5], (fptype)v[6], (fptype)v[7] };
#endif
    return out;
  }

  //--------------------------------------------------------------------------

  inline fptype_v
  fpvsplit1( const fptype2_v& v )
  {
    /*
    fptype_v out = {}; // see #594
    for( int ieppV = 0; ieppV < neppV; ieppV++ )
    {
      out[ieppV] = v[ieppV+neppV];
    }
    */
#if MGONGPU_CPPSIMD == 2
    fptype_v out =
      { (fptype)v[2], (fptype)v[3] };
#elif MGONGPU_CPPSIMD == 4
    fptype_v out =
      { (fptype)v[4], (fptype)v[5], (fptype)v[6], (fptype)v[7] };
#elif MGONGPU_CPPSIMD == 8
    fptype_v out =
      { (fptype)v[8], (fptype)v[9], (fptype)v[10], (fptype)v[11], (fptype)v[12], (fptype)v[13], (fptype)v[14], (fptype)v[15] };
#endif
    return out;
  }

  //--------------------------------------------------------------------------
}
#endif

#endif // MGONGPUVECTORSSPLITMERGE_H
