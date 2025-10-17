// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: Z. Wettersten (Oct 2024) for the MG5aMC CUDACPP plugin.

#include "Bridge.h"
#include "CPPProcess.h"
#include "GpuRuntime.h"

#ifndef _FBRIDGE_H_
#define _FBRIDGE_H_

extern "C"
{
#ifdef MGONGPUCPP_GPUIMPL
  using namespace mg5amcGpu;
#else
  using namespace mg5amcCpu;
#endif

  using FORTRANFPTYPE = double;

  void fbridgecreate_( CppObjectInFortran** ppbridge, const int* pnevtF, const int* pnparF, const int* pnp4F );

  void fbridgedelete_( CppObjectInFortran** ppbridge );

  void fbridgesequence_( CppObjectInFortran** ppbridge,
                         const FORTRANFPTYPE* momenta,
                         const FORTRANFPTYPE* gs,
                         const FORTRANFPTYPE* rndhel,
                         const FORTRANFPTYPE* rndcol,
                         const unsigned int* channelIds,
                         FORTRANFPTYPE* mes,
                         int* selhel,
                         int* selcol,
                         const bool* pgoodHelOnly );

  void fbridgesequence_nomultichannel_( CppObjectInFortran** ppbridge,
                                        const FORTRANFPTYPE* momenta,
                                        const FORTRANFPTYPE* gs,
                                        const FORTRANFPTYPE* rndhel,
                                        const FORTRANFPTYPE* rndcol,
                                        FORTRANFPTYPE* mes,
                                        int* selhel,
                                        int* selcol,
                                        const bool* pgoodHelOnly );

  void fbridgegetngoodhel_( CppObjectInFortran** ppbridge, unsigned int* pngoodhel, unsigned int* pntothel );
}
#endif // _FBRIDGE_H_