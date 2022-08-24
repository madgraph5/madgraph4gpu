//==========================================================================
// This file has been automatically generated for %(output_name)s by
%(info_lines)s
//==========================================================================

#ifndef HelAmps_%(model_name)s_H
#define HelAmps_%(model_name)s_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuVectors.h"

#include "Parameters_%(model_name)s.h"

//#include <cmath>
//#include <cstdlib>
//#include <iomanip>
//#include <iostream>

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
%(function_declarations)s
