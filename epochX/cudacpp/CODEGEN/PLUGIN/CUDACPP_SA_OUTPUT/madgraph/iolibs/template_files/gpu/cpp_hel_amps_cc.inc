//==========================================================================
// This file has been automatically generated for %(output_name)s by
%(info_lines)s
//==========================================================================

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"
#include "HelAmps_sm.h"

//#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

mgDebugDeclare();

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

%(function_definitions)s
} // end namespace
