//==========================================================================
// This file has been automatically generated for C++ Standalone by
%(info_lines)s
//==========================================================================

%(helamps_h)s

#ifndef MG5_%(process_file_name)s_H
#define MG5_%(process_file_name)s_H

#include <vector>
#include <cassert>
#include <iostream>

#include "mgOnGpuConfig.h"
#include "mgOnGpuVectors.h"
#include "mgOnGpuTypes.h"

#include "Kokkos_Core.hpp"

#include "Parameters_%(model_name)s.h"

namespace mg5amcGpu
{

%(process_class_definitions)s


%(process_function_definitions)s

} // end namespace mg5amcGpu

#endif // MG5_%(process_file_name)s_H