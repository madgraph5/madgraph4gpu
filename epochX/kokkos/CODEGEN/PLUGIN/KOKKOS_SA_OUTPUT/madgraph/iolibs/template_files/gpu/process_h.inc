//==========================================================================
// This file has been automatically generated for C++ Standalone by
%(info_lines)s
//==========================================================================

// %%(helamps_h)s

#ifndef MG5_%(process_file_name)s_H
#define MG5_%(process_file_name)s_H

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "Kokkos_Core.hpp"
#include "Parameters_%(model_name)s.h"

//--------------------------------------------------------------------------

namespace dependentCouplings = Parameters_%(model_name)s_dependentCouplings;
namespace independentCouplings = Parameters_%(model_name)s_independentCouplings;

template <typename T>
%(all_helicities)s
};

%(process_class_definitions)s


%(process_function_definitions)s


#endif // MG5_%(process_file_name)s_H