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

#include "mgKokkosConfig.h"
#include "mgKokkosTypes.h"
#include "Kokkos_Core.hpp"

#include "Parameters_%(model_name)s.h"

namespace Proc
{

%(process_class_definitions)s


%(process_function_definitions)s

} // end namespace Proc

#endif // MG5_%(process_file_name)s_H