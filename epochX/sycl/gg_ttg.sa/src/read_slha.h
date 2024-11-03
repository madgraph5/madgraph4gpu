// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Sep 2010) for the MG5aMC CPP backend.
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

#ifndef READ_SLHA_H
#define READ_SLHA_H

#include <map>
#include <string>
#include <sstream>
#include <vector>

// We haven't checked which filesystem to include yet
#ifndef INCLUDE_STD_FILESYSTEM_EXPERIMENTAL

// Check for feature test macro for <filesystem>
#   if defined(__cpp_lib_filesystem)
#       define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 0

// Check for feature test macro for <experimental/filesystem>
#   elif defined(__cpp_lib_experimental_filesystem)
#       define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 1

// We can't check if headers exist...
// Let's assume experimental to be safe
#   elif !defined(__has_include)
#       define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 1

// Check if the header "<filesystem>" exists
#   elif __has_include(<filesystem>)

// If we're compiling on Visual Studio and are not compiling with C++17, we need to use experimental
#       ifdef _MSC_VER

// Check and include header that defines "_HAS_CXX17"
#           if __has_include(<yvals_core.h>)
#               include <yvals_core.h>

// Check for enabled C++17 support
#               if defined(_HAS_CXX17) && _HAS_CXX17
// We're using C++17, so let's use the normal version
#                   define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 0
#               endif
#           endif

// If the marco isn't defined yet, that means any of the other VS specific checks failed, so we need to use experimental
#           ifndef INCLUDE_STD_FILESYSTEM_EXPERIMENTAL
#               define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 1
#           endif

// Not on Visual Studio. Let's use the normal version
#       else // #ifdef _MSC_VER
#           define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 0
#       endif

// Check if the header "<filesystem>" exists
#   elif __has_include(<experimental/filesystem>)
#       define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 1

// Fail if neither header is available with a nice error message
#   else
#       error Could not find system header "<filesystem>" or "<experimental/filesystem>"
#   endif

// We priously determined that we need the exprimental version
#   if INCLUDE_STD_FILESYSTEM_EXPERIMENTAL
// Include it
#       include <experimental/filesystem>

// We need the alias from std::experimental::filesystem to std::filesystem
namespace filsys = std::experimental::filesystem;

// We have a decent compiler and can use the normal version
#   else
// Include it
#       include <filesystem>
namespace filsys = std::filesystem;
#   endif

#endif // #ifndef INCLUDE_STD_FILESYSTEM_EXPERIMENTAL

class SLHABlock
{
public:
  SLHABlock(std::string name = ""){_name = name;}
  ~SLHABlock(){}

  void set_entry(std::vector<int> indices, double value);
  double get_entry(std::vector<int> indices, double def_val = 0);
  void set_name(std::string name) {_name = name;}
  std::string get_name(){return _name;}
  int get_indices() { return _indices;}

private:
  std::string _name;
  std::map<std::vector<int>, double> _entries;
  unsigned int _indices;
};

class SLHAReader
{
public:
  SLHAReader(std::string file_name = "", bool verbose=true)
  {if(file_name != "") read_slha_file(file_name, verbose);}

  void read_slha_file(std::string file_name, bool verbose);
  double get_block_entry(std::string block_name, std::vector<int> indices,
                         double def_val = 0);
  double get_block_entry(std::string block_name, int index,
                         double def_val = 0);
  void set_block_entry(std::string block_name, std::vector<int> indices,
                       double value);
  void set_block_entry(std::string block_name, int index,
                       double value);
private:
  std::map<std::string, SLHABlock> _blocks;
};

#endif
