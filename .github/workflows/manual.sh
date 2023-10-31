#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Oct 2023) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2023) for the MG5aMC CUDACPP plugin.

# Verbose script
###set -x

# Automatic exit on error
set -e

# Usage
function usage() {
  echo "Usage: $(basename $0) codegen|tput"
  exit 1
}

# Check input arguments
if [ "$1" == "codegen" ]; then
  stage=$1
  shift
elif [ "$1" == "tput" ]; then
  stage=$1
  shift
else
  usage
fi

# Exit status
status=0

# Initialise
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "[manual.sh] $stage starting at $(date)"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

echo
echo "Current directory is $(pwd)"
echo "Current git commit is $(git log --oneline -n1 | cut -d' ' -f1)"
topdir=$(pwd)

#echo
#echo "Contents of . (start)"
#ls
#echo "Contents of . (end)"

#echo
#echo "Contents of MG5aMC/mg5amcnlo (start)"
#ls MG5aMC/mg5amcnlo
#echo "Contents of MG5aMC/mg5amcnlo (end)"

# Code generation
if [ "$stage" == "codegen" ]; then
  cd ${topdir}/epochX/cudacpp
  ###processes="$(git ls-tree --name-only HEAD *.mad *.sa)"
  processes="gg_tt.mad" # FOR QUICK TESTS
  for proc in $processes; do
    echo 
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    echo "Code generation for ${proc}"
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    if [ "${proc%.mad}" != "${proc}" ]; then
      # Generate code and check clang formatting
      ./CODEGEN/generateAndCompare.sh -q ${proc%.mad} --mad
    elif [ "${proc%.sa}" != "${proc}" ]; then
      # Generate code and check clang formatting
      ./CODEGEN/generateAndCompare.sh -q ${proc%.sa}
    else
      echo "WARNING! SKIP process directory '${proc}' because it does not end in .mad or .sa"
    fi
    # Check if there are any differences to the current repo
    git checkout HEAD ${proc}/CODEGEN*.txt
    if [ "${proc%.mad}" != "${proc}" ]; then
      git checkout HEAD ${proc}/Cards/me5_configuration.txt
      ###sed -i 's/DEFAULT_F2PY_COMPILER=f2py.*/DEFAULT_F2PY_COMPILER=f2py3/' ${proc}/Source/make_opts
      git checkout HEAD ${proc}/Source/make_opts
    fi
    echo
    echo "git diff (start)"
    git diff --exit-code
    echo "git diff (end)"
  done
fi

# Tput tests (build, runTest.exe, check.exe, gcheck.exe)
if [ "$stage" == "tput" ]; then
  cd ${topdir}/epochX/cudacpp
  # Install ccache
  echo
  echo "wget -q https://github.com/ccache/ccache/releases/download/v4.8.3/ccache-4.8.3-linux-x86_64.tar.xz"
  wget -q https://github.com/ccache/ccache/releases/download/v4.8.3/ccache-4.8.3-linux-x86_64.tar.xz
  echo
  echo "tar -xvf ccache-4.8.3-linux-x86_64.tar.xz"
  tar -xvf ccache-4.8.3-linux-x86_64.tar.xz
  export PATH=$(pwd)/ccache-4.8.3-linux-x86_64:$PATH
  # Configure ccache
  mkdir $(pwd)/CCACHE_DIR
  export CCACHE_DIR=$(pwd)/CCACHE_DIR
  echo
  echo "ccache --version | head -1"
  ccache --version | head -1
  # Dump ccache status before the builds
  echo
  echo "ccache -s (before the builds)"
  ccache -s
  # Configure CXX
  echo
  echo "g++ --version"
  g++ --version
  export CXX=g++
  # Build and test
  export USECCACHE=1 # enable ccache in madgraph4gpu builds
  echo "./tput/teeThroughputX.sh -makej -ggtt -makeclean"
  ./tput/teeThroughputX.sh -makej -ggtt -makeclean
  # Dump ccache status after the builds
  echo
  echo "ccache -s (after the builds)"
  ccache -s
fi

# Finalise
echo
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
if [ $status -eq 0 ]; then # keep this just in case, but 'set -e' causes immediate exit on any error
  echo "[manual.sh] $stage finished with status=$status (OK) at $(date)"
else
  echo "[manual.sh] $stage finished with status=$status (NOT OK) at $(date)"
fi
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
exit $status
