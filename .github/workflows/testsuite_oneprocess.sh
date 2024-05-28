#!/bin/bash
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Oct 2023) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2023-2024) for the MG5aMC CUDACPP plugin.

# Verbose script
###set -x

# Automatic exit on error
###set -e

# Path to the top directory of madgraphgpu
# In the CI this would be simply $(pwd), but allow the script to be run also outside the CI
echo "Executing $0 $*"; echo
topdir=$(cd $(dirname $0)/../..; pwd)

#----------------------------------------------------------------------------------------------------------------------------------

# Code generation stage
function codegen() {
  if [ "$1" == "" ] || [ "$2" != "" ]; then echo "Usage: $(basename $0) <process>"; exit 1; fi
  proc=$1
  if [ "${proc%.mad}" == "${proc}" ] && [ "${proc%.sa}" == "${proc}" ]; then echo "Usage: $(basename $0) <process.mad|process.sa>"; exit 1; fi
  # Generate code and check clang formatting
  cd ${topdir}/epochX/cudacpp
  echo "Current directory is $(pwd)"
  echo 
  echo "*******************************************************************************"
  echo "*** code generation for ${proc}"
  echo "*******************************************************************************"
  if [ "${proc%.mad}" != "${proc}" ]; then
    ./CODEGEN/generateAndCompare.sh -q ${proc%.mad} --mad
  else
    ./CODEGEN/generateAndCompare.sh -q ${proc%.sa}
  fi
  # Check if there are any differences to the current repo
  ###compare=1 # enable comparison to current git repo
  compare=0 # disable comparison to current git repo
  if [ "${compare}" != "0" ] && [ "$(git ls-tree --name-only HEAD ${proc})" != "" ]; then
    echo
    echo "Compare newly generated code for ${proc} to that in the madgraph4gpu github repository"
    git checkout HEAD ${proc}/CODEGEN*.txt
    if [ "${proc%.mad}" != "${proc}" ]; then
      git checkout HEAD ${proc}/Cards/me5_configuration.txt
      git checkout HEAD ${proc}/Source/make_opts
    fi
    echo "git diff (start)"
    git diff --exit-code
    echo "git diff (end)"
  else
    echo
    echo "(SKIP comparison of newly generated code for ${proc} to that in the madgraph4gpu github repository)"
  fi
}

#----------------------------------------------------------------------------------------------------------------------------------

function setup_ccache {
  # Set up ccache environment
  export PATH=${topdir}/BIN:$PATH
  export CCACHE_DIR=${topdir}/CCACHE_DIR
}

#----------------------------------------------------------------------------------------------------------------------------------

# Before-build stage (analyse data retrieved from cache, download ccache executable and googletest if not retrieved from cache)
function before_build() {
  # Install and configure ccache
  if [ -d ${topdir}/DOWNLOADS ]; then
    echo "Directory ${topdir}/DOWNLOADS already exists (retrieved from cache)"
  else
    echo "Directory ${topdir}/DOWNLOADS does not exist: create it"
    mkdir ${topdir}/DOWNLOADS
    cd ${topdir}/DOWNLOADS
    echo "Current directory is $(pwd)"
    echo
    echo "wget -q https://github.com/ccache/ccache/releases/download/v4.8.3/ccache-4.8.3-linux-x86_64.tar.xz"
    wget -q https://github.com/ccache/ccache/releases/download/v4.8.3/ccache-4.8.3-linux-x86_64.tar.xz
    echo
    echo "tar -xvf ccache-4.8.3-linux-x86_64.tar.xz"
    tar -xvf ccache-4.8.3-linux-x86_64.tar.xz
  fi
  mkdir ${topdir}/BIN
  cd ${topdir}/BIN
  ln -sf ${topdir}/DOWNLOADS/ccache-4.8.3-linux-x86_64/ccache .
  # Set up ccache environment
  setup_ccache
  # Create the CCACHE_DIR directory if it was not retrieved from the cache
  echo
  if [ -d ${CCACHE_DIR} ]; then
    echo "Directory CCACHE_DIR=${CCACHE_DIR} already exists (retrieved from cache)"
  else
    echo "Directory CCACHE_DIR=${CCACHE_DIR} does not exist: create it"
    mkdir ${CCACHE_DIR}
  fi
  # Dump ccache status before the builds
  echo
  echo "ccache --version | head -1"
  ccache --version | head -1
  echo
  echo "CCACHE_DIR=${CCACHE_DIR}"
  echo "du -sm ${CCACHE_DIR}"
  du -sm ${CCACHE_DIR}
  echo
  echo "ccache -s (before the builds)"
  ccache -s
  # Check if googletest has already been installed and configured
  echo
  if [ -d ${topdir}/test/googletest ]; then
    echo "Directory ${topdir}/test/googletest already exists (retrieved from cache)"
    echo "ls ${topdir}/test/googletest (start)"
    ls ${topdir}/test/googletest
    echo "ls ${topdir}/test/googletest (end)"
  else
    echo "Directory ${topdir}/test/googletest does not exist: it will be created during the build"
  fi
}

#----------------------------------------------------------------------------------------------------------------------------------

# Build stage
function build() {
  if [ "$1" == "" ] || [ "$2" != "" ]; then echo "Usage: $(basename $0) <process>"; exit 1; fi
  proc=$1
  if [ "${proc%.mad}" == "${proc}" ] && [ "${proc%.sa}" == "${proc}" ]; then echo "Usage: $(basename $0) <process.mad|process.sa>"; exit 1; fi
  # Set up build environment
  setup_ccache
  export USECCACHE=1 # enable ccache in madgraph4gpu builds
  export CXX=g++ # set up CXX that is needed by cudacpp.mk
  ###echo; echo "$CXX --version"; $CXX --version
  export USEBUILDDIR=1
  # Iterate over P* directories and build
  cd ${topdir}/epochX/cudacpp/${proc}
  echo "Current directory is $(pwd)"
  gtestlibs=0
  pdirs="$(ls -d SubProcesses/P*_*)"
  for pdir in ${pdirs}; do
    pushd $pdir >& /dev/null
    echo
    echo "*******************************************************************************"
    echo "*** build ${proc} ($(basename $(pwd)))"
    echo "*******************************************************************************"
    echo
    echo "Building in $(pwd)"
    if [ "${gtestlibs}" == "0" ]; then
      # Build googletest once and for all to avoid issues in parallel builds
      gtestlibs=1
      make -f cudacpp.mk gtestlibs
    fi
    # NB: 'make bldall' internally checks if 'which nvcc' and 'which hipcc' succeed before attempting to build cuda and hip
    make -j bldall
    popd >& /dev/null
  done
}

#----------------------------------------------------------------------------------------------------------------------------------

# After-build stage (analyse data to be saved in updated cache)
function after_build() {
  # Set up ccache environment
  setup_ccache
  # Dump ccache status after the builds
  echo
  echo "CCACHE_DIR=${CCACHE_DIR}"
  echo "du -sm ${CCACHE_DIR}"
  du -sm ${CCACHE_DIR}
  echo
  echo "ccache -s (after the builds)"
  ccache -s
  # Check contents of googletest
  echo
  echo "ls ${topdir}/test/googletest (start)"
  ls ${topdir}/test/googletest
  echo "ls ${topdir}/test/googletest (end)"
  # Check contents of build directories
  echo
  echo "ls -d ${topdir}/epochX/cudacpp/*.*/SubProcesses/P*_*/build* (start)"
  ls -d ${topdir}/epochX/cudacpp/*.*/SubProcesses/P*_*/build*
  echo "ls -d ${topdir}/epochX/cudacpp/*.*/SubProcesses/P*_*/build* (end)"
}

#----------------------------------------------------------------------------------------------------------------------------------

function runExe() {
  echo
  echo "Execute $*"
  if [ -f $1 ]; then $*; else echo "(SKIP missing $1)"; fi
}

# Tput-test stage (runTest.exe, check.exe, gcheck.exe)
function tput_test() {
  if [ "$1" == "" ] || [ "$2" != "" ]; then echo "Usage: $(basename $0) <process>"; exit 1; fi
  proc=$1
  if [ "${proc%.mad}" == "${proc}" ] && [ "${proc%.sa}" == "${proc}" ]; then echo "Usage: $(basename $0) <process.mad|process.sa>"; exit 1; fi
  # Iterate over P* directories and run tests
  cd ${topdir}/epochX/cudacpp/${proc}
  echo "Current directory is $(pwd)"
  pdirs="$(ls -d SubProcesses/P*_*)"
  for pdir in ${pdirs}; do
    pushd $pdir >& /dev/null
    echo
    echo "*******************************************************************************"
    echo "*** tput-test ${proc} ($(basename $(pwd)))"
    echo "*******************************************************************************"
    echo
    echo "Testing in $(pwd)"
    # FIXME1: this is just a quick test, eventually port here tput tests from throughputX.sh
    # (could move some throughputX.sh functions to a separate script included both here and there)
    # FIXME2: enable FPEs
    # FIXME3: handle all d/f/m, inl0/1, hrd0/1 etc...
    unamep=$(uname -p)
    unames=$(uname -s)
    for backend in cuda hip cppnone cppsse4 cppavx2 cpp512y cpp512z; do
      # Skip GPU tests for NVidia and AMD unless nvcc and hipcc, respectively, are in PATH
      if ! nvcc --version &> /dev/null; then
        if [ "${backend}" == "cuda" ]; then echo; echo "(SKIP ${backend} because nvcc is missing on this node)"; continue; fi
      elif ! hipcc --version &> /dev/null; then
        if [ "${backend}" == "hip" ]; then echo; echo "(SKIP ${backend} because hipcc is missing on this node)"; continue; fi
      fi
      # Skip C++ tests for unsupported simd modes as done in tput tests (prevent illegal instruction crashes #791)
      if [ "${unamep}" != "x86_64" ]; then
        if [ "${backend}" == "cppavx2" ]; then echo; echo "(SKIP ${backend} which is not supported on ${unamep})"; continue; fi
        if [ "${backend}" == "cpp512y" ]; then echo; echo "(SKIP ${backend} which is not supported on ${unamep})"; continue; fi
        if [ "${backend}" == "cpp512z" ]; then echo; echo "(SKIP ${backend} which is not supported on ${unamep})"; continue; fi
      elif [ "${unames}" == "Darwin" ]; then
        if [ "${backend}" == "cpp512y" ]; then echo; echo "(SKIP ${backend} which is not supported on ${unames})"; continue; fi
        if [ "${backend}" == "cpp512z" ]; then echo; echo "(SKIP ${backend} which is not supported on ${unames})"; continue; fi
      elif [ "$(grep -m1 -c avx512vl /proc/cpuinfo)" != "1" ]; then
        if [ "${backend}" == "cpp512y" ]; then echo; echo "(SKIP ${backend} which is not supported - no avx512vl in /proc/cpuinfo)"; continue; fi
        if [ "${backend}" == "cpp512z" ]; then echo; echo "(SKIP ${backend} which is not supported - no avx512vl in /proc/cpuinfo)"; continue; fi
      fi
      if ls -d build.${backend}* > /dev/null 2>&1; then
        bdirs="$(ls -d build.${backend}*)"
        for bdir in ${bdirs}; do
          runExe ${bdir}/runTest.exe
          if [ -f ${bdir}/check.exe ]; then
            runExe ${bdir}/check.exe -p 1 32 1
          elif [ -f ${bdir}/gcheck.exe ]; then
            runExe ${bdir}/gcheck.exe -p 1 32 1
          else
            echo "ERROR! Neither ${bdir}/check.exe nor ${bdir}/gcheck.exe was found?"; exit 1
          fi
        done
      fi
    done
    popd >& /dev/null
  done
}

#----------------------------------------------------------------------------------------------------------------------------------

# Usage
function usage() {
  echo "Usage: $(basename $0) <${stages// /|}> <proc.sa|proc.mad>"
  exit 1
}

#----------------------------------------------------------------------------------------------------------------------------------

# Valid stages
stages="codegen before_build build after_build tput_test"

# Check input arguments
for astage in $stages; do
  if [ "$1" == "$astage" ]; then
    stage=$1; proc=$2; shift; shift; break
  fi
done
if [ "$stage" == "" ] || [ "$proc" == "" ] || [ "$1" != "" ]; then usage; fi

# Start
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "[testsuite_oneprocess.sh] $stage ($proc) starting at $(date)"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

# Execute stage
( set -e; $stage $proc) # execute this within a subprocess and fail immediately on error
status=$?

# Finish
echo
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
if [ $status -eq 0 ]; then
  echo "[testsuite_oneprocess.sh] $stage ($proc) finished with status=$status (OK) at $(date)"
else
  echo "[testsuite_oneprocess.sh] $stage ($proc) finished with status=$status (NOT OK) at $(date)"
fi
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
exit $status
