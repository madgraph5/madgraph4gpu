#!/bin/bash

#
#   __  __               _    ____                          _       _  _      ____   ____    _   _ 
#  |  \/  |   __ _    __| |  / ___|  _ __    __ _   _ __   | |__   | || |    / ___| |  _ \  | | | |
#  | |\/| |  / _` |  / _` | | |  _  | '__|  / _` | | '_ \  | '_ \  | || |_  | |  _  | |_) | | | | |
#  | |  | | | (_| | | (_| | | |_| | | |    | (_| | | |_) | | | | | |__   _| | |_| | |  __/  | |_| |
#  |_|  |_|  \__,_|  \__,_|  \____| |_|     \__,_| | .__/  |_| |_|    |_|    \____| |_|      \___/ 
#                                                  |_|                                             
#
#
#   Bash script for compiling and executing physics processes using the MadGraph5_aMC@NLO GPU development framework
#   using CUDA/HIP
#
#   Author: Jorgen Teig, CERN 2023
#

helpFunction()
{
    echo ""
    echo "Usage: $0 -n gg_ttgg -b 1024 -t 128 -i 10"
    echo -e "\t-n Name of the physics process being built and run"
    echo -e "\t-b Blocks per grid"
    echo -e "\t-t Threads per block"
    echo -e "\t-i Iterations"
    echo -e "\t-r Branch"
    echo -e "\t-c CUDA or HIP compilation"
    echo -e "\t-m Makefile arguments"
    exit 1 # Exit script after printing help
}

while getopts "n:b:t:i:r:m:a:" opt
do
    case "$opt" in
        n ) MG_PROC="$OPTARG" ;; #process to target
        b ) blocksPerGrid="$OPTARG" ;;
        t ) threadsPerBlock="$OPTARG" ;;
        i ) iterations="$OPTARG" ;;
        r ) branch="$OPTARG" ;;
        a ) gpuCompiler="$OPTARG" ;;
        m ) makeArgs="$OPTARG" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
    esac
done

# Print helpFunction in case parameters are empty
if [ -z "${MG_PROC}" ] || [ -z "${blocksPerGrid}" ] || [ -z "${threadsPerBlock}" ] || [ -z "${iterations}" ]
then
    echo "Some or all of the parameters are empty";
    helpFunction
fi

# Begin script in case all parameters are correct

##################################################################

# Set variables for later use

# CUDA
if [[ -z "${gpuCompiler}" ]] || [[ "${gpuCompiler,,}" == "cuda" ]]; then
    if [[ -z "$CUDA_HOME" ]]; then
        # Check if CUDA_HOME has not been set from the outside, usefull in CI/CD
        export CUDA_HOME="`which nvcc 2>/dev/null`" && while [ -L "$compiler" ]; do compiler=`readlink "$compiler"`; done && echo "$$compiler"
    fi
    # Sets CUDA in PATH
    export PATH=$CUDA_HOME:$PATH

# HIP
else if [[ "${gpuCompiler,,}" == "hip" ]]
    if [[ -z "$HIP_HOME" ]]; then
        # Check if HIP_HOME has not been set from the outside, usefull in CI/CD
        export HIP_HOME="`which hipcc 2>/dev/null`" && while [ -L "$compiler" ]; do compiler=`readlink "$compiler"`; done && echo "$$compiler"
    fi
    # Sets HIP to PATH
    export PATH=$HIP_HOME:$PATH
fi

# Prefix for saving the JSON files in workspace folder in the tools/profiling directory
prefix=$(dirname "$0")

export USEBUILDDIR=1
export NTPBMAX=1024
export CXX=`which g++`
export FC=`which gfortran`

export MG_EXE="./gcheck.exe" #GPU
#export MG_EXE="./check.exe" #CPU

export WORKSPACE=$prefix/workspace_mg4gpu

REPORT_FOLDER="${WORKSPACE}/$(date +"%y-%m-%d")_${CUDA_NAME_PREFIX}_${branch}"

mkdir $WORKSPACE 2>/dev/null; true
mkdir $REPORT_FOLDER 2>/dev/null; true

export MG_PROC_DIR=$prefix/../../epochX/cudacpp/$MG_PROC
export MG_SP_DIR=$MG_PROC_DIR/SubProcesses/P1_Sigma_sm_*
export MG5AMC_CARD_PATH=$MG_PROC_DIR/Cards

# Build executable

cd $MG_SP_DIR
make -j $makeArgs

# Run executable

cd build.${makeArgs:3}*
mkdir -p perf/data/ 2>/dev/null; true
$MG_EXE -j $blocksPerGrid $threadsPerBlock $iterations

echo "${MG_EXE} -j ${blocksPerGrid} ${threadsPerBlock} ${iterations}"

cd perf/data/
mv 0-perf-test-run0.json ${REPORT_FOLDER}/test_${MG_PROC}_${CUDA_NAME_PREFIX}_${blocksPerGrid}_${threadsPerBlock}_${iterations}.json