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
#   using oneAPI/SYCL
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
    echo -e "\t-d Flag for setting device id"
    exit 1 # Exit script after printing help
}

while getopts "n:b:t:i:r:d:" opt
do
    case "$opt" in
        n ) MG_PROC="$OPTARG" ;; #process to target
        b ) blocksPerGrid="$OPTARG" ;;
        t ) threadsPerBlock="$OPTARG" ;;
        i ) iterations="$OPTARG" ;;
        r ) branch="$OPTARG" ;;
        d ) DEVICE_ID="$OPTARG" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
    esac
done

# Print helpFunction in case parameters are empty
if [ -z "${MG_PROC}" ] || [ -z "${blocksPerGrid}" ] || [ -z "${threadsPerBlock}" ] || [ -z "${iterations}" ]
then
    echo "Some or all of the parameters are empty";
    helpFunction
fi

# Added check if the SYCL_NAME_PREFIX variable are not set
if [ -z "$SYCL_NAME_PREFIX" ]; then
    echo "WARNING: SYCL_NAME_PREFIX is not set. Cannot append system info to JSON file names!"
fi

##################################################################

# Assign correct SM level for NVIDIA GPUs

# Check if nvidia-smi command exists
if command -v nvidia-smi > /dev/null 2>&1; then

    # Get the name of the GPU
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader)

    # GPU (DEVICE_ID=2 for oneAPI toolkit runs on GPUs, else DEVICE_ID=0 with LLVM compiler)
    export DEVICE_ID=2
    # CPU
    #export DEVICE_ID=1
else
    echo "nvidia-smi non existent on system, Nvidia GPU possibly not present!"
    exit
fi

case $GPU_NAME in
    *V100S* ) export SM_LEVEL="sm_70" ;;
    *A100* ) export SM_LEVEL="sm_80" ;;
esac

##################################################################

# Begin script in case all parameters and GPU specific settings are set

##################################################################

# Set variables for later use

# Assumes that this is run from profiling directory in the repo
prefix=$(pwd)

export USEBUILDDIR=1
export NTPBMAX=1024
export CUDA_PATH=/usr/local/cuda-12.0/
export WORKSPACE=$prefix/workspace_mg4gpu

export CXTYPE="thrust"

# Old SYCLFLAGS
# export SYCLFLAGS="-fsycl -fsycl-targets=nvptx64-nvidia-cuda -Xsycl-target-backend '--cuda-gpu-arch=$SM_LEVEL' -fgpu-rdc --cuda-path=$CUDA_PATH"

export SYCLFLAGS="-fsycl -fsycl-targets=nvptx64-nvidia-cuda -Xcuda-ptxas --maxrregcount=255 -Xcuda-ptxas --verbose -Xsycl-target-backend --cuda-gpu-arch=$SM_LEVEL"

# Compilation using OneAPI Toolkit through CVMFS
#export CXX=/cvmfs/projects.cern.ch/intelsw/oneAPI/linux/x86_64/2023/compiler/2023.0.0/linux/bin-llvm/clang++

# Compilation with LLVM DPC++ compiler
export DPCPP_HOME=/afs/cern.ch/work/j/jteig/sycl_workspace
export CXX=$DPCPP_HOME/llvm/llvm-20230418-fea99cc9ad67-gcc-11.2.1-cuda-12.0/bin/clang++

# Sets CUDA in PATH
export PATH=$CUDA_HOME:$PATH

# Branch should be enviroment variable in main script and then passed down if none then it is not displayed in prefix
REPORT_FOLDER="${WORKSPACE}/$(date +"%y-%m-%d")_${SYCL_NAME_PREFIX}_${branch}"

mkdir -p $WORKSPACE/mg4gpu/lib 2>/dev/null; true
mkdir -p $WORKSPACE/mg4gpu/bin 2>/dev/null; true
mkdir $REPORT_FOLDER 2>/dev/null; true

export MG4GPU_LIB=$WORKSPACE/mg4gpu/lib
export MG4GPU_BIN=$WORKSPACE/mg4gpu/bin

export MG_PROC_DIR=$prefix/../../epochX/sycl/$MG_PROC
export MG_SP_DIR=$MG_PROC_DIR/SubProcesses/P1_*

export MG_LIBS_DIR="${MG4GPU_LIB}/build_${MG_PROC}_${SYCL_NAME_PREFIX}"

if [[ -z "${DPCPP_HOME}" ]]; then
    export MG_LIBS="$MG_LIBS_DIR"
else
    export MG_LIBS="$DPCPP_HOME/llvm/build/lib:$MG_LIBS_DIR"
fi

export MG_EXE_DIR="${MG4GPU_BIN}/build_${MG_PROC}_${SYCL_NAME_PREFIX}"
export MG_EXE="$MG_EXE_DIR/check.exe"
export MG5AMC_CARD_PATH=$MG_PROC_DIR/Cards

# Build executable
cd $MG_SP_DIR
make -j build.d_inl0_hrd1/check.exe
mv -f ../../lib/build.*/ $MG_LIBS_DIR #2>/dev/null; true
mv -f build.*/ $MG_EXE_DIR

# Run executable
cd $WORKSPACE

if [ $DEVICE_ID == "info" ]; then
    # Add MG Libs to linker library path and display the devices
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MG_LIBS $MG_EXE --param_card $MG5AMC_CARD_PATH/param_card.dat --device_info 32 32 10

else
    # Add MG Libs to linker library path and run the executable
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MG_LIBS $MG_EXE -j --json_file ${REPORT_FOLDER}/test_${MG_PROC}_${SYCL_NAME_PREFIX}_${blocksPerGrid}_${threadsPerBlock}_${iterations}.json --param_card $MG5AMC_CARD_PATH/param_card.dat --device_id $DEVICE_ID $blocksPerGrid $threadsPerBlock $iterations
fi
