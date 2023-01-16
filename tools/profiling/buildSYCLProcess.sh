#!/bin/bash

helpFunction()
{
    echo ""
    echo "Usage: $0 -n gg_ttgg -b 1024 -t 128 -i 10"
    echo -e "\t-n Name of the physics process being built and run"
    echo -e "\t-b Blocks per grid"
    echo -e "\t-t Threads per block"
    echo -e "\t-i Iterations"
    echo -e "\t-r Branch"
    exit 1 # Exit script after printing help
}

while getopts "n:b:t:i:r:" opt
do
    case "$opt" in
        n ) MG_PROC="$OPTARG" ;; #process to target
        b ) blocksPerGrid="$OPTARG" ;;
        t ) threadsPerBlock="$OPTARG" ;;
        i ) iterations="$OPTARG" ;;
        r ) branch="$OPTARG" ;;
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

# Assign correct SM level for NVIDIA GPUs

# Check if nvidia-smi command exists
if command -v nvidia-smi > /dev/null 2>&1; then

    # Get the name of the GPU
    GPU_NAME=$(lshw -C display | grep -i "product:" | awk -F'[][]' '{print $2}')
else
    echo "nvidia-smi non existent on system, Nvidia GPU not present!"
    exit

case $GPU_NAME in
    Tesla V100S PCIe 32GB )
        export SM_LEVEL="sm_70"

        # GPU
        export DEVICE_ID=0
        # CPU
        #export DEVICE_ID=1
        ;;
    A100 PCIe 40GB )
        export SM_LEVEL="sm_80"

        # GPU
        export DEVICE_ID=2
        # CPU
        #export DEVICE_ID=1
        ;;
esac

##################################################################

# Set variables for later use

# Assumes that this is run from profiling directory in the repo
prefix=$(pwd)

#export DPCPP_HOME=/p/project/prpb109/sycl_workspace
export DPCPP_HOME=/afs/cern.ch/work/j/jteig/sycl_workspace
export USEBUILDDIR=1
export NTPBMAX=1024
export CXX=$DPCPP_HOME/llvm/build/bin/clang++
export CUDA_PATH=/usr/local/cuda-11.8/
export SYCLFLAGS="-fsycl -fsycl-targets=nvptx64-nvidia-cuda -Xsycl-target-backend '--cuda-gpu-arch=$SM_LEVEL' -fgpu-rdc --cuda-path=$CUDA_PATH"
export WORKSPACE=$prefix/workspace_mg4gpu
#export NAME_PREFIX="sycl_v100s_cuda_11.6.2_gcc_11.3"
#export NAME_PREFIX="sycl_Xeon-Silver-4216_a100s_cuda-11.6.2_gcc-11.3"

# Branch should be enviroment variable in main script and then passed down if none then it is not displayed in prefix
REPORT_FOLDER="${WORKSPACE}/$(date +"%y-%m-%d")_${SYCL_NAME_PREFIX}_${branch}"

# If unknown set at the run step after running LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MG_LIBS $MG_EXE --param_card $MG5AMC_CARD_PATH/param_card.dat --device_info 1024 128 10
# GPU
export DEVICE_ID=0
# CPU
#export DEVICE_ID=1

# Finds correct subprocess
case $MG_PROC in
    ee_mumu ) export MG_SUBPROC="P1_Sigma_sm_epem_mupmum" ;;
    ee_mumu.sa ) export MG_SUBPROC="P1_Sigma_sm_epem_mupmum" ;;
    gg_tt ) export MG_SUBPROC="P1_Sigma_sm_gg_ttx" ;;
    gg_tt.sa ) export MG_SUBPROC="P1_Sigma_sm_gg_ttx" ;;
    gg_ttg ) export MG_SUBPROC="P1_Sigma_sm_gg_ttxg" ;;
    gg_ttg.sa ) export MG_SUBPROC="P1_Sigma_sm_gg_ttxg" ;;
    gg_ttgg ) export MG_SUBPROC="P1_Sigma_sm_gg_ttxgg" ;;
    gg_ttgg.sa ) export MG_SUBPROC="P1_Sigma_sm_gg_ttxgg" ;;
    gg_ttggg ) export MG_SUBPROC="P1_Sigma_sm_gg_ttxggg" ;;
    gg_ttggg.sa ) export MG_SUBPROC="P1_Sigma_sm_gg_ttxggg" ;;
esac

mkdir -p $WORKSPACE/mg4gpu/lib 2>/dev/null; true
mkdir -p $WORKSPACE/mg4gpu/bin 2>/dev/null; true
mkdir $REPORT_FOLDER 2>/dev/null; true

export MG4GPU_LIB=$WORKSPACE/mg4gpu/lib
export MG4GPU_BIN=$WORKSPACE/mg4gpu/bin

export MG_PROC_DIR=$prefix/../../epochX/sycl/$MG_PROC
export MG_SP_DIR=$MG_PROC_DIR/SubProcesses/$MG_SUBPROC

export MG_LIBS_DIR="${MG4GPU_LIB}/build_${MG_PROC}_${SYCL_NAME_PREFIX}"
export MG_LIBS="$DPCPP_HOME/llvm/build/lib:$MG_LIBS_DIR"

export MG_EXE_DIR="${MG4GPU_BIN}/build_${MG_PROC}_${SYCL_NAME_PREFIX}"
export MG_EXE="$MG_EXE_DIR/check.exe"
export MG5AMC_CARD_PATH=$MG_PROC_DIR/Cards

# Build executable
cd $MG_SP_DIR
make -j
mv ../../lib/build.d_inl0*/ $MG_LIBS_DIR #2>/dev/null; true
mv build.d_inl0*/ $MG_EXE_DIR #2>/dev/null; true

# Run executable
cd $WORKSPACE

# Display the devices
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MG_LIBS $MG_EXE --param_card $MG5AMC_CARD_PATH/param_card.dat --device_info 32 32 10

# Add MG Libs to linker library path and run the executable
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MG_LIBS $MG_EXE -j --json_file ${REPORT_FOLDER}/test_${MG_PROC}_${SYCL_NAME_PREFIX}_${blocksPerGrid}_${threadsPerBlock}_${iterations}.json --param_card $MG5AMC_CARD_PATH/param_card.dat --device_id $DEVICE_ID $blocksPerGrid $threadsPerBlock $iterations

# View output
#nano $REPORT_FOLDER/test_${SYCL_NAME_PREFIX}_${MG_PROC}_${blocksPerGrid}_${threadsPerBlock}_${iterations}.json-+