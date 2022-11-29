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

# Set variables for later use

# Assumes that this is run from profiling directory in the repo
prefix=$(pwd)

export USEBUILDDIR=1
export NTPBMAX=1024
export CXX=/cvmfs/sft.cern.ch/lcg/releases/gcc/11.3.0-ad0f5/x86_64-centos8/bin/g++
export MG_EXE="./gcheck.exe" #GPU
#export MG_EXE="./check.exe" #CPU
export CUDA_HOME=/usr/local/cuda-11.6/
export FC=`which gfortran`
export WORKSPACE=$prefix/workspace_mg4gpu
#export NAME_PREFIX="cudacpp_v100s_cuda_11.6.2_gcc_11.3"

REPORT_FOLDER="${WORKSPACE}/$(date +"%y-%m-%d")_${CUDA_NAME_PREFIX}_${branch}"

# Sets CUDA in PATH
export PATH=$CUDA_HOME:$PATH

mkdir $WORKSPACE 2>/dev/null; true
mkdir $REPORT_FOLDER 2>/dev/null; true

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

export MG_PROC_DIR=$prefix/../../epochX/cudacpp/$MG_PROC
export MG_SP_DIR=$MG_PROC_DIR/SubProcesses/$MG_SUBPROC
export MG5AMC_CARD_PATH=$MG_PROC_DIR/Cards

# Build executable

cd $MG_SP_DIR
make

# Run executable

cd build*
mkdir -p perf/data/ 2>/dev/null; true
$MG_EXE -j $blocksPerGrid $threadsPerBlock $iterations

cd perf/data/
mv 0-perf-test-run0.json ${REPORT_FOLDER}/test_${MG_PROC}_${CUDA_NAME_PREFIX}_${blocksPerGrid}_${threadsPerBlock}_${iterations}.json