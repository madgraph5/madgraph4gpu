#!/bin/bash

omp=0
if [ "$1" == "-omp" ]; then
  omp=1
  shift
fi

if [ "$1" != "" ]; then
  echo "Usage: $0 [-omp]"
  exit 1
fi

exes=
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.none/check.exe"
###exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.none/gcheck.exe"
###exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.none/hcheck.exe"
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.avx512/check.exe"
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.avx512/gcheck.exe"
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.avx512/hcheck.exe"
exes="$exes ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/check.exe"
exes="$exes ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/gcheck.exe"

pushd ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum >& /dev/null
  pwd
  make AVX=avx512
  make AVX=none
popd >& /dev/null

pushd ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum >& /dev/null
  pwd
  make
popd >& /dev/null

for exe in $exes; do
  echo "-------------------------------------------------------------------------"
  unset OMP_NUM_THREADS
  # For TIMEFORMAT see https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
  TIMEFORMAT=$'real\t%3lR' && time $exe -p 2048 256 12 2>&1 | egrep '(Process|fptype_sv|OMP threads|EvtsPerSec\[Matrix|MeanMatrix|TOTAL       :)'
  if [ "${omp}" == "1" ] && [ "${exe%%/check*}" != "${exe}" ]; then 
    echo "-------------------------------------------------------------------------"
    export OMP_NUM_THREADS=$(nproc --all)
    TIMEFORMAT=$'real\t%3lR' && time $exe -p 2048 256 12 2>&1 | egrep '(Process|fptype_sv|OMP threads|EvtsPerSec\[Matrix|MeanMatrix|TOTAL       :)'
  elif [ "${exe%%/gcheck*}" != "${exe}" ]; then 
    sudo LD_LIBRARY_PATH=${LD_LIBRARY_PATH} $(which ncu) --metrics launch__registers_per_thread --target-processes all -k "sigmaKinE" --kernel-regex-base mangled --print-kernel-base mangled $exe -p 2048 256 1 | egrep '(sigmaKin|registers)' | tr "\n" " " | awk '{print $1, $2, $3, $15, $17}'
  fi
done
echo "-------------------------------------------------------------------------"
