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
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/check.exe"
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/gcheck.exe"
exes="$exes ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/check.exe"
exes="$exes ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/gcheck.exe"

pattern="Process|OMP threads|EvtsPerSec\[Matrix|MeanMatrix|FP precision|TOTAL       :"
# Optionally add other patterns here for some specific configurations (e.g. clang)
pattern="(${pattern})"

for exe in $exes; do
  if [ ! -f $exe ]; then continue; fi
  echo "-------------------------------------------------------------------------"
  unset OMP_NUM_THREADS
  # For TIMEFORMAT see https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
  TIMEFORMAT=$'real\t%3lR' && time $exe -p 2048 256 12 2>&1 | egrep "$pattern"
  if [ "${omp}" == "1" ] && [ "${exe%%/check*}" != "${exe}" ]; then 
    echo "-------------------------------------------------------------------------"
    export OMP_NUM_THREADS=$(nproc --all)
    TIMEFORMAT=$'real\t%3lR' && time $exe -p 2048 256 12 2>&1 | egrep "$pattern"
  elif [ "${exe%%/gcheck*}" != "${exe}" ]; then 
    ###sudo LD_LIBRARY_PATH=${LD_LIBRARY_PATH} $(which ncu) --metrics launch__registers_per_thread --target-processes all -k "sigmaKinE" --kernel-regex-base mangled --print-kernel-base mangled $exe -p 2048 256 1 | egrep '(sigmaKin|registers)' | tr "\n" " " | awk '{print $1, $2, $3, $15, $17}'
    $(which ncu) --metrics launch__registers_per_thread --target-processes all --kernel-id "::sigmaKin:" --print-kernel-base mangled $exe -p 2048 256 1 | egrep '(sigmaKin|registers)' | tr "\n" " " | awk '{print $1, $2, $3, $15, $17}'
  fi
done
echo "-------------------------------------------------------------------------"
