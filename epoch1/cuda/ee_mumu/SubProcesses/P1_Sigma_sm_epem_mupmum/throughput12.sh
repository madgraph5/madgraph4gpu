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
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/hcheck.exe"
exes="$exes ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/check.exe"
exes="$exes ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/gcheck.exe"

for exe in $exes; do
  unset OMP_NUM_THREADS
  # For TIMEFORMAT see https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
  TIMEFORMAT=$'real\t%3lR' && time $exe -p 2048 256 12 2>&1 | egrep -i '(process|OMP threads|EvtsPerSec\[Matrix|MeanMatrix|TOTAL       :)'
  echo "-------------------------------------------------------------------------"
  if [ "${omp}" == "1" ] && [ "${exe%%/check*}" != "${exe}" ]; then 
    export OMP_NUM_THREADS=$(nproc --all)
    TIMEFORMAT=$'real\t%3lR' && time $exe -p 2048 256 12 2>&1 | egrep -i '(process|OMP threads|EvtsPerSec\[Matrix|MeanMatrix|TOTAL       :)'
    echo "-------------------------------------------------------------------------"
  fi
done
