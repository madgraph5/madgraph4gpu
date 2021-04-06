#!/bin/bash

exes=
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.none/check.exe"
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.none/gcheck.exe"
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.avx512/check.exe"
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.avx512/gcheck.exe"
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
  unset OMP_NUM_THREADS
  echo "-------------------------------------------------------------------------"
  # For TIMEFORMAT see https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
  TIMEFORMAT=$'real\t%3lR' && time $exe -p 2048 256 12 2>&1 | egrep '(Process|fptype_sv|OMP threads|EvtsPerSec\[Matrix|MeanMatrix|TOTAL       :)'
  if [ "${exe%%/check*}" != "${exe}" ]; then 
    export OMP_NUM_THREADS=$(nproc --all)
    echo "-------------------------------------------------------------------------"
    TIMEFORMAT=$'real\t%3lR' && time $exe -p 2048 256 12 2>&1 | egrep '(Process|fptype_sv|OMP threads|EvtsPerSec\[Matrix|MeanMatrix|TOTAL       :)'
  fi
done
echo "-------------------------------------------------------------------------"
