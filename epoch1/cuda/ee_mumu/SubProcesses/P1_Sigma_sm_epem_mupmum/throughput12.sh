#!/bin/bash

exes=
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/check.exe"
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/gcheck.exe"
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/hcheck.exe"
exes="$exes ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/check.exe"
exes="$exes ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/gcheck.exe"

for exe in $exes; do
  # For TIMEFORMAT see https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
  TIMEFORMAT=$'real\t%3lR' && time $exe -p 2048 256 12 2>&1 | egrep -i '(process|OMP threads|EvtsPerSec\[Matrix|MeanMatrix|TOTAL       :)'
  echo "-------------------------------------------------------------------------"
done
