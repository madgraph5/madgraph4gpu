#!/bin/bash

omp=1 # new default: OMP only for epoch1
if [ "$1" == "-omp" ]; then
  omp=2
  shift
elif [ "$1" == "-noomp" ]; then
  omp=0
  shift
fi

avxall=0
if [ "$1" == "-avxall" ]; then
  avxall=1
  shift
fi

if [ "$1" != "" ]; then
  echo "Usage: $0 [-omp|-noomp] [-avxall]"
  exit 1
fi

exes=
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.none/check.exe"
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.none/gcheck.exe"
###exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.none/hcheck.exe"
if [ "${avxall}" == "1" ]; then 
  exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.sse4/check.exe"
  exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.avx2/check.exe"
fi
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.512y/check.exe"
###exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.512y/gcheck.exe"
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.512y/hcheck.exe"
if [ "${avxall}" == "1" ]; then 
  exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.512z/check.exe"
fi
exes="$exes ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/check.exe"
exes="$exes ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/gcheck.exe"

export USEBUILDDIR=1
pushd ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum >& /dev/null
  pwd
  make AVX=none
  if [ "${avxall}" == "1" ]; then make AVX=sse4; fi
  if [ "${avxall}" == "1" ]; then make AVX=avx2; fi
  make AVX=512y
  if [ "${avxall}" == "1" ]; then make AVX=512z; fi
popd >& /dev/null

pushd ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum >& /dev/null
  pwd
  make
popd >& /dev/null

function runExe() {
  exe=$1
  ###echo "runExe $exe OMP=$OMP_NUM_THREADS"
  # For TIMEFORMAT see https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
  if [ "${exe%%/hcheck*}" != "${exe}" ]; then 
    TIMEFORMAT=$'real\t%3lR' && time $exe -p 2048 256 12 2>&1 | grep -v NaN | egrep '(Process|fptype_sv|OMP threads|EvtsPerSec\[Matrix|MeanMatrix|TOTAL       :|TotalEventsComputed)' | sort -k"2.1,2.6" -r | uniq 
  else
    TIMEFORMAT=$'real\t%3lR' && time $exe -p 2048 256 12 2>&1 | egrep '(Process|fptype_sv|OMP threads|EvtsPerSec\[Matrix|MeanMatrix|TOTAL       :)'
  fi
}

function runNcu() {
  exe=$1
  ###echo "runExe $exe OMP=$OM_NUM_THREADS (NCU)"
  ###sudo LD_LIBRARY_PATH=${LD_LIBRARY_PATH} $(which ncu) --metrics launch__registers_per_thread --target-processes all -k "sigmaKinE" --kernel-regex-base mangled --print-kernel-base mangled $exe -p 2048 256 1 | egrep '(sigmaKin|registers)' | tr "\n" " " | awk '{print $1, $2, $3, $15, $17}'
  $(which ncu) --metrics launch__registers_per_thread --target-processes all -k "sigmaKinE" --kernel-regex-base mangled --print-kernel-base mangled $exe -p 2048 256 1 | egrep '(sigmaKin|registers)' | tr "\n" " " | awk '{print $1, $2, $3, $15, $17}'
}

for exe in $exes; do
  echo "-------------------------------------------------------------------------"
  unset OMP_NUM_THREADS
  if [ "${exe%%/hcheck*}" != "${exe}" ]; then 
    export OMP_NUM_THREADS=$(nproc --all)
  fi
  runExe $exe
  if [ "${exe%%/check*}" != "${exe}" ]; then 
    obj=${exe%%/check*}/CPPProcess.o; ./simdSymSummary.sh -stripdir ${obj}
    if [ "${omp}" != "0" ] && [ "${exe%%/epoch1*}" != "${exe}" ]; then 
      echo "-------------------------------------------------------------------------"
      export OMP_NUM_THREADS=$(nproc --all)
      runExe $exe
    elif [ "${omp}" == "2" ] && [ "${exe%%/epoch2*}" != "${exe}" ]; then 
      echo "-------------------------------------------------------------------------"
      export OMP_NUM_THREADS=$(nproc --all)
      runExe $exe
    fi
  elif [ "${exe%%/gcheck*}" != "${exe}" ]; then 
    runNcu $exe
  fi
done
echo "-------------------------------------------------------------------------"
