#!/bin/bash

omp=0
avxall=0
ep2=0
cpp=1

function usage()
{
  echo "Usage: $0 [-omp] [-avxall] [-ep2]"
  echo "Usage: $0 [-nocpp] [-ep2]"
  exit 1
}

while [ "$1" != "" ]; do
  if [ "$1" == "-omp" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options -omp and -nocpp are incompatible"; usage; fi
    omp=1
    shift
  elif [ "$1" == "-avxall" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options -avxall and -nocpp are incompatible"; usage; fi
    avxall=1
    shift
  elif [ "$1" == "-nocpp" ]; then
    if [ "${avxall}" == "1" ]; then echo "ERROR! Options -avxall and -nocpp are incompatible"; usage; fi
    if [ "${omp}" == "1" ]; then echo "ERROR! Options -avxall and -nocpp are incompatible"; usage; fi
    cpp=0
    shift
  elif [ "$1" == "-ep2" ]; then
    ep2=1
    shift
  else
    usage
  fi
done

exes=
exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.none/gcheck.exe"
if [ "${cpp}" == "1" ]; then 
  exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.none/check.exe"
fi
if [ "${avxall}" == "1" ]; then 
  exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.sse4/check.exe"
  exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.avx2/check.exe"
fi
if [ "${cpp}" == "1" ]; then 
  exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.512y/check.exe"
fi
if [ "${avxall}" == "1" ]; then 
  exes="$exes ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/build.512z/check.exe"
fi
if [ "${ep2}" == "1" ]; then 
  exes="$exes ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/gcheck.exe"
  if [ "${cpp}" == "1" ]; then 
    exes="$exes ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/check.exe"
  fi
fi

export USEBUILDDIR=1
pushd ../../../../../epoch1/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum >& /dev/null
pwd
make AVX=none
if [ "${avxall}" == "1" ]; then make AVX=sse4; fi
if [ "${avxall}" == "1" ]; then make AVX=avx2; fi
make AVX=512y # always consider 512y as the reference, even if for clang avx2 is slightly faster...
if [ "${avxall}" == "1" ]; then make AVX=512z; fi
popd >& /dev/null

if [ "${ep2}" == "1" ]; then 
  pushd ../../../../../epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum >& /dev/null
  pwd
  make
  popd >& /dev/null
fi

function runExe() {
  exe=$1
  ###echo "runExe $exe OMP=$OMP_NUM_THREADS"
  pattern="Process|fptype_sv|OMP threads|EvtsPerSec\[Matrix|MeanMatrix|FP precision|TOTAL       :"
  # Optionally add other patterns here for some specific configurations (e.g. clang)
  pattern="${pattern}|CUCOMPLEX"
  pattern="${pattern}|COMMON RANDOM"
  if perf --version >& /dev/null; then
    # -- Newer version using perf stat
    pattern="${pattern}|instructions|cycles"
    pattern="${pattern}|elapsed"
    perf stat $exe -p 2048 256 12 2>&1 | egrep "(${pattern})" | grep -v "Performance counter stats"
  else
    # -- Older version using time
    # For TIMEFORMAT see https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
    TIMEFORMAT=$'real\t%3lR' && time $exe -p 2048 256 12 2>&1 | egrep "(${pattern})"
  fi
}

function runNcu() {
  exe=$1
  ###echo "runExe $exe OMP=$OM_NUM_THREADS (NCU)"
  $(which ncu) --metrics launch__registers_per_thread --target-processes all --kernel-id "::sigmaKin:" --print-kernel-base mangled $exe -p 2048 256 1 | egrep '(sigmaKin|registers)' | tr "\n" " " | awk '{print $1, $2, $3, $15, $17}'
}

echo -e "\nOn $HOSTNAME ($(nvidia-smi -L | awk '{print $5}')):"
for exe in $exes; do
  if [ ! -f $exe ]; then continue; fi
  echo "-------------------------------------------------------------------------"
  unset OMP_NUM_THREADS
  runExe $exe
  if [ "${exe%%/check*}" != "${exe}" ]; then 
    obj=${exe%%/check*}/CPPProcess.o; ./simdSymSummary.sh -stripdir ${obj}
    if [ "${omp}" == "1" ]; then 
      echo "-------------------------------------------------------------------------"
      export OMP_NUM_THREADS=$(nproc --all)
      runExe $exe
    fi
  elif [ "${exe%%/gcheck*}" != "${exe}" ]; then 
    runNcu $exe
  fi
done
echo "-------------------------------------------------------------------------"
