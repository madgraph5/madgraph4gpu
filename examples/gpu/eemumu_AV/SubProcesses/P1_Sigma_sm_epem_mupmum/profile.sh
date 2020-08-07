#!/bin/bash

usage(){
  echo "Usage: $0 [--cc] [-p #blocks #threads #iterations]"
  exit 1
}

# Profile cuda or cpp?
tag=cu
args="2048 256 12" # DEFAULT 20.08.06 AND BEFORE
###args="32 256 12" # FOR LOCAL/GLOBAL/SHARED COMPARISON 20.08.07
if [ "$1" == "--cc" ]; then
  tag=cc
  args="64 64 12"
  shift
fi

# Override blocks/threads/iterations
if [ "$1" == "-p" ]; then
  if [ "$4" != "" ]; then
    args="$2 $3 $4"    
    shift 4
  else
    usage
  fi
fi

if [ "$1" != "" ]; then usage; fi

if [ "$tag" == "cu" ]; then
  cmd="./gcheck.exe -p $args"
else
  cmd="./check.exe -p $args"
fi

arg1=$(echo $args | cut -d' ' -f1)
arg2=$(echo $args | cut -d' ' -f2)
arg3=$(echo $args | cut -d' ' -f3)

host=$(hostname)
if [ "${host%%raplab*}" != "${host}" ]; then
  logs=logs_raplab
elif [ "${host%%cern.ch}" != "${host}" ] && [ "${host##b}" != "${host}" ]; then
  logs=logs_lxbatch
else
  logs=logs
fi
trace=$logs/eemumuAV_${tag}_`date +%m%d_%H%M`_b${arg1}_t${arg2}_i${arg3}
( time ${cmd} ) 2>&1 | tee ${trace}.txt
if [ "${host%%cern.ch}" != "${host}" ] && [ "${host##b}" != "${host}" ]; then
  if [ "$tag" == "cu" ]; then
    /usr/local/cuda-11.0/bin/ncu --set full -o ${trace} ${cmd}
  fi
  ###/usr/local/cuda-10.1/bin/nsys profile -o ${trace} ${cmd}
  ###/usr/local/cuda-10.2/bin/nsys profile -o ${trace} ${cmd}
  ###/cvmfs/sft.cern.ch/lcg/releases/cuda/10.2-9d877/x86_64-centos7-gcc62-opt/bin/nsys profile -o ${trace} ${cmd}
  /cvmfs/sft.cern.ch/lcg/releases/cuda/11.0RC-d9c38/x86_64-centos7-gcc62-opt/bin/nsys profile -o ${trace} ${cmd}
  echo ""
  echo "TO ANALYSE TRACE FILES:"
  ###echo "  ncu-ui &"
  ###echo "  nsight-sys &"
  echo "  Launch the Nsight Compute or Nsight System GUI from Windows"
else
  if [ "$tag" == "cu" ]; then
    ncu --set full -o ${trace} ${cmd}
  fi
  nsys profile -o ${trace} ${cmd}
  echo ""
  echo "TO ANALYSE TRACE FILES:"
  echo "  ncu-ui &"
  echo "  nsight-sys &"
fi

