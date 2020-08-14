#!/bin/bash

usage(){
  echo "Usage: $0 [--cc] [-p #blocks #threads #iterations] [-l label]"
  exit 1
}

# Default options
tag=cu
cuargs="16384 32 12" # NEW DEFAULT 20.08.10 (faster on local, and allows comparison to global and shared memory)
ccargs="  256 32 12" # Similar to cuda config, but faster than using "16384 32 12"
args=
label=

# Command line arguments
while [ "$1" != "" ]; do
  # Profile C++ instead of cuda
  if [ "$1" == "--cc" ]; then
    tag=cc
  # Override blocks/threads/iterations
  # (NB do not exceed 12 iterations: profiling overhead per iteration is huge)
  elif [ "$1" == "-p" ]; then
    if [ "$4" != "" ]; then
      args="$2 $3 $4"    
      shift 4
    else
      usage
    fi
  # Label
  elif [ "$1" == "-l" ]; then
    if [ "$2" != "" ]; then
      label="$2"
      shift 2
    else
      usage
    fi
  # Invalid arguments
  else
    usage
  fi
done

if [ "$label" == "" ]; then
  echo "ERROR! You must specify a label"
  usage
fi

if [ "$tag" == "cu" ]; then
  if [ "$args" == "" ]; then args=$cuargs; fi
  cmd="./gcheck.exe -p $args"
else
  if [ "$args" == "" ]; then args=$ccargs; fi
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
if [ "$label" != "" ]; then trace=${trace}_${label}; fi

echo
echo "PROFILE output: ${trace}.*"
echo

\rm -f ${trace}.*

hostname > ${trace}.txt
echo "nproc=$(nproc)" >> ${trace}.txt
echo >> ${trace}.txt
( time ${cmd} ) 2>&1 | tee -a ${trace}.txt
nvidia-smi -q -d CLOCK >> ${trace}.txt

if [ "${host%%cern.ch}" != "${host}" ] && [ "${host##b}" != "${host}" ]; then
  if [ "$tag" == "cu" ]; then
    /usr/local/cuda-11.0/bin/ncu --set full -o ${trace} ${cmd}
  fi
  ###/usr/local/cuda-10.1/bin/nsys profile -o ${trace} ${cmd}
  ###/usr/local/cuda-10.2/bin/nsys profile -o ${trace} ${cmd}
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

