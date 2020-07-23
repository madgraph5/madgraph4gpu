#!/bin/bash

if [ "x${1}" = "x" ] || [ "x{$2}" = "x" ]; then 
  echo error
  exit 1;
fi
 
blocks=$1
threads=$2

iters=20
if [ "x${3}" != "x" ]; then iters=$3; fi

skips=10
if [ "x${4}" != "x" ]; then skips=$4; fi

ncu --set full -f -s ${skips} -o profile_`date +%Y%m%d_%H%M`_${blocks}_${threads}_${iters}_s${skips} ./check.exe ${blocks} ${threads} ${iters}
