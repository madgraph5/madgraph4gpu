#!/bin/bash
host=$(hostname)
if [ "${host%%raplab*}" != "${host}" ]; then
  logs=logs_raplab
else
  logs=logs
fi
cmd="./check.exe -p 2048 256 12"
trace=$logs/eemumuAV_`date +%m%d_%H%M`
time ${cmd} | tee ${trace}.txt
if [ "${host%%raplab*}" != "${host}" ]; then
  ncu -o ${trace} ${cmd}
  nsys profile -o ${trace} ${cmd}
else
  /usr/local/cuda-11.0/bin/ncu -o ${trace} ${cmd}
  /usr/local/cuda-10.1/bin/nsys profile -o ${trace} ${cmd}
  echo ""
  echo "TO ANALYSE TRACE FILES:"
  echo "  /usr/local/cuda-11.0/bin/ncu-ui &"
  echo "  /usr/local/cuda-10.1/bin/nsight-sys &"
fi

