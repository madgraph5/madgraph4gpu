#!/bin/bash
host=$(hostname)
if [ "${host%%raplab*}" != "${host}" ]; then
  logs=logs_raplab
elif [ "${host%%cern.ch}" != "${host}" ] && [ "${host##b}" != "${host}" ]; then
  logs=logs_lxbatch
else
  logs=logs
fi
cmd="./gcheck.exe -p 2048 256 12"
trace=$logs/eemumuAV_`date +%m%d_%H%M`
( time ${cmd} ) 2>&1 | tee ${trace}.txt
if [ "${host%%cern.ch}" != "${host}" ] && [ "${host##b}" != "${host}" ]; then
  /usr/local/cuda-11.0/bin/ncu --set full -o ${trace} ${cmd}
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
  ncu --set full -o ${trace} ${cmd}
  nsys profile -o ${trace} ${cmd}
  echo ""
  echo "TO ANALYSE TRACE FILES:"
  echo "  ncu-ui &"
  echo "  nsight-sys &"
fi

