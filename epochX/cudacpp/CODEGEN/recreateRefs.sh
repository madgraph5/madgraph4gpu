#!/bin/bash

cudacppdir=$(cd $(dirname $0)/..; pwd)
echo ${cudacppdir}

cd ${cudacppdir}
#for d in $(git ls-tree --name-only HEAD *.mad/SubProcesses/P* | grep 'P1_gg_tt'); do
for d in gg_tt.mad/SubProcesses/P1_gg_ttx gg_ttg.mad/SubProcesses/P1_gg_ttxg; do
  echo; echo "================================================================================="
  dd=${cudacppdir}/$d
  cd $dd
  pwd
  ###make cleanall # overkill (executed from every P* subdirectory...)
  make -f cudacpp.mk -j BACKEND=cuda
  CUDACPP_RUNTEST_DUMPEVENTS=1 ./runTest_cuda.exe
  \cp ../../test/ref/dump* ../../../CODEGEN/PLUGIN/CUDACPP_SA_OUTPUT/test/ref  
done
