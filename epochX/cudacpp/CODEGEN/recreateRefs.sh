#!/bin/bash

cudacppdir=$(cd $(dirname $0)/..; pwd)
echo ${cudacppdir}

cd ${cudacppdir}
for d in $(git ls-tree --name-only HEAD *.mad/SubProcesses/P*); do
  echo; echo "================================================================================="
  dd=${cudacppdir}/$d
  cd $dd
  pwd
  ###make cleanall # overkill (executed from every P* subdirectory...)
  make -f cudacpp.mk -j BACKEND=cuda USEBUILDDIR=1
  CUDACPP_RUNTEST_DUMPEVENTS=1 ./build.cuda_d_inl0_hrd0/runTest_cuda.exe
  \cp ../../test/ref/dump* ../../../CODEGEN/PLUGIN/CUDACPP_SA_OUTPUT/test/ref  
done
