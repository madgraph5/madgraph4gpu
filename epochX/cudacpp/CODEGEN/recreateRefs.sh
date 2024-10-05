#!/bin/bash

cudacppdir=$(cd $(dirname $0)/..; pwd)
echo ${cudacppdir}
cd ${cudacppdir}

function usage()
{
  echo "Usage: $0 [<processdirectory>]"
  exit 1
}

procs=$(git ls-tree --name-only HEAD *.mad)
if [ "$2" != "" ]; then
  usage
elif [ "$1" != "" ]; then
  if [ -d $1 ]; then
    procs=$1
  else
    echo "ERROR! Directory not found: $1"
    usage
  fi
fi

for p in ${procs}; do
  for d in $(ls -d ${p}/SubProcesses/P*); do
    echo; echo "================================================================================="
    dd=${cudacppdir}/$d
    pushd $dd >& /dev/null
    pwd
    ###make cleanall # overkill (executed from every P* subdirectory...)
    make -f cudacpp.mk -j BACKEND=cuda USEBUILDDIR=1
    CUDACPP_RUNTEST_DUMPEVENTS=1 ./build.cuda_d_inl0_hrd0/runTest_cuda.exe
    \cp ../../test/ref/dump* ../../../CODEGEN/PLUGIN/CUDACPP_SA_OUTPUT/test/ref  
    popd >& /dev/null
  done
done
