#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Jun 2023) for the MG5aMC CUDACPP plugin.

set -e # fail on error

function usage()
{
  echo "Usage: $0 <-FORTRAN|-CUDA|-CPP>"
  exit 1
}

bckend=
while [ "$1" != "" ]; do
  if [ "$1" == "-FORTRAN" ] || [ "$1" == "-CUDA" ] || [ "$1" == "-CPP" ]; then
    if [ "$bckend" == "" ]; then bckend=${1/-/}; else echo "ERROR! Backend already set"; usage; fi
    shift
  else
    echo "ERROR! Invalid option '$1'"
    usage
  fi
done
if [ "$bckend" == "" ]; then echo "ERROR! No backend was specified"; usage; fi

proc=gg_tt.mad

cd $(dirname $0)/..
echo "Execute $(basename $0) for process ${proc} in directory $(pwd)"
procdir=$(pwd)/${proc}
cd $procdir

function lauX_cleanup()
{
  rm -f crossx.html index.html
  rm -f SubProcesses/results.dat
  rm -rf Events HTML; mkdir Events HTML; touch Events/.keep HTML/.keep
  for d in SubProcesses/P*; do cd $d; rm -rf gensym input_app.txt symfact.dat G[0-9]* ajob[0-9]*; cd -; done
}

# Cleanup before launch
lauX_cleanup
rm -f SubProcesses/ME5_debug
cp SubProcesses/randinit SubProcesses/randinit.BKP # save the initial randinit
cp Cards/run_card.dat Cards/run_card.dat.BKP # save the initial run_card.dat

# Set the backend in run_card.dat
sed -i "s/FORTRAN = cudacpp_backend/${bckend} = cudacpp_backend/" Cards/run_card.dat

# Launch (generate_events)
###set -x # verbose
MG5AMC_CARD_PATH=$(pwd)/Cards ./bin/generate_events -f # (BUG #683: this does not return an error code even if it fails)
###set +x # not verbose

# Cleanup after launch
lauX_cleanup
mv SubProcesses/randinit.BKP SubProcesses/randinit # restore the initial randinit
mv Cards/run_card.dat.BKP Cards/run_card.dat # restore the initial run_card.dat
