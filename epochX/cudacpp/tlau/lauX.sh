#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Jun 2023) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2023) for the MG5aMC CUDACPP plugin.

set -e # fail on error

scrdir=$(cd $(dirname $0); pwd)

function usage()
{
  echo "Usage:   $0 <-FORTRAN|-CUDA|-CPP> <proc>"
  echo "Example: $0 -CPP gg_tt"
  exit 1
}

bckend=
proc=
suff=.mad
while [ "$1" != "" ]; do
  if [ "$1" == "-FORTRAN" ] || [ "$1" == "-CUDA" ] || [ "$1" == "-CPP" ]; then
    if [ "$bckend" == "" ]; then bckend=${1/-/}; else echo "ERROR! Backend already set"; usage; fi
  elif [ "$proc" == "" ]; then
    proc=${1}
  else
    echo "ERROR! Invalid option '$1': process directory already set to '${proc}'"
    usage
  fi
  shift
done
if [ "$bckend" == "" ]; then echo "ERROR! No backend was specified"; usage; fi
if [ "$proc" == "" ]; then echo "ERROR! No process directory was specified"; usage; fi
if [ ! -d ${proc}${suff} ]; then echo "ERROR! Process directory '${proc}' does not exist"; usage; fi

cd $(dirname $0)/..
echo "Execute $(basename $0) for process ${proc} and backend ${bckend} in directory $(pwd)"
procdir=$(pwd)/${proc}${suff}
cd $procdir
resultsdir=${scrdir}/logs_${proc/_}_${bckend/}

function getnevt()
{
  if [ "${proc}" == "gg_tt" ]; then
    nevt=10000
  elif [ "${proc}" == "gg_ttg" ]; then
    nevt=10000
  elif [ "${proc}" == "gg_ttgg" ]; then
    nevt=1000
  elif [ "${proc}" == "gg_ttggg" ]; then
    nevt=100
  else
    echo "WARNING! Unknown process ${proc}" > /dev/stderr
    nevt=100
  fi
  echo $nevt
}

function lauX_makeclean()
{
  for d in SubProcesses/P*; do cd $d; make cleanall; cd -; break; done
}

function lauX_cleanup()
{
  rm -f crossx.html index.html
  rm -f SubProcesses/results.dat
  rm -rf Events HTML; mkdir Events HTML; touch Events/.keep HTML/.keep
  for d in SubProcesses/P*; do cd $d; rm -rf gensym input_app.txt symfact.dat G[0-9]* ajob[0-9]*; cd - > /dev/null; done
}

# Clean builds before launch
lauX_makeclean >& /dev/null

# Clean config before launch
rm -rf ${resultsdir}; mkdir ${resultsdir}
lauX_cleanup
rm -f SubProcesses/ME5_debug
echo "r=21" > SubProcesses/randinit # just in case a previous test was not cleaned up
cp SubProcesses/randinit SubProcesses/randinit.BKP # save the initial file
sed -i "s/.* = nevents/  10000 = nevents/" Cards/run_card.dat # just in case
sed -i "s/.* = cudacpp_backend/CPP = cudacpp_backend/" Cards/run_card.dat # just in case
cp Cards/run_card.dat Cards/run_card.dat.BKP # save the initial file
sed -i "s/      NEVENTS = .*/      NEVENTS = 10000/" Source/run_card.inc # just in case
cp Source/run_card.inc Source/run_card.inc.BKP # save the initial file
sed -i "s/8192 1 1/%(event)s         %(maxiter)s           %(miniter)s/" bin/internal/gen_ximprove.py # just in case
cp bin/internal/gen_ximprove.py bin/internal/gen_ximprove.py.BKP # save the initial file
sed -i "s/'int', 8192,'Number of points/'int', 1000,'Number of points/" bin/internal/madevent_interface.py # just in case
sed -i "s/'int', 1, 'Number of iterations'/'int', 5, 'Number of iterations'/" bin/internal/madevent_interface.py # just in case
cp bin/internal/madevent_interface.py bin/internal/madevent_interface.py.BKP # save the initial file

# Set the number of events and iterations in the survey step
sed -i "s/'int', 1000,'Number of points/'int', 8192,'Number of points/" bin/internal/madevent_interface.py
sed -i "s/'int', 5, 'Number of iterations'/'int', 1, 'Number of iterations'/" bin/internal/madevent_interface.py

# Set the number of events and iterations in the refine step
sed -i "s/%(event)s         %(maxiter)s           %(miniter)s/8192 1 1/" bin/internal/gen_ximprove.py

# Set the number of unweighted events in run_card.dat
nevt=$(getnevt)
sed -i "s/ 10000 = nevents/ ${nevt} = nevents/" Cards/run_card.dat

# Set the backend in run_card.dat
sed -i "s/CPP = cudacpp_backend/${bckend} = cudacpp_backend/" Cards/run_card.dat

# Launch (generate_events)
# (BUG #683: generate_events does not return an error code even if it fails)
###set -x # verbose
START=$(date +%s)
echo "START: $(date)"  |& tee ${resultsdir}/output.txt
MG5AMC_CARD_PATH=$(pwd)/Cards time ./bin/generate_events -f |& tee -a ${resultsdir}/output.txt
echo "END: $(date)" |& tee -a ${resultsdir}/output.txt
END=$(date +%s)
echo "ELAPSED: $((END-START)) seconds" |& tee -a ${resultsdir}/output.txt
###set +x # not verbose

# Process and keep results
\rm HTML/results.pkl
mv Events ${resultsdir}; mv HTML ${resultsdir}
gunzip ${resultsdir}/Events/run_01/unweighted_events.lhe.gz

# Clean config after launch
lauX_cleanup
mv SubProcesses/randinit.BKP SubProcesses/randinit # restore the initial file
mv Cards/run_card.dat.BKP Cards/run_card.dat # restore the initial file
mv Source/run_card.inc.BKP Source/run_card.inc # restore the initial file
mv bin/internal/gen_ximprove.py.BKP bin/internal/gen_ximprove.py # restore the initial file
mv bin/internal/madevent_interface.py.BKP bin/internal/madevent_interface.py # restore the initial file

# Add an 80-character separator
echo ""
echo "********************************************************************************"
echo ""
