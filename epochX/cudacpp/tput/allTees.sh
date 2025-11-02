#!/bin/bash
# Copyright (C) 2020-2025 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Apr 2022) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2022-2025) for the MG5aMC CUDACPP plugin.

scrdir=$(cd $(dirname $0); pwd)

# By default, use the madevent+cudacpp version of code and tee scripts (use -sa to use the standalone version instead)
# By default, build and run all tests (use -makeonly to only build all tests)
opts=
suff=".mad"

# New (Oct 2025): systematically run 'make clean' (unless -nomakeclean is specified)
# As troughputX.sh is now using 'make -f cudacpp.mk', rebuilds are fast as they are fully available from ccache
makeclean=-makeclean

# By default, build and run all backends
# (AV private config: on itgold91 build and run only the C++ backends)
bblds=
if [ "$(hostname)" == "itgold91.cern.ch" ]; then bblds=-cpponly; fi

# Usage
function usage()
{
  echo "Usage (1): $0 [-e] [-sa] [-makeonly] [-nomakeclean] [-hip|-nocuda|-cpponly] [-bsmonly|-nobsm|-scalingonly|-blasonly|-blasandscalingonly|-graphsonly]"
  echo "Run tests and check all logs"
  echo ""
  echo "Usage (2): $0 -checkonly"
  echo "Check existing logs without running any tests"
  exit 1
}

# Parse command line arguments
checkonly=0
rndhst=-curhst
sm=1
bsm=1
scaling=1
blas=1
graphsonly=0
if [ "$1" == "-checkonly" ]; then
  # Check existing logs without running any tests?
  checkonly=1
  shift
  if [ "$1" != "" ]; then usage; fi
fi
while [ "${checkonly}" == "0" ] && [ "$1" != "" ]; do
  if [ "$1" == "-e" ]; then
    # Fail on error?
    set -e
    shift
  elif [ "$1" == "-sa" ]; then
    # Use standalone_cudacpp builds instead of madevent+cudacpp?
    opts+=" -sa"
    suff=".sa"
    shift
  elif [ "$1" == "-makeonly" ]; then
    # Only build all tests instead of building and running them?
    opts+=" -makeonly"
    shift
  elif [ "$1" == "-nomakeclean" ]; then
    # Skip -makeclean (e.g. for brand new generated/downloaded code)
    makeclean=$1
    shift
  elif [ "$1" == "-hip" ]; then
    if [ "${bblds}" != "" ] && [ "${bblds}" != "$1" ]; then echo "ERROR! Incompatible option $1: backend builds are already defined as '$bblds'"; usage; fi
    bblds="$1"
    shift
  elif [ "$1" == "-nocuda" ]; then
    if [ "${bblds}" != "" ] && [ "${bblds}" != "$1" ]; then echo "ERROR! Incompatible option $1: backend builds are already defined as '$bblds'"; usage; fi
    bblds="$1"
    shift
  elif [ "$1" == "-cpponly" ]; then
    if [ "${bblds}" != "" ] && [ "${bblds}" != "$1" ]; then echo "ERROR! Incompatible option $1: backend builds are already defined as '$bblds'"; usage; fi
    bblds="$1"
    shift
  elif [ "$1" == "-bsmonly" ] && [ "${sm}${scaling}${bsm}${blas}" == "1111" ]; then
    sm=0
    bsm=1
    scaling=0
    blas=0
    shift
  elif [ "$1" == "-nobsm" ] && [ "${sm}${scaling}${bsm}${blas}" == "1111" ]; then
    sm=1
    bsm=0
    scaling=1
    blas=1
    shift
  elif [ "$1" == "-scalingonly" ] && [ "${sm}${scaling}${bsm}${blas}" == "1111" ]; then
    sm=0
    bsm=0
    scaling=1
    blas=0
    shift
  elif [ "$1" == "-blasonly" ] && [ "${blas}${scaling}${bsm}${blas}" == "1111" ]; then
    sm=0
    bsm=0
    scaling=0
    blas=1
    shift
  elif [ "$1" == "-blasandscalingonly" ] && [ "${blas}${scaling}${bsm}${blas}" == "1111" ]; then
    sm=0
    bsm=0
    scaling=1
    blas=1
    shift
  elif [ "$1" == "-graphsonly" ]; then
    graphsonly=1
    shift
  else
    usage
  fi
done

# Graphs only?
if [ "${graphsonly}" == "1" ]; then
  if [ "${blas}${scaling}${bsm}${blas}" == "1111" ]; then
    sm=0
    bsm=0
    scaling=0
    blas=0
  else
    usage
  fi
fi

# Check logs
function checklogs()
{
  cd $scrdir/..
  # Print out any errors in the logs (exclude scaling logs)
  if ! egrep -i '(error|fault|failed)' ./tput/logs_*/*.txt; then echo "No errors found in logs"; fi
  # Print out any FPEs or '{ }' in the logs
  echo
  if ! egrep '(^Floating Point Exception|{ })' tput/logs* -r; then echo "No FPEs or '{ }' found in logs"; fi
  # Print out any aborts in the logs (exclude scaling logs)
  echo
  txt=$(grep Abort ./tput/logs_*/*.txt | sed "s|\:.*SubProcesses/P|: P|")
  if [ "${txt}" == "" ]; then
    echo "No aborts found in logs"
  else
    echo "${txt}"
  fi
  # Print out any asserts/aborts in scaling logs
  echo
  txt=$(egrep -i '(abort|assert)' ./tput/logs_*/*.scaling | sed "s|\:.*SubProcesses/P|: P|" | sort -u)
  if [ "${txt}" == "" ]; then
    echo "No aborts or asserts found in scaling logs"
  else
    echo "${txt}"
  fi

  # Print out the MEK channelid debugging output (except for '{ }')
  echo
  \grep MEK ${scrdir}/logs_*/* | sed "s|${scrdir}/logs_||" | grep -v '{ }' | sed 's|_mad.*DEBUG:||' | sort -u
}
if [ "${checkonly}" != "0" ]; then checklogs; exit 0; fi

# Define builds
if [ "$bblds" == "-nocuda" ]; then
  # Random numbers use common instead of curand
  rndhst=-common
  opts+=" -nocuda"
elif [ "$bblds" == "-cpponly" ]; then
  # Random numbers use common instead of curand
  rndhst=-common
  opts+=" -cpponly"
elif [ "$bblds" == "-hip" ]; then # NB: currently (Sep 2024) this is identical to -nocuda
  #### Random numbers use hiprand instead of curand?
  #### This needs ROCm 6.2 (see https://github.com/ROCm/hipRAND/issues/76)
  ###rndhst=-hirhst
  # Random numbers use common (not hiprand) instead of curand
  rndhst=-common
  opts+=" -nocuda"
fi

# This is a script to launch in one go all tests for the main physics processes in this repository
# (Initially it reproduced the logs in tput at the time of commit c0c276840, just before the large alphas PR #434)

cd $scrdir/..
started="STARTED  AT $(date)"

# (+36: 36/165) Six logs (double/mixed/float x hrd0/hrd1 x inl0) in each of the six SM processes (always force rebuild)
\rm -rf gg_ttggg${suff}/lib/build.none_*
cmd="./tput/teeThroughputX.sh -dmf -hrd -makej -eemumu -ggtt -ggttg -ggttgg -ggttggg -gqttq -makeclean ${opts}"
tmp1=$(mktemp)
if [ "${sm}" == "1" ]; then
  $cmd; status=$?
  ls -ltr ee_mumu${suff}/lib/build.none_*_inl0_hrd* gg_tt${suff}/lib/build.none_*_inl0_hrd* gg_tt*g${suff}/lib/build.none_*_inl0_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp1
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended1="$cmd\nENDED(1) AT $(date) [Status=$status]"

# (+18: 54/165) Three scaling logs (double/mixed/float x hrd0 x inl0) in each of the six SM processes
cmd="./tput/teeThroughputX.sh -dmf -makej -eemumu -ggtt -ggttg -ggttgg -ggttggg -gqttq -scaling ${makeclean} ${opts}"
if [ "${scaling}" == "1" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended1sc="$cmd\nENDED(1-scaling) AT $(date) [Status=$status]"

# (+6: 60/165) Three extra logs (double/mixed/float x hrd0 x inl0 + blasOn) only in two of the six SM processes
cmd="./tput/teeThroughputX.sh -ggtt -ggttgg -dmf -blasOn ${makeclean} ${opts}" # no rebuild needed
if [ "${blas}" == "1" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended2="$cmd\nENDED(2) AT $(date) [Status=$status]"

# (+12: 72/165) Three scaling logs (double/mixed/float x hrd0 x inl0 + blasOn) only in four of the six SM processes
cmd="./tput/teeThroughputX.sh -ggtt -ggttg -ggttgg -ggttggg -dmf -blasOn -scaling ${makeclean} ${opts}"
if [ "${blas}" == "1" ] || [ "${scaling}" == "1" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended2sc="$cmd\nENDED(2-scaling) AT $(date) [Status=$status]"

# (+6: 78/165) Three extra logs (double/mixed/float x hrd0 x inl0 + graphs) only in two of the six SM processes
cmd="./tput/teeThroughputX.sh -ggttgg -ggttggg -dmf -useGraphs ${makeclean} ${opts}"
if [ "${sm}" == "1" ] || [ "${graphsonly}" == "1" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended3="$cmd\nENDED(3) AT $(date) [Status=$status]"

# (+12: 90/165) Three scaling logs (double/mixed/float x hrd0 x inl0 + graphs) only in four of the six SM processes
cmd="./tput/teeThroughputX.sh -ggtt -ggttg -ggttgg -ggttggg -dmf -useGraphs -scaling ${makeclean} ${opts}"
if [ "${scaling}" == "1" ] || [ "${graphsonly}" == "1" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended3sc="$cmd\nENDED(3-scaling) AT $(date) [Status=$status]"

# (+12: 102/165) Four extra logs (double/float x hrd0/hrd1 x inl1) only in three of the six SM processes
\rm -rf gg_ttg${suff}/lib/build.none_*
\rm -rf gg_ttggg${suff}/lib/build.none_*
cmd="./tput/teeThroughputX.sh -d_f -hrd -makej -eemumu -ggtt -ggttgg -inlonly ${makeclean} ${opts}"
tmp4=$(mktemp)
if [ "${sm}" == "1" ]; then
  $cmd; status=$?
  ls -ltr ee_mumu${suff}/lib/build.none_*_inl1_hrd* gg_tt${suff}/lib/build.none_*_inl1_hrd* gg_tt*g${suff}/lib/build.none_*_inl1_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp4
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended4="$cmd\nENDED(4) AT $(date) [Status=$status]"

# (+12: 114/165) Two extra logs (double/float x hrd0 x inl0 + bridge) in all six SM processes
cmd="./tput/teeThroughputX.sh -makej -eemumu -ggtt -ggttg -ggttgg -ggttggg -gqttq -d_f -bridge ${makeclean} ${opts}"
if [ "${sm}" == "1" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended5="$cmd\nENDED(5) AT $(date) [Status=$status]"

# (+6: 120/165) Two extra logs (double/float x hrd0 x inl0 + rmbhst) only in three of the six SM processes
cmd="./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -d_f -rmbhst ${makeclean} ${opts}"
if [ "${sm}" == "1" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended6="$cmd\nENDED(6) AT $(date) [Status=$status]"

# (+6: 126/165) Two extra logs (double/float x hrd0 x inl0 + rndhst) only in three of the six SM processes
cmd="./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -d_f ${rndhst} ${makeclean} ${opts}"
if [ "${sm}" == "1" ] && [ "${rndhst}" != "-common" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended7="$cmd\nENDED(7) AT $(date) [Status=$status]"

# (+6: 132/165) Two extra logs (double/float x hrd0 x inl0 + common) only in three of the six SM processes
cmd="./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -d_f -common ${makeclean} ${opts}"
if [ "${sm}" == "1" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended8="$cmd\nENDED(8) AT $(date) [Status=$status]"

# (+6: 138/165) Three extra logs (double/float x hrd0 x inl0 + noBlas) only in two of the six SM processes (always force rebuild)
cmd="./tput/teeThroughputX.sh -ggtt -ggttgg -dmf -noBlas -makeclean ${opts}"
if [ "${blas}" == "1" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended9="$cmd\nENDED(9) AT $(date) [Status=$status]"

# (+24: 162/165) Six extra logs (double/mixed/float x hrd0/hrd1 x inl0) only in the four BSM processes (always force rebuild)
cmd="./tput/teeThroughputX.sh -dmf -hrd -makej -susyggtt -susyggt1t1 -smeftggtttt -heftggbb -makeclean ${opts}"
tmp10=$(mktemp)
if [ "${bsm}" == "1" ]; then
  $cmd; status=$?
  ls -ltr susy_gg_tt${suff}/lib/build.none_*_inl0_hrd* susy_gg_t1t1${suff}/lib/build.none_*_inl0_hrd* smeft_gg_tttt${suff}/lib/build.none_*_inl0_hrd* heft_gg_bb${suff}/lib/build.none_*_inl0_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp10
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended10="$cmd\nENDED(10) AT $(date) [Status=$status]"

# (+1: 164/165) Two logs (mixed x hrd0 x inl0 x dcd0/dcd1) in ggttg5 (always force rebuild)
cmd="./tput/teeThroughputX.sh -makej -ggttg5 -dcd -makeclean ${opts}"
if [ "${sm}" == "1" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended11="$cmd\nENDED(11) AT $(date) [Status=$status]"

# (+1: 165/165) One log (mixed x hrd0 x inl0 + graphs) in ggttg5 (no need to rebuild)
cmd="./tput/teeThroughputX.sh -makej -ggttg5 -useGraphs ${opts}"
if [ "${sm}" == "1" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended12="$cmd\nENDED(12) AT $(date) [Status=$status]"

echo
echo "Build(1):"
cat $tmp1
echo
echo "Build(4):"
cat $tmp4
echo
echo "Build(10):"
cat $tmp10
echo
echo -e "$started"
echo -e "$ended1"
echo -e "$ended1sc"
echo -e "$ended2"
echo -e "$ended2sc"
echo -e "$ended3"
echo -e "$ended3sc"
echo -e "$ended4"
echo -e "$ended5"
echo -e "$ended6"
echo -e "$ended7"
echo -e "$ended8"
echo -e "$ended9"
echo -e "$ended10"
echo -e "$ended11"
echo -e "$ended12"

# Finally check all logs
echo
checklogs
