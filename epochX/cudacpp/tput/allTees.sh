#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Apr 2022) for the MG5aMC CUDACPP plugin.

scrdir=$(cd $(dirname $0); pwd)

# By default, use the madevent+cudacpp version of code and tee scripts (use -sa to use the standalone version instead)
# By default, build and run all tests (use -makeonly to only build all tests)
opts=
suff=".mad"

# Parse command line arguments
ggttggg=-ggttggg
rndhst=-curhst
bsm=
while [ "$1" != "" ]; do
  if [ "$1" == "-short" ]; then
    # Short (no ggttggg) or long version?
    ggttggg=
    shift
  elif [ "$1" == "-e" ]; then
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
  elif [ "$1" == "-hip" ]; then
    #### Random numbers use hiprand instead of curand?
    ###rndhst=-hirhst
    # See https://github.com/ROCm/hipRAND/issues/76
    # Random numbers use common (not hiprand) instead of curand?
    rndhst=-common
    shift
  elif [ "$1" == "-bsmonly" ] && [ "$bsm" != "-nobsm" ]; then
    bsm=$1
    shift
  elif [ "$1" == "-nobsm" ] && [ "$bsm" != "-bsmonly" ]; then
    bsm=$1
    shift
  else
    echo "Usage: $0 [-short] [-e] [-sa] [-makeonly] [-hip] [-bsmonly|-nobsm]"
    exit 1
  fi
done

# This is a script to launch in one go all tests for the (4 or) 5 main processes in this repository
# It reproduces the logs in tput at the time of commit c0c276840654575d9fa0c3f3c4a0088e57764dbc
# This is the commit just before the large alphas PR #434

cd $scrdir/..
started="STARTED  AT $(date)"

# (36/102) Six logs (double/float/mixed x hrd0/hrd1 x inl0) in each of the six SM processes
\rm -rf gg_ttggg${suff}/lib/build.none_*
cmd="./tput/teeThroughputX.sh -mix -hrd -makej -eemumu -ggtt -ggttg -ggttgg -gqttq $ggttggg -makeclean ${opts}"
tmp1=$(mktemp)
if [ "${bsm}" != "-bsmonly" ]; then
  $cmd; status=$?
  ls -ltr ee_mumu${suff}/lib/build.none_*_inl0_hrd* gg_tt${suff}/lib/build.none_*_inl0_hrd* gg_tt*g${suff}/lib/build.none_*_inl0_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp1
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended1="$cmd\nENDED(1) AT $(date) [Status=$status]"

# (48/102) Four extra logs (double/float x hrd0/hrd1 x inl1) only in three of the six SM processes
\rm -rf gg_ttg${suff}/lib/build.none_*
\rm -rf gg_ttggg${suff}/lib/build.none_*
cmd="./tput/teeThroughputX.sh -flt -hrd -makej -eemumu -ggtt -ggttgg -inlonly -makeclean ${opts}"
tmp2=$(mktemp)
if [ "${bsm}" != "-bsmonly" ]; then
  $cmd; status=$?
  ls -ltr ee_mumu${suff}/lib/build.none_*_inl1_hrd* gg_tt${suff}/lib/build.none_*_inl1_hrd* gg_tt*g${suff}/lib/build.none_*_inl1_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp2
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended2="$cmd\nENDED(2) AT $(date) [Status=$status]"

# (60/102) Two extra logs (double/float x hrd0 x inl0 + bridge) in all six SM processes (rebuild from cache)
cmd="./tput/teeThroughputX.sh -makej -eemumu -ggtt -ggttg -gqttq -ggttgg $ggttggg -flt -bridge -makeclean ${opts}"
if [ "${bsm}" != "-bsmonly" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended3="$cmd\nENDED(3) AT $(date) [Status=$status]"

# (66/102) Two extra logs (double/float x hrd0 x inl0 + rmbhst) only in three of the six SM processes (no rebuild needed)
cmd="./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt -rmbhst ${opts}"
if [ "${bsm}" != "-bsmonly" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended4="$cmd\nENDED(4) AT $(date) [Status=$status]"

# (72/102) Two extra logs (double/float x hrd0 x inl0 + rndhst) only in three of the six SM processes (no rebuild needed)
cmd="./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt ${rndhst} ${opts}"
if [ "${bsm}" != "-bsmonly" ] && [ "${rndhst}" != "-common" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended5="$cmd\nENDED(5) AT $(date) [Status=$status]"

# (78/102) Two extra logs (double/float x hrd0 x inl0 + common) only in three of the six SM processes (no rebuild needed)
cmd="./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt -common ${opts}"
if [ "${bsm}" != "-bsmonly" ]; then
  $cmd; status=$?
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended6="$cmd\nENDED(6) AT $(date) [Status=$status]"

# (102/102) Six extra logs (double/float/mixed x hrd0/hrd1 x inl0) only in the four BSM processes
cmd="./tput/teeThroughputX.sh -mix -hrd -makej -susyggtt -susyggt1t1 -smeftggtttt -heftggbb -makeclean ${opts}"
tmp3=$(mktemp)
if [ "${bsm}" != "-nobsm" ]; then
  $cmd; status=$?
  ls -ltr susy_gg_tt${suff}/lib/build.none_*_inl0_hrd* susy_gg_t1t1${suff}/lib/build.none_*_inl0_hrd* smeft_gg_tttt${suff}/lib/build.none_*_inl0_hrd* heft_gg_bb${suff}/lib/build.none_*_inl0_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp2
else
  cmd="SKIP '$cmd'"; echo $cmd; status=$?
fi
ended7="$cmd\nENDED(7) AT $(date) [Status=$status]"

echo
echo "Build(1):"
cat $tmp1
echo
echo "Build(2):"
cat $tmp2
echo
echo -e "$started"
echo -e "$ended1"
echo -e "$ended2"
echo -e "$ended3"
echo -e "$ended4"
echo -e "$ended5"
echo -e "$ended6"
echo -e "$ended7"

if [ "$ggttggg" == "" ]; then
  echo
  echo "To complete the test for ggttggg type:"
  echo "  ./tput/teeThroughputX.sh -flt -hrd -makej -ggttggg -makeclean ${opts}"
  echo "  ./tput/teeThroughputX.sh -makej -ggttggg -flt -bridge -makeclean ${opts}"
fi

# Print out any errors in the logs
echo
if ! egrep -i 'error' ./tput/logs_* -r; then echo "No errors found in logs"; fi
