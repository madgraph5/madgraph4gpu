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
  else
    echo "Usage: $0 [-short] [-e] [-sa] [-makeonly]"
    exit 1
  fi
done

# This is a script to launch in one go all tests for the (4 or) 5 main processes in this repository
# It reproduces the logs in tput at the time of commit c0c276840654575d9fa0c3f3c4a0088e57764dbc
# This is the commit just before the large alphas PR #434

cd $scrdir/..
started="STARTED  AT $(date)"

# (36/78) Six logs (double/float/mixed x hrd0/hrd1 x inl0) in each of the six processes
\rm -rf gg_ttggg${suff}/lib/build.none_*
cmd="./tput/teeThroughputX.sh -mix -hrd -makej -eemumu -ggtt -ggttg -ggttgg -gqttq $ggttggg -makeclean ${opts}"
$cmd; status=$?
ended1="$cmd\nENDED(1) AT $(date) [Status=$status]"
tmp1=$(mktemp)
ls -ltr ee_mumu${suff}/lib/build.none_*_inl0_hrd* gg_tt${suff}/lib/build.none_*_inl0_hrd* gg_tt*g${suff}/lib/build.none_*_inl0_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp1

# (48/78) Four extra logs (double/float x hrd0/hrd1 x inl1) only in three of the six processes
\rm -rf gg_ttg${suff}/lib/build.none_*
\rm -rf gg_ttggg${suff}/lib/build.none_*
cmd="./tput/teeThroughputX.sh -flt -hrd -makej -eemumu -ggtt -ggttgg -inlonly -makeclean ${opts}"
$cmd; status=$?
ended2="$cmd\nENDED(2) AT $(date) [Status=$status]"
tmp2=$(mktemp)
ls -ltr ee_mumu${suff}/lib/build.none_*_inl1_hrd* gg_tt${suff}/lib/build.none_*_inl1_hrd* gg_tt*g${suff}/lib/build.none_*_inl1_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp2

# (60/78) Two extra logs (double/float x hrd0 x inl0 + bridge) in all six processes (rebuild from cache)
cmd="./tput/teeThroughputX.sh -makej -eemumu -ggtt -ggttg -gqttq -ggttgg $ggttggg -flt -bridge -makeclean ${opts}"
$cmd; status=$?
ended3="$cmd\nENDED(3) AT $(date) [Status=$status]"

# (66/78) Two extra logs (double/float x hrd0 x inl0 + rmbhst) only in three of the six processes (no rebuild needed)
cmd="./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt -rmbhst ${opts}"
$cmd; status=$?
ended4="$cmd\nENDED(4) AT $(date) [Status=$status]"

# (72/78) Two extra logs (double/float x hrd0 x inl0 + curhst) only in three of the six processes (no rebuild needed)
cmd="./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt -curhst ${opts}"
$cmd; status=$?
ended5="$cmd\nENDED(5) AT $(date) [Status=$status]"

# (78/78) Two extra logs (double/float x hrd0 x inl0 + common) only in three of the six processes (no rebuild needed)
cmd="./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt -common ${opts}"
$cmd; status=$?
ended6="$cmd\nENDED(6) AT $(date) [Status=$status]"

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

if [ "$ggttggg" == "" ]; then
  echo
  echo "To complete the test for ggttggg type:"
  echo "  ./tput/teeThroughputX.sh -flt -hrd -makej -ggttggg -makeclean ${opts}"
  echo "  ./tput/teeThroughputX.sh -makej -ggttggg -flt -bridge -makeclean ${opts}"
fi

# Print out any errors in the logs
echo
if ! egrep -i 'error' ./tput/logs_* -r; then echo "No errors found in logs"; fi
