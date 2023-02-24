#!/bin/bash

scrdir=$(cd $(dirname $0); pwd)

# By default, use the madevent+cudacpp version of code and tee scripts
sa=
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
    sa=-sa
    suff=".sa"
    shift
  else
    echo "Usage: $0 [-short] [-e] [-sa]"
    exit 1
  fi
done

# This is a script to launch in one go all tests for the (4 or) 5 main processes in this repository
# It reproduces the logs in tput at the time of commit c0c276840654575d9fa0c3f3c4a0088e57764dbc
# This is the commit just before the large alphas PR #434

cd $scrdir/..
started="STARTED  AT $(date)"

# (30/70) Four logs (double/float/mixed x hrd0/hrd1 x inl0) in each of the five processes
\rm -rf gg_ttggg${suff}/lib/build.none_*
cmd="./tput/teeThroughputX.sh -mix -hrd -makej -eemumu -ggtt -ggttg -ggttgg $ggttggg -makeclean ${sa}"
$cmd; status=$?
ended1="$cmd\nENDED(1) AT $(date) [Status=$status]"
tmp1=$(mktemp)
ls -ltr ee_mumu${suff}/lib/build.none_*_inl0_hrd* gg_tt${suff}/lib/build.none_*_inl0_hrd* gg_tt*g${suff}/lib/build.none_*_inl0_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp1

# (42/70) Four extra logs (double/float x hrd0/hrd1 x inl1) only in three of the five processes
\rm -rf gg_ttg${suff}/lib/build.none_*
\rm -rf gg_ttggg${suff}/lib/build.none_*
cmd="./tput/teeThroughputX.sh -flt -hrd -makej -eemumu -ggtt -ggttgg -inlonly -makeclean ${sa}"
$cmd; status=$?
ended2="$cmd\nENDED(2) AT $(date) [Status=$status]"
tmp2=$(mktemp)
ls -ltr ee_mumu${suff}/lib/build.none_*_inl1_hrd* gg_tt${suff}/lib/build.none_*_inl1_hrd* gg_tt*g${suff}/lib/build.none_*_inl1_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp2

# (52/70) Two extra logs (double/float x hrd0 x inl0 + bridge) in all five processes (rebuild from cache)
cmd="./tput/teeThroughputX.sh -makej -eemumu -ggtt -ggttg -ggttgg $ggttggg -flt -bridge -makeclean ${sa}"
$cmd; status=$?
ended3="$cmd\nENDED(3) AT $(date) [Status=$status]"

# (58/70) Two extra logs (double/float x hrd0 x inl0 + rmbhst) only in three of the five processes (no rebuild needed)
cmd="./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt -rmbhst ${sa}"
$cmd; status=$?
ended4="$cmd\nENDED(4) AT $(date) [Status=$status]"

# (64/70) Two extra logs (double/float x hrd0 x inl0 + curhst) only in three of the five processes (no rebuild needed)
cmd="./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt -curhst ${sa}"
$cmd; status=$?
ended5="$cmd\nENDED(5) AT $(date) [Status=$status]"

# (70/70) Two extra logs (double/float x hrd0 x inl0 + common) only in three of the five processes (no rebuild needed)
cmd="./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt -common ${sa}"
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
  echo "  ./tput/teeThroughputX.sh -flt -hrd -makej -ggttggg -makeclean ${sa}"
  echo "  ./tput/teeThroughputX.sh -makej -ggttggg -flt -bridge -makeclean ${sa}"
fi

# Print out any errors in the logs
echo
if ! egrep -i 'error' ./tput/logs_* -r; then echo "No errors found in logs"; fi
