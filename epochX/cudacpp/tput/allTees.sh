#!/bin/bash

scrdir=$(cd $(dirname $0); pwd)

# By default, use the madevent+cudacpp version of code and tee scripts
sa=

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

# (20/56) Four logs (double/float x hrd0/hrd1 x inl0) in each of the five processes
\rm -rf gg_ttggg/lib/build.none_*
cmd="./tput/teeThroughputX.sh ${sa} -flt -hrd -makej -makeclean -eemumu -ggtt -ggttg -ggttgg $ggttggg"
$cmd; status=$?
ended1="$cmd\nENDED(1) AT $(date) [Status=$status]"
tmp1=$(mktemp)
ls -ltr ee_mumu/lib/build.none_*_inl0_hrd* gg_tt/lib/build.none_*_inl0_hrd* gg_tt*g/lib/build.none_*_inl0_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp1

# (32/56) Four extra logs (double/float x hrd0/hrd1 x inl1) only in three of the five processes
\rm -rf gg_ttg/lib/build.none_*
\rm -rf gg_ttggg/lib/build.none_*
cmd="./tput/teeThroughputX.sh ${sa} -flt -hrd -makej -makeclean -eemumu -ggtt -ggttgg -inlonly"
$cmd; status=$?
ended2="$cmd\nENDED(2) AT $(date) [Status=$status]"
tmp2=$(mktemp)
ls -ltr ee_mumu/lib/build.none_*_inl1_hrd* gg_tt/lib/build.none_*_inl1_hrd* gg_tt*g/lib/build.none_*_inl1_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp2

# (38/56) Two extra logs (double/float x hrd0 x inl0 + bridge) only in three of the five processes (rebuild from cache)
cmd="./tput/teeThroughputX.sh ${sa} -makej -makeclean -eemumu -ggtt -ggttgg -flt -bridge"
$cmd; status=$?
ended3="$cmd\nENDED(3) AT $(date) [Status=$status]"

# (44/56) Two extra logs (double/float x hrd0 x inl0 + rmbhst) only in three of the five processes (no rebuild needed)
cmd="./tput/teeThroughputX.sh ${sa} -eemumu -ggtt -ggttgg -flt -rmbhst"
$cmd; status=$?
ended4="$cmd\nENDED(4) AT $(date) [Status=$status]"

# (50/56) Two extra logs (double/float x hrd0 x inl0 + curhst) only in three of the five processes (no rebuild needed)
cmd="./tput/teeThroughputX.sh ${sa} -eemumu -ggtt -ggttgg -flt -curhst"
$cmd; status=$?
ended5="$cmd\nENDED(5) AT $(date) [Status=$status]"

# (56/56) Two extra logs (double/float x hrd0 x inl0 + common) only in three of the five processes (no rebuild needed)
cmd="./tput/teeThroughputX.sh ${sa} -eemumu -ggtt -ggttgg -flt -common"
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
  echo "  ./tput/teeThroughputX.sh ${sa} -flt -hrd -makej -makeclean -ggttggg"
fi
