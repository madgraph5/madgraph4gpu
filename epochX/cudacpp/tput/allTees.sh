#!/bin/bash

scrdir=$(cd $(dirname $0); pwd)

# Use the -mad version of code and tee scripts
mad=-mad

# Short (no ggttggg) or long version?
if [ "$1" == "-short" ]; then
  ggttggg=
  shift
elif [ "$1" == "" ]; then
  ggttggg=-ggttggg
  shift
else
  echo "Usage: $0 [-short]"
  exit 1
fi

# This is a script to launch in one go all tests for the (4 or) 5 main processes in this repository
# It reproduces the logs in tput at the time of commit c0c276840654575d9fa0c3f3c4a0088e57764dbc
# This is the commit just before the large alphas PR #434

cd $scrdir/..
started="STARTED  AT $(date)"

# Four logs (double/float x hrd0/hrd1 x inl0) in each of the five processes
\rm -rf gg_ttggg/lib/build.none_*
./tput/teeThroughputX.sh ${mad} -flt -hrd -makej -makeclean -eemumu -ggtt -ggttg -ggttgg $ggttggg
ended1="ENDED(1) AT $(date)"
tmp1=$(mktemp)
ls -ltr ee_mumu/lib/build.none_*_inl0_hrd* gg_tt/lib/build.none_*_inl0_hrd* gg_tt*g/lib/build.none_*_inl0_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp1

# Four extra logs (double/float x hrd0/hrd1 x inl1) only in three of the five processes
\rm -rf gg_ttg/lib/build.none_*
\rm -rf gg_ttggg/lib/build.none_*
./tput/teeThroughputX.sh ${mad} -flt -hrd -makej -makeclean -eemumu -ggtt -ggttgg -inlonly
ended2="ENDED(2) AT $(date)"
tmp2=$(mktemp)
ls -ltr ee_mumu/lib/build.none_*_inl1_hrd* gg_tt/lib/build.none_*_inl1_hrd* gg_tt*g/lib/build.none_*_inl1_hrd* | egrep -v '(total|\./|\.build|_common|^$)' > $tmp2

# Two extra logs (double/float x hrd0 x inl0 + bridge) only in three of the five processes (no rebuild needed)
./tput/teeThroughputX.sh ${mad} -eemumu -ggtt -ggttgg -flt -bridge
ended3="ENDED(3) AT $(date)"

# Two extra logs (double/float x hrd0 x inl0 + rmbhst) only in three of the five processes (no rebuild needed)
./tput/teeThroughputX.sh ${mad} -eemumu -ggtt -ggttgg -flt -rmbhst
ended4="ENDED(4) AT $(date)"

# Two extra logs (double/float x hrd0 x inl0 + curhst) only in three of the five processes (no rebuild needed)
./tput/teeThroughputX.sh ${mad} -eemumu -ggtt -ggttgg -flt -curhst
ended5="ENDED(5) AT $(date)"

# Two extra logs (double/float x hrd0 x inl0 + common) only in three of the five processes (no rebuild needed)
./tput/teeThroughputX.sh ${mad} -eemumu -ggtt -ggttgg -flt -common
ended6="ENDED(6) AT $(date)"

echo
echo "Build(1):"
cat $tmp1
echo
echo "Build(2):"
cat $tmp2
echo
echo "$started"
echo "$ended1"
echo "$ended2"
echo "$ended3"
echo "$ended4"
echo "$ended5"

if [ "$ggttggg" == "" ]; then
  echo
  echo "To complete the test for ggttggg type:"
  echo "  ./tput/teeThroughputX.sh ${mad} -flt -hrd -makej -makeclean -ggttggg"
fi
