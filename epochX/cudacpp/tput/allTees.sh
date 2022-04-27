#!/bin/bash

scrdir=$(cd $(dirname $0); pwd)

# This is a script to launch in one go all tests for the 5 main processes in this repository
# It reproduces the logs in tput at the time of commit c0c276840654575d9fa0c3f3c4a0088e57764dbc
# This is the commit just before the large alphas PR #434

cd $scrdir/..
started="STARTED AT $(date)"

# Four logs (double/float x hrd0/hrd1 x inl0) in each of the five processes
./tput/teeThroughputX.sh -flt -hrd -makej -makeclean -eemumu -ggtt -ggttg -ggttgg -ggttggg
ended1="ENDED(1) AT $(date)"

# Four extra logs (double/float x hrd0/hrd1 x inl1) only in three of the five processes
./tput/teeThroughputX.sh -flt -hrd -makej -makeclean -eemumu -ggtt -ggttgg -inlonly
ended2="ENDED(2) AT $(date)"

# Two extra logs (double/float x hrd0 x inl0 + bridge) only in three of the five processes (no rebuild needed)
./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt -bridge
ended3="ENDED(3) AT $(date)"

# Two extra logs (double/float x hrd0 x inl0 + rmbhst) only in three of the five processes (no rebuild needed)
./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt -rmbhst
ended4="ENDED(4) AT $(date)"

# Two extra logs (double/float x hrd0 x inl0 + curhst) only in three of the five processes (no rebuild needed)
./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt -curhst
ended5="ENDED(5) AT $(date)"

# Two extra logs (double/float x hrd0 x inl0 + common) only in three of the five processes (no rebuild needed)
./tput/teeThroughputX.sh -eemumu -ggtt -ggttgg -flt -common
ended6="ENDED(6) AT $(date)"

echo
echo "$started"
echo "$ended1"
echo "$ended2"
echo "$ended3"
echo "$ended4"
echo "$ended5"
