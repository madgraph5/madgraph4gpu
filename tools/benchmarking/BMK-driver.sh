#!/bin/bash

### TO KILL A TEST:
###kill $(ps -ef | egrep '(mg5amc-madgraph4gpu-2022-bmk.sh|throughputX.sh|check.exe)' | grep -v grep | awk '{print $2}')

# Node-specific configuration
image=oras://registry.cern.ch/hep-workloads/mg5amc-madgraph4gpu-2022-bmk:v0.6
if [ "$(hostname)" == "pmpe04.cern.ch" ]; then
  export SINGULARITY_CACHEDIR=/data/BMK2021/SINGULARITY_CACHEDIR
  unset SINGULARITY_TMPDIR
  resDir=/data/BMK2021/TMP_RESULTS
  tstDir=BMK-pmpe04
  mkdir -p ${tstDir}
  jts=""; for j in 1 2 4 8 16 32 48; do for t in 1; do if [ $((j*t)) -le 64 ]; then jts="$jts [$j,$t]"; fi; done; done 
elif [ "$(hostname)" == "itscrd70.cern.ch" ]; then
  export SINGULARITY_TMPDIR=/scratch/TMP_AVALASSI/
  export SINGULARITY_CACHEDIR=/scratch/SINGULARITY_CACHEDIR
  resDir=/scratch/TMP_RESULTS
  tstDir=BMK-itscrd70
  mkdir -p ${tstDir}
  jts=""; for j in 1 2 4 6; do for t in 1; do if [ $((j*t)) -le 8 ]; then jts="$jts [$j,$t]"; fi; done; done 
else
  echo "ERROR! Unknown host $(hostname)"; exit 1
fi

# Loop over tests
###events=1 # NOT ENOUGH FOR GGTTG AND ALSO GGTT (BMK-1056)
events=10
for jt in $jts; do
  jt=${jt:1:-1}; t=${jt#*,}; j=${jt%,*}
  ###echo "jobs=$j, threads=$t"
  wDir="${tstDir}/sa-cpp$(printf '%s%03i%s%03i%s%03i' '-j' ${j} '-t' ${t} '-e' ${events})"
  ###echo ${resDir}/${wDir}
  \rm -rf ${resDir}/${wDir}
  singularity run -B ${resDir}:/results ${image} --extra-args '-eemumu -ggtt -ggttg -ggttgg -dbl -flt -inl0 -inl1 --cpu' -c${j} -t${t} -e${events} -w /results/${wDir} -W
  ###ls -l ${resDir}/${wDir}
done
