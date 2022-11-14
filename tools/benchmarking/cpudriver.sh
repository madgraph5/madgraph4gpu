#!/bin/bash

### TO KILL A TEST:
###kill $(ps -ef | egrep '(mg5amc-madgraph4gpu-2022-bmk.sh|throughputX.sh|check.exe)' | grep -v grep | awk '{print $2}')

startdate=$(date)

# Defaults for all nodes (may be overridden)
image=oras://registry.cern.ch/hep-workloads/mg5amc-madgraph4gpu-2022-bmk:v0.6
extraargs="-eemumu -ggtt -ggttg -ggttgg -dbl -flt -inl0 -inl1 --cpu"
###events=1 # NOT ENOUGH FOR GGTTG AND ALSO GGTT (BMK-1056)
events=10

# Node-specific configuration
if [ "$(hostname)" == "pmpe04.cern.ch" ]; then
  export SINGULARITY_CACHEDIR=/data/BMK2021/SINGULARITY_CACHEDIR
  unset SINGULARITY_TMPDIR
  resDir=/data/BMK2021/TMP_RESULTS
  tstDir=BMK-pmpe04
  jts=""; for j in 1 2 4 8 16 32 48; do for t in 1; do if [ $((j*t)) -le 64 ]; then jts="$jts [$j,$t]"; fi; done; done 
elif [ "$(hostname)" == "itscrd70.cern.ch" ]; then
  export SINGULARITY_TMPDIR=/scratch/TMP_AVALASSI/
  export SINGULARITY_CACHEDIR=/scratch/SINGULARITY_CACHEDIR
  resDir=/scratch/TMP_RESULTS
  tstDir=BMK-itscrd70
  jts=""; for j in 1 2 4 6; do for t in 1; do if [ $((j*t)) -le 8 ]; then jts="$jts [$j,$t]"; fi; done; done 
elif [ "$(hostname)" == "jwlogin08.juwels" ]; then
  export APPTAINER_TMPDIR=/p/scratch/prpb109/avalassi/TMP_AVALASSI/
  export APPTAINER_CACHEDIR=/p/scratch/prpb109/avalassi/SINGULARITY_CACHEDIR
  export SINGULARITY_CACHEDIR=${APPTAINER_CACHEDIR}
  resDir=/p/scratch/prpb109/avalassi/TMP_RESULTS
  tstDir=BMK-jwlogin08
  jts=""; for j in 1 2 5 10 20 40 80; do for t in 1; do if [ $((j*t)) -le 80 ]; then jts="$jts [$j,$t]"; fi; done; done 
  extraargs="-ggttgg -dbl -flt -inl0 --cpu" # SHORTER TESTS ON JUWELS
  events=1 # SHORTER TESTS ON JUWELS
elif [ "$(hostname)" == "bmk-ironic-0731f1ce3b.cern.ch" ]; then
  export SINGULARITY_TMPDIR=/var/benchmark/avalassi/TMP_AVALASSI/
  export SINGULARITY_CACHEDIR=/var/benchmark/avalassi/SINGULARITY_CACHEDIR
  resDir=/var/benchmark/avalassi/TMP_RESULTS
  tstDir=BMK-bmk6130
  jts=""; for j in 1 2 4 8 16 32 64 96; do for t in 1; do if [ $((j*t)) -le 128 ]; then jts="$jts [$j,$t]"; fi; done; done 
else
  echo "ERROR! Unknown host $(hostname)"; exit 1
fi

# Loop over tests
for jt in $jts; do
  jt=${jt:1:-1}; t=${jt#*,}; j=${jt%,*}
  ###echo "jobs=$j, threads=$t"
  wDir="${tstDir}/sa-cpp$(printf '%s%03i%s%03i%s%03i' '-j' ${j} '-t' ${t} '-e' ${events})"
  ###echo ${resDir}/${wDir}
  \rm -rf ${resDir}/${wDir}
  singularity run -B ${resDir}:/results ${image} --extra-args "${extraargs}" -c${j} -t${t} -e${events} -w /results/${wDir} -W
  ###ls -l ${resDir}/${wDir}
done

echo "Started at ${startdate}"
echo "Completed at $(date)"
