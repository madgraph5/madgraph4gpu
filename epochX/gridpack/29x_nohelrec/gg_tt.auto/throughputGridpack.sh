#!/bin/bash

function usage()
{
  echo "Usage: $0 [-cleanbuild|-skiprun] [nevt(def=1000) [rndm(def=1)]]"
  exit 1
}

cleanbuild=0
skiprun=0
if [ "$1" == "-cleanbuild" ]; then
  cleanbuild=1
  shift
elif [ "$1" == "-skiprun" ]; then
  skiprun=1
  shift
fi
nevt=1000
rndm=1

re='^[0-9]+$'
if [ "$1" != "" ] && [[ $1 =~ $re ]]; then nevt=$1; shift; fi
if [ "$1" != "" ] && [[ $1 =~ $re ]]; then rndm=$1; shift; fi
if [ "$1" != "" ]; then usage; fi
###echo cleanbuild=$cleanbuild
###echo skiprun=$skiprun

dir=$(cd $(dirname $0); pwd)
cd $dir

if [ "$cleanbuild" == "1" ]; then
  \rm -rf madevent
  git checkout madevent
  ./rebuild.sh
fi

if [ "$skiprun" == "0" ]; then
  \rm -f events.lhe.gz
  ./run.sh $nevt $rndm
fi

mkdir -p ../tput
out=$(cd ../tput; pwd)/log_$(echo $(basename $dir) | sed 's/.auto/_auto/').txt
echo "Output: $out"

\rm -f $out
touch $out

echo

for Gdir in madevent/SubProcesses/P1_*/G*; do
  log=$Gdir/counters_log.txt
  if [ -f $log ]; then printf "%-5s " $(basename $Gdir) | tee -a $out; cat $log | grep '^PROGRAM' | tee -a $out; fi
done
cat madevent/SubProcesses/P1_*/G*/counters_log.txt | grep '^PROGRAM' | awk 'BEGIN{stot=0;ctot=0};{stot+=substr($3,0,length($3)-1)};END{printf "TOTAL PROGRAM    : %9.4fs\n",stot}' | tee -a $out

echo | tee -a $out

for Gdir in madevent/SubProcesses/P1_*/G*; do
  log=$Gdir/counters_log.txt
  if [ -f $log ]; then printf "%-5s " $(basename $Gdir) | tee -a $out; cat $log | grep '^MATRIX1(a)' | tee -a $out; fi
done
cat madevent/SubProcesses/P1_*/G*/counters_log.txt | grep '^MATRIX1(a)' | awk 'BEGIN{stot=0;ctot=0};{stot+=substr($3,0,length($3)-1);ctot+=$5};END{printf "TOTAL MATRIX1(a) : %9.4fs for %8d MATRIX1 calls  => throughput is %8.2E calls/s\n",stot,ctot,ctot/stot}' | tee -a $out

echo | tee -a $out

for Gdir in madevent/SubProcesses/P1_*/G*; do
  log=$Gdir/counters_log.txt
  if [ -f $log ]; then printf "%-5s " $(basename $Gdir) | tee -a $out; cat $log | grep '^MATRIX1(b)' | tee -a $out; fi
done
cat madevent/SubProcesses/P1_*/G*/counters_log.txt | grep '^MATRIX1(b)' | awk 'BEGIN{stot=0;ctot=0};{stot+=substr($3,0,length($3)-1);ctot+=$5};END{printf "TOTAL MATRIX1(b) : %9.4fs for %8d SMATRIX1 calls => throughput is %8.2E calls/s\n",stot,ctot,ctot/stot}' | tee -a $out

echo | tee -a $out

for Gdir in madevent/SubProcesses/P1_*/G*; do
  log=$Gdir/counters_log.txt
  if [ -f $log ]; then printf "%-5s " $(basename $Gdir) | tee -a $out; cat $log | grep '^SMATRIX1' | tee -a $out; fi
done
cat madevent/SubProcesses/P1_*/G*/counters_log.txt | grep ^SMATRIX1 | awk 'BEGIN{stot=0;ctot=0};{stot+=substr($3,0,length($3)-1);ctot+=$5};END{printf "TOTAL SMATRIX1   : %9.4fs for %8d SMATRIX1 calls => throughput is %8.2E calls/s\n",stot,ctot,ctot/stot}' | tee -a $out

echo

###exit 0 # comment out to skip final cleanup
\rm -f events.lhe.gz
\rm -rf madevent/Events/GridRun_1/
\rm -f madevent/SubProcesses/P1*/ajob1
\rm -f madevent/SubProcesses/P1*/G*/input_sg.txt
\rm -f madevent/SubProcesses/P1*/G*/moffset.dat
git checkout madevent/index.html madevent/Cards/grid_card.dat madevent/SubProcesses/randinit madevent/SubProcesses/P1*/G* |& grep -v "Updated 0 paths"
###\rm -f madevent/SubProcesses/P1*/madevent
