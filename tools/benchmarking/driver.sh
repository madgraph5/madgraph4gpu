#!/bin/bash

tstName="check-test"
exe=./build.none/check.exe 

if ! calc -v >& /dev/null; then echo "ERROR! please install calc"; exit 1; fi

#--- [lhcb-knl-scripts] from kill.sh ----

function killExe() {
  exeName=check.exe
  procs=`ps -efH | grep ${USER} | grep $exeName | grep -v grep | awk '{ print $2 }'`
  if [ "$procs" != "" ]; then
    echo Killing $exeName processes $procs
    echo "$procs" | xargs kill -9
    return 1
  else
    echo No $exeName processes to kill
    return 0
  fi
}

function killAndWait() {
  for i in `seq 1 5`; do
    if ! killExe; then
      echo "WARNING! ${tstName} killed: sleep and try again ($i/5)"
      sleep 5
    else
      break
    fi
  done
}

killAndWait

#--- [lhcb-knl-scripts] from createTEST.sh ----

tstDir=$(pwd)/BMKTST
###\rm -rf ${tstDir}
mkdir -p ${tstDir}

jts="[1,1] [1,4]"

\rm -f ${tstDir}/alljts
space=""
for jt in $jts; do
  jt=${jt:1:-1}; t=${jt#*,}; j=${jt%,*}
  echo -n "$space"$(printf '%03i/%03i' $j $t) >> ${tstDir}/alljts
  space=" "
done

for jt in $jts; do
  jt=${jt:1:-1}; t=${jt#*,}; j=${jt%,*}
  echo "jobs=$j, threads=$t"
  ###if [ $(($j * $t)) -gt 300 ]; then continue; fi
  wDir="${tstDir}/${tstName}$(printf '.j%03i.t%03i' $j $t)"
  if [ ! -d $wDir ]; then
    mkdir $wDir
    echo -n $j > $wDir/njobs.txt
  fi
done

#--- [lhcb-knl-scripts] from run.sh ----

function runjob {
  echo c0 ${1} `pwd -P`
  jobno=${1}
  if [ ! -d ${jobno} ]; then mkdir ${jobno}; fi
  cd ${jobno}
  echo c1 ${1} `pwd -P`
  ###$exe -p 2048 256 1 | egrep '(OMP threads|TotalEventsComputed|TOTAL   \(3\))'
  echo "DUMMY TEST!"
  echo c2 ${1} `pwd -P`
}

function runandwait {
  echo b0 `pwd -P`
  echo `date` STARTING jobs  
  for x in `seq \`cat njobs.txt\``; do
    runjob ${x} &
  done
  sleep 1
  touch WAITING
  echo `date`: START WAITING
  ls -l `pwd`/WAITING
  sleep 5
  wait
  \rm -f WAITING
  echo `date` FINISHED jobs  
  echo b1 `pwd -P`
}

cd ${tstDir}

testfile=tests.txt
\rm -f ${testfile}
touch ${testfile}

for jt in `cat alljts`; do
  j=${jt%/*}
  t=${jt#*/}
  echo
  echo `date` starting jobs=${j}, threads=${t}
  wDir=${tstName}.j${j}.t${t}
  if [ -e ${wDir} ]; then
    echo a0 `pwd -P`
    echo j${j}.t${t} >> ${testfile}
    cd ${wDir}
    echo a1 `pwd -P`
    \rm -f WAITING
    ( runandwait ) &
    waiterPid=$! # See https://stackoverflow.com/a/10028986
    while [ ! -f WAITING ]; do sleep 1; done
    timeout=1 # 1 minute
    echo `date`": INFO: TIMEOUT $timeout minutes"
    for i in `seq 1 $timeout`; do
      if [ -f WAITING ]; then
        memuse=$(free -m | grep ^Mem: | awk '{printf "%6d",$3}')
        memtot=$(free -m | grep ^Mem: | awk '{printf "%6d", $2}')
        swpuse=$(free -m | grep ^Swap: | awk '{printf "%6d", $3}')
        swptot=$(free -m | grep ^Swap: | awk '{printf "%6d", $2}')
        mempct=$(calc -d "printf('%5.3f\n',$memuse/$memtot)")
        swppct=$(calc -d "printf('%5.3f\n',$swpuse/$swptot)")
        memuse2=$(calc -d "printf('%6d',$memuse+$swpuse)")
        mempct2=$(calc -d "printf('%5.3f\n',$memuse2/$memtot)")
        echo `date`": INFO: MEM  USED $mempct ($memuse/$memtot)"
        echo `date`": INFO: SWAP USED $swppct ($swpuse/$swptot)"
        echo `date`": INFO: M+SW USED $mempct2 ($memuse2/$memtot)"
        if [[ $mempct2 > 0.95 ]]; then
          echo `date`": ERROR: M+SW USED > 0.95"
          break
        fi
        echo `date`": WAITING ($i/$timeout)"
        sleep 10
      fi
      for i2 in `seq 1 5`; do
	if [ -f WAITING ]; then
	  sleep 10
	fi
      done
    done
    echo `date`: WAITING NO LONGER
    if [ -f WAITING ]; then
      echo `date`: STOP WAITING
      killAndWait
      if ps $waiterPid > /dev/null; then kill -9 $waiterPid; fi
    fi
    echo a2 `pwd -P`
    cd ${tstDir}
    echo a3 `pwd -P`
  fi
  echo `date` finished jobs=${j}, thrs=${t}
done

