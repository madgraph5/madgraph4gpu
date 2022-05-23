#!/bin/bash

set +x # not verbose
set -e # fail on error

scrdir=$(cd $(dirname $0); pwd)
bckend=$(basename $(cd $scrdir; cd ..; pwd)) # cudacpp or alpaka
topdir=$(cd $scrdir; cd ../../..; pwd)

# Can events be generated? (i.e. has issue #14 been fixed?)
genevt=0 # present implementation (workaround for issue #14)
###genevt=1 # test if issue #14 has been fixed

function usage()
{
  if [ "$genevt" == "0" ]; then
    echo "Usage: $0 <processes [-eemumu][-ggtt][-ggttg][-ggttgg][-ggttggg]> [-d] [-makeonly|-makeclean|-makecleanonly] [-keeprdat]" > /dev/stderr
  else
    echo "Usage: $0 <processes [-eemumu][-ggtt][-ggttg][-ggttgg][-ggttggg]> [-d] [-makeonly|-makeclean|-makecleanonly]" > /dev/stderr
  fi
  exit 1
}

##########################################################################
# PART 0 - decode command line arguments
##########################################################################

debug=0

eemumu=0
ggtt=0
ggttg=0
ggttgg=0
ggttggg=0

maketype=
###makej=

if [ "$genevt" == "1" ]; then keeprdat=0; fi

while [ "$1" != "" ]; do
  if [ "$1" == "-d" ]; then
    debug=1
    shift
  elif [ "$1" == "-eemumu" ]; then
    eemumu=1
    shift
  elif [ "$1" == "-ggtt" ]; then
    ggtt=1
    shift
  elif [ "$1" == "-ggttg" ]; then
    ggttg=1
    shift
  elif [ "$1" == "-ggttgg" ]; then
    ggttgg=1
    shift
  elif [ "$1" == "-ggttggg" ]; then
    ggttggg=1
    shift
  elif [ "$1" == "-makeonly" ] || [ "$1" == "-makeclean" ] || [ "$1" == "-makecleanonly" ]; then
    if [ "${maketype}" != "" ] && [ "${maketype}" != "$1" ]; then
      echo "ERROR! Options -makeonly, -makeclean and -makecleanonly are incompatible"; usage
    fi
    maketype="$1"
    shift
  elif [ "$genevt" == "0" ] && [ "$1" == "-keeprdat" ]; then
    keeprdat=1
    shift
  else
    usage
  fi
done
###exit 1

# Check that at least one process has been selected
if [ "${eemumu}" == "0" ] && [ "${ggtt}" == "0" ] && [ "${ggttg}" == "0" ] && [ "${ggttgg}" == "0" ] && [ "${ggttggg}" == "0" ]; then usage; fi

# Always test only the .mad/ directories (hardcoded)
suffs=".mad/"

# Determine the working directory below topdir based on suff, bckend and <process>
function showdir()
{
  if [ "${suff}" == ".mad/" ]; then
    if [ "${eemumu}" == "1" ]; then 
      dir=$topdir/epochX/${bckend}/ee_mumu${suff}SubProcesses/P1_ll_ll
    elif [ "${ggtt}" == "1" ]; then 
      dir=$topdir/epochX/${bckend}/gg_tt${suff}SubProcesses/P1_gg_ttx
    elif [ "${ggttg}" == "1" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttg${suff}SubProcesses/P1_gg_ttxg
    elif [ "${ggttgg}" == "1" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttgg${suff}SubProcesses/P1_gg_ttxgg
    elif [ "${ggttggg}" == "1" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttggg${suff}SubProcesses/P1_gg_ttxggg
    fi
  else
    if [ "${eemumu}" == "1" ]; then 
      dir=$topdir/epochX/${bckend}/ee_mumu${suff}SubProcesses/P1_Sigma_sm_epem_mupmum
    elif [ "${ggtt}" == "1" ]; then 
      dir=$topdir/epochX/${bckend}/gg_tt${suff}SubProcesses/P1_Sigma_sm_gg_ttx
    elif [ "${ggttg}" == "1" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttg${suff}SubProcesses/P1_Sigma_sm_gg_ttxg
    elif [ "${ggttgg}" == "1" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttgg${suff}SubProcesses/P1_Sigma_sm_gg_ttxgg
    elif [ "${ggttggg}" == "1" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttggg${suff}SubProcesses/P1_Sigma_sm_gg_ttxggg
    fi
  fi
  echo $dir
}

# Determine the appropriate number of events for the specific process
function getnevt()
{
  if [ "${eemumu}" == "1" ]; then # computed xsec is zero (https://github.com/oliviermattelaer/mg5amc_test/issues/13)
    ###nevt=16384 # computes 524320 MEs in 7.4s
    nevt=2048 # computes 65558 MEs in 1.1s
  elif [ "${ggtt}" == "1" ]; then 
    nevt=16384 # computes 16416 MEs (writes to file 368 events) in 1.2s
  elif [ "${ggttg}" == "1" ]; then # event generation fails (https://github.com/oliviermattelaer/mg5amc_test/issues/14)
    nevt=4096 # computes 4128 MEs in 1.2s
  elif [ "${ggttgg}" == "1" ]; then # event generation fails (https://github.com/oliviermattelaer/mg5amc_test/issues/14)
    ###nevt=1024 # computes 1056 MEs in 2.6s
    nevt=512 # computes 544 MEs in 1.5s
  elif [ "${ggttggg}" == "1" ]; then # event generation fails (https://github.com/oliviermattelaer/mg5amc_test/issues/14)
    nevt=64 # computes 96 MEs in 5.3s
  else
    echo "ERROR! Unknown process" > /dev/stderr; usage
  fi
  echo $nevt
}

# Create an input file that is appropriate for the specific process
function getinputfile()
{
  nevt=$(getnevt)
  tmp=$(mktemp)
  if [ "$1" == "-fortran" ]; then
    mv ${tmp} ${tmp}_fortran
    tmp=${tmp}_fortran
  elif [ "$1" == "-cuda" ]; then
    mv ${tmp} ${tmp}_cuda
    tmp=${tmp}_cuda
    nloop=32768
    while [ $nloop -gt $nevt ]; do (( nloop = nloop / 2 )); done
    echo "${nloop} ! Number of events in a single CUDA iteration (nb_page_loop)" >> ${tmp}
  elif [ "$1" == "-cpp" ]; then
    mv ${tmp} ${tmp}_cuda
    tmp=${tmp}_cuda
    echo "32 ! Number of events in a single C++ iteration (nb_page_loop)" >> ${tmp}
  else
    echo "Usage: getinputfile <backend [-fortran][-cuda]-cpp]>"
    exit 1
  fi
  cat << EOF >> ${tmp}
${nevt} 1 1 ! Number of events and max and min iterations
0.000001 ! Accuracy (ignored because max iterations = min iterations)
0 ! Grid Adjustment 0=none, 2=adjust (NB if = 0, ftn26 will still be used if present)
0 ! Suppress Amplitude 1=yes (i.e. use MadEvent single-diagram enhancement)
0 ! Helicity Sum/event 0=exact
2 ! Channel number for single-diagram enhancement multi-channel (IGNORED as suppress amplitude is 0?)
EOF
  echo ${tmp}
}

# Run check.exe or gcheck.exe (depending on $1) and parse its output
function runcheck()
{
  if [ "$1" == "" ] || [ "$2" != "" ]; then echo "Usage: runcheck <check/gcheck executable>"; exit 1; fi
  if [ "${1/gcheck}" != "$1" ]; then
    txt=GCHECK
  elif [ "${1/check}" != "$1" ]; then
    txt=CHECK
  else
    echo "ERROR! Unknown check executable '$1'"; exit 1
  fi
  nevt=$(getnevt)
  nthr=32
  while [ $nthr -gt $nevt ]; do (( nthr = nthr / 2 )); done
  (( nblk = nevt/nthr )) # NB integer division
  (( nevt2 = nblk*nthr ))
  if [ "$nevt" != "$nevt2" ]; then echo "ERROR! nevt($nevt) != nevt2($nevt2)=nthr($nthr)*nblk($nblk)"; exit 1; fi
  pattern="Process|Workflow|EvtsPerSec\[MECalc"
  echo -e "\n*** EXECUTE $txt -p $nblk $nthr 1 --bridge ***"
  $1 -p $nblk $nthr 1 --bridge | egrep "(${pattern})"
  echo -e "\n*** EXECUTE $txt -p $nblk $nthr 1 ***"
  $1 -p $nblk $nthr 1 | egrep "(${pattern})"
}

# Run madevent (or cmadevent or gmadevent, depending on $1) and parse its output
function runmadevent()
{
  if [ "$1" == "" ] || [ "$2" != "" ]; then echo "Usage: runmadevent <madevent executable>"; exit 1; fi
  if [ "${1/cmadevent}" != "$1" ]; then
    tmpin=$(getinputfile -cpp)
  elif [ "${1/gmadevent}" != "$1" ]; then
    tmpin=$(getinputfile -cuda)
  else
    tmpin=$(getinputfile -fortran)
  fi
  if [ ! -f $tmpin ]; then echo "ERROR! Missing input file $tmpin"; exit 1; fi
  tmp=$(mktemp)
  set +e # do not fail on error
  if [ "${debug}" == "1" ]; then
    echo "--------------------"; cat ${tmpin}; echo "--------------------"
    echo "Executing '$timecmd $1 < ${tmpin} > ${tmp}'"
  fi
  $timecmd $1 < ${tmpin} > ${tmp}
  if [ "$?" != "0" ]; then echo "ERROR! '$timecmd $1 < ${tmpin} > ${tmp}' failed"; tail -10 $tmp; exit 1; fi
  xsec=$(cat ${tmp} | grep --binary-files=text 'Cross sec =' | awk '{print 0+$NF}')
  if [ "${xsec}" != "" ]; then
    echo " [XSECTION] Cross section = ${xsec}"
  else
    echo -e " [XSECTION] ERROR! No cross section in log file:\n   $tmp\n   ..."
    tail -10 $tmp
    exit 1
  fi
  mch=$(cat ${tmp} | grep --binary-files=text 'MULTI_CHANNEL =' | awk '{print $NF}')
  chid=$(cat ${tmp} | grep --binary-files=text 'CHANNEL_ID =' | awk '{print $NF}')
  echo " [XSECTION] MultiChannel = ${mch}"
  echo " [XSECTION] ChannelId = ${chid}"
  evtf=$(cat ${tmp} | grep --binary-files=text 'events.' | grep 'Found' | awk '{print $2}')
  evtw=$(cat ${tmp} | grep --binary-files=text 'events.' | grep 'Wrote' | awk '{print $2}')
  if [ "${evtf}" != "" ] && [ "${evtw}" != "" ]; then
    echo " [UNWEIGHT] Wrote ${evtw} events (found ${evtf} events)"  
  fi
  if [ "${1/cmadevent}" != "$1" ] || [ "${1/gmadevent}" != "$1" ]; then
    # Hack: use awk to convert Fortran's 0.42E-01 into 4.20e-02
    cat ${tmp} | grep --binary-files=text MERATIOS \
      | awk -v sep=" 1 - " '{i=index($0,sep); if(i>0){print substr($0,0,i-1) sep 0+substr($0,i+length(sep))} else print $0}' \
      | awk -v sep=" 1 + " '{i=index($0,sep); if(i>0){print substr($0,0,i-1) sep 0+substr($0,i+length(sep))} else print $0}' \
      | awk -v sep1=" AVG = " -v sep2=" +- " '{i1=index($0,sep1); i2=index($0,sep2); if(i1>0 && i2>0){print substr($0,0,i1-1) sep1 0+substr($0,i1+length(sep1),i2-i1) sep2 0+substr($0,i2+length(sep2))} else print $0}'
  fi
  cat ${tmp} | grep --binary-files=text COUNTERS
  set -e # fail on error
}

##########################################################################
# PART 1 - build madevent
##########################################################################

for suff in $suffs; do

  dir=$(showdir)
  if [ ! -d $dir ]; then echo "WARNING! Skip missing directory $dir"; continue; fi
  echo "Working directory: $dir"
  cd $dir

  if [ "${maketype}" == "-makeclean" ]; then make cleanall; echo; fi
  if [ "${maketype}" == "-makecleanonly" ]; then make cleanall; echo; continue; fi
  make -j

done

if [ "${maketype}" == "-makecleanonly" ]; then printf "\nMAKE CLEANALL COMPLETED\n"; exit 0; fi
if [ "${maketype}" == "-makeonly" ]; then printf "\nMAKE COMPLETED\n"; exit 0; fi

##########################################################################
# PART 2 - run madevent
##########################################################################

printf "\nDATE: $(date '+%Y-%m-%d_%H:%M:%S')\n\n"

for suff in $suffs; do

  dir=$(showdir)
  if [ ! -d $dir ]; then echo "WARNING! Skip missing directory $dir"; continue; fi
  echo "Working directory: $dir"
  cd $dir

  # Disable OpenMP multithreading in Fortran
  ###export OMP_NUM_THREADS=1 # not needed in .mad directories (OpenMP MT disabled in the code)

  # Use the time command?
  ###timecmd=time
  timecmd=

  if [ "$genevt" == "1" ]; then

    # FUTURE IMPLEMENTATION? Currently fails with https://github.com/oliviermattelaer/mg5amc_test/issues/14
    # First execution: compute xsec (create results.dat)
    cd $dir
    \rm -f ftn26
    if [ "${keeprdat}" == "0" ]; then \rm -f results.dat; fi
    if [ ! -f results.dat ]; then
      echo -e "\n*** EXECUTE MADEVENT (create results.dat) ***"
      runmadevent ./madevent
    fi
    # Second execution: compute xsec and generate events (read results.dat and create events.lhe)
    echo -e "\n*** EXECUTE MADEVENT (create events.lhe) ***"
    runmadevent ./madevent
    echo -e "\n*** EXECUTE CMADEVENT_CUDACPP (create events.lhe) ***"
    runmadevent ./cmadevent_cudacpp
    echo -e "\n*** EXECUTE GMADEVENT_CUDACPP (create events.lhe) ***"
    runmadevent ./gmadevent_cudacpp

  else

    # CURRENT (TEMPORARY?) IMPLEMENTATION! Work around https://github.com/oliviermattelaer/mg5amc_test/issues/14
    # First execution ONLY (no event generation): compute xsec (create results.dat)
    cd $dir
    \rm -f ftn26 results.dat
    echo -e "\n*** EXECUTE MADEVENT (create results.dat) ***"
    runmadevent ./madevent
    \rm -f ftn26 results.dat
    echo -e "\n*** EXECUTE CMADEVENT_CUDACPP (create results.dat) ***"
    runmadevent ./cmadevent_cudacpp
    runcheck ./check.exe
    \rm -f ftn26 results.dat
    echo -e "\n*** EXECUTE GMADEVENT_CUDACPP (create results.dat) ***"
    runmadevent ./gmadevent_cudacpp
    runcheck ./gcheck.exe

  fi

done
printf "\nTEST COMPLETED\n"
