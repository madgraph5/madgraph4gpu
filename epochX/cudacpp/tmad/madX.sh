#!/bin/bash

set +x # not verbose
set -e # fail on error

scrdir=$(cd $(dirname $0); pwd)
bckend=$(basename $(cd $scrdir; cd ..; pwd)) # cudacpp or alpaka
topdir=$(cd $scrdir; cd ../../..; pwd)

function usage()
{
  echo "Usage: $0 <processes [-eemumu][-ggtt][-ggttg][-ggttgg][-ggttggg]> [-d] [-makeonly|-makeclean|-makecleanonly] [-keeprdat]" > /dev/stderr
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

keeprdat=0

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
  elif [ "$1" == "-keeprdat" ]; then
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
  tmpdir=/tmp/$USER
  mkdir -p $tmpdir
  if [ "${eemumu}" == "1" ]; then 
    tmp=$tmpdir/input_eemumu
  elif [ "${ggtt}" == "1" ]; then 
    tmp=$tmpdir/input_ggtt
  elif [ "${ggttg}" == "1" ]; then 
    tmp=$tmpdir/input_ggttg
  elif [ "${ggttgg}" == "1" ]; then 
    tmp=$tmpdir/input_ggttgg
  elif [ "${ggttggg}" == "1" ]; then 
    tmp=$tmpdir/input_ggttggg
  else
    echo "ERROR! cannot determine input file name"; exit 1
  fi
  \rm -f ${tmp}; touch ${tmp}
  if [ "$1" == "-fortran" ]; then
    mv ${tmp} ${tmp}_fortran
    tmp=${tmp}_fortran
  elif [ "$1" == "-cuda" ]; then
    mv ${tmp} ${tmp}_cuda
    tmp=${tmp}_cuda
    echo "+1 ! Fortran bridge mode (CppOnly=1, FortranOnly=0, BothQuiet=-1, BothDebug=-2)" >> ${tmp}
    nloop=32768
    while [ $nloop -gt $nevt ]; do (( nloop = nloop / 2 )); done
    echo "${nloop} ! Number of events in a single CUDA iteration (nb_page_loop)" >> ${tmp}
  elif [ "$1" == "-cpp" ]; then
    mv ${tmp} ${tmp}_cpp
    tmp=${tmp}_cpp
    echo "+1 ! Fortran bridge mode (CppOnly=1, FortranOnly=0, BothQuiet=-1, BothDebug=-2)" >> ${tmp}
    echo "32 ! Number of events in a single C++ or CUDA iteration (nb_page_loop)" >> ${tmp}
  else
    echo "Usage: getinputfile <backend [-fortran][-cuda]-cpp]>"
    exit 1
  fi
  cat << EOF >> ${tmp}
${nevt} 1 1 ! Number of events and max and min iterations
0.000001 ! Accuracy (ignored because max iterations = min iterations)
0 ! Grid Adjustment 0=none, 2=adjust (NB if = 0, ftn26 will still be used if present)
1 ! Suppress Amplitude 1=yes (i.e. use MadEvent single-diagram enhancement)
0 ! Helicity Sum/event 0=exact
1 ! Channel number (1-N) for single-diagram enhancement multi-channel (NB used even if suppress amplitude is 0!)
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
  cmd=$1
  if [ "${cmd/cmadevent}" != "$cmd" ]; then
    tmpin=$(getinputfile -cpp)
    cmd=${cmd/.\//.\/build.${avx}_d_inl0_hrd0\/}
  elif [ "${cmd/gmadevent2}" != "$cmd" ]; then
    cmd=${cmd/gmadevent2/gmadevent} # hack: run cuda gmadevent with cpp input file
    cmd=${cmd/.\//.\/build.none_d_inl0_hrd0\/}
    tmpin=$(getinputfile -cpp)
  elif [ "${cmd/gmadevent}" != "$cmd" ]; then
    cmd=${cmd/.\//.\/build.none_d_inl0_hrd0\/}
    tmpin=$(getinputfile -cuda)
  else # assume this is madevent (do not check)
    tmpin=$(getinputfile -fortran)
  fi
  if [ ! -f $tmpin ]; then echo "ERROR! Missing input file $tmpin"; exit 1; fi
  tmp=${tmpin/input/output}
  \rm -f ${tmp}; touch ${tmp}  
  set +e # do not fail on error
  if [ "${debug}" == "1" ]; then
    echo "--------------------"; cat ${tmpin}; echo "--------------------"
    echo "Executing '$timecmd $cmd < ${tmpin} > ${tmp}'"
  fi
  $timecmd $cmd < ${tmpin} > ${tmp}
  if [ "$?" != "0" ]; then echo "ERROR! '$timecmd $cmd < ${tmpin} > ${tmp}' failed"; tail -10 $tmp; exit 1; fi
  fbm=$(cat ${tmp} | grep --binary-files=text 'FBRIDGE_MODE =' | awk '{print $NF}')
  nbp=$(cat ${tmp} | grep --binary-files=text 'NB_PAGE_LOOP =' | awk '{print $NF}')
  mch=$(cat ${tmp} | grep --binary-files=text 'MULTI_CHANNEL =' | awk '{print $NF}')
  conf=$(cat ${tmp} | grep --binary-files=text 'Running Configuration Number:' | awk '{print $NF}')
  chid=$(cat ${tmp} | grep --binary-files=text 'CHANNEL_ID =' | awk '{print $NF}')
  echo " [XSECTION] nb_page_loop = ${nbp}"
  echo " [XSECTION] MultiChannel = ${mch}"
  echo " [XSECTION] Configuration = ${conf}"
  echo " [XSECTION] ChannelId = ${chid}"
  xsec=$(cat ${tmp} | grep --binary-files=text 'Cross sec =' | awk '{print 0+$NF}')
  xsec2=$(cat ${tmp} | grep --binary-files=text 'Actual xsec' | awk '{print $NF}')
  if [ "${fbm}" != "" ]; then
    echo " [XSECTION] Cross section = ${xsec} [${xsec2}] fbridge_mode=${fbm}"
  elif [ "${xsec2}" != "" ]; then
    echo " [XSECTION] Cross section = ${xsec} [${xsec2}]"
  elif [ "${xsec}" != "" ]; then
    echo " [XSECTION] Cross section = ${xsec}"
  else
    echo -e " [XSECTION] ERROR! No cross section in log file:\n   $tmp\n   ..."
    tail -10 $tmp
    exit 1
  fi
  evtf=$(cat ${tmp} | grep --binary-files=text 'events.' | grep 'Found' | awk '{print $2}')
  evtw=$(cat ${tmp} | grep --binary-files=text 'events.' | grep 'Wrote' | awk '{print $2}')
  if [ "${evtf}" != "" ] && [ "${evtw}" != "" ]; then
    echo " [UNWEIGHT] Wrote ${evtw} events (found ${evtf} events)"  
  fi
  if [ "${cmd/cmadevent}" != "$cmd" ] || [ "${cmd/gmadevent}" != "$cmd" ]; then
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
  echo "Working directory (build): $dir"
  cd $dir

  if [ "${maketype}" == "-makeclean" ]; then make cleanall; echo; fi
  if [ "${maketype}" == "-makecleanonly" ]; then make cleanall; echo; continue; fi
  make -j avxall

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
  echo "Working directory (run): $dir"
  cd $dir

  # Disable OpenMP multithreading in Fortran
  ###export OMP_NUM_THREADS=1 # not needed in .mad directories (OpenMP MT disabled in the code)

  # Use the time command?
  ###timecmd=time
  timecmd=

  # Show results.dat?
  ###rdatcmd="stat results.dat"
  rdatcmd="echo"

  # DEFAULT IMPLEMENTATION : compute cross section and then generate events
  cd $dir

  # (1) MADEVENT
  \rm -f results.dat # ALWAYS remove results.dat before the first madevent execution
  if [ ! -f results.dat ]; then
    echo -e "\n*** (1) EXECUTE MADEVENT (create results.dat) ***"
    \rm -f ftn26
    runmadevent ./madevent
    \cp -p results.dat results.dat.ref
  fi
  echo -e "\n*** (1) EXECUTE MADEVENT (create events.lhe) ***"
  ${rdatcmd} | grep Modify | sed 's/Modify/results.dat /'
  \rm -f ftn26
  runmadevent ./madevent
  ${scrdir}/dummyColor.sh events.lhe events.lhe.ref
  ${scrdir}/dummyHelicities.sh events.lhe.ref events.lhe.ref2
  \mv events.lhe.ref2 events.lhe.ref
  
  # (2) CMADEVENT_CUDACPP
  for avx in none sse4 avx2 512y 512z; do
    if [ "${keeprdat}" == "1" ]; then \cp -p results.dat.ref results.dat; else \rm -f results.dat; fi  
    if [ ! -f results.dat ]; then
      echo -e "\n*** (2-$avx) EXECUTE CMADEVENT_CUDACPP (create results.dat) ***"
      \rm -f ftn26
      runmadevent ./cmadevent_cudacpp
    fi
    echo -e "\n*** (2-$avx) EXECUTE CMADEVENT_CUDACPP (create events.lhe) ***"
    ${rdatcmd} | grep Modify | sed 's/Modify/results.dat /'
    \rm -f ftn26
    runmadevent ./cmadevent_cudacpp
    runcheck ./check.exe
    \mv events.lhe events.lhe.cpp
    echo -e "\n*** (2-$avx) Compare CMADEVENT_CUDACPP events.lhe to MADEVENT events.lhe reference (with dummy colors and helicities) ***"
    if ! diff events.lhe.cpp events.lhe.ref; then echo "ERROR! events.lhe.cpp and events.lhe.ref differ!"; exit 1; else echo -e "\nOK! events.lhe.cpp and events.lhe.ref are identical"; fi
  done

  # (3) GMADEVENT_CUDACPP
  if [ "${keeprdat}" == "1" ]; then \cp -p results.dat.ref results.dat; else \rm -f results.dat; fi  
  if [ ! -f results.dat ]; then
    echo -e "\n*** (3) EXECUTE GMADEVENT_CUDACPP (create results.dat) ***"
    \rm -f ftn26
    runmadevent ./gmadevent2_cudacpp # hack: run cuda gmadevent with cpp input file
  fi
  echo -e "\n*** (3) EXECUTE GMADEVENT_CUDACPP (create events.lhe) ***"
  ${rdatcmd} | grep Modify | sed 's/Modify/results.dat /'
  \rm -f ftn26
  runmadevent ./gmadevent2_cudacpp # hack: run cuda gmadevent with cpp input file
  runcheck ./gcheck.exe
  \mv events.lhe events.lhe.cuda
  echo -e "\n*** (3) Compare GMADEVENT_CUDACPP events.lhe to MADEVENT events.lhe reference (with dummy colors and helicities) ***"
  if ! diff events.lhe.cuda events.lhe.ref; then echo "ERROR! events.lhe.cuda and events.lhe.ref differ!"; exit 1; else echo -e "\nOK! events.lhe.cuda and events.lhe.ref are identical"; fi

  # (3bis) GMADEVENT_CUDACPP
  if [ "${keeprdat}" == "1" ]; then \cp -p results.dat.ref results.dat; else \rm -f results.dat; fi  
  if [ ! -f results.dat ]; then
    echo -e "\n*** (3bis) EXECUTE GMADEVENT_CUDACPP (create results.dat) ***"
    \rm -f ftn26
    runmadevent ./gmadevent_cudacpp
  fi
  echo -e "\n*** (3bis) EXECUTE GMADEVENT_CUDACPP (create events.lhe) ***"
  ${rdatcmd} | grep Modify | sed 's/Modify/results.dat /'
  \rm -f ftn26
  runmadevent ./gmadevent_cudacpp
  runcheck ./gcheck.exe

  # Cleanup
  \rm results.dat
  \rm results.dat.ref
  \rm events.lhe
  \rm events.lhe.*

done
printf "\nTEST COMPLETED\n"
