#!/bin/bash

set +x # not verbose
set -e # fail on error

scrdir=$(cd $(dirname $0); pwd)
bckend=$(basename $(cd $scrdir; cd ..; pwd)) # cudacpp or alpaka
topdir=$(cd $scrdir; cd ../../..; pwd)

function usage()
{
  echo "Usage: $0 <processes [-eemumu][-ggtt][-ggttg][-ggttgg][-ggttggg]> [-makeclean] [-runclean]" > /dev/stderr
  exit 1
}

##########################################################################
# PART 0 - decode command line arguments
##########################################################################

eemumu=0
ggtt=0
ggttg=0
ggttgg=0
ggttggg=0

makeclean=0
runclean=0

while [ "$1" != "" ]; do
  if [ "$1" == "-eemumu" ]; then
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
  elif [ "$1" == "-makeclean" ]; then
    makeclean=1
    shift
  elif [ "$1" == "-runclean" ]; then
    runclean=1
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

# Create an input file that is appropriate for the specific process
function inputfile()
{
  if [ "${eemumu}" == "1" ]; then
    nevt=16384 # 16384 unweighted events require 524320 MEs
  elif [ "${ggtt}" == "1" ]; then 
    nevt=16384 # 16384 unweighted events require 524320 MEs
  #  dir=$topdir/epochX/${bckend}/gg_tt${suff}SubProcesses/P1_gg_ttx
  #elif [ "${ggttg}" == "1" ]; then 
  #  dir=$topdir/epochX/${bckend}/gg_ttg${suff}SubProcesses/P1_gg_ttxg
  #elif [ "${ggttgg}" == "1" ]; then 
  #  dir=$topdir/epochX/${bckend}/gg_ttgg${suff}SubProcesses/P1_gg_ttxgg
  #elif [ "${ggttggg}" == "1" ]; then 
  #  dir=$topdir/epochX/${bckend}/gg_ttggg${suff}SubProcesses/P1_gg_ttxggg
  else
    echo "ERROR! Unknown process" > /dev/stderr; usage
  fi
  tmp=$(mktemp)
  cat << EOF > ${tmp}
${nevt} 1 1 ! Number of events and max and min iterations
0.000001 ! Accuracy (ignored because max iterations = min iterations)
0 ! Grid Adjustment 0=none, 2=adjust (NB if = 0, ftn26 will still be used if present)
0 ! Suppress Amplitude 1=yes (i.e. use MadEvent single-diagram enhancement)
0 ! Helicity Sum/event 0=exact
1 ! Channel number for single-diagram enhancement multi-channel (IGNORED as suppress amplitude is 0?)
EOF
  echo ${tmp}
}

# Run madevent (or cmadevent or gmadevent, depending on $1) and parse its output
function runmadevent()
{
  if [ "$1" == "" ] || [ "$2" != "" ]; then echo "Usage: runmadevent <madevent executable>"; exit 1; fi
  tmpin=$(inputfile)
  if [ ! -f $tmpin ]; then echo "ERROR! Missing input file $tmpin"; exit 1; fi
  tmp=$(mktemp)
  if [ ! -f results.dat ]; then grepA=4; else grepA=12; fi
  $timecmd $1 < ${tmpin} > ${tmp}
  ###echo $tmp
  cat ${tmp} | grep --binary-files=text -B1 -A${grepA} Accumulated
  if [ "${1/cmadevent}" != "$1" ] || [ "${1/gmadevent}" != "$1" ]; then
    # Hack: use awk to convert Fortran's 0.42E-01 into 4.20e-02
    cat ${tmp} | grep --binary-files=text MERATIOS \
      | awk -v sep=" 1 - " '{i=index($0,sep); if(i>0){print substr($0,0,i-1) sep 0+substr($0,i+length(sep))} else print $0}' \
      | awk -v sep=" 1 + " '{i=index($0,sep); if(i>0){print substr($0,0,i-1) sep 0+substr($0,i+length(sep))} else print $0}' \
      | awk -v sep1=" AVG = " -v sep2=" +- " '{i1=index($0,sep1); i2=index($0,sep2); if(i1>0 && i2>0){print substr($0,0,i1-1) sep1 0+substr($0,i1+length(sep1),i2-i1) sep2 0+substr($0,i2+length(sep2))} else print $0}'
  fi
  cat ${tmp} | grep --binary-files=text COUNTERS
}

##########################################################################
# PART 1 - build and run madevent
##########################################################################

for suff in $suffs; do

  dir=$(showdir)
  if [ ! -d $dir ]; then echo "WARNING! Skip missing directory $dir"; continue; fi

  ##########################################################################
  # PART 1a - build
  ##########################################################################

  echo "Working directory: $dir"
  cd $dir
  if [ "${makeclean}" == "1" ]; then make cleanall; echo; fi
  make -j

  ##########################################################################
  # PART 1b - run madevent
  ##########################################################################

  # Disable OpenMP multithreading in Fortran
  ###export OMP_NUM_THREADS=1 # not needed in .mad directories (OpenMP MT disabled in the code)

  # Use the time command?
  ###timecmd=time
  timecmd=

  # First execution: this will create results.dat
  cd $dir
  if [ "${runclean}" == "1" ]; then \rm -f results.dat; fi
  if [ ! -f results.dat ]; then
    echo -e "\n*** EXECUTE MADEVENT (create results.dat) ***"
    runmadevent ./madevent
  fi
  # Second execution: this will read results.dat and create events.lhe 
  echo -e "\n*** EXECUTE MADEVENT (create events.lhe) ***"
  runmadevent ./madevent
  echo -e "\n*** EXECUTE CMADEVENT_CUDACPP (create events.lhe) ***"
  runmadevent ./cmadevent_cudacpp
  echo -e "\n*** EXECUTE GMADEVENT_CUDACPP (create events.lhe) ***"
  runmadevent ./gmadevent_cudacpp
  
done

