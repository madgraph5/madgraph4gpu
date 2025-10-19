#!/bin/bash
# Copyright (C) 2020-2025 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (May 2022) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2022-2025) for the MG5aMC CUDACPP plugin.

scrdir=$(cd $(dirname $0); pwd)

host=$(hostname)
if [ "${host/juwels}" != "${host}" ]; then ${scrdir}/juwelspatch.sh; fi # workaround for #498

# Usage
function usage()
{
  echo "Usage (1): $0 [-short|-ggttggg] [-bsmonly|-nobsm] [-makeclean] [+10x] [-hip]"
  echo "Run tests and check all logs"
  echo ""
  echo "Usage (2): $0 -checkonly"
  echo "Check existing logs without running any tests"
  exit 1
}

# Parse command line arguments
checkonly=0
short=0
bsm=
flts=-dmf # "d m f" (alternative: -d_f i.e. "d f")
makeclean=-makeclean
rmrdat=
add10x=
hip=
if [ "$1" == "-checkonly" ]; then
  # Check existing logs without running any tests?
  checkonly=1
  shift
  if [ "$1" != "" ]; then usage; fi
fi
while [ "${checkonly}" == "0" ] && [ "$1" != "" ]; do
  if [ "$1" == "-short" ]; then
    short=1 # all (possibly including bsm) but ggttggg
    shift
  elif [ "$1" == "-ggttggg" ]; then
    short=-1 # only ggttggg (implies no bsm!)
    shift
  elif [ "$1" == "-nomakeclean" ]; then
    makeclean=$1
    shift
  elif [ "$1" == "+10x" ]; then
    add10x=$1
    shift
  elif [ "$1" == "-bsmonly" ] && [ "$bsm" != "-nobsm" ]; then
    bsm=$1
    shift
  elif [ "$1" == "-nobsm" ] && [ "$bsm" != "-bsmonly" ]; then
    bsm=$1
    shift
  elif [ "$1" == "-hip" ]; then
    hip=$1
    shift
  else
    usage
  fi
done

# Run all tests
if [ "${checkonly}" == "0" ]; then
  started="STARTED  AT $(date)"
  # SM tests
  cmd1="(none)"
  if [ "${bsm}" != "-bsmonly" ]; then
    if [ "$short" == "1" ]; then
      cmd1="${scrdir}/teeMadX.sh -eemumu -ggtt -ggttg -ggttgg -gqttq $flts $makeclean $rmrdat $add10x $hip"
    elif [ "$short" == "-1" ]; then
      cmd1="${scrdir}/teeMadX.sh -ggttggg $flts $makeclean $rmrdat $add10x $hip"
    else
      cmd1="${scrdir}/teeMadX.sh -eemumu -ggtt -ggttg -ggttgg -gqttq -ggttggg $flts $makeclean $rmrdat $add10x $hip"
    fi
    ${cmd1}
  fi
  status=$?
  ended1="(SM tests)\nENDED(1) AT $(date) [Status=$status]"
  # BSM tests
  cmd2="(none)"
  if [ "${bsm}" != "-nobsm" ]; then
    if [ "$short" != "-1" ]; then
      cmd2="${scrdir}/teeMadX.sh -heftggbb -susyggtt -susyggt1t1 -smeftggtttt $flts $makeclean $rmrdat $add10x $hip"
    fi
    ${cmd2}
  fi
  status=$?
  ended2="(BSM tests)\nENDED(1) AT $(date) [Status=$status]"
  # Timing information
  echo
  printf "\n%80s\n" |tr " " "#"
  echo
  echo -e "$started"
  echo -e "$cmd1"
  echo -e "$ended1"
  echo -e "$cmd2"
  echo -e "$ended2"
  echo
fi

# Print out the number of "OK!"s in each log (expect 24)
for f in ${scrdir}/logs_*_mad/log_*; do echo $(cat $f | grep OK  | wc -l) $f; done # expect 24

# Print out any errors or aborts in the logs
echo
txt=$(egrep -i '(error|abort)' tmad/logs* -r | sed 's/:0:rocdevice.cpp.*Aborting.*/rocdevice.cpp: Aborting/' | sed "s/Gpu.*Assert/Assert/")
if [ "${txt}" == "" ]; then
  echo "No errors or aborts found in logs"
else
  echo "${txt}"
fi

# Print out any asserts in the logs
echo
txt=$(grep assert tmad/logs* -r | sed "s/Gpu.*Assert/Assert/")
if [ "${txt}" == "" ]; then
  echo "No asserts found in logs"
else
  echo "${txt}"
fi
  
# Print out any segfaults in the logs
echo
txt=$(grep -i segmentation tmad/logs* -r | sed "s/Gpu.*Assert/Assert/")
if [ "${txt}" == "" ]; then
  echo "No segmentation fault found in logs"
else
  echo "${txt}"
fi

# Print out the MEK channelid debugging output
echo
\grep MEK ${scrdir}/logs_*/* | sed "s|${scrdir}/logs_||" | sed 's|_mad.*DEBUG:||' | sort -u
