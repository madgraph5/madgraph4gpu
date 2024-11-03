#!/bin/bash
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (May 2022) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

scrdir=$(cd $(dirname $0); pwd)

host=$(hostname)
if [ "${host/juwels}" != "${host}" ]; then ${scrdir}/juwelspatch.sh; fi # workaround for #498

short=0
bsm=
flts=-dmf # "d m f" (alternative: -d_f i.e. "d f")
makeclean=
rmrdat=
add10x="+10x"
hip=

while [ "$1" != "" ]; do
  if [ "$1" == "-short" ]; then
    short=1 # all (possibly including bsm) but ggttggg
    shift
  elif [ "$1" == "-ggttggg" ]; then
    short=-1 # only ggttggg (implies no bsm!)
    shift
  elif [ "$1" == "-makeclean" ]; then
    makeclean=$1
    shift
  elif [ "$1" == "-no10x" ]; then
    add10x=""
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
    echo "Usage: $0 [-short|-ggttggg] [-bsmonly|-nobsm] [-makeclean] [-no10x] [-hip]"
    exit 1
  fi
done

started="STARTED  AT $(date)"

if [ "${bsm}" != "-bsmonly" ]; then
  if [ "$short" == "1" ]; then
    ${scrdir}/teeMadX.sh -eemumu -ggtt -ggttg -ggttgg -gqttq $flts $makeclean $rmrdat $add10x $hip
  elif [ "$short" == "-1" ]; then
    ${scrdir}/teeMadX.sh -ggttggg $flts $makeclean $rmrdat $add10x $hip
  else
    ${scrdir}/teeMadX.sh -eemumu -ggtt -ggttg -ggttgg -gqttq -ggttggg $flts $makeclean $rmrdat $add10x $hip
  fi
fi
status=$?
ended1="(SM tests)\nENDED(1) AT $(date) [Status=$status]"

if [ "${bsm}" != "-nobsm" ]; then
  if [ "$short" != "-1" ]; then
    ${scrdir}/teeMadX.sh -heftggbb -susyggtt -susyggt1t1 -smeftggtttt $flts $makeclean $rmrdat $add10x $hip
  fi
fi
status=$?
ended2="(BSM tests)\nENDED(1) AT $(date) [Status=$status]"

# Print out the number of "OK!"s in each log (expect 24)
echo
printf "\n%80s\n" |tr " " "#"
echo
echo -e "$started"
echo -e "$ended1"
echo -e "$ended2"
echo
for f in ${scrdir}/logs_*_mad/log_*; do echo $(cat $f | grep OK  | wc -l) $f; done # expect 24

# Print out the MEK channelid debugging output
echo
\grep MEK ${scrdir}/logs_*/* | sed "s|${scrdir}/logs_||" | sed 's|_mad.*DEBUG:||' | sort -u
