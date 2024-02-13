#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (May 2022) for the MG5aMC CUDACPP plugin.

scrdir=$(cd $(dirname $0); pwd)

host=$(hostname)
if [ "${host/juwels}" != "${host}" ]; then ${scrdir}/juwelspatch.sh; fi # workaround for #498

short=0
flts=-mix # "d f m" (alternative: -flt i.e. "d f")
makeclean=
rmrdat=
add10x="+10x"

while [ "$1" != "" ]; do
  if [ "$1" == "-short" ]; then
    short=1 # all but ggttggg
    shift
  elif [ "$1" == "-ggttggg" ]; then
    short=-1 # only ggttggg
    shift
  elif [ "$1" == "-makeclean" ]; then
    makeclean=$1
    shift
  elif [ "$1" == "-no10x" ]; then
    add10x=""
    shift
  else
    echo "Usage: $0 [-short|-ggttggg] [-makeclean] [-no10x]"
    exit 1
  fi
done

if [ "$short" == "1" ]; then
  ${scrdir}/teeMadX.sh -eemumu -ggtt -ggttg -ggttgg -gqttq $flts $makeclean $rmrdat $add10x
elif [ "$short" == "-1" ]; then
  ${scrdir}/teeMadX.sh -ggttggg $flts $makeclean $rmrdat $add10x
else
  ${scrdir}/teeMadX.sh -eemumu -ggtt -ggttg -ggttgg -gqttq -ggttggg $flts $makeclean $rmrdat $add10x
fi


# Print out the number of "OK!"s in each log (expect 24)
echo
for f in ${scrdir}/logs_*_mad/log_*; do echo $(cat $f | grep OK  | wc -l) $f; done # expect 24
