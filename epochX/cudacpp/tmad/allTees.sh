#!/bin/bash

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
  else
    echo "Usage: $0 [-short|-ggttggg] [-makeclean]"
    exit 1
  fi
done

if [ "$short" == "1" ]; then
  ${scrdir}/teeMadX.sh -eemumu -ggtt -ggttg -ggttgg $flts $makeclean $rmrdat $add10x
elif [ "$short" == "-1" ]; then
  ${scrdir}/teeMadX.sh -ggttggg $flts $makeclean $rmrdat $add10x
else
  ${scrdir}/teeMadX.sh -eemumu -ggtt -ggttg -ggttgg -ggttggg $flts $makeclean $rmrdat $add10x
fi
