#!/bin/bash

scrdir=$(cd $(dirname $0); pwd)

host=$(hostname)
if [ "${host/juwels}" != "${host}" ]; then ${scrdir}/juwelspatch.sh; fi # workaround for #498

short=0
makeclean=
rmrdat=
add10x="+10x +100x"

while [ "$1" != "" ]; do
  if [ "$1" == "-short" ]; then
    short=1
    shift
  elif [ "$1" == "-makeclean" ]; then
    makeclean=$1
    shift
  elif [ "$1" == "-no100x" ]; then
    add10x="+10x"
    shift
  else
    echo "Usage: $0 [-short] [-makeclean] [-no100x]"
    exit 1
  fi
done

if [ "$short" == "1" ]; then
  ${scrdir}/teeMadX.sh -eemumu -ggtt -ggttg -ggttgg $makeclean $rmrdat $add10x
else
  ${scrdir}/teeMadX.sh -eemumu -ggtt -ggttg -ggttgg -ggttggg $makeclean $rmrdat $add10x
fi
