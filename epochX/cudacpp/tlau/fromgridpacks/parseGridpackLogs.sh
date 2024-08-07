#!/bin/sh

if [ "$1" == "" ] || [ "$2" != "" ]; then
  echo "Usage:   $0 <proc.mad>"
  echo "Example: $0 pp_dy012j.mad"
  exit 1
fi

procdir=$1
if [ ! -d $procdir ]; then
  echo "ERROR! Directory not found $procdir"
  exit 1
fi

for outfile in $procdir/*/output.txt; do
  echo $outfile
  cat $outfile | grep "__CUDACPP_DEBUG: GridPackCmd.launch finished" | sed 's/__CUDACPP_DEBUG: /----> /'
  for msg in "PROGRAM TOTAL   " "Fortran Overhead" "Fortran MEs     "; do
    cat $outfile | grep "\[COUNTERS\]" | grep "${msg}" | sed 's/s for.*//' | sed 's/s$//' \
      | awk -vmsg="${msg}" -vttot=0 '{jtot=$NF; ttot += jtot}; END{print "[madevent COUNTERS]", msg, ttot}'
  done
  echo "--------------------------------------------------------------------------------"
done
