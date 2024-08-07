#!/bin/sh

if [ "$1" == "" ] || [ "$2" != "" ]; then
  echo "Usage:   $0 <output.txt>"
  echo "Example: $0 pp_dy012j.mad/fortran/output.txt"
  exit 1
fi

outfile=$1
if [ ! -f $outfile ]; then
  echo "ERROR! File not found $outfile"
  exit 1
fi

cat $outfile | grep "__CUDACPP_DEBUG: GridPackCmd.launch finished" | sed 's/__CUDACPP_DEBUG: /----> /'
for msg in "PROGRAM TOTAL   " "Fortran Overhead" "Fortran MEs     "; do
  cat $outfile | grep "\[COUNTERS\]" | grep "${msg}" | sed 's/s for.*//' | sed 's/s$//' \
    | awk -vmsg="${msg}" -vttot=0 '{jtot=$NF; ttot += jtot}; END{print "[madevent COUNTERS]", msg, ttot}'
done
