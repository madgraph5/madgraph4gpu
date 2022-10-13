#!/bin/bash

# Temmporary workaround for #403
# Replace event-by-event helicity information in <events.lhe.in> by a dummy helicity "0." (floating point) in <events.lhe.out>
if [ "$1" == "$2" ] || [ "$2" == "" ] || [ "$3" != "" ]; then
  echo "Usage: $0 <events.lhe.in> <events.lhe.out>"
  exit 1
fi
infile=$1
outfile=$2

if [ ! -f ${infile} ]; then
  echo "ERROR! File not found: ${infile}"
  exit 1
fi

cat ${infile} | awk \
  'BEGIN{line=1000; npar=0}
   {if ($1=="<event>"){line=0; print $0}
    else {line=line+1;
          if (line==1){$3=sprintf("%.2E",$3); print $0} # keep 3 digits "1.23" in XWGTUP (first one from sprintf is non-0)
          else print $0}}' > ${outfile}
