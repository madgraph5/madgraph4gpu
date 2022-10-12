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
          if (line==1){npar=$1; print $0}
          else if (line<=1+npar){split("AAA"$0,a,FS,seps); # see https://stackoverflow.com/a/40456663 (add AAA to keep leading whitespaces)
                                 a[NF+1]=sprintf("%"length($NF)-1"d.",0); # helicity is field $NF (but a[NF+1] because of AAA)
                                 printf("%s",seps[1]); # print only leading whitespaces (not AAA)
                                 for(i=2;i<=NF+1;i++) printf("%s%s",a[i],seps[i]); # NF+1 because of AAA
                                 print ""} # print newline
          else print $0}}' > ${outfile}
