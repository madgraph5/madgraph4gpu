#!/bin/bash

# Temmporary workaround for #402
# Replace event-by-event color information in <events.lhe.in> by a dummy color "0" (integer) in <events.lhe.out>
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
                                 a[6]=sprintf("%"length($5)"d",0); # 1st color is field 5 (but a[6] because of AAA)
                                 a[7]=sprintf("%"length($6)"d",0); # 2nd color is field 6 (but a[7] because of AAA)
                                 printf("%s",seps[1]); # print only leading whitespaces (not AAA)
                                 for(i=2;i<=NF+1;i++) printf("%s%s",a[i],seps[i]); # NF+1 because of AAA
                                 print ""} # print newline
          else print $0}}' > ${outfile}
