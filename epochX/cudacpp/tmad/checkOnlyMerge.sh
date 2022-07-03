#!/bin/bash

if [ "$1" == "" ] || [ "$2" != "" ]; then
  echo "Usage: $0 <logfile.txt>"
  exit 1
fi

logfileold=${1}
logfiletmp=${1}.TMP
###echo logfileold=${logfileold}
###echo logfiletmp=${logfiletmp}

if [ ! -f ${logfileold} ]; then echo "ERROR! ${logfileold} not found"; exit 1; fi
if [ ! -f ${logfiletmp} ]; then echo "ERROR! ${logfiletmp} not found"; exit 1; fi

logfilenew=${1}.NEW
\rm -f ${logfilenew}

cat $logfileold | awk -v tmp=$logfiletmp '
  {if ($2=="EXECUTE") {do {getline < tmp} while ($2!="EXECUTE"); print $0;
                       for (i=0; i<3; i++) {getline; getline < tmp; print $0}}
   else print $0}
' > $logfilenew

###diff $logfilenew $logfileold
\mv $logfilenew $logfileold
