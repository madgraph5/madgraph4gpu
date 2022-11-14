#!/bin/bash

scrdir=$(cd $(dirname $0); pwd)
topdir=${scrdir%/*/*/*} # The directory with this script is 3 levels below the git topdir
###echo $scrdir
###echo $topdir

reldir=${PWD#$topdir/} # Relative path from the git topdir to the current directory
###echo $reldir

function usage()
{
  echo "Usage:   $0 <path-to-file> [-n <nold max (default=10)>] [-t <text to grep>] [-log]"
  echo "Example: $0 tput/logs_eemumu_auto/log_eemumu_auto_d_inl0.txt -t '(EvtsPerSec\[MECalcOnly\]|nvcc|registers)' -log -n 5"
  exit 1
}

if [ "$1" == "" ]; then usage; fi
file=$1; shift

nold=10
text=
log=
while [ "$1" != "" ]; do  
  if [ "$1" == "-n" ] && [ "$2" != "" ]; then
    nold=$2; shift; shift
  elif [ "$1" == "-t" ] && [ "$2" != "" ]; then
    text=$2; shift; shift
  elif [ "$1" == "-log" ]; then
    log=1; shift
  else
    usage
  fi
done
###echo file: $file
###echo nold: $nold
###echo text: $text

if [ ! -f "$file" ]; then
  echo "ERROR! File not found: $file"
  exit 1
fi

revs=$(git log  --oneline $file | cut -d' ' -f1)
iold=0
for rev in $revs; do
  if [[ $iold -ge $nold ]]; then break; fi
  iold=$((iold+1))
  if [ "$log" == "1" ]; then
    git log -1 $rev; echo
  else
    echo REVISION=$rev
  fi
  if [ "$text" != "" ]; then
    git show $rev:$reldir/$file | egrep "$text"
  else
    git show $rev:$reldir/$file | more -1000
  fi
  echo
done

