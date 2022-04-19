#!/bin/bash

status=0

scrdir=$(cd $(dirname $0); pwd)

if [ "$1" == "" ] || [ "$2" != "" ]; then
  echo "Usage: $0 <dir>"
  exit 1 
fi
dir=$1

if [ ! -d ${dir} ]; then echo "ERROR! Directory ${dir} does not exist"; exit 1; fi
if [ ! -d ${dir}.BKP ]; then echo "ERROR! Directory ${dir}.BKP does not exist"; exit 1; fi

set -x

mv ${dir}.BKP ${dir}.BKP.tmp
mv ${dir} ${dir}.BKP
mv ${dir}.BKP.tmp ${dir}
