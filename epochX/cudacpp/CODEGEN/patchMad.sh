#!/bin/bash

status=0

scrdir=$(cd $(dirname $0); pwd)

if [ "${1%.madonly}" == "$1" ] || [ "$2" != "" ]; then
  echo "Usage: $0 <process.madonly>"
  exit 1 
fi
dir=$1
###echo "Current dir: $pwd"
###echo "Input dir to patch: $dir"

if [ ! -e ${dir} ]; then echo "ERROR! Directory $dir does not exist"; exit 1; fi

# These two steps are part of "cd Source; make" but they actually are code-generating steps
${dir}/bin/madevent treatcards run
${dir}/bin/madevent treatcards param

# Cleanup
\rm -f ${dir}/crossx.html
\rm -f ${dir}/index.html
\rm -f ${dir}/madevent.tar.gz
\rm -rf ${dir}/bin/internal/__pycache__
\rm -rf ${dir}/bin/internal/ufomodel/__pycache__
echo -e "index.html" > ${dir}/.gitignore
touch ${dir}/Events/.keepme
