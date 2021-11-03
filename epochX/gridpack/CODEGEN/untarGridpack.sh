#!/bin/bash

if [ "${1%.tar.gz}" == "$1" ] || [ "$2" == "" ] || [ "$3" != "" ]; then
  echo "Usage: $0 <input gridpack.tar.gz> <output directory>"
  exit 1 
fi
tar=$1
dir=$2
###echo "Input file: $tar"
###echo "Output dir: $dir"

if [ ! -e ${tar} ]; then echo "ERROR! File $tar does not exist"; exit 1; fi
if [ -d ${dir} ]; then echo "ERROR! Directory $dir already exists"; exit 1; fi
if [ -e ${dir} ]; then echo "ERROR! File $dir already exists"; exit 1; fi
if ! mkdir -p ${dir}; then echo "ERROR! Directory $dir could not be created"; exit 1; fi

cd ${dir}

# Copy the tarfile to the local directory, untar it and remove it
cp ${tar} .
tar -xzf $(basename ${tar})
rm -f $(basename ${tar})

# Restore data (see https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/GridDevelopment)
pushd madevent >& /dev/null
./bin/internal/restore_data default >& /dev/null
popd >& /dev/null

# Basic cleanup
pushd madevent/Source >& /dev/null
make clean
popd >& /dev/null
\rm madevent/bin/madevent
\rm madevent/lib/*.a
\rm madevent/Source/BIAS/dummy/*.o
\rm $(find . -name *.pyc)
touch madevent/Events/.keepme

# Dump the final contents of the local directory
echo "In $(pwd):"
ls -l .
