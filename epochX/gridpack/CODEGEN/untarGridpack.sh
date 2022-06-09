#!/bin/bash

status=0

scrdir=$(cd $(dirname $0); pwd)

if [ "${1%.tar.gz}" == "$1" ] || [ "$2" != "" ]; then
  echo "Usage: $0 <gridpack.tar.gz>"
  exit 1 
fi
tar=$1
dir=${1%.tar.gz}
###echo "Input file: $tar"
###echo "Output dir: $dir"

if [ ! -e ${tar} ]; then echo "ERROR! File $tar does not exist"; exit 1; fi
if [ -e ${dir} ] && ! rm -rf ${dir}; then echo "ERROR! $dir exists and could not be removed"; exit 1; fi
if ! mkdir -p ${dir}; then echo "ERROR! Directory $dir could not be created"; exit 1; fi

# Copy the tarfile to the local directory, untar it and remove it
cd ${dir}
echo "Untarring gridpack in $(pwd)"
cp ${tar} .
tar -xzf $(basename ${tar})
rm -f $(basename ${tar})

# Restore data (see https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/GridDevelopment)
pushd madevent > /dev/null
./bin/internal/restore_data default > /dev/null
popd > /dev/null

# Basic cleanup
pushd madevent/Source > /dev/null
make clean > /dev/null
popd > /dev/null
\rm madevent/bin/madevent
\rm madevent/lib/*.a
\rm madevent/Source/BIAS/dummy/*.o
\rm $(find . -name *.pyc)
touch madevent/Events/.keepme

# Copy or replace some custom scripts
\cp -dpr ${scrdir}/MG5aMC_patches/rebuild.sh .
\cp -dpr ${scrdir}/MG5aMC_patches/throughputGridpack.sh .

# Fix runtime errors from run.sh
cat madevent/Cards/me5_configuration.txt | sed 's/mg5_path/#mg5_path/' > madevent/Cards/me5_configuration.txt.new
\mv madevent/Cards/me5_configuration.txt.new madevent/Cards/me5_configuration.txt

# Inject C++ counters into the Fortran code
###exit 0
for p1dir in madevent/SubProcesses/P1_*; do
  cd $p1dir
  echo "Applying patches in $(pwd)"  
  \cp -dpr ${scrdir}/MG5aMC_patches/timer.h .
  \cp -dpr ${scrdir}/MG5aMC_patches/counters.cpp .
  if ! patch -i ${scrdir}/MG5aMC_patches/patch.driver.f; then status=1; fi
  if ! patch -i ${scrdir}/MG5aMC_patches/patch.matrix1_optim.f; then status=1; fi
  \rm -f matrix1_optim.f.orig
  cd - > /dev/null
done
cd madevent/SubProcesses
if ! patch -i ${scrdir}/MG5aMC_patches/patch.makefile; then status=1; fi
cd - > /dev/null

# Replace "-O" by "-O3 -ffast-math" globally
cat madevent/Source/make_opts | sed "s/GLOBAL_FLAG=-O /GLOBAL_FLAG=-O3 -ffast-math /" > madevent/Source/make_opts.new
\mv madevent/Source/make_opts.new madevent/Source/make_opts

# Dump the final contents of the local directory
echo "In $(pwd):"
ls -l .

exit $status
