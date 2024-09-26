#!/bin/bash
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Sep 2024) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

# Path to the top directory of madgraphgpu
# In the CI this would be simply $(pwd), but allow the script to be run also outside the CI
echo "Executing $0 $*"
topdir=$(cd $(dirname $0)/../..; pwd)

# Copy all relevant plugin files to a temporary directory
tmpdir=$(mktemp -d)
outdir=${tmpdir}/CUDACPP_OUTPUT
mkdir ${outdir}
cd ${topdir}/epochX/cudacpp/CODEGEN/PLUGIN/CUDACPP_SA_OUTPUT
for file in $(git ls-tree --name-only HEAD -r); do
  if [ "${file/acceptance_tests}" != "${file}" ]; then continue; fi # acceptance_tests are not needed for code generation
  mkdir -p ${outdir}/$(dirname ${file})
  cp -dp ${file} ${outdir}/${file} # preserve symlinks for AUTHORS, COPYING, COPYING.LESSER and COPYRIGHT
done

# Create the tgz archive
outtgz=cudacpp.tar.gz
cd ${tmpdir}
tar -czf ${outtgz} CUDACPP_OUTPUT
mv ${outtgz} ${topdir}
echo "Archive available on ${topdir}/${outtgz}"
