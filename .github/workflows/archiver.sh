#!/bin/bash
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Sep 2024) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

# Path to the top directory of madgraphgpu
# In the CI this would be simply $(pwd), but allow the script to be run also outside the CI
echo "Executing $0 $*"
topdir=$(cd $(dirname $0)/../..; pwd)

# Create a temporary directory and a VERSION.txt file
tmpdir=$(mktemp -d)
outdir=${tmpdir}/CUDACPP_OUTPUT
mkdir ${outdir}
cd ${topdir}/epochX/cudacpp/CODEGEN/PLUGIN/CUDACPP_SA_OUTPUT
outfile=${outdir}/VERSION.txt
touch ${outfile}
dateformat='%Y-%m-%d_%H:%M:%S UTC'
echo "cudacpp_version              = $(cat __init__.py | awk '/__version__/{print $3}' | sed 's/(//' | sed 's/)//' | sed 's/,/./g')" >> ${outfile}
echo "mg5_version_current          = $(cat ../../../../../MG5aMC/mg5amcnlo/VERSION | awk '/version =/{print $3}' | sed -r 's/[^0-9.]//g')" >> ${outfile}
echo "" >> ${outfile}
echo "mg5_version_minimal          = $(cat __init__.py | awk '/minimal_mg5amcnlo_version/{print $3}'  | sed 's/(//' | sed 's/)//' | sed 's/,/./g')" >> ${outfile}
echo "mg5_version_latest_validated = $(cat __init__.py | awk '/latest_validated_version/{print $3}'  | sed 's/(//' | sed 's/)//' | sed 's/,/./g')" >> ${outfile}
echo "" >> ${outfile}
echo "TARBALL DATE:  $(date -u +"${dateformat}")" >> ${outfile}
echo "" >> ${outfile}
TZ=UTC git --no-pager log -n1 --date=format-local:"${dateformat}" --pretty=format:'commit         %h%nAuthor:        %an%nAuthorDate:    %ad%nCommitter:     %cn%nCommitterDate: %cd%nMessage:       "%s"%n' >> ${outfile}
python -c 'print("="*132)'; cat ${outfile}; python -c 'print("="*132)'
cp ${outfile} ${topdir}
echo "VERSION.txt file available on ${topdir}/$(basename ${outfile})"

# Copy all relevant plugin files to the temporary directory
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
