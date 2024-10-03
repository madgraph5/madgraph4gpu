#!/bin/bash
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Sep 2024) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

# Path to the top directory of madgraphgpu
# In the CI this would be simply $(pwd), but allow the script to be run also outside the CI
echo "Executing $0 $*"
topdir=$(cd $(dirname $0)/../..; pwd)

# Check that all git submodules have been updated
cd ${topdir}
if ! git submodule status | grep '^ ' > /dev/null; then
  echo "ERROR! There are git submodules that need to be updated"
  git submodule status
  exit 1
fi
mg5_commit_current=$(git submodule status | awk '/ MG5aMC\/mg5amcnlo /{print substr($1,0,7)}')

# Create a temporary directory and a VERSION.txt file
cd ${topdir}/epochX/cudacpp/CODEGEN/PLUGIN/CUDACPP_SA_OUTPUT
tmpdir=$(mktemp -d)
outdir=${tmpdir}/CUDACPP_OUTPUT
mkdir ${outdir}
outfile=${outdir}/VERSION.txt
touch ${outfile}
dateformat='%Y-%m-%d_%H:%M:%S UTC'
echo "(From CUDACPP_OUTPUT/__init__.py)" >> ${outfile}
echo "cudacpp_version              = $(cat __init__.py | awk '/__version__/{print $3}' | sed 's/(//' | sed 's/)//' | sed 's/,/./g')" >> ${outfile}
echo "mg5_version_minimal          = $(cat __init__.py | awk '/minimal_mg5amcnlo_version/{print $3}'  | sed 's/(//' | sed 's/)//' | sed 's/,/./g')" >> ${outfile}
echo "mg5_version_latest_validated = $(cat __init__.py | awk '/latest_validated_version/{print $3}'  | sed 's/(//' | sed 's/)//' | sed 's/,/./g')" >> ${outfile}
echo "" >> ${outfile}
echo "(From MG5AMC/mg5amcnlo)" >> ${outfile}
echo "mg5_version_current          = $(cat ../../../../../MG5aMC/mg5amcnlo/VERSION | awk '/version =/{print $3}' | sed -r 's/[^0-9.]//g')" >> ${outfile}
echo "mg5_commit_current           = ${mg5_commit_current}" >> ${outfile}
echo "" >> ${outfile}
echo "TARBALL DATE:  $(date -u +"${dateformat}")" >> ${outfile}
echo "" >> ${outfile}
TZ=UTC git --no-pager log -n1 --date=format-local:"${dateformat}" --pretty=format:'commit         %h%nAuthor:        %an%nAuthorDate:    %ad%nCommitter:     %cn%nCommitterDate: %cd%nMessage:       "%s"%n' >> ${outfile}
python3 -c 'print("="*132)'; cat ${outfile}; python3 -c 'print("="*132)'
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
