#!/bin/bash
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Sep 2024) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

set -e # fail on error

topdir=$(cd $(dirname $0)/../..; pwd)

# The remote repo used for tagging
REMOTE=origin

# The branch name used for tagging
BRANCH=${REMOTE}/install # temporary! it should be master eventually...

# Determine mg5_version (as in HEPToolInstaller.py)
fVERSION=${topdir}/MG5aMC/mg5amcnlo/VERSION 
if [ ! -f ${fVERSION} ]; then echo "ERROR! File ${fVERSION} not found"; exit 1; fi 
VERSION=$(python3 -c "import re; print(re.findall(r'version\s*=\s*([\.\d]*)','$(cat ${fVERSION} | grep ^version)')[0])")
###echo $VERSION

# Determine UTC date
DATE=$(date -u +on%y%m%d_at%H%M%S)
###echo $DATE

# Build the unique and running tag names
PREFIX=TEST_cudacpp
uniqTAG=${PREFIX}_for${VERSION}_${DATE}
runnTAG=${PREFIX}_for${VERSION}_latest
echo ${uniqTAG}
echo ${runnTAG}

# Create the unique (~permanent) tag
git tag ${uniqTAG} ${BRANCH} -m "Unique tag ${uniqTAG}"

# Create the running tag
# (use '-f' to replace any previously existing tag with the same name)
git tag -f ${runnTAG} ${uniqTAG} -m "Running tag ${runnTAG} (linked to unique tag ${uniqTAG})"

# Push the tags to the remote repository
# (use '-f' to replace any previously existing tag with the same name)
git push -f ${REMOTE} --tags
