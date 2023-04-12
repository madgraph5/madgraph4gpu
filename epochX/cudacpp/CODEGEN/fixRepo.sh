#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Apr 2023) for the MG5aMC CUDACPP plugin.

SCRDIR=$(cd $(dirname $0); pwd)
TOPDIR=$(dirname $SCRDIR) # e.g. epochX/cudacpp if $SCRDIR=epochX/cudacpp/CODEGEN

if [ "$1" == "" ]; then
  echo "Usage: $0 <procdir1> [<procdir2>...]"
  exit 1
fi
procdirs=$@

for procdir in $procdirs; do
  if [ ! -d $TOPDIR/$procdir ]; then echo "ERROR! Directory not found $TOPDIR/$procdir"; exit 1; fi
done

# Add to the git repo the copyright files in <procdir>/src if they exist
for procdir in $procdirs; do
  echo "=== $TOPDIR/$procdir"
  cd $TOPDIR/$procdir
  for fcopy in COPYRIGHT COPYING COPYING.LESSER AUTHORS; do
    if [ -f src/${fcopy} ]; then git add src/${fcopy}; fi
  done
done

