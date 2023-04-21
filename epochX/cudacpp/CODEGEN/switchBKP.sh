#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Apr 2022) for the MG5aMC CUDACPP plugin.

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
