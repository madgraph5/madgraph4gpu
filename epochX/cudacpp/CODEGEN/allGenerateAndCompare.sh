#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Oct 2023) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.

set -e # fail on error

sa=1
mad=1
if [ "$1" == "--sa" ]; then
  mad=0; shift
elif [ "$1" == "--mad" ]; then
  sa=0; shift
fi
if [ "$1" != "" ]; then echo "Usage: $0 [--sa|--mad]"; exit 1; fi

cd $(dirname $0)/..

if [ "$mad" == "1" ]; then
  ./CODEGEN/generateAndCompare.sh -q ee_mumu --mad
  ./CODEGEN/generateAndCompare.sh -q gg_tt --mad
  ./CODEGEN/generateAndCompare.sh -q gg_ttg --mad
  ./CODEGEN/generateAndCompare.sh -q gg_ttgg --mad
  ./CODEGEN/generateAndCompare.sh -q gg_ttggg --mad
  ./CODEGEN/generateAndCompare.sh -q gq_ttq --mad
  ./CODEGEN/generateAndCompare.sh -q gg_tt01g --mad
  ./CODEGEN/generateAndCompare.sh -q pp_tt012j --mad
fi

if [ "$sa" == "1" ]; then
  ./CODEGEN/generateAndCompare.sh -q ee_mumu
  ./CODEGEN/generateAndCompare.sh -q gg_tt
  ./CODEGEN/generateAndCompare.sh -q gg_ttg
  ./CODEGEN/generateAndCompare.sh -q gg_ttgg
  ./CODEGEN/generateAndCompare.sh -q gg_ttggg
  ./CODEGEN/generateAndCompare.sh -q gq_ttq
  ./CODEGEN/generateAndCompare.sh -q heft_gg_h
fi
