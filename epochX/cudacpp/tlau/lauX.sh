#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Jun 2023) for the MG5aMC CUDACPP plugin.

cd $(dirname $0)/..
proc=gg_tt.mad
echo "Execute $(basename $0) for process ${proc} in directory $(pwd)"

set -x # verbose

rm -rf $proc; git checkout $proc; cd $proc
MG5AMC_CARD_PATH=$(pwd)/Cards ./bin/generate_events -f

