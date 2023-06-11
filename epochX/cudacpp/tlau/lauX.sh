#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Jun 2023) for the MG5aMC CUDACPP plugin.

cd $(dirname $0)/..
proc=gg_tt.mad
echo "Execute $(basename $0) for process ${proc} in directory $(pwd)"
procdir=$(pwd)/gg_tt.mad
cd $procdir

# Cleanup before launch
rm -rf Events HTML; mkdir Events HTML; touch Events/.keep HTML/.keep
rm -f SubProcesses/ME5_debug
cp SubProcesses/randinit SubProcesses/randinit.BKP # save the initial randinit

# Launch (generate_events)
set -x # verbose
MG5AMC_CARD_PATH=$(pwd)/Cards ./bin/generate_events -f
set +x # not verbose

# Cleanup after launch
mv SubProcesses/randinit.BKP SubProcesses/randinit # restore the initial randinit
