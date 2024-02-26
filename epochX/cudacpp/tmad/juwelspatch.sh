#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Jul 2022) for the MG5aMC CUDACPP plugin.

topdir=$(cd $(dirname $0)/..; pwd)

procs="ee_mumu gg_tt gg_ttg gg_ttgg gg_ttggg gq_ttq"
for proc in $procs; do
  cd $topdir/${proc}.mad/Source
  cat vector.inc | sed 's/16384/32/' > vector.inc.NEW
  \mv vector.inc.NEW vector.inc
  cat coupl.inc | sed 's/16384/32/' > coupl.inc.NEW
  \mv coupl.inc.NEW coupl.inc
done
