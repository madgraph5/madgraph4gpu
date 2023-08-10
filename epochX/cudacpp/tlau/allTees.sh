#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (June 2023) for the MG5aMC CUDACPP plugin.

for proc in gg_tt gg_ttg gg_ttgg gg_ttggg; do
  ./tlau/lauX.sh -CUDA ${proc}
  ./tlau/lauX.sh -FORTRAN ${proc}
  ./tlau/lauX.sh -CPP ${proc}
done

