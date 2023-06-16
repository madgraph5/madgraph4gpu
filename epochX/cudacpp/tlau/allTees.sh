#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (June 2023) for the MG5aMC CUDACPP plugin.

for dir in gg_tt.mad gg_ttg.mad gg_ttgg.mad gg_ttggg.mad; do
  ./tlau/lauX.sh -CUDA ${dir}
  ./tlau/lauX.sh -FORTRAN ${dir}
  ./tlau/lauX.sh -CPP ${dir}
done

