#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.

diff --no-dereference -x '*log.txt' -x '*.o' -x '*.o.*' -x '*.a' -x '*.exe' -x 'lib' -x 'build.*' -x '.build.*' -x '*~' $* 
