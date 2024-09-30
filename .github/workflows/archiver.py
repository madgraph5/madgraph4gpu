#!/bin/bash
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: O. Mattelaer (Sep 2024) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

import subprocess
import sys

def get_all_tags():
    out = subprocess.check_output(['git', 'tag']).decode()
    return out.split('\n')

def get_supported_versions(tags):
    PREFIX = 'cudacpp_for'
    SUFFIX = '_latest'
    versions = [ t[len(PREFIX):-len(SUFFIX)] for t in tags if t.startswith(PREFIX) and t.endswith(SUFFIX)]
    versions = set(versions)
    return versions

def create_info_file(path, versions):
    line = "%(version)s https://github.com/valassi/madgraph4gpu/releases/download/TEST_cudacpp_for%(version)s_latest/cudacpp.tar.gz\n"
    fsock = open(path, 'w')
    for v in versions:
        fsock.write(line%{'version':v})

if "__main__" == __name__:
    if len(sys.argv) != 2:
        print('Usage: python3 %s <filename>'%sys.argv[0])
        sys.exit(1)
    tags = get_all_tags()
    print('Tags:', tags)
    versions = get_supported_versions(tags)
    print('Versions:', versions)
    create_info_file(sys.argv[1], versions)
