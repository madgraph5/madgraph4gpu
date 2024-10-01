#!/bin/bash
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: O. Mattelaer (Sep 2024) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

import subprocess
import sys

def get_repo():
    out = subprocess.check_output(['git', 'config', '--get', 'remote.origin.url']).decode()
    return out.split(':')[1].split('.')[0].split('/') # returns (repo_owner, repo_name)

def get_all_tags():
    out = subprocess.check_output(['git', 'tag']).decode()
    return out.split('\n')

def get_supported_versions(tags):
    versions = [ t[len(PREFIX):-len(SUFFIX)] for t in tags if t.startswith(PREFIX) and t.endswith(SUFFIX)]
    versions = set(versions)
    return versions

def create_infodat_file(path, versions):
    line = "%(version)s https://github.com/%(repo)s/releases/download/%(prefix)s%(version)s%(suffix)s/cudacpp.tar.gz\n"
    with open(path, 'w') as fsock:
        for v in versions:
            fsock.write(line%{'repo':GITHUB_REPO, 'prefix':PREFIX, 'version':v, 'suffix':SUFFIX})

if "__main__" == __name__:
    if len(sys.argv) != 2:
        print('Usage: python3 %s <infodat>'%sys.argv[0])
        sys.exit(1)
    repo_owner, repo_name = get_repo()
    ###print('Repo owner:', repo_owner)
    ###print('Repo name:', repo_name)
    GITHUB_REPO = '%s/%s/'%(repo_owner, repo_name)
    PREFIX = 'cudacpp_for'
    if repo_owner != 'madgraph5' : PREFIX = repo_owner + "_" + PREFIX # TEMPORARY! this will change eventually...
    if repo_name != 'madgraph4gpu' : raise Exception('Invalid repo_name (expect "madgraph4gpu")') # TEMPORARY! this will change eventually...
    SUFFIX = '_latest'
    tags = get_all_tags()
    ###print('Tags:', tags)
    versions = get_supported_versions(tags)
    ###print('Supported versions:', versions)
    create_infodat_file(sys.argv[1], versions)
