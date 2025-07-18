#!/bin/bash
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: O. Mattelaer (Sep 2024) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

import subprocess
import sys
import re

def get_all_tags():
    out = subprocess.check_output(['git', 'tag']).decode()
    return out.split('\n')

def get_supported_versions(tags):
    # versions = [ t[len(PREFIX):-len(SUFFIX)] for t in tags if t.startswith(PREFIX) and t.endswith(SUFFIX)]
    versions = re.findall(rf"{PREFIX}(\d+.\d+.\d+)_v(\d+.\d+.\d+)$","\n".join(tags), re.M)
    versions = set(versions)
    return versions

def create_infodat_file(path, mg_cudacpp_versions):
    line = "%(version)s https://github.com/%(repo)s/releases/download/%(prefix)s%(version)s_v%(suffix)s/cudacpp.tar.gz\n"
    with open(path, 'w+') as fsock:
        for v, c in mg_cudacpp_versions:
            fsock.write(line%{'repo':GITHUB_REPO, 'prefix':PREFIX, 'version':v, 'suffix':c})

if "__main__" == __name__:
    if len(sys.argv) != 3:
        print('Usage: python3 %s <repoowner/reponame> <infodat>'%sys.argv[0])
        sys.exit(1)
    print('Executing: python3 %s "%s" "%s"'%( sys.argv[0],sys.argv[1],sys.argv[2]))
    repo = sys.argv[1]
    infodat = sys.argv[2]
    repo_owner, repo_name = repo.split('/')
    ###print('Repo owner:', repo_owner)
    ###print('Repo name:', repo_name)
    GITHUB_REPO = '%s/%s/'%(repo_owner, repo_name)
    PREFIX = 'cudacpp_for'
    if repo_owner != 'madgraph5' : PREFIX = repo_owner + "_" + PREFIX # TEMPORARY! this will change eventually...
    if repo_name != 'madgraph4gpu' : raise Exception('Invalid repo_name "%s" (expect "madgraph4gpu")'%repo_name) # TEMPORARY! this will change eventually...
    tags = get_all_tags()
    # SUFFIX = '_latest'
    ###print('Tags:', tags)
    versions = get_supported_versions(tags)
    ###print('Supported versions:', versions)
    create_infodat_file(infodat, versions)
