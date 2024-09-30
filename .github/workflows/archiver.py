#!/bin/bash
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: O. Mattelaer (Sep 2024) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

import subprocess
import sys

GITHUB_REPO = 'valassi/madgraph4gpu/'

PREFIX = 'cudacpp_for'
SUFFIX = '_latest'

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

def create_infohtm_file(path, versions):
    from datetime import datetime
    line = "<li><a href=https://github.com/%(repo)s/releases/tag/%(prefix)s%(version)s%(suffix)s>%(version)s</a>\n"
    with open(path, 'w') as fsock:
        fsock.write("""<!DOCTYPE html>
<html>
<body>
<h3>List of cudacpp releases</h3>
<ul>
""")
        for v in versions:
            fsock.write(line%{'repo':GITHUB_REPO, 'prefix':PREFIX, 'version':v, 'suffix':SUFFIX})
        fsock.write("""</ul>
<hr>
<i>Last modified on %s</i>.
</body>
</html>"""%datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"))

if "__main__" == __name__:
    if len(sys.argv) != 3:
        print('Usage: python3 %s <infodat> <infohtml>'%sys.argv[0])
        sys.exit(1)
    tags = get_all_tags()
    ###print('Tags:', tags)
    versions = get_supported_versions(tags)
    print('Supported versions:', versions)
    create_infodat_file(sys.argv[1], versions)
    create_infohtm_file(sys.argv[2], versions)
