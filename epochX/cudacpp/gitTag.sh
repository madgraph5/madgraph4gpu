#!/bin/bash
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Sep 2024) for the MG5aMC CUDACPP plugin (based on earlier work for the HEP_OSlibs project).
# Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

set -e # fail on error

cd $(dirname $0)
topdir=$(cd ../..; pwd) # the top directory in this local repo
scr=$(basename $(cd ..; pwd))/$(basename $(pwd))/$(basename $0) # the name of this script in the git repo

skipFetch=0
###skipFetch=1 # FOR DEBUGGING!

# The tag prefix used for all tags handled by this script
PREFIX=TEST_cudacpp

# Usage
function usage()
{
  echo "Usage (1): $0 [-f] <tagsuffix>"
  echo "Creates a new version tag and pushes it to the remote repository"
  echo "The github CI will then create also a running tag with '_latest' suffix"
  echo "Valid formats for <tagsuffix> are 'n1.n2.n3' or 'n1.n2.n3_txt' where txt only contains letters or digits"
  echo "Version number 'n1.n2.n3' must match that in the CUDACPP_OUTPUT/__init__.py file"
  echo "Use the -f option to delete and recreate a version tag that already exists"
  echo ""
  echo "Usage (2): $0 -l"
  echo "Shows all available git tags (these are not necessarily github releases)"
  exit 1
}

# Command line arguments
force=""
if [ "$1" == "-f" ]; then
  force=$1
  shift
fi
if [ "$1" == "" ] || [ "$2" != "" ]; then
  usage
elif [ "$1" == "-l" ]; then
  tagsuffix=""
elif [ "${1/-}" != "${1}" ]; then # exclude immediately tagsuffix beginning with "-" (print usage on '-h')
  usage
else
  tagsuffix="$1"
  shift
fi

# Determine the local git branch
brn=`git branch | grep \* | cut -d ' ' -f2`
echo "INFO: git branch is  '${brn}'"
if [ "${brn}" != "install" ]; then # TEMPORARY! this should be branch 'master' eventually...
  echo "ERROR! Invalid local branch ${brn}"
  exit 1
fi

# Determine the remote git origin and the corresponding project name
# Stop if the remote git origin is not at github.com
url=`git config --get remote.origin.url`
echo "INFO: git remote is  '${url}'"
url=${url%.git}
prj=${url##*github.com:}
if [ "${url}" == "${prj}" ]; then
  echo "ERROR! git remote does not seem to be at github.com"
  exit 1
fi
echo "INFO: git project is '${prj}'"
if [ "${prj}" != "valassi/madgraph4gpu" ]; then # TEMPORARY! this will change eventually...
  echo "ERROR! Invalid project ${prj}"
  exit 1
fi

# Check that all changes are committed to git (except at most for this script only in the '-l' mode)
if [ "$(git status --porcelain | grep -v ^??)" != "" ]; then
  if [ "$tagsuffix" != "" ] || [ "$(git status --porcelain | grep -v ^??)" != " M ${scr}" ]; then
    echo "ERROR! There are local changes not committed to git yet"
    exit 1
  fi
fi

# Check that the local and remote branch are in sync
# See https://stackoverflow.com/a/3278427, https://stackoverflow.com/a/17938274
if [ "${skipFetch}" == "0" ]; then
  git fetch
fi
if [ "$(git rev-parse @{0})" != "$(git rev-parse @{u})" ]; then
  echo "ERROR! You need to git push or git pull"
  git status
  exit 1
fi

# Remove local git tags that are no longer on the remote repo
# See https://stackoverflow.com/a/27254532
if [ "${skipFetch}" == "0" ]; then
  git fetch --prune origin +refs/tags/*:refs/tags/*
fi
echo ""

# Retrieve the list of existing tags
existing_tags=$(git tag -l | grep ^${PREFIX} || true)

#
# OPTION 1 - LIST TAGS
#
if [ "$tagsuffix" == "" ]; then

  # See https://stackoverflow.com/questions/13208734
  echo "INFO: list existing tags starting with '${PREFIX}' (detailed list)"
  ###git --no-pager log --oneline --decorate --tags --no-walk
  if [ "${existing_tags}" == "" ]; then echo "[None]"; fi
  for tag in ${existing_tags}; do
    echo "--------------------------------------------------------------------------------"
    echo "$tag"
    echo "--------------------------------------------------------------------------------"
    ###git --no-pager for-each-ref --format="%(refname)" refs/tags/${tag}
    git --no-pager for-each-ref --format="TagDate:       %(creatordate)" refs/tags/${tag}
    git --no-pager log -n1 ${tag} --pretty=format:'commit         %h%nAuthor:        %an%nAuthorDate:    %ad%nMessage:       "%s"%n'
    echo ""
  done

#
# OPTION 2 - CREATE NEW TAGS
#
else

  # Determine cudacpp_version (as in archiver.sh)
  echo "INFO: determine cudacpp and mg5amc versions"
  cudacpp_version=$(cat ${topdir}/epochX/cudacpp/CODEGEN/PLUGIN/CUDACPP_SA_OUTPUT/__init__.py | awk '/__version__/{print $3}' | sed 's/(//' | sed 's/)//' | sed 's/,/./g')
  echo "> cudacpp_version = $cudacpp_version"

  # Determine mg5_version (as in HEPToolInstaller.py)
  mg5VERSION=${topdir}/MG5aMC/mg5amcnlo/VERSION 
  if [ ! -f ${mg5VERSION} ]; then echo "ERROR! File ${mg5VERSION} not found"; exit 1; fi 
  mg5_version=$(python3 -c "import re; print(re.findall(r'version\s*=\s*([\.\d]*)','$(cat ${mg5VERSION} | grep ^version)')[0])")
  echo "> mg5_version_current = $mg5_version"
  echo ""

  # Validate tagsuffix
  # See https://stackoverflow.com/questions/229551
  # See https://stackoverflow.com/questions/55486225
  echo "INFO: validate tag suffix '${tagsuffix}'"
  if [[ $tagsuffix == *" "* ]]; then
    echo "ERROR! Invalid tag suffix '${tagsuffix}' (no spaces allowed)"; exit 1
  fi
  if [ `python3 -c "import re; print(re.match('^[0-9]+\.[0-9]+\.[0-9]+(|_[0-9a-z]+)$','${tagsuffix}') is not None)"` == "False" ]; then
    echo "ERROR! Invalid tag suffix '${tagsuffix}' (valid formats are 'n1.n2.n3' or 'n1.n2.n3_txt' where txt only contains letters or digits)"; exit 1
  fi
  if [ "${tagsuffix/${cudacpp_version}}" == "${tagsuffix}" ]; then
    echo "ERROR! Invalid tag suffix '${tagsuffix}' (the version number does not match the cudacpp_version)"; exit 1
  fi
  echo ""

  # List the existing tags
  echo "INFO: list existing tags starting with '${PREFIX}' (short list)"
  if [ "${existing_tags}" == "" ]; then echo "> [None]"; fi
  for tag in ${existing_tags}; do
    echo "> $tag"
  done
  echo ""

  # Build the version tag and running tag names
  versTAG=${PREFIX}_for${mg5_version}_v${tagsuffix}
  runnTAG=${PREFIX}_for${mg5_version}_latest
  echo "INFO: will create version tag ${versTAG}"
  echo "(the CI will create running tag ${runnTAG})"

  # Check if the version tag already exists
  for tag in ${existing_tags}; do
    if [ "${versTAG}" == "${tag}" ] && [ "${force}" != "-f" ]; then
      echo "ERROR! Tag ${tag} already exists: use the -f option to recreate it"
      exit 1
    fi
  done
  echo ""

  # Create the version (~permanent) tag
  # (optionally use '-f' to replace any previously existing tag with the same name)
  echo "INFO: create version tag ${versTAG}"
  git tag ${force} ${versTAG} -m "Version tag ${versTAG}" -m "Tag created on $(date)"

  # Push the tag to the remote repository
  # (use '-f' to replace any previously existing tag with the same name)
  echo "INFO: push tag to the remote repository"
  git push -f --tags

fi
