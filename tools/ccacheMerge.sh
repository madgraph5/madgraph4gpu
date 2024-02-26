#!/bin/bash

# Merge the contents of srcCcacheDir into dstCcacheDir
# See https://lists.samba.org/archive/ccache/2018q1/001531.html
if [ "$2" == "" ] || [ "$3" != "" ]; then
  echo "Usage: $0 srcCcacheDir dstCcacheDir"
  exit 1
fi
src=$1
dst=$2

if [ ! -d $src ]; then echo "ERROR! $src directory not found"; exit 1; fi
if [ ! -d $dst ]; then echo "ERROR! $dst directory not found"; exit 1; fi

if ! ccache --version | grep 'ccache version'; then echo "ERROR! ccache executable not found"; exit 1; fi

echo
date

echo
echo "=== Before the merge:"
echo "Source ccache directory:"
echo $src
echo "  Disk usage: $(du -sb $src | awk '{print $1}') bytes"
echo "  Files on disk: $(find $src -type f | wc -l)"
ccache -d $src -s | grep 'Cache size'
echo "Destination ccache directory:"
echo $dst
echo "  Disk usage: $(du -sb $dst | awk '{print $1}') bytes"
echo "  Files on disk: $(find $dst -type f | wc -l)"
ccache -d $dst -s | grep 'Cache size'

echo
echo "=== Merging..."
# Note1: ? is a glob for single-character directories 0 1 2 3 4 5 6 7 8 9 a b c d e f
# Note2: cp uses the default settings and overwrites by default
# Note3: any files already present in the dst directory but not in the source directory are kept
###cd $src; cp -a --parents ? $dst # original implementation do-nothing 80s for 9000 files
# Note1: -a == -rlptgoD (recursive; keep links/perms/times/group/owner/special), use this
# Note2: -P == --partial --progress (verbose progress display), use this
# Note3: -u == --update (skip files that are newer on the receiver), use this
# Note4: --delete (delete files that only exist on the receiver), do NOT use this
# Note5: trailing / in src, no trailing slash in dst!
rsync -aPu $src/ $dst # new faster implementation do-nothing 1s for 9000 files
ccache -d $dst -c

echo
echo "=== After the merge:"
echo "Source ccache directory:"
echo $src
echo "  Disk usage: $(du -sb $src | awk '{print $1}') bytes"
echo "  Files on disk: $(find $src -type f | wc -l)"
ccache -d $src -s | grep 'Cache size'
echo "Destination ccache directory:"
echo $dst
echo "  Disk usage: $(du -sb $dst | awk '{print $1}') bytes"
echo "  Files on disk: $(find $dst -type f | wc -l)"
ccache -d $dst -s | grep 'Cache size'

echo
date
