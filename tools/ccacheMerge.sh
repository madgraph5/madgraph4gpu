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
cd $src
# Note1: ? is a glob for single-character directories 0 1 2 3 4 5 6 7 8 9 a b c d e f
# Note2: cp uses the default settings and overwrites by default
# Note3: any files already present in the dst directory but not in the source directory are kept
cp -a --parents ? $dst
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
