#!/bin/bash

if [ "$1" == "-short" ]; then
  ./tmad/teeMadX.sh -eemumu -ggtt -ggttg -ggttgg -makeclean
elif [ "$1" == "" ]; then
  ./tmad/teeMadX.sh -eemumu -ggtt -ggttg -ggttgg -ggttggg -makeclean
else
  echo "Usage: $0 [-short]"
  exit 1
fi
