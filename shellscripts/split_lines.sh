#!/bin/bash

if ! [ -d $1 ]; then
  echo "Must pass a directory name to this script"
  exit
fi

function linenum {
  IFS=.
  read -a fnarr <<< "$1"
  echo ${fnarr[2]}
}

for f in $1/*.nc; do
  base=$(basename $f)
  dirname="lines/$(linenum $base)"
  echo ${dirname}
  if ! [ -d $dirname ]; then
    mkdir -p $dirname
  fi
  cp $f $dirname
done 
