#!/bin/bash

# wrapper script to determine the code version
# usage: getversion.sh <CPMD root directory>

# check if git is available
if [ $# -ne 1 ]
then 
  echo "usage: getversion.sh <CPMD root directory>"
  exit 1
fi
# default version from VERSION file
if [ -f $1/VERSION ]
then
  VERSION="$(cat $1/VERSION)"
else
  echo "VERSION file not found. Exiting" 1>&2
  exit 2
fi
# check for git command
if !command -v git &> /dev/null
then
   echo "git command not available. Using VERSION from file" 1>&2
else
# check if we have a git repository
  if git -C $1 rev-parse --is-inside-work-tree &> /dev/null
  then
     VERSION=$(git -C $1 describe --tags --abbrev=1 --always)
  else
     echo "CPMD-root not a git repository. Using VERSION from file" 1>&2
  fi 
fi

echo $VERSION
