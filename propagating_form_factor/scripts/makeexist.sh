#! /bin/sh

# Arguments
case $# in
[0-1]) echo "Usage: $0 <file> <cfg>"
   exit 1
esac

file=$1
cfg=$2

if [ ! -e $file ]
then
  make -f Make_select $file "CFG=${cfg}"
fi
