#! /bin/sh

# For loosening error tolerances.  Inspect the differences first!
# Assumes that the file out.test.* has been generated already
# Run this script from the application directory

target=$1
if [ -z "$target" ]
then
  echo "Usage $0 <target>"
  exit 1
fi

perl ../../trainerrfile.pl $target.test-tmp $target.sample-tmp $target.errtol
perl ../../diffn3.pl $target.test-tmp $target.sample-tmp $target.errtol
