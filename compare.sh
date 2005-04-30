#! /bin/sh

# For comparing test and sample output by hand
# Run this script from the application directory

PROJ=$1
PATTERNS=`grep PATTERNS Make_test | awk '{print $(NF-1),$NF}'`

../headtail.pl $PATTERNS < out.test.$PROJ > out.test.$PROJ.tmp
../headtail.pl $PATTERNS < out.sample.$PROJ > out.sample.$PROJ.tmp
../diffn3.pl out.test.$PROJ.tmp out.sample.$PROJ.tmp out.errtol.$PROJ

