#! /bin/sh

# For comparing test and sample output by hand
# Run this script from the application directory

# Example of usage

#   cd ks_spectrum
#   ../compare.sh ks_spectrum_asqtad.fpi.1

PROJ=$1
PATTERNS=`grep ${PROJ} test/checklist | awk '{print $(NF-1),$NF}'`

perl ../headtail.pl $PATTERNS < test/$PROJ.test-out > test/$PROJ.test-tmp
perl ../headtail.pl $PATTERNS < test/$PROJ.sample-out > test/$PROJ.sample-tmp
perl ../diffn3.pl test/$PROJ.test-tmp test/$PROJ.sample-tmp test/$PROJ.errtol

