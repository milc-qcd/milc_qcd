#! /bin/sh

# It is assumed that the full test suite has been run with
# make -f Make_test_work all
# and the results are in */maketest.log

# This script searches the application directories for any C code file
# that was not mentioned in maketest.log

# (It can't tell whether the code was actually called from another one)
# We do not test use of the library routines

MILCHOME=$HOME/milc/v6   # Home directory for MILC code

# List of application directories

APPDIRS=" \
  arb_dirac_eigen \
  arb_dirac_invert \
  clover_dynamical \
  clover_hybrids \
  clover_invert \
  h_dibaryon \
  heavy \
  hqet_heavy_to_light \
  hvy_qpot \
  ks_dynamical \
  ks_eigen \
  ks_hybrids2 \
  ks_imp_dyn1 \
  ks_imp_dyn2 \
  propagating_form_factor \
  pure_gauge \
  schroed_cl_inv \
  schroed_ks_dyn \
  schroed_pg \
  smooth_inst \
  string_break \
  symanzik_sl32 \
  wilson_dynamical \
  wilson_hybrids \
  wilson_invert \
  wilson_static"

# List of supporting directories to be checked

GENDIRS="
  generic \
  generic_clover \
  generic_form \
  generic_ks \
  generic_pg \
  generic_schroed \
  generic_wilson"

echo "--------------------------------------"
echo "supporting directories"
echo "--------------------------------------"
for dir in $GENDIRS
do
    cd $MILCHOME/$dir
    cfiles=`ls -1 *.c`
    for cfile in $cfiles
    do
	grep -qw $cfile $MILCHOME/*/maketest.log
	if [ $? -eq 1 ]
	then
	    # Different message if file already marked unmaintained
	    grep -q "NOT MAINTAINED" $cfile
	    if [ $? -eq 1 ]
            then
		echo $dir/$cfile NOT COMPILED
 	    else
		echo $dir/$cfile unmaintained
	    fi
	fi
    done
done

for dir in $APPDIRS
do
    echo "--------------------------------------"
    echo $dir
    echo "--------------------------------------"

    cd $MILCHOME/$dir
    cfiles=`ls -1 *.c`
    for cfile in $cfiles
    do
	grep -qw $cfile $MILCHOME/*/maketest.log
	if [ $? -eq 1 ]
	then
	    # Different message if file already marked unmaintained
	    grep -q "NOT MAINTAINED" $cfile
	    if [ $? -eq 1 ]
            then
		echo $dir/$cfile NOT COMPILED
 	    else
		echo $dir/$cfile unmaintained
	    fi
	fi
    done
done


