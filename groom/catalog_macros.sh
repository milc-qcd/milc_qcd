#! /bin/sh

# This script scans the code for #ifdef's and #if defined's
# and lists the macros in use

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

LIBDIR=libraries
INCLUDEDIR=include
FILEUTILITYDIR=file_utilities

for dir in $GENDIRS $APPDIRS $LIBDIR $INCLUDEDIR $FILEUTILITYDIR
do
    echo "--------------------------------------"
    echo $dir
    echo "--------------------------------------"
    
    cd $MILCHOME/$dir
    scratch=/tmp/macro_catalog.$dir
    cat /dev/null > $scratch
    cfiles=`ls -1 *.c`
    hfiles=`ls -1 *.h`
    for file in $cfiles $hfiles
    do
	grep -E "^#if|^#els" $file |\
        awk -f $MILCHOME/groom/macros.awk >> $scratch
    done

    # Process raw scratch file
    
    sort $scratch | uniq
    /bin/rm $scratch
done
	