#! /bin/sh

# Makes 2pt file lists based on specified getlist and hitlist
# Also calls "make" to fetch files from Unitree

# Usage

#    setup_filelist.sh <getlist> <hitlist> <kappa_spect>

# kappa_spect should be given without the leading 0, as in ".1245"

# Arguments
case $# in
[0-2]) echo "Usage: $0 <getlist> <hitlist> <kappa_spect>"
       exit 1
esac

getlist=$1
hitlist=$2
kappa_spect=$3

# Parameters for this script
# Fragments of file names

HH2_GLname=HH2_GL
HL2_GEname=HL2_GE_0${kappa_spect}
HL2_GGname=HL2_GG_0${kappa_spect}
HL2_GLname=HL2_GL_0${kappa_spect}
LL2_GGname=LL2_GG_0${kappa_spect}

# There is a separate subdirectory for each configuration
#datadir=/scratch/scphdeta/prop_form/b560m01/data
datadir=/work/u2219/prop_form/b560m01

cat /dev/null > HH2_GLlist
cat /dev/null > HL2_GElist
cat /dev/null > HL2_GGlist
cat /dev/null > HL2_GLlist
cat /dev/null > LL2_GGlist

for cfg in `cat ${getlist}`
do
  if ! grep ${cfg} ${hitlist} > /dev/null
  then

# Two pt functions alone

    name2list="$HH2_GLname $HL2_GEname $HL2_GGname $HL2_GLname $LL2_GGname"
    for name2 in $name2list
    do
      file2="${datadir}/${cfg}/${name2}.${cfg}"

      # Name of list file is derived from the 1st 6 letters of the file name
      list2=`echo $name2 | awk '{print substr($1,1,6)}'`
      list2=${list2}list

      makeexist.sh $file2 $cfg
      echo ${file2} >> ${list2}
    done

  fi
done



