#! /bin/sh

# Makes 3pt file lists based on supplied getlist and hitlist
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

HL3_p0name=HL3_p0_0${kappa_spect}
HL3_p1name=HL3_p1_0${kappa_spect}
HL3_p2name=HL3_p2_0${kappa_spect}

HH3_p0name=HH3_p0_0${kappa_spect}
HH3_p1name=HH3_p1_0${kappa_spect}
HH3_p2name=HH3_p2_0${kappa_spect}

# There is a separate subdirectory for each configuration
#datadir=/scratch/scphdeta/prop_form/b560m01/data
datadir=/work/u2219/prop_form/b560m01

cat /dev/null > HL3_p0list
cat /dev/null > HL3_p1list
cat /dev/null > HL3_p2list

cat /dev/null > HH3_p0list
cat /dev/null > HH3_p1list
cat /dev/null > HH3_p2list

for cfg in `cat ${getlist}`
do
  if ! grep ${cfg} ${hitlist} > /dev/null
  then

#   Heavy-light 3 pt function

    for name3 in $HL3_p0name $HL3_p1name $HL3_p2name
    do
      file3="${datadir}/${cfg}/${name3}.${cfg}"
      makeexist.sh $file3 $cfg
 
      # Name of list file is derived from the 1st 6 letters of the file name
      list3=`echo $name3 | awk '{print substr($1,1,6)}'`
      list3=${list3}list

      echo ${file3} >> ${list3}
      echo "   ${file2}" >> ${list3}
      echo "   ${file1}" >> ${list3}

    done

#   Heavy-heavy 3 pt function

    for name3 in $HH3_p0name $HH3_p1name $HH3_p2name
    do
      file3="${datadir}/${cfg}/${name3}.${cfg}"
#      makeexist.sh $file3 $cfg
 
      # Name of list file is derived from the 1st 6 letters of the file name
      list3=`echo $name3 | awk '{print substr($1,1,6)}'`
      list3=${list3}list

      echo ${file3} >> ${list3}
      echo "   ${file2}" >> ${list3}
      echo "   ${file1}" >> ${list3}

    done
  fi
done



