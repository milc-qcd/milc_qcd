#! /bin/sh

# Makes 3pt and 2pt file lists based on getlist and hitlist
# Also calls "make" to fetch files from Unitree

# Usage

#    setup_filelist.sh <getlist> <hitlist> <kappa_spect> <spect_index>

# kappa_spect should be given without the leading 0, as in ".1245"

# Arguments
case $# in
[0-3]) echo "Usage: $0 <getlist> <hitlist> <kappa_spect (.nnnn)> <spect_index>"
       exit 1
esac

getlist=$1
hitlist=$2
kappa_spect=$3
sp=$4

# Parameters for this script
# Fragments of file names

HH2_GLname=HH2_GL
HL2_GEname=HL2_GE_0${kappa_spect}
HL2_GGname=HL2_GG_0${kappa_spect}
HL2_GLname=HL2_GL_0${kappa_spect}
LL2_GGname=LL2_GG_0${kappa_spect}

HL3_p0name=HL3_p0_0${kappa_spect}
HL3_p1name=HL3_p1_0${kappa_spect}
HL3_p2name=HL3_p2_0${kappa_spect}

HH3_p0name=HH3_p0_0${kappa_spect}
HH3_p1name=HH3_p1_0${kappa_spect}
HH3_p2name=HH3_p2_0${kappa_spect}

# There is a separate subdirectory for each configuration
#datadir=/scratch/scphdeta/prop_form/b560m01/data
datadir=/work/u2219/prop_form/b560m01

cat /dev/null > HH2_GL_sp${sp}list
cat /dev/null > HL2_GE_sp${sp}list
cat /dev/null > HL2_GG_sp${sp}list
cat /dev/null > HL2_GL_sp${sp}list
cat /dev/null > LL2_GG_sp${sp}list

cat /dev/null > HL3_p0_sp${sp}list
cat /dev/null > HL3_p1_sp${sp}list
cat /dev/null > HL3_p2_sp${sp}list

cat /dev/null > HH3_p0_sp${sp}list
cat /dev/null > HH3_p1_sp${sp}list
cat /dev/null > HH3_p2_sp${sp}list

for cfg in `cat ${getlist}`
do
  if ! grep ${cfg} ${hitlist} > /dev/null
  then

#   Heavy-light 3 pt function

    file1="${datadir}/${cfg}/${HL2_GEname}.${cfg}"
    file2="${datadir}/${cfg}/${LL2_GGname}.${cfg}"

    makeexist.sh $file1 $cfg
    makeexist.sh $file2 $cfg

    for name3 in $HL3_p0name $HL3_p1name $HL3_p2name
    do
      file3="${datadir}/${cfg}/${name3}.${cfg}"
      makeexist.sh $file3 $cfg
 
      # Name of list file is derived from the 1st 6 letters of the file name
      list3=`echo $name3 | awk '{print substr($1,1,6)}'`
      list3=${list3}_sp${sp}list

      echo ${file3} >> ${list3}
      echo "   ${file2}" >> ${list3}
      echo "   ${file1}" >> ${list3}

    done

#   Heavy-heavy 3 pt function

    file1="${datadir}/${cfg}/${HL2_GEname}.${cfg}"
    file2="${datadir}/${cfg}/${HL2_GGname}.${cfg}"

    makeexist.sh $file1 $cfg
    makeexist.sh $file2 $cfg

    for name3 in $HH3_p0name $HH3_p1name $HH3_p2name
    do
      file3="${datadir}/${cfg}/${name3}.${cfg}"
      makeexist.sh $file3 $cfg
 
      # Name of list file is derived from the 1st 6 letters of the file name
      list3=`echo $name3 | awk '{print substr($1,1,6)}'`
      list3=${list3}_sp${sp}list

      echo ${file3} >> ${list3}
      echo "   ${file2}" >> ${list3}
      echo "   ${file1}" >> ${list3}

    done

# Two pt functions alone

    name2list="$HH2_GLname $HL2_GEname $HL2_GGname $HL2_GLname $LL2_GGname"
    for name2 in $name2list
    do
      file2="${datadir}/${cfg}/${name2}.${cfg}"

      # Name of list file is derived from the 1st 6 letters of the file name
      list2=`echo $name2 | awk '{print substr($1,1,6)}'`
      list2=${list2}_sp${sp}list

      makeexist.sh $file2 $cfg
      echo ${file2} >> ${list2}
    done

  fi
done



