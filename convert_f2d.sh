#! /bin/sh

# For use with the MILC code v6.21aug02

# ancient D. Toussaint - Created
# 8/15/02 C. DeTar - Integrated with v6.31aug02

# Converts the entire code from single precision to double precision

# Must be run from the top level directory, i.e. the directory
# containing this script.

# Converts "float" to "double" for all *.c and *.h files in all
# subdirectories with some exceptions
# Makes miscellaneous other changes

# After conversion, be sure to build libraries with Make_vanilla.
# until we have double precision assembly language code.

# Usage...

#  convert_f2d.sh

echo "WARNING: This script changes ALL code in the current directory"
echo "and ALL of its subdirectories."
echo "ARE YOU SURE YOU WANT TO DO THIS (yes/no)?"
read reply
if [ "$reply" != "yes" ]
then
  exit 1
fi

echo "List of edited files"

# Edit all *.c and *.h files
for f in `find . -name "*.[ch]" -print`
do
# Files to leave intact:
if [ $f != "./generic/io_lat4_double.c" -a $f != "./generic/io_wb3_double.c" -a $f != "./generic_form/load_smearing_double.c" -a $f != "./file_utilities/check_gauge.c" -a $f != "./file_utilities/check_prop.c" ]
then
  echo $f

  ex - $f << EOF
g/floatsum/s//DUMMY1/g
g/floatmax/s//DUMMY3/g
g/broadcast_float/s//DUMMY2/g
g/fcomplex/s/float/FLOAT/g
g/float/s//Real/g
g/MPI_FLOAT/s//MILC_MPI_REAL/g
g/FLOAT/s//MILC_REAL/g
g/DUMMY1/s//floatsum/g
g/DUMMY3/s//floatmax/g
g/DUMMY2/s//broadcast_float/g
g/scanf/s/%e/%leHELP/g
g/scanf/s/%f/%lfHELP/g
g/rsum/s//sumReal/g
g/gssum/s//gsumHELP/g
g/gshigh/s//gshighHELP/g
g/pkcplx/s//pkcplxHELP/g
g/
x
EOF

fi
done

# Edit all Make_template files to plug in the double versions of file I/O 
for f in `find . -name Make_template -print`
do
if [ $f != "./generic/Make_template" -a $f != "./generic_form/Make_template" ]
then
  echo $f

  ex - $f << EOF
g/io_lat4\.o/s//${IO_LAT_REAL}/g
g/io_wb3\.o/s//${IO_PROP_REAL}/g
g/load_smearing\.o/s//${LOAD_SMEARING_REAL}/g
g/
x
EOF
fi

done

