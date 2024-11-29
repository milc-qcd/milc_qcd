#!/bin/bash

# Run as:
#	bash compile_qioqmp.sh

echo "module list:"
module list

echo "LD_LIBRARY_PATH:"
echo $LD_LIBRARY_PATH

# Build location
install=`pwd`

# Download and prepare directory structure
cd ${install}
if [ -d qmp ]
then
  cd qmp
  git pull
else
  git clone https://github.com/usqcd-software/qmp.git
fi

cd ${install}
if [ -d qio ]
then
  cd qio
  git pull
else
  git clone --recursive https://github.com/usqcd-software/qio.git
  cd qio
  bash autogen.sh
fi

cd ${install}
mkdir -p build/qmp build/qio

# Build QMP
echo -e "\nConfiguring and building QMP..."
cd ${install}/build/qmp
../../qmp/configure \
     --prefix=${install}/build/qmp/qmp \
     --with-qmp-comms-type=MPI \
     CFLAGS="-O3"

make
make install

pwd
echo $QMP_LIBS

# Build QIO
echo -e "\nConfiguring and building QIO..."
cd ${install}/build/qio
../../qio/configure \
     --prefix=${install}/build/qio/qio \
     CFLAGS="-O3" \
     --with-qmp=${install}/build/qmp/qmp \
     --enable-qmp-route

make
make install

