#!/bin/bash

# Run as:
#	bash compile_quda.sh


module reset

echo "module list:"
module list

echo "LD_LIBRARY_PATH:"
echo $LD_LIBRARY_PATH


if [ -d quda ]
then
  cd quda
  git pull
  git checkout develop
else
  git clone https://github.com/lattice/quda
  cd quda
  git checkout develop
fi
cd ..

if [ -d build ]
then
  cd build
else
  mkdir build
  cd build
fi

cmake ../quda -DCMAKE_BUILD_TYPE=RELEASE \
	-DCMAKE_INSTALL_PREFIX=`pwd`/usqcd \
	-DQUDA_GPU_ARCH=sm_80 -DQUDA_DIRAC_DEFAULT_OFF=ON -DQUDA_DIRAC_STAGGERED=ON \
	-DQUDA_QMP=ON -DQUDA_QIO=ON \
	-DQUDA_MULTIGRID=OFF \
	-DQUDA_SMEAR_GAUSS_TWOLINK=ON \
	-DQUDA_DOWNLOAD_USQCD=ON -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpiCC \

make -j 32 >& make_quda.log
make install
