#!/bin/bash

# Run as:
#	bash compile_quda.sh

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

export CRAY_ACCEL_TARGET=nvidia80

CUBLAS_LIB=/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/math_libs/lib64/libcublas.so
CUFFT_LIB=/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/math_libs/lib64/libcufft.so

cmake ../quda -DCMAKE_BUILD_TYPE=RELEASE \
	-DCMAKE_INSTALL_PREFIX=`pwd`/usqcd \
	-DQUDA_GPU_ARCH=sm_80 -DQUDA_DIRAC_DEFAULT_OFF=ON -DQUDA_DIRAC_STAGGERED=ON \
	-DQUDA_QMP=ON -DQUDA_QIO=ON \
	-DQUDA_MULTIGRID=OFF \
	-DQUDA_SMEAR_GAUSS_TWOLINK=ON \
	-DQUDA_DOWNLOAD_USQCD=ON -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC \
	-DCUDA_cublas_LIBRARY=$CUBLAS_LIB \
	-DCUDA_cufft_LIBRARY=$CUFFT_LIB \

make -j 32 >& make_quda.log
make install
