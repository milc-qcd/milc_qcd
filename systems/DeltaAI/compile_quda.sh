#!/bin/bash

# Run as:
#	bash compile_quda.sh


module load gcc-native/12.3

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

## now explictly adding this variable
export CRAY_ACCEL_TARGET=nvidia90

CUBLAS=/opt/nvidia/hpc_sdk/Linux_aarch64/24.3/math_libs/12.3/lib64/libcublas.so
CUFFT=/opt/nvidia/hpc_sdk/Linux_aarch64/24.3/math_libs/12.3/lib64/libcufft.so

cmake ../quda -DCMAKE_BUILD_TYPE=RELEASE \
	-DCMAKE_INSTALL_PREFIX=`pwd`/usqcd \
	-DQUDA_GPU_ARCH=sm_90 -DQUDA_DIRAC_DEFAULT_OFF=ON -DQUDA_DIRAC_STAGGERED=ON \
	-DQUDA_QMP=ON -DQUDA_QIO=ON \
	-DQUDA_MULTIGRID=OFF \
	-DQUDA_SMEAR_GAUSS_TWOLINK=ON \
	-DQUDA_DOWNLOAD_USQCD=ON -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC \
	-DCUDA_cublas_LIBRARY=$CUBLAS \
	-DCUDA_cufft_LIBRARY=$CUFFT \

make -j 32 >& make_quda.log
make install
