#!/bin/bash

# Run as:
#	bash compile_quda.sh

################### Setup environment for Frontier

module reset
module load PrgEnv-amd amd/7.1.1 rocm/7.1.1
module load craype-accel-amd-gfx90a
module load cmake
module load ninja
module list

# Per OLCF docs, set HIPFLAGS when compiling without Cray compiler wrappers
HIPFLAGS="--offload-arch=gfx90a "

MY_CFLAGS="-I${MPICH_DIR}/include --offload-arch=gfx90a -g -pg"
MY_LDFLAGS="-Wl,-rpath=${MPICH_DIR}/lib -L${MPICH_DIR}/lib -lmpi ${CRAY_XPMEM_POST_LINK_OPTS} -lxpmem ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a} --offload-arch=gfx90a -g -pg"

################### Get QUDA and prepare build/ directory

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

################### Compile

cmake ../quda \
    -G "Ninja" \
    -DQUDA_TARGET_TYPE=HIP \
    -DQUDA_GPU_ARCH=gfx90a \
    -DROCM_PATH=${ROCM_PATH} \
    -DCMAKE_INSTALL_PREFIX=`pwd`/usqcd \
    -DCMAKE_BUILD_TYPE=RELEASE \
    -DQUDA_DIRAC_DEFAULT_OFF=ON \
    -DQUDA_DIRAC_STAGGERED=ON \
    -DQUDA_BUILD_SHAREDLIB=ON \
    -DQUDA_QMP=ON \
    -DQUDA_QIO=ON \
    -DQUDA_DOWNLOAD_USQCD=ON \
    -DQUDA_MULTIGRID=OFF \
    -DCMAKE_BUILD_TYPE="DEVEL" \
    -DCMAKE_CXX_COMPILER="hipcc" \
    -DCMAKE_C_COMPILER="hipcc" \
    -DBUILD_SHARED_LIBS=ON \
    -DQUDA_BUILD_SHAREDLIB=ON \
    -DQUDA_BUILD_ALL_TESTS=ON \
    -DQUDA_CTEST_DISABLE_BENCHMARKS=ON \
    -DCMAKE_C_STANDARD=99 \
    -DCMAKE_CXX_FLAGS="${MY_CFLAGS}" \
    -DCMAKE_C_FLAGS="${MY_CFLAGS}" \
          -DCMAKE_HIP_FLAGS="${MY_CFLAGS}" \
    -DCMAKE_SHARED_LINKER_FLAGS="${MY_LDFLAGS}" \
    -DCMAKE_EXE_LINKER_FLAGS="${MY_LDFLAGS}" \

cmake --build . -j 16 -v 2>&1 | tee -a cmake_quda.log
echo "" | tee -a cmake_quda.out
cmake --install . 2>&1 | tee -a cmake_quda.log

