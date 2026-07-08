#!/bin/bash

# Run as:
#	bash compile_quda.sh

module reset
module load cmake
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
  git checkout feature/sycl
fi
cd ..

if [ -d build ]
then
  cd build
else
  mkdir build
  cd build
fi

export QUDA_SYCL_TARGETS="spir64_gen"
export SYCL_LINK_FLAGS="$SYCL_LINK_FLAGS -Xs \"-device pvc\""
export SYCL_LINK_FLAGS="$SYCL_LINK_FLAGS -fsycl-device-code-split=per_kernel"
export SYCL_LINK_FLAGS="$SYCL_LINK_FLAGS -fsycl-max-parallel-link-jobs=32"
export SYCL_LINK_FLAGS="$SYCL_LINK_FLAGS -flink-huge-device-code"

export QUDA_WARP_SIZE=16
#export QUDA_WARP_SIZE=32
export QUDA_MAX_BLOCK_SIZE=1024
export QUDA_MAX_ARGUMENT_SIZE=2048
export QUDA_TEST_NUMPROCS=1

LARGE_REGS=1    #turn off for now
if [ "X$LARGE_REGS" = "X1" ]; then
  echo "Using large register file"
  #export CXXFLAGS="$CXXFLAGS -Xs \"-options -ze-opt-large-register-file\""
  #export LDFLAGS="$LDFLAGS -Xs \"-options -ze-opt-large-register-file\""
  #export SYCL_FLAGS="$SYCL_FLAGS -Xs \"-options -ze-opt-large-register-file\""
  export SYCL_LINK_FLAGS="$SYCL_LINK_FLAGS -Xs \"-options -ze-opt-large-register-file\""
  export QUDA_MAX_BLOCK_SIZE=512
fi

echo "QUDA_TARGET:"
echo $QUDA_TARGET
echo "QUDA_SYCL_TARGETS:"
echo $QUDA_SYCL_TARGETS
echo "SYCL_FLAGS:"
echo $SYCL_FLAGS
echo "SYCL_LINK_FLAGS:"
echo $SYCL_LINK_FLAGS
echo "CXXFLAGS:"
echo $CXXFLAGS
echo "LDFLAGS:"
echo $LDFLAGS

cmake ../quda -DCMAKE_BUILD_TYPE=RELEASE \
        -DCMAKE_INSTALL_PREFIX=`pwd`/usqcd \
        -DQUDA_BUILD_SHAREDLIB=ON \
        -DQUDA_DIRAC_DEFAULT_OFF=ON \
        -DQUDA_DIRAC_STAGGERED=ON \
        -DQUDA_DOWNLOAD_USQCD=ON -DQUDA_QMP=ON -DQUDA_QIO=ON \
        -DQUDA_MULTIGRID=OFF \
        -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
        -DQUDA_TARGET_TYPE=SYCL
        
cmake --build . -j 16 -v
cmake --install . 
