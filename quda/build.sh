#! /bin/bash

# Build quda
BRANCH=develop
ARCH=sm_80 # V100 sm_70; A100 sm_80
QUDA_HOME=`pwd`
QUDA_INSTALL=${QUDA_HOME}/install

if [ -d quda ]
then
  pushd quda
  git pull
  git checkout ${BRANCH}
  popd
else
  git clone -b ${BRANCH} https://github.com/lattice/quda
fi


mkdir -p build
pushd build

QUDA_BUILD=`pwd`
#MULTIGRID_FLAG="-DQUDA_MULTIGRID=ON -DQUDA_MULTIGRID_NVEC_LIST=24,64,96"
MULTIGRID_FLAG=""

cmake \
    -DCMAKE_BUILD_TYPE=RELEASE \
    -DCMAKE_INSTALL_PREFIX=${QUDA_INSTALL} \
    -DQUDA_BUILD_SHAREDLIB=ON \
    -DQUDA_COVDEV=ON \
    -DQUDA_DIRAC_DEFAULT_OFF=ON \
    -DQUDA_DIRAC_STAGGERED=ON \
    -DQUDA_DOWNLOAD_USQCD=ON -DQUDA_QIO=ON -DQUDA_QMP=ON \
    -DQUDA_GPU_ARCH=${ARCH} \
    ${MULTIGRID_FLAG} \
    ${QUDA_HOME}/quda

make -k -j16 install

popd


