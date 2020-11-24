#! /bin/bash

ARCH=$1     # Choices: scalar, avx2, avx512
PK_CC=$2
PK_CXX=$3
GIT_BRANCH=devel

if [ -z "${ARCH}" ]
then
   echo "Usage $0 <scalar|avx2|avx512> <CC> <CXX>"
   exit 1
fi

MAKE="make -j4"

TOPDIR=`pwd`
SRCDIR=${TOPDIR}/qphix
INSTALL_PREFIX=${TOPDIR}/install

QMP_DIR=${HOME}/scidac/install/qmp

if [ ! -d ${SRCDIR} ]
then
  echo "Fetching ${GIT_BRANCH} branch of package from github"
  git clone https://github.com/JeffersonLab/qphix/ -b ${GIT_BRANCH}
fi

pushd build

# Invoke cmake, but preset some of the cmake parameters
cmake \
-D isa:STRING=${ARCH} \
-D CMAKE_INSTALL_PREFIX:PATH=${INSTALL_PREFIX} \
-D parallel_arch:STRING=parscalar \
-D mpi_comms:BOOL=OFF \
-D QMP_DIR:PATH=${QMP_DIR} \
${SRCDIR}

# It may be necessary at this point to edit some of the cmake parameters by hand
# before running make
# ccmake ${SRCDIR}

${MAKE}
${MAKE} install

popd
