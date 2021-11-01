#! /bin/bash

# Usage

#  build.sh <scalar|knl|hsw|skx> <CC> <CXX>

# where CC is the C compiler (currently ignored)
#       CXX is the C++ compiler

# You may also need to edit milc-qphix/Makefile_qphixlib

ARCH=$1     # Choices: scalar, knl, hsw
PK_CC=$2
PK_CXX=$3
GIT_BRANCH=master

if [ -z "${ARCH}" ]
then
   echo "Usage $0 <scalar|knl|skx|hsw> <CC> <CXX>"
   exit 1
fi

if [ ${ARCH} = "scalar" ]
then
  modecmd="mode=scalar"
  TARGET=scalar
elif [ ${ARCH} = "knl" ]
then
  modecmd="mode=mic"
  TARGET=avx512
elif [ ${ARCH} = "skx" ]
then
  modecmd="mode=mic"
  TARGET=avx512
elif [ ${ARCH} = "hsw" ]
then
  modecmd="mode=hsw"
  TARGET=avx512
else
  modecmd="mode=mic"
  TARGET=avx2
fi
   
###############################
#  Build milc-qphix-codegen
###############################

dir=milc-qphix-codegen
MAKE="make -j4"

if [ ! -d ${dir} ]
then
  echo "Fetching ${GIT_BRANCH} branch of package from github"
  git clone https://github.com/JeffersonLab/${dir} -b ${GIT_BRANCH}
fi

pushd ${dir}

${MAKE} ${TARGET} ${modecmd}

popd


###############################
#  Build milc-qphix
###############################

dir=milc-qphix
MAKE="make -j4"

if [ ! -d ${dir} ]
then
  echo "Fetching ${GIT_BRANCH} branch of package from github"
  git clone https://github.com/JeffersonLab/${dir} -b ${GIT_BRANCH}
fi

pushd ${dir}

${MAKE} -f Makefile_qphixlib ${modecmd} "ARCH=${ARCH}" "PK_CXX=${PK_CXX}"

popd

