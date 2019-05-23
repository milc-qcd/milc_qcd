#! /bin/bash

# Usage

#  build.sh <scalar|knl|hsw|skx> <CXX>

# where CXX is the C++ compiler

# You may also need to edit milc-qphix/Makefile

ARCH=$1     # Choices: scalar, knl, hsw
PK_CXX=$2

#GIT_BRANCH=gauge_force
GIT_BRANCH=feature/fermion-force

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

# HACK: Need to copy files from milc-qphix/avx512 to milc-qphix-codegen/avx512
# Note, These files work only for double precision and avx512
#/bin/cp milc-qphix/avx512/ff_* ${dir}/avx512

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
mkdir -p lib

${MAKE} "ARCH=${ARCH}" "PK_CXX=${PK_CXX}"

popd

