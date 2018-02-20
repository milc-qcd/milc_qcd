#! /bin/bash

ARCH=$1     # Choices: scalar, avx2, avx512
SOALEN=$2   # Choices: 4, 8 or 1 (for scalar)
PK_CC=$3
PK_CXX=$4
GIT_BRANCH=devel

MAKE="make -j4"

TOPDIR=`pwd`
SRCDIR=${TOPDIR}/qphix

ARCH=$1     # Choices: scalar, avx2, avx512
SOALEN=$2   # Choices: 4, 8 or 1 (for scalar)
PK_CC=$3
PK_CXX=$4
GIT_BRANCH=devel

MAKE="make -j4"

TOPDIR=`pwd`
SRCDIR=${TOPDIR}/qphix

if [ ! -d ${SRCDIR} ]
then
  echo "Fetching ${GIT_BRANCH} branch of package from github"
  git clone https://github.com/JeffersonLab/qphix/ -b ${GIT_BRANCH}
fi

pushd build
cmake ${SRCDIR}

sed -e 's/isa:STRING=.*/isa:STRING='${ARCH}'/ CMakeCache.txt > foo
mv foo CMakeCache.txt

${MAKE}
${MAKE} install

popd
