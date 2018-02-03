#! /bin/bash

ARCH=$1     # Choices: scalar, avx2, avx512
SOALEN=$2   # Choices: 4, 8 or 1 (for scalar)
PK_CC=$3
PK_CXX=$4
GIT_BRANCH=devel

MAKE="make -j4"

TOPDIR=`pwd`
SRCDIR=${TOPDIR}/qphix

if [ ! -d qphix ]
then
  echo "Fetching ${GIT_BRANCH} branch of package from github"
  git clone https://github.com/JeffersonLab/qphix/ -b ${GIT_BRANCH}
  rm -f qphix/configure qphix/Makefile  # (Will be rebuilt)

  cd qphix
  autoreconf -vif
  cd ..
fi

mkdir -p build/dslash-${ARCH}-s${SOALEN} install/dslash-${ARCH}-s${SOALEN}
pushd build/dslash-${ARCH}-s${SOALEN}

# Configure only if not already configured
if [ ! -f Makefile ]
then

  case ${ARCH} in
  
    scalar)
  
  ${SRCDIR}/configure \
  	--prefix=${TOPDIR}/install/dslash-scalar-s${SOALEN} \
          --with-qmp=${HOME}/scidac/install/qmp-single \
  	--enable-proc= \
  	--enable-soalen=1 \
  	--enable-clover \
          --enable-openmp \
  	--enable-cean \
  	--enable-mm-malloc \
  	--enable-parallel-arch=scalar \
          CXXFLAGS="${PK_CXXFLAGS}" \
          CFLAGS="${PK_CFLAGS}" \
  	CXX="${PK_CXX}" \
  	CC="${PK_CC}" \
  	--host=x86_64-linux-gnu --build=none \
  
        ;;
  
    avx2)
  
  PK_CXXFLAGS="-openmp -g -O2 -finline-functions -fno-alias -std=c++11 -xCORE-AVX2 -vec-report -restrict"
  PK_CFLAGS="-openmp -g  -O2 -fno-alias -std=c99 -xCORE-AVX2 -vec-report -restrict"
  
  ${SRCDIR}/configure \
  	--prefix=${TOPDIR}/install/dslash-avx2-s${SOALEN} \
          --with-qmp=${HOME}/scidac/install/qmp-cori-icc \
  	--enable-proc=AVX2 \
  	--enable-soalen=${SOALEN} \
  	--enable-clover \
          --enable-openmp \
  	--enable-cean \
  	--enable-mm-malloc \
  	--enable-parallel-arch=parscalar \
          CXXFLAGS="${PK_CXXFLAGS}" \
          CFLAGS="${PK_CFLAGS}" \
  	CXX="${PK_CXX}" \
  	CC="${PK_CC}" \
  	--host=x86_64-linux-gnu --build=none \
  
        ;;
  
    avx512)
  
  PK_CXXFLAGS="-openmp -g -O2 -finline-functions -fno-alias -std=c++11 -xMIC-AVX512 -vec-report -restrict"
  PK_CFLAGS="-openmp -g  -O2 -fno-alias -std=c99 -xMIC-AVX512 -vec-report -restrict" \
  ${SRCDIR}/configure \
  	--prefix=${TOPDIR}/install/dslash-avx512-s${SOALEN} \
          --with-qmp=${HOME}/scidac/install/qmp-cori-omp-knl-icc \
  	--enable-proc=AVX512 \
  	--enable-soalen=${SOALEN} \
  	--enable-clover \
          --enable-openmp \
  	--enable-cean \
  	--enable-mm-malloc \
  	--enable-parallel-arch=parscalar \
          CXXFLAGS="${PK_CXXFLAGS}" \
          CFLAGS="${PK_CFLAGS}" \
  	CXX="${PK_CXX}" \
  	CC="${PK_CC}" \
  	--host=x86_64-linux-gnu --build=none \
  
        ;;
  
    *)
  	echo "Unsupported ARCH ${ARCH}"
          exit 1;
  esac

fi # If ! -f Makefile

${MAKE}
${MAKE} install

popd
