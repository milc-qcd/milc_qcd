#! /bin/bash

ARCH=$1
PK_CC=$2
PK_CXX=$3
GIT_BRANCH=develop

if [ -z ${PK_CXX} ]
then
  echo "Usage $0 <scalar|avx512|avx2> <CC> <CXX>"
  exit 1
fi

case ${ARCH} in
    scalar|avx512|avx2)
      ;;
    *)
      echo "Unsupported ARCH"
      echo "Usage $0 <scalar|avx512|avx2> <CC> <CXX>"
      exit 1
esac

MAKE="make -j4"

TOPDIR=`pwd`
SRCDIR=${TOPDIR}/Grid
BUILDDIR=${TOPDIR}/build-${ARCH}
INSTALLDIR=${TOPDIR}/install-${ARCH}

MAKE="make -j4"

if [ ! -d ${SRCDIR} ]
then
  echo "Fetching ${GIT_BRANCH} branch of Grid package from github"
  git clone https://github.com/paboyle/Grid -b ${GIT_BRANCH}
fi

# Fetch Eigen package, set up Make.inc files and create Grid configure
pushd ${SRCDIR}
./bootstrap.sh
popd

# Configure only if not already configured                                          
mkdir -p ${BUILDDIR}
pushd ${BUILDDIR}
if [ ! -f Makefile ]
then
  echo "Configuring Grid for ${ARCH} in ${BUILDDIR}"

  case ${ARCH} in

    scalar)

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --disable-openmp \
            --enable-precision=double \
            --enable-simd=GEN \
            --enable-comms=none \
            CXX="${PK_CXX}" \
            CC="${PK_CC}" \
             ;;

    avx2)

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --enable-mkl=yes \
            --enable-precision=double \
            --enable-simd=GEN \
            --enable-comms=mpi \
            CXX="${PK_CXX}" \
            CXXFLAGS="-std=c++11 -xMIC-AVX512" \

             ;;
    avx512)

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --enable-precision=double \
            --enable-simd=GEN \
            --enable-comms=mpi \
            --host=x86_64-unknown-linux-gnu \
	    --with-hdf5=/opt/cray/pe/hdf5/1.10.0/INTEL/15.0 \
	    --with-lime=${HOME}/scidac/install/qio-cori-omp-knl-icc/lib \
            CXX="${PK_CXX}" \
            CXXFLAGS="-std=c++11 -xMIC-AVX512" \


	    #            --enable-mkl=yes \

             ;;
    *)
    echo "Unsupported ARCH ${ARCH}"
          exit 1;
  esac

  echo "Building in ${BUILDDIR}"
  ${MAKE}

  echo "Installing in ${INSTALLDIR}"
  ${MAKE} install

fi     
popd




