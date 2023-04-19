#! /bin/bash

ARCH=$1
PK_CC=$2
PK_CXX=$3
#GIT_REPO=https://github.com/aportelli/Hadrons
GIT_REPO=https://github.com/milc-qcd/Hadrons
#GIT_BRANCH=develop
GIT_BRANCH=feature/staggered-a2a-ml

if [ -z ${PK_CXX} ]
then
  echo "Usage $0 <scalar|avx2|avx512-knl|avx512-skx|gpu-cuda|gpu-hip|gpu-sycl> <PK_CC> <PK_CXX>"
  exit 1
fi

case ${ARCH} in
    scalar|avx2|avx512-knl|avx512-skx|gpu-cuda|gpu-hip|gpu-sycl)
      ;;
    *)
      echo "Unsupported ARCH"
      echo "Usage $0 <scalar|avx2|avx512-knl|avx512-skx|gpu-cuda|gpu-hip|gpu-sycl> <PK_CC> <PK_CXX>"
      exit 1
esac

TOPDIR=`pwd`
SRCDIR=${TOPDIR}/Hadrons
BUILDDIR=${TOPDIR}/build-hadrons-${ARCH}
INSTALLDIR=${TOPDIR}/install-hadrons-${ARCH}
GRIDINSTALLDIR=${TOPDIR}/install-grid-${ARCH}

MAKE=make

if [ ! -d ${SRCDIR} ]
then
  echo "Fetching ${GIT_BRANCH} branch of Hadrons package from github"
  git clone ${GIT_REPO} -b ${GIT_BRANCH}
fi

# Initialization
pushd ${SRCDIR}
./bootstrap.sh
popd

# Configure only if not already configured                                          
mkdir -p ${BUILDDIR}
pushd ${BUILDDIR}
if [ ! -f Makefile ]
then
  echo "Configuring Hadrons for ${ARCH} in ${BUILDDIR}"

  case ${ARCH} in

    scalar)

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --with-grid=${GRIDINSTALLDIR} \
            CXX="${PK_CXX}" CC="${PK_CC}" \
            CXXFLAGS="-std=gnu++17 -O0 -g -fpermissive -Wno-psabi" \

       status=$?
             ;;

    avx2)

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --with-grid=${GRIDINSTALLDIR} \
            CXX="${PK_CXX}" CC="${PK_CC}" \
            CXXFLAGS="-std=c++17 -xCORE-AVX2" \

       status=$?
             ;;
    avx512-knl)

       INCMKL="-I/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include"
       LIBMKL="-L/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin"

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --with-grid=${GRIDINSTALLDIR} \
            --host=x86_64-unknown-linux-gnu \
            CXX="${PK_CXX}" CC="${PK_CC}" \
            CXXFLAGS="-std=c++17 -xMIC-AVX512 -DHAVE_LIME -O2 -g -simd -qopenmp" \

       status=$?
       echo "Configure exit status $status"
       ;;
    avx512-skx)

       INCMKL="-I/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include"
       LIBMKL="-L/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin"

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --with-grid=${GRIDINSTALLDIR} \
            --host=x86_64-unknown-linux-gnu \
            CXX="${PK_CXX}" CC="${PK_CC}" \
            CXXFLAGS="-std=c++17 -xCORE-AVX512 -O2 -g -simd -qopenmp" \

       status=$?
       echo "Configure exit status $status"
       ;;

    gpu-cuda)
	# Cori: salloc -C gpu -t 60 -N 1 -c 10 --gres=gpu:1 -A m1759
	${SRCDIR}/configure \
        --prefix ${INSTALLDIR}      \
        --with-grid=${GRIDINSTALLDIR} \
        --host=x86_64-unknown-linux-gnu

        status=$?
        echo "Configure exit status $status"
	;;

    gpu-hip)

	source ${TOPDIR}/env.sh
	# export PATH=/opt/rocm/bin:${PATH}

	${SRCDIR}/configure \
            --prefix ${INSTALLDIR}      \
            --with-grid=${GRIDINSTALLDIR} \
            --host=x86_64-unknown-linux-gnu \
	    CXXFLAGS="-std=c++17"

        status=$?
        echo "Configure exit status $status"
	;;

    gpu-sycl)

	${SRCDIR}/configure \
        --prefix ${INSTALLDIR}      \
        --with-grid=${GRIDINSTALLDIR} \
        --host=x86_64-unknown-linux-gnu

        status=$?
        echo "Configure exit status $status"
	;;

    *)
    echo "Unsupported ARCH ${ARCH}"
          exit 1;
  esac

  if [ $status -ne 0 ]
  then
      echo "Quitting because of configure errors"
  fi

fi

echo "Building in ${BUILDDIR}"
${MAKE} -k -j16

echo "Installing in ${INSTALLDIR}"
${MAKE} install

popd




