#! /bin/bash

ARCH=$1
PK_CC=$2
PK_CXX=$3
GIT_REPO=https://github.com/milc-qcd/Grid
GIT_BRANCH=develop

if [ -z ${PK_CXX} ]
then
  echo "Usage $0 <scalar|avx2|avx512-knl|avx512-skx|gpu-cuda> <PK_CC> <PK_CXX>"
  exit 1
fi

case ${ARCH} in
    scalar|avx512-knl|avx512-skx|avx2|gpu-cuda)
      ;;
    *)
      echo "Unsupported ARCH"
      echo "Usage $0 <scalar|avx2|avx512-knl|avx512-skx|gpu-cuda> <PK_CC> <PK_CXX>"
      exit 1
esac

TOPDIR=`pwd`
SRCDIR=${TOPDIR}/Grid
BUILDDIR=${TOPDIR}/build-grid-${ARCH}
INSTALLDIR=${TOPDIR}/install-grid-${ARCH}

MAKE="make V=1"

if [ ! -d ${SRCDIR} ]
then
  echo "Fetching ${GIT_BRANCH} branch of Grid package from github"
  git clone ${GIT_REPO} -b ${GIT_BRANCH}
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
            --enable-simd=GEN \
            --enable-comms=none \
	    --with-lime=${HOME}/scidac/install/qio-single \
	    --with-fftw=${HOME}/fftw/build-gcc \
            --with-openssl=/global/common/cori/software/openssl/1.1.0a/hsw \
            CXX="${PK_CXX}" \
            CXXFLAGS="-std=gnu++17 -Wno-psabi" \

# 	    --with-hdf5=/opt/cray/pe/hdf5/1.10.0/INTEL/15.0 \

       status=$?
             ;;

    avx2)

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --enable-mkl=yes \
            --enable-simd=GEN \
            --enable-comms=mpi \
	    --with-lime=${HOME}/scidac/install/qio-skx \
            --with-openssl=/global/common/cori/software/openssl/1.1.0a/hsw \
	    --with-hdf5=/opt/cray/pe/hdf5/1.12.0.0/INTEL/19.1 \
            CXX="${PK_CXX}" CC="${PK_CC}" \
            CXXFLAGS="-std=c++11 -xCORE-AVX2" \

       status=$?
             ;;
    avx512-knl)

       INCMKL="-I/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include"
       LIBMKL="-L/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin"

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --enable-simd=KNL \
            --enable-comms=mpi \
            --host=x86_64-unknown-linux-gnu \
	    --with-lime=${HOME}/scidac/install/qio-cori-extend-omp-knl-icc \
	    --with-hdf5=/opt/cray/pe/hdf5/1.12.0.0/INTEL/19.1 \
            --with-openssl=/global/common/cori/software/openssl/1.1.0a/hsw \
            CXX="${PK_CXX}" CC="${PK_CC}" \
            CXXFLAGS="-std=c++17 -xMIC-AVX512 -O2 -g -simd -qopenmp" \


       status=$?
       echo "Configure exit status $status"
       ;;
    avx512-skx)

       INCMKL="-I/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include"
       LIBMKL="-L/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin"

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --enable-simd=KNL \
            --enable-comms=mpi \
            --host=x86_64-unknown-linux-gnu \
	    --with-lime=${HOME}/scidac/install/qio-cori-omp-knl-icc \
            CXX="${PK_CXX}" CC="${PK_CC}" \
            CXXFLAGS="-std=c++17 -xCORE-AVX512 -O2 -g -simd -qopenmp" \

	    # --with-hdf5=/opt/cray/pe/hdf5/1.10.0.3/INTEL/16.0 \
            # --with-openssl=/global/common/cori/software/openssl/1.1.0a/hsw \

       status=$?
       echo "Configure exit status $status"
       ;;
    gpu-cuda)
	# Cori: salloc -C gpu -t 60 -N 1 -c 10 --gres=gpu:1 -A m1759
	${SRCDIR}/configure \
             --prefix ${INSTALLDIR}      \
	     --enable-comms=mpi          \
             --host=x86_64-unknown-linux-gnu \
             CXX=nvcc                    \
             LDFLAGS=-L$HOME/prefix/lib/ \
             CXXFLAGS="-ccbin ${PK_CXX} -gencode arch=compute_70,code=sm_70 -I$HOME/prefix/include/ -std=c++11" 
        status=$?
        echo "Configure exit status $status"
	;;

    #              --enable-simd=GEN           \

    *)
    echo "Unsupported ARCH ${ARCH}"
          exit 1;
  esac

  if [ $status -ne 0 ]
  then
      echo "Quitting because of configure errors"
  else
    echo "Building in ${BUILDDIR}"
    ${MAKE} -k -j20

    echo "Installing in ${INSTALLDIR}"
    ${MAKE} install
  fi

fi     
popd

# Might need to do these by hand...

# CayleyFermion5DInstantiationZWilsonImplF.cc
# g++-8 -DHAVE_CONFIG_H -I. -I/u/inscc/detar/milc_qcd/Grid/Grid/Grid    -I/home/falco/detar/milc/milc_qcd/Grid/Grid  -I/global/common/cori/software/openssl/1.1.0a/hsw/include -I/u/inscc/detar/scidac/install/qio-single/include -I/u/inscc/detar/fftw/build-gcc/include -fopenmp  -O3 -std=gnu++17 -Wno-psabi  -fno-strict-aliasing -c -o qcd/action/fermion/instantiation/ZWilsonImplF/CayleyFermion5DInstantiationZWilsonImplF.o /u/inscc/detar/milc_qcd/Grid/Grid/Grid/qcd/action/fermion/instantiation/ZWilsonImplF/CayleyFermion5DInstantiationZWilsonImplF.cc

# WilsonKernelsInstantiationWilsonImplDF.cc
# g++-8 -DHAVE_CONFIG_H -I. -I/u/inscc/detar/milc_qcd/Grid/Grid/Grid    -I/home/falco/detar/milc/milc_qcd/Grid/Grid  -I/global/common/cori/software/openssl/1.1.0a/hsw/include -I/u/inscc/detar/scidac/install/qio-single/include -I/u/inscc/detar/fftw/build-gcc/include -fopenmp  -O3 -std=gnu++17 -Wno-psabi  -fno-strict-aliasing -c -o qcd/action/fermion/instantiation/WilsonImplDF/WilsonKernelsInstantiationWilsonImplDF.o /u/inscc/detar/milc_qcd/Grid/Grid/Grid/qcd/action/fermion/instantiation/WilsonImplDF/WilsonKernelsInstantiationWilsonImplDF.cc

# CayleyFermion5DInstantiationGparityWilsonImplD.cc 
# g++-8 -DHAVE_CONFIG_H -I. -I/u/inscc/detar/milc_qcd/Grid/Grid/Grid    -I/home/falco/detar/milc/milc_qcd/Grid/Grid  -I/global/common/cori/software/openssl/1.1.0a/hsw/include -I/u/inscc/detar/scidac/install/qio-single/include -I/u/inscc/detar/fftw/build-gcc/include -fopenmp  -O3 -std=gnu++17 -Wno-psabi  -fno-strict-aliasing -c -o qcd/action/fermion/instantiation/GparityWilsonImplD/CayleyFermion5DInstantiationGparityWilsonImplD.o /u/inscc/detar/milc_qcd/Grid/Grid/Grid/qcd/action/fermion/instantiation/GparityWilsonImplD/CayleyFermion5DInstantiationGparityWilsonImplD.cc 

# CayleyFermion5DInstantiationZWilsonImplFH.cc
# g++-8 -DHAVE_CONFIG_H -I. -I/u/inscc/detar/milc_qcd/Grid/Grid/Grid    -I/home/falco/detar/milc/milc_qcd/Grid/Grid  -I/global/common/cori/software/openssl/1.1.0a/hsw/include -I/u/inscc/detar/scidac/install/qio-single/include -I/u/inscc/detar/fftw/build-gcc/include -fopenmp  -O3 -std=gnu++17 -Wno-psabi  -fno-strict-aliasing -c -o qcd/action/fermion/instantiation/ZWilsonImplFH/CayleyFermion5DInstantiationZWilsonImplFH.o /u/inscc/detar/milc_qcd/Grid/Grid/Grid/qcd/action/fermion/instantiation/ZWilsonImplFH/CayleyFermion5DInstantiationZWilsonImplFH.cc

# g++-8 -DHAVE_CONFIG_H -I. -I/u/inscc/detar/milc_qcd/Grid/Grid/Grid    -I/home/falco/detar/milc/milc_qcd/Grid/Grid  -I/global/common/cori/software/openssl/1.1.0a/hsw/include -I/u/inscc/detar/scidac/install/qio-single/include -I/u/inscc/detar/fftw/build-gcc/include -fopenmp  -O3 -std=gnu++17 -Wno-psabi  -fno-strict-aliasing  -c -o qcd/action/fermion/instantiation/WilsonImplF/CayleyFermion5DInstantiationWilsonImplF.o /u/inscc/detar/milc_qcd/Grid/Grid/Grid/qcd/action/fermion/instantiation/WilsonImplF/CayleyFermion5DInstantiationWilsonImplF.cc 

# g++-8 -DHAVE_CONFIG_H -I. -I/u/inscc/detar/milc_qcd/Grid/Grid/Grid    -I/home/falco/detar/milc/milc_qcd/Grid/Grid  -I/global/common/cori/software/openssl/1.1.0a/hsw/include -I/u/inscc/detar/scidac/install/qio-single/include -I/u/inscc/detar/fftw/build-gcc/include -fopenmp  -O3 -std=gnu++17 -Wno-psabi  -fno-strict-aliasing -c -o qcd/action/fermion/instantiation/GparityWilsonImplFH/CayleyFermion5DInstantiationGparityWilsonImplFH.o /u/inscc/detar/milc_qcd/Grid/Grid/Grid/qcd/action/fermion/instantiation/GparityWilsonImplFH/CayleyFermion5DInstantiationGparityWilsonImplFH.cc 

# g++-8 -DHAVE_CONFIG_H -I. -I/u/inscc/detar/milc_qcd/Grid/Grid/Grid    -I/home/falco/detar/milc/milc_qcd/Grid/Grid  -I/global/common/cori/software/openssl/1.1.0a/hsw/include -I/u/inscc/detar/scidac/install/qio-single/include -I/u/inscc/detar/fftw/build-gcc/include -fopenmp  -O3 -std=gnu++17 -Wno-psabi  -fno-strict-aliasing -c -o qcd/action/fermion/instantiation/GparityWilsonImplF/CayleyFermion5DInstantiationGparityWilsonImplF.o /u/inscc/detar/milc_qcd/Grid/Grid/Grid/qcd/action/fermion/instantiation/GparityWilsonImplF/CayleyFermion5DInstantiationGparityWilsonImplF.cc
