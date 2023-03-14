#! /bin/bash

ARCH=$1
PK_CC=$2
PK_CXX=$3
GIT_REPO=https://github.com/milc-qcd/Grid
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
SRCDIR=${TOPDIR}/Grid
BUILDDIR=${TOPDIR}/build-grid-${ARCH}
INSTALLDIR=${TOPDIR}/install-grid-${ARCH}

MAKE=make

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
	    --with-lime=${HOME}/crusher/quda/install/qio \
	    --with-fftw=${HOME}/fftw/build-gcc \
            --with-mpfr=${HOME}/mpfr \
            CXX="${PK_CXX}" \
            CXXFLAGS="-std=gnu++17 -O0 -g -Wno-psabi" \

#            --with-openssl=/global/common/cori/software/openssl/1.1.0a/hsw \
# 	    --with-hdf5=/opt/cray/pe/hdf5/1.10.0/INTEL/15.0 \
#            --disable-gparity \
#	    --disable-zmobius \
#	    --disable-fermion-reps \

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
	# Summit: ./build-Grid.sh gpu-cuda mpicc mpiCC
	# Perlmutter ./build-Grid.sh gpu-cuda cc CC
	${SRCDIR}/configure \
             --prefix ${INSTALLDIR}       \
	     --enable-comms=mpi           \
	     --enable-comms-threads       \
	     --enable-simd=GPU            \
	     --enable-shm=shmnone         \
             --enable-gen-simd-width=64   \
	     --enable-accelerator=cuda    \
	     --disable-fermion-reps       \
	     --enable-unified             \
	     --disable-gparity            \
             --host=x86_64-unknown-linux-gnu \
	     --with-mpfr=${HOME}/perlmutter/mpfr \
	     --with-hdf5=${HOME}/perlmutter/hdf5 \
	     --with-lime=${HOME}/perlmutter/build/usqcd \
             CXX="nvcc"                \
	     LDFLAGS="-cudart shared " \
             CXXFLAGS="-ccbin ${PK_CXX} -gencode arch=compute_80,code=sm_80 -std=c++14 -cudart shared" \

        status=$?
        echo "Configure exit status $status"
	;;


    gpu-hip)

	source ${TOPDIR}/env.sh

	${SRCDIR}/configure \
	     --prefix ${INSTALLDIR}      \
	     --disable-fermion-reps      \
             --disable-gparity \
             --disable-zmobius \
	     --enable-accelerator=hip \
	     --enable-comms=mpi3-auto \
	     --enable-gen-simd-width=64 \
	     --enable-shm=nvlink \
	     --enable-simd=GPU \
	     --enable-tracing=timer \
             --enable-unified=no \
	     --with-fftw=${FFTW_DIR}/.. \
	     --with-gmp=${OLCF_GMP_ROOT} \
	     --with-hdf5=${OLCF_HDF5_ROOT} \
 	     --with-lime=${INSTALLROOT}/qio \
	     --with-mpfr=/opt/cray/pe/gcc/mpfr/3.1.4/ \
	     CXX=hipcc    CXXLD=hipcc \
	     MPICXX=${MPICH_DIR}/bin/mpicxx \
	     CXXFLAGS="${MPI_CFLAGS} -I${ROCM_PATH}/include -std=c++14 -O3 -fPIC -fopenmp --amdgpu-target=gfx90a" \
	     LDFLAGS="-L/lib64 -L${ROCM_PATH}/lib -lamdhip64 ${MPI_LDFLAGS}" \

        status=$?
        echo "Configure exit status $status"
	cp ${BUILDDIR}/grid-config ${INSTALLDIR}/bin

	;;

    gpu-sycl)

	# ./build-Grid.sh gpu-sycl dpcpp dpcpp


	${SRCDIR}/configure \
	 --prefix ${INSTALLDIR}      \
	 --enable-simd=GPU \
	 --enable-comms=mpi \
	 --enable-gen-simd-width=64  \
         --disable-gparity \
         --disable-zmobius \
         --disable-fermion-reps \
         --enable-accelerator=sycl   \
	 --enable-unified=yes \
	 CXXCPP="/soft/packaging/spack-builds/linux-opensuse_leap15-x86_64/gcc-10.2.0/gcc-10.2.0-yudlyezca7twgd5o3wkkraur7wdbngdn/bin/cpp" \
         CXX="${PK_CXX}" CC="${PK_CC}" \
	 CXXFLAGS="-cxx=dpcpp -fsycl-unnamed-lambda -fsycl -no-fma -std=c++17 -O0 -g" \
	 LDFLAGS="-fsycl-device-code-split=per_kernel -fsycl-device-lib=all" \

	 
#	 CXXFLAGS="-cxx=dpcpp -fsycl-unnamed-lambda -fsycl -no-fma -std=c++17" \

	 #	     --enable-comms=mpi          \
#	     --with-lime=${HOME}/scidac/install/qio-gcc \

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
${MAKE} -k -j20

echo "Installing in ${INSTALLDIR}"
${MAKE} install

popd


