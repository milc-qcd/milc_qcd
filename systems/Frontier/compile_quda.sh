#!/bin/bash

# Run as:
#	bash compile_quda.sh

################### Setup environment for Frontier

module reset
module load PrgEnv-amd amd/5.3.0 rocm/5.3.0
module load craype-accel-amd-gfx90a
module load cray-mpich/8.1.28
module load cmake
module load perftools
module load ninja
module list

export PK_BUILD_TYPE="Release"
export MPICH_ROOT=${CRAY_MPICH_ROOTDIR}
export GTL_ROOT=${MPICH_ROOT}/gtl/lib
export MPICH_DIR=${MPICH_ROOT}/ofi/rocm-compiler/5.0
export PATH=${ROCM_PATH}/bin:${ROCM_PATH}/llvm/bin:${PATH}
export LD_LIBRARY_PATH=${INSTALLDIR}/lib:${ROCM_PATH}/llvm/lib64:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${ROCM_PATH}/llvm/lib:${MPICH_DIR}/lib:${GTL_ROOT}:${LD_LIBRARY_PATH}

MPI_CFLAGS="-I${MPICH_DIR}/include -g "
# Note the flags needed to enable XPMEM support when compiling with hipcc:
MPI_LDFLAGS="-g -Wl,-rpath=${MPICH_DIR}/lib -L${MPICH_DIR}/lib -lmpi -L${GTL_ROOT} -Wl,-rpath=${GTL_ROOT} -lmpi_gtl_hsa ${CRAY_XPMEM_POST_LINK_OPTS} -lxpmem ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a}"

MY_CFLAGS="$(pat_opts include hipcc gpu) $(pat_opts pre_compile hipcc gpu) ${MPI_CFLAGS} --offload-arch=gfx90a $(pat_opts post_compile hipcc gpu) -g -pg"
HIPFLAGS="$(pat_opts include hipcc gpu) $(pat_opts pre_compile hipcc gpu) --offload-arch=gfx90a $(pat_opts post_compile hipcc gpu)"
MY_LDFLAGS="$(pat_opts pre_link hipcc gpu) ${MPI_LDFLAGS} --offload-arch=gfx90a  $(pat_opts post_link hipcc gpu) -g -pg"

################### Get QUDA and prepare build/ directory

if [ -d quda ]
then
  cd quda
  git pull
  git checkout develop
else
  git clone https://github.com/lattice/quda
  cd quda
  git checkout develop
fi
cd ..

if [ -d build ]
then
  cd build
else
  mkdir build
  cd build
fi

################### Compile

cmake ../quda \
    -G "Ninja" \
    -DQUDA_TARGET_TYPE=HIP \
    -DQUDA_GPU_ARCH=gfx90a \
    -DROCM_PATH=${ROCM_PATH} \
    -DCMAKE_INSTALL_PREFIX=`pwd`/usqcd \
    -DCMAKE_BUILD_TYPE=RELEASE \
    -DQUDA_DIRAC_DEFAULT_OFF=ON \
    -DQUDA_DIRAC_STAGGERED=ON \
    -DQUDA_BUILD_SHAREDLIB=ON \
    -DQUDA_QMP=ON \
    -DQUDA_QIO=ON \
    -DQUDA_DOWNLOAD_USQCD=ON \
    -DQUDA_MULTIGRID=OFF \
    -DQUDA_SMEAR_GAUSS_TWOLINK=ON \
    -DCMAKE_BUILD_TYPE="DEVEL" \
    -DCMAKE_CXX_COMPILER="hipcc" \
    -DCMAKE_C_COMPILER="hipcc" \
    -DBUILD_SHARED_LIBS=ON \
    -DQUDA_BUILD_SHAREDLIB=ON \
    -DQUDA_BUILD_ALL_TESTS=ON \
    -DQUDA_CTEST_DISABLE_BENCHMARKS=ON \
    -DCMAKE_C_STANDARD=99 \
    -DCMAKE_CXX_FLAGS="${MY_CFLAGS}" \
    -DCMAKE_C_FLAGS="${MY_CFLAGS}" \
          -DCMAKE_HIP_FLAGS="${MY_CFLAGS}" \
    -DCMAKE_SHARED_LINKER_FLAGS="${MY_LDFLAGS}" \
    -DCMAKE_EXE_LINKER_FLAGS="${MY_LDFLAGS}" \

cmake --build . -j 16 -v 2>&1 | tee -a cmake_quda.log
echo "" | tee -a cmake_quda.out
cmake --install . 2>&1 | tee -a cmake_quda.log

