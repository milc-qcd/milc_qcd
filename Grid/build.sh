#! /bin/bash

ARCH=$1
PK_CC=$2
PK_CXX=$3
GIT_BRANCH=feature/block-cg

MAKE="make -j4"

TOPDIR=`pwd` 
dir=Grid
SRCDIR=${TOPDIR}/${dir}

MAKE="make -j4"

if [ ! -d ${dir} ]
then
  echo "Fetching ${GIT_BRANCH} branch of package from github"
  git clone https://github.com/paboyle/Grid -b ${GIT_BRANCH}
  rm -f ${dir}/Makefile  # (Will be rebuilt)                                   

  cd ${dir}
  autoreconf -vif
  cd ..
fi

# Configure only if not already configured                                          
pushd ${dir}
if [ ! -f Makefile ]
then

  case ${ARCH} in

    scalar)

       ${SRCDIR}/configure \
            --prefix=`pwd` \
            --disable-openmp \
            CXX="${PK_CXX}" \
            CC="${PK_CC}" \
             ;;

    avx2)

       ${SRCDIR}/configure \
            --prefix=`pwd` \
            --enable-mkl=yes \
            CXX="${PK_CXX}" \
            CC="${PK_CC}" \
             ;;
    avx512)

       ${SRCDIR}/configure \
            --prefix=`pwd` \
            --enable-mkl=yes \
            CXX="${PK_CXX}" \
            CC="${PK_CC}" \
             ;;
    *)
    echo "Unsupported ARCH ${ARCH}"
          exit 1;
  esac

fi     
popd




