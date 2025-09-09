#!/bin/sh

# Run after running compile_qioqmp.sh
# Note that the ks_spectrum application can be run without QIO and QMP. If so,
# remove the relevant lines from the make command.
# Run as:
#	bash compile_ks_spectrum_hisq.sh

pushd .

module reset
module load gompi
module list

echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

# Paths to QMP and QIO installations
QIOPAR=`pwd`/build/qio/qio
QMPPAR=`pwd`/build/qmp/qmp

# Download/update MILC
if [ -d milc_qcd ]
then
  cd milc_qcd/milc_qcd
  git checkout develop
  git pull
else
  mkdir milc_qcd
  cd milc_qcd
  git clone https://github.com/milc-qcd/milc_qcd.git
  cd milc_qcd
  git checkout develop
fi

############ Make ks_spectrum_hisq ##################
cd ks_spectrum
cp ../Makefile .
make clean

MY_CC=cc \
MY_CXX=CC \
ARCH="" \
COMPILER="gnu" \
OPT="-O3 -Ofast -g" \
LDFLAGS="-g -lmpi" \
PRECISION=2 \
MPP=true \
OMP=true \
WANTQIO=true \
WANTQMP=true \
QIOPAR=${QIOPAR} \
QMPPAR=${QMPPAR} \
CGEOM="-DFIX_NODE_GEOM -DFIX_IONODE_GEOM" \
KSCGMULTI="-DKS_MULTICG=HYBRID " \
CTIME="-DNERSC_TIME -DCGTIME -DFFTIME -DFLTIME -DGFTIME -DREMAP -DPRTIME -DIOTIME -DGS_TIME" \
make -j 1 ks_spectrum_hisq >& make_ks_spectrum_hisq.log

popd

echo ""
echo ""
echo "Check that compilation was successful by viewing make_ks_baryon.log:"
echo ""
tail milc_qcd/milc_qcd/ks_spectrum/make_ks_spectrum_hisq.log

