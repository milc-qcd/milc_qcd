#  MIMD version 7
#  Standard application makefile.
#  Replaces Make_vanilla, Make_linux_mpi, Make_linux_qdp, Make_qcdoc_gcc
#  Do not use for making libraries
#  Copy this file into the application directory and run make there.
#

# This file
MAKEFILE = Makefile

#----------------------------------------------------------------------
#  User choices - edit to suit 
#----------------------------------------------------------------------
# 0. Shell (if it really matters)

#SHELL = /bin/bash

#----------------------------------------------------------------------
# 1. Architecture

# Compiling for a parallel machine?  blank for a scalar machine
#MPP = true

# Cross-compiling for the QCDOC?  blank if we are not.
# For the QCDOC, be sure to source the QOS setup script before running make,
# and be sure to set MPP = true
#QCDOC = true

#----------------------------------------------------------------------
# 2. Compiler
CC = gcc   # ( mpicc cc gcc pgcc g++ )

#CC = /usr/local/mpich/bin/mpicc -cc=${GCC_DIR}/bin/gcc  # FNAL
#CC = env GCC_EXEC_PREFIX=$(GCC_EXEC_PREFIX) powerpc-gnu-elf-gcc # QCDOC

#----------------------------------------------------------------------
# 3. Compiler optimization level
OPT              = -O3  # ( -g -O, etc )

#----------------------------------------------------------------------
# 4. Other compiler optimization flags.  Uncomment stanza to suit.
#-------------- Gnu C -------------------------------------
#OCFLAGS = -Wall   # ( -Wall, etc )

#OCFLAGS = -fexpensive-optimizations -fpeephole -fstrength-reduce -march=i586  # Simone's pick for PIII/gcc version 2.95.2.1 19991024 (release)
#OCFLAGS = -fexpensive-optimizations -funroll-loops -fpeephole -fstrength-reduce -fschedule-insns2 -march=i586 # works best for matrix x vector
#OCFLAGS =  -march=pentium4 -mfpmath=sse -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer # J. Osborn 10/20/04
#OCFLAGS =  -march=pentium4 -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer # J. Osborn 10/24/04
#-------------- Intel icc/ecc -----------------------------------
#OCFLAGS = -tpp2 -static

#-------------- Portland Group ----------------------------
#OCFLAGS = -tp p6 -Munroll=c:4,n:4
#OCFLAGS= -mpentiumpro -march=pentiumpro -funroll-all-loops -malign-double -D_REENTRANT  # Pentium pro

#-------------- TCS alpha -------------------------------------
#OCFLAGS = -float -arch=ev68 # ( -arch=ev67 )

#-------------- SUN SPARC ---------------------------------
#OCFLAGS= -fast -dalign -xlibmil -fsimple=2 -fns  #Ultra

#-------------- SGI Origin -------------------------------------
#OCFLAGS = -64 -mips4 -OPT:IEEE_arithmetic=3:roundoff=3:alias=restrict -TENV:X=1
#OCFLAGS = -64 -mips4 -r10000 -OPT:IEEE_arithmetic=3:roundoff=3:alias=restrict -TENV:X=1

#-------------- Blue Horizon -------------------------------------
#CARCH = -qarch=pwr3 -qtune=pwr3   # Architecture: ( ppc pwr2 )
# For Blue Horizon (memory flags for 1 GByte)
#OCFLAGS= ${CARCH} -Q=500 -qmaxmem=-1 -bmaxdata:0x40000000 -bmaxstack:0x8000000

#----------------------------------------------------------------------
# 5. Choose large file support
#CLFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE64 # Large files gcc only
CLFS =             # Not researched for others
#CLFS = -D_LARGE_FILES   # Blue Horizon

#----------------------------------------------------------------------
# 6. Installation-specific MPI includes and libraries
#    Not needed if using mpicc or on single processor

#----------------- MPICH/GM ---------------------------------------------
#IMPI = #-I/usr/local/mpich/include
#IMPI = -I/usr/include   # Pittsburgh TCS. Really!

#LMPI = #-L/usr/local/mpich/lib/shared -L/usr/local/mpich/lib -lmpich -L/opt/gm/lib/ -lgm
#LMPI = --lfmpi -lmpi -lelan -lelan3 -lrmscall -lmach # Pittsburgh TCS

#----------------- MPIPRO ---------------------------------------------
#IMPI =
# With Redhat, for MPI PRO see rpm -ql mpipro | more
#LMPI = -lmpipro_tv -lpthread
#LMPI = -lmpipro -lvipl -lpthread

#----------------- MVICH ----------------------------------------------
#IMPI = -I/uufs/icebox/sys/src/mpich/1.2.0-via/include  # MVICH
#LMPI = -L/uufs/icebox/sys/pkg/mpich/1.2.0-via/lib -lmpi -lvipl -lpthread # MVICH

#----------------------------------------------------------------------
# 7. I/O routines
# Both io_nonansi and io_ansi should work on a scalar machine
# Solaris 2.6 gave "bad file number" errors with io_ansi.  CD
MACHINE_DEP_IO   = io_ansi.o # (io_ansi.o io_nonansi.o)

#----------------------------------------------------------------------
# 8. QDP/C options

# Edit these "wants"

# Choose "true" or blank. Implies HAVEQIO and HAVEQMP.
WANTQDP =

# Choose "true" or "". Implies HAVEQMP.
WANTQIO =

# Choose "true" or "".
WANTQMP =

#  Edit these locations for the SciDAC packages
# It is assumed that these are the parents of "include" and "lib"

SCIDAC = ${HOME}/scidac
QIO = $(SCIDAC)/qio
# Parallel version
QMPPAR = ${SCIDAC}/qmp
# Single processor version
QMPSNG = ${SCIDAC}/qmp-single
QDP = ${SCIDAC}/qdp
QLA = ${SCIDAC}/qla

# Do not change these conditionals

HAVEQDP = ${WANTQDP}

ifeq ($(strip ${HAVEQDP}),true)
   HAVEQIO = true #
else
   HAVEQIO = ${WANTQIO}
endif

ifeq ($(strip ${HAVEQIO}),true)
   HAVEQMP = true
else
   HAVEQMP = ${WANTQMP}
endif

# Shouldn't need to change the rest of this stanza (except see "QCDOC
# nonstandard")

ifeq ($(strip ${MPP}),true)
   QMP = ${QMPPAR}
else
   QMP = ${QMPSNG}
endif

IQMP = -I${QMP}/include
LQMP = -L$(QMP)/lib -lqmp
#LQMP = -L$(QMP)/lib -lqcd_api     # QCDOC nonstandard
#LQMP = -L$(QMP)/lib -lqcd_api_nn -lscu -lscu_devel # QCDOC nonstandard

IQIO = -I${QIO}/include
LQIO = -L${QIO}/lib -lqio -llime 

IQLA = -I${QLA}/include
LQLA = -L${QLA}/lib -lqla_int -lqla_f -lqla_f3 -lqla_df -lqla_d3 -lqla_df3 \
          -lqla_dq3 -lqla_q3 -lqla_d -lqla_dq -lqla_q -lqla_random

IQDP = -I${QDP}/include
LQDP = -L${QDP}/lib -lqdp_int -lqdp_f -lqdp_f3 -lqdp_d -lqdp_d3 -lqdp_common

ifeq ($(strip ${HAVEQMP}),true)
  LIBQMP = ${LQMP}
  INCQMP = ${IQMP}
endif

ifeq ($(strip ${HAVEQIO}),true)
  LIBQIO = ${LQIO}
  INCQIO = ${IQIO}
endif

ifeq ($(strip ${HAVEQDP}),true)
  LIBQDP = ${LQDP} ${LQLA}
  INCQDP = ${IQDP} ${IQLA}
endif

LIBSCIDAC = ${LIBQDP} ${LIBQIO} ${LIBQMP}
INCSCIDAC = ${INCQDP} ${INCQIO} ${INCQMP}

#----------------------------------------------------------------------
# 9. Linker
LD               = ${CC}

#----------------------------------------------------------------------
# 10. Extra linker flags
LDFLAGS          =   # most
#LDFLAGS          = -fast     # Sun SPARC
#LDFLAGS          =   -L$(QOS)/quser/gcc-lib-user/// -Xlinker   # QCDOC
#LDFLAGS          = -64 -L/usr/lib64 # SGIPC

#----------------------------------------------------------------------
# 11. Extra libraries
LIBADD = 

#----------------------------------------------------------------------
# 12. Precision 

# 1 = single precision; 2 = double
PRECISION = 2

#----------------------------------------------------------------------
# 13. SSE

# SSE for P3 or P4 or Athlon Thunderbird and compliant compilers
# -DSSE says we have the SSE package.  It does not do global inlining
# but it allows selective inlining through explicit inline macro calls.
# -DSSE_INLINE together with -DSSE gives automatic global inlining.
# See also the libraries Make_SSE_nasm for building non-inline SSE
# Don't use -DSSE together with Make_SSE_nasm.

SSEOPT = -DSSE # -DSSE_INLINE

#----------------------------------------------------------------------
# 14. Other miscellaneous macros you want for all of your compilations

CODETYPE  = # -DQDP_PROFILE
# Choices include -DPREFETCH (not recommended)

#----------------------------------------------------------------------
# 15. Choose MILC library make file in libraries directory.  
#    CHECK IT FOR FURTHER OPTIONS!

MAKELIBRARIES = Make_vanilla  # or Make_SSE_nasm (don't use -DSSE with this)
# Other choices Make_RS6K  Make_alpha

#----------------------------------------------------------------------
# End of user choices.  Please, also, check choices in include/config.h.
#----------------------------------------------------------------------

ifeq ($(strip ${MPP}),true)
  ifeq ($(strip ${HAVEQMP}),true)
     COMMTYPE = QMP
     COMMPKG = com_qmp.o
  else
     COMMTYPE = MPI
     COMMPKG = com_mpi.o
  endif
else
  ifeq ($(strip ${HAVEQMP}),true)
     COMMTYPE = SINGLE
     COMMPKG = com_qmp.o
  else
     COMMTYPE = SINGLE
     COMMPKG = com_vanilla.o
  endif
endif

ifeq ($(strip ${HAVEQDP}),true)
  PREC = -DPRECISION=${PRECISION} -DQDP_Precision=${PRECISION}
else
  PREC = -DPRECISION=${PRECISION}
endif

# Complete set of compiler flags - do not change
CFLAGS = ${OPT} ${OCFLAGS} -D${COMMTYPE} ${CODETYPE} ${SSEOPT} \
	${PREC} ${CLFS} ${INCSCIDAC} -I${MYINCLUDEDIR} ${DEFINES}

ILIB = ${LIBSCIDAC} ${LMPI} ${LIBADD}

include Make_template
