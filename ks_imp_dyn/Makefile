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
# 3. Precision 

# 1 = single precision; 2 = double
PRECISION = 1
#----------------------------------------------------------------------
# 4. Compiler
# Choices include mpicc cc gcc pgcc g++
CC = gcc

#CC = /usr/local/mpich/bin/mpicc -cc=${GCC_DIR}/bin/gcc  # FNAL
#CC = env GCC_EXEC_PREFIX=$(GCC_EXEC_PREFIX) powerpc-gnu-elf-gcc # QCDOC

# C++ compiler (needed only for QCDOC Level 3 wrapper)
CXX = env GCC_EXEC_PREFIX=$(GCC_EXEC_PREFIX) powerpc-gnu-elf-g++

#----------------------------------------------------------------------
# 5. Compiler optimization level
# Choices include -g -O, etc
OPT              = -O3

#----------------------------------------------------------------------
# 6. Other compiler optimization flags.  Uncomment stanza to suit.
#-------------- Gnu C -------------------------------------
#OCFLAGS = -Wall # ( -Wall, etc )

#OCFLAGS = -fexpensive-optimizations -funroll-loops -fpeephole -fstrength-reduce -fschedule-insns2 -fprefetch-loop-arrays # QCDOC
#OCFLAGS = -fexpensive-optimizations -fpeephole -fstrength-reduce -march=i586  # Simone's pick for PIII/gcc version 2.95.2.1 19991024 (release)
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
# 7. Choose large file support
#CLFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE64 # Large files gcc only
CLFS = # Not researched for others
#CLFS = -D_LARGE_FILES   # Blue Horizon

#----------------------------------------------------------------------
# 8. Installation-specific MPI includes and libraries
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
# 9. I/O routines
# Both io_nonansi and io_ansi should work on a scalar machine
# Solaris 2.6 gave "bad file number" errors with io_ansi.  CD
MACHINE_DEP_IO   = io_ansi.o # (io_ansi.o io_nonansi.o)

#----------------------------------------------------------------------
# 10. QDP/C options

# Edit these "wants"

# Choose "true" or blank. Does not imply any other SciDAC package
WANTQOP =

# Choose "true" or blank. Implies HAVEQIO and HAVEQMP.
WANTQDP = 

# Choose "true" or "". Implies HAVEQMP.
WANTQIO = true

# Choose "true" or "".
WANTQMP = true

#  Edit these locations for the SciDAC packages
# It is assumed that these are the parents of "include" and "lib"

SCIDAC = ${HOME}/scidac
QIOSNG = $(SCIDAC)/qio-single
QIOPAR = $(SCIDAC)/qio
# Parallel version
QMPPAR = ${SCIDAC}/qmp
# Single processor version
QMPSNG = ${SCIDAC}/qmp-single
QDP = ${SCIDAC}/qdp
QLA = ${SCIDAC}/qla
# Level 3
QOP = /home/chulwoo/SciDAC/Asqtad_L3

# Make_template_scidac defines these macros:
# HAVEQDP HAVEQIO HAVEQMP (Says what we are compiling with)
# LIBSCIDAC INCSCIDAC (The -L and -I compiler and linker lists)
# SCIDAC_LIBRARIES SCIDAC_HEADERS  (Lists for make dependencies)

include ../Make_template_scidac

#----------------------------------------------------------------------
# 11. Linker
LD               = ${CC}

#----------------------------------------------------------------------
# 12. Extra linker flags
#LDFLAGS = -Xlinker   # QCDOC
#LDFLAGS          = -fast     # Sun SPARC
#LDFLAGS          = -64 -L/usr/lib64 # SGIPC

#----------------------------------------------------------------------
# 13. Extra libraries
LIBADD =

#----------------------------------------------------------------------
# 14. Inlining choices

# SSE ASM and explicit C inlining is available for some of the library
# functions.

# Use SSE for P3 or P4 or Athlon Thunderbird and compilers, such
# as gcc that accept ASM macros

# While both SSE and C inline macros can be invoked by changing
# the function call to the macro name, with SSE_INLINE or C_INLINE
# these macros are automatically invoked where available

# See also the libraries Make_SSE_nasm for building non-inline SSE
# Some compilers don't like -DSSE_INLINE with the debugging -g option.

# Choose nothing, C_INLINE or SSE_INLINE, or both
INLINEOPT = -DC_INLINE # -DSSE_INLINE

#----------------------------------------------------------------------
# 15. Other miscellaneous macros you want for all of your compilations

CODETYPE =# -DQDP_PROFILE
# Choices include -DPREFETCH (not recommended)

#----------------------------------------------------------------------
# 16. Choose MILC library make file in libraries directory.  
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
CFLAGS = ${OPT} ${OCFLAGS} -D${COMMTYPE} ${CODETYPE} ${INLINEOPT} \
	${PREC} ${CLFS} ${INCSCIDAC} -I${MYINCLUDEDIR} ${DEFINES} ${DARCH}

ILIB = ${LIBSCIDAC} ${LMPI} ${LIBADD}

check:
	make -f Make_test check

test_clean:
	make -f Make_test test_clean

include Make_template
