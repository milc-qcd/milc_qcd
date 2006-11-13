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

OPT              = -O3 -Wall

#----------------------------------------------------------------------
# 6. Other compiler optimization flags.  Uncomment stanza to suit.
#-------------- Gnu C -------------------------------------
#OCFLAGS = -Wall # ( -Wall, etc )

#OCFLAGS = -fexpensive-optimizations -funroll-loops -fpeephole -fstrength-reduce -fschedule-insns2 -fprefetch-loop-arrays # QCDOC
#OCFLAGS = -fexpensive-optimizations -fpeephole -fstrength-reduce -march=i586  # Simone's pick for PIII/gcc version 2.95.2.1 19991024 (release)
#OCFLAGS = -fexpensive-optimizations -funroll-loops -fpeephole -fstrength-reduce -fschedule-insns2 -march=i586 # works best for matrix x vector
#OCFLAGS =  -march=pentium4 -mfpmath=sse -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer # J. Osborn 10/20/04
#OCFLAGS =  -march=pentium4 -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer # J. Osborn 10/24/04

#------------------------ BlueGene -----------------------------------
# OCFLAGS = -qarch=440d -qtune=440

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
#CLFS = -D_LARGE_FILES   # AIX

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

# Choose the QOP package.  
# Choices are 
#    QCDOC for Level 3 QCDCOC
#    QDP   for Level 3 QOPQDP
#    MILC (nonoptimized MILC implementation for testing)
#    blank if you don't want QOP
WANTQOP = 

# Choose "true" or blank. Implies HAVEQIO and HAVEQMP.
WANTQDP = 

# Choose "true" or "". Implies HAVEQMP.
WANTQIO = 

# Choose "true" or "".
WANTQMP = 

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
QOP = ${SCIDAC}/qop
QOPQDP = ${SCIDAC}/qopqdp

# Make_template_qop defines these macros:
# HAVEQOP
# LIBQOP INCQOP
# INCDEPQOP LIBDEPQOP
# GENERICQOP ASQINVERTQOP ASQFORCEQOP

include ../Make_template_qop

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

# USE INLINE SSE WITH EXTREME CAUTION!  IT MAY GIVE WRONG RESULTS.

# SSE ASM and explicit C inlining is available for some of the library
# functions.

# Use SSE for P3 or P4 or Athlon Thunderbird and compilers, such
# as gcc that accept ASM macros

# There is also a little-tested alternate single precision SSE package in the
# directory sse_opteron.  To use it, you need to rename the
# sse directory to sse_p4 (or something) and then rename the
# sse_opteron directory to "sse".

# Both SSE and C inline macros can be invoked selectively by defining
# SSE_INLINE and C_INLINE and changing the function call to the macro
# name.

# To invoke the inline macros globally (where available) without
# changing the code, define, instead, SSE_GLOBAL_INLINE or
# C_GLOBAL_INLINE.

# You may use SSE and C inlining at the same time.  The SSE version
# takes precedence when both are available.

# See also the libraries Make_SSE_nasm for building non-inline SSE
# Some compilers don't like -DSSE_INLINE with the debugging -g option.

# Choose nothing, C_INLINE or SSE_INLINE, or both
INLINEOPT = -DC_GLOBAL_INLINE # -DSSE_GLOBAL_INLINE -DC_INLINE

#----------------------------------------------------------------------
# 15. Miscellaneous macros for performance control and metric

#     Define them with a -D prefix.

#------------------------------
# Print timing statistics.

# Use any combination of these
# CGTIME CG Solver
# FFTIME Fermion force
# LLTIME Link fattening
# GFTIME Gauge forde

CTIME =#

#------------------------------
# Profiling

# QDP_PROFILE         Generates a report for all QDP routines

CPROF =#

#------------------------------
# Layout

# These are currently selected only by editing Make_template.
#  Choices
#    layout_hyper_prime.o          Standard hypercubic
#    layout_timeslices.o           Puts full timeslices on each node.
#    layout_squares.o              Puts 2D slices on each node.
#    layout_hyper_tstretch.o       Rarely used.  Fewer timeslices on last node.
#                                  In case node no is incommensurate with nt
#    layout_timeslices_2.o         Untested.  Supposedly more forgiving.

#------------------------------
# Improved staggered CG inverter and Dslash

#  Note, some options still require editing Make_template
#  Choices 
#   dslash_fn.o                   Overlaps computation with backward gathers.
#   dslash_fn2.o                  Does all gathers before computation
#                                   but has fewer FORALLSITES loops
#   dslash_fn_dblstore.o          Double store version of dslash_fn.o
#                                   Supports GATHER13

# DBLSTORE_FN  Copies backward links.  Requires more memory.
#              You must also call for dslash_fn_dblstore.o for now.
# GATHER13     Combine third and next neighbor gathers in Dslash.
#              For now, works only with dslash_fn_dblstore.o
# FEWSUMS      Fewer CG reductions

KSCGSTORE =#

#------------------------------
# Staggered fermion force routines

# These are currently selected only by editing Make_template

# Choices
#  fermion_force_asqtad3.o   Optimized for the Asqtad action
#  fermion_force_general.o   Takes any quark action

#------------------------------
# Prefetching

# PREFETCH (Not working yet)

CPREFETCH = #

#------------------------------
# Multimass improved KS CG solvers

# KS_MULTICG=OFFSET  The basic multicg solver.
# KS_MULTICG=HYBRID  Solve with multicg and polish off with single mass CG.
# KS_MULTICG=FAKE    Iterate the single mass solver.
# KS_MULTICG=REVERSE Iterate in reverse order
# KS_MULTICG=REVHYB  Same as HYBRID but with vectors in reverse order.

KSCGMULTI =#

#------------------------------
# Multifermion force routines

# KSMULTIFF=FNMAT    Construct matrix parallel transporters
#                    and use with sum of outer products of sources 
# KSMULTIFF=FNMATREV Older version of FNMAT.  Traverses path
#                    in reverse order
# KSMULTIFF=ASVEC    Use improved Asqtad, parallel transporting
#                    groups of source vectors.  See VECLENGTH.
# VECLENGTH=n        Number of source vectors to process in one group.
#                    Applies only to the ASVEC option

KSFFMULTI =#

#------------------------------
# Summary

CODETYPE = ${CTIME} ${CPROF} ${KSCGSTORE} ${CPREFETCH} ${KSCGMULTI}\
 ${KSFFMULTI}

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
  QDPPREC = -DQDP_Precision=${PRECISION}
endif

ifeq ($(strip ${HAVEQOP}),true)
  QOPPREC = -DQOP_Precision=${PRECISION}
endif

PREC = -DPRECISION=${PRECISION} ${QDPPREC} ${QOPPREC}

# Complete set of compiler flags - do not change
CFLAGS = ${OPT} ${OCFLAGS} -D${COMMTYPE} ${CODETYPE} ${INLINEOPT} \
	${PREC} ${CLFS} ${INCSCIDAC} -I${MYINCLUDEDIR} ${DEFINES} ${DARCH}

ILIB = ${LIBSCIDAC} ${LMPI} ${LIBADD}

check:
	make -f Make_test check

test_clean:
	make -f Make_test test_clean

include Make_template
