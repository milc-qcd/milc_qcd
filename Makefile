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
# 1. Shell (if it really matters)

#SHELL = /bin/bash

#----------------------------------------------------------------------
# 2. Architecture

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

ifeq ($(strip ${MPP}),true)
  CC = /usr/local/mvapich/bin/mpicc
else
  CC = gcc
endif

#CC = /usr/local/mvapich/bin/mpicc  # FNAL
#CC = powerpc-gnu-elf-gcc           # QCDOC
#----------------------------------------------------------------------
# 5. Compiler optimization level
# Choices include -g -O, etc

OPT              = -O3

#----------------------------------------------------------------------
# 6. Other compiler optimization flags.  Uncomment stanza to suit.
#-------------- Gnu C -------------------------------------
OCFLAGS = -Wall # ( -Wall, etc )

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
# 7. Choose large file support.

CLFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE # Large files gcc only
#CLFS = # Not researched for others
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

MACHINE_DEP_IO   = io_ansi.o # (io_ansi.o io_nonansi.o io_dcap.o)

# Uncomment if you have and want to support dcache I/O
# (Forces use of io_dcap.o)

# WANTDCAP = true

# The location of the installed dcap libraries. Uncomment and enter if
# it is not already defined as a system environment variable.

# DCAP_DIR = /usr/local/develop/dcache

# Choose the appropriate library path according to the addressing
# size of the machine

# DCAPLIB  = lib64 # (lib64 lib)

#----------------------------------------------------------------------
# 10. SciDAC package options

# Edit these "wants"

# Choose the QOP package.  
# Choices are 
#    QCDOC for Level 3 QCDOC
#    QDP   for Level 3 QOPQDP
#    MILC (nonoptimized MILC implementation for testing)
#    blank if you don't want QOP

WANTQOP = #QDP (or QCDOC or MILC)

# Backward compatibility for QOPQDP:
# As of version qopqdp 0.9.0 the normalization convention for the
# staggered inverter changed.  If you are using a version of QOPQDP
# with the old convention, define this macro:
CQOPQDP_NORM = #-DOLD_QOPQDP_NORM

# Choose "true" or blank. Implies HAVEQIO and HAVEQMP.
WANTQDP = 

# Choose "true" or "". Implies HAVEQMP.
WANTQIO = 

# Choose "true" or "".
WANTQMP = 


# Edit these locations for the installed SciDAC packages
# It is assumed that these are the parents of "include" and "lib"

SCIDAC = ${HOME}/scidac
# Parallel versions
QMPPAR = ${SCIDAC}/qmp
QIOPAR = $(SCIDAC)/qio
# Single processor versions
QMPSNG = ${SCIDAC}/qmp-single
QIOSNG = $(SCIDAC)/qio-single
QLA = ${SCIDAC}/qla
# Either version
QDP = ${SCIDAC}/qdp-single
QOPQDP = ${SCIDAC}/qopqdp-single

QOP = ${QOPQDP}

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

# Choose nothing or
#  [ C_INLINE | C_GLOBAL_INLINE ] [ SSE_INLINE | SSE_GLOBAL_INLINE ]
INLINEOPT = -DC_GLOBAL_INLINE -DSSE_GLOBAL_INLINE #-DC_INLINE

# There are special single-precision macros for the AMD Opteron
# To get them, uncomment the next line

#INLINEOPT += -DSSEOPTERON

#----------------------------------------------------------------------
# 15. Miscellaneous macros for performance control and metric

#     Define them with a -D prefix.

#------------------------------
# Print timing statistics.
# Applications: many

# Use any combination of these
# CGTIME CG Solver
# FFTIME Fermion force
# LLTIME Link fattening
# GFTIME Gauge force

# REMAP  report remapping time for QDP, QOP in conjunction with above

CTIME =# -DCGTIME -DFFTIME -DLLTIME -DGFTIME -DREMAP

#------------------------------
# Profiling
# Applications:  QDP

# QDP_PROFILE         Generates a report for all QDP routines

CPROF =#

#------------------------------
# Troubleshooting
# Applications: All

# COM_CRC            Message passing test.  Checksums on all gathers.

#------------------------------
# Debugging
# Applications:  all

# CHECK_MALLOC        Report malloc/free activity.
#                     Process output using check_malloc.pl
# CG_DEBUG            Print debugging information for the inverters.
#
# MILC_GLOBAL_DEBUG   So far applies only to ks_imp_rhmc HISQ code

CDEBUG =#

#------------------------------
# Layout
# Applications: all

# These are currently selected only by editing Make_template.
#  Choices
#    layout_hyper_prime.o          Standard hypercubic
#    layout_timeslices.o           Puts full timeslices on each node.
#    layout_squares.o              Puts 2D slices on each node.
#    layout_hyper_tstretch.o       Rarely used.  Fewer timeslices on last node.
#                                  In case node no is incommensurate with nt
#    layout_timeslices_2.o         Untested.  Supposedly more forgiving.

#------------------------------
# Compute node grid layout

#     layout_hyper_prime views the machine as a grid and selects a layout
#     that distributes the lattice uniformly across the nodes and
#     minimizes the surface to volume ratio.  On fixed grid machines
#     such as the IBM Bluegene it may be better to control the
#     geometry.  In that case some of our applications support
#     the compiler macro FIX_NODE_GEOM.  You must then provide
#     and extra list of dimensions in the parameter input file.
#     See e.g. ks_imp_rhmc.

CGEOM =#-DFIX_NODE_GEOM

#------------------------------
# I/O node grid layout

#     For some applications and architectures it is more efficient to
#     split large files and have a subset of processors handle the I/O.
#     We define the I/O node partitions by hypercubes in the same manner 
#     as the compute node geometry above.  The I/O node geometry must
#     be commensurate with the compute node geometry.
#     For applications that support it, the I/O node geometry is specified
#     by the macro FIX_IONODE_GEOM.  Then the parameter input file
#     includes a list of dimensions.

CGEOM +=# -DFIX_IONODE_GEOM

#------------------------------
# Improved staggered CG inverter and Dslash
# Applications: ks_imp_dyn ks_imp_rhmc ks_imp_invert_multi ks_hl_spectrum

#  Note, some options still require editing Make_template
#  Choices 
#   dslash_fn.o                   Overlaps computation with backward gathers.
#   dslash_fn2.o                  Does all gathers before computation
#                                   but has fewer FORALLSITES loops
#   dslash_fn_dblstore.o          Double store version of dslash_fn.o
#                                   Supports GATHER13

# DBLSTORE_FN    Copies backward links.  Requires more memory.
#                You must also call for dslash_fn_dblstore.o for now.
# D_FN_GATHER13  Combine third and next neighbor gathers in Dslash.
#                For now, works only with dslash_fn_dblstore.o
# FEWSUMS        Fewer CG reduction calls

KSCGSTORE = -DDBLSTORE_FN -DD_FN_GATHER13 -DFEWSUMS

#------------------------------
# Staggered fermion force routines
# Applications: ks_imp_dyn ks_imp_rhmc

# These are currently selected only by editing Make_template

# Choices
#  fermion_force_asqtad.o    Optimized for the Asqtad action
#  fermion_force_general.o   Takes any quark action

#------------------------------
# Prefetching
# Applications: all

# PREFETCH (Not working yet)

CPREFETCH = #

#------------------------------
# Multimass improved KS CG solvers
# Applications: ks_imp_rhmc ks_imp_invert_multi

# Choices
# KS_MULTICG=OFFSET  The basic multicg solver.
# KS_MULTICG=HYBRID  Solve with multicg and polish off with single mass CG.
# KS_MULTICG=FAKE    Iterate the single mass solver.
# KS_MULTICG=REVERSE Iterate in reverse order
# KS_MULTICG=REVHYB  Same as HYBRID but with vectors in reverse order.

KSCGMULTI = -DKS_MULTICG=HYBRID

#------------------------------
# Multifermion force routines
# Applications: ks_imp_rhmc

# Choices
# KS_MULTIFF=FNMAT    Construct matrix parallel transporters
#                     and use with sum of outer products of sources 
# KS_MULTIFF=FNMATREV Older version of FNMAT.  Traverses path
#                     in reverse order
# KS_MULTIFF=ASVEC    Use improved Asqtad, parallel transporting
#                     groups of source vectors.  See VECLENGTH.

# Additional options
# VECLENGTH=n        Number of source vectors to process in one group.
#                    Applies only to the ASVEC option

KSFFMULTI = -DKS_MULTIFF=FNMAT


#------------------------------
# RHMC molecular dynamics algorithm
# Applications: ks_imp_rhmc

# Choices
# (Not needed if the make target has a preset value in Make_template)

# INT_ALG=INT_LEAPFROG
# INT_ALG=INT_OMELYAN
# INT_ALG=INT_2EPS_3TO1
# INT_ALG=INT_2EPS_2TO1
# INT_ALG=INT_2G1F
# INT_ALG=INT_3G1F
# INT_ALG=INT_4MN4FP
# INT_ALG=INT_4MN5FV
# INT_ALG=INT_FOURSTEP
# INT_ALG=INT_PLAY

KSRHMCINT =#

#------------------------------
# Clover inverter choice
# Applications: clover_invert
#

# CL_CG=BICG  Biconjugate gradient
# CL_CG=CG    Standard CG
# CL_CG=MR    Minimum residue
# CL_CG=HOP   Hopping

CLCG = #-DCL_CG=BICG 

#------------------------------
# Summary

CODETYPE = ${CTIME} ${CPROF} ${CDEBUG} ${CGEOM} ${KSCGSTORE} ${CPREFETCH} \
 ${KSCGMULTI} ${KSFFMULTI} ${KSRHMCINT} ${CLCG} ${CQOP} ${CQOPQDP_NORM}

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

ifeq ($(strip ${WANTDCAP}),true)
   MACHINE_DEP_IO = io_dcap.o
   OCFLAGS += -I${DCAP_DIR}/include
   LDFLAGS += -L${DCAP_DIR}/${DCAPLIB} -Wl,--rpath,${DCAP_DIR}/${DCAPLIB} -ldcap
endif

CPREC = -DPRECISION=${PRECISION} ${QDPPREC} ${QOPPREC}

# Complete set of compiler flags - do not change
CFLAGS = ${OPT} ${OCFLAGS} -D${COMMTYPE} ${CODETYPE} ${INLINEOPT} \
	${CPREC} ${CLFS} ${INCSCIDAC} -I${MYINCLUDEDIR} ${DEFINES} ${DARCH} \
	${IMPI}

ILIB = ${LIBSCIDAC} ${LMPI} ${LIBADD}

check:
	make -f Make_test check

test_clean:
	make -f Make_test test_clean

include Make_template
