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

#----------------------------------------------------------------------
# 3. Precision 

# 1 = single precision; 2 = double
PRECISION = 1

#----------------------------------------------------------------------
# 4. Compiler
# Choices include mpicc cc gcc pgcc g++
# Note: If you are linking with QUDA, you need the C++ linker,
# so you may as well compile and link everything, including the libraries, with C++

ifeq ($(strip ${MPP}),true)
  CC = mpicc
#  CC = mpixlc_r # BG/P
else
  CC = gcc
endif

#CC = /usr/local/mvapich/bin/mpicc  # FNAL
#----------------------------------------------------------------------
# 5. Compiler optimization level
# Choices include -g -O, etc

OPT              = -O3

#----------------------------------------------------------------------
# 6. Other compiler optimization flags.  Uncomment stanza to suit.
#-------------- Gnu C -------------------------------------
OCFLAGS = -Wall # ( -Wall, etc )

#OCFLAGS = -fexpensive-optimizations -fpeephole -fstrength-reduce -march=i586  # Simone's pick for PIII/gcc version 2.95.2.1 19991024 (release)
#OCFLAGS = -fexpensive-optimizations -funroll-loops -fpeephole -fstrength-reduce -fschedule-insns2 -march=i586 # works best for matrix x vector
#OCFLAGS =  -march=pentium4 -mfpmath=sse -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer # J. Osborn 10/20/04
#OCFLAGS =  -march=pentium4 -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer # J. Osborn 10/24/04

#------------------------ BlueGene -----------------------------------
# OCFLAGS = -qarch=450 -qlanglvl=stdc99 # BG/P
# OCFLAGS = -qarch=440d # BG/L

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

WANTQOP = # true # or blank. Implies HAVEQDP, HAVEQOP, HAVEQMP.

WANTQIO = true # or blank.  Implies HAVEQMP.

WANTQMP = # true or blank.

# Edit these locations for the installed SciDAC packages
# It is assumed that these are the parents of "include" and "lib"

SCIDAC = ${HOME}/scidac/install
# Parallel versions
QMPPAR = ${SCIDAC}/qmp
QIOPAR = $(SCIDAC)/qio
# Single processor versions
QMPSNG = ${SCIDAC}/qmp-single
QIOSNG = $(SCIDAC)/qio-single
QLA = ${SCIDAC}/qla
# Either version
QDP = ${SCIDAC}/qdp
QOPQDP = ${SCIDAC}/qopqdp
#QOPQDP = ${SCIDAC}/qopqdp-lapack # BG/P

QOP = ${QOPQDP}

# Make_template_scidac defines these macros:
# HAVEQOP HAVEQDP HAVEQIO HAVEQMP (Says what we are compiling with)
# LIBSCIDAC INCSCIDAC (The -L and -I compiler and linker lists)
# SCIDAC_LIBRARIES SCIDAC_HEADERS  (Lists for make dependencies)
# CSCIDAC (List of compiler macros for SciDAC modules)

include ../Make_template_scidac

#----------------------------------------------------------------------
# 11. FFTW3 Options

WANTFFTW = #true

ifeq ($(strip ${WANTFFTW}),true)
FFTW=${HOME}/fftw/build-gcc
INCFFTW = -I${FFTW}/include
LIBFFTW = -L${FFTW}/lib

ifeq ($(strip ${PRECISION}),1)
  LIBFFTW += -lfftw3f
else
  LIBFFTW += -lfftw3
endif
endif

#----------------------------------------------------------------------
# 12. LAPACK Options (for qopqdp-lapack and arb_overlap )

#LIBLAPACK = -L/opt/ibmcmp/xlf/bg/11.1/lib /soft/apps/LAPACK/liblapack_bgp.a /soft/apps/LIBGOTO/libgoto.a -lxlf90 -lxlsmp # LAPACK on BG/P

LIBLAPACK = # -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -llapack -lblas -lgfortran

#----------------------------------------------------------------------
# 13. PRIMME Options (for arb_overlap).  REQUIRES LAPACK AS WELL.

WANTPRIMME = #true

ifeq ($(strip ${WANTPRIMME}),true)
  LIBPRIMME = -L${HOME}/milc/install/PRIMME -lzprimme
endif

#----------------------------------------------------------------------
# 14. GPU/QUDA Options

WANTQUDA    = #true
WANT_CL_BCG_GPU = #true
WANT_FN_CG_GPU = #true
WANT_FL_GPU = #true
WANT_FF_GPU = #true
WANT_GF_GPU = #true

ifeq ($(strip ${WANTQUDA}),true)

  QUDA_HOME = ${HOME}/quda
  INCQUDA = -I${QUDA_HOME}/include -I/lib -I${QUDA_HOME}/tests
  LIBQUDA = -L${QUDA_HOME}/lib -lquda
  QUDA_LIBRARIES = ${QUDA_HOME}/lib

  CUDA_HOME = /usr/local/cuda
  INCQUDA += -I${CUDA_HOME}/include
  LIBQUDA += -L${CUDA_HOME}/lib64 -lcudart

# Definitions of compiler macros -- don't change.  Could go into a Make_template_QUDA

  ifeq ($(strip ${WANT_CL_BCG_GPU}),true)
    HAVE_CL_GPU = true
    CGPU += -DUSE_CL_GPU
  endif

  ifeq ($(strip ${WANT_FN_CG_GPU}),true)
    HAVE_FN_CG_GPU = true
    CGPU += -DUSE_CG_GPU
  endif

  ifeq ($(strip ${WANT_GF_GPU}),true)
    HAVE_GF_GPU = true
    CGPU += -DUSE_GF_GPU
  endif

  ifeq ($(strip ${WANT_FL_GPU}),true)
    HAVE_FL_GPU = true
    CGPU += -DUSE_FL_GPU
  endif

  ifeq ($(strip ${WANT_FF_GPU}),true)
    HAVE_FF_GPU = true
    CGPU += -DUSE_FF_GPU
  endif

endif

#----------------------------------------------------------------------
# 15. Linker
LD               = ${CC}

#----------------------------------------------------------------------
# 16. Extra ld flags

#LDFLAGS          = -fast     # Sun SPARC
#LDFLAGS          = -64 -L/usr/lib64 # SGIPC

#----------------------------------------------------------------------
# 17. Extra include paths
INCADD = ${INCFFTW} ${INCQUDA} 

#----------------------------------------------------------------------
# 18. Extra libraries
LIBADD = ${LIBFFTW} ${LIBPRIMME} ${LIBLAPACK} ${LIBQUDA}

#----------------------------------------------------------------------
# 19. Inlining choices

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
INLINEOPT = -DC_GLOBAL_INLINE # -DSSE_GLOBAL_INLINE #-DC_INLINE

# There are special single-precision macros for the AMD Opteron
# To get them, uncomment the next line

#INLINEOPT += -DSSEOPTERON

#----------------------------------------------------------------------
# 20. Miscellaneous macros for performance control and metric

#     Define them with a -D prefix.

#------------------------------
# Print timing statistics.
# Applications: many

# Use any combination of these
# CGTIME CG Solver
# FFTIME Fermion force
# FLTIME Link fattening
# GFTIME Gauge force
# IOTIME I/O timing
# PRTIME print time (clover_invert2)

# REMAP  report remapping time for QDP, QOP in conjunction with above

CTIME = # -DCGTIME -DFFTIME -DFLTIME -DGFTIME -DREMAP -DPRTIME -DIOTIME

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
# Debugging and diagnostics
# Applications:  all

# CHECK_MALLOC        Report malloc/free activity.
#                     (Then process stdout using check_malloc.pl)
# CG_DEBUG            Print debugging information for the inverters.
# CG_OK               Print inverter convergence information even when OK
# REMAP_STDIO_APPEND  All nodes append to stdout.
#
# HISQ_SVD_VALUES_INFO Print HISQ SVD diagnostics
# HISQ_SVD_COUNTER    Print summary count of SVD uses
# HISQ_FORCE_FILTER_COUNTER Print summary count of force filter applications.

CDEBUG = -DCG_OK # -DCHECK_MALLOC -DREMAP_STDIO_APPEND

#------------------------------
# Backward compatibility

# As of version qopqdp 0.9.0 the normalization convention for the
# staggered inverter changed.  If you are using a version of QOPQDP
# with the old convention, define this macro:
CCOMPAT += #-DOLD_QOPQDP_NORM

# Prior to version 7.7.2 the conversion from staggeredd to naive was peculiar.
CCOMPAT += #-DOLD_STAGGERED2NAIVE

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
# Applications: ks_imp_rhmc ks_measure ks_spectrum

# Choices
# KS_MULTICG=OFFSET  The basic multicg solver.
# KS_MULTICG=HYBRID  Solve with multicg and polish off with single mass CG.
# KS_MULTICG=FAKE    Iterate the single mass solver.
# KS_MULTICG=REVERSE Iterate in reverse order
# KS_MULTICG=REVHYB  Same as HYBRID but with vectors in reverse order.

# HALF_MIXED         If PRECISION=2, do multimass solve in single precision
#                    and single-mass refinements in double
# NO_REFINE          No refinements except for masses with nonzero Naik eps

KSCGMULTI = -DKS_MULTICG=HYBRID -DHALF_MIXED # -DNO_REFINE

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

# Staggered spin-taste operator shift for applications that construct
# nonlocal interpolating operators.  By default the shift operator is
# symmetric (forward and backward).  This macro makes it one-sided
# (forward).

# The Fat-Naik variants of the nonlocal operators are currently
# constructed only from the forward fat links.

KSSHIFT = # -DONE_SIDED_SHIFT

#------------------------------
# Clover inverter choice
# Applications: clover_invert, clover_invert2
#

# CL_CG=BICG  Biconjugate gradient
# CL_CG=CG    Standard CG
# CL_CG=MR    Minimum residue
# CL_CG=HOP   Hopping

# HALF_MIXED  Do double-precision inversion with single, or single with half (if supported)
# MAX_MIXED   Do double-precision inversion with half-precision (if supported)
# SCALE_PROP  Do rescaling for the clover propagator

CLCG = -DCL_CG=BICG 

#------------------------------
# Propagator storage
# Applications: clover_invert2
#

# CLOV_LEAN   Write intermediate propagators to disk, saving memory
#             but increasing run time somewhat

CLMEM = #-DCLOV_LEAN

#------------------------------
# Summary

CODETYPE = ${CTIME} ${CPROF} ${CDEBUG} ${CGEOM} ${KSCGSTORE} ${CPREFETCH} \
 ${KSCGMULTI} ${KSFFMULTI} ${KSRHMCINT} ${KSSHIFT} ${CLCG} ${CLMEM} ${CQOP} \
 ${CCOMPAT}

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

ifeq ($(strip ${WANTFFTW}),true)
  HAVEFFTW = true
endif

ifeq ($(strip ${WANTPRIMME}),true)
  HAVEPRIMME = true
endif

# Make_template_combos defines convenience macros for interdependent
# groups of compilation units.  They are used to specify build lists.

include ../Make_template_combos

CPREC = -DPRECISION=${PRECISION} ${QDPPREC} ${QOPPREC}
DARCH = ${CSCIDAC} ${CGPU}

# Complete set of compiler flags - do not change
CFLAGS = ${OPT} ${OCFLAGS} -D${COMMTYPE} ${CODETYPE} ${INLINEOPT} \
	${CPREC} ${CLFS} ${INCSCIDAC} -I${MYINCLUDEDIR} ${DARCH} \
	${DEFINES} ${ADDDEFINES} ${IMPI} ${INCADD}

ILIB = ${LIBSCIDAC} ${LMPI} ${LIBADD}

.PHONY: time check test_clean
time:
	make -f Make_time time

check: test_clean
	cd test ; perl ../../check.pl ${EXEC} ${CASE} ${PREC} < checklist

test_clean:
	cd test ; make test_clean

include Make_template
