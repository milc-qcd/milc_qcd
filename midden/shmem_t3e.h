/* USMID @(#)libsma/include/shmem.h	20.5	05/16/96 08:52:25 */

/*
 *      (C) COPYRIGHT CRAY RESEARCH, INC.
 *      UNPUBLISHED PROPRIETARY INFORMATION.
 *      ALL RIGHTS RESERVED.
 */

/*
 * 	This file contains user definitions for the shmem routines.
 */

#ifndef _MPP_SHMEM_H
#define _MPP_SHMEM_H

#include <sys/cdefs.h>

int	__shmem_h_version_1;		/* identify version of header file */
#ifdef _CRAY1
#pragma _CRI common __shmem_h_version_1
#endif

#ifdef	_CRAYMPP
#include <mpp/mpphw.h>
#endif

#ifdef	_CRAY1
#include <mpt.h>	/* bring in defs for _my_pe() and _num_pes() */
#endif


#if	defined(_RELEASE) || defined(_CRAYC)	/* CRI SCC or C++ compiler */

#ifdef __cplusplus				/* C++ compiler */
#include <intrinsics.h>
#endif 

/**#include "../include/complex.h"**/
#endif

#ifdef __cplusplus				/* C++ compiler */
/**#define	_SHMEM_COMPLEXD	**/	/* double_complex */ complex
#elif	defined(_RELEASE) || defined(_CRAYC)	/* CRI SCC compiler */
/**#define	_SHMEM_COMPLEXF		Real complex **/
/**#define	_SHMEM_COMPLEXD		double complex**/
#endif


/*
 *	CONSTANTS
 */

#ifdef _CRAYT3E

#define _CACHELINE 8
#define _LOG2MAXPE 15			/* assume NPES <= 32768 */

#define _SHMEM_BCAST_SYNC_SIZE  	((2+_LOG2MAXPE)*_CACHELINE)
#define _SHMEM_BARRIER_SYNC_SIZE        ((1+_LOG2MAXPE)*3)
#define _SHMEM_REDUCE_SYNC_SIZE		((2+_LOG2MAXPE)*_CACHELINE)
#define _SHMEM_COLLECT_SYNC_SIZE	((2+_LOG2MAXPE)*_CACHELINE)
#define _SHMEM_SYNC_VALUE       	(-1)
#define _SHMEM_REDUCE_MIN_WRKDATA_SIZE  8

#else

#define _SHMEM_BCAST_SYNC_SIZE          128
#define _SHMEM_BARRIER_SYNC_SIZE        36
#define _SHMEM_REDUCE_SYNC_SIZE         36
#define _SHMEM_COLLECT_SYNC_SIZE        36
#define _SHMEM_SYNC_VALUE               (-1)
#define _SHMEM_REDUCE_MIN_WRKDATA_SIZE  8

#endif

#define	_SHMEM_CMP_EQ	0
#define	_SHMEM_CMP_NE	1
#define	_SHMEM_CMP_GT	2
#define	_SHMEM_CMP_LE	3
#define	_SHMEM_CMP_LT	4
#define	_SHMEM_CMP_GE	5

/*
 *	PROTOTYPES: USER FUNCTIONS
 */
__BEGIN_DECLS 

extern void shmem_get(void *, void *, int, int);
extern void shmem_put(void *, void *, int, int);
extern void shmem_iget(void *, void *, int, int, int, int);
extern void shmem_iput(void *, void *, int, int, int, int);
extern void shmem_ixget(void *, void *, long *, int, int);
extern void shmem_ixput(void *, void *, long *, int, int);

extern void shmem_get32(void *, void *, int, int);
extern void shmem_put32(void *, void *, int, int);
extern void shmem_iget32(void *, void *, int, int, int, int);
extern void shmem_iput32(void *, void *, int, int, int, int);
extern void shmem_ixget32(void *, void *, short *, int, int);
extern void shmem_ixput32(void *, void *, short *, int, int);

extern void shmem_short_p(short *addr, short value, int pe);
extern void shmem_int_p(int *addr, int value, int pe);
extern void shmem_long_p(long *addr, long value, int pe);
extern void shmem_Real_p(Real *addr, Real value, int pe);
extern void shmem_double_p(double *addr, double value, int pe);

extern short shmem_short_g(short *addr, int pe);
extern int shmem_int_g(int *addr, int pe);
extern long shmem_long_g(long *addr, int pe);
extern Real shmem_Real_g(Real *addr, int pe);
extern double shmem_double_g(double *addr, int pe);

extern long   shmem_swap       (long   *addr, long   value, int pe);
extern short  shmem_short_swap (short  *addr, short  value, int pe);
extern int    shmem_int_swap   (int    *addr, int    value, int pe);
extern long   shmem_long_swap  (long   *addr, long   value, int pe);
extern Real  shmem_Real_swap (Real  *addr, Real  value, int pe);
extern double shmem_double_swap(double *addr, double value, int pe);

extern short  shmem_short_mswap (short *addr, short mask, short value, int pe);
extern int    shmem_int_mswap   (int   *addr, int   mask, int   value, int pe);
extern long   shmem_long_mswap  (long  *addr, long  mask, long  value, int pe);

extern short  shmem_short_cswap (short *addr, short match, short value, int pe);
extern int    shmem_int_cswap   (int   *addr, int   match, int   value, int pe);
extern long   shmem_long_cswap  (long  *addr, long  match, long  value, int pe);

extern short  shmem_short_finc (short  *addr, int pe);
extern short  shmem_short_fadd (short  *addr, short value, int pe);

extern void shmem_set_lock(long *);
extern void shmem_clear_lock(long *);
extern int shmem_test_lock(long *);

extern void shmem_barrier(int, int, int, long *);
extern void shmem_wait(long *, long);
extern void shmem_wait_until(long *, int, long);

extern void shmem_set_cache_inv(void);
extern void shmem_set_cache_line_inv(long *);
extern void shmem_clear_cache_inv(void);
extern void shmem_udcflush(void);
extern void shmem_udcflush_line(long *);

extern int shmem_n_pes(void);
extern int shmem_my_pe(void);
extern void shmem_stack(void*);
extern void shmem_fence(void);
extern void shmem_quiet(void);

extern void shmem_int_sum_to_all(int *, int *, int, int, int, int, int *,
	long *);
extern void shmem_int_max_to_all(int *, int *, int, int, int, int, int *,
	long *);
extern void shmem_int_min_to_all(int *, int *, int, int, int, int, int *,
	long *);
extern void shmem_int_prod_to_all(int *, int *, int, int, int, int, int *,
	long *);
extern void shmem_int_and_to_all(int *, int *, int, int, int, int, int *,
	long *);
extern void shmem_int_or_to_all(int *, int *, int, int, int, int, int *,
	long *);
extern void shmem_int_xor_to_all(int *, int *, int, int, int, int, int *,
	long *);

extern void shmem_short_sum_to_all(short *, short *, int, int, int, int,
	short *, long *);
extern void shmem_short_max_to_all(short *, short *, int, int, int, int,
	short *, long *);
extern void shmem_short_min_to_all(short *, short *, int, int, int, int,
	short *, long *);
extern void shmem_short_prod_to_all(short *, short *, int, int, int, int,
	short *, long *);
extern void shmem_short_and_to_all(short *, short *, int, int, int, int,
	short *, long *);
extern void shmem_short_or_to_all(short *, short *, int, int, int, int,
	short *, long *);
extern void shmem_short_xor_to_all(short *, short *, int, int, int, int,
	short *, long *);

extern void shmem_double_sum_to_all(double *, double *, int, int, int, int,
	double *, long *);
extern void shmem_double_max_to_all(double *, double *, int, int, int, int,
	double *, long *);
extern void shmem_double_min_to_all(double *, double *, int, int, int, int,
	double *, long *);
extern void shmem_double_prod_to_all(double *, double *, int, int, int, int,
	double *, long *);

extern void shmem_Real_sum_to_all(Real *, Real *, int, int, int, int,
	Real *, long *);
extern void shmem_Real_max_to_all(Real *, Real *, int, int, int, int,
	Real *, long *);
extern void shmem_Real_min_to_all(Real *, Real *, int, int, int, int,
	Real *, long *);
extern void shmem_Real_prod_to_all(Real *, Real *, int, int, int, int,
	Real *, long *);

#ifdef _SHMEM_COMPLEXF
extern void shmem_complexf_sum_to_all(_SHMEM_COMPLEXF *, _SHMEM_COMPLEXF *,
	int, int, int, int, _SHMEM_COMPLEXF *, long *);
extern void shmem_complexf_prod_to_all(_SHMEM_COMPLEXF *, _SHMEM_COMPLEXF *,
	int, int, int, int, _SHMEM_COMPLEXF *, long *);
#endif	/* _SHMEM_COMPLEXF */

#ifdef _SHMEM_COMPLEXD
extern void shmem_complexd_sum_to_all(_SHMEM_COMPLEXD *, _SHMEM_COMPLEXD *,
	int, int, int, int, _SHMEM_COMPLEXD *, long *);
extern void shmem_complexd_prod_to_all(_SHMEM_COMPLEXD *, _SHMEM_COMPLEXD *,
	int, int, int, int, _SHMEM_COMPLEXD *, long *);
#endif	/* _SHMEM_COMPLEXD */

extern void shmem_collect(long *, long *, int, int, int, int, long *);
extern void shmem_fcollect(long *, long *, int, int, int, int, long *);

extern void shmem_collect32(short *, short *, int, int, int, int, long *);
extern void shmem_fcollect32(short *, short *, int, int, int, int, long *);

extern void shmem_broadcast(void *, void *, int, int, int, int, int,
	long *);
extern void shmem_broadcast32(void *, void *, int, int, int, int, int,
	long *);


/*
 *	PROTOTYPES: CRI-INTERNAL
 *
 *	Here are prototypes for CRI-internal versions of the above.  The
 *	names which begin with an underscrore provide a means to avoid 
 *	forbidden collisions with the C user namespace.
 */
extern void _shmem_get(void *, void *, int, int);
extern void _F_shmem_get(void *, void *, int, int);
extern void _shmem_put(void *, void *, int, int);
extern void _F_shmem_put(void *, void *, int, int);
extern void _shmem_iget(void *, void *, int, int, int, int);
extern void _shmem_iput(void *, void *, int, int, int, int);
extern void _shmem_ixget(void *, void *, long *, int, int);
extern void _shmem_ixput(void *, void *, long *, int, int);
 
extern void _shmem_get32(void *, void *, int, int);
extern void _shmem_put32(void *, void *, int, int);
extern void _shmem_iget32(void *, void *, int, int, int, int);
extern void _shmem_iput32(void *, void *, int, int, int, int);
extern void _shmem_ixget32(void *, void *, short *, int, int);
extern void _shmem_ixput32(void *, void *, short *, int, int);
 
extern void _shmem_short_p(short *addr, short value, int pe);
extern void _shmem_int_p(int *addr, int value, int pe);
extern void _shmem_long_p(long *addr, long value, int pe);
extern void _shmem_Real_p(Real *addr, Real value, int pe);
extern void _shmem_double_p(double *addr, double value, int pe);

extern short _shmem_short_g(short *addr, int pe);
extern int _shmem_int_g(int *addr, int pe);
extern long _shmem_long_g(long *addr, int pe);
extern Real _shmem_Real_g(Real *addr, int pe);
extern double _shmem_double_g(double *addr, int pe);

extern long   _shmem_swap       (long   *addr, long   value, int pe);
extern short  _shmem_short_swap (short  *addr, short  value, int pe);
extern int    _shmem_int_swap   (int    *addr, int    value, int pe);
extern long   _shmem_long_swap  (long   *addr, long   value, int pe);
extern Real  _shmem_Real_swap (Real  *addr, Real  value, int pe);
extern double _shmem_double_swap(double *addr, double value, int pe);

extern short  _shmem_short_finc (short  *addr, int pe);
extern short  _shmem_short_fadd (short  *addr, short value, int pe);

extern void _shmem_set_lock(long *);
extern void _shmem_clear_lock(long *);
extern int _shmem_test_lock(long *);
 
extern void _shmem_barrier(int, int, int, long *);
extern void _shmem_wait(long *, long);
extern void _shmem_wait_until(long *, int, long);
extern long _shmem_swap(long *, long, int);
 
extern void _shmem_set_cache_inv(void);
extern void _shmem_set_cache_line_inv(long *);
extern void _shmem_clear_cache_inv(void);
extern void _shmem_udcflush(void);
extern void _shmem_udcflush_line(long *);
 
extern int _shmem_n_pes(void);
extern int _shmem_my_pe(void);
extern void _shmem_stack(void*);
extern void _shmem_fence(void);
extern void _shmem_quiet(void);


__END_DECLS

/*
 *	EXTERNAL DATA DECLARATIONS
 *
 *      External data.
 */
#ifdef _CRAY1
extern int _shmem_pe_tc_offset[];
#pragma _CRI common _shmem_pe_tc_offset
#endif


/*
 *	UNDEF SECTION
 *
 *	This section has an #undef for every SHMEM routine.
 *	This allows the source file which includes <mpp/shmem.h> to
 *	#define a SHMEM function name to some other value, thus avoiding
 *	the prototype definition above.  This, in turn, allows the name
 *	to appear in a "#pragma _CRI duplicate" directive.   By convention,
 *	the libsma source file which defines foo and _foo defines foo and
 *	addes a "#pragma _CRI duplicate foo as _foo" directive.  Therefore,
 *	only the _foo name (with and underscore) need be #undef'ed.
 *
 */
#undef _shmem_get
#undef _F_shmem_get
#undef _shmem_put
#undef _F_shmem_put
#undef _shmem_iget
#undef _shmem_iput
#undef _shmem_ixget
#undef _shmem_ixput
 
#undef _shmem_get32
#undef _shmem_put32
#undef _shmem_iget32
#undef _shmem_iput32
#undef _shmem_ixget32
#undef _shmem_ixput32
 
#undef _shmem_short_p
#undef _shmem_int_p
#undef _shmem_long_p
#undef _shmem_Real_p
#undef _shmem_double_p

#undef _shmem_short_g
#undef _shmem_int_g
#undef _shmem_long_g
#undef _shmem_Real_g
#undef _shmem_double_g

#undef _shmem_set_lock
#undef _shmem_clear_lock
#undef _shmem_test_lock
 
#undef _shmem_barrier
#undef _shmem_wait
#undef _shmem_wait_until
#undef _shmem_swap
 
#undef _shmem_set_cache_inv
#undef _shmem_set_cache_line_inv
#undef _shmem_clear_cache_inv
#undef _shmem_udcflush
#undef _shmem_udcflush_line
 
#undef _shmem_n_pes
#undef _shmem_my_pe
#undef _shmem_stack
#undef _shmem_fence
#undef _shmem_quiet



/*
 *   OVERVIEW
 *
 *	Some functions also have macro equivalents for better performance.
 *
 *	There are two ways to disable the use of SHMEM C macros:
 *
 *		1) Add an "#undef name" CPP directive in your source file
 *		   after the inclusion of the <mpp/shmem.h> header file.
 *		2) Set the _SHMEM_MACRO_OPT macro to 0 prior to inclusion
 *	           of the <mpp/shmem.h> header file.
 *
 *	The default setting for _SHMEM_MACRO_OPT is 2, which results in the
 *	use of some macros and/or inline functions to handle various
 *	SHMEM requests.  Currently defined settings for _SHMEM_MACRO_OPT are:
 *
 *		0 	No macros are used to activate inline SHMEM requests.
 *		1 	Only macros which do not depend on CRI C compiler
 *			features (intrinsics and built-in macros) are defined.
 *		2	Moderately aggressive use of macros.  (default)
 *		3	Very aggressive use of macros.
 *
 *
 *   CONVENTIONS
 *
 *	Any inlined macro has the same name as the corresponding true 
 *	function.  It may be #undef'ed by a user who wishes to use 
 *	the true function.  Every SHMEM function has an equivalent function
 *	name preceded by an underscore.  This name should be used by 
 *	CRI libraries to avoid collisions with the user namespace.
 *	Therefore, an inlined macro name with the leading underscore must 
 *	also exist for internal library use.  It is desirable for all users
 *	to be able to access SHMEM functions the same way whether or not 
 *	the function is inlined via macro.
 * 
 *	Any inlined SHMEM function has a version of the name beginning with 
 *	"_I_".  This name is always available.  It is used by library code 
 *	which implements the true (non-inlined) function.
 */

#ifndef _SHMEM_MACRO_OPT
#if	defined(_RELEASE) || defined(_CRAYC)	/* if CRI SCC or C++ compiler */
#define _SHMEM_MACRO_OPT 2	/* default CRI macro optimization level */
#else
#define _SHMEM_MACRO_OPT 1	/* default non-CRI macro optimization level */
#endif
#endif	/* ! _SHMEM_MACRO_OPT */

#if	!defined(_LINT_) 


/*****************************************************************************
 *
 *	Macros at level 1 do not depend on CRI C compiler features.
 */
#if	_SHMEM_MACRO_OPT >= 1

/*	
 *	shmem_fence
 */
#if	defined(_CRAYT3D) 

#define shmem_fence             _I_shmem_fence
#define _shmem_fence            _I_shmem_fence
#define _I_shmem_fence()        ((void) 0)

#elif  defined(_CRAYT3E) 
 
#define shmem_fence             _I_shmem_fence
#define _shmem_fence            _I_shmem_fence
#define _I_shmem_fence()        (_shmem_quiet())
 
#elif  defined(_CRAY1)
 
#define shmem_fence             _I_shmem_fence
#define _shmem_fence            _I_shmem_fence
#define _I_shmem_fence()        (_shmem_quiet())
 
#endif

/*	
 *	shmem_stack
 */
#if	defined(_CRAYT3E) 
#define	shmem_stack			_I_shmem_stack
#define	_shmem_stack			_I_shmem_stack
#define	_I_shmem_stack(x)		((void) 0)
#endif
 
/*
 *	The following cache management functions are no-ops on T3E and 
 *	CRAY C90 systems. 
 */
#if	defined(_CRAYT3E) || \
	(defined(_CRAY1) && _MAXVL == 128 && !defined(_ADDR64))

#define	shmem_udcflush			_I_shmem_udcflush
#define	_shmem_udcflush			_I_shmem_udcflush
#define	_I_shmem_udcflush()		((void) 0)

#define	shmem_udcflush_line		_I_shmem_udcflush_line
#define	_shmem_udcflush_line		_I_shmem_udcflush_line
#define	_I_shmem_udcflush_line(x)	((void) 0)

#define	shmem_clear_cache_inv		_I_shmem_clear_cache_inv
#define	_shmem_clear_cache_inv		_I_shmem_clear_cache_inv
#define	_I_shmem_clear_cache_inv()	((void) 0)

#define	shmem_set_cache_inv		_I_shmem_set_cache_inv
#define	_shmem_set_cache_inv		_I_shmem_set_cache_inv
#define	_I_shmem_set_cache_inv()	((void) 0)

#define	shmem_set_cache_line_inv	_I_shmem_set_cache_line_inv
#define	_shmem_set_cache_line_inv	_I_shmem_set_cache_line_inv
#define	_I_shmem_set_cache_line_inv(x)	((void) 0)

#endif	/* T3E or C90 */

/*
 *	Inline some cache management functions on CRAY J90 and CRAY YMP 
 *	systems.
 */
#if	defined(_CRAY1) && _MAXVL == 64

#define	shmem_udcflush			_I_shmem_udcflush
#define	_shmem_udcflush			_I_shmem_udcflush
static void
_I_shmem_udcflush(void)
{
	_semts(1);
	_semclr(1);
}
#pragma _CRI inline _I_shmem_udcflush

#define	shmem_udcflush_line		_I_shmem_udcflush_line
#define	_shmem_udcflush_line		_I_shmem_udcflush_line
#define	_I_shmem_udcflush_line(x)	_I_shmem_udcflush()

#endif	/* J90 (or YMP) */


#endif	/* _SHMEM_MACRO_OPT >= 1 */

/*****************************************************************************
 *
 *	Macros at level 2 and higher might depend on CRI C compiler features.
 */
#if	_SHMEM_MACRO_OPT >= 2

#if	defined(_CRAYT3D) || defined(_CRAYT3E)

#define shmem_quiet		_I_shmem_quiet
#define _shmem_quiet		_I_shmem_quiet
#define _I_shmem_quiet()	_remote_write_barrier()

#elif	defined(_CRAY1)

#define shmem_quiet		_I_shmem_quiet
#define _shmem_quiet		_I_shmem_quiet
#define _I_shmem_quiet()	_cmr()

#endif


/*
 *	One-word PUT inline functions.  
 *
 *	shmem_short_p
 *	shmem_int_p
 *	shmem_long_p
 *	shmem_Real_p
 *	shmem_double_p
 */ 

#ifdef	_CRAYT3D
/*
 *	One-word PUT functions are translated into _shmem_put() calls on T3D
 *	systems.  
 */
#define shmem_short_p _I_shmem_short_p
#define _shmem_short_p _I_shmem_short_p
static void
_I_shmem_short_p(short *addr, short value, int pe)
{
	_shmem_put32(addr, &value, 1, pe);
}
#pragma _CRI inline _I_shmem_short_p

#define shmem_int_p _I_shmem_int_p
#define _shmem_int_p _I_shmem_int_p
static void
_I_shmem_int_p(int *addr, int value, int pe)
{
	_shmem_put(addr, &value, 1, pe);
}
#pragma _CRI inline _I_shmem_int_p

#define shmem_long_p _I_shmem_long_p
#define _shmem_long_p _I_shmem_long_p
#define _I_shmem_long_p(a, value, pe) \
	(_I_shmem_int_p((int*)(a), (int)(value), (pe)))

#define shmem_Real_p _I_shmem_Real_p
#define _shmem_Real_p _I_shmem_Real_p
static void
_I_shmem_Real_p(Real *addr, Real value, int pe)
{
	_shmem_put32(addr, &value, 1, pe);
}
#pragma _CRI inline _I_shmem_Real_p

#define shmem_double_p _I_shmem_double_p
#define _shmem_double_p _I_shmem_double_p
static void
_I_shmem_double_p(double *addr, double value, int pe)
{
	_shmem_put(addr, &value, 1, pe);
}
#pragma _CRI inline _I_shmem_double_p

#endif	/* _CRAYT3D */


#ifdef	_CRAYT3E
/*
 *	One-word PUT functions are inlined on T3E systems.  
 */

#define shmem_short_p _I_shmem_short_p
#define _shmem_short_p _I_shmem_short_p
static void
_I_shmem_short_p(short *addr, short value, int pe)
{
        volatile long * const Ecmd =
                (volatile long *)(_PUT(_MPC_E_REG_SADE) | _MPC_EOM_32BIT);
	volatile short * const sade =
		(volatile short*)(_MPC_E_REG_BASE + 8*_MPC_E_REG_SADE);

#pragma sade 1
	*sade = value;
        _write_memory_barrier();    /* SADE store must precede PUT */
        *Ecmd =                         /* issue the PUT E register command */
                (_MPC_E_REG_STRIDE1 << _MPC_BS_EDATA_MOBE) |
                (              (pe) << _MPC_BS_DFLTCENTPE) |
                (      (long)(addr)                    );
        _memory_barrier();    /* PUT must precede load of EREG_PENDING */
}
#pragma _CRI inline _I_shmem_short_p

#define shmem_int_p _I_shmem_int_p
#define _shmem_int_p _I_shmem_int_p
static void
_I_shmem_int_p(int *addr, int value, int pe)
{
        volatile long * const Ecmd =
                (volatile long *)(_PUT(_MPC_E_REG_SADE));
	volatile int * const sade =
		((volatile int*)_MPC_E_REG_BASE) + _MPC_E_REG_SADE;

#pragma sade 1
	*sade = value;
        _write_memory_barrier();    /* SADE store must precede PUT */
        *Ecmd =                         /* issue the PUT E register command */
                (_MPC_E_REG_STRIDE1 << _MPC_BS_EDATA_MOBE) |
                (              (pe) << _MPC_BS_DFLTCENTPE) |
                (      (long)(addr)                    );
        _memory_barrier();    /* PUT must precede load of EREG_PENDING */
}
#pragma _CRI inline _I_shmem_int_p
 
#define shmem_long_p _I_shmem_long_p
#define _shmem_long_p _I_shmem_long_p
#define _I_shmem_long_p(a, value, pe) \
	(_I_shmem_int_p((int*)(a), (int)(value), (pe)))

#define shmem_Real_p _I_shmem_Real_p
#define _shmem_Real_p _I_shmem_Real_p
static void
_I_shmem_Real_p(Real *addr, Real value, int pe)
{
        volatile long * const Ecmd =
                (volatile long *)(_PUT(_MPC_E_REG_SADE) | _MPC_EOM_32BIT);
	volatile Real * const sade =
		(volatile Real*)(_MPC_E_REG_BASE + 8*_MPC_E_REG_SADE);

#pragma sade 1
	*sade = value;
        _write_memory_barrier();    /* SADE store must precede PUT */
        *Ecmd =                         /* issue the PUT E register command */
                (_MPC_E_REG_STRIDE1 << _MPC_BS_EDATA_MOBE) |
                (              (pe) << _MPC_BS_DFLTCENTPE) |
                (      (long)(addr)                    );
        _memory_barrier();    /* PUT must precede load of EREG_PENDING */
}
#pragma _CRI inline _I_shmem_Real_p

#define shmem_double_p _I_shmem_double_p
#define _shmem_double_p _I_shmem_double_p
static void
_I_shmem_double_p(double *addr, double value, int pe)
{
        volatile long * const Ecmd =
                (volatile long *)(_PUT(_MPC_E_REG_SADE));
	volatile double * const sade =
		((volatile double*)_MPC_E_REG_BASE) + _MPC_E_REG_SADE;

#pragma sade 1
	*sade = value;
        _write_memory_barrier();    /* SADE store must precede PUT */
        *Ecmd =                         /* issue the PUT E register command */
                (_MPC_E_REG_STRIDE1 << _MPC_BS_EDATA_MOBE) |
                (              (pe) << _MPC_BS_DFLTCENTPE) |
                (      (long)(addr)                    );
        _memory_barrier();    /* PUT must precede load of EREG_PENDING */
}
#pragma _CRI inline _I_shmem_double_p

#endif	/* _CRAYT3E */


#ifdef	_CRAY1
/*
 *	One-word PUT functions are inlined on PVP systems.
 */
#define shmem_long_p _I_shmem_long_p
#define _shmem_long_p _I_shmem_long_p
static void
_I_shmem_long_p(long *addr, long value, int pe)
{
        int bias = _shmem_pe_tc_offset[pe] - _shmem_pe_tc_offset[_pe];
	*(addr + bias) = value;
}
#pragma _CRI inline _I_shmem_long_p

#define shmem_short_p _I_shmem_short_p
#define _shmem_short_p _I_shmem_short_p
#define _I_shmem_short_p(a, value, pe) \
	(_I_shmem_long_p((long*)(a), (long)(value), (pe)))

#define shmem_int_p _I_shmem_int_p
#define _shmem_int_p _I_shmem_int_p
#define _I_shmem_int_p(a, value, pe) \
	(_I_shmem_long_p((long*)(a), (long)(value), (pe)))

#define shmem_Real_p _I_shmem_Real_p
#define _shmem_Real_p _I_shmem_Real_p
static void
_I_shmem_Real_p(Real *addr, Real value, int pe)
{
        int bias = _shmem_pe_tc_offset[pe] - _shmem_pe_tc_offset[_pe];
	*(addr + bias) = value;
}
#pragma _CRI inline _I_shmem_Real_p

#define shmem_double_p _I_shmem_double_p
#define _shmem_double_p _I_shmem_double_p
#define _I_shmem_double_p(a, value, pe) \
	(_I_shmem_Real_p((Real*)(a), (Real)(value), (pe)))

#endif	/* _CRAY1 */


/*
 *	One-word GET inline functions.  
 *
 *	shmem_short_p
 *	shmem_int_g
 *	shmem_long_g
 *	shmem_Real_g
 *	shmem_double_g
 */ 

#ifdef	_CRAYT3D
/*
 *	One-word GET functions are translated into _shmem_get() calls on T3D
 *	systems.  
 */
#define shmem_short_g _I_shmem_short_g
#define _shmem_short_g _I_shmem_short_g
static short
_I_shmem_short_g(short *addr, int pe)
{
	short target;
	_shmem_get32(&target, addr, 1, pe);
	return (target);
}
#pragma _CRI inline _I_shmem_short_g

#define shmem_int_g _I_shmem_int_g
#define _shmem_int_g _I_shmem_int_g
static int
_I_shmem_int_g(int *addr, int pe)
{
	int target;
	_shmem_get(&target, addr, 1, pe);
	return (target);
}
#pragma _CRI inline _I_shmem_int_g

#define shmem_long_g _I_shmem_long_g
#define _shmem_long_g _I_shmem_long_g
#define _I_shmem_long_g(a, pe) ((long)_I_shmem_int_g((int*)(a), (pe)))

#define shmem_Real_g _I_shmem_Real_g
#define _shmem_Real_g _I_shmem_Real_g
static Real
_I_shmem_Real_g(Real *addr, int pe)
{
	Real target;
	_shmem_get32(&target, addr, 1, pe);
	return (target);
}
#pragma _CRI inline _I_shmem_Real_g

#define shmem_double_g _I_shmem_double_g
#define _shmem_double_g _I_shmem_double_g
static double
_I_shmem_double_g(double *addr, int pe)
{
	double target;
	_shmem_get(&target, addr, 1, pe);
	return (target);
}
#pragma _CRI inline _I_shmem_double_g

#endif	/* _CRAYT3D */


#ifdef	_CRAYT3E
/*
 *	One-word GET functions are inlined on T3E systems.  
 */

#define shmem_short_g _I_shmem_short_g
#define _shmem_short_g _I_shmem_short_g
static short
_I_shmem_short_g(short *addr, int pe)
{
        volatile long * const Ecmd =
                (volatile long *)(_GET(_MPC_E_REG_SADE) | _MPC_EOM_32BIT);
	volatile short * const sade =
		(volatile short*)(_MPC_E_REG_BASE + 8*_MPC_E_REG_SADE);

#pragma sade 1
        _write_memory_barrier();    	/* prior stores must complete */
        *Ecmd =                         /* issue the GET E register command */
                (_MPC_E_REG_STRIDE1 << _MPC_BS_EDATA_MOBE) |
                (              (pe) << _MPC_BS_DFLTCENTPE) |
                (      (long)(addr)                    );
        _memory_barrier();    /* GET must precede load of SADE */
	return (*sade);
}
#pragma _CRI inline _I_shmem_short_g

#define shmem_int_g _I_shmem_int_g
#define _shmem_int_g _I_shmem_int_g
static int
_I_shmem_int_g(int *addr, int pe)
{
        volatile long * const Ecmd =
                (volatile long *)(_GET(_MPC_E_REG_SADE));
	volatile int * const sade =
		((volatile int*)_MPC_E_REG_BASE) + _MPC_E_REG_SADE;

#pragma sade 1
        _write_memory_barrier();    	/* prior stores must complete */
        *Ecmd =                         /* issue the GET E register command */
                (_MPC_E_REG_STRIDE1 << _MPC_BS_EDATA_MOBE) |
                (              (pe) << _MPC_BS_DFLTCENTPE) |
                (      (long)(addr)                    );
        _memory_barrier();    /* GET must precede load of SADE */
	return (*sade);
}
#pragma _CRI inline _I_shmem_int_g
 
#define shmem_long_g _I_shmem_long_g
#define _shmem_long_g _I_shmem_long_g
#define _I_shmem_long_g(a, pe) ((long)_I_shmem_int_g((int*)(a), (pe)))


#define shmem_Real_g _I_shmem_Real_g
#define _shmem_Real_g _I_shmem_Real_g
static Real
_I_shmem_Real_g(Real *addr, int pe)
{
        volatile long * const Ecmd =
                (volatile long *)(_GET(_MPC_E_REG_SADE) | _MPC_EOM_32BIT);
	volatile Real * const sade =
		(volatile Real*)(_MPC_E_REG_BASE + 8*_MPC_E_REG_SADE);

#pragma sade 1
        _write_memory_barrier();    	/* prior stores must complete */
        *Ecmd =                         /* issue the GET E register command */
                (_MPC_E_REG_STRIDE1 << _MPC_BS_EDATA_MOBE) |
                (              (pe) << _MPC_BS_DFLTCENTPE) |
                (      (long)(addr)                    );
        _memory_barrier();    /* GET must precede load of SADE */
	return (*sade);
}
#pragma _CRI inline _I_shmem_Real_g

#define shmem_double_g _I_shmem_double_g
#define _shmem_double_g _I_shmem_double_g
static double
_I_shmem_double_g(double *addr, int pe)
{
        volatile long * const Ecmd =
                (volatile long *)(_GET(_MPC_E_REG_SADE));
	volatile double * const sade =
		((volatile double*)_MPC_E_REG_BASE) + _MPC_E_REG_SADE;

#pragma sade 1
        _write_memory_barrier();    	/* prior stores must complete */
        *Ecmd =				/* issue the GET E register command */
                (_MPC_E_REG_STRIDE1 << _MPC_BS_EDATA_MOBE) |
                (              (pe) << _MPC_BS_DFLTCENTPE) |
                (      (long)(addr)                    );
        _memory_barrier();    /* GET must precede load of SADE */
	return (*sade);
}
#pragma _CRI inline _I_shmem_double_g

#endif	/* _CRAYT3E */


#ifdef	_CRAY1

#define shmem_long_g _I_shmem_long_g
#define _shmem_long_g _I_shmem_long_g
static long
_I_shmem_long_g(long * restrict addr, int pe)
{
        int bias = _shmem_pe_tc_offset[pe] - _shmem_pe_tc_offset[_pe];
	return *(addr + bias);
}
#pragma _CRI inline _I_shmem_long_g

#define shmem_short_g _I_shmem_short_g
#define _shmem_short_g _I_shmem_short_g
#define _I_shmem_short_g(a, pe) ((short)_I_shmem_long_g((long*)(a), (pe)))

#define shmem_int_g _I_shmem_int_g
#define _shmem_int_g _I_shmem_int_g
#define _I_shmem_int_g(a, pe) ((int)_I_shmem_long_g((long*)(a), (pe)))

#define shmem_Real_g _I_shmem_Real_g
#define _shmem_Real_g _I_shmem_Real_g
static Real
_I_shmem_Real_g(Real * restrict addr, int pe)
{
        int bias = _shmem_pe_tc_offset[pe] - _shmem_pe_tc_offset[_pe];
	return *(addr + bias);
}
#pragma _CRI inline _I_shmem_Real_g

#define shmem_double_g _I_shmem_double_g
#define _shmem_double_g _I_shmem_double_g
#define _I_shmem_double_g(a, pe) ((double)_I_shmem_Real_g((Real*)(a), (pe)))

#endif	/* _CRAY1 */


/*
 *	shmem_put and shmem_get are inlined or partially inlined on some
 *	architectures.
 */

#ifdef	_CRAYT3E

#define shmem_put _I_shmem_put
#define _shmem_put _I_shmem_put
static void
_I_shmem_put(
void * itarg,
void * isrc,
int len,
int remotepe)
{
	if (len == 1)
		_I_shmem_long_p(itarg, *(long*)isrc, remotepe);
	else
		_F_shmem_put(itarg, isrc, len, remotepe);
}
#pragma _CRI inline _I_shmem_put

#endif	/* _CRAYT3E */


#ifdef	_CRAY1

#ifndef _SHMEM_PUTGET_INLINE_THRESHOLD
#ifdef _ADDR64
#define _SHMEM_PUTGET_INLINE_THRESHOLD 512	/* TS */
#elif _MAXVL == 128
#define _SHMEM_PUTGET_INLINE_THRESHOLD 512	/* C90 */
#else
#define _SHMEM_PUTGET_INLINE_THRESHOLD 128	/* J90, YMP */
#endif
#endif  /* ! _SHMEM_PUTGET_INLINE_THRESHOLD */

#define shmem_put _I_shmem_put
#define _shmem_put _I_shmem_put
static void
_I_shmem_put(
void * restrict itarg,
void * restrict isrc,
int len,
int remotepe)
{
	if (len > _SHMEM_PUTGET_INLINE_THRESHOLD)
		_F_shmem_put(itarg, isrc, len, remotepe);
	else {
		int bias, i;
		bias = _shmem_pe_tc_offset[remotepe] - _shmem_pe_tc_offset[_pe];
#pragma _CRI ivdep
		for (i=0; i<len;  i++)
		    ((int * )itarg)[i+bias] = ((int * )isrc)[i];
	}
}
#pragma _CRI inline _I_shmem_put

#define shmem_get _I_shmem_get
#define _shmem_get _I_shmem_get
static void
shmem_get(
void * restrict itarg,
void * restrict isrc,
int len,
int remotepe)
{
	if (len > _SHMEM_PUTGET_INLINE_THRESHOLD)
		_F_shmem_get(itarg,isrc,len,remotepe);
	else {
		int bias, i;
		bias = _shmem_pe_tc_offset[remotepe] - _shmem_pe_tc_offset[_pe];
#pragma _CRI ivdep
		for (i=0; i<len;  i++)
		    ((int * )itarg)[i] = ((int * )isrc)[i+bias];
	}
}
#pragma _CRI inline _I_shmem_get

#endif	/* _CRAY1 */


/*
 *	shmem_wait is inlined on T3E.
 */

#ifdef	_CRAYT3E
#define shmem_wait _I_shmem_wait
#define _shmem_wait _I_shmem_wait
static void
_I_shmem_wait(long *addr, long value)
{
	do {} while (*(volatile long*)addr == value);
}
#pragma _CRI inline _I_shmem_wait

#endif	/* _CRAYT3E */


#endif	/* _SHMEM_MACRO_OPT >= 2 */

/*****************************************************************************
 *
 *	Macros at level 3 are not defined by default.
 */
#if	_SHMEM_MACRO_OPT >= 3

#endif	/* _SHMEM_MACRO_OPT >= 3 */

#endif	/* !_LINT_ */

#endif /* !_MPP_SHMEM_H */
