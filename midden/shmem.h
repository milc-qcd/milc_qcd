/* USMID @(#)libsma/shmem.h	11.6	10/06/94 09:23:38 */

/*
 *      (C) COPYRIGHT CRAY RESEARCH, INC.
 *      UNPUBLISHED PROPRIETARY INFORMATION.
 *      ALL RIGHTS RESERVED.
 */

/*
 * This file contains user definitions for the shmem routines.
 */

#ifndef _MPP_SHMEM_H
#define _MPP_SHMEM_H

/**#include "../include/complex.h"**/
#include <sys/cdefs.h>

/************************ Constant definitions  *****************************/

#define _SHMEM_BCAST_SYNC_SIZE  	128
#define _SHMEM_BARRIER_SYNC_SIZE        36
#define _SHMEM_REDUCE_SYNC_SIZE		36
#define _SHMEM_COLLECT_SYNC_SIZE	36
#define _SHMEM_SYNC_VALUE       	(-1)
#define _SHMEM_REDUCE_MIN_WRKDATA_SIZE  8

/************************ Function prototypes *******************************/

__BEGIN_DECLS

extern void shmem_get __((long *, long *, int, int));
extern void shmem_put __((long *, long *, int, int));
extern void shmem_iget __((long *, long *, int, int, int, int));
extern void shmem_iput __((long *, long *, int, int, int, int));
extern void shmem_ixget __((long *, long *, long *, int, int));
extern void shmem_ixput __((long *, long *, long *, int, int));

extern void shmem_get32 __((short *, short *, int, int));
extern void shmem_put32 __((short *, short *, int, int));
extern void shmem_iget32 __((short *, short *, int, int, int, int));
extern void shmem_iput32 __((short *, short *, int, int, int, int));
extern void shmem_ixget32 __((short *, short *, short *, int, int));
extern void shmem_ixput32 __((short *, short *, short *, int, int));

extern void shmem_int_sum_to_all __((int *, int *, int, int, int, int, int *,
	long *));
extern void shmem_int_max_to_all __((int *, int *, int, int, int, int, int *,
	long *));
extern void shmem_int_min_to_all __((int *, int *, int, int, int, int, int *,
	long *));
extern void shmem_int_prod_to_all __((int *, int *, int, int, int, int, int *,
	long *));
extern void shmem_int_and_to_all __((int *, int *, int, int, int, int, int *,
	long *));
extern void shmem_int_or_to_all __((int *, int *, int, int, int, int, int *,
	long *));
extern void shmem_int_xor_to_all __((int *, int *, int, int, int, int, int *,
	long *));

extern void shmem_short_sum_to_all __((short *, short *, int, int, int, int,
	short *, long *));
extern void shmem_short_max_to_all __((short *, short *, int, int, int, int,
	short *, long *));
extern void shmem_short_min_to_all __((short *, short *, int, int, int, int,
	short *, long *));
extern void shmem_short_prod_to_all __((short *, short *, int, int, int, int,
	short *, long *));
extern void shmem_short_and_to_all __((short *, short *, int, int, int, int,
	short *, long *));
extern void shmem_short_or_to_all __((short *, short *, int, int, int, int,
	short *, long *));
extern void shmem_short_xor_to_all __((short *, short *, int, int, int, int,
	short *, long *));

extern void shmem_double_sum_to_all __((double *, double *, int, int, int, int,
	double *, long *));
extern void shmem_double_max_to_all __((double *, double *, int, int, int, int,
	double *, long *));
extern void shmem_double_min_to_all __((double *, double *, int, int, int, int,
	double *, long *));
extern void shmem_double_prod_to_all __((double *, double *, int, int, int, int,
	double *, long *));

extern void shmem_Real_sum_to_all __((Real *, Real *, int, int, int, int,
	Real *, long *));
extern void shmem_Real_max_to_all __((Real *, Real *, int, int, int, int,
	Real *, long *));
extern void shmem_Real_min_to_all __((Real *, Real *, int, int, int, int,
	Real *, long *));
extern void shmem_Real_prod_to_all __((Real *, Real *, int, int, int, int,
	Real *, long *));

/**extern void shmem_complexf_sum_to_all __((complex Real *, complex Real *,
	int, int, int, int, complex Real *, long *));
extern void shmem_complexf_prod_to_all __((complex Real *, complex Real *,
	int, int, int, int, complex Real *, long *));

extern void shmem_complexd_sum_to_all __((complex double *, complex double *,
	int, int, int, int, complex double *, long *));
extern void shmem_complexd_prod_to_all __((complex double *, complex double *,
	int, int, int, int, complex double *, long *));**/

extern void shmem_collect __((long *, long *, int, int, int, int, long *));
extern void shmem_fcollect __((long *, long *, int, int, int, int, long *));

extern void shmem_collect32 __((short *, short *, int, int, int, int, long *));
extern void shmem_fcollect32 __((short *, short *, int, int, int, int, long *));

extern void shmem_broadcast __((long *, long *, int, int, int, int, int,
	long *));
extern void shmem_broadcast32 __((short *, short *, int, int, int, int, int,
	long *));

extern void shmem_barrier __((int, int, int, long *));
extern void shmem_wait __((long *, int));
extern long shmem_swap __((long *, long, int));

extern void shmem_set_cache_inv __((void));
extern void shmem_set_cache_line_inv __((long *));
extern void shmem_clear_cache_inv __((void));
extern void shmem_udcflush __((void));
extern void shmem_udcflush_line __((long *));

__END_DECLS

#endif /* !_MPP_SHMEM_H */
