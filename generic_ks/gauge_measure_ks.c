/*********************** gauge_measure_ks.c  -- ****************************/
/* MIMD version 7 */
/* This version was introduced so KS code could support
   resident GPU operation with phases always in on the CPU.  */

#include "generic_ks_includes.h"        /* definitions files and prototypes */

/* this is needed to query ANISOTROPY, NREPS */
#define GAUGE_ACTION_PART1
/* defines NREPS NLOOP MAX_LENGTH MAX_NUM */
#include <gauge_action.h>
#undef GAUGE_ACTION_PART1

void g_measure_ks( ) {
#if defined (HAVE_QUDA) && defined(USE_GA_GPU) && !defined(ANISOTROPY) && !defined(BPCORR) && NREPS == 1
    g_measure_gpu();
#else
    rephase( OFF );
    g_measure();
    rephase( ON );
#endif
}

/* gauge_measure_ks.c */
