/*********************** ploop3_ks.c  -- ****************************/
/* MIMD version 7 */
/* This version was introduced so KS code could support
   resident GPU operation with phases always in on the CPU.  */

#include "generic_ks_includes.h"        /* definitions files and prototypes */

complex ploop3_ks(void) {
    complex p_loop;
#if defined (HAVE_QUDA) && defined(USE_GA_GPU) && !defined(BPCORR)
    p_loop = ploop_gpu();
#else
    p_loop = ploop_cpu();
#endif
    return p_loop;
}

/* ploop3_ks.c */
