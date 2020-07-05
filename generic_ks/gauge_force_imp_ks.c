/*********************** gauge_force_imp_ks.c  -- ****************************/
/* MIMD version 7 */
/* This version was introduced so KS code could support
   resident GPU operation with phases always in on the CPU.  */

#include "generic_ks_includes.h"	/* definitions files and prototypes */

void imp_gauge_force_ks( Real eps, field_offset mom_off ){
#ifdef USE_GF_GPU
  imp_gauge_force_gpu(eps, mom_off);
#elif USE_GF_QPHIX
  rephase(OFF);
  imp_gauge_force_qphix(eps, mom_off);
  rephase(ON);
#else
  rephase(OFF);
  imp_gauge_force_cpu(eps, mom_off);
  rephase(ON);
#endif
}

/* End of gauge_force_imp_ks.c */
