/*********************** gauge_action_imp_ks.c  -- ****************************/
/* MIMD version 7 */
/* This version was introduced so KS code could support
   resident GPU operation with phases always in on the CPU.  */

#include "generic_ks_includes.h"        /* definitions files and prototypes */

double imp_gauge_action_ks(void) {

#if defined(USE_GA_GPU) && defined(HAVE_GRID) && !defined(ANISOTROPY)
  return imp_gauge_action_grid();
#elif defined(USE_GA_GPU) && defined(HAVE_GRID) && !defined(ANISOTROPY)
  return imp_gauge_action_quda();
#else
  rephase(OFF);
  double g_action = imp_gauge_action();
  rephase(ON);
  return g_action;
#endif
}

/* End of gauge_action_imp_ks.c */
