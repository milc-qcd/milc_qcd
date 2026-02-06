/*********************** reunitarize2.c ***************************/
/* MIMD version 7 */

/* reunitarize the link matrices */
/* This version expects KS phases to be in */

#include "generic_ks_includes.h"

void reunitarize_ks() {

  /* Use QUDA if gauge-force is enabled for GPU, but fallback to CPU
     if Schroedinger functional boundary conditions are enabled */

#ifdef SCHROED_FUN
  node0_printf("%s not supported on GPU, using CPU fallback\n", __func__);
#endif

#if defined(USE_GA_GPU) && defined(HAVE_QUDA) && !defined(SCHROED_FUN)
  reunitarize_gpu();
#else
  rephase(OFF);
  reunitarize_cpu();
  rephase(ON);
#endif

}
