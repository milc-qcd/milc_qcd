/*********************** reunitarize2.c ***************************/
/* MIMD version 7 */

/* reunitarize the link matrices */
/* This version expects KS phases to be in */

#include "generic_ks_includes.h"

void reunitarize_ks() {

#ifdef USE_GF_GPU // temporarily disable

  /* Use QUDA if gauge-force is enabled for GPU, but fallback to CPU
     if Schroedinger functional boundary conditions are enabled */
#ifdef SCHROED_FUN
  node0_printf("%s not supported on GPU, using CPU fallback\n", __func__);
  rephase(OFF);
  reunitarize_cpu();
  rephase(ON);
#else
  reunitarize_gpu();
#endif

#else
  rephase(OFF);
  reunitarize_cpu();
  rephase(ON);
#endif

}
