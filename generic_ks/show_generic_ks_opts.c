/********************** show_generic_ks_opts.c *************************************/
/* MIMD Version 7 */

/* List options selected in the compilation */
/* These options pertain to link construction and the inverter only */

#include "generic_ks_includes.h"

void show_generic_ks_opts( void ){

#ifdef DBLSTORE_FN
  node0_printf("DBLSTORE_FN\n");
#endif

#ifdef D_FN_GATHER13
  node0_printf("D_FN_GATHER13\n");
#endif

#ifdef FEWSUMS
  node0_printf("FEWSUMS\n");
#endif

#ifdef PREFETCH
  node0_printf("PREFETCH\n");
#endif

#ifdef KS_MULTICG
  node0_printf("KS_MULTICG=%s\n",ks_multicg_opt_chr());
#endif

}
