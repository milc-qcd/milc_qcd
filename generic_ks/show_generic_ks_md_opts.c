/********************** show_generic_ks_md_opts.c *************************************/
/* MIMD Version 7 */

/* List options selected in the compilation */
/* These options pertain to the fermion force */

#include "generic_ks_includes.h"

void show_generic_ks_md_opts( void ){

#ifdef KS_MULTIFF
  node0_printf("KS_MULTIFF=%s\n",ks_multiff_opt_chr());

#ifdef VECLENGTH
  node0_printf("VECLENGTH=%d\n",VECLENGTH);
#endif
#endif

}
