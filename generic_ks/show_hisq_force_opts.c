/***************** show_hisq_force_opts.c *********************************/
/* MIMD Version 7 */

/* List options selected in the compilation */

#include "generic_ks_includes.h"

/* Print compiler option macros */
void 
show_hisq_force_opts(void){

#ifdef HISQ_FORCE_FILTER
  node0_printf("HISQ_FORCE_FILTER = %g\n",HISQ_FORCE_FILTER);
#endif

#ifdef HISQ_FF_MULTI_WRAPPER
  node0_printf("HISQ_FF_MULTI_WRAPPER is ON\n");
#endif

#ifdef HISQ_FF_DEBUG
  node0_printf("HISQ_FF_DEBUG is ON\n");
#endif

}
