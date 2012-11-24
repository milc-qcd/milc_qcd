/***************** show_hypisq_force_opts.c *********************************/

/* CD PLACEHOLDER COPIED FROM show_hisq_force_opts.c NEEDS DEVELOPMENT */

/* MIMD Version 7 */

/* List options selected in the compilation */

#include "generic_ks_includes.h"

/* Print compiler option macros */
void 
show_hypisq_force_opts(void){

#ifdef HYPISQ_FORCE_FILTER
  node0_printf("HYPISQ_FORCE_FILTER = %g\n",HYPISQ_FORCE_FILTER);
#endif

#ifdef HYPISQ_FF_MULTI_WRAPPER
  node0_printf("HYPISQ_FF_MULTI_WRAPPER is ON\n");
#endif

#ifdef HYPISQ_FF_DEBUG
  node0_printf("HYPISQ_FF_DEBUG is ON\n");
#endif

}
