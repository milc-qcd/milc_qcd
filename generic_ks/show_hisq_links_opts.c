/****************** show_hisq_links_opts.c **************************/
/* MIMD Version 7 */

/* List options selected in the compilation */

#include "generic_ks_includes.h"

/* Print compiler option macros */
void
show_hisq_links_opts(void){

#ifdef HISQ_REUNIT_ALLOW_SVD
  node0_printf("HISQ_REUNIT_ALLOW_SVD\n");
#endif

#ifdef HISQ_REUNIT_SVD_ONLY
  node0_printf("HISQ_REUNIT_SVD_ONLY (used together with HISQ_REUNIT_ALLOW_SVD)\n");
#endif
  
#ifdef HISQ_REUNIT_SVD_REL_ERROR
  node0_printf("HISQ_REUNIT_SVD_REL_ERROR = %g\n",HISQ_REUNIT_SVD_REL_ERROR);
#endif

#ifdef HISQ_REUNIT_SVD_ABS_ERROR
  node0_printf("HISQ_REUNIT_SVD_ABS_ERROR = %g\n",HISQ_REUNIT_SVD_ABS_ERROR);
#endif

#ifdef HISQ_SVD_VALUES_INFO
  node0_printf("HISQ_SVD_VALUES_INFO\n");
#endif

#ifdef SU3_UNIT_ANALYTIC_FOLLOW_PREC
  node0_printf("SU3_UNIT_ANALYTIC_FOLLOW_PREC is defined but ignored.\n");
  node0_printf("QOPQDP always uses double precision\n");
#endif

}

