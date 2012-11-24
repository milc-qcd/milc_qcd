/****************** show_hypisq_links_opts.c **************************/

/* CD PLACEHOLDER COPIED FROM show_hisq_links_opts.c NEEDS DEVELOPMENT */

/* MIMD Version 7 */

/* List options selected in the compilation */

#include "generic_ks_includes.h"

/* Print compiler option macros */
void
show_hypisq_links_opts(void){

#ifdef HYPISQ_REUNIT_ALLOW_SVD
  node0_printf("HYPISQ_REUNIT_ALLOW_SVD\n");
#endif

#ifdef HYPISQ_REUNIT_SVD_ONLY
  node0_printf("HYPISQ_REUNIT_SVD_ONLY (used together with HYPISQ_REUNIT_ALLOW_SVD)\n");
#endif
  
#ifdef HYPISQ_REUNIT_SVD_REL_ERROR
  node0_printf("HYPISQ_REUNIT_SVD_REL_ERROR = %g\n",HYPISQ_REUNIT_SVD_REL_ERROR);
#endif

#ifdef HYPISQ_REUNIT_SVD_ABS_ERROR
  node0_printf("HYPISQ_REUNIT_SVD_ABS_ERROR = %g\n",HYPISQ_REUNIT_SVD_ABS_ERROR);
#endif

#ifdef HYPISQ_SVD_VALUES_INFO
  node0_printf("HYPISQ_SVD_VALUES_INFO\n");
#endif

#ifdef SU3_UNIT_ANALYTIC_FOLLOW_PREC
  node0_printf("SU3_UNIT_ANALYTIC_FOLLOW_PREC is defined but ignored.\n");
  node0_printf("QOPQDP always uses double precision\n");
#endif

}

