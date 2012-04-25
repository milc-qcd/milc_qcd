/******************** report_invert_status.c *****************************/
/* MIMD version 7 */

#include "generic_includes.h"

/* Report inversion status */
void report_status(quark_invert_control *qic){

  if(this_node != 0)return;
  if((qic->resid > 0 && qic->final_rsq <= qic->resid * qic->resid ) ||
     (qic->relresid > 0 && qic->final_relrsq <= qic->relresid * qic->relresid ))
#ifdef CG_OK
    printf(" OK converged final_rsq= %.2g (cf %.2g) rel = %.2g (cf %.2g) restarts = %d iters= %d\n",
	   qic->final_rsq, qic->resid * qic->resid, 
	   qic->final_relrsq, qic->relresid * qic->relresid, 
	   qic->final_restart, qic->final_iters );
#else
  ;
#endif
  else
    printf(" NOT converged final_rsq= %.2g (cf %.2g) rel = %.2g (cf %.2g) restarts = %d iters= %d \n",
	   qic->final_rsq, qic->resid * qic->resid, 
	   qic->final_relrsq, qic->relresid * qic->relresid, 
	   qic->final_restart, qic->final_iters );
} /* report_status */

