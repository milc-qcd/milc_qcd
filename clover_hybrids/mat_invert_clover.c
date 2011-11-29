/************************ mat_invert_clover.c *******************/
/* MIMD version 7 */
/**   Wrapper for the calculation of clover
 **   action quark propagators
 **   7/11/98 Updated for new v5 inverters CD
 **/

#include "cl_hyb_includes.h"

int mat_invert( field_offset src, field_offset dest ){
  register int i;
  register site *s;
  int iters;

  if( src != F_OFFSET(CHI) ){
    FORALLSITES(i,s) s->CHI = *(wilson_vector *)F_PT(s,src);
  }


  /* Load inversion control structure */
  qic.prec = PRECISION;
  qic.min = 0;
  qic.max = niter;
  qic.nrestart = 5;
  qic.parity = EVENANDODD;
  qic.start_flag = 0;
  qic.nsrc = 1;
  qic.resid = sqrt(rsqprop);
  qic.relresid = 0;

  /* Load Dirac matrix parameters */
  dcp.Kappa = kappa;
  dcp.Clov_c = clov_c;
  dcp.U0 = u0;
  
#ifdef BI
  iters = 
    wilson_invert_site(F_OFFSET(CHI),dest,
		       bicgilu_cl_site,&qic,(void *)&dcp);
#else
  
  iters = 
    wilson_invert_site(F_OFFSET(CHI),dest,
		       cgilu_cl_site,&qic,(void *)&dcp);
#endif
  return(iters);
}
