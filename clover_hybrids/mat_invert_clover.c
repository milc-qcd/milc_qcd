/************************ mat_invert_clover.c *******************/
/* MIMD version 6 */
/**   Wrapper for the calculation of clover
 **   action quark propagators
 **   7/11/98 Updated for new v5 inverters CD
 **/

#include "cl_hyb_includes.h"

int mat_invert( field_offset src, field_offset dest ){
  register int i;
  register site *s;
  int iters;

  if( src != F_OFFSET(chi) ){
    FORALLSITES(i,s) s->chi = *(wilson_vector *)F_PT(s,src);
  }

  /* Load inversion control structure */
  qic.min = 0;
  qic.max = niter;
  qic.nrestart = 5;
  qic.resid = rsqprop;
  qic.start_flag = 0;
  qic.wv1 = F_OFFSET(tmp);
  qic.wv2 = F_OFFSET(mp);

  /* Load Dirac matrix parameters */
  dcp.Kappa = kappa;
  dcp.Clov_c = clov_c;
  dcp.U0 = u0;
  dcp.work_f_mn = F_OFFSET(f_mn);
  
#ifdef BI
  /* Load temporaries specific to inverter */
  qic.wv3 = F_OFFSET(tmpb);  /* Called rv in bicg */
  qic.wv4 = F_OFFSET(sss);
  
  iters = 
    wilson_invert(F_OFFSET(chi),dest,F_OFFSET(quark_save),
			      bicgilu_cl,&qic,(void *)&dcp);
#else
  
  iters = 
    wilson_invert(F_OFFSET(chi),dest,F_OFFSET(quark_save),
			      cgilu_cl,&qic,(void *)&dcp);
#endif
  return(iters);
}
