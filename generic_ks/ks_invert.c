/********************** ks_invert.c *********************************/
/* MIMD version 6 */

/*  Generic wrapper for KS quark inverter
 *
 *  The source routine must be called outside before this
 *  inversion code is called.
 *
 *  This wrapper works for multiple sources as well as a single source
 *
 */

#include "generic_ks_includes.h"

int ks_invert( /* Return value is number of iterations taken */
    field_offset src,   /* type su3_vector (preloaded source) */
    field_offset dest,  /* type su3_vector (answer and initial guess) */
    int (*invert_func)(field_offset src, field_offset dest,
			quark_invert_control *qic,
			void *dmp),
    quark_invert_control *qic, /* inverter control */
    void *dmp                 /* Passthrough Dirac matrix parameters */
    )
{
  int tot_iters, irestart, Minsav;
  int i ; 
  register site *s; 

  Minsav = qic->min;

  /* Inversion with restart (restart is appropriate for CG and BiCG) */

  for(tot_iters = 0,irestart = 0; irestart < qic->nrestart; irestart++)
    {
      
      tot_iters += invert_func(src,dest,qic,dmp);

      if(qic->converged)break;

      qic->start_flag = 1;
      /* No minimum when restarting */
      qic->min = 0;
    }

  if(this_node==0)
    {
      if(qic->converged)
	printf(" NOT converged size_r= %.2g iters= %d\n",
	       qic->size_r, tot_iters);
      else
	printf(" OK converged size_r= %.2g iters= %d\n",
	       qic->size_r, tot_iters);
    }
  
  qic->min = Minsav;
  return tot_iters;
} /* wilson_invert.c */
