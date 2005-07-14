/********************** wilson_invert_lean.c *********************************/
/* MIMD version 7 */

/*  Generic wrapper for quark inverter
 *
 *  This wrapper works for WILSON or CLOVER fermions.
 *
 *  This lean version recreates the source before each restart
 *
 *
 */

/* Modifications

   4/30/00 Set min = 0 on restart CD
   4/26/98 Wrapper made generic after introducing structures for
           inversion parameters   CD
   CD 14 Nov 1997
*/

#include "generic_wilson_includes.h"

int wilson_invert_lean( /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (where source is to be created)*/
    field_offset dest,  /* type wilson_vector (answer and initial guess) */
    void (*source_func)(field_offset src, 
			wilson_quark_source *wqs),  /* source function */
    wilson_quark_source *wqs, /* source parameters */
    int (*invert_func)(field_offset src, field_offset dest,
			quark_invert_control *qic,void *dcp),
    quark_invert_control *qic, /* inverter control */
    void *dmp              /* Passthrough Dirac matrix parameters */
    )
{
  int tot_iters, irestart, Minsav;

  Minsav = qic->min;

  /* Inversion with restart (restart is appropriate for CG and BiCG) */
  /* Source in src is overwritten by inverter so we recreate each time */
  for(tot_iters = 0,irestart = 0; irestart < qic->nrestart; irestart++)
    {
      /* Make the source */
      source_func(src,wqs);
      
      /* Do the inversion */
      tot_iters += invert_func(src,dest,qic,dmp);

      /* Check for convergence */
      if(qic->size_r < qic->resid)break;

      /* Restart trial solution is not zero any more */
      qic->start_flag = 1;

      /* No minimum when restarting */
      qic->min = 0;
    }
  if(this_node==0)
    {
      if(qic->size_r > qic->resid)
	printf(" NOT converged size_r= %.2g iters= %d\n",
	       qic->size_r, tot_iters);
      else
	printf(" OK converged size_r= %.2g iters= %d\n",
	       qic->size_r, tot_iters);
    }
  
  /* Restore minimum */
  qic->min = Minsav;
  return tot_iters;
}

/* wilson_invert_lean.c */
