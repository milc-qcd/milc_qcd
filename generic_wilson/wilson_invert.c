/********************** wilson_invert.c *********************************/
/* MIMD version 6 */

/*  Generic wrapper for quark inverter
 *
 *  This wrapper works for WILSON or CLOVER fermions.
 *
 *  The source routine must be called outside before this
 *  inversion code is called.
 *
 *
 */

/* Modifications

   4/26/98 Wrapper made generic after introducing structures for
           inversion parameters   CD
*/

#include "generic_wilson_includes.h"

int wilson_invert( /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (where source is to be created)*/
    field_offset dest,  /* type wilson_vector (answer and initial guess) */
    field_offset sav,   /* type wilson_vector (for saving source) */
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

  /*** store the source for future use *****/
  FORALLSITES(i,s){
    copy_wvec( (wilson_vector *)F_PT(s,src)  , 
	       (wilson_vector *)F_PT(s,sav)  ) ; 
  }

  Minsav = qic->min;

  /* Inversion with restart (restart is appropriate for CG and BiCG) */

  for(tot_iters = 0,irestart = 0; irestart < qic->nrestart; irestart++)
    {
      
      /* Do the inversion */
      tot_iters += invert_func(src,dest,qic,dmp);

      /* Check for convergence */
      if(qic->size_r < qic->resid)break;

      /* Restore the source for restart (but not for exit) */
      FORALLSITES(i,s){
	copy_wvec( (wilson_vector *)F_PT(s,sav)  , 
		   (wilson_vector *)F_PT(s,src)  ) ; 
      }
      
      /* Restart trial solution is not zero any more */
      qic->start_flag = 1;

      /* No minimum when restarting */
      qic->min = 0;
    }

  /* Report inversion status */
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
} /* wilson_invert.c */
