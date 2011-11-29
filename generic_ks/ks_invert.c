/********************** ks_invert.c *********************************/
/* MIMD version 7 */

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
		       Real mass, imp_ferm_links_t *fn),
    quark_invert_control *qic, /* inverter control */
    Real mass,
    imp_ferm_links_t *fn
    )
{
  int tot_iters, irestart, Minsav;
  Minsav = qic->min;

  /* Inversion with restart (restart is appropriate for CG and BiCG) */

  for(tot_iters = 0,irestart = 0; irestart < qic->nrestart; irestart++)
    {
      
      tot_iters += invert_func(src,dest,qic,mass,fn);

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
#ifdef CG_OK
      else
	printf(" OK converged size_r= %.2g iters= %d\n",
	       qic->size_r, tot_iters);
#endif
    }
  
  qic->min = Minsav;
  return tot_iters;
} /* ks_invert.c */

/* This variant builds the source, based on the ksqs specification */

int ks_invert_ksqs( /* Return value is number of iterations taken */
    quark_source *ksqs, /* source parameters */
    int (*source_func_field)(su3_vector *src, 
			      quark_source *ksqs),  /* source function */
    su3_vector *dest,  /* answer and initial guess */
    int (*invert_func)(su3_vector *src, su3_vector *dest,
		       quark_invert_control *qic, Real mass, 
		       imp_ferm_links_t *fn),
    quark_invert_control *qic, /* inverter control */
    Real mass,
    imp_ferm_links_t *fn
    )
{
  int tot_iters;
  su3_vector *src;
  char myname[] = "ks_invert_ksqs";

  /* Make the source */
  /* PAD may be used to avoid cache trashing */
#define PAD 0
  src = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
  if(src == NULL){
    printf("%s(%d): Can't allocate src\n",myname,this_node);
    terminate(1);
  }
  if(source_func_field(src,ksqs)){
    printf("%s(%d): error getting source\n",myname,this_node);
    terminate(1);
  };

  /* Do the inversion */
  tot_iters = invert_func(src,dest,qic,mass,fn);

  report_status(qic);

  free(src);
  return tot_iters;
} /* ks_invert_ksqs */


