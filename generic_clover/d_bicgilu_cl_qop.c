/******* d_bicgilu_cl_qop.c - BiCGstab-ILU interface for QOP ****/
/* MIMD version 7 */

#include "generic_clover_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_clover_qop.h"

/* Standard MILC interface for the clover inverter */

int bicgilu_cl_field(    /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    )
{
  int iterations_used = 0;

  /* Set initial guess at solution */
  if(qic->start_flag == START_ZERO_GUESS)
    clear_wv_field(dest);

  if(qic->prec == 1)
    iterations_used = 
      bicgilu_cl_milc2qop_F( src, dest, qic, dmp );
  else if(qic->prec == 2)
    iterations_used = 
      bicgilu_cl_milc2qop_D( src, dest, qic, dmp );
  else{
    printf("bicgilu_cl_field(%d): Bad precision selection %d\n",
	   this_node,qic->prec);
    terminate(1);
  }
  
  total_iters += iterations_used;
  return iterations_used;
}

int bicgilu_cl_site(    /* Return value is number of iterations taken */
    field_offset src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    )
{
  int i;
  site *s;
  wilson_vector *t_src, *t_dest;
  int iters;

#define PAD 0   /* In case we need to align for cache efficiency */
  t_src  = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  t_dest = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));

  if(t_src == NULL || t_dest == NULL){
    printf("bicgilu_cl_site(%d): Can't allocate src and dest\n",this_node);
    terminate(1);
  }

  /* copy src and dest to temporary */
  FORALLSITES(i,s) {
    t_src[i] = *(wilson_vector *)F_PT(s,src);
    t_dest[i] = *(wilson_vector *)F_PT(s,dest);
  }

  iters = bicgilu_cl_field(t_src, t_dest, qic, dmp);

  /* copy dest back */
  FORALLSITES(i,s) {
    *(wilson_vector *)F_PT(s,dest) = t_dest[i];
  }

  free(t_dest);
  free(t_src);

  return iters;
} /* bicgilu_cl_site */

/* Plain Wilson variant */

int bicgilu_w_field(     /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the WILSON Dirac matrix */
    )
{
  /* Input parameters are for Dirac Wilson! Need remapping. */
  dirac_clover_param dcp;

  map_dwp_to_dcp(&dcp, (dirac_wilson_param *)dmp);
  return bicgilu_cl_field(src, dest, qic, (void *)&dcp);
}


/* Site-based - kept for backward compatibility */
int bicgilu_w_site(      /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    )
{
  /* Input parameters are for Dirac Wilson! Need remapping. */
  dirac_clover_param dcp;

  map_dwp_to_dcp(&dcp, (dirac_wilson_param *)dmp);
  return bicgilu_cl_site(src, dest, qic, (void *)&dcp);
}

