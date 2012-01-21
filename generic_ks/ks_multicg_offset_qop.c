/******* ks_multicg_offset_qop.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* Multi-mass CG inverter for staggered fermions */

/* These are wrappers for branching to the appropriate precision */

/* 12/06/06 C. DeTar created */

#include "generic_ks_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"

/* Standard MILC interface for the Asqtad multimass inverter 
   single source, multiple masses.  Works from the prevailing MILC precision */

/* Offsets are 4 * mass * mass and must be positive */
int 
ks_multicg_offset_field_cpu(	/* Return value is number of iterations taken */
   su3_vector *src,	/* source vector (type su3_vector) */
   su3_vector **psim,	/* solution vectors */
   ks_param *ksp,	/* the offsets */
   int num_offsets,	/* number of offsets */
   quark_invert_control *qic,  /* inversion parameters */
   imp_ferm_links_t *fn       /* Storage for fat and Naik links */
   )
{

  int i, iters;

  if(qic->prec == 1)
    iters = ks_multicg_offset_field_F(src, psim, ksp, num_offsets, qic, fn);
  else
    iters = ks_multicg_offset_field_D(src, psim, ksp, num_offsets, qic, fn);

  /* Copy multimass qic results to the rest of the structures */
  for(i = 1; i < num_offsets; i++){
    qic[i].final_rsq     = qic[0].final_rsq;
    qic[i].final_relrsq  = qic[0].final_relrsq;
    qic[i].size_r        = qic[0].size_r;   
    qic[i].size_relr     = qic[0].size_relr;
    qic[i].final_iters   = qic[0].final_iters;
    qic[i].final_restart = qic[0].final_restart;
    qic[i].converged     = qic[0].converged;  
  }

  return iters;
}

