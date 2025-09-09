/******* d_congrad5_fn_qop.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */
/* It invokes an inverter with the appropriate precision */

#include "generic_ks_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"

int
ks_congrad_parity_cpu( su3_vector *src, su3_vector *sol, 
		       quark_invert_control *qic, Real mass,
		       imp_ferm_links_t *fn){

  int iterations_used;

  if(qic->prec == 1)
    iterations_used = 
      ks_congrad_milcfield2qop_F( src, sol, qic, mass, fn );
  else
    iterations_used = 
      ks_congrad_milcfield2qop_D( src, sol, qic, mass, fn );
  
  total_iters += iterations_used;
  return iterations_used;
}

int ks_congrad_block_parity_cpu(int nsrc, su3_vector **t_src, su3_vector **t_dest, 
				quark_invert_control *qic, Real mass,
				imp_ferm_links_t *fn)
{

  /* FAKE version for now */
  int iters = 0;
  for(int i = 0; i < nsrc; i++)
    iters += ks_congrad_parity_cpu(t_src[i], t_dest[i], qic, mass, fn);
  return iters;
}
