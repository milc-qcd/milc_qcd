/******* d_congrad5_fn_qdp.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC QDP inverter */
/* It invokes an inverter with the appropriate precision */

#include "generic_ks_includes.h"
#include "../include/generic_qdp.h"
#include "../include/generic_ks_qdp.h"

/* New API for site arguments */

int 
ks_congrad_site( field_offset milc_src, field_offset milc_sol, 
		 quark_invert_control *qic, Real mass,
		 ferm_links_t *fn)
{
  int iterations_used;

  if(qic->prec == 1)
    iterations_used = 
      ks_congrad_milc2qdp_F( milc_src, milc_sol, qic, mass, fn);
  else
    iterations_used = 
      ks_congrad_milc2qdp_D( milc_src, milc_sol, qic, mass, fn);
  
  total_iters += iterations_used;
  return iterations_used;
}

/* New API for field arguments */

int 
ks_congrad_field( su3_vector *milc_src, su3_vector *milc_sol, 
		  quark_invert_control *qic, Real mass,
		  ferm_links_t *fn )
{
  int iterations_used;

  if(qic->prec == 1)
    iterations_used = 
      ks_congrad_milcfield2qdp_F( milc_src, milc_sol, qic, mass, fn);
  else
    iterations_used = 
      ks_congrad_milcfield2qdp_D( milc_src, milc_sol, qic, mass, fn);
  
  total_iters += iterations_used;
  return iterations_used;
}


/* Traditional MILC API for site arguments and no relative residual test */

int
ks_congrad(field_offset f_src, field_offset f_dest, Real mass,
	   int niter, int nrestart, Real rsqmin, int prec, int milc_parity, 
	   Real *final_rsq_ptr, ferm_links_t *fn )
{
  int iteration;
  quark_invert_control qic;

  /* Pack structure */
  qic.prec      = prec;
  qic.parity    = milc_parity;
  qic.max       = niter;
  qic.nrestart  = nrestart;
  qic.resid     = rsqmin;
  qic.relresid  = 0;     /* Suppresses this test */

  if(prec == 1)
    iteration = 
      ks_congrad_milc2qdp_F(f_src, f_dest, &qic, mass, fn );
  else
    iteration = 
      ks_congrad_milc2qdp_D(f_src, f_dest, &qic, mass, fn );
  
  *final_rsq_ptr = qic.final_rsq;
  total_iters += iteration;
  return iteration;
}

