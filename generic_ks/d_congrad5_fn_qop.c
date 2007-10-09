/******* d_congrad5_fn_qop.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */
/* It invokes an inverter with the appropriate precision */

#include "generic_ks_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"

/* New API for site arguments */

int ks_congrad_site( field_offset milc_src, field_offset milc_sol, 
		     quark_invert_control *qic, Real mass,
		     fn_links_t *fn, ks_action_paths *ap )
{
  int iterations_used;

  if(qic->prec == 1)
    iterations_used = 
      ks_congrad_milc2qop_F( milc_src, milc_sol, qic, mass, fn, ap );
  else
    iterations_used = 
      ks_congrad_milc2qop_D( milc_src, milc_sol, qic, mass, fn, ap );
  
  total_iters += iterations_used;
  return iterations_used;
}

/* New API for field arguments */

int ks_congrad_field( su3_vector *milc_src, su3_vector *milc_sol, 
		      quark_invert_control *qic, Real mass,
		     fn_links_t *fn, ks_action_paths *ap )
{
  int iterations_used;

  if(qic->prec == 1)
    iterations_used = 
      ks_congrad_milcfield2qop_F( milc_src, milc_sol, qic, mass, fn, ap );
  else
    iterations_used = 
      ks_congrad_milcfield2qop_D( milc_src, milc_sol, qic, mass, fn, ap );
  
  total_iters += iterations_used;
  return iterations_used;
}

/* Traditional MILC API for site arguments and no relative residual test */

int ks_congrad( field_offset milc_src, field_offset milc_sol, Real mass,
	        int niter, int nrestart, Real rsqmin, int prec, 
		int milc_parity, Real* final_rsq_ptr,
		fn_links_t *fn, ks_action_paths *ap )
{
  int iterations_used;
  quark_invert_control qic;

  /* Pack structure */
  qic.prec      = prec;
  qic.parity    = milc_parity;
  qic.max       = niter;
  qic.nrestart  = nrestart;
  qic.resid     = rsqmin;
  qic.relresid  = 0;     /* Suppresses this test */

  if(prec == 1)
    iterations_used = 
      ks_congrad_milc2qop_F( milc_src, milc_sol, &qic, mass, fn, ap );
  else
    iterations_used = 
      ks_congrad_milc2qop_D( milc_src, milc_sol, &qic, mass, fn, ap );
  
  *final_rsq_ptr = qic.final_rsq;
  total_iters += iterations_used;
  return iterations_used;
}
