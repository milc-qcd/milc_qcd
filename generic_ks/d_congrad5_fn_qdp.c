/******* d_congrad5_fn_qdp.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC QDP inverter */
/* It invokes an inverter with the appropriate precision */

#include "generic_ks_includes.h"
#include "../include/generic_qdp.h"
#include "../include/generic_ks_qdp.h"

/* Standard MILC interface for the Asqtad inverter */

int
ks_congrad(field_offset f_src, field_offset f_dest, Real mass,
	   int niter, int nrestart, Real rsqmin, int prec, int parity, 
	   Real *final_rsq_ptr)
{
  int iteration;

  if(prec == 1)
    iteration = 
      ks_congrad_milc2qdp_F(f_src, f_dest, mass, niter, nrestart, 
			    rsqmin, parity, final_rsq_ptr);
  else
    iteration = 
      ks_congrad_milc2qdp_D(f_src, f_dest, mass, niter, nrestart, 
			    rsqmin, parity, final_rsq_ptr);
  
  return iteration;
}
