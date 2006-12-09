/******* d_congrad5_fn_qop.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */
/* It invokes an inverter with the appropriate precision */

#include "generic_ks_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"

/* Standard MILC interface for the Asqtad inverter */

int ks_congrad( field_offset milc_src, field_offset milc_sol, Real mass,
	        int niter, int nrestart, Real rsqmin, int prec, 
		int milc_parity, Real* final_rsq_ptr )
{
  int iterations_used;

  if(prec == 1)
    iterations_used = 
      ks_congrad_milc2qop_F( milc_src, milc_sol, mass, niter, nrestart, 
			     rsqmin, milc_parity, final_rsq_ptr);
  else
    iterations_used = 
      ks_congrad_milc2qop_D( milc_src, milc_sol, mass, niter, nrestart, 
			     rsqmin, milc_parity, final_rsq_ptr);
  
  total_iters += iterations_used;
  return iterations_used;
}
