/******* ks_multicg_qop.c - multi-mass CG for SU3/fermions ****/
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
ks_multicg_offset(	/* Return value is number of iterations taken */
   field_offset src,	/* source vector (type su3_vector) */
   su3_vector **psim,	/* solution vectors */
   Real *offsets,	/* the offsets */
   int num_offsets,	/* number of offsets */
   int niter,		/* maximal number of CG interations */
   Real rsqmin,	        /* desired residue squared */
   int prec,            /* internal precision for inversion */
   int parity,		/* parity to be worked on */
   Real *final_rsq_ptr	/* final residue squared */
   )
{
  if(prec == 1)
    return ks_multicg_offset_F(src, psim, offsets, num_offsets, niter, 
			       rsqmin, parity, final_rsq_ptr);
  else
    return ks_multicg_offset_D(src, psim, offsets, num_offsets, niter, 
			       rsqmin, parity, final_rsq_ptr);
}

/* Standard MILC interface for the Asqtad multimass inverter 
   single source, multiple masses.  Uses the prevailing MILC precision */

int
ks_multicg_mass(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors (preallocated) */
    Real *masses,	/* the masses */
    int num_masses,	/* number of masses */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* internal precision for inversion */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
			)
{
  if(prec == 1)
    return ks_multicg_mass_F(src, psim, masses, num_masses, niter, 
			     rsqmin, parity, final_rsq_ptr);
  else
    return ks_multicg_mass_D(src, psim, masses, num_masses, niter, 
			     rsqmin, parity, final_rsq_ptr);
}


