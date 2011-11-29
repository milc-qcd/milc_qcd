OBSOLETE
/******* ks_multicg_offset_qdp.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* Multi-mass CG inverter for staggered fermions */

/* These are wrappers for branching to the appropriate precision */

/* 12/06/06 C. DeTar created */

#include "generic_ks_includes.h"
#include "../include/generic_qdp.h"
#include "../include/generic_ks_qdp.h"

/* Standard MILC interface for the Asqtad multimass inverter 
   single source, multiple masses.  Works from the prevailing MILC precision */

/* Offsets are 4 * mass * mass and must be positive */
int 
ks_multicg_offset_field(	/* Return value is number of iterations taken */
   su3_vector *src,	/* source vector (type su3_vector) */
   su3_vector **psim,	/* solution vectors */
   Real *offsets,	/* the offsets */
   int num_offsets,	/* number of offsets */
   quark_invert_control qic;
   Real *final_rsq_ptr,	/* final residue squared */
   ferm_links_t *fn     /* Storage for fat and Naik links */
   )
{
  if(prec == 1)
    return ks_multicg_offset_field_F(src, psim, offsets, num_offsets, fn);
  else
    return ks_multicg_offset_field_D(src, psim, offsets, num_offsets, fn);
}

#if 0
/* Standard MILC interface for the Asqtad multimass inverter 
   single source, multiple masses.  Uses the prevailing MILC precision */

int
ks_multicg_mass_field(	 /* Return value is number of iterations taken */
    su3_vector *src,	 /* source vector (type su3_vector) */
    su3_vector **psim,	 /* solution vectors (preallocated) */
    ks_param *masses,	 /* KS parameters, including the masses */
    int num_masses,	 /* number of masses */
    quark_invert_control qic;
    ferm_links_t *fn     /* Storage for fat and Naik links */
			)
{
  if(prec == 1)
    return ks_multicg_mass_field_F(src, psim, ksp, num_masses, qic, fn);
  else
    return ks_multicg_mass_field_D(src, psim, ksp, num_masses, qic, fn);
}

#endif
