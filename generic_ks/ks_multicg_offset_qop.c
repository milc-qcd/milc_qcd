/******* ks_multicg_qop.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */

/* 6/2006 C. DeTar Created */

/*
 * $Log: ks_multicg_offset_qop.c,v $
 * Revision 1.2  2006/08/22 21:16:40  detar
 * Move functionality from ks_multicg_qop.c to ks_multicg_offset_qop.c
 * Remove extraneous globals from ks_multicg_offset.c
 * Remove make rule for ks_multicg_qop.c from Make_template
 *
 * Revision 1.1  2006/08/13 15:01:50  detar
 * Realign procedures to accommodate ks_imp_rhmc code
 * Add Level 3 wrappers and MILC dummy Level 3 implementation
 *
 */

#include "generic_ks_includes.h"
#include "../include/loopend.h"
#include <qop.h>

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/ks_multicg_offset_qop.c,v 1.2 2006/08/22 21:16:40 detar Exp $";

/* Standard MILC interface for the Asqtad multimass inverter 
   single source, multiple masses.  Uses the prevailing precision */

/* Offsets are 4 * mass * mass and must be positive */
int ks_multicg_offset(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
    )
{
  int num_masses = num_offsets;
  int i,j;
  Real *masses;
  int iterations_used;
  Real *masses2[1];
  int nmass[1], nsrc;
  field_offset milc_srcs[1];
  su3_vector **milc_sols[1];

  /* Generate mass table */
  masses = (Real *)malloc(sizeof(Real)*num_masses);
  if(masses == NULL){
    printf("ks_multicg_mass: No room for masses\n");
    terminate(1);
  }
  for(i = 0; i < num_masses; i++){
    if(offsets[i] < 0){
      printf("ks_multicg_mass(%d): called with negative offset %e\n",
		   this_node,offsets[i]);
      terminate(1);
    }
    masses[i] = sqrt(offsets[i]/4.0);
  }

  /* Set up general source and solution pointers for one mass, one source */
  nsrc = 1;
  milc_srcs[0] = src;

  nmass[0] = num_masses;
  masses2[0] = masses;

  /* Require zero initial guess for multicg */
  for(i = 0; i < num_masses; i++)
    for(j = 0; j < sites_on_node; j++)
      clearvec(psim[i]+j);
	
  /* Just set pointer for source 1 array of solutions */
  milc_sols[0] =  psim;

  iterations_used = ks_congrad_qop_site2field( niter, rsqmin, 
					       masses2, nmass, milc_srcs,
					       milc_sols, nsrc, final_rsq_ptr,
					       parity );

  free(masses);
  return iterations_used;
}


/* Standard MILC interface for the Asqtad multimass inverter 
   single source, multiple masses.  Uses the prevailing precision */

int ks_multicg_mass(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors (preallocated) */
    Real *masses,	/* the masses */
    int num_masses,	/* number of masses */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
    )
{

  int iterations_used;
  Real *masses2[1];
  int nmass[1], nsrc;
  field_offset milc_srcs[1];
  su3_vector **milc_sols[1];

  /* Set up general source and solution pointers for one mass, one source */
  nsrc = 1;
  milc_srcs[0] = src;

  nmass[0] = num_masses;
  masses2[0] = masses;

  milc_sols[0] =  psim;

  iterations_used = ks_congrad_qop_site2field( niter, rsqmin, 
					       masses2, nmass, milc_srcs,
					       milc_sols, nsrc, final_rsq_ptr,
					       parity );

  total_iters += iterations_used;
  return( iterations_used );
}


