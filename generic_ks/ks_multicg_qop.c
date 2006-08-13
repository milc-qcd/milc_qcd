/******* ks_multicg_qop.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */

/* 6/2006 C. DeTar Created */

/*
 * $Log: ks_multicg_qop.c,v $
 * Revision 1.1  2006/08/13 15:01:50  detar
 * Realign procedures to accommodate ks_imp_rhmc code
 * Add Level 3 wrappers and MILC dummy Level 3 implementation
 *
 */

#include "generic_ks_includes.h"
#include "../include/loopend.h"
#include <qop.h>

int ks_congrad_qop_site2field(int niter, Real rsqmin, 
			      Real *masses[], int nmass[], 
			      field_offset milc_srcs[], 
			      su3_vector **milc_sols[],
			      int nsrc, Real* final_rsq_ptr, int milc_parity );

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/Attic/ks_multicg_qop.c,v 1.1 2006/08/13 15:01:50 detar Exp $";

/* Set the KS multicg inverter flavor depending on the macro KS_MULTICG */

#define OFFSET  0
#define HYBRID  1
#define FAKE    2
#define REVERSE 3
#define REVHYB  4

ks_multicg_t ks_multicg_init(){
#if (KS_MULTICG == OFFSET)
  return ks_multicg_offset;
#elif (KS_MULTICG == HYBRID)
  return ks_multicg_hybrid;
#elif (KS_MULTICG == FAKE)
  return ks_multicg_fake;
#elif (KS_MULTICG == REVERSE)
  return ks_multicg_reverse;
#elif (KS_MULTICG == REVHYB)
  return ks_multicg_revhyb;
#elif defined(KS_MULTICG)
  node0_print ("ks_multicg_init: unknown or missing KS_MULTICG macro\n");
  return NULL;
#else
  return ks_multicg_offset; /* Default when KS_MULTICG is not defined */
#endif
}

#undef OFFSET
#undef HYBRID
#undef FAKE  
#undef REVERSE
#undef REVHYB

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
  int i;
  Real *masses;
  int status;
  int num_masses = num_offsets;

  masses = (Real *)malloc(sizeof(Real)*num_offsets);
  if(masses == NULL){
    printf("ks_multicg_mass: No room for masses\n");
    terminate(1);
  }
  for(i = 0; i < num_offsets; i++){
    if(offsets[i] < 0){
      printf("ks_multicg_mass(%d): called with negative offset %e\n",
		   this_node,offsets[i]);
      terminate(1);
    }
    masses[i] = sqrt(offsets[i]/4.0);
  }
  status = ks_multicg_mass(src, psim, masses, num_masses, niter, rsqmin,
			   parity, final_rsq_ptr);
  free(masses);
  return status;
}

// mock up multicg by repeated calls to ordinary cg
int ks_multicg_fake(	/* Return value is number of iterations taken */
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
    int i,j,iters; site *s;
#ifdef ONEMASS
    field_offset tmp = F_OFFSET(xxx);
#else
    field_offset tmp = F_OFFSET(xxx1);
#endif

    iters=0;
    for(i=0;i<num_offsets;i++){
       iters += ks_congrad( src, tmp, 0.5*sqrt(offsets[i]), niter, rsqmin, parity, final_rsq_ptr );
       FORALLSITES(j,s){ psim[i][j] = *((su3_vector *)F_PT(s,tmp)); }
    }
    return(iters);
}

// Do a multimass CG followed by calls to individual CG's
// to finish off.
int ks_multicg_hybrid(	/* Return value is number of iterations taken */
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
    int i,j,iters; site *s;
#ifdef ONEMASS
    field_offset tmp = F_OFFSET(xxx);
#else
    field_offset tmp = F_OFFSET(xxx1);
#endif

    ks_multicg_offset( src, psim, offsets, num_offsets, niter, rsqmin, parity, final_rsq_ptr);
    for(i=0;i<num_offsets;i++){
       FORSOMEPARITY(j,s,parity){ *((su3_vector *)F_PT(s,tmp)) = psim[i][j]; } END_LOOP
       iters += ks_congrad( src, tmp, 0.5*sqrt(offsets[i]), niter/5, rsqmin, parity, final_rsq_ptr );
       FORSOMEPARITY(j,s,parity){ psim[i][j] = *((su3_vector *)F_PT(s,tmp)); } END_LOOP
    }
    return(iters);
}
