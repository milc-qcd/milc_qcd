/******* ks_multicg.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* Wrappers for multi-mass CG inverter for staggered fermions

   8/12 C. DeTar added macros for selecting the multicg inverter option
*/


#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/loopend.h"

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



