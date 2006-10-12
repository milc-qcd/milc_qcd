/******* ks_multicg.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* Wrappers for multi-mass CG inverter for staggered fermions

   8/12 C. DeTar added macros for selecting the multicg inverter option
*/


#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/loopend.h"
#include <string.h>

/* Set the KS multicg inverter flavor depending on the macro KS_MULTICG */

enum ks_multicg_opt_t { KS_MULTICG_OFFSET, KS_MULTICG_HYBRID, KS_MULTICG_FAKE, 
			KS_MULTICG_REVERSE, KS_MULTICG_REVHYB };
static enum ks_multicg_opt_t ks_multicg_opt = KS_MULTICG_HYBRID;   /* Default */

/**********************************************************************/
/*   Set optimization choice                                          */
/**********************************************************************/
/* returns 1 for error and 0 for success */

int ks_multicg_set_opt(char opt_string[]){
  if(strcmp(opt_string,"OFFSET") == 0)
    ks_multicg_opt = KS_MULTICG_OFFSET;
  else if(strcmp(opt_string,"HYBRID") == 0)
    ks_multicg_opt = KS_MULTICG_HYBRID;
  else if(strcmp(opt_string,"FAKE") == 0)
    ks_multicg_opt = KS_MULTICG_FAKE;
  else if(strcmp(opt_string,"REVERSE") == 0)
    ks_multicg_opt = KS_MULTICG_REVERSE;
  else if(strcmp(opt_string,"REVHYB") == 0)
    ks_multicg_opt = KS_MULTICG_REVHYB;
  else{
    printf("ks_multicg_set_opt: Unrecognized type %s\n",opt_string);
    return 1;
  }
  /*printf("ks_multicg_set_opt: set opt to %d\n",ks_multicg_opt);*/
  return 0;
}

/**********************************************************************/
/*   Wrapper for the multimass inverter with multiple sources         */
/**********************************************************************/
int ks_multicg(	        /* Return value is number of iterations taken */
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

  if(ks_multicg_opt == KS_MULTICG_OFFSET)
    return ks_multicg_offset( src, psim, offsets, num_offsets, 
			      niter, rsqmin, parity, final_rsq_ptr);
  else if(ks_multicg_opt == KS_MULTICG_HYBRID)
    return ks_multicg_hybrid( src, psim, offsets, num_offsets, 
			      niter, rsqmin, parity, final_rsq_ptr);
  else if(ks_multicg_opt == KS_MULTICG_FAKE)
    return ks_multicg_fake( src, psim, offsets, num_offsets, 
			      niter, rsqmin, parity, final_rsq_ptr);
  else if(ks_multicg_opt == KS_MULTICG_REVERSE)
    return ks_multicg_reverse( src, psim, offsets, num_offsets, 
			      niter, rsqmin, parity, final_rsq_ptr);
  else if(ks_multicg_opt == KS_MULTICG_REVHYB)
    return ks_multicg_revhyb( src, psim, offsets, num_offsets, 
			      niter, rsqmin, parity, final_rsq_ptr);
  else
    return 0;
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
    int i,j,iters=0; site *s;
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



