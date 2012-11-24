/******* ks_multicg_offset_qop_P.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */

/* NOTE: This code is actually an include file for ks_multicg_qop_F.c
   and ks_multicg_qop_D.c, so any edits should be consistent with this
   purpose. */

/* Entry points (must be redefined to precision-specific names)

   KS_MULTICG_OFFSET_FIELD

   Externals (must be redefined to precision-specific names)

   KS_CONGRAD_QOP_FIELD2FIELD   

*/


/* 6/2006 C. DeTar created */

/*
 * $Log: ks_multicg_offset_qop_P.c,v $
 * Revision 1.8  2012/11/24 00:02:53  detar
 * Add placeholders for HYPISQ action.  Support HISQ action within ks_imp_dyn.
 *
 * Revision 1.7  2011/12/03 22:30:09  detar
 * Cosmetic: Fix file name in top comment line
 *
 * Revision 1.6  2011/11/29 20:45:57  detar
 * Support new fermion links scheme
 *
 * Revision 1.5  2007/11/16 04:07:15  detar
 * Add parity to QDP "set" utilities
 *
 * Revision 1.4  2007/11/09 16:42:41  detar
 * Pull FN link calculation out of inverters
 *
 * Revision 1.3  2007/10/09 20:10:14  detar
 * Add ferm_links_t and ks_action_paths structures and pass them as params
 *
 * Revision 1.2  2007/05/21 05:06:50  detar
 * Change stopping condition to true residual.
 *
 * Revision 1.1  2006/12/09 13:52:40  detar
 * Add mixed precision capability for KS inverter in QOP and QDP
 *
 * Revision 1.3  2006/11/04 23:41:21  detar
 * Add QOP and QDP support for FN fermion links
 * Create QDP version of fermion_links_fn_multi
 * Add nrestart parameter for ks_congrad
 *
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

#if ( QOP_PrecisionInt == 1 )

#define KS_MULTICG_OFFSET_FIELD         ks_multicg_offset_field_F
#define KS_CONGRAD_QOP_FIELD2FIELD      ks_congrad_qop_F_field2field
#define MYREAL                          float

#else

#define KS_MULTICG_OFFSET_FIELD         ks_multicg_offset_field_D
#define KS_CONGRAD_QOP_FIELD2FIELD      ks_congrad_qop_D_field2field
#define MYREAL                          double

#endif

#include "generic_ks_includes.h"
#include "../include/generic_ks_qop.h"
#include "../include/loopend.h"

//static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/ks_multicg_offset_qop_P.c,v 1.8 2012/11/24 00:02:53 detar Exp $";

/* Standard MILC interface for the Asqtad multimass inverter 
   single source, multiple masses.  Uses the prevailing precision */

/* Offsets are 4 * mass * mass and must be positive */
int KS_MULTICG_OFFSET_FIELD(	      /* Return value is number of iterations taken */
    su3_vector *src,	      /* source vector (type su3_vector) */
    su3_vector **psim,	      /* solution vectors */
    ks_param *ksp,	      /* the offsets */
    int num_offsets,	      /* number of offsets */
    quark_invert_control *qic,/* inversion parameters */
    imp_ferm_links_t *fn      /* Storage for fat and Naik links */
    )
{
  int num_masses = num_offsets;
  int i,j;
  site *s;
  MYREAL *masses;
  int iterations_used;
  MYREAL *masses2[1];
  int nmass[1], nsrc;
  int parity = qic[0].parity;     /* MILC parity */
  su3_vector *milc_srcs[1];
  su3_vector **milc_sols[1];
  char myname[] = "ks_multicg_offset_field";
  
  /* Generate mass table */
  masses = (MYREAL *)malloc(sizeof(MYREAL)*num_masses);
  if(masses == NULL){
    printf("%s(%d): No room for masses\n",myname,this_node);
    terminate(1);
  }
  for(i = 0; i < num_masses; i++){
    if(ksp[i].offset <= 0){
      printf("%s(%d): called with nonpositive offset %e\n",
	     myname, this_node,ksp[i].offset);
      terminate(1);
    }
    masses[i] = sqrt(ksp[i].offset/4.0);
  }
  
  /* Set up general source and solution pointers for one mass, one source */
  nsrc = 1;
  milc_srcs[0] = src;
  
  nmass[0] = num_masses;
  masses2[0] = masses;
  
  /* Require zero initial guess for multicg */
  for(j = 0; j < num_masses; j++){
    FORSOMEPARITY(i,s,parity){
      clearvec(psim[j]+i);
    }
  } END_LOOP
  
  /* Just set pointer for source 1 array of solutions */
  milc_sols[0] =  psim;
  
  iterations_used = KS_CONGRAD_QOP_FIELD2FIELD( qic, masses2, nmass, milc_srcs,
						milc_sols, nsrc, fn );
  
  free(masses);
  return iterations_used;
}

