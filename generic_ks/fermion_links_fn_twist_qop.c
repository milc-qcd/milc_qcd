/************** fermion_links_fn_twist_qop.c **************************/
/* MIMD version 7 */

/* Apply (or undo) a boundary twist to the fat and long links in the
   fn_links_t structure and shift the KS phases and antiperiodic BC.to
   a new origin.
*/

/* The resulting transformation is equivalent to applying it to the
   underlying gauge field and then recomputing the fat and long links
   from the resulting twisted gauge field.

   The twist is specified by four real numbers, bdry_phase[mu]. The
   phase is defined as exp(i Pi bdry_phase[mu]).  Thus bdry_phase[TUP]
   = 1 has the effect of switching from antiperiodic to periodic
   boundary conditions in the time direction.

   The boundary of the lattice is defined relative to the origin
   coordinate vectior r0.  When the twist is applied, we also shift
   the KS phases and antiperiodic boundary condition to correspond to
   the new origin r0.

   The convention adopted here is that the untwisted fields follow
   the standard antiperiodic BC with the KS phases in and defined
   relative to the natural origin (0,0,0,0).  The twisted fields
   obey the twisted BC with the KS phases in but defined relative to
   the shifted origin r0.

   The link phase information, origin shift, and twist state is kept
   in the fermion links structure.

   Entry points:

   create_link_phase_info
   set_boundary_twist_fn
   copy_link_phase_info
   destroy_link_phase_info
   boundary_twist_fn

*/

#include "generic_ks_includes.h"
#include "../include/generic_ks_qop.h"

/*----------------------------------------------------------*/
/* Create, set, copy, destroy the link phase info structure */
/*----------------------------------------------------------*/

link_phase_info_t *
create_link_phase_info(void){
  
  link_phase_info_t *lp;
  char myname[] = "create_link_phase_info";
  int dir;

  lp = (link_phase_info_t *)malloc(sizeof(link_phase_info_t));
  if(lp == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  lp->twist_in = OFF;

  FORALLUPDIR(dir){
    lp->bdry_phase[dir] = 0.;
    lp->r0[dir] = 0;
  }

  return lp;
}

void 
set_boundary_twist_fn(fn_links_qop_t *fn, Real bdry_phase[4], int r0[4]){
  int dir;
  
  FORALLUPDIR(dir){
    fn->phase->bdry_phase[dir] = bdry_phase[dir];
    fn->phase->r0[dir] = r0[dir];
  }
}

void 
copy_link_phase_info(link_phase_info_t *dest, link_phase_info_t *src){
  int dir;

  if(src == NULL || dest == NULL)return;
  
  FORALLUPDIR(dir){
    dest->bdry_phase[dir] = src->bdry_phase[dir];
  }

  dest->twist_in = src->twist_in;
}

void 
destroy_link_phase_info(link_phase_info_t *lp){
  if(lp == NULL)return;
  free(lp);
}

/*--------------------------------------------------------------------*/
/* Shift the time antiperiodic boundary condition and KS phases       */
/*--------------------------------------------------------------------*/

/* STANDARD MILC STAGGERED PHASES */
/*    eta_t = 1
      eta_x = (-1)^t
      eta_y = (-1)^(t+x)
      eta_z = (-1)^(t+x+y)
*/

/* QOP signmask values for MILC staggered phases:
   There is one four-bit mask value for each direction.
   The directions are numbered x = 0, y = 1, z = 2, t = 3.
   The four bits in a given mask value are indexed tzyx.  
   When the bit in a given mask is 1, the indexed coordinate appears 
   in the exponent of (-1) for the eta in the above expressions.
   The mask values are given as octal numbers */


#define MILC_KS_PHASES  { 010, 011, 012, 000 }

static void 
KS_rephase(fn_links_qop_t *fn, int r0[]){
  QOP_F3_FermionLinksAsqtad *al_F = fn->al_F;
  QOP_D3_FermionLinksAsqtad *al_D = fn->al_D;
  QOP_bc_t cp;
  int mask[4] = MILC_KS_PHASES;
  QOP_staggered_sign_t ks_phases = { mask };
  int dir;

  cp.phase = (QOP_Complex *)malloc(sizeof(QOP_Complex)*4);

  /* Conventional antiperiodic bc */
  FORALLUPDIR(dir){
    cp.phase[dir].re = 1.0;
    cp.phase[dir].im = 0.0;
  }
  cp.phase[TUP].re = -1.0;

  if(al_F != NULL){
    QOP_F3_asqtad_rephase_L(al_F, r0, &cp, &ks_phases);
  }
  if(al_D != NULL){
    QOP_D3_asqtad_rephase_L(al_D, r0, &cp, &ks_phases);
  }
      
  free(cp.phase);
}

/*--------------------------------------------------------------------*/
/* Do the momentum twist (AKA boundary twist)                         */
/*--------------------------------------------------------------------*/

static void
apply_twist(fn_links_qop_t *fn, int r0[], complex cphase[4]){
  
  QOP_F3_FermionLinksAsqtad *al_F = fn->al_F;
  QOP_D3_FermionLinksAsqtad *al_D = fn->al_D;
  QOP_bc_t cp;
  int dir;

  cp.phase = (QOP_Complex *)malloc(sizeof(QOP_Complex)*4);
  FORALLUPDIR(dir){
    cp.phase[dir].re = cphase[dir].real;
    cp.phase[dir].im = cphase[dir].imag;
  }
  
  /* QOP links */
  if(al_F != NULL){
    QOP_F3_asqtad_rephase_L(al_F, r0, &cp, 0);
  }
  if(al_D != NULL){
    QOP_D3_asqtad_rephase_L(al_D, r0, &cp, 0);
  }

  free(cp.phase);
}

void 
boundary_twist_fn(fn_links_qop_t *fn, int flag) {
  char myname[] = "boundary_twist_fn";

  link_phase_info_t *lp = fn->phase;
  int *status_now = &lp->twist_in;
  Real *bdry_phase = lp->bdry_phase;
  int *r0 = lp->r0;
  int dir;
  complex cphase[4];
  int no_twist;
  int zero[4] = {0, 0, 0, 0};

#ifndef FN
  node0_printf("Boundary twists are supported only for FN actions\n");
  terminate(1);
#endif

  /* Check that we have links defined */
  if(fn==NULL){
    node0_printf("%s: Attempting to twist undefined FN links\n",myname);
    terminate(1);
  }

  /* Announce action */

  if(flag == ON){
    node0_printf("Turning ON boundary phases %g %g %g %g to FN links r0 %d %d %d %d\n",
		 bdry_phase[XUP], bdry_phase[YUP], bdry_phase[ZUP], bdry_phase[TUP],
		 r0[XUP], r0[YUP], r0[ZUP], r0[TUP]);
  } else {
    node0_printf("Turning OFF boundary phases %g %g %g %g to FN links r0 %d %d %d %d\n",
		 bdry_phase[XUP], bdry_phase[YUP], bdry_phase[ZUP], bdry_phase[TUP],
		 r0[XUP], r0[YUP], r0[ZUP], r0[TUP]);
  }

  /* Check to make sure we are going in expected direction */
  if( !( flag==ON && *status_now==OFF )  && !( flag==OFF && *status_now==ON ) ){
    node0_printf("%s: DUMMY: you fouled up the phases\n",myname);
    terminate(1);
  }

  /* Remove KS phases and APBC if they are changing */
  /* (No %2 on r0[TUP] because it shifts the APBC) */
  if(r0[XUP]%2 != 0 || r0[YUP]%2 != 0 || r0[ZUP]%2 != 0 || r0[TUP] != 0){

    if(*status_now == OFF){
      /* Switch off the conventional KS phases and antiperiodic boundary condition */
      KS_rephase(fn, zero);
    } else {
      /* Switch off the unconventional KS phases and antiperiodic boundary condition */
      KS_rephase(fn, r0);
    }
  }
  
  /* No twist needed if boundary phases are all zero */

  no_twist = 1;
  FORALLUPDIR(dir){
    if(bdry_phase[dir] != 0.)no_twist = 0;
  }

  if(!no_twist){

    /* No double wraparound supported */
    
    if(nx < 3 || ny < 3 ||  nz < 3 || nt < 3){
      node0_printf("%s: nonsense for dimensions less than 3.\n", myname);
    }
    

    if(*status_now == OFF){
      FORALLUPDIR(dir){
	cphase[dir] = ce_itheta(PI*bdry_phase[dir]);
      }
    }
    else { /* *status_now == ON */
      FORALLUPDIR(dir){
	cphase[dir] = ce_itheta(-PI*bdry_phase[dir]);
      }
    }

    apply_twist(fn, r0, cphase);

  }

  if(r0[XUP]%2 != 0 || r0[YUP]%2 != 0 || r0[ZUP]%2 != 0 || r0[TUP] != 0){
    
    if(*status_now == OFF){
      /* Switch on the unconventional KS phases and antiperiodic boundary condition */
      KS_rephase(fn, r0);
    } else {
      /* Switch on the conventional KS phases and antiperiodic boundary condition */
      KS_rephase(fn, zero);
    }
  }
  
  *status_now = flag;
}

