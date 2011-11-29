/************** fermion_links_fn_twist_milc.c **************************/
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
set_boundary_twist_fn(fn_links_t *fn, Real bdry_phase[4], int r0[4]){
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

/* Change sign of su3_matrix */
static void 
negate_su3_matrix(su3_matrix *link){
  int j,k;
  for(j=0;j<3;j++)for(k=0;k<3;k++){
      CNEGATE(link->e[j][k],link->e[j][k]);
    }
}

/* Relocate the time antiperiodic boundary condition from the old
   origin t0_old to the new origin t0_new */

#define BOOL_NEQ(a,b) ((((a) && !(b))) || ((!(a) && (b))))

static void 
switch_time_apbc(fn_links_t *fn, int t0_old, int t0_new){
  su3_matrix *fat = get_fatlinks(fn);
  su3_matrix *lng = get_lnglinks(fn);
  su3_matrix *fatback = get_fatbacklinks(fn);
  su3_matrix *lngback = get_lngbacklinks(fn);
  int t_old, t_new;
  int i;
  site *s;
  
  /* Nothing to do if the times are the same */
  if(t0_old == t0_new)return;

  FORALLSITES(i,s){
    t_old = (s->t + nt - t0_old) % nt;
    t_new = (s->t + nt - t0_new) % nt;
    if(BOOL_NEQ( t_old == nt-1, t_new == nt-1 ))
      negate_su3_matrix(fat+4*i+TUP);
    if(lng != NULL && 
       BOOL_NEQ( t_old >= nt-3, t_new >= nt-3 ))
      negate_su3_matrix(lng+4*i+TUP);
    if(fatback != NULL && 
       BOOL_NEQ( t_old == 0, t_new == 0 ))
      negate_su3_matrix(fatback+4*i+TUP);
    if(lngback != NULL && 
       BOOL_NEQ( t_old <= 2, t_new <= 2))
      negate_su3_matrix(lngback+4*i+TUP);
  }
}

/* Compute the hypercube coordinate relative to an offset.  We assume
   that all lattice dimensions are even, as they should be for
   staggered fermions! */
static short *
hyp_coord(site *s, int r0[]){
  static short h[4];
  h[XUP] = (s->x - r0[XUP]) & 0x1;
  h[YUP] = (s->y - r0[YUP]) & 0x1;
  h[ZUP] = (s->z - r0[ZUP]) & 0x1;
  h[TUP] = (s->t - r0[TUP]) & 0x1;
  return h;
}

/* STANDARD MILC STAGGERED PHASES */
/*    eta_t = 1
      eta_x = (-1)^t
      eta_y = (-1)^(t+x)
      eta_z = (-1)^(t+x+y)
*/
/*	phase of link(i,mu) = prod(nu<mu) { -1^i[nu] }		*/
/*	all t phases for t=nt-1 time slice get extra minus sign	*/
/*	   to give antiperiodic boundary conditions		*/
       
void 
alpha_offset(short phase[], site *sit, int r0[]){
  short *h = hyp_coord(sit, r0);

  phase[TUP] = 1;
  if( h[TUP]%2 == 1) phase[XUP] = -1;  else phase[XUP] = 1;
  if( h[XUP]%2 == 1) phase[YUP] = -phase[XUP]; else phase[YUP] = phase[XUP];
  if( h[YUP]%2 == 1) phase[ZUP] = -phase[YUP]; else phase[ZUP] = phase[YUP];
}


/* Switch KS phases from origin r0_old to origin r0_new */

static void 
switch_KS_phases(fn_links_t *fn, int r0_old[], int r0_new[]){
  su3_matrix *fat = get_fatlinks(fn);
  su3_matrix *lng = get_lnglinks(fn);
  su3_matrix *fatback = get_fatbacklinks(fn);
  su3_matrix *lngback = get_lngbacklinks(fn);
  int i, dir;
  short p_old[4], p_new[4];
  site *s;

  /* Nothing to do if we are just shifting 2^4 hypercubes */

  if((r0_old[0]-r0_new[0])%2 == 0 && 
     (r0_old[1]-r0_new[1])%2 == 0 && 
     (r0_old[2]-r0_new[2])%2 == 0 &&
     (r0_old[3]-r0_new[3])%2 == 0)
    return;
  
  FORALLSITES(i,s){
    alpha_offset(p_old, s, r0_old);
    alpha_offset(p_new, s, r0_new);
    FORALLUPDIR(dir){
      if(p_old[dir]*p_new[dir] == -1){
	negate_su3_matrix(fat+4*i+dir);
	if(lng != NULL) negate_su3_matrix(lng+4*i+dir);
	if(fatback != NULL) negate_su3_matrix(fatback+4*i+dir);
	if(lngback != NULL) negate_su3_matrix(lngback+4*i+dir);
      }
    }
  }
}

/*--------------------------------------------------------------------*/
/* Do the momentum twist on the FN links (AKA boundary twist)         */
/*--------------------------------------------------------------------*/

/* Apply complex phase */
static void 
phase_mult_su3_matrix( su3_matrix *a, complex c ){
  int i,j;
  complex z;
  for(i=0;i<3;i++)for(j=0;j<3;j++){
      CMUL(a->e[i][j], c, z);
      a->e[i][j] = z;
    }
}

/* Remove complex phase */
static void 
phase_div_su3_matrix( su3_matrix *a, complex c ){
  int i,j;
  complex z;
  for(i=0;i<3;i++)for(j=0;j<3;j++){
      CMUL_J(a->e[i][j], c, z);
      a->e[i][j] = z;
    }
}

static void 
apply_twist(fn_links_t *fn, int r0[], complex cphase[]){

  su3_matrix *fat = get_fatlinks(fn);
  su3_matrix *lng = get_lnglinks(fn);
  su3_matrix *fatback = get_fatbacklinks(fn);
  su3_matrix *lngback = get_lngbacklinks(fn);
  int xp, yp, zp, tp;
  int i;
  site *s;

  FORALLSITES(i,s){
    xp = (s->x + nx - r0[XUP]) % nx;
    yp = (s->y + ny - r0[YUP]) % ny;
    zp = (s->z + nz - r0[ZUP]) % nz;
    tp = (s->t + nt - r0[TUP]) % nt;

    /* fat links: multiply only on surface */
    if( xp == nx - 1 )
      phase_mult_su3_matrix(fat+4*i+XUP, cphase[XUP]);
    if( yp == ny - 1 )
      phase_mult_su3_matrix(fat+4*i+YUP, cphase[YUP]);
    if( zp == nz - 1 )
      phase_mult_su3_matrix(fat+4*i+ZUP, cphase[ZUP]);
    if( tp == nt - 1 )
      phase_mult_su3_matrix(fat+4*i+TUP, cphase[TUP]);
    
    /* long links: multiply three deep on surface */
    if( xp >= nx-3 )
      phase_mult_su3_matrix(lng+4*i+XUP, cphase[XUP]);
    if( yp >= ny-3 )
      phase_mult_su3_matrix(lng+4*i+YUP, cphase[YUP]);
    if( zp >= nz-3 )
      phase_mult_su3_matrix(lng+4*i+ZUP, cphase[ZUP]);
    if( tp >= nt-3 )
      phase_mult_su3_matrix(lng+4*i+TUP, cphase[TUP]);
    
    /* Backward links are stored as the adjoint */

    /* fatback links */
    if(fatback != NULL){
      if( xp == 0 )
	phase_div_su3_matrix(fatback+4*i+XUP, cphase[XUP]);
      if( yp == 0 )
	phase_div_su3_matrix(fatback+4*i+YUP, cphase[YUP]);
      if( zp == 0 )
	phase_div_su3_matrix(fatback+4*i+ZUP, cphase[ZUP]);
      if( tp == 0 )
	phase_div_su3_matrix(fatback+4*i+TUP, cphase[TUP]);
    }
    
    /* longback links */
    if(lngback != NULL){
      if( xp <= 2)
	phase_div_su3_matrix(lngback+4*i+XUP, cphase[XUP]);
      if( yp <= 2)
	phase_div_su3_matrix(lngback+4*i+YUP, cphase[YUP]);
      if( zp <= 2)
	phase_div_su3_matrix(lngback+4*i+ZUP, cphase[ZUP]);
      if( tp <= 2)
	phase_div_su3_matrix(lngback+4*i+TUP, cphase[TUP]);
    }
  }
}

void 
boundary_twist_fn(fn_links_t *fn, int flag) {
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

  /* Shift the KS phases and antiperiodic BC, depending on the offset */

  if(*status_now == OFF){
    /* Shift the origin from zero to r0 */
    switch_time_apbc(fn, 0, r0[3]);
    switch_KS_phases(fn, zero, r0);
  } else {
    /* Shift the origin from r0 back to zero */
    switch_time_apbc(fn, r0[3], 0);
    switch_KS_phases(fn, r0, zero);
  }


  /* No twist needed if boundary phases are all zero */

  no_twist = 1;
  FORALLUPDIR(dir){
    if(bdry_phase[dir] != 0.)no_twist = 0;
  }

  if(no_twist){
    *status_now = flag;
    return;
  }

  /* No double wraparound supported */

  if(nx < 3 || ny < 3 ||  nz < 3 || nt < 3){
    node0_printf("%s: twisting unsupported for dimensions less than 3.\n", myname);
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

  /* Now apply twist to fat and lng members */
  /* and fatback AND longback members if used */

  apply_twist(fn, r0, cphase);

  *status_now = flag;
}
