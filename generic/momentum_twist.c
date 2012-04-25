/***************** momentum_twist.c *******************************/
/* MIMD version 7*/

#include "generic_includes.h"

/* Apply phase change to gauge field. */

/* All link matrices in direction mu are multiplied by the same
   phase exp(i Pi bdry_phase[mu]/n[mu]) where n[mu] is the lattice
   extent in the mu direction */

/* This transformation is U(1) gauge equivalent to multiplying the
   last links in direction mu by exp(i Pi bdry_phase[mu]) */

/* Multiply an SU(3) matrix by a complex number */
static void 
phase_mult_su3_matrix( su3_matrix *a, complex c ){
  int i,j;
  complex z;
  for(i=0;i<3;i++)for(j=0;j<3;j++){
      CMUL(a->e[i][j], c, z);
      a->e[i][j] = z;
    }
}

/* Apply phase change to gauge field in site structure */
/* Phase is applied to ALL LINKS, divided equally in each direction */

void 
momentum_twist_site(Real bdry_phase[4], int sign) {
  //char myname[] = "momentum_twist_site";
  int i,dir;
  site *s;
  complex cphase[4];
  int no_twist;
  int nmu[4] = {nx, ny, nz, nt};

  /* The phases */

  no_twist = 1;
  FORALLUPDIR(dir){
    if(bdry_phase[dir] != 0.)no_twist = 0;
    cphase[dir] = ce_itheta(sign*PI*bdry_phase[dir]/nmu[dir]);
  }

  /* No twist needed if momentum phases are all zero */

  if(no_twist)
    return;

  /* Now apply twist to link member of the site structure */
  
  node0_printf("Applying momentum phases %g %g %g %g to gauge field\n",
	       sign*bdry_phase[0], sign*bdry_phase[1], 
	       sign*bdry_phase[2], sign*bdry_phase[3]);
  
  FORALLSITES(i,s){
    phase_mult_su3_matrix(&(s->link[XUP]), cphase[XUP]);
    phase_mult_su3_matrix(&(s->link[YUP]), cphase[YUP]);
    phase_mult_su3_matrix(&(s->link[ZUP]), cphase[ZUP]);
    phase_mult_su3_matrix(&(s->link[TUP]), cphase[TUP]);
  }
  
} /* momentum_twist_site */

/* Apply phase change to gauge field in site structure */
/* Phase is applied to the BOUNDARY LINKS */

void 
boundary_twist_site(Real bdry_phase[4], int r0[4], int sign) {
  //char myname[] = "boundary_twist_site";
  int i,dir;
  site *s;
  int bc_coord[4];
  complex cphase[4];
  int no_twist;

  /* The phases */

  no_twist = 1;
  FORALLUPDIR(dir){
    if(bdry_phase[dir] != 0.)no_twist = 0;
    cphase[dir] = ce_itheta(sign*PI*bdry_phase[dir]);
  }

  /* No twist needed if boundary phases are all zero */

  if(no_twist)
    return;

  /* Now apply the twist to the link member of the site structure */
  
  node0_printf("Applying boundary phases %g %g %g %g to gauge field\n",
	       sign*bdry_phase[0], sign*bdry_phase[1], 
	       sign*bdry_phase[2], sign*bdry_phase[3]);
  
  bc_coord[XUP] = (r0[XUP] + nx - 1) % nx;
  bc_coord[YUP] = (r0[YUP] + ny - 1) % ny;
  bc_coord[ZUP] = (r0[ZUP] + nz - 1) % nz;
  bc_coord[TUP] = (r0[TUP] + nt - 1) % nt;

  FORALLSITES(i,s){
    if( s->x == bc_coord[XUP])
      phase_mult_su3_matrix(&s->link[XUP], cphase[XUP]);
    if( s->y == bc_coord[YUP])
      phase_mult_su3_matrix(&s->link[YUP], cphase[YUP]);
    if( s->z == bc_coord[ZUP])
      phase_mult_su3_matrix(&s->link[ZUP], cphase[ZUP]);
    if( s->t == bc_coord[TUP])
      phase_mult_su3_matrix(&s->link[TUP], cphase[TUP]);
  }

  FORALLUPDIR(dir){
    boundary_phase[dir] = bdry_phase[dir];
  }
  
} /* boundary_twist_site */

/*--------------------------------------------------------------------*/
/* Apply progressive phase change to color vector field */
/* Note, this is equivalent to the gauge field transformation above */

/* Apply Fourier phase to an SU(3) vector */
static void 
phase_mult_su3_vector( su3_vector *a, complex c ){
  int i;
  complex z;

  for(i=0;i<3;i++){
    CMUL(a->c[i], c, z);
    a->c[i] = z;
  }
}

void 
rephase_v_field(su3_vector *v, Real bdry_phase[4], int r0[4], int sign) {

  int i;
  site *s;
  complex c;
  Real p;

  /* No transformation if phases are all zero */
  if(bdry_phase[XUP] == 0. && bdry_phase[YUP] == 0. &&
     bdry_phase[ZUP] == 0. && bdry_phase[TUP] == 0.)
    return;

  node0_printf("Applying momentum phases %g %g %g %g sign %d to vector field\n",
	       bdry_phase[0], bdry_phase[1], bdry_phase[2], bdry_phase[3], sign);

  FORALLSITES(i,s){
    p = bdry_phase[XUP]*((int)(s->x + nx - r0[XUP]) % nx)/(Real)nx 
      + bdry_phase[YUP]*((int)(s->y + ny - r0[YUP]) % ny)/(Real)ny 
      + bdry_phase[ZUP]*((int)(s->z + nz - r0[ZUP]) % nz)/(Real)nz 
      + bdry_phase[TUP]*((int)(s->t + nt - r0[TUP]) % nt)/(Real)nt;
    c = ce_itheta(sign*PI*p);
    phase_mult_su3_vector(v+i, c);
  }  
}

/*--------------------------------------------------------------------*/
/* Apply progressive phase change to a Wilson vector field */
/* Note, this is equivalent to the gauge field transformation above */

/* Apply Fourier phase to a Wilson vector */
static void 
phase_mult_wilson_vector( wilson_vector *a, complex c ){
  int i,j;
  complex z;

  for(j=0;j<4;j++)
    for(i=0;i<3;i++){
      CMUL(a->d[j].c[i], c, z);
      a->d[j].c[i] = z;
    }
}

void 
rephase_wv_field(wilson_vector *wv, Real bdry_phase[4], int r0[4], int sign) {

  int i;
  site *s;
  complex c;
  Real p;

  /* No transformation if phases are all zero */
  if(bdry_phase[XUP] == 0. && bdry_phase[YUP] == 0. &&
     bdry_phase[ZUP] == 0. && bdry_phase[TUP] == 0.)
    return;

  node0_printf("Applying momentum phases %g %g %g %g sign %d to Wilson vector field\n",
	       bdry_phase[0], bdry_phase[1], bdry_phase[2], bdry_phase[3], sign);

  FORALLSITES(i,s){
    p = bdry_phase[XUP]*((int)(s->x + nx - r0[XUP]) % nx)/(Real)nx 
      + bdry_phase[YUP]*((int)(s->y + ny - r0[YUP]) % ny)/(Real)ny 
      + bdry_phase[ZUP]*((int)(s->z + nz - r0[ZUP]) % nz)/(Real)nz 
      + bdry_phase[TUP]*((int)(s->t + nt - r0[TUP]) % nt)/(Real)nt;
    c = ce_itheta(sign*PI*p);
    phase_mult_wilson_vector(wv+i, c);
  }  
}

