/********************* rephase.c *******************************/
/* MIMD version 7*/

#include "generic_ks_includes.h"
#include "../include/openmp_defs.h"

/********* phaseset() - set up KS phase vectors **********/
/* ANTIPERIODIC bc's in t and PERIODIC in x,y,z */

/* The default phase set has standard KS phases with periodic BC in
   the three spatial directions and antiperiodic or periodic in
   time. */

/* Compute the hypercube coordinate relative to an offset.  We assume
   that all lattice dimensions are even, as they should be for
   staggered fermions! */
static void
hyp_coord(short h[], const site *s, const int r0[]){
  h[XUP] = (s->x - r0[XUP]) & 0x1;
  h[YUP] = (s->y - r0[YUP]) & 0x1;
  h[ZUP] = (s->z - r0[ZUP]) & 0x1;
  h[TUP] = (s->t - r0[TUP]) & 0x1;
}

/* STANDARD MILC STAGGERED PHASES */
/*    eta_t = 1
      eta_x = (-1)^t
      eta_y = (-1)^(t+x)
      eta_z = (-1)^(t+x+y)
*/
/*      with mu = 0,1,2,3 for t,x,y,z                           */
/*	phase of link(i,mu) = prod(nu<mu) { -1^i[nu] }		*/
/*	all t phases for t=nt-1 time slice get extra minus sign	*/
/*	   to give antiperiodic boundary conditions		*/

static void
alpha_offset(short h[], short phase[], const site *sit, const int r0[]){

  hyp_coord(h, sit, r0);

  phase[TUP] = 1;
  if( h[TUP]%2 == 1) phase[XUP] = -1;  else phase[XUP] = 1;
  if( h[XUP]%2 == 1) phase[YUP] = -phase[XUP]; else phase[YUP] = phase[XUP];
  if( h[YUP]%2 == 1) phase[ZUP] = -phase[YUP]; else phase[ZUP] = phase[YUP];
}

#if 0
static void
alpha_apb_offset(short phase[], site *sit, int r0[]){

  short h[4];

  alpha_offset(h, phase, sit, r0);

  /* antiperiodic boundary conditions in Euclidean time */
  if( ((sit->t - r0[TUP] + nt ) % nt) == nt-1) {
    phase[TUP] = -phase[TUP];
  }

  return phase;
}
#endif

void phaseset() {
  register site *sit; /* pointer to current site */
  register int i;

  /* This choice can be superceded by the boundary_twist parameter */
  //#ifdef PERIODICBC
  //  node0_printf("with periodic boundary conditions in time\n");
  //#endif

  /* STANDARD MILC PHASES */
  /*    eta_t = 1
	eta_x = (-1)^t
	eta_y = (-1)^(t+x)
	eta_z = (-1)^(t+x+y)
  */
  /*	phase of link(i,mu) = prod(nu<mu) { -1^i[nu] }		*/
  /*	all t phases for t=nt-1 time slice get extra minus sign	*/
  /*	   to give antiperiodic boundary conditions		*/
  
  FORALLSITES_OMP(i,sit,default(shared)){
    sit->phase[TUP] = 1.0;
    if( (sit->t)%2 == 1) sit->phase[XUP] = -1.0;
    else sit->phase[XUP] = 1.0;
    if( (sit->x)%2 == 1) sit->phase[YUP] = -sit->phase[XUP];
    else sit->phase[YUP] = sit->phase[XUP];
    if( (sit->y)%2 == 1) sit->phase[ZUP] = -sit->phase[YUP];
    else sit->phase[ZUP] = sit->phase[YUP];

//    /* ELVIRA'S CODE PHASES */
//
//    /*    eta_x = 1
//	  eta_y = (-1)^x
//	  eta_z = (-1)^(x+y)
//	  eta_t = (-1)^(x+y+z)
//    */
//
//  FORALLSITES(i,sit){
//    sit->phase[XUP] = 1.0;
//    if( (sit->x)%2 == 1) sit->phase[YUP] = -1.0;
//    else sit->phase[YUP] = 1.0;
//    if( (sit->y)%2 == 1) sit->phase[ZUP] = -sit->phase[YUP];
//    else sit->phase[ZUP] = sit->phase[YUP];
//    if( (sit->z)%2 == 1) sit->phase[TUP] = -sit->phase[ZUP];
//    else sit->phase[TUP] = sit->phase[ZUP];
    
//#ifndef PERIODICBC
    if( sit->t == nt-1) {
      /* antiperiodic boundary conditions in Euclidean time */
      sit->phase[TUP] = -sit->phase[TUP];
    }
    //#endif
  } END_LOOP_OMP
}


/* put Kogut-Sussind and boundary condition phase factors into or
   out of lattice */

void rephase( int flag ){
  register int i,j,k,dir;
  register site *s;
  /* Check to make sure we are going in expected direction */
  if( !( (flag==ON && phases_in==OFF) || (flag==OFF && phases_in==ON) ) ){
    node0_printf("rephase: DUMMY: you fouled up the phases\n");
    terminate(1);
  }
  FORALLSITES_OMP(i,s,private(dir,j,k)){
    for(dir=XUP;dir<=TUP;dir++){
      for(j=0;j<3;j++)
	for(k=0;k<3;k++){
	  lattice[i].link[dir].e[j][k].real *= lattice[i].phase[dir];
	  lattice[i].link[dir].e[j][k].imag *= lattice[i].phase[dir];
	}
    }
  } END_LOOP_OMP;

  phases_in = flag;
} /* rephase */

/* put Kogut-Sussind and boundary condition phase factors into or
   out of a field consisting of four SU(3) matrices per site.  The
   phases are defined relative to the offset origin r0. 
*/
void rephase_field_offset( su3_matrix *internal_links, int flag, 
			   int* status_now, const int r0[] ){
  register int i,j,k,dir;
  register site *s;
  short h[4] __attribute__ ((aligned (8)));
  short p[4] __attribute__ ((aligned (8)));

  /* Check to make sure we are going in expected direction */
  if( status_now != NULL)
    if( !( flag==ON && *status_now==OFF )  && 
	!( flag==OFF && *status_now==ON ) ){
      node0_printf("rephase_field: DUMMY: you fouled up the phases\n");
      terminate(1);
    }
  FORALLSITES_OMP(i,s,private(dir,j,k,h,p)){
    alpha_offset(h, p, s, r0);
    for(dir=XUP;dir<=TUP;dir++){
      for(j=0;j<3;j++)for(k=0;k<3;k++){
        (internal_links[4*i+dir].e[j][k]).real *= p[dir];
        (internal_links[4*i+dir].e[j][k]).imag *= p[dir];
      }
    }
  }
  END_LOOP_OMP;

  if(status_now != NULL)
    *status_now = flag;

} /* rephase_field_offset */

/* Remove any persistent phase */
/* Just avoids the error message about botched phases */
void rephase_field_unset(su3_matrix *links, int *status_now,
			 int r0_now[], int *changed){
  if(*status_now == OFF)return;
  node0_printf("Removing phases\n"); fflush(stdout);
  rephase_field_offset( links, OFF, status_now, r0_now );
  for(int i = 0; i < 4; i++)r0_now[i] = 0;
  *changed = 1;
}

/* Support persistent KS phases in a link field 
   Here we support a protocol that inserts phases and keeps them there
   unless they need to be changed.
*/

/* TODO: Need a structure containing the link field and its phase state */
void rephase_field_set(su3_matrix *links, const int flag, int *status_now,
		       int r0_now[], int *changed, const int r0[]){
  /* Skip setting if the phases are ON with a compatible offset */
  if(flag == OFF){
    rephase_field_unset(links, status_now, r0_now, changed);
  } else {
    if(*status_now == ON)
      {
	if(r0[0] % 2 != r0_now[0] % 2 ||
	   r0[1] % 2 != r0_now[1] % 2 ||
	   r0[2] % 2 != r0_now[2] % 2 ||
	   r0[3] % 2 != r0_now[3] % 2){
	  /* Incomptible: Redo the phase offset */
	  node0_printf("Replacing phases\n"); fflush(stdout);
	  rephase_field_offset( links, OFF, status_now, r0_now );
	  rephase_field_offset( links, ON, status_now, r0 );
	  for(int i = 0; i < 4; i++)r0_now[i] = r0[i];
	  *changed = 1;
	}
      } else {
      /* Insert the phases */
      node0_printf("Inserting phases\n"); fflush(stdout);
      rephase_field_offset( links, ON, status_now, r0 );
      for(int i = 0; i < 4; i++)r0_now[i] = r0[i];
      *changed = 1;
    }
  }
}

/* conventional antiperiodic boundary conditions in Euclidean time */
/* Do not use for long links! */
void apply_apbc( su3_matrix *links, int r0t ){

  int i;
  site *s;

  FORALLSITES_OMP(i,s,default(shared)){
    if( s->t == nt-1 - r0t % nt){
      scalar_mult_su3_matrix( links + 4*i + TUP, -1., links + 4*i + TUP );
    }
  }
  END_LOOP_OMP;
}
