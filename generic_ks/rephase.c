/********************* rephase.c *******************************/
/* MIMD version 7*/

#include "generic_ks_includes.h"

/********* phaseset() - set up KS phase vectors **********/
/* ANTIPERIODIC bc's in t and PERIODIC in x,y,z */

/* The default phase set has standard KS phases with periodic BC in
   the three spatial directions and antiperiodic or periodic in
   time. */

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
/*      with mu = 0,1,2,3 for t,x,y,z                           */
/*	phase of link(i,mu) = prod(nu<mu) { -1^i[nu] }		*/
/*	all t phases for t=nt-1 time slice get extra minus sign	*/
/*	   to give antiperiodic boundary conditions		*/

static short *
alpha_offset(site *sit, int r0[]){
  static short phase[4];
  short *h = hyp_coord(sit, r0);

  phase[TUP] = 1;
  if( h[TUP]%2 == 1) phase[XUP] = -1;  else phase[XUP] = 1;
  if( h[XUP]%2 == 1) phase[YUP] = -phase[XUP]; else phase[YUP] = phase[XUP];
  if( h[YUP]%2 == 1) phase[ZUP] = -phase[YUP]; else phase[ZUP] = phase[YUP];

  return phase;
}

#if 0
static short *
alpha_apb_offset(site *sit, int r0[]){

  short *phase = alpha_offset(sit, r0);

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
  
  FORALLSITES(i,sit){
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
  }
}


/************************** rephase() ******************************/
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
  FORALLSITES(i,s){
    for(dir=XUP;dir<=TUP;dir++){
      for(j=0;j<3;j++)for(k=0;k<3;k++){
	  s->link[dir].e[j][k].real *= s->phase[dir];
	  s->link[dir].e[j][k].imag *= s->phase[dir];
	}
    }
  }

  phases_in = flag;
} /* rephase.c */

/************************** rephase_field() ******************************/
/* put Kogut-Sussind and boundary condition phase factors into or
   out of a field consisting of four SU(3) matrices per site.  The
   phases are defined relative to the offset origin r0. */
void rephase_field_offset( su3_matrix *internal_links, int flag, 
			   int* status_now, int r0[] ){
  register int i,j,k,dir;
  register site *s;

  /* Check to make sure we are going in expected direction */
  if( status_now != NULL)
    if( !( flag==ON && *status_now==OFF )  && 
	!( flag==OFF && *status_now==ON ) ){
      node0_printf("rephase_field: DUMMY: you fouled up the phases\n");
      terminate(1);
    }
  FORALLSITES(i,s){
    short *p = alpha_offset(s, r0);
    for(dir=XUP;dir<=TUP;dir++){
      for(j=0;j<3;j++)for(k=0;k<3;k++){
        (internal_links[4*i+dir].e[j][k]).real *= p[dir];
        (internal_links[4*i+dir].e[j][k]).imag *= p[dir];
      }
    }
  }

  if(status_now != NULL)
    *status_now = flag;

} /* rephase_field_offset.c */

// /************************** rephase_offset() ******************************/
// /* put Kogut-Sussind and anti-periodic boundary condition phase factors into or
//    out of lattice.  Phases are defined relative to the origin r0.*/
// 
// void 
// rephase_offset( int flag, int r0[] ){
//   register int i,j,k,dir;
//   register site *s;
//   char myname[] = "rephase_offset";
//   
//   /* Check to make sure we are going in expected direction */
//   if( !( (flag==ON && phases_in==OFF) || 
// 	 (flag==OFF && phases_in==ON) ) ){
//     node0_printf("%s: you fouled up the phases\n",myname);
//     terminate(1);
//   }
//   FORALLSITES(i,s){
//     short *p = alpha_apb_offset(s, r0);
//     for(dir=XUP;dir<=TUP;dir++){
//       for(j=0;j<3;j++)for(k=0;k<3;k++){
// 	  s->link[dir].e[j][k].real *= p[dir];
// 	  s->link[dir].e[j][k].imag *= p[dir];
// 	}
//     }
//   }
//   
//   phases_in = flag;
// } /* rephase_offset.c */

