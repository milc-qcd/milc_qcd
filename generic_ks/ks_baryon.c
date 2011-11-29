/********************** ks_baryon.c *********************************/
/* MIMD version 7 */

/*  Tie together staggered propagators to make baryon correlators
 *  (nucleon and delta)
 *
 */

/* Delta propagator code taken from nl_spectrum.c */
/* The method for the delta is described in M.F.L. Golterman and */
/* J. Smit, Nucl. Phys. B 255, 328 (1985)                        */
/* Equation (6.3) defines the sink operator for the delta        */

/* C.E.D.  7/91  first verion
 * D.T. modified 5/97 for improved action program, changed source
 *   slice logic, use gen_pt[8-15] instead of gen_pt2[].
 * For measuring propagation IN THE T DIRECTION
 * With correct t offset in E_MES_PRO.  Now reported as E_PI_PRO. C.D.7/17/91
 * modified 2/99 for new output format, DT
 *
 * Spectrum for Kogut-Susskind nonlocal hadrons, including delta
 * pion, and a few local ones for checking.  Uses E and O wall sources
 * This version does arbitrary number of wall sources 
 * This version DOES NOT FIX THE GAUGE.  MUST BE DONE BEFORE CALLING!!
 */

#include "generic_ks_includes.h"
#include <string.h>

/*------------------------------------------------------------------*/
/* Symbolic names for propagators.  prop_SOURCE_SINK */
/* Name by KS flavor content e.g. Goldstone pion = pion5 */
/* operators:
   pion5:	local 0-+:  (flavor)gamma_5     partner=0+-  phase=(1)
   pion05:	local 0-+:  gamma_0 gamma_5     partner=0++  phase=(-1)^(x+y+z+t)
   pioni5:	one-link 0-+:  gamma_i gamma_5	partner=0+-
   pionij:	one-link 0-+:  gamma_i gamma_j  partner=0++
   pioni:	two-link 0-+:  gamma_i          partner=0++
   pioni0:	two-link 0-+:  gamma_i gamma_0	partner=0+-
   pions:	three-link 0-+:  1 ("singlet")  partner=0++
   pion0:	three-link 0-+:  gamma_0	partner=0+-

   rhoi:	local 1--: gamma_i              partner=1+-  phase=(-1)^(dir) (VT)
   rhoi0:	local 1--: gamma_i gamma_0      partner=1++  phase=(-1)^(x+y+z+t+dir) (PV)
   rhos:	one-link 1--: 1 ("singlet")     partner=1+-
   rho0:	one-link 1--: gamma_0           partner=1++
*/

enum baryon_type_type {
  nucleon,
  delta,
  MAX_BARYON_TYPE
};

static char *baryon_type_string[MAX_BARYON_TYPE]  = { 
    "nucleon",
    "delta",
};

/*------------------------------------------------------------------*/
/* Map a label to the baryon-type index */

int baryon_type_index(char *label){
  int i;
  for(i = 0; i < MAX_BARYON_TYPE; i++){
    if(strcmp(label,baryon_type_string[i]) == 0)return i;
  }
  return -1;  /* Error condition */
}

/* Map an index to the label */

char *baryon_type_label(int index){
  return baryon_type_string[index];
}


static void accum_delta_prop(int i, int c[3], su3_matrix *tmp)
     /* Load a determinant for the calculation of the delta propagator */
     /* Each column of the determinant is the symmetrically shifted   */
     /* "q" propagator.  Each column has a differenct shift direction */
     /* Each column also corresponds to a different source wall color */
     /* as specified by the color vector c                            */
     /* Result is put in tempmat1 */
{
  int j=-1,dir;
  su3_matrix *af, *ab;
  
  FORALLUPDIRBUT(TUP,dir)
    {
      /* For timeward propagation, this step is pedantic */
      /* since we could simply take j = dir */
      /* but for spacelike propagation, it is necessary */
      switch (dir)
	{
	case XUP: j = 0; break;
	case YUP: j = 1; break;
	case ZUP: j = 2;
	}

      /* gen_pt points to the "q" propagator on the next neighbor    */

      af = (su3_matrix *)gen_pt[dir][i];
      ab = (su3_matrix *)gen_pt[dir+4][i];
      add_su3_vector((su3_vector *)(ab->e[c[j]]),
		     (su3_vector *)(af->e[c[j]]),
		     (su3_vector *)(tmp->e[j]) );

    }
}

static void delta_prop(int c0,int c1,int c2,int perm,int t,Real *delprop)

/* Calculates the contribution to the delta propagator at time t */
/* With source colors c0, c1, c2 */
/* perm specifies the sign of the contribution to be calculated */
{
  register int i;
  register complex cc;
  int x,y,z;
  int c[3];
  su3_matrix tempmat;

  c[0] = c0;
  c[1] = c1;
  c[2] = c2;

  /*  Calculate for the cube origin only */
  for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
	if( node_number(x,y,z,t) != mynode() )continue;
	i=node_index(x,y,z,t);
	
	accum_delta_prop(i,c,&tempmat);
	
	cc = det_su3( &tempmat );
	
	if(perm>0)*delprop += cc.real;
	else      *delprop -= cc.real;  
	
      }
}

#include <assert.h>

/* The input ks_prop_field must come from an even_and_odd wall source */

static void ks_delta( ks_prop_field *qp, complex prop[]){

  site* s;
  int i,c0,c1,dir,t;
  su3_matrix *propmat;
  msg_tag *mtag[16];

  if(qp->nc != 3){
    node0_printf("ks_delta: baryon propagators must all have only three colors\n");
    return;
  }

  /* Required gen_pt dimension */
  assert(N_POINTERS >= 8);

  /* Map ks_prop_field to temporary */
  propmat = create_m_field();

  FORALLSITES(i,s){
    for(c0 = 0; c0 < 3; c0++)
      for(c1 = 0; c1 < 3; c1++)
	propmat[i].e[c0][c1] = qp->v[c0][i].c[c1];
  }

  /* Do forward and backward shifts needed for symmetric shifts Dq */
      
  FORALLUPDIRBUT(TUP,dir) {
    
    mtag[dir] = start_gather_field(propmat, sizeof(su3_matrix), dir, 
				   EVENANDODD, gen_pt[dir]);
	  
    mtag[dir+4] = start_gather_field(propmat, sizeof(su3_matrix), OPP_DIR(dir), 
				     EVENANDODD, gen_pt[dir+4]);
    wait_gather(mtag[dir]);
    wait_gather(mtag[dir+4]);
    
  }
      
  
  /* Calculate and dump delta propagator */
  for(t=0; t<nt; t++) {

    prop[t].imag = 0;

    /* Calculate contribution for each permutation of source color */
    delta_prop (0,1,2, 1, t, &prop[t].real);
    delta_prop (1,2,0, 1, t, &prop[t].real);
    delta_prop (2,0,1, 1, t, &prop[t].real);
    delta_prop (1,0,2,-1, t, &prop[t].real);
    delta_prop (0,2,1,-1, t, &prop[t].real);
    delta_prop (2,1,0,-1, t, &prop[t].real);

  }
      
  /* Clean up gathers */
  FORALLUPDIRBUT(TUP,dir) {
    cleanup_gather(mtag[dir]);
    cleanup_gather(mtag[dir+4]);
  }
  
  destroy_m_field(propmat);
  
} /* ks_delta */

#if 0
static void ks_nucleon(ks_prop_field *qp, complex prop[]){
  
  int i, c0, c1;
  int x,y,z,t;
  complex cc;
  su3_matrix propmat;

  if(qp->nc != 3){
    node0_printf("ks_nucleon: baryon propagators must all have only three colors\n");
    return;
  }

  for(t=0; t<nt; t++) {
    for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
	  if( node_number(x,y,z,t) != this_node )continue;
	  i=node_index(x,y,z,t);
	  for(c0 = 0; c0 < 3; c0++)
	    for(c1 = 0; c1 < 3; c1++)
	      propmat.e[c0][c1] = qp->v[c0][i].c[c1];
	  cc = det_su3( &propmat );
	  CSUM(prop[t], cc);
	}
  }
}
#endif

static void ks_nucleon_nd(ks_prop_field *qp0, ks_prop_field *qp1, ks_prop_field *qp2,
			  complex prop[]){

  int i, c;
  int x,y,z,t;
  complex cc;
  su3_matrix propmat;

  if(qp0->nc != 3 || qp1->nc != 3 || qp2->nc != 3){
    node0_printf("ks_nucleon_nd: baryon propagators must all have only three colors\n");
    return;
  }

  for(t=0; t<nt; t++) {
    for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
	  if( node_number(x,y,z,t) != this_node )continue;
	  i=node_index(x,y,z,t);
	  for(c = 0; c < 3; c++){
	    propmat.e[0][c] = qp0->v[0][i].c[c];
	    propmat.e[1][c] = qp1->v[1][i].c[c];
	    propmat.e[2][c] = qp2->v[2][i].c[c];
	  }
	  cc = det_su3( &propmat );
	  CSUM(prop[t], cc);
	}
  }
}

/*---------------------------------------------------------------------*/

static void 
norm_corr(int phase, Real fact, complex prop[])
{
  int k;
  complex z = {0.,0.};
  
  for(k=0; k<nt; k++){
    switch(phase){
    case 0:
      z =            prop[k];
      break;
    case 1:
      TIMESPLUSI(    prop[k], z);
      break;
    case 2:
      TIMESMINUSONE( prop[k], z);
      break;
    case 3:
      TIMESMINUSI(   prop[k], z);
    }
    CMULREAL(z,fact,prop[k]);
  }
} /* norm */

/*---------------------------------------------------------------------*/

void ks_baryon_nd(complex *prop[],
		  ks_prop_field *qp0, ks_prop_field *qp1, ks_prop_field *qp2,
		  int num_corr_b, int baryon_type_snk[], int phase[], Real fact[]){
  int i;

  for(i = 0; i < num_corr_b; i++){
    if(baryon_type_snk[i] == nucleon)
      ks_nucleon_nd(qp0, qp1, qp2, prop[i]);

    else if(baryon_type_snk[i] == delta){
      if(qp0 != qp1 || qp0 != qp2){
	printf("WARNING: No provision for non-degenerate deltas!\n");
      }
      ks_delta(qp0, prop[i]);
    }
    else{
      printf("ERROR: Bad baryon sink type %d\n",baryon_type_snk[i]);
    }
    norm_corr( phase[i], fact[i], prop[i] );
  }
}
