/******** spectrum_multimom.c *************/
/* MIMD version 7*/
/* DT 4/01
   Rho masses for various quark masses and momenta.  Do many valence masses,
   one of which is presumably the dynamical mass.
   Also need pseudoscalars with one valence quark and one dynamical
   quark.  ("kaons")
*/
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.  Assumes "dslash_site" has
   been defined to be the appropriate "dslash_fn_site" or "dslash_eo_site"
*/

/* Spectrum for Kogut-Susskind hadrons.
   Any antiperiodic boundary conditions (time direction usually)
   are absorbed in the link matrices (ie rephase(ON) )

    Gauge fixing should be done before calling this function
    (Coulomb gauge, probably).
*/


/* Symbolic names for propagators.  prop__SINK */
/* operators:
   pion:	qbar gamma_5 x gamma_5 q:  sign = +1
   pion2:	qbar gamma_0 x gamma_5 q> sign =  (-1)^(x+y+z+t)
   rho:		1--: qbar gamma_i gamma_0 q
		(-1)^(dir) (VT)
   rho2:	1--: qbar gamma_i q
		(-1)^(x+y+z+t+dir) (PV)
*/
enum prop_name { 
    prop_kaon_000,	/* zero momentum kaon */
    prop_kaon_001,	/* one unit of momentum */
    prop_kaon_011,	/* momentum(1,1,0) */
    prop_rho_000,	/* zero momentum rho */
    prop_rho_001perp,	/* one unit of momentum, polarization perp. to mom */
    prop_rho_001para,	/* one unit of momentum, polarization parallel to mom */
    prop_rho_011perp,	/* one unit of momentum, polarization perp. to both mom */
    prop_rho_011para,	/* one unit of momentum, polarization parallel to one mom */

    nprops		/* nprops = number of propagators */
};

#include "generic_ks_includes.h"
#include "../include/dslash_ks_redefine.h"

/* Various meson operators.  Source code later in this file */
/* All of these operators contain an extra gamma_5, entering when
    antiquark propagators are computed from the same source as quark
    propagators */
   /* quark-antiquark operators*/
void mult_pion_mom_temp( int fb, int px, int py, int pz,
    su3_vector *src, field_offset dest );
void mult_rho_mom_temp( int fb, int pdir, int px, int py, int pz,
    su3_vector *src, field_offset dest );
int test_converge(int t_source);

int spectrum_multimom( Real dyn_mass, Real low_mass, Real mass_inc, int nmasses, Real tol, ferm_links_t *fn){
  /* arguments are dynamical mass, lowest mass, distance between masses, number of masses,
     tolerance for inverter check.
     return C.G. iteration number */

  int cgn;
  register int i,j,pdir,mdir;
  register site* s;
  register complex cc;
  register int t_source;
  int j_dyn;	/* index of dynamical mass */
  int color;	/* color for source */
  int src_count; /* number of source time slices used */
  Real *masses; /* array of masses */
  Real finalrsq;
  su3_vector **quark_props;
  complex **props;	/* arrays of propagators */
  int prec = PRECISION;  /* make CG precision the same as the
			    prevailing precision */

  cgn=0; /* number of CG iterations */

  /* array of masses */
  j_dyn=-1;
  masses = (Real *)malloc(nmasses*sizeof(Real));
  for(i=0; i<nmasses; i++){
    masses[i]=low_mass+i*mass_inc;
    if( fabs(masses[i]-dyn_mass) < 0.0001 )j_dyn=i;
  }
  if(j_dyn== -1){node0_printf("DUMMY: one mass should be %e\n",dyn_mass);}

  /* allocate space for quark and antiquark propagators */
  quark_props = (su3_vector **)malloc(nmasses*sizeof(su3_vector *));
  for(i=0; i<nmasses; i++){
      quark_props[i] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  }
  /* allocate space for meson propagators (one number per timeslice) */
  /* for meson number n, mass number m, use props[nmasses*n+m] */
  props = (complex **)malloc( nprops*nmasses*sizeof(complex *) );
  props[0] = (complex *)malloc( nprops*nmasses*nt*sizeof(complex) );
  for(i=1; i<nprops*nmasses; i++) props[i] = props[i-1] + nt;

  /* set propagators to zero */
  for(cc.real=cc.imag=0.0,i=0;i<nprops*nmasses;i++)for(j=0;j<nt;j++){
    props[i][j]=cc;
  }

  /* loop over "source" time slice */
  for(src_count=0,t_source=source_start; t_source<nt && src_count<n_sources;
    t_source += source_inc,src_count++){
    /* Corner wall source */
    /* Use quark_source for quark source */
    if(this_node==0)printf("spectrum_multimom(): source time = %d\n",t_source);

    for(color=0;color<3;color++){
	FORALLSITES(i,s){
	    clearvec( &(s->quark_source) );
if( t_source%2 != 0 ){node0_printf("Even sources only!\n"); terminate(0);}
	    if( s->t==t_source && (s->x%2)==(s->t)%2 && (s->y%2)==0 && (s->z%2)==0 ){
		 s->quark_source.c[color].real  = 1.0;
		 cc = ce_itheta(  2*PI*(s->x)/((Real)nx) ) ; CSUM( s->quark_source.c[color], cc );
		 cc = ce_itheta(  2*PI*(s->y)/((Real)ny) ) ; CSUM( s->quark_source.c[color], cc );
		 cc = ce_itheta(  2*PI*(s->z)/((Real)nz) ) ; CSUM( s->quark_source.c[color], cc );
		/**
		s->quark_source.c[color].real += cos( 2*PI*(s->x)/((Real)nx) );
		s->quark_source.c[color].real += cos( 2*PI*(s->y)/((Real)ny) );
		s->quark_source.c[color].real += cos( 2*PI*(s->z)/((Real)nz) );
		**/
	    }
	}

	/* compute M^-1 * quark_source */
	cgn += ks_multicg_mass( F_OFFSET(quark_source), quark_props, masses, 
				nmasses, niter, rsqprop, prec, EVEN, 
				&finalrsq, fn);
	/* Multiply by Madjoint. Note this assumes source on even sites only */
	/****** NEW CODE **/
	for(j=0;j<nmasses;j++){
	    dslash_field( quark_props[j], quark_props[j], ODD, fn );
	    FOREVENSITES(i,s){
		scalar_mult_su3_vector( &(quark_props[j][i]), 2.0*masses[j], &(quark_props[j][i]) );
	    }
	    FORODDSITES(i,s){
		scalar_mult_su3_vector( &(quark_props[j][i]), -1.0, &(quark_props[j][i]) );
	    }

	    FORALLSITES(i,s) s->ttt = quark_props[j][i];
	    check_invert( F_OFFSET(ttt), F_OFFSET(quark_source), masses[j],tol,
			  fn);
	} /* j=masses */
	/**** OLD CODE *
	for(j=0;j<nmasses;j++){
	    FOREVENSITES(i,s) s->ttt = quark_props[j][i];
	    dslash_site( F_OFFSET(ttt), F_OFFSET(ttt), ODD, fn );
	    scalar_mult_latvec( F_OFFSET(ttt),  2.0*masses[j], F_OFFSET(ttt), EVEN );
	    scalar_mult_latvec( F_OFFSET(ttt), -1.0, F_OFFSET(ttt), ODD );
	    FORALLSITES(i,s) quark_props[j][i] = s->ttt;

	    check_invert( F_OFFSET(ttt), F_OFFSET(quark_source), masses[j],tol, fn);
	} ** j=masses **
	*END OLD CODE **/


	/* 0-+ (kaon) PROPAGATORS ***************************************/
	/* pdir = polarization direction, mdir = momentum direction */
	for(j=0;j<nmasses;j++){
	        mult_pion_mom_temp( FORWARDS,0,0,0,
	            quark_props[j], F_OFFSET(g_rand) );
	        FORALLSITES(i,s){
	            cc = su3_dot( &(quark_props[j_dyn][i]), &(s->g_rand) );
	            CSUM( props[nmasses*prop_kaon_000+j][(s->t+nt-t_source)%nt], cc );
	        }
    
	        for( mdir=XUP; mdir<=ZUP; mdir++ ){
		    switch(mdir){ /* here mdir is momentum direction */
		        case XUP: mult_pion_mom_temp( FORWARDS,1,0,0,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case YUP: mult_pion_mom_temp( FORWARDS,0,1,0,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case ZUP: mult_pion_mom_temp( FORWARDS,0,0,1,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		    }
	            FORALLSITES(i,s){
	                cc = su3_dot( &(quark_props[j_dyn][i]), &(s->g_rand) );
	                CSUM( props[nmasses*prop_kaon_001+j][(s->t+nt-t_source)%nt], cc );
	            }
		    switch(mdir){ /* here mdir is momentum direction */
		        case XUP: mult_pion_mom_temp( FORWARDS,-1,0,0,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case YUP: mult_pion_mom_temp( FORWARDS,0,-1,0,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case ZUP: mult_pion_mom_temp( FORWARDS,0,0,-1,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		    }
	            FORALLSITES(i,s){
	                cc = su3_dot( &(quark_props[j_dyn][i]), &(s->g_rand) );
	                CSUM( props[nmasses*prop_kaon_001+j][(s->t+nt-t_source)%nt], cc );
	            }
		    switch(mdir){  /*  here mdir is the direction NOT used */
		        case XUP: mult_pion_mom_temp( FORWARDS,0,1,-1,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case YUP: mult_pion_mom_temp( FORWARDS,-1,0,1,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case ZUP: mult_pion_mom_temp( FORWARDS,1,-1,0,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		    }
	            FORALLSITES(i,s){
	                cc = su3_dot( &(quark_props[j_dyn][i]), &(s->g_rand) );
	                CSUM( props[nmasses*prop_kaon_011+j][(s->t+nt-t_source)%nt], cc );
	            }
		    switch(mdir){  /*  here mdir is the direction NOT used */
		        case XUP: mult_pion_mom_temp( FORWARDS,0,-1,1,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case YUP: mult_pion_mom_temp( FORWARDS,1,0,-1,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case ZUP: mult_pion_mom_temp( FORWARDS,-1,1,0,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		    }
	            FORALLSITES(i,s){
	                cc = su3_dot( &(quark_props[j_dyn][i]), &(s->g_rand) );
	                CSUM( props[nmasses*prop_kaon_011+j][(s->t+nt-t_source)%nt], cc );
	            }
	        } /* mdir */
	} /* j (masses)*/

	/* 1-- (rho) PROPAGATORS ***************************************/
	/* pdir = polarization direction, mdir = momentum direction */
	for(j=0;j<nmasses;j++){
	    for( pdir=XUP; pdir<=ZUP; pdir++){
	        mult_rho_mom_temp( FORWARDS,pdir,0,0,0,
	            quark_props[j], F_OFFSET(g_rand) );
	        FORALLSITES(i,s){
	            cc = su3_dot( &(quark_props[j][i]), &(s->g_rand) );
	            CSUM( props[nmasses*prop_rho_000+j][(s->t+nt-t_source)%nt], cc );
	        }
    
	        for( mdir=XUP; mdir<=ZUP; mdir++ ){
		    switch(mdir){ /* here mdir is momentum direction */
		        case XUP: mult_rho_mom_temp( FORWARDS,pdir,1,0,0,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case YUP: mult_rho_mom_temp( FORWARDS,pdir,0,1,0,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case ZUP: mult_rho_mom_temp( FORWARDS,pdir,0,0,1,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		    }
	            FORALLSITES(i,s){
	                cc = su3_dot( &(quark_props[j][i]), &(s->g_rand) );
		        if( mdir==pdir ){
	                    CSUM( props[nmasses*prop_rho_001para+j][(s->t+nt-t_source)%nt], cc );
		        }
		        else {
	                    CSUM( props[nmasses*prop_rho_001perp+j][(s->t+nt-t_source)%nt], cc );
		        }
	            }
		    switch(mdir){ /* here mdir is momentum direction */
		        case XUP: mult_rho_mom_temp( FORWARDS,pdir,-1,0,0,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case YUP: mult_rho_mom_temp( FORWARDS,pdir,0,-1,0,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case ZUP: mult_rho_mom_temp( FORWARDS,pdir,0,0,-1,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		    }
	            FORALLSITES(i,s){
	                cc = su3_dot( &(quark_props[j][i]), &(s->g_rand) );
		        if( mdir==pdir ){
	                    CSUM( props[nmasses*prop_rho_001para+j][(s->t+nt-t_source)%nt], cc );
		        }
		        else {
	                    CSUM( props[nmasses*prop_rho_001perp+j][(s->t+nt-t_source)%nt], cc );
		        }
	            }
		    switch(mdir){  /*  here mdir is the direction NOT used */
		        case XUP: mult_rho_mom_temp( FORWARDS,pdir,0,1,-1,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case YUP: mult_rho_mom_temp( FORWARDS,pdir,-1,0,1,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case ZUP: mult_rho_mom_temp( FORWARDS,pdir,1,-1,0,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		    }
	            FORALLSITES(i,s){
	                cc = su3_dot( &(quark_props[j][i]), &(s->g_rand) );
		        if( mdir==pdir ){
	                    CSUM( props[nmasses*prop_rho_011perp+j][(s->t+nt-t_source)%nt], cc );
		        }
		        else {
	                    CSUM( props[nmasses*prop_rho_011para+j][(s->t+nt-t_source)%nt], cc );
		        }
	            }
		    switch(mdir){  /*  here mdir is the direction NOT used */
		        case XUP: mult_rho_mom_temp( FORWARDS,pdir,0,-1,1,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case YUP: mult_rho_mom_temp( FORWARDS,pdir,1,0,-1,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		        case ZUP: mult_rho_mom_temp( FORWARDS,pdir,-1,1,0,
                	    quark_props[j], F_OFFSET(g_rand) ); break;
		    }
	            FORALLSITES(i,s){
	                cc = su3_dot( &(quark_props[j][i]), &(s->g_rand) );
		        if( mdir==pdir ){
	                    CSUM( props[nmasses*prop_rho_011perp+j][(s->t+nt-t_source)%nt], cc );
		        }
		        else {
	                    CSUM( props[nmasses*prop_rho_011para+j][(s->t+nt-t_source)%nt], cc );
		        }
	            }
	        } /* mdir */
	    } /* pdir */
	} /* j (masses)*/

    }
  } /* end loop on t_source */


  /* Sum propagator arrays over nodes */
  /* print out propagators */
  for(i=0;i<nmasses*nprops;i++){
      g_veccomplexsum( props[i] , nt );
      for(j=0;j<nt;j++){
          CDIVREAL(props[i][j],n_sources*nx*ny*nz,props[i][j]);
      }
  }
  if(this_node==0){


    for(i=0;i<nmasses;i++){
        printf("STARTPROP\n");
        printf("MASSES:  %.5e   %.5e\n",masses[i],masses[j_dyn]);
        printf("SOURCE: CORNER_MOM\n");
        printf("SINKS: KAON_000 KAON_001  KAON_011\n");
        for(j=0;j<nt;j++){
	    printf("%d %e %e %e %e %e %e\n",j,
	      props[nmasses*prop_kaon_000+i][j].real, props[nmasses*prop_kaon_000+i][j].imag,
	      props[nmasses*prop_kaon_001+i][j].real, props[nmasses*prop_kaon_001+i][j].imag,
	      props[nmasses*prop_kaon_011+i][j].real, props[nmasses*prop_kaon_011+i][j].imag );
        }
        printf("ENDPROP\n");
        printf("STARTPROP\n");
        printf("MASSES:  %.5e   %.5e\n",masses[i],masses[i]);
        printf("SOURCE: CORNER_MOM\n");
        printf("SINKS: RHO_000 RHO_001PERP  RHO_001PARA  RHO_011PERP RHO_011PARA\n");
        for(j=0;j<nt;j++){
	    printf("%d %e %e %e %e %e %e %e %e %e %e\n",j,
	      props[nmasses*prop_rho_000+i][j].real, props[nmasses*prop_rho_000+i][j].imag,
	      props[nmasses*prop_rho_001perp+i][j].real, props[nmasses*prop_rho_001perp+i][j].imag,
	      props[nmasses*prop_rho_001para+i][j].real, props[nmasses*prop_rho_001para+i][j].imag,
	      props[nmasses*prop_rho_011perp+i][j].real, props[nmasses*prop_rho_011perp+i][j].imag,
	      props[nmasses*prop_rho_011para+i][j].real, props[nmasses*prop_rho_011para+i][j].imag );
        }
        printf("ENDPROP\n");
    } /* i (masses) */

    fflush(stdout);
  } /* end if(this_node==0) */

  /* free arrays */
  free(props[0]); free(props);
  for(i=0; i<nmasses; i++){
      free(quark_props[i]);
  }
  free(quark_props);
  
  return(cgn);
} /* spectrum_mom */



/* "Multiply by" the quark-antiquark local pion operator, mom=px,py,pz */
void mult_pion_mom_temp( int fb, int px, int py, int pz,
    su3_vector *src, field_offset dest ){
   /* operator is gamma_pdir, another (-1)^(x+y+z+t) for antiquark */
    register int i;
    register site *s;
    complex phase;
    FORALLSITES(i,s){
	phase = ce_itheta( ((double)(px*s->x)/(double)(nx) +
	    (double)(py*s->y)/(double)(ny) +
	    (double)(pz*s->z)/(double)(nz))*2.0*PI );
	if( fb==BACKWARDS )phase.imag = -phase.imag;
	c_scalar_mult_su3vec( &(src[i]), &phase,
		(su3_vector *)F_PT(s,dest) );
    }
}


/* "Multiply by" the quark-antiquark local rho operator, mom=px,py,pz */
void mult_rho_mom_temp( int fb, int pdir, int px, int py, int pz,
    su3_vector *src, field_offset dest ){
   /* operator is gamma_pdir, another (-1)^(x+y+z+t) for antiquark */
    register int i;
    register site *s;
    complex phase;
    FORALLSITES(i,s){
	phase = ce_itheta( ((double)(px*s->x)/(double)(nx) +
	    (double)(py*s->y)/(double)(ny) +
	    (double)(pz*s->z)/(double)(nz))*2.0*PI );
	if( fb==BACKWARDS )phase.imag = -phase.imag;
	if( ((((short *)&(s->x))[pdir]) & 0x1) != 0 ){
	     CMULREAL( phase, -1.0, phase );
	}
	c_scalar_mult_su3vec( &(src[i]), &phase,
		(su3_vector *)F_PT(s,dest) );
    }
}
