/******** spectrum_mom.c *************/
/* MIMD version 7*/
/* DT 6/97
* start from ../ks_hybrids2/spectrum_hybrids4.c
* Need hadrons with various momenta, flexibility to change what we look at.
* DT 7/97 change -1^t in pi_2 operator to agree with old spectrum program.
* DT 9/97 added some more different momenta for quenched Naik project
* DT 3/99 masses are arguments.
*/
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.
*/
#define mat_invert mat_invert_uml
/**#define mat_invert mat_invert_cg**/

/**#define PION_000_SOURCE**/
#define PION_001_SOURCE
#define PION_011_SOURCE
/**#define PION_111_SOURCE**/
/**#define PION_002_SOURCE**/
/**#define PION_022_SOURCE**/
/**#define PION_222_SOURCE**/
/**#define PION2_000_SOURCE**/
/**#define PION2_001_SOURCE**/
/**#define RHO_000_SOURCE**/
#define RHO_001_SOURCE
/**#define RHO_010_SOURCE**/
/**#define RHO_011_SOURCE**/
/**#define RHO_110_SOURCE**/
/**#define RHO2_000_SOURCE**/
/**#define RHO2_001_SOURCE**/

/* Spectrum for Kogut-Susskind hadrons.
   Any antiperiodic boundary conditions (time direction usually)
   are absorbed in the link matrices (ie rephase(ON) )

    Gauge fixing should be done before calling this function
    (Coulomb gauge, probably).
*/


/* Symbolic names for propagators.  prop_SOURCE_SINK */
/* operators:
   pion:	qbar gamma_5 x gamma_5 q:  sign = +1
   pion2:	qbar gamma_0 x gamma_5 q> sign =  (-1)^(x+y+z+t)
   rho:		1--: qbar gamma_i gamma_0 q
		(1)^(dir) (VT)
   rho2:	1--: qbar gamma_i q
		(1)^(x+y+z+t+dir) (PV)
*/
enum prop_name { 
    prop_pi_000,	/* pion at momentum 0,0,0 */
    prop_pi_001,	/* pion at momentum 0,0,1 */
    prop_pi_011,	/* pion at momentum 0,1,0 */
    prop_pi_111,
    prop_pi_002,
    prop_pi_022,
    prop_pi_222,
    prop_pi2_000,
    prop_pi2_001,

    prop_rho_000,
    prop_rho_001,
    prop_rho_010,
    prop_rho_011,
    prop_rho_110,
    prop_rho2_000,
    prop_rho2_001,

    nprops		/* nprops = number of propagators */
};

#include "generic_ks_includes.h"

/* Various meson operators.  Source code later in this file */
/* All of these operators contain an extra gamma_5, entering when
    antiquark propagators are computed from the same source as quark
    propagators */
   /* quark-antiquark operators*/
void mult_pi_mom( int fb, int px, int py, int pz,
    field_offset src, field_offset dest );
void mult_pi2_mom( int fb, int px, int py, int pz,
    field_offset src, field_offset dest);
void mult_rho_mom( int fb, int pdir, int px, int py, int pz,
     field_offset src, field_offset dest );
void mult_rho2_mom( int fb, int pdir, int px, int py, int pz,
     field_offset src, field_offset dest );
int test_converge(int t_source);

int spectrum_mom( Real qmass, Real amass, field_offset temp, Real tol,
		  ferm_links_t *fn)
{
  /* arguments are quark and antiquark masses, return C.G. iteration number */

  int cgn;
  register int i,j;
  register site* s;
  register complex cc;
  register int t_source;
  int color;	/* color for source */
  int src_count; /* number of source time slices used */
  complex **props;	/* arrays of propagators */

  cgn=0; /* number of CG iterations */

  /* allocate arrays to accumulate propagators */
  props = (complex **)malloc(nprops*sizeof(complex *));
  props[0] = (complex *)malloc(nprops*nt*sizeof(complex));
  for(i=1;i<nprops;i++)props[i]=props[i-1]+nt;

  /* set propagators to zero */
  for(cc.real=cc.imag=0.0,i=0;i<nprops;i++)for(j=0;j<nt;j++){
    props[i][j]=cc;
  }

  /* loop over "source" time slice */
  for(src_count=0,t_source=source_start; t_source<nt && src_count<n_sources;
    t_source += source_inc,src_count++){
    /* Wall source */
    /* Use quark_source for quark source */
    if(this_node==0)printf("spectrum_mom(): source time = %d\n",t_source);

    for(color=0;color<3;color++){
	FORALLSITES(i,s){
	    clearvec( &(s->quark_source) );
/**if( (s->x%2)==0 && (s->y%2)==0 &&  (s->z%2)==0 )**/
	    if(s->t==t_source) s->quark_source.c[color].real=1.0;
	}

	/* compute M^-1 * quark_source */
	cgn += mat_invert( F_OFFSET(quark_source), F_OFFSET(quark_prop), 
			   temp, qmass, PRECISION, fn );
	/*if(t_source==0)test_converge(t_source);*/ /*TEMP*/
	/* TEMP: test inversion, */
	check_invert( F_OFFSET(quark_prop), F_OFFSET(quark_source), qmass,
		      tol, fn);

	/* Begin with the pion source and sink operator */
	/* Use mat_invert to construct anti_prop */
	/* ACTUALLY, nothing to do for pion - anti_prop = quark_prop */
	/* tie propagators together at sink end to project out desired mesons */
	/* add into propagator arrays at distance from source */

#ifdef PION_000_SOURCE
	/* 0-+ (pion) PROPAGATORS ***************************************/
	mult_pi_mom( FORWARDS,0,0,0, F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_pi_mom( FORWARDS,0,0,0, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pi_000][(s->t+nt-t_source)%nt], cc );
	}
#endif /*PION_000_SOURCE*/

#ifdef PION_001_SOURCE
	/* now the pion source at momentum 001 */
	mult_pi_mom( FORWARDS,0,0,1, F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_pi_mom( FORWARDS,0,0,1, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pi_001][(s->t+nt-t_source)%nt], cc );
	}
#endif /*PION_001_SOURCE*/

#ifdef PION_011_SOURCE
	/* now the pion source at momentum 011 */
	mult_pi_mom( FORWARDS,0,1,1, F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_pi_mom( FORWARDS,0,1,1, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pi_011][(s->t+nt-t_source)%nt], cc );
	}
#endif /*PION_011_SOURCE*/

#ifdef PION_111_SOURCE
	/* now the pion source at momentum 111 */
	mult_pi_mom( FORWARDS,1,1,1, F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_pi_mom( FORWARDS,1,1,1, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pi_111][(s->t+nt-t_source)%nt], cc );
	}
#endif /*PION_111_SOURCE*/

#ifdef PION_002_SOURCE
	/* now the pion source at momentum 002 */
	mult_pi_mom( FORWARDS,0,0,2, F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_pi_mom( FORWARDS,0,0,2, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pi_002][(s->t+nt-t_source)%nt], cc );
	}
#endif /*PION_002_SOURCE*/

#ifdef PION_022_SOURCE
	/* now the pion source at momentum 022 */
	mult_pi_mom( FORWARDS,0,2,2, F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_pi_mom( FORWARDS,0,2,2, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pi_022][(s->t+nt-t_source)%nt], cc );
	}
#endif /*PION_022_SOURCE*/

#ifdef PION_222_SOURCE
	/* now the pion source at momentum 222 */
	mult_pi_mom( FORWARDS,2,2,2, F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_pi_mom( FORWARDS,2,2,2, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pi_222][(s->t+nt-t_source)%nt], cc );
	}
#endif /*PION_222_SOURCE*/

#ifdef PION2_000_SOURCE
	/* now the pion2 source */
	mult_pi2_mom( FORWARDS,0,0,0,F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_pi2_mom( FORWARDS,0,0,0, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pi2_000][(s->t+nt-t_source)%nt], cc );
	}
#endif /*PION2_000_SOURCE*/

#ifdef PION2_001_SOURCE
	/* now the pion2 source at momentum 001 */
	mult_pi2_mom( FORWARDS,0,0,1,F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_pi2_mom( FORWARDS,0,0,1, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pi2_001][(s->t+nt-t_source)%nt], cc );
	}
#endif /*PION2_001_SOURCE*/


	/* 1-- (rho) PROPAGATORS ***************************************/
#ifdef RHO_000_SOURCE
	/* Now the rho source */
	mult_rho_mom( FORWARDS,ZUP,0,0,0,
	    F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_rho_mom( FORWARDS,ZUP,0,0,0,
	    F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_rho_000][(s->t+nt-t_source)%nt], cc );
	}
#endif /*RHO_000_SOURCE*/

#ifdef RHO_001_SOURCE
	mult_rho_mom( FORWARDS, ZUP, 0,0,1,
	    F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_rho_mom( FORWARDS, ZUP, 0,0,1,
	    F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_rho_001][(s->t+nt-t_source)%nt], cc );
	}
#endif /*RHO_001_SOURCE*/

#ifdef RHO_010_SOURCE
	mult_rho_mom( FORWARDS, ZUP, 0,1,0,
	    F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_rho_mom( FORWARDS, ZUP, 0,1,0,
	     F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_rho_010][(s->t+nt-t_source)%nt], cc );
	}
#endif /*RHO_010_SOURCE*/

#ifdef RHO_011_SOURCE
	mult_rho_mom( FORWARDS, ZUP, 0,1,1,
	     F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_rho_mom( FORWARDS, ZUP, 0,1,1,
	     F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_rho_011][(s->t+nt-t_source)%nt], cc );
	}
#endif /*RHO_011_SOURCE*/

#ifdef RHO_110_SOURCE
	mult_rho_mom( FORWARDS, ZUP, 1,1,0,
	    F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_rho_mom( FORWARDS, ZUP, 1,1,0,
	    F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_rho_110][(s->t+nt-t_source)%nt], cc );
	}
#endif /*RHO_110_SOURCE*/

#ifdef RHO2_000_SOURCE
	/* Now the rho2 source */
	mult_rho2_mom( FORWARDS, ZUP, 0,0,0,
	    F_OFFSET(quark_source), F_OFFSET(g_rand) );
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_rho2_mom( FORWARDS, ZUP, 0,0,0,
	    F_OFFSET(quark_prop), F_OFFSET(g_rand) );
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_rho2_000][(s->t+nt-t_source)%nt], cc );
	}
#endif /*RHO2_000_SOURCE*/

#ifdef RHO2_001_SOURCE
	mult_rho2_mom( FORWARDS, ZUP, 0,0,1,
	    F_OFFSET(quark_source), F_OFFSET(g_rand));
	cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );
	mult_rho2_mom( FORWARDS, ZUP, 0,0,1,
	    F_OFFSET(quark_prop), F_OFFSET(g_rand));
	FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_rho2_001][(s->t+nt-t_source)%nt], cc );
	}
#endif /*RHO2_000_SOURCE*/

    }
  } /* end loop on t_source */


  /* Sum propagator arrays over nodes */
  /* print out propagators */
  g_veccomplexsum( props[0] , nprops*nt );
  for(i=0;i<nprops;i++)for(j=0;j<nt;j++){
    CDIVREAL(props[i][j],n_sources*nx*ny*nz,props[i][j]);
  }
  if(this_node==0){

#ifdef PION_000_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: PION_000\n");
    printf("SINKS: PION_000\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	props[prop_pi_000][j].real, props[prop_pi_000][j].imag );
    }
    printf("ENDPROP\n");
#endif /*PION_000_SOURCE*/

#ifdef PION_001_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: PION_001\n");
    printf("SINKS: PION_001\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_pi_001][j].real, props[prop_pi_001][j].imag );
    }
    printf("ENDPROP\n");
#endif /*PION_001_SOURCE*/

#ifdef PION_011_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: PION_011\n");
    printf("SINKS: PION_011\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_pi_011][j].real, props[prop_pi_011][j].imag );
    }
    printf("ENDPROP\n");
#endif /*PION_011_SOURCE*/

#ifdef PION_111_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: PION_111\n");
    printf("SINKS: PION_111\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_pi_111][j].real, props[prop_pi_111][j].imag );
    }
    printf("ENDPROP\n");
#endif /*PION_111_SOURCE*/

#ifdef PION_002_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: PION_002\n");
    printf("SINKS: PION_002\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_pi_002][j].real, props[prop_pi_002][j].imag );
    }
    printf("ENDPROP\n");
#endif /*PION_002_SOURCE*/

#ifdef PION_022_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: PION_022\n");
    printf("SINKS: PION_022\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_pi_022][j].real, props[prop_pi_022][j].imag );
    }
    printf("ENDPROP\n");
#endif /*PION_022_SOURCE*/

#ifdef PION_222_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: PION_222\n");
    printf("SINKS: PION_222\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_pi_222][j].real, props[prop_pi_222][j].imag );
    }
    printf("ENDPROP\n");
#endif /*PION_222_SOURCE*/

#ifdef PION2_000_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: PION2_000\n");
    printf("SINKS: PION2_000\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_pi2_000][j].real, props[prop_pi2_000][j].imag );
    }
    printf("ENDPROP\n");
#endif /*PION2_000_SOURCE*/

#ifdef PION2_001_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: PION2_001\n");
    printf("SINKS: PION2_001\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_pi2_001][j].real, props[prop_pi2_001][j].imag );
    }
    printf("ENDPROP\n");
#endif /*PION2_001_SOURCE*/

#ifdef RHO_000_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: RHO_000\n");
    printf("SINKS: RHO_000\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_rho_000][j].real, props[prop_rho_000][j].imag );
    }
    printf("ENDPROP\n");
#endif /*RHO_000_SOURCE*/

#ifdef RHO_001_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: RHO_001\n");
    printf("SINKS: RHO_001\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_rho_001][j].real, props[prop_rho_001][j].imag );
    }
    printf("ENDPROP\n");
#endif /*RHO_001_SOURCE*/

#ifdef RHO_010_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: RHO_010\n");
    printf("SINKS: RHO_010\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_rho_010][j].real, props[prop_rho_010][j].imag );
    }
    printf("ENDPROP\n");
#endif /*RHO_010_SOURCE*/

#ifdef RHO_011_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: RHO_011\n");
    printf("SINKS: RHO_011\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_rho_011][j].real, props[prop_rho_011][j].imag );
    }
    printf("ENDPROP\n");
#endif /*RHO_011_SOURCE*/

#ifdef RHO_110_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: RHO_110\n");
    printf("SINKS: RHO_110\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_rho_110][j].real, props[prop_rho_110][j].imag );
    }
    printf("ENDPROP\n");
#endif /*RHO_110_SOURCE*/

#ifdef RHO2_000_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: RHO2_000\n");
    printf("SINKS: RHO2_000\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_rho2_000][j].real, props[prop_rho2_000][j].imag );
    }
    printf("ENDPROP\n");
#endif /*RHO2_000_SOURCE*/

#ifdef RHO2_001_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: RHO2_001\n");
    printf("SINKS: RHO2_001\n");
    for(j=0;j<nt;j++){
	printf("%d %e %e\n",j,
	  props[prop_rho2_001][j].real, props[prop_rho2_001][j].imag );
    }
    printf("ENDPROP\n");
#endif /*RHO2_001_SOURCE*/

    fflush(stdout);
  } /* end if(this_node==0) */

  /* free arrays */
  free(props[0]); free(props);
  
  return(cgn);
} /* spectrum_mom */


/* "Multiply by" the quark-antiquark local pion operator at momentum mom */
/* argument "fb" = "forwards/backwards" decides exp( +- k.x ) */
void mult_pi_mom( int fb, int px, int py, int pz,
    field_offset src, field_offset dest ){
   /* operator is (-1)^(x+y+z+t), another (-1)^(x+y+z+t) for antiquark,
	so just multiply by phase */
    register int i;
    register site *s;
    complex phase;
    FORALLSITES(i,s){
	phase = ce_itheta( ((double)(px*s->x)/(double)(nx) +
	    (double)(py*s->y)/(double)(ny) +
	    (double)(pz*s->z)/(double)(nz))*2.0*PI );
	if( fb==BACKWARDS )phase.imag = -phase.imag;
	c_scalar_mult_su3vec( (su3_vector *)F_PT(s,src), &phase,
	    (su3_vector *)F_PT(s,dest) );
    }
}

/* "Multiply by" the second quark-antiquark local pion operator, mom=px,py,pz */
void mult_pi2_mom( int fb, int px, int py, int pz, 
    field_offset src, field_offset dest ){
   /* operator is gamma_0 (1), another (-1)^(x+y+z+t) for antiquark, 
	so just multiply by (-1)^(x+y+z+t) */
    register int i;
    register site *s;
    complex phase;
    FORALLSITES(i,s){
	phase = ce_itheta( ((double)(px*s->x)/(double)(nx) +
	    (double)(py*s->y)/(double)(ny) +
	    (double)(pz*s->z)/(double)(nz))*2.0*PI );
	if( fb==BACKWARDS )phase.imag = -phase.imag;
	if( ((s->x+s->y+s->z)&0x1) != 0 ) CMULREAL( phase, -1.0, phase );
	c_scalar_mult_su3vec( (su3_vector *)F_PT(s,src), &phase,
	    (su3_vector *)F_PT(s,dest) );
    }
}

/* "Multiply by" the quark-antiquark local rho operator, mom=px,py,pz */
void mult_rho_mom( int fb, int pdir, int px, int py, int pz,
    field_offset src, field_offset dest ){
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
	c_scalar_mult_su3vec( (su3_vector *)F_PT(s,src), &phase,
		(su3_vector *)F_PT(s,dest) );
    }
}

/* "Multiply by" the second quark-antiquark local rho operator, mom=px,py,pz */
void mult_rho2_mom( int fb, int pdir, int px, int py, int pz,
   field_offset src, field_offset dest ){
   /* operator is gamma_0 gamma_pdir, another (-1)^(x+y+z+t) for antiquark */
    register int i;
    register site *s;
    complex phase;
    FORALLSITES(i,s){
	phase = ce_itheta( ((double)(px*s->x)/(double)(nx) +
	    (double)(py*s->y)/(double)(ny) +
	    (double)(pz*s->z)/(double)(nz))*2.0*PI );
	if( fb==BACKWARDS )phase.imag = -phase.imag;
	if( ((((short *)&(s->x))[pdir] + s->x+s->y+s->z) & 0x1) != 0 ){
	     CMULREAL( phase, -1.0, phase );
	}
	c_scalar_mult_su3vec( (su3_vector *)F_PT(s,src), &phase,
		(su3_vector *)F_PT(s,dest) );
    }
}


