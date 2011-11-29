/******** spectrum_hybrids5.c *************/
/* MIMD version 7*/
/* DT 11/21/95 started */
/* DT 2/21/96 version2.  Make it easy to add operators.  Keep collection
   of "mult_by_1mps" routines around.  Select only operators that look
   interesting in previous tests. 
   DT 7/10/96 multiple sink operators for each source.
   New output format:  "SINK:" line lists sink operators, 
    "operator_SRC:" lines give propagators from each source.
   DT 10/4/96 add 0+- "baryon number" operator
   DT 1/21/97 use defines to select desired source operators, add
	1-+ 4 quark operator QQQQ
   DT 4/4/97  convert to KS quarks
   CD 11/09/97 added parentheses to A & B == C -> changed results
   DT 9/01 fixed many bugs
   DT 1/02 need flavor structures that are singlets under spatial rotations
   DT 1/24/02 symmetrize - average over field at quark and antiquark
   CD 8/02, moved all nonhybrid q qbar code to spectrum_nlpi2.c
*/
#define mat_invert mat_invert_uml

/**#define ONEMP_0_SOURCE**/
#define ONEMP_S_SOURCE

/* Spectrum for Kogut-Susskind hybrid mesons. 
   Any antiperiodic boundary conditions (time direction usually)
   are absorbed in the link matrices (ie phases_in=ON )

   Also does spectrum for a few conventional mesons, as reference.

    Gauge fixing should be done before calling this function
    (Coulomb gauge, probably).
*/


/* Symbolic names for propagators.  prop_SOURCE_SINK */
/* operators:
   pion:	conventional 0-+:  qbar gamma_5 q
		(1)
   pion2:	conventional 0-+:  qbar gamma_0 gamma_5 q
		(-1)^(x+y+z+t)
   0mp:		0-+ hybrid
   rho:		conventional 1--: qbar gamma_i gamma_0 q
		(1)^(dir) (VT)
   rho2:	conventional 1--: qbar gamma_i q
		(1)^(x+y+z+t+dir) (PV)
   1mm:		1-- hybrid
   a1:		conventional 1++: qbar gamma_5 gamma_i q
   a1P:		P wave a_1: qbar eps_ijk gamma_j partial_k q
   1pp:		1++ hybrid
   0pm:		0+- hybrid
   0pmP		0+- P-wave hybrid
   0pmB		0+- bilinear, qbar gamma_0 q
   0mm		0-- hybrid
   0mmP:	0-- P-wave hybrid
   1mp0:	1-+ hybrid, qbar gamma_i q eps_ijk B_j, flavor gamma_0
   1mps:	1-+ hybrid, qbar gamma_i q eps_ijk B_j, flavor 1
   1mps:	1-+ hybrid, qbar gamma_0 q E_i
   qqqq	1-+, a_1 plus pion
*/
enum prop_name { 

    prop_1mp0_1mp0,
    prop_1mps_1mps,
    prop_qqqq_1mp0,
    prop_qqqq_1mps,

    prop_2pm_2pm,

    nprops		/* nprops = number of propagators */
};

#include "generic_ks_includes.h"
#include "../include/field_strength.h"

/* Various meson operators.  Source code later in this file */
/* All of these operators contain an extra gamma_5, entering when
    antiquark propagators are computed from the same source as quark
    propagators */

  /* exotic hybrids, and one exotic bilinear */
void mult_zero_pm( field_offset src, field_offset dest );
void mult_zero_mm( field_offset src, field_offset dest );
void mult_1mp0( int pdir, field_offset src, field_offset dest );
void mult_1mps( int pdir, field_offset src, field_offset dest );
void mult_1mps_rev( int pdir, field_offset src, field_offset dest );
void mult_1mpE( int pdir, field_offset src, field_offset dest );
void mult_qqqq1mp0( int src_t, int pdir, field_offset src, field_offset dest,
		    field_offset work, field_offset temp, Real mass,
		    ferm_links_t *fn);
   /* non-exotic hybrids*/
void mult_zero_mp( field_offset src, field_offset dest );
void mult_one_mm( int pdir, field_offset src, field_offset dest );
void mult_one_pp( int pdir, field_offset src, field_offset dest );

void mult_by_field_strength( int dir1, int dir2,
    field_offset src, field_offset dest );
int test_converge(int t_source);

int spectrum_hybrids( Real mass, field_offset temp, Real tol, ferm_links_t *fn )
{ /* return the C.G. iteration number */

  int cgn;
  register int i,j;
  register site* s;
  register complex cc;
  register int t_source;
  int dir;	/* direction in lattice */
  int color;	/* color for source */
  int src_count; /* number of source time slices used */
  complex **props;	/* arrays of propagators */

  node0_printf("Spectrum_hybrids: mass = %e\n",mass);
  cgn=0; /* number of CG iterations */

  /* allocate arrays to accumulate propagators */
  props = (complex **)malloc(nprops*sizeof(complex *));
  props[0] = (complex *)malloc(nprops*nt*sizeof(complex));
  for(i=1;i<nprops;i++)props[i]=props[i-1]+nt;

  /* compute the field strength tensor, uses smeared links */
  /* links may already be smeared once for "fat link" inversions. */
  /* uses field_strength[] for temporary space */

#define SMEAR 32
#define staple_weight 0.25

/* Temporary - running in old version6 */
#define APE_PROJECT

#ifdef APE_PROJECT
#define NHIT 10
#else
#define NHIT 0
#endif

  rephase(OFF);
  if(this_node==0)printf("SMEARING IS ON, level %d\n",SMEAR);
  if(SMEAR>=1){
    node0_printf( "Smearing: staple_weight = %e\n", staple_weight);
#ifdef APE_PROJECT
    node0_printf( "APE projection is on, NHIT = %d\n",NHIT);
#endif
    ape_smear( F_OFFSET(link[0]), F_OFFSET(tmplink[0]), staple_weight, u0, 1, 3*NHIT, 0.);
  }
  else FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
    s->tmplink[dir] = s->link[dir];
  }
  for(j=2;j<=SMEAR;j++){
    FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
      s->field_strength[dir] = s->tmplink[dir];
    }
    /**    node0_printf( "Smearing: staple_weight = %e\n", staple_weight); **/
#ifdef APE_PROJECT
    /**    node0_printf( "APE projection is on, NHIT = %d\n",NHIT); **/
#endif
    ape_smear( F_OFFSET(field_strength[0]), F_OFFSET(tmplink[0]), staple_weight, u0, 1, 3*NHIT, 0.);
  }
  make_field_strength(F_OFFSET(tmplink[0]), F_OFFSET(field_strength[0]));
  rephase(ON);
  /**temp_test();**/

  /* set propagators to zero */
  for(cc.real=cc.imag=0.0,i=0;i<nprops;i++)for(j=0;j<nt;j++){
    props[i][j]=cc;
  }

  /* loop over "source" time slice */
  for(src_count=0,t_source=source_start; t_source<nt && src_count<n_sources;
    t_source += source_inc,src_count++){
    /* Wall source */
    /* Use quark_source for quark source */
    if(this_node==0)printf("spectrum_hybrids(): source time = %d\n",t_source);

    for(color=0;color<3;color++){
        FORALLSITES(i,s){
	    clearvec( &(s->quark_source) );
            if(s->t==t_source) s->quark_source.c[color].real=1.0;
	}

        /* compute M^-1 * quark_source */
        cgn += mat_invert( F_OFFSET(quark_source), F_OFFSET(quark_prop),
			   temp, mass, PRECISION, fn );
	/*if(t_source==0)test_converge(t_source);*/ /*TEMP*/
	/* TEMP: test inversion, */
	check_invert( F_OFFSET(quark_prop), F_OFFSET(quark_source), mass, tol,
		      fn);

        /* Begin with the pion source and sink operator */
        /* Use mat_invert to construct anti_prop */
	/* ACTUALLY, nothing to do for pion - anti_prop = quark_prop */
        /* tie propagators together at sink end to project out desired mesons */
        /* add into propagator arrays at distance from source */


#ifdef ZEROPM_SOURCE
	/* 0+- (exotic) PROPAGATORS ***************************************/
        /* Now the 0+-_0+- propagator */
        mult_zero_pm( F_OFFSET(quark_source), F_OFFSET(g_rand) );
        cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop),
			   temp, mass, PRECISION, fn );
        mult_zero_pm( F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_0pm_0pm][(s->t+nt-t_source)%nt], cc );
        }
#endif /*ZEROPM_SOURCE*/

#ifdef ONEMP_0_SOURCE
	/* 1-+ (exotic) PROPAGATORS ***************************************/
        /* Now do the 1-+_1-+ source.  For the moment, Z component only */
        mult_1mp0( ZUP, F_OFFSET(quark_source), F_OFFSET(g_rand) );
        cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop),
			   temp, mass, PRECISION, fn );
        mult_1mp0( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_1mp0_1mp0][(s->t+nt-t_source)%nt], cc );
        }
#endif /*ONEMP_0_SOURCE*/

#ifdef ONEMP_S_SOURCE
	/* 1-+ (exotic) PROPAGATORS ***************************************/
        /* Now do the 1-+_1-+ source.  For the moment, Z component only */
        mult_1mps( ZUP, F_OFFSET(quark_source), F_OFFSET(g_rand) );
        cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop),
			   temp, mass, PRECISION, fn );
        mult_1mps( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
/**mult_1mps_rev( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );**/
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_1mps_1mps][(s->t+nt-t_source)%nt], cc );
        }
#endif /*ONEMP_S_SOURCE*/

#ifdef ONEMP_E_SOURCE
        /* Now the 1-+2_1+-2 source.  For the moment, Z component only */
        mult_1mpE( ZUP, F_OFFSET(quark_source), F_OFFSET(g_rand) );
        cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop),
			   temp, mass, PRECISION, fn );
        mult_1mpE( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_1mps_1mps][(s->t+nt-t_source)%nt], cc );
        }
        mult_1mp0( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_1mps_1mp0][(s->t+nt-t_source)%nt], cc );
        }
#endif /*ONEMP_E_SOURCE*/

#ifdef QQQQ_SOURCE
        /* Now the 4-quark (a1-pi) source.  For the moment, Z component only */
        mult_qqqq1mp0( t_source,  ZUP,
		       F_OFFSET(quark_source), F_OFFSET(g_rand), 
		       F_OFFSET(anti_prop), F_OFFSET(phi1), mass1, fn );
        cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop),
			   temp, mass, PRECISION, fn );

        mult_1mpE( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_qqqq_1mps][(s->t+nt-t_source)%nt], cc );
        }
        mult_1mp0( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_qqqq_1mp0][(s->t+nt-t_source)%nt], cc );
        }
#endif /*QQQQ_SOURCE*/

	
	/* 2+- (exotic) PROPAGATORS ***************************************/
    }
  } /* end loop on t_source */

  /* Sum propagator arrays over nodes */
  /* print out propagators */
  g_veccomplexsum( props[0] , nprops*nt );
  for(i=0;i<nprops;i++)for(j=0;j<nt;j++){
    CDIVREAL(props[i][j],nx*ny*nz*n_sources,props[i][j]);
  }
  if(this_node==0){

#ifdef ZEROPM_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",mass,mass );
    printf("SOURCE: 0PM\n");
    printf("SINKS: 0PM 0PMP 0PMB\n");
    for(j=0;j<nt;j++){
        printf("%d %e %e %e %e %e %e\n",j,
          props[prop_0pm_0pm][j].real, props[prop_0pm_0pm][j].imag,
          props[prop_0pm_0pmP][j].real, props[prop_0pm_0pmP][j].imag,
          props[prop_0pm_0pmB][j].real, props[prop_0pm_0pmB][j].imag );
    }
    printf("ENDPROP\n");
#endif /*ZEROPM_SOURCE*/

#ifdef ONEMP_0_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",mass,mass );
    printf("SOURCE: ONEMP_0\n");
    printf("SINKS: ONEMP_0\n");
    for(j=0;j<nt;j++){
        printf("%d %e %e\n",j,
          props[prop_1mp0_1mp0][j].real, props[prop_1mp0_1mp0][j].imag );
    }
    printf("ENDPROP\n");
#endif /*ONEMP_0_SOURCE*/

#ifdef ONEMP_S_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",mass,mass );
    printf("SOURCE: ONEMP_S\n");
    printf("SINKS: ONEMP_S\n");
    for(j=0;j<nt;j++){
        printf("%d %e %e\n",j,
          props[prop_1mps_1mps][j].real, props[prop_1mps_1mps][j].imag );
    }
    printf("ENDPROP\n");
#endif /*ONEMP_S_SOURCE*/

#ifdef QQQQ_SOURCE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",mass,mass );
    printf("SOURCE: QQQQ\n");
    printf("SINKS: ONEMP_0 ONEMP_S\n");
    for(j=0;j<nt;j++){
        printf("%d %e %e %e %e\n",j,
          props[prop_qqqq_1mp0][j].real, props[prop_qqqq_1mp0][j].imag,
          props[prop_qqqq_1mps][j].real, props[prop_qqqq_1mps][j].imag );
    }
    printf("ENDPROP\n");

#ifdef NONSENSE
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",mass,mass );
    printf("SOURCE: QQQQ\n");
    printf("SINKS: A1 B1\n");
    for(j=0;j<nt;j++){
        printf("%d %e %e %e %e\n",j,
          props[prop_qqqq_junk][j].real, props[prop_qqqq_junk][j].imag,
          props[prop_qqqq_junk2][j].real, props[prop_qqqq_junk2][j].imag );
    }
    printf("ENDPROP\n");
#endif /*NONSENSE*/
#endif /*QQQQ_SOURCE*/

    fflush(stdout);
  } /* end if(this_node==0) */

  /* free arrays */
  free(props[0]); free(props);
  
  return(cgn);
} /* spectrum_hybrids */


/* "Multiply by" the "S-wave" zero-plus-minus operator 
	quark operator is a1 = (-1)^(x+y+z+t) gamma_i = 1++,
	gluon operator is B_i = 1+-,  take dot product */
void mult_zero_pm( field_offset src, field_offset dest ){
    register int dir,i;
    register site *s;

if(this_node==0)printf("OOPS, mult_zero_pm NOT WRITTEN\n");
exit(0);
    /* set destination to zero */
    FORALLSITES(i,s){
	clearvec( (su3_vector *)F_PT(s,dest) );
    }

    /* Loop over spatial directions */
    for(dir=XUP;dir<=ZUP;dir++){
        /* Multiply by dir+1,dir+2 component of magnetic field,
           multiply by gamma_dir (-1)^(x+y+z+t)(for a1) and gamma_five (for
           antiquark propagator ) */
        mult_by_field_strength( (dir+1)%3, (dir+2)%3, src, F_OFFSET(cg_p) );
        FORALLSITES(i,s){
        }
    } /* end loop on dir */
}


/* "Multiply by" the zero-minus-minus operator 
	quark operator is a1 = (-1)^(x+y+z+t) gamma_i = 1++,
	gluon operator is E_i = 1--,  take dot product */
void mult_zero_mm( field_offset src, field_offset dest ){
    register int dir,i;
    register site *s;

if(this_node==0)printf("OOPS, mult_zero_mm NOT WRITTEN\n");
exit(0);
    /* set destination to zero */
    FORALLSITES(i,s){
	clearvec( (su3_vector *)F_PT(s,dest) );
    }

    /* Loop over spatial directions */
    for(dir=XUP;dir<=ZUP;dir++){
        /* Multiply by dir component of electric field,
           multiply by gamma_dir (-1)^(x+y+z+t)(for a1) and gamma_five (for
           antiquark propagator ) */
        mult_by_field_strength( dir, TUP, src, F_OFFSET(cg_p) );
        FORALLSITES(i,s){
        }
    } /* end loop on dir */
}



/* "Multiply by" a one-minus-plus operator 
    quark operator is "rho0", gluon operator is magnetic field B.
    Cross product is spin one.
    Z component = "rho_X*F_XZ + rho_Y*F_YZ.
   "pdir" is the polarization direction of the meson */
void mult_1mp0( int pdir, field_offset src, field_offset dest ){
    /* use cg_p and ttt as temporary storage */
    register int dir,i;
    register site *s;

    /* set destination to zero */
    FORALLSITES(i,s){
	clearvec( (su3_vector *)F_PT(s,dest) );
    }

    /* Loop over spatial directions orthogonal to pdir */
    for(dir=XUP;dir<=ZUP;dir++)if(dir != pdir){
	/* Multiply by dir,pdir component of magnetic field, 
	   multiply by flavor gamma_0 rho (includes gamma_5 for antiquark
	   propagator) */
	mult_by_field_strength( dir, pdir, src, F_OFFSET(cg_p) );
	mult_rho0( dir, F_OFFSET(cg_p),F_OFFSET(ttt) );
        FORALLSITES(i,s){
            add_su3_vector( (su3_vector *)F_PT(s,dest),
                (su3_vector *)&(s->ttt), (su3_vector *)F_PT(s,dest) );
        }
	mult_rho0( dir, src, F_OFFSET(cg_p) );
	mult_by_field_strength( dir, pdir, F_OFFSET(cg_p), F_OFFSET(ttt)  );
        FORALLSITES(i,s){
            add_su3_vector( (su3_vector *)F_PT(s,dest),
                (su3_vector *)&(s->ttt), (su3_vector *)F_PT(s,dest) );
        }
    } /* end loop on dir */
} /* end mult_1mp0 */

/* "Multiply by" a one-minus-plus operator 
    quark operator is "rhos", gluon operator is magnetic field B.
    Cross product is spin one.
    Z component = "rho_X*F_XZ + rho_Y*F_YZ.
   "pdir" is the polarization direction of the meson */
void mult_1mps( int pdir, field_offset src, field_offset dest ){
    /* use cg_p and ttt as temporary storage */
    register int dir,i;
    register site *s;

    /* set destination to zero */
    FORALLSITES(i,s){
	clearvec( (su3_vector *)F_PT(s,dest) );
    }

    /* Loop over spatial directions orthogonal to pdir */
    for(dir=XUP;dir<=ZUP;dir++)if(dir != pdir){
	/* Multiply by dir,pdir component of magnetic field, 
	   multiply by flavor singlet rho (includes gamma_5 for antiquark
	   propagator) */
	mult_by_field_strength( dir, pdir, src, F_OFFSET(cg_p) );
	mult_rhois( dir, F_OFFSET(cg_p),F_OFFSET(ttt) );
        FORALLSITES(i,s){
            add_su3_vector( (su3_vector *)F_PT(s,dest),
                (su3_vector *)&(s->ttt), (su3_vector *)F_PT(s,dest) );
        }
	mult_rhois( dir, src, F_OFFSET(cg_p) );
	mult_by_field_strength( dir, pdir, F_OFFSET(cg_p), F_OFFSET(ttt)  );
        FORALLSITES(i,s){
            add_su3_vector( (su3_vector *)F_PT(s,dest),
                (su3_vector *)&(s->ttt), (su3_vector *)F_PT(s,dest) );
        }
    } /* end loop on dir */
} /* end mult_1mps */
void mult_1mps_rev( int pdir, field_offset src, field_offset dest ){
    /* use cg_p and ttt as temporary storage */
    register int dir,i;
    register site *s;

    /* set destination to zero */
    FORALLSITES(i,s){
	clearvec( (su3_vector *)F_PT(s,dest) );
    }

    /* Loop over spatial directions orthogonal to pdir */
    for(dir=XUP;dir<=ZUP;dir++)if(dir != pdir){
	/* Multiply by dir,pdir component of magnetic field, 
	   multiply by flavor gamma_0 rho (includes gamma_5 for antiquark
	   propagator) */
	mult_rhois( dir, src, F_OFFSET(cg_p) );
	mult_by_field_strength( dir, pdir, F_OFFSET(cg_p), F_OFFSET(ttt)  );
        FORALLSITES(i,s){
            add_su3_vector( (su3_vector *)F_PT(s,dest),
                (su3_vector *)&(s->ttt), (su3_vector *)F_PT(s,dest) );
        }
    } /* end loop on dir */
} /* end mult_1mps_rev */

/* "Multiply by" another one-minus-plus operator 
    quark operator is zero-plus-minus bilinear (parity partner of pion!)
    gluon operator is electric field E.
    Z component = "psibar gamma_0 psi * F_0Z
   "pdir" is the polarization direction of the meson */
void mult_1mpE( int pdir, field_offset src, field_offset dest ){
    /* use cg_p as temporary storage */
    register int i;
    register site *s;

    /* Multiply by pdir component of electric field, 
    multiply by gamma_0 and gamma_five (for antiquark propagator ) */
    mult_by_field_strength( TUP, pdir, src, F_OFFSET(cg_p) );
    FORALLSITES(i,s){
	if( (s->t & 0x1) == 0 ){
	    *(su3_vector *)F_PT(s,dest) = *((su3_vector *)&(s->cg_p));
	} else {
	    scalar_mult_su3_vector( (su3_vector *)&(s->cg_p), -1.0,
		(su3_vector *)F_PT(s,dest) );
	}
    }
} /* end mult_1mpE */

/* "Multiply by" another one-minus-plus operator 
   a 4-quark operator, a_1 plus pion.  All quark lines different
   flavors, so don't worry about multiple contractions.
   "src_t" is the operator time slice
   "pdir" is the polarization direction of the meson */
void mult_qqqq1mp0( int src_t, int pdir, field_offset src, field_offset dest,
		    field_offset work, field_offset temp, Real mass,
		    ferm_links_t *fn){
    register int i;
    register site *s;
    /* use "work" and "temp" as temporary storage */
if(this_node==0)printf("OOPS, mult_qqqq1mp0 NOT WRITTEN\n");
exit(0);

    /* Multiply by a_1 operator, extra (-1)^(x+y+z+t) for direction of
	quark line */
    /* compute propagator, back to same time slice */
    /* multiply by pion operator */
    FORALLSITES(i,s){
	/* two (-1)^(x+y+z+t)'s = 1 */
    }
    mat_invert( dest, work, temp, mass, PRECISION, fn );
    FORALLSITES(i,s){
	if(s->t== src_t ){
	    /* really three (-1)^(x+y+z+t)'s - both propagators have one end here */
	}
	else{
	}
    }
} /* end mult_qqqq1mp0 */

/* "Multiply by" a one-minus-minus hybrid operator 
    quark operator is "pion", gluon operator is magnetic field.
   "pdir" is the polarization direction of the meson */
void mult_one_mm( int pdir, field_offset src, field_offset dest ){
    /* use cg_p as temporary storage */

    /* multiply by magnetic field.  (-1)^(x+y+z+t) for antiquark cancels
	(-1)^(x+y+z+t) in pion operator. */
    mult_by_field_strength( (pdir+1)%3, (pdir+2)%3, src, dest );
} /* end mult_one_mm */


/* "Multiply by" a one-plus-plus hybrid operator 
    quark operator is "rho", gluon operator is electric field.
    1++ = epsilon_ijk \psibar gamma_j \psi F_{0,k}
   "pdir" is the polarization direction of the meson */
void mult_one_pp( int pdir, field_offset src, field_offset dest ){
    /* use cg_p as temporary storage */
    register int i,j,k;
    register site *s;
    su3_vector tvec;
    register Real sign;

    /* set destination to zero */
    FORALLSITES(i,s){
	clearvec( (su3_vector *)F_PT(s,dest) );
    }
    /* loop over directions, sign is sign of epsilon tensor */
    for(j=XUP;j<=ZUP;j++)for(k=XUP;k<=ZUP;k++){
	if( j==k || j==pdir || k== pdir )continue;
	else if ( j == ((pdir+1)%3) )sign = 1.0;
	else sign = -1.0;

        mult_by_field_strength( TUP, k, src, F_OFFSET(cg_p) );
	FORALLSITES(i,s){
	    if( ( (((short *)&(s->x))[pdir]) & 0x1 ) == 0 ) sign *= 1.0;
	    else sign *= -1.0;
	    scalar_mult_add_su3_vector( (su3_vector *)F_PT(s,dest), &tvec, sign,
		(su3_vector *)F_PT(s,dest) );
	}

    }/* end loop over j,k (directions) */

} /* end mult_one_pp */

/* "Multiply by" a zero-minus-plus hybrid operator 
    quark operator is "rho", gluon operator is magnetic field.
    0-+ = epsilon_ijk \psibar gamma_i \psi F_{j,k} */
void mult_zero_mp( field_offset src, field_offset dest ){
    /* use cg_p as temporary storage */
    register int i,j,k,in;
    register site *s;
    register Real x;

    /* set destination to zero */
    FORALLSITES(i,s){
	clearvec( (su3_vector *)F_PT(s,dest) );
    }
    /* loop over directions, */
    for(i=XUP;i<=ZUP;i++){
	/* antisymmetry of epsilon and F_{jk} means no need to sum all terms */
	j=(i+1)%3; k = (i+2)%3;

        mult_by_field_strength( j, k, src, F_OFFSET(cg_p) );
	FORALLSITES(in,s){
	    if( ( (((short *)&(s->x))[i]) & 0x1 ) == 0 ) x = 1.0;
	    else x = -1.0;
	    scalar_mult_add_su3_vector( (su3_vector *)F_PT(s,dest), &(s->cg_p),
		x, (su3_vector *)F_PT(s,dest) );
	}

    }/* end loop over i,j,k (directions) */

} /* end mult_zero_mp */


/* Multiply by the field strength.  Arguments are two indices on
  field strength tensor, source and dest vectors */
void mult_by_field_strength( int dir1, int dir2,
    field_offset src, field_offset dest ){
    su3_vector tvec;

    register site *s;
    register int i;
    int fs_dir=-99;	/* index in field strength tensor, FS_XY to FS_ZT */
    Real sign;	/* minus one if indices in "wrong" order */

    /* figure out what fs_dir is */
    if(dir1==dir2){
	if(this_node==0)printf("Bad call to mult_by_field_strength\n");
	terminate(0);
    }
    if(dir1>dir2){ 
	sign= -1.0;
	i=dir1; dir1=dir2; dir2=i;	/* switch dirs */
    }
    else sign= 1.0;
    switch( dir1+10*dir2 ){
	case XUP+10*YUP: fs_dir = FS_XY; break;
	case XUP+10*ZUP: fs_dir = FS_XZ; break;
	case XUP+10*TUP: fs_dir = FS_XT; break;
	case YUP+10*ZUP: fs_dir = FS_YZ; break;
	case YUP+10*TUP: fs_dir = FS_YT; break;
	case ZUP+10*TUP: fs_dir = FS_ZT; break;
	default: 
            if(this_node==0)printf("Bad call to mult_by_field_strength\n");
            terminate(0); break;
    }

    /* multiply by matrix */
    FORALLSITES(i,s){
	mult_su3_mat_vec( &(s->field_strength[fs_dir]), 
	    (su3_vector *)F_PT(s,src), &tvec );
	scalar_mult_su3_vector( &tvec, sign, (su3_vector *)F_PT(s,dest) );
    }
}

/*TEMP TEST change gauge, see if things change appropriately*
temp_test(){
int dir,fs_dir;
complex cc;
site *s;
s = &(lattice[1]); 
   int nvector = 3;
   field_offset vector_offset[3] = { F_OFFSET(g_rand), F_OFFSET(phi), 
        F_OFFSET(xxx) };
   int vector_parity[3] = { EVENANDODD, EVEN, EVEN };
   int nantiherm = 4;
   field_offset antiherm_offset[4] = { F_OFFSET(mom[0]), F_OFFSET(mom[1]),
       F_OFFSET(mom[2]), F_OFFSET(mom[3]) };
   field_offset antiherm_parity[4] = { EVENANDODD, EVENANDODD, EVENANDODD,
       EVENANDODD }

    for(dir=XUP;dir<=TUP;dir++){
	cc = det_su3( &(s->tmplink[dir]) );
	printf("BEFORE: %d %d %d %d dir= %d  Det = %f  %f\n",
	s->x,s->y,s->z,s->t,dir,cc.real,cc.imag);
    }
    for(fs_dir=FS_XY;fs_dir<=FS_ZT;fs_dir++){
	cc = det_su3( &(s->field_strength[fs_dir]) );
	printf("BEFORE: %d %d %d %d dir= %d  Det = %f  %f\n",
	s->x,s->y,s->z,s->t,fs_dir,cc.real,cc.imag);
    }

   rephase( OFF );
   gaugefix_combo(XUP,(Real)1.8,500,(Real)GAUGE_FIX_TOL,
       nvector,vector_offset,vector_parity,
       nantiherm,antiherm_offset,antiherm_parity);
#ifdef APE_PROJECT
    node0_printf( "APE projection is on, NHIT = %d\n",NHIT);
#endif
    ape_smear(F_OFFSET(link[0]), F_OFFSET(tmplink[0]), staple_weight, u0, 1, 3*NHIT, 0.);
   make_field_strength(F_OFFSET(link[0]), F_OFFSET(field_strength[0]));
   rephase( ON );
    for(dir=XUP;dir<=TUP;dir++){
	cc = det_su3( &(s->tmplink[dir]) );
	printf("AFTER: %d %d %d %d dir= %d  Det = %f  %f\n",
	s->x,s->y,s->z,s->t,dir,cc.real,cc.imag);
    }
    for(fs_dir=FS_XY;fs_dir<=FS_ZT;fs_dir++){
	cc = det_su3( &(s->field_strength[fs_dir]) );
	printf("AFTER: %d %d %d %d dir= %d  Det = %f  %f\n",
	s->x,s->y,s->z,s->t,fs_dir,cc.real,cc.imag);
    }
}
*END TEMP TEST */

