/******** spectrum_nlpi2.c *************/
/* MIMD version 7*/
/* DT 2/25/98 started
   All KS pions, and some rhos
   DT 3/99, masses are arguments
   CD 8/02, collected all q qbar operators here and removed duplicates
*/
#define mat_invert mat_invert_uml

/* Spectrum for Kogut-Susskind hybrid mesons. 
   Any antiperiodic boundary conditions (time direction usually)
   are absorbed in the link matrices (ie phases_in=ON )

   Also does spectrum for a few conventional mesons, as reference.

    Gauge fixing should be done before calling this function
    (Coulomb gauge, probably).
*/


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
enum prop_name { 
    prop_pion5_pion5,
    prop_pion05_pion05,
    prop_pioni5_pioni5,
    prop_pionij_pionij,
    prop_pioni_pioni,
    prop_pioni0_pioni0,
    prop_pions_pions,
    prop_pion0_pion0,
    prop_rhoi_rhoi,
    prop_rhoi0_rhoi0,
    prop_rhos_rhos,
    prop_rho0_rho0,
    nprops		/* nprops = number of propagators */
};

#include "generic_ks_includes.h"

/* Various meson operators.  Source code later in this file */
/* All of these operators contain an extra gamma_5, entering when
    antiquark propagators are computed from the same source as quark
    propagators */
   /* quark-antiquark operators*/
void mult_pion5( field_offset src, field_offset dest ) ;
void mult_pion05( field_offset src, field_offset dest ) ;
void mult_pioni5( int flavor_dir, field_offset src, field_offset dest ) ;
void mult_pionij( int flavor_dir, field_offset src, field_offset dest ) ;
void mult_pioni( int fdir, field_offset src, field_offset dest ) ;
void mult_pioni0( int fdir, field_offset src, field_offset dest ) ;
void mult_pions(field_offset src, field_offset dest ) ;
void mult_pion0(field_offset src, field_offset dest ) ;
void mult_rhoi( int pdir, field_offset src, field_offset dest ) ;
void mult_rhoi0( int pdir, field_offset src, field_offset dest ) ;
void mult_rhos( int fdir,  field_offset src, field_offset dest ) ;
void mult_rho0( int fdir,  field_offset src, field_offset dest ) ;
void mult_a1( int pdir,  field_offset src, field_offset dest ) ;
void mult_b1( int pdir,  field_offset src, field_offset dest ) ;

int test_converge(int t_source);

int spectrum_nlpi2( Real qmass, Real amass, field_offset temp, 
		    Real tol, ferm_links_t *fn ){
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
    if(this_node==0)printf("spectrum_nlpi2(): source time = %d\n",t_source);


    for(color=0;color<3;color++){
      FORALLSITES(i,s){
	clearvec( &(s->quark_source) );
/**if(s->parity==EVEN)**/
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_SPECTRUM_DEBUG_POINT_SOURCE
    if(i==node_index(0,0,0,t_source)) { /* source at (0,0,0,0) */
#endif /* HISQ_SPECTRUM_DEBUG_POINT_SOURCE */
#endif /* MILC_GLOBAL_DEBUG */
	if(s->t==t_source) 
	  s->quark_source.c[color].real=1.0; /* Even Wall source */
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_SPECTRUM_DEBUG_POINT_SOURCE
    }
#endif /* HISQ_SPECTRUM_DEBUG_POINT_SOURCE */
#endif /* MILC_GLOBAL_DEBUG */
      }

        /* compute M^-1 * quark_source, M^-1 * antiquark_source */
        cgn += mat_invert( F_OFFSET(quark_source), F_OFFSET(quark_prop), 
			   temp, qmass, PRECISION, fn );
	/*if(t_source==0)test_converge(t_source);*/ /*TEMP*/
	/* TEMP: test inversion, */
	check_invert( F_OFFSET(quark_prop), F_OFFSET(quark_source), qmass,
		      tol, fn);

	/* make antiquark source by summing over desired operators
	   times quark source */
	/* First source couples to pion5, pioni5, pioni, pions, rhoi, rhos */
	mult_pion5( F_OFFSET(quark_source), F_OFFSET(g_rand) );
	mult_pioni( ZUP, F_OFFSET(quark_source), F_OFFSET(anti_prop) );
	FORALLSITES(i,s){ add_su3_vector( &(s->g_rand), &(s->anti_prop),
	    &(s->g_rand) ); }
	mult_rhoi( ZUP, F_OFFSET(quark_source), F_OFFSET(anti_prop) );
	FORALLSITES(i,s){ add_su3_vector( &(s->g_rand), &(s->anti_prop),
	    &(s->g_rand) ); }
        cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );

        /* Loop over all desired sink operators */
        /*tie propagators together at sink end to project out desired mesons */
        /* add into propagator arrays at distance from source */

        mult_pion5( F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pion5_pion5][(s->t+nt-t_source)%nt], cc );
        }

        mult_pioni5( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pioni5_pioni5][(s->t+nt-t_source)%nt], cc );
        }

        mult_pioni( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pioni_pioni][(s->t+nt-t_source)%nt], cc );
        }

        mult_pions( F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pions_pions][(s->t+nt-t_source)%nt], cc );
        }

	/* 1-- (rho) PROPAGATORS ***************************************/
        mult_rhoi( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_rhoi_rhoi][(s->t+nt-t_source)%nt], cc );
        }

        mult_rhos( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_rhos_rhos][(s->t+nt-t_source)%nt], cc );
        }

	/* 2nd source for pion05, pionij, pioni0, pion0, rhoi0, rho0 */
	mult_pion05( F_OFFSET(quark_source), F_OFFSET(g_rand) );
	mult_pioni0( ZUP, F_OFFSET(quark_source), F_OFFSET(anti_prop) );
	FORALLSITES(i,s){ add_su3_vector( &(s->g_rand), &(s->anti_prop),
	    &(s->g_rand) ); }
	mult_rhoi0( ZUP, F_OFFSET(quark_source), F_OFFSET(anti_prop) );
	FORALLSITES(i,s){ add_su3_vector( &(s->g_rand), &(s->anti_prop),
	    &(s->g_rand) ); }
        cgn += mat_invert( F_OFFSET(g_rand), F_OFFSET(anti_prop), 
			   temp, amass, PRECISION, fn );

        mult_pion05( F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pion05_pion05][(s->t+nt-t_source)%nt], cc );
        }

        mult_pionij( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pionij_pionij][(s->t+nt-t_source)%nt], cc );
        }

        mult_pioni0( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pioni0_pioni0][(s->t+nt-t_source)%nt], cc );
        }

        mult_pion0( F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_pion0_pion0][(s->t+nt-t_source)%nt], cc );
        }

        mult_rhoi0( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_rhoi0_rhoi0][(s->t+nt-t_source)%nt], cc );
        }

        mult_rho0( ZUP, F_OFFSET(quark_prop), F_OFFSET(g_rand) );
        FORALLSITES(i,s){
	    cc = su3_dot( &(s->anti_prop), &(s->g_rand) );
	    CSUM( props[prop_rho0_rho0][(s->t+nt-t_source)%nt], cc );
        }

    }
  } /* end loop on t_source */

  /* Sum propagator arrays over nodes */
  /* print out propagators */
  g_veccomplexsum( props[0] , nprops*nt );
  for(i=0;i<nprops;i++)for(j=0;j<nt;j++){
    CDIVREAL(props[i][j],nx*ny*nz*n_sources,props[i][j]);
  }
  if(this_node==0){
  
    /* First source couples to pion5, pioni5, pioni, pions, rhoi, rhos */
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: FUNNYWALL1\n");
    printf("SINKS: PION_5 PION_i5 PION_i PION_s RHO_i RHO_s \n");
    for(j=0;j<nt;j++){
      printf("%d %e %e %e %e %e %e %e %e %e %e %e %e\n",j,
        props[prop_pion5_pion5][j].real, props[prop_pion5_pion5][j].imag,
        props[prop_pioni5_pioni5][j].real, props[prop_pioni5_pioni5][j].imag,
        props[prop_pioni_pioni][j].real, props[prop_pioni_pioni][j].imag,
        props[prop_pions_pions][j].real, props[prop_pions_pions][j].imag,
        props[prop_rhoi_rhoi][j].real, props[prop_rhoi_rhoi][j].imag,
        props[prop_rhos_rhos][j].real, props[prop_rhos_rhos][j].imag);
    }
    printf("ENDPROP\n");
  
    /* 2nd source for pion05, pionij, pioni0, pion0, rhoi0, rho0 */
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",qmass,amass);
    printf("SOURCE: FUNNYWALL2\n");
    printf("SINKS: PION_05 PION_ij PION_i0 PION_0 RHO_i0 RHO_0 \n");
    for(j=0;j<nt;j++){
      printf("%d %e %e %e %e %e %e %e %e %e %e %e %e\n",j,
        props[prop_pion05_pion05][j].real, props[prop_pion05_pion05][j].imag,
        props[prop_pionij_pionij][j].real, props[prop_pionij_pionij][j].imag,
        props[prop_pioni0_pioni0][j].real, props[prop_pioni0_pioni0][j].imag,
        props[prop_pion0_pion0][j].real, props[prop_pion0_pion0][j].imag,
        props[prop_rhoi0_rhoi0][j].real, props[prop_rhoi0_rhoi0][j].imag,
        props[prop_rho0_rho0][j].real, props[prop_rho0_rho0][j].imag);
    }
    printf("ENDPROP\n");

    fflush(stdout);
  } /* end if(this_node==0) */

  /* free arrays */
  free(props[0]); free(props);
  
  return(cgn);
} /* spectrum_nlpi2 */


/* "Multiply by" the quark-antiquark local pion operator */
void mult_pion5( field_offset src, field_offset dest ){
   /* operator is (-1)^(x+y+z+t), another (-1)^(x+y+z+t) for antiquark,
	so nothing to do */
    register int i;
    register site *s;
    FORALLSITES(i,s){
	*(su3_vector *)F_PT(s,dest) = *(su3_vector *)F_PT(s,src);
    }
}

/* "Multiply by" the second quark-antiquark local pion operator */
void mult_pion05( field_offset src, field_offset dest ){
   /* operator is gamma_0 (1), another (-1)^(x+y+z+t) for antiquark, 
	so just multiply by (-1)^(x+y+z+t) */
    register int i;
    register site *s;
    FORALLSITES(i,s){
	if( ((s->x+s->y+s->z+s->t)&0x1) == 0 ){
	    *(su3_vector *)F_PT(s,dest) = *(su3_vector *)F_PT(s,src);
	} else {
	    scalar_mult_su3_vector( (su3_vector *)F_PT(s,src), -1.0,
		(su3_vector *)F_PT(s,dest) );
	}
    }
}

/* "Multiply by" the one link pion operator.
   "pi_i5" or gamma_5 x gamma_i gamma_5 */
void mult_pioni5( int fdir, field_offset src, field_offset dest ){
   /* operator is gamma_0 (1), another (-1)^(x+y+z+t) for antiquark, 
	so just multiply by (-1)^(x+y+z+t) */
    register int i;
    register site *s;
    short coords[4];
    Real sign;

    /* apply the symmetric shift opperator */
    sym_shift(fdir, src, dest);
    FORALLSITES(i,s){
	coords[XUP] = s->x;
	coords[YUP] = s->y;
	coords[ZUP] = s->z;
	coords[TUP] = s->t;
	sign = +1.0;
	if( (coords[fdir]%2)==1) sign = -sign;
	if( s->parity==ODD ) sign = -sign; /* the extra gamma_5 */
	scalar_mult_su3_vector( (su3_vector *)F_PT(s,dest), sign, 
				(su3_vector *)F_PT(s,dest) );
    }
}

/* "Multiply by" the one link pion operator.
   "pi_ij" or gamma_0 gamma_5 x gamma_i gamma_j */
void mult_pionij( int fdir, field_offset src, field_offset dest ){
   /* operator is gamma_0 (1), another (-1)^(x+y+z+t) for antiquark, 
	so just multiply by (-1)^(x+y+z+t) */
    register int i;
    register site *s;
    Real sign;
    short coords[4];

    /* apply the symmetric shift opperator */
    sym_shift(fdir, src, dest);
    FORALLSITES(i,s){
	coords[XUP] = s->x;
	coords[YUP] = s->y;
	coords[ZUP] = s->z;
	coords[TUP] = s->t;
	sign = +1.0;
	if( (coords[fdir]%2)==1) sign = -sign;
	/* two gamma_5's = nothing */
	scalar_mult_su3_vector( (su3_vector *)F_PT(s,dest), sign, 
				(su3_vector *)F_PT(s,dest) );
    }
}

/* "Multiply by" the two link pion operator. *
 * "pi_4"                                    */
void mult_pioni( int fdir, field_offset src, field_offset dest )
{
  register int i,d1,d2;
  register site *s;
  Real sign, epsilon_sign;
  short coords[4];
  short dir[2] ;
  
  switch(fdir)
    {
    case 0: dir[0]=1 ; dir[1] = 2 ; break ;
    case 1: dir[0]=2 ; dir[1] = 0 ; break ;
    case 2: dir[0]=0 ; dir[1] = 1 ; break ;
    default:   node0_printf("ERROR! invalid direction %i\n",fdir);
    }
  
  /*clean up dest */
  /**/
  FORALLSITES(i,s){
    clearvec((su3_vector *)F_PT(s,dest));
  }    
  /**/
  epsilon_sign=.5 ; /* devide by two for averaging (1,2), (2,1) paths */
  for(d1=0,d2=1;d1<2;d1++,d2--)
    {
      /*      printf("In mult_pioni: (d1,d2): (%i,%i)\n",dir[d1],dir[d2]); */
      /* apply the symmetric shift opperator in dir2 */
      sym_shift(dir[d2], src, F_OFFSET(tempvec[0]));
      /* multiply by \zeta_dir2 */
      FORALLSITES(i,s){
	coords[XUP] = s->x;
	coords[YUP] = s->y;
	coords[ZUP] = s->z;
	coords[TUP] = s->t;
	/*   because the phases are on we multiply by  * 
	 *  \zeta * \eta = \epsilon * (-1)^coord[dir2] */
	if (s->parity==EVEN) 
	  sign =  1.0 ;
	else
	  sign = -1.0 ;
	if(coords[dir[d2]]%2==1) sign = -sign ;
	scalar_mult_su3_vector( &(s->tempvec[0]), sign, &(s->tempvec[0]) );
      }
      /* apply the symmetric shift opperator in dir1 */
      sym_shift(dir[d1], F_OFFSET(tempvec[0]),F_OFFSET(tempvec[1]));
      /* multiply by \zeta_dir1 */
      FORALLSITES(i,s){
	coords[XUP] = s->x;
	coords[YUP] = s->y;
	coords[ZUP] = s->z;
	coords[TUP] = s->t;
	/*   because the phases are on we multiply by  * 
	 *  \zeta * \eta = \epsilon * (-1)^coord[dir1] *
	 *  Here we have to multiply with the extra    *
	 *  \epsilon for the anti-quark propagator.    *
	 * so we are left with only  (-1)^coord[dir1]  */
	sign=1.0 ;
	if(coords[dir[d1]]%2==1) sign = -sign ;
	sign *= epsilon_sign ;
	scalar_mult_sum_su3_vector( (su3_vector *)F_PT(s,dest),
				    &(s->tempvec[1]), sign );
      }
      epsilon_sign = -epsilon_sign ;
    }
}


/* "Multiply by" the two link pion operator. */
void mult_pioni0( int fdir, field_offset src, field_offset dest )
{
  register int i;
  register site *s;

  /* We only need to multiply the pi_4 by \zeta_4\eta_4=(-1)^(x+y+z) */
  mult_pioni( fdir, src, dest ) ;
  FORALLSITES(i,s){
    if((s->x+s->y+s->z)%2==1)
      scalar_mult_su3_vector( (su3_vector *)F_PT(s,dest), -1.0, 
			      (su3_vector *)F_PT(s,dest) );
  }
}


/* "Multiply by" the three link pion operator.  */
void mult_pions(field_offset src, field_offset dest )
{
  register int i;
  register site *s;
  int c ;
  /* All permutation with appropriate sign */
  struct {
    int d[3];
    Real sign ;
  } p[6]={{{0,1,2},+1.0/6.0},
	  {{1,2,0},+1.0/6.0},
	  {{2,0,1},+1.0/6.0},
	  {{0,2,1},-1.0/6.0},
	  {{1,0,2},-1.0/6.0},
	  {{2,1,0},-1.0/6.0}}; /* The factor of 6 accounts for the *
				* multiplicity of the permutations */
  
  /*clean up dest */
  FORALLSITES(i,s){
    clearvec((su3_vector *)F_PT(s,dest));
  }    
  for(c=0;c<6;c++)
    {
      zeta_shift(3,p[c].d,src,F_OFFSET(tempvec[1])) ;
      FORALLSITES(i,s){
	scalar_mult_sum_su3_vector((su3_vector *)F_PT(s,dest), 
				   &(s->tempvec[1]), p[c].sign );
      }
    }
  /* multiply by \epsilon for the anti-quark */
  FORODDSITES(i,s)
    scalar_mult_su3_vector((su3_vector *)F_PT(s,dest), -1.0, 
			   (su3_vector *)F_PT(s,dest) );
}

/* "Multiply by" the three link pion operator.  */
void mult_pion0(field_offset src, field_offset dest )
{
  register int i;
  register site *s;
  
  /* We only need to multiply the pi_5 by \zeta_4\eta_4=(-1)^(x+y+z) */
  mult_pions( src, dest ) ;
  FORALLSITES(i,s){
    if((s->x+s->y+s->z)%2==1)
      scalar_mult_su3_vector( (su3_vector *)F_PT(s,dest), -1.0, 
			      (su3_vector *)F_PT(s,dest) );
  }
}

/* "Multiply by" the quark-antiquark local rho operator */
void mult_rhoi( int pdir,  field_offset src, field_offset dest ){
   /* operator is gamma_pdir, another (-1)^(x+y+z+t) for antiquark */
    register int i;
    register site *s;
    FORALLSITES(i,s){
	if( ((((short *)&(s->x))[pdir]) & 0x1) == 0 ){
	    *(su3_vector *)F_PT(s,dest) = *(su3_vector *)F_PT(s,src);
	} else {
	    scalar_mult_su3_vector( (su3_vector *)F_PT(s,src), -1.0,
		(su3_vector *)F_PT(s,dest) );
	}
    }
}

/* "Multiply by" the second quark-antiquark local rho operator */
void mult_rhoi0( int pdir,  field_offset src, field_offset dest ){
   /* operator is gamma_0 gamma_pdir, another (-1)^(x+y+z+t) for antiquark */
    register int i;
    register site *s;
    FORALLSITES(i,s){
	if( ((((short *)&(s->x))[pdir] + s->x+s->y+s->z) & 0x1) == 0 ){
	    *(su3_vector *)F_PT(s,dest) = *(su3_vector *)F_PT(s,src);
	} else {
	    scalar_mult_su3_vector( (su3_vector *)F_PT(s,src), -1.0,
		(su3_vector *)F_PT(s,dest) );
	}
    }
}

/* "Multiply by" the quark-antiquark one link rho operator */
void mult_rhos( int fdir,  field_offset src, field_offset dest )
{
  register int i;
  register site *s;  

  /* apply the symmetric shift opperator */
  sym_shift(fdir, src, dest);
  FORALLSITES(i,s){
    /* \eta_k already in the phases                   * 
     * only the \epsilon for the anti-quark is needed */
    if(s->parity==ODD)
      scalar_mult_su3_vector( (su3_vector *)F_PT(s,dest), -1.0,
			      (su3_vector *)F_PT(s,dest) );
  }
}

/* "Multiply by" the second quark-antiquark one link rho operator */
void mult_rho0( int fdir,  field_offset src, field_offset dest )
{
  register int i;
  register site *s;  

  /* apply the symmetric shift opperator */
  sym_shift(fdir, src, dest);
  FORALLSITES(i,s){
    /* \eta_k already in the phases                                     * 
     * only the \epsilon for the anti-quark is needed times (-1)^(x+y+z)*
     * this is equal to (-1)^t                                          */
    if((s->t)%2==1)
      scalar_mult_su3_vector( (su3_vector *)F_PT(s,dest), -1.0,
			      (su3_vector *)F_PT(s,dest) );
  }
}

/* "Multiply by" the quark-antiquark local a1 operator */
void mult_a1( int pdir,  field_offset src, field_offset dest ){
   /* operator is gamma_pdir, (-1)^(x+y+z+t), another (-1)^(x+y+z+t)
	for antiquark */
    register int i;
    register site *s;
if(this_node==0)printf("OOPS, mult_a1 NOT WRITTEN\n");
exit(0);
    FORALLSITES(i,s){
    }
}

/* "Multiply by" the quark-antiquark local b1 operator */
void mult_b1( int pdir,  field_offset src, field_offset dest ){
   /* operator is gamma_pdir, gamma_0, (-1)^(x+y+z+t), another (-1)^(x+y+z+t) for antiquark */
    register int i;
    register site *s;
if(this_node==0)printf("OOPS, mult_b1 NOT WRITTEN\n");
exit(0);
    FORALLSITES(i,s){
    }
}

