/******** spectrum_pwave.c *************/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE */
/* DT 11/21/95 started */

/* Spectrum for Wilson hybrid mesons with exotic quantum
   numbers.
   Any antiperiodic boundary conditions (time direction usually)
   are absorbed in the link matrices (ie boundary_flip(MINUS) )

    Gauge fixing should be done before calling this function
    (Coulomb gauge, probably).

    This version has P wave sources for the a_1 quark-antiquark
    operator used in constructing 0+- and 0-- sources.

    The a_1 in this version is a P-wave operator:
    The 0-+ in this version is a hybrid operator: rho for the
	quark-antiquark and color magnetic field.
	1-- dot 1+- = 0-+
	(expect to overlap the pion)
    The 1-- in this version is a hybrid operator: pion for the
	quark-antiquark and color magnetic field.
	0-+ cross 1+- = 1--
	(expect to overlap the rho)
    The 1++ in this version is a hybrid operator: rho for the
	quark-antiquark and color electric field.
	1-- cross 1-- =  ... 1++  ...
	(expect to overlap the a_1)

   BICONGRAD must be defined, since I use vtmp and sss as temporary
   storage.
   Smearing and field strength computation are assumed already done,
	most likely by spectrum_hybrids().

   multiply by field strength and matrix inverters are in spectrum_hybrids.c
*/

#include "cl_hyb_includes.h"

void mult_zero_pm_P( field_offset src, field_offset dest );
void mult_zero_mm_P( field_offset src, field_offset dest );
void mult_a1_P( int pdir, field_offset src, field_offset dest );
void mult_zero_mp( field_offset src, field_offset dest );
void mult_one_mm( int pdir, field_offset src, field_offset dest );
void mult_one_pp( int pdir, field_offset src, field_offset dest );

void mult_by_field_strength( int dir1, int dir2,
    field_offset src, field_offset dest );
int mat_invert( field_offset src, field_offset dest );
void check_invert( field_offset src, field_offset dest );
int test_converge(int t_source);

int spectrum_pwave(){ /* return the C.G. iteration number */

  int cgn;
  register int i,j;
  register site* s;
  register complex cc;
  register Real phase;
  Real finalrsq;
  register int t_source;
  int dir;	/* direction in lattice */
  int spin;	/* spin for source */
  int color;	/* color for source */
  int src_count; /* number of source time slices used */
  complex *prop_a1_P, *prop_0pm_P, *prop_0mm_P, *prop_0mp, *prop_1mm, *prop_1pp;

  cgn=0; /* number of CG iterations */

  /* allocate arrays to accumulate propagators, set them to zero */
  prop_a1_P=(complex *)malloc(nt*sizeof(complex));
  prop_0pm_P=(complex *)malloc(nt*sizeof(complex));
  prop_0mm_P=(complex *)malloc(nt*sizeof(complex));
  prop_0mp=(complex *)malloc(nt*sizeof(complex));
  prop_1mm=(complex *)malloc(nt*sizeof(complex));
  prop_1pp=(complex *)malloc(nt*sizeof(complex));
  for(cc.real=cc.imag=0.0,i=0;i<nt;i++){
    prop_a1_P[i]=cc;
    prop_0pm_P[i]=prop_0mm_P[i]=prop_0mp[i]=prop_1mm[i]=prop_1pp[i]=cc;
  }

  /* loop over "source" time slice */
  for(src_count=0,t_source=source_start; t_source<nt && src_count<n_sources;
    t_source += source_inc,src_count++){
    /* Wall source */
    /* Use quark_source for quark source */
    if(this_node==0)printf("spectrum_pwave(): source time = %d\n",t_source);

    for(spin=0;spin<4;spin++)for(color=0;color<3;color++){
        FORALLSITES(i,s){
	    clear_wvec( &(s->quark_source) );
            if(s->t==t_source) s->quark_source.d[spin].c[color].real=1.0;
	}

        /* compute M^-1 * quark_source */
        cgn += mat_invert( F_OFFSET(quark_source), F_OFFSET(quark_prop) );

        /* First the a1 operator (1++, epsilon_ijk gamma_j deriv_k ) */
        mult_a1_P( ZUP, F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_a1_P( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( prop_a1_P[(s->t+nt-t_source)%nt], cc );
        }

        /* Now the 0+- P wave operator */
        /* Source for antiquark propagator - mult_zero_pm_P() includes
	    gamma matrix for the antiquark propagator */
        mult_zero_pm_P( F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_zero_pm_P( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( prop_0pm_P[(s->t+nt-t_source)%nt], cc );
        }

        /* Now do the 0-- operator */
        /* Source for antiquark propagator (mult_zero_mm_P includes
		gamma_5 for antiquark propagator) */
        mult_zero_mm_P( F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_zero_mm_P( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( prop_0mm_P[(s->t+nt-t_source)%nt], cc );
        }

        /* Now do the 0-+ hybrid operator.  */
        /* Source for antiquark propagator */
        mult_zero_mp( F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_zero_mp( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( prop_0mp[(s->t+nt-t_source)%nt], cc );
        }

        /* Now do the 1-- hybrid operator.  For the moment, Z component only */
        /* Source for antiquark propagator */
        mult_one_mm( ZUP, F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_one_mm( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( prop_1mm[(s->t+nt-t_source)%nt], cc );
        }

        /* Now do the 1++ hybrid operator.  For the moment, Z component only */
        /* Source for antiquark propagator */
        mult_one_pp( ZUP, F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_one_pp( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( prop_1pp[(s->t+nt-t_source)%nt], cc );
        }

    }
  } /* end loop on t_source */

  /* Sum propagator arrays over nodes */
  /* print out propagators */
  g_veccomplexsum( prop_a1_P , nt );
  g_veccomplexsum( prop_0pm_P , nt );
  g_veccomplexsum( prop_0mm_P , nt );
  g_veccomplexsum( prop_0mp , nt );
  g_veccomplexsum( prop_1mm , nt );
  g_veccomplexsum( prop_1pp , nt );
  for(i=0;i<nt;i++){
    CDIVREAL(prop_a1_P[i]  ,nx*ny*nz*src_count,prop_a1_P[i]);
    CDIVREAL(prop_0pm_P[i] ,nx*ny*nz*src_count,prop_0pm_P[i]);
    CDIVREAL(prop_0mm_P[i] ,nx*ny*nz*src_count,prop_0mm_P[i]);
    CDIVREAL(prop_0mp[i] ,nx*ny*nz*src_count,prop_0mp[i]);
    CDIVREAL(prop_1mm[i] ,nx*ny*nz*src_count,prop_1mm[i]);
    CDIVREAL(prop_1pp[i] ,nx*ny*nz*src_count,prop_1pp[i]);
    if(this_node==0){
	printf("HYBRI2_PROP: %d %e %e %e %e %e %e\n",i,
          prop_a1_P[i].real, prop_a1_P[i].imag,
	  prop_0mm_P[i].real, prop_0mm_P[i].imag,
          prop_0pm_P[i].real, prop_0pm_P[i].imag );
	printf("HYBRI3_PROP: %d %e %e %e %e %e %e\n",i,
          prop_1mm[i].real,prop_1mm[i].imag,
          prop_1pp[i].real,prop_1pp[i].imag,
          prop_0mp[i].real, prop_0mp[i].imag );
    }
  }

  /* free arrays */
  free(prop_a1_P); free(prop_0pm_P); free(prop_0mm_P);
  free(prop_0mp); free(prop_1mm); free(prop_1pp);
  
  return(cgn);
} /* spectrum_pwave */


/* "Multiply by" the zero-plus-minus P wave operator 
	quark operator is a1 = epsilon_ijk gamma_j deriv_k = 1++,
	deriv is forward_deriv - backward_deriv
	gluon operator is B_i = 1+-,  take dot product 
	see "sources.tex" for details. */
    /* Uses MP, vtmp and sss for temporary storage */
void mult_zero_pm_P( field_offset src, field_offset dest ){
    register int dir1,dir2,i;
    register site *s;
    wilson_vector tvec1,tvec2;
    msg_tag *tag0,*tag1;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( F_PT(s,dest) );
    }

    /* Loop over spatial directions */
    for(dir1=XUP;dir1<=ZUP;dir1++)for(dir2=XUP;dir2<=ZUP;dir2++)if(dir1!=dir2){

        /* multiply by gamma_dir1 (for a1) and gamma_five (for
           antiquark propagator ) */
        FORALLSITES(i,s){
            mult_by_gamma( (wilson_vector *)F_PT(s,src), &tvec1, dir1 );
            mult_by_gamma( &tvec1, &(s->sss), GAMMAFIVE );
        }

	/* parallel transport sss from positive dir2 direction */
	/* parallel transport sss, from negative dir2 direction */
        tag0=start_gather_site( F_OFFSET(sss), sizeof(wilson_vector),
            dir2, EVENANDODD, gen_pt[0] );
	FORALLSITES(i,s){
	    mult_adj_mat_wilson_vec( &(s->link[dir2]), &(s->sss),  &(s->vtmp) );
	}
        tag1=start_gather_site( F_OFFSET(vtmp), sizeof(wilson_vector),
            OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
	wait_gather(tag0);
	wait_gather(tag1);
	FORALLSITES(i,s){
	    mult_mat_wilson_vec( &(s->link[dir2]),
		(wilson_vector *)gen_pt[0][i], &tvec2 );
            sub_wilson_vector( &tvec2, (wilson_vector *)gen_pt[1][i], &(s->MP));
	}
	cleanup_gather(tag0);
	cleanup_gather(tag1);
        mult_by_field_strength( dir1, dir2, F_OFFSET(MP), F_OFFSET(vtmp) );
	FORALLSITES(i,s){
	    add_wilson_vector( F_PT(s,dest), &(s->vtmp), F_PT(s,dest) );
	}

        /* Multiply sss by dir1,dir2 component of magnetic field, */
	/* keep in temp. vector mp */
        mult_by_field_strength( dir1, dir2, F_OFFSET(sss), F_OFFSET(MP) );

	/* parallel transport MP from positive dir2 direction */
	/* parallel transport MP from negative dir2 direction */
        tag0=start_gather_site( F_OFFSET(MP), sizeof(wilson_vector),
            dir2, EVENANDODD, gen_pt[0] );
	FORALLSITES(i,s){mult_adj_mat_wilson_vec( &(s->link[dir2]), &(s->MP), &(s->vtmp) );
	}
        tag1=start_gather_site( F_OFFSET(vtmp), sizeof(wilson_vector),
            OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
	wait_gather(tag0);
	wait_gather(tag1);
	FORALLSITES(i,s){
	    mult_mat_wilson_vec( &(s->link[dir2]),
		(wilson_vector *)gen_pt[0][i], &tvec1 );
            sub_wilson_vector( &tvec1, (wilson_vector *)gen_pt[1][i], &tvec1 );
            add_wilson_vector( (wilson_vector *)F_PT(s,dest), &tvec1,
                (wilson_vector *)F_PT(s,dest) );
	}
	cleanup_gather(tag0);
	cleanup_gather(tag1);

    } /* end loops on dir1 and dir2 */
}


/* "Multiply by" the zero-minus-minus P wave operator 
	quark operator is a1 = epsilon_ijk gamma_j  deriv_k = 1++,
	gluon operator is E_i = 1--,  take dot product */
void mult_zero_mm_P( field_offset src, field_offset dest ){
    register int in,i,j,k;
    register site *s;
    wilson_vector tvec1,tvec2;
    msg_tag *tag0,*tag1;
    register Real sign;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( F_PT(s,dest) );
    }
    /* loop over directions, sign is sign of epsilon tensor */
    for(i=XUP;i<=ZUP;i++)for(j=XUP;j<=ZUP;j++)for(k=XUP;k<=ZUP;k++){
	if( j==k || j==i || k== i )continue;
	else if ( j == ((i+1)%3) )sign = 1.0;
	else sign = -1.0;

        /* multiply by gamma_j (for a1) and gamma_five (for
           antiquark propagator ) */
        FORALLSITES(in,s){
            mult_by_gamma( (wilson_vector *)F_PT(s,src), &tvec1, j );
            mult_by_gamma( &tvec1, &(s->sss), GAMMAFIVE );
        }

	/* parallel transport sss from positive k direction */
	/* parallel transport sss, from negative k direction */
	/* subtract, and multiply by field strength at site */
        tag0=start_gather_site( F_OFFSET(sss), sizeof(wilson_vector),
            k, EVENANDODD, gen_pt[0] );
	FORALLSITES(in,s){
	    mult_adj_mat_wilson_vec( &(s->link[k]), &(s->sss),  &(s->vtmp) );
	}
        tag1=start_gather_site( F_OFFSET(vtmp), sizeof(wilson_vector),
            OPP_DIR(k), EVENANDODD, gen_pt[1] );
	wait_gather(tag0);
	wait_gather(tag1);
	FORALLSITES(in,s){
	    mult_mat_wilson_vec( &(s->link[k]),
		(wilson_vector *)gen_pt[0][in], &tvec2 );
            sub_wilson_vector( &tvec2, (wilson_vector *)gen_pt[1][in],&(s->MP));
	}
	cleanup_gather(tag0);
	cleanup_gather(tag1);
        mult_by_field_strength( TUP, i, F_OFFSET(MP), F_OFFSET(vtmp) );
	FORALLSITES(in,s){
	    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest), &(s->vtmp),
		sign, (wilson_vector *)F_PT(s,dest) );
	}

        /* Multiply sss by 0,i component of magnetic field, */
	/* keep in temp. vector MP */
        mult_by_field_strength( TUP, i, F_OFFSET(sss), F_OFFSET(MP) );

	/* parallel transport MP from positive k direction */
	/* parallel transport MP from negative k direction */
        tag0=start_gather_site( F_OFFSET(MP), sizeof(wilson_vector),
            k, EVENANDODD, gen_pt[0] );
	FORALLSITES(in,s){
	    mult_adj_mat_wilson_vec( &(s->link[k]), &(s->MP), &(s->vtmp) );
	}
        tag1=start_gather_site( F_OFFSET(vtmp), sizeof(wilson_vector),
            OPP_DIR(k), EVENANDODD, gen_pt[1] );
	wait_gather(tag0);
	wait_gather(tag1);
	FORALLSITES(in,s){
	    mult_mat_wilson_vec( &(s->link[k]),
		(wilson_vector *)gen_pt[0][in], &tvec1 );
            sub_wilson_vector( &tvec1, (wilson_vector *)gen_pt[1][in], &tvec1 );
            scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest), &tvec1, sign,
                (wilson_vector *)F_PT(s,dest) );
	}
	cleanup_gather(tag0);
	cleanup_gather(tag1);
    } /* end loops on directions i,j,k */

}


/* "Multiply by" the a1 P wave operator */
void mult_a1_P( int pdir,  field_offset src, field_offset dest ){
    register int i,j,k;
    register site *s;
    wilson_vector tvec1;
    Real sign;
    msg_tag *tag0,*tag1;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( F_PT(s,dest) );
    }
    /* loop over directions, sign is sign of epsilon tensor */
    for(j=XUP;j<=ZUP;j++)for(k=XUP;k<=ZUP;k++){
	if( j==k || j==pdir || k== pdir )continue;
	else if ( j == ((pdir+1)%3) )sign = 1.0;
	else sign = -1.0;

        /* multiply by gamma_j (for a1) and gamma_five (for
           antiquark propagator ) */
        FORALLSITES(i,s){
            mult_by_gamma( (wilson_vector *)F_PT(s,src), &tvec1, j );
            mult_by_gamma( &tvec1, &(s->MP), GAMMAFIVE );
        }
	/* parallel transport MP from forwards and backwards. */
        tag0=start_gather_site( F_OFFSET(MP), sizeof(wilson_vector),
            k, EVENANDODD, gen_pt[0] );
	FORALLSITES(i,s){
	    mult_adj_mat_wilson_vec( &(s->link[k]), &(s->MP), &(s->vtmp) );
	}
        tag1=start_gather_site( F_OFFSET(vtmp), sizeof(wilson_vector),
            OPP_DIR(k), EVENANDODD, gen_pt[1] );
	wait_gather(tag0);
	wait_gather(tag1);
	FORALLSITES(i,s){
	    mult_mat_wilson_vec( &(s->link[k]),
		(wilson_vector *)gen_pt[0][i], &tvec1 );
	    sub_wilson_vector( &tvec1, (wilson_vector *)gen_pt[1][i], &tvec1);
            scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest),
		&tvec1, sign, (wilson_vector *)F_PT(s,dest) );
	}
	cleanup_gather(tag0);
	cleanup_gather(tag1);

    } /* end loops on j and k (directions) */
}

/* "Multiply by" the one-minus-minus operator 
    quark operator is "pion", gluon operator is magnetic field.
   "pdir" is the polarization direction of the meson */
void mult_one_mm( int pdir, field_offset src, field_offset dest ){
    /* use mp as temporary storage */
    register int dir,i;
    register site *s;
    wilson_vector tvec1,tvec2;

    /* multiply by magnetic field.  gamma_5 for antiquark cancels
	gamma_5 in pion operator. */
    mult_by_field_strength( (pdir+1)%3, (pdir+2)%3, src, dest );
} /* end mult_one_mm */


/* "Multiply by" the one-plus-plus operator 
    quark operator is "rho", gluon operator is electric field.
    1++ = epsilon_ijk \psibar gamma_j \psi F_{0,k}
   "pdir" is the polarization direction of the meson */
void mult_one_pp( int pdir, field_offset src, field_offset dest ){
    /* use MP as temporary storage */
    register int i,j,k;
    register site *s;
    wilson_vector tvec;
    register Real sign;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( F_PT(s,dest) );
    }
    /* loop over directions, sign is sign of epsilon tensor */
    for(j=XUP;j<=ZUP;j++)for(k=XUP;k<=ZUP;k++){
	if( j==k || j==pdir || k== pdir )continue;
	else if ( j == ((pdir+1)%3) )sign = 1.0;
	else sign = -1.0;

        mult_by_field_strength( TUP, k, src, F_OFFSET(MP) );
	FORALLSITES(i,s){
            mult_by_gamma( &(s->MP), &tvec, j );
	    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest), &tvec, sign,
		(wilson_vector *)F_PT(s,dest) );
	}

    }/* end loop over j,k (directions) */

    /* multiply everything by gamma_5 for antiquarks */
    FORALLSITES(i,s){
        mult_by_gamma( F_PT(s,dest), &tvec, GAMMAFIVE );
	*(wilson_vector *)F_PT(s,dest) = tvec;
    }
} /* end mult_one_pp */

/* "Multiply by" the zero-minum-plus operator 
    quark operator is "rho", gluon operator is magnetic field.
    0-+ = epsilon_ijk \psibar gamma_i \psi F_{j,k} */
void mult_zero_mp( field_offset src, field_offset dest ){
    /* use MP as temporary storage */
    register int i,j,k,in;
    register site *s;
    wilson_vector tvec;
    register Real sign;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( F_PT(s,dest) );
    }
    /* loop over directions, sign is sign of epsilon tensor */
    for(i=XUP;i<=ZUP;i++){
	/* antisymmetry of epsilon and F_{jk} means no need to sum all terms */
	j=(i+1)%3; k = (i+2)%3;

        mult_by_field_strength( j, k, src, F_OFFSET(MP) );
	FORALLSITES(in,s){
            mult_by_gamma( &(s->MP), &tvec, i );
	    add_wilson_vector( (wilson_vector *)F_PT(s,dest), &tvec,
		(wilson_vector *)F_PT(s,dest) );
	}

    }/* end loop over i,j,k (directions) */

    /* multiply everything by gamma_5 for antiquarks */
    FORALLSITES(in,s){
        mult_by_gamma( F_PT(s,dest), &tvec, GAMMAFIVE );
	*(wilson_vector *)F_PT(s,dest) = tvec;
    }
} /* end mult_zero_mp */
