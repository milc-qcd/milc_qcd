/******** mat_invert.c *************/
/* MIMD version 7*/
/* DT 6/97
* separated from spectrum_mom.c 12/97
* Modify check_invert() to compute magnitude of error 11/98 DT
*/
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.  Assumes "dslash" has
   been defined to be the appropriate "dslash_fn" or "dslash_eo"
*/
#include "generic_ks_includes.h"
#include "../include/dslash_ks_redefine.h"

/* dst = M src. With parity selection */

void ks_dirac_op( su3_vector *src, su3_vector *dst, Real mass, 
		  int parity, imp_ferm_links_t *fn){
    register int i;
    register site *s;

    dslash_field( src, dst, parity, fn);
    FORSOMEPARITY(i,s,parity){
      scalar_mult_add_su3_vector( dst+i, src+i, +2.0*mass, dst+i);
    }
}

/* dst = Madj src. With parity selection */

void ks_dirac_adj_op( su3_vector *src, su3_vector *dst, Real mass,
		      int parity, imp_ferm_links_t *fn){
    register int i;
    register site *s;

    dslash_field( src, dst, parity, fn);
    FORSOMEPARITY(i,s,parity){
      scalar_mult_su3_vector( dst+i, -1.0, dst+i);
      scalar_mult_add_su3_vector( dst+i, src+i, 2.0*mass, dst+i);
    }
}

/* dst = Madj dst. With parity selection */

void ks_dirac_adj_op_inplace( su3_vector *dst, Real mass,
			      int parity, imp_ferm_links_t *fn){
    register int i;
    register site *s;
    su3_vector *tvec = create_v_field();

    dslash_field( dst, tvec, parity, fn);
    FORSOMEPARITY(i,s,parity){
      scalar_mult_su3_vector( dst+i, 2.0*mass, dst+i);
      scalar_mult_add_su3_vector( dst+i, tvec+i, -1.0, dst+i);
    }
    destroy_v_field(tvec);
}

/* dst = Madj M src with parity selection */

void ks_dirac_opsq( su3_vector *src, su3_vector *dst, Real mass, int parity,
		    imp_ferm_links_t *fn){
    register int i;
    register site *s;
    int otherparity = 0;
    Real msq_x4 = 4.0*mass*mass;
    su3_vector *tmp;

    tmp = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
    if(tmp==NULL){
      printf("ks_dirac_opsq(%d): no room for tmp\n",this_node);
      terminate(1);
    }

    switch(parity){
    case(EVEN): otherparity=ODD; break;
    case(ODD):  otherparity=EVEN; break;
    case(EVENANDODD): otherparity=EVENANDODD;
    }

    dslash_field( src, tmp, otherparity, fn);
    dslash_field( tmp, dst, parity, fn);
    FORSOMEPARITY(i,s,parity){
      scalar_mult_su3_vector( dst+i, -1.0, dst+i);
      scalar_mult_add_su3_vector( dst+i, src+i, msq_x4, dst+i);
    }

    free(tmp);
}

/* This algorithm solves even and odd sites separately */

int mat_invert_cg_field(su3_vector *src, su3_vector *dst, 
			 quark_invert_control *qic,
			 Real mass, imp_ferm_links_t *fn ){
    int cgn;
    su3_vector *tmp;

    tmp = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
    if(tmp==NULL){
      printf("mat_invert_cg_field(%d): no room for tmp\n",this_node);
      terminate(1);
    }

    /* "Precondition" both even and odd sites */
    /* temp <- M_adj * src */
    /* The following operation is done in the prevailing
       precision.  The algorithm needs to be fixed! */
    ks_dirac_adj_op( src, tmp, mass, EVENANDODD, fn);

    /* We don't call with EVENANDODD anymore because we are
       transitioning to the QOP/QDP standard */

    /* dst_e <- (M_adj M)^-1 temp_e  (even sites only) */
    qic->parity = EVEN;
    cgn = ks_congrad_field( tmp, dst, qic, mass, fn );

    /* dst_o <- (M_adj M)^-1 temp_o  (odd sites only) */
    qic->parity = ODD;
    cgn += ks_congrad_field( tmp, dst, qic, mass, fn );

    free(tmp);

    //    check_invert_field2( dst, tmp, mass, 1e-6, fn);
    //    check_invert_field( dst, src, mass, 1e-6, fn);
    return cgn;
}

/* Compute M^-1 * phi, answer in dest
  Uses phi, ttt, resid, xxx, and cg_p as workspace */
int mat_invert_cg( field_offset src, field_offset dest, field_offset temp,
		   Real mass, int prec, imp_ferm_links_t *fn ){
    int cgn;
    quark_invert_control qic;
    int i; site *s;

    qic.prec       = prec;
    qic.min        = 0;
    qic.max        = niter;
    qic.nrestart   = nrestart;
    qic.parity     = EVENANDODD;
    qic.start_flag = 1;
    qic.nsrc = 1;
    qic.resid      = sqrt(rsqprop);
    qic.relresid   = 0;

    su3_vector *tsrc = create_v_field_from_site_member(src);
    su3_vector *tdest = create_v_field_from_site_member(dest);

    cgn = mat_invert_cg_field( tsrc, tdest, &qic, mass, fn );

    FORALLSITES(i,s){
      su3vec_copy(tdest+i, (su3_vector *)F_PT(s,dest));
    }

    destroy_v_field(tdest);
    destroy_v_field(tsrc);

    return cgn;
}

/* Invert using Leo's UML trick */

/* Our M is     (  2m		D_eo   )
		( -D_eo^adj	2m )

    Note D_oe = -D_eo^adj 

   define U =	( 2m		-D_eo )
		(  0		2m    )

	  L =   ( 2m		0  )
		( D_eo^adj	2m )

   observe UML =    2m ( 4m^2+D_eo D_eo^adj	0    )
		       (  0			4m^2 )
   and the upper left corner is what our conjugate gradient(even)
   inverts, except for a funny minus sign in our cong. grad. and a
   factor of 2m

   Use M^-1 phi = L L^-1 M^-1 U^-1 U phi
		= L  (UML)^-1 U phi
		= L  -(1/2m)(congrad) U phi

   Note that

      M^-1 = ( 2m A          - A D_eo )
             ( - B D_oe        2m B   )

where  A = (4m^2+D_eo D_eo^adj)^-1 and B = (4m^2+D_oe^adj D_oe)^-1

*/
         
/* This algorithm solves even sites, reconstructs odd and then polishes
   to compensate for loss of significance in the reconstruction
*/

int mat_invert_uml_field(su3_vector *src, su3_vector *dst, 
			 quark_invert_control *qic,
			 Real mass, imp_ferm_links_t *fn ){
    int cgn;
    register int i;
    register site *s;
    su3_vector *tmp = create_v_field();
    su3_vector *ttt = create_v_field();
    int even_iters;

    /* "Precondition" both even and odd sites */
    /* temp <- M_adj * src */

    ks_dirac_adj_op( src, tmp, mass, EVENANDODD, fn );

    /* dst_e <- (M_adj M)^-1 tmp_e  (even sites only) */
    qic->parity     = EVEN;
    cgn = ks_congrad_field( tmp, dst, qic, mass, fn );
    even_iters = qic->final_iters;

    /* reconstruct odd site solution */
    /* dst_o <-  1/2m (Dslash_oe*dst_e + src_o) */
    dslash_field( dst, ttt, ODD, fn );
    FORODDSITES(i,s){
      sub_su3_vector( src+i, ttt+i, dst+i);
      scalar_mult_su3_vector( dst+i, 1.0/(2.0*mass), dst+i );
    }

    /* Polish off odd sites to correct for possible roundoff error */
    /* dst_o <- (M_adj M)^-1 temp_o  (odd sites only) */
    qic->parity = ODD;
    cgn = ks_congrad_field( tmp, dst, qic, mass, fn );
    qic->final_iters += even_iters;

    //    check_invert_field( dst, src, mass, 1e-6, fn);
    destroy_v_field(tmp);
    destroy_v_field(ttt);

    return cgn;
}

int mat_invert_uml(field_offset src, field_offset dest, field_offset temp,
		   Real mass, int prec, imp_ferm_links_t *fn ){
    int cgn;
    register int i;
    register site *s;
    quark_invert_control qic;
    su3_vector *tsrc = create_v_field_from_site_member(src);
    su3_vector *tdest = create_v_field_from_site_member(dest);

    qic.prec       = prec;
    qic.min        = 0;
    qic.max        = niter;
    qic.nrestart   = nrestart;
    qic.start_flag = 0;
    qic.nsrc       = 1;
    qic.resid      = sqrt(rsqprop);
    qic.relresid   = 0;

    cgn = mat_invert_uml_field(tsrc, tdest, &qic, mass, fn);

    FORALLSITES(i,s){
      su3vec_copy(tdest+i, (su3_vector *)F_PT(s,dest));
    }

    // check_invert( dest, src, mass, 1e-6, fn);

    destroy_v_field(tdest);
    destroy_v_field(tsrc);
    return cgn;
}

/* FOR TESTING: multiply src by matrix and check against dest */
void check_invert_field( su3_vector *src, su3_vector *dest, Real mass,
			 Real tol, imp_ferm_links_t *fn){
    register int i,k,flag;
    register site *s;
    Real r_diff, i_diff;
    double sum,sum2,dflag,dmaxerr,derr;
    su3_vector *tmp;

    tmp = create_v_field();

    /* Compute tmp = M src */
    ks_dirac_op( src, tmp, mass, EVENANDODD, fn);

    sum2=sum=0.0;
    dmaxerr=0;
    flag = 0;
    FORALLSITES(i,s){
	for(k=0;k<3;k++){
	    r_diff = dest[i].c[k].real - tmp[i].c[k].real;
	    i_diff = dest[i].c[k].imag - tmp[i].c[k].imag;
	    if( fabs(r_diff) > tol || fabs(i_diff) > tol ){
	      printf("site %d color %d  expected ( %.4e , %.4e ) got ( %.4e , %.4e )\n",
		     i,k,
		     dest[i].c[k].real, dest[i].c[k].imag,
		     tmp[i].c[k].real, tmp[i].c[k].imag);
	      flag++;
	    }
	    derr = r_diff*r_diff + i_diff*i_diff;
	    if(derr>dmaxerr)dmaxerr=derr;
 	    sum += derr;
	}
	sum2 += magsq_su3vec( dest+i );
    }
    g_doublesum( &sum );
    g_doublesum( &sum2 );
    dflag=flag;
    g_doublesum( &dflag );
    g_doublemax( &dmaxerr );
    if(this_node==0){
      if(sum2 > 0)
	printf("Inversion checked, frac. error = %e\n",sqrt(sum/sum2));
      printf("Flagged comparisons = %d\n",(int)dflag);
      printf("Max err. = %e",sqrt(dmaxerr));
      if(sum2 > 0)
	printf(" frac. = %e",sqrt(dmaxerr*volume/sum2));
      printf("\n");
      fflush(stdout);
    }
    destroy_v_field(tmp);
}

/* FOR TESTING: multiply src by matrix and check against dest */
void check_invert( field_offset src, field_offset dest, Real mass,
		   Real tol, imp_ferm_links_t *fn){

  su3_vector *tsrc = create_v_field_from_site_member(src);
  su3_vector *tdest = create_v_field_from_site_member(dest);
  
  check_invert_field( tsrc, tdest, mass, tol, fn );
  
  destroy_v_field(tdest);
  destroy_v_field(tsrc);
}

/* FOR TESTING: multiply src by Madj M and check against dest */
void check_invert_field2( su3_vector *src, su3_vector *dest, Real mass,
			  Real tol, imp_ferm_links_t *fn){
    register int i,k,flag;
    register site *s;
    Real r_diff, i_diff;
    double sum,sum2,dflag,dmaxerr,derr;
    su3_vector *tmp;

    tmp = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
    if(tmp==NULL){
      printf("check_invert_field2(%d): no room for tmp\n",this_node);
      terminate(1);
    }

    /* Compute tmp = (Madj M) src */
    ks_dirac_opsq( src, tmp, mass, EVENANDODD, fn);

    sum2=sum=0.0;
    dmaxerr=0;
    flag = 0;
    FORALLSITES(i,s){
	for(k=0;k<3;k++){
	    r_diff = dest[i].c[k].real - tmp[i].c[k].real;
	    i_diff = dest[i].c[k].imag - tmp[i].c[k].imag;
	    if( fabs(r_diff) > tol || fabs(i_diff) > tol ){
	      printf("site %d color %d  expected ( %.4e , %.4e ) got ( %.4e , %.4e )\n",
		     i,k,
		     dest[i].c[k].real, dest[i].c[k].imag,
		     tmp[i].c[k].real, tmp[i].c[k].imag);
	      flag++;
	    }
	    derr = r_diff*r_diff + i_diff*i_diff;
	    if(derr>dmaxerr)dmaxerr=derr;
 	    sum += derr;
	}
	sum2 += magsq_su3vec( dest+i );
    }
    g_doublesum( &sum );
    g_doublesum( &sum2 );
    dflag=flag;
    g_doublesum( &dflag );
    g_doublemax( &dmaxerr );
    if(this_node==0){
      printf("Inversion checked, frac. error = %e\n",sqrt(sum/sum2));
      printf("Flagged comparisons = %d\n",(int)dflag);
      printf("Max err. = %e frac. = %e\n",sqrt(dmaxerr),
	     sqrt(dmaxerr*volume/sum2));
      fflush(stdout);
    }
    free(tmp);
}

