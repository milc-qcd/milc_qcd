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

/* dst = M src. No parity selection here */

void ks_dirac_op( su3_vector *src, su3_vector *dst, Real mass, 
		  ferm_links_t *fn){
    register int i;
    register site *s;

    dslash_field( src, dst, EVENANDODD, fn);
    FORALLSITES(i,s){
      scalar_mult_add_su3_vector( dst+i, src+i, +2.0*mass, dst+i);
    }
}

/* dst = Madj src. No parity selection here */

void ks_dirac_adj_op( su3_vector *src, su3_vector *dst, Real mass,
		      ferm_links_t *fn){
    register int i;
    register site *s;

    dslash_field( src, dst, EVENANDODD, fn);
    FORALLSITES(i,s){
      scalar_mult_su3_vector( dst+i, -1.0, dst+i);
      scalar_mult_add_su3_vector( dst+i, src+i, 2.0*mass, dst+i);
    }
}

/* dst = Madj M src with parity selection */

void ks_dirac_opsq( su3_vector *src, su3_vector *dst, Real mass, int parity,
		    ferm_links_t *fn){
    register int i;
    register site *s;
    int otherparity = 0;
    Real msq_x4 = 4.0*mass*mass;
    su3_vector *tmp;

    tmp = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
    if(tmp==NULL){
      printf("ks_dirac_op(%d): no room for tmp\n",this_node);
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

/* Compute M^-1 * phi, answer in dest
  Uses phi, ttt, resid, xxx, and cg_p as workspace */
int mat_invert_cg( field_offset src, field_offset dest, field_offset temp,
		   Real mass, int prec, ferm_links_t *fn ){
    int cgn;
    quark_invert_control qic;

    qic.prec       = prec;
    qic.max        = niter;
    qic.nrestart   = nrestart;
    qic.resid      = rsqprop;
    qic.relresid   = 0;
    qic.parity     = EVENANDODD;

    //    clear_latvec( dest, EVENANDODD );
    /* Multiply SOURCE by Madjoint */
    dslash_site( src, F_OFFSET(ttt), EVENANDODD, fn );
    scalar_mult_add_latvec( F_OFFSET(ttt), src,
       -2.0*mass, F_OFFSET(ttt), EVENANDODD);
    scalar_mult_latvec( F_OFFSET(ttt), -1.0, F_OFFSET(ttt), EVENANDODD );
    cgn = ks_congrad_site( F_OFFSET(ttt), dest, &qic, mass, fn );
    return cgn;
}


/* This algorithm solves even and odd sites separately */

int mat_invert_cg_field(su3_vector *src, su3_vector *dst, 
			 quark_invert_control *qic,
			 Real mass, ferm_links_t *fn ){
    int cgn;
    su3_vector *tmp;

    tmp = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
    if(tmp==NULL){
      printf("mat_invert_cg_field(%d): no room for tmp\n",this_node);
      terminate(1);
    }

    /* "Precondition" both even and odd sites */
    /* temp <- M_adj * src */
    ks_dirac_adj_op( src, tmp, mass, fn);

    /* We don't call with EVENANDODD anymore because we are
       transitioning to the QOP/QDP standard */

    /* dst_e <- (M_adj M)^-1 temp_e  (even sites only) */
    qic->parity = EVEN;
    cgn = ks_congrad_field( tmp, dst, qic, mass, fn );

    /* dst_o <- (M_adj M)^-1 temp_o  (odd sites only) */
    qic->parity = ODD;
    cgn = ks_congrad_field( tmp, dst, qic, mass, fn );

    free(tmp);

    // check_invert_field2( dst, tmp, mass, 1e-6, fn);
    // check_invert_field( dst, src, mass, 1e-6, fn);
    return cgn;
}

/* Compute M^-1 * phi, answer in dest
  Uses phi, ttt, resid, xxx, and cg_p as workspace */
int mat_invert_cg_odd( field_offset src, field_offset dest, field_offset temp,
		   Real mass, int prec, ferm_links_t *fn ){
    int cgn;
    quark_invert_control qic;

    qic.prec       = prec;
    qic.max        = niter;
    qic.nrestart   = nrestart;
    qic.resid      = rsqprop;
    qic.relresid   = 0;
    qic.parity     = ODD;

    //    clear_latvec( dest, EVENANDODD );
    /* Multiply SOURCE by Madjoint ODD sites only */
    dslash_site( src, F_OFFSET(ttt), ODD, fn);
    scalar_mult_add_latvec( F_OFFSET(ttt), src,
       -2.0*mass, F_OFFSET(ttt), ODD);
    scalar_mult_latvec( F_OFFSET(ttt), -1.0, F_OFFSET(ttt), ODD );
    cgn = ks_congrad_site( F_OFFSET(ttt), dest, &qic, mass, fn);
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
*/
         
int mat_invert_uml_old(field_offset src, field_offset dest, field_offset temp,
		       Real mass, int prec, ferm_links_t *fn )
{
    int cgn;
    Real finalrsq;
    register int i;
    register site *s;

    if( src==temp ){
	printf("BOTCH\n"); exit(0);
    }
    /* multiply by U - even sites only */
    dslash_site( src, F_OFFSET(ttt), EVEN, fn);
    scalar_mult_add_latvec( F_OFFSET(ttt), src,
       -2.0*mass, temp, EVEN);
    scalar_mult_latvec( temp, -1.0, temp, EVEN);
    /* invert with M_adj M even */
    cgn = ks_congrad( temp, dest, mass, niter, nrestart, rsqprop,
		      prec, EVEN, &finalrsq, fn );
    /* multiply by (1/2m)L, does nothing to even sites */
    /* fix up odd sites , 1/2m (Dslash_oe*dest_e + phi_odd) */
    dslash_site( dest, F_OFFSET(ttt), ODD, fn );
    FORODDSITES(i,s){
	sub_su3_vector( (su3_vector *)F_PT(s,src), &(s->ttt), 
	    (su3_vector *)F_PT(s,dest) );
	scalar_mult_su3_vector( (su3_vector *)F_PT(s,dest), 1.0/(2.0*mass),
	    (su3_vector *)F_PT(s,dest) );
    }
    return cgn;
}

/* This algorithm solves even sites, reconstructs odd and then polishes
   to compensate for loss of significance in the reconstruction
*/

int mat_invert_uml_field(su3_vector *src, su3_vector *dst, 
			 quark_invert_control *qic,
			 Real mass, ferm_links_t *fn ){
    int cgn;
    register int i;
    register site *s;
    su3_vector *tmp, *ttt;

    tmp = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
    if(tmp==NULL){
      printf("mat_invert_uml_field(%d): no room for tmp\n",this_node);
      terminate(1);
    }
    ttt = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
    if(ttt==NULL){
      printf("mat_invert_uml_field(%d): no room for ttt\n",this_node);
      terminate(1);
    }
    /* "Precondition" both even and odd sites */
    /* temp <- M_adj * src */
    dslash_field( src, ttt, EVENANDODD, fn);
    FORALLSITES(i,s){
      scalar_mult_add_su3_vector( ttt+i, src+i,-2.0*mass, tmp+i);
      scalar_mult_su3_vector( tmp+i, -1.0, tmp+i);
    }

    /* dst_e <- (M_adj M)^-1 temp_e  (even sites only) */
    qic->parity     = EVEN;
    cgn = ks_congrad_field( tmp, dst, qic, mass, fn );

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

    // check_invert_field( dst, src, mass, 1e-6, fn);
    free(tmp);
    free(ttt);

    return cgn;
}

int mat_invert_uml(field_offset src, field_offset dest, field_offset temp,
		   Real mass, int prec, ferm_links_t *fn ){
    int cgn;
    register int i;
    register site *s;
    quark_invert_control qic;

    qic.prec       = prec;
    qic.max        = niter;
    qic.nrestart   = nrestart;
    qic.resid      = rsqprop;
    qic.relresid   = 0;

    if( src==temp ){
	printf("BOTCH\n"); exit(0);
    }
    /* "Precondition" both even and odd sites */
    /* temp <- M_adj * src */
    dslash_site( src, F_OFFSET(ttt), EVENANDODD, fn);
    scalar_mult_add_latvec( F_OFFSET(ttt), src,
       -2.0*mass, temp, EVENANDODD);
    scalar_mult_latvec( temp, -1.0, temp, EVENANDODD);

    /* dest_e <- (M_adj M)^-1 temp_e  (even sites only) */
    qic.parity     = EVEN;
    cgn = ks_congrad_site( temp, dest, &qic, mass, fn );

    /* reconstruct odd site solution */
    /* dest_o <-  1/2m (Dslash_oe*dest_e + src_o) */
    dslash_site( dest, F_OFFSET(ttt), ODD, fn );
    FORODDSITES(i,s){
	sub_su3_vector( (su3_vector *)F_PT(s,src), &(s->ttt), 
	    (su3_vector *)F_PT(s,dest) );
	scalar_mult_su3_vector( (su3_vector *)F_PT(s,dest), 1.0/(2.0*mass),
	    (su3_vector *)F_PT(s,dest) );
    }

    /* Polish off odd sites to correct for possible roundoff error */
    /* dest_o <- (M_adj M)^-1 temp_o  (odd sites only) */
    qic.parity = ODD;
    cgn = ks_congrad_site( temp, dest, &qic, mass, fn );

    // check_invert( dest, src, mass, 1e-6, fn);

    return cgn;
}

/* FOR TESTING: multiply src by matrix and check against dest */
void check_invert( field_offset src, field_offset dest, Real mass,
		   Real tol, ferm_links_t *fn){
    register int i,k,flag;
    register site *s;
    Real r_diff, i_diff;
    double sum,sum2,dflag,dmaxerr,derr;
    dslash_site( src, F_OFFSET(cg_p), EVENANDODD, fn);
    FORALLSITES(i,s){
	scalar_mult_add_su3_vector( &(s->cg_p), (su3_vector *)F_PT(s,src),
	    +2.0*mass, &(s->cg_p) );
    }
    sum2=sum=0.0;
    dmaxerr=0;
    flag = 0;
    FORALLSITES(i,s){
	for(k=0;k<3;k++){
	    r_diff = ((su3_vector *)F_PT(s,dest))->c[k].real
		    - s->cg_p.c[k].real;
	    i_diff = ((su3_vector *)F_PT(s,dest))->c[k].imag
		    - s->cg_p.c[k].imag;
	    if( fabs(r_diff) > tol || fabs(i_diff) > tol ){
	      printf("site %d color %d  expected ( %.4e , %.4e ) got ( %.4e , %.4e )\n",
		     i,k,
		     ((su3_vector *)F_PT(s,dest))->c[k].real,
		     ((su3_vector *)F_PT(s,dest))->c[k].imag,
		     s->cg_p.c[k].real,s->cg_p.c[k].imag);
	      flag++;
	    }
	    derr = r_diff*r_diff + i_diff*i_diff;
	    if(derr>dmaxerr)dmaxerr=derr;
 	    sum += derr;
	}
	sum2 += magsq_su3vec( (su3_vector *)F_PT(s,dest) );
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
}

/* FOR TESTING: multiply src by matrix and check against dest */
void check_invert_field( su3_vector *src, su3_vector *dest, Real mass,
			 Real tol, ferm_links_t *fn){
    register int i,k,flag;
    register site *s;
    Real r_diff, i_diff;
    double sum,sum2,dflag,dmaxerr,derr;
    su3_vector *tmp;

    tmp = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
    if(tmp==NULL){
      printf("check_invert_field(%d): no room for tmp\n",this_node);
      terminate(1);
    }

    /* Compute tmp = M src */
    ks_dirac_op( src, tmp, mass, fn);

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

/* FOR TESTING: multiply src by Madj M and check against dest */
void check_invert_field2( su3_vector *src, su3_vector *dest, Real mass,
			  Real tol, ferm_links_t *fn){
    register int i,k,flag;
    register site *s;
    Real r_diff, i_diff;
    double sum,sum2,dflag,dmaxerr,derr;
    su3_vector *tmp;

    tmp = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
    if(tmp==NULL){
      printf("check_invert_field(%d): no room for tmp\n",this_node);
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

