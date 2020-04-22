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
#include "../include/openmp_defs.h"

/*****************************************************************************/
/* dst = M src. With parity selection */

void ks_dirac_op( su3_vector *src, su3_vector *dst, Real mass, 
		  int parity, imp_ferm_links_t *fn){
    register int i;
    register site *s;

    dslash_field( src, dst, parity, fn);
    FORSOMEPARITYDOMAIN_OMP(i,s,parity,){
      scalar_mult_add_su3_vector( dst+i, src+i, +2.0*mass, dst+i);
    } END_LOOP_OMP;
}

/*****************************************************************************/
/* dst = Madj src. With parity selection */

void ks_dirac_adj_op( su3_vector *src, su3_vector *dst, Real mass,
		      int parity, imp_ferm_links_t *fn){
    register int i;
    register site *s;

    dslash_field( src, dst, parity, fn);
    FORSOMEPARITYDOMAIN_OMP(i,s,parity,){
      scalar_mult_su3_vector( dst+i, -1.0, dst+i);
      scalar_mult_add_su3_vector( dst+i, src+i, 2.0*mass, dst+i);
    } END_LOOP_OMP;
}

/*****************************************************************************/
/* dst = Madj dst. With parity selection */

void ks_dirac_adj_op_inplace( su3_vector *dst, Real mass,
			      int parity, imp_ferm_links_t *fn){
    register int i;
    register site *s;
    su3_vector *tvec = create_v_field();

    dslash_field( dst, tvec, parity, fn);
    FORSOMEFIELDPARITY_OMP(i,parity,){
      scalar_mult_su3_vector( dst+i, 2.0*mass, dst+i);
      scalar_mult_add_su3_vector( dst+i, tvec+i, -1.0, dst+i);
    } END_LOOP_OMP;
    destroy_v_field(tvec);
}

/*****************************************************************************/
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
    FORSOMEPARITYDOMAIN_OMP(i,s,parity,){
      scalar_mult_su3_vector( dst+i, -1.0, dst+i);
      scalar_mult_sum_su3_vector( dst+i, src+i, msq_x4 );
    } END_LOOP_OMP;

    free(tmp);
}

/*****************************************************************************/
/* Returns the dot product of two fermion vectors */
static void dot_product(su3_vector *vec1, su3_vector *vec2, 
			double_complex *dot, int parity) {
  register double re,im ;
  register  int i;
  
  re=im=0.0;
  FORSOMEFIELDPARITY_OMP(i,parity,reduction(+:re,im)){
    complex cc = su3_dot( &(vec1[i]), &(vec2[i]) );
    re += cc.real ;
    im += cc.imag ;
  } END_LOOP_OMP;
  dot->real = re ; 
  dot->imag = im ;
  g_dcomplexsum(dot);
}

/*****************************************************************************/
/* Returns vec2 = vec2 - cc*vec1   cc is a double complex   */
static void complex_vec_mult_sub(double_complex *cc, su3_vector *vec1, 
			  su3_vector *vec2, int parity){

  register  int i;
  complex sc ;
  
  sc.real= (Real)(cc->real) ; 
  sc.imag= (Real)(cc->imag) ;

  FORSOMEFIELDPARITY_OMP(i,parity,){
    c_scalar_mult_sub_su3vec(&(vec2[i]), (&sc), &(vec1[i])) ;
  } END_LOOP_OMP;
}

/*****************************************************************************/
/*  Projects out the eigVec from the vec. Num is the Number of vectors  *
 *  The vectors are assumed to be orthonormal.                             */
static void project_out(su3_vector *vec, int Num, int parity){
  register int i ;
  double_complex cc ;
  double ptime = -dclock();

  for(i=Num-1;i>-1;i--){
    dot_product(eigVec[i], vec, &cc, parity) ;
    complex_vec_mult_sub(&cc, eigVec[i], vec, parity);
  }

  ptime += dclock();
  node0_printf("Time to project out %d modes %g sec\n", Num, ptime);
}

/*****************************************************************************/
/* Construct the exact solution in the space spanned by eigVec 
   keep the trial solution in the complementary space */

static void deflate(su3_vector *dst, su3_vector *src, Real mass, int Num, int parity){
  int i, j;
  double_complex *c;

  /* Remove the Num low eigenmodes from the trial solution, leaving the
     high eigenmode trial solution */

  project_out(dst, Num, parity);

  /* Then add the exact solution back */
  /* dst_eo <- sum_j ((eigVec_eo[j].src_eo)/(eigVal[j]+4*mass*mass)) eigVec_eo[j] */

  c = (double_complex *)malloc(Num*sizeof(double_complex));
  for( j = 0; j < Num; j++){
    double re = 0.0;
    double im = 0.0;
    FORSOMEFIELDPARITY_OMP(i,parity,reduction(+:re,im)){
      complex cc = su3_dot( eigVec[j]+i, src+i );
      re += cc.real;
      im += cc.imag;
    } END_LOOP_OMP;
    c[j] = dcmplx(re, im);
  }
  g_vecdcomplexsum( c, Num );
  for( j = 0; j < Num; j++){
    CDIVREAL( c[j], eigVal[j]+4.0*mass*mass, c[j] );
    complex ctmp = cmplx(c[j].real, c[j].imag);
    FORSOMEFIELDPARITY_OMP(i,parity,){
      c_scalar_mult_add_su3vec( dst+i, &ctmp, eigVec[j]+i );
    } END_LOOP_OMP;
  }
  free(c);
}

/*****************************************************************************/
/* This algorithm solves even and odd sites separately */

int mat_invert_cg_field(su3_vector *src, su3_vector *dst, 
			quark_invert_control *qic,
			Real mass, imp_ferm_links_t *fn ){
    int cgn;
    su3_vector *tmp;
    double dtime;

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

    /* Do deflation if we have eigenvectors and the deflate parameter is true */

    if(param.eigen_param.Nvecs > 0 && qic->deflate){

      dtime = - dclock();
#ifdef CGTIME
      node0_printf("deflating on even sites for mass %g with %d eigenvec\n", mass, param.eigen_param.Nvecs);
#endif
      
      deflate(dst, tmp, mass, param.eigen_param.Nvecs, EVEN);
      
      dtime += dclock();
#ifdef CGTIME
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
#endif
    }
      
    /* dst_e <- (M_adj M)^-1 temp_e  (even sites only) */
    qic->parity = EVEN;
    cgn = ks_congrad_field( tmp, dst, qic, mass, fn );

    if(param.eigen_param.Nvecs > 0 && qic->deflate){

      dtime = - dclock();
#ifdef CGTIME
      node0_printf("deflating on odd sites for mass %g with %d eigenvec\n", mass, param.eigen_param.Nvecs);
#endif
      
      deflate(dst, tmp, mass, param.eigen_param.Nvecs, ODD);
      
      dtime += dclock();
#ifdef CGTIME
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
#endif
    }

    /* dst_o <- (M_adj M)^-1 temp_o  (odd sites only) */
    qic->parity = ODD;
    cgn += ks_congrad_field( tmp, dst, qic, mass, fn );

    free(tmp);

    //    check_invert_field2( dst, tmp, mass, 1e-6, fn, EVENANDODD);
    //    check_invert_field( dst, src, mass, 1e-6, fn, EVENANDODD);
    return cgn;
}

/*****************************************************************************/
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
    qic.nsrc = 1;
    qic.resid      = sqrt(rsqprop);
    qic.relresid   = 0;

    su3_vector *tsrc = create_v_field_from_site_member(src);
    su3_vector *tdest = create_v_field_from_site_member(dest);

    cgn = mat_invert_cg_field( tsrc, tdest, &qic, mass, fn );

    FORALLSITES_OMP(i,s,){
      su3vec_copy(tdest+i, (su3_vector *)F_PT(s,dest));
    } END_LOOP_OMP;

    destroy_v_field(tdest);
    destroy_v_field(tsrc);

    return cgn;
}

/*****************************************************************************/
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

where  A = (4m^2+D_eo D_eo^adj)^-1 and B = (4m^2+D_eo^adj D_eo)^-1

    Note: -D_oe = D_eo^adj and  B D_eo^adj = D_eo^adj A

*/
         
/*****************************************************************************/
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
    double dtime;

    /* "Precondition" both even and odd sites */
    /* temp <- M_adj * src */

    ks_dirac_adj_op( src, tmp, mass, EVENANDODD, fn );

#if EIGMODE != EIGCG
    if(param.eigen_param.Nvecs > 0 && qic->deflate){

      dtime = - dclock();
      node0_printf("deflating on even sites for mass %g with %d eigenvec\n", mass, param.eigen_param.Nvecs);
      
      deflate(dst, tmp, mass, param.eigen_param.Nvecs, EVEN);

      dtime += dclock();
#ifdef CGTIME
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
#endif
    }
#endif

    /* dst_e <- (M_adj M)^-1 tmp_e  (even sites only) */
    qic->parity     = EVEN;
#if EIGMODE == EIGCG
    cgn = ks_inc_eigCG_parity(tmp, dst, eigVal, eigVec, &param.eigcgp, qic, mass, fn);
    report_status(qic);
#else
    cgn = ks_congrad_field( tmp, dst, qic, mass, fn );
#endif
    even_iters = qic->final_iters;

    /* reconstruct odd site solution */
    /* dst_o <-  1/2m (Dslash_oe*dst_e + src_o) */
    dslash_field( dst, ttt, ODD, fn );
    FORODDFIELDSITES_OMP(i,){
      sub_su3_vector( src+i, ttt+i, dst+i);
      scalar_mult_su3_vector( dst+i, 1.0/(2.0*mass), dst+i );
    } END_LOOP_OMP;

#if EIGMODE != EIGCG
    if(param.eigen_param.Nvecs > 0 && qic->deflate){

      dtime = - dclock();
      node0_printf("deflating on odd sites for mass %g with %d eigenvec\n", mass, param.eigen_param.Nvecs);
      
      deflate(dst, tmp, mass, param.eigen_param.Nvecs, ODD);
      
      dtime += dclock();
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
    }
#endif

    /* Polish off odd sites to correct for possible roundoff error */
    /* dst_o <- (M_adj M)^-1 temp_o  (odd sites only) */
    qic->parity = ODD;
    cgn += ks_congrad_field( tmp, dst, qic, mass, fn );
    qic->final_iters += even_iters;

    //    check_invert_field( dst, src, mass, 1e-6, fn, EVENANDODD);
    destroy_v_field(tmp);
    destroy_v_field(ttt);

    return cgn;
}

/*****************************************************************************/
/* This algorithm solves even sites, reconstructs odd and then polishes
   to compensate for loss of significance in the reconstruction
   BLOCKCG version
*/

int mat_invert_block_uml(su3_vector **src, su3_vector **dst, 
			 Real mass, int nsrc, quark_invert_control *qic,
			 imp_ferm_links_t *fn){
  
  int cgn;
  register int i, is;
  register site *s;
  su3_vector **tmp = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  su3_vector *ttt = create_v_field();
  double dtime;
  int even_iters;
  
  for(is = 0; is < nsrc; is++)
    tmp[is] = create_v_field();
  
  /* "Precondition" both even and odd sites */
  /* tmp <- M_adj * src */
  
  for(is = 0; is < nsrc; is++){
    ks_dirac_adj_op( src[is], tmp[is], mass, EVENANDODD, fn );
    
    if(param.eigen_param.Nvecs > 0 && qic->deflate){
      dtime = - dclock();
      node0_printf("deflating on even sites for mass %g with %d eigenvec\n",
		   mass, param.eigen_param.Nvecs);
      
      deflate(dst[is], tmp[is], mass, param.eigen_param.Nvecs, EVEN);
      
      dtime += dclock();
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
    }
  }
  
  /* dst_e <- (M_adj M)^-1 tmp_e  (even sites only) */
  qic->parity     = EVEN;
  cgn = ks_congrad_block_field(nsrc, tmp, dst, qic, mass, fn );
  even_iters = qic->final_iters;
  
  /* reconstruct odd site solution */
  /* dst_o <-  1/2m (Dslash_oe*dst_e + src_o) */
  for(is = 0; is < nsrc; is++){
    dslash_field( dst[is], ttt, ODD, fn );
    FORODDFIELDSITES_OMP(i,){
      sub_su3_vector( src[is]+i, ttt+i, dst[is]+i);
      scalar_mult_su3_vector( dst[is]+i, 1.0/(2.0*mass), dst[is]+i );
    } END_LOOP_OMP;

    if(param.eigen_param.Nvecs > 0 && qic->deflate){
      dtime = - dclock();
      node0_printf("deflating on odd sites for mass %g with %d eigenvec\n",
		   mass, param.eigen_param.Nvecs);
      
      deflate(dst[is], tmp[is], mass, param.eigen_param.Nvecs, ODD);
      
      dtime += dclock();
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
    }
  }
  
  /* Polish off odd sites to correct for possible roundoff error */
  /* dst_o <- (M_adj M)^-1 temp_o  (odd sites only) */
  qic->parity = ODD;
  cgn += ks_congrad_block_field(nsrc, tmp, dst, qic, mass, fn );
  qic->final_iters += even_iters;
  
  //    check_invert_field( dst, src, mass, 1e-6, fn, EVENANDODD);
  for(is = 0; is < nsrc; is++)
    destroy_v_field(tmp[is]);
  free(tmp);
  destroy_v_field(ttt);
  
  return cgn;
}

/*****************************************************************************/
/* This algorithm solves the Dirac equation for both parities using
   staggered multigrid */

int mat_invert_block_mg(su3_vector **src, su3_vector **dst, 
			Real mass, int nsrc, quark_invert_control *qic,
			imp_ferm_links_t *fn){
  
  int cgn = 0;
  register int is;
  
  /* Temporary until there is multi-rhs support for multigrid */
  for(is = 0; is < nsrc; is++)
    cgn += mat_invert_mg_field_gpu(src[is], dst[is], qic, mass, fn );
  
  return cgn;
}

/*****************************************************************************/
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
    qic.nsrc       = 1;
    qic.resid      = sqrt(rsqprop);
    qic.relresid   = 0;

    cgn = mat_invert_uml_field(tsrc, tdest, &qic, mass, fn);

    FORALLSITES_OMP(i,s,){
      su3vec_copy(tdest+i, (su3_vector *)F_PT(s,dest));
    } END_LOOP_OMP;

    // check_invert( dest, src, mass, 1e-6, fn);

    destroy_v_field(tdest);
    destroy_v_field(tsrc);
    return cgn;
}

/*****************************************************************************/
/* Generic inverter entry point for a single mass and source */

int mat_invert_field(su3_vector *src, su3_vector *dst, 
		     quark_invert_control *qic,
		     Real mass, imp_ferm_links_t *fn, int use_precond){

  int cgn = 0;

  if(qic->inv_type == CGTYPE){
    if(use_precond)
      /* Preconditioned inversion */
      cgn = mat_invert_uml_field(src, dst, qic, mass, fn );
    else
      /* Unpreconditioned inversion */
      cgn = mat_invert_cg_field(src, dst, qic, mass, fn );
  } else {
    /* inv_type == MGTYPE */
#ifdef USE_CG_GPU
    /* Currently only available through QUDA on GPUs */
    cgn = mat_invert_mg_field_gpu(src, dst, qic, mass, fn );
#else
    node0_printf("mat_invert_field: ERROR. Multigrid is available only with GPU compilation\n");
    terminate(1);
#endif
  }
  return cgn;
}

/*****************************************************************************/
/* Generic multi-rhs inversion                                               */
/*****************************************************************************/

int mat_invert_block(su3_vector **src, su3_vector **dst, 
		     Real mass, int nsrc, quark_invert_control *qic,
		     imp_ferm_links_t *fn){
  int cgn;
  if(qic->inv_type == CGTYPE){
    cgn = mat_invert_block_uml(src, dst, mass, nsrc, qic, fn);
  } else {
    /* inv_type == MGTYPE */
#ifdef USE_CG_GPU
    /* Currently only available through QUDA on GPUs */
    cgn = mat_invert_block_mg(src, dst, mass, nsrc, qic, fn);
#else
    node0_printf("mat_invert_block: ERROR. Multigrid is available only with GPU compilation\n");
    terminate(1);
#endif
  }
  return cgn;
}
  
/*****************************************************************************/
/* FOR TESTING: multiply src by matrix and check against dest */
void check_invert_field( su3_vector *src, su3_vector *dest, Real mass,
			 Real tol, imp_ferm_links_t *fn, int parity){
    register int i,k,flag;
    register site *s;
    Real r_diff, i_diff;
    double sum,sum2,dflag,dmaxerr,derr;
    su3_vector *tmp;

    tmp = create_v_field();

    /* Compute tmp = M src */
    node0_printf("check_invert_field: calling ks_dirac_op\n"); fflush(stdout);
    ks_dirac_op( src, tmp, mass, parity, fn);

    sum2=sum=0.0;
    dmaxerr=0;
    flag = 0;
    node0_printf("check_invert_field: checking diffs\n"); fflush(stdout);
    FORSOMEFIELDPARITY(i,parity){
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
	if(flag > 200)break;  // Don't write too many lines for debugging
    } END_LOOP;
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

/*****************************************************************************/
/* FOR TESTING: multiply src by matrix and check against dest */
void check_invert( field_offset src, field_offset dest, Real mass,
		   Real tol, imp_ferm_links_t *fn){

  su3_vector *tsrc = create_v_field_from_site_member(src);
  su3_vector *tdest = create_v_field_from_site_member(dest);
  
  check_invert_field( tsrc, tdest, mass, tol, fn, EVENANDODD );
  
  destroy_v_field(tdest);
  destroy_v_field(tsrc);
}

/*****************************************************************************/
/* FOR TESTING: multiply src by Madj M and check against dest */
void check_invert_field2( su3_vector *src, su3_vector *dest, Real mass,
			  Real tol, imp_ferm_links_t *fn, int parity){
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
    ks_dirac_opsq( src, tmp, mass, parity, fn);

    sum2=sum=0.0;
    dmaxerr=0;
    flag = 0;
    FORSOMEFIELDPARITY(i,parity){
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
    } END_LOOP;
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

/*****************************************************************************/
/* Creates an array of vectors for the block-cg solver */

static su3_vector **create_su3_vector_array(int n){
  su3_vector **a;
  int i;

  a = (su3_vector **)malloc(n*sizeof(su3_vector *));
  if(a == NULL){
    printf("f_meas: No room for array\n");
    terminate(1);
  }
  for(i = 0; i < n; i++) a[i] = create_v_field();
  return a;
}

/*****************************************************************************/
/* Destroys an array of vectors */

static void destroy_su3_vector_array(su3_vector **a, int n){
  int i;

  if(a == NULL)return;
  for(i = 0; i < n; i++)
    if(a[i] != NULL)
      destroy_v_field(a[i]);
}
