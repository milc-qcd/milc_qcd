/******** mat_invert.c *************/
/* MIMD version 7*/

/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.  Assumes "dslash" has
   been defined to be the appropriate "dslash_fn" or "dslash_eo"
*/

#include "generic_ks_includes.h"
#include "../include/dslash_ks_redefine.h"
#include "../include/openmp_defs.h"
#ifdef HAVE_QUDA
#include <quda_milc_interface.h>
#include "../include/generic_quda.h"
#endif

#ifdef CGTIME
static const char *prec_label[2] = {"F", "D"};
#endif

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

    switch(parity){
    case(EVEN): otherparity=ODD; break;
    case(ODD):  otherparity=EVEN; break;
    case(EVENANDODD): otherparity=EVENANDODD;
    }

    dslash_fn_field( src, dst, otherparity, fn);
    dslash_fn_field( dst, dst, parity, fn);
    FORSOMEFIELDPARITY_OMP(i,parity,){
      scalar_mult_add_su3_vector( dst+i, src+i, -msq_x4, dst+i );
    } END_LOOP_OMP;

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

static void deflate(su3_vector *dst, su3_vector *src, Real mass,
		    int Num, int parity, int zerodst){
  int i, j;
  double_complex *c;

  /* Remove the Num low eigenmodes from the trial solution, leaving
     the high eigenmode trial solution. Skip this step if dst is
     supposed to be zero */

  if(!zerodst)
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

/* dst = (M^adj M)^{-1} M^adj src */

int mat_invert_cg_field(su3_vector *src, su3_vector *dst, 
			quark_invert_control *qic,
			Real mass, imp_ferm_links_t *fn ){
    su3_vector *tmp;
    double dtime;

    node0_printf("Plain CG inversion from given initial guess for mass %f\n", mass);

    tmp = create_v_field();

    /* "Precondition" both even and odd sites */
    /* tmp <- M_adj * src */
    /* The following operation is done in the prevailing
       precision.  The algorithm needs to be fixed! */
    ks_dirac_adj_op( src, tmp, mass, EVENANDODD, fn);

    /* Do deflation if we have eigenvectors and the deflate parameter is true */
    /* Skip MILC CPU deflation if using QUDA deflation */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
    if(param.eigen_param.Nvecs > 0 && qic->deflate){

      dtime = - dclock();
#ifdef CGTIME
      node0_printf("deflating on even sites for mass %g with %d eigenvec\n",
		   mass, param.eigen_param.Nvecs);
#endif
      
      /* dst is replaced by the exact low mode solution plus the provided high mode */
      deflate(dst, tmp, mass, param.eigen_param.Nvecs, EVEN, 0);
      
      dtime += dclock();
#ifdef CGTIME
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
#endif
    }
#endif
      
    /* dst_e <- (M_adj M)^-1 tmp_e  (even sites only) */
    qic->parity = EVEN;
    int cgn = ks_congrad_field( tmp, dst, qic, mass, fn );
    int even_iters = qic->final_iters;

    /* Skip MILC CPU deflation if using QUDA deflation */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
    if(param.eigen_param.Nvecs > 0 && qic->deflate){

      dtime = - dclock();
#ifdef CGTIME
      node0_printf("deflating on odd sites for mass %g with %d eigenvec\n", mass, param.eigen_param.Nvecs);
#endif
      
      deflate(dst, tmp, mass, param.eigen_param.Nvecs, ODD, 0);
      
      dtime += dclock();
#ifdef CGTIME
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
#endif
    }
#endif

    /* dst_o <- (M_adj M)^-1 tmp_o  (odd sites only) */
    qic->parity = ODD;
    cgn += ks_congrad_field( tmp, dst, qic, mass, fn );
    qic->final_iters += even_iters;

    destroy_v_field(tmp);

    //    node0_printf("Entering check_invert_field in mat_invert_cg_field\n");
    //    fflush(stdout);
    //    check_invert_field( dst, src, mass, 1e-6, fn, EVENANDODD);
    
    return cgn;
}

/*****************************************************************************/
/* This algorithm solves even and odd sites separately with zero initial guess */
/* This method may reduce the error in low-mode components for a given residual tolerance */
/* It is also helpful when the source is purely even or purely odd  */

/* dst = M_adj * (MM^adj)^{-1} src */

int mat_invert_cgz_field(su3_vector *src, su3_vector *dst, 
			 quark_invert_control *qic,
			 Real mass, imp_ferm_links_t *fn ){
    su3_vector *tmp;
    double dtime;

    node0_printf("Plain CG inversion from zero initial guess for mass %f\n", mass);

    /* Create a null tmp as zero initial guess */
    tmp = create_v_field();

    /* Put "exact" low-mode even-site solution in tmp if deflate parameter is true */

    /* Skip MILC CPU deflation if using QUDA deflation */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
    if(param.eigen_param.Nvecs > 0 && qic->deflate){

      dtime = - dclock();
#ifdef CGTIME
      node0_printf("deflating on even sites for mass %g with %d eigenvec\n", mass, param.eigen_param.Nvecs);
#endif
      
      /* tmp is replaced by the exact low mode solution plus zero high mode */
     deflate(tmp, src, mass, param.eigen_param.Nvecs, EVEN, 1);
      
     dtime += dclock();
#ifdef CGTIME
     node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
#endif
    }
#endif
      
    /* Solve for all modes using tmp as an initial guess */
    /* tmp_e <- (M_adj M)^-1 src_e  (even sites only) */
    qic->parity = EVEN;
    int cgn = ks_congrad_field( src, tmp, qic, mass, fn );
    int even_iters = qic->final_iters;

    /* Put "exact" low-mode odd-site solution in tmp if deflate parameter is true */

    /* Skip MILC CPU deflation if using QUDA deflation */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
    if(param.eigen_param.Nvecs > 0 && qic->deflate){

      dtime = - dclock();
#ifdef CGTIME
      node0_printf("deflating on odd sites for mass %g with %d eigenvec\n",
		   mass, param.eigen_param.Nvecs);
#endif
      
      deflate(tmp, src, mass, param.eigen_param.Nvecs, ODD, 1);
      
      dtime += dclock();
#ifdef CGTIME
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
#endif
    }
#endif

    /* Solve for all modes using tmp as an initial guess */
    /* tmp_o <- (M_adj M)^-1 src_o  (odd sites only) */
    qic->parity = ODD;
    cgn += ks_congrad_field( src, tmp, qic, mass, fn );
    qic->final_iters += even_iters;

    /* Convert solution to normal equation into solution to M dst = src */
    /* dst <- M_adj * tmp */
    /* The following operation is done in the prevailing
       precision.  The algorithm needs to be fixed! */
    ks_dirac_adj_op( tmp, dst, mass, EVENANDODD, fn);

    destroy_v_field(tmp);

    //    node0_printf("Entering check_invert_field in mat_invert_cg_field\n");
    //    fflush(stdout);
    //    check_invert_field( dst, src, mass, 1e-6, fn, EVENANDODD);
    
    return cgn;
}

/*****************************************************************************/
/* Compute M^-1 * phi, answer in dest
  Uses phi, ttt, resid, xxx, and cg_p as workspace */
int mat_invert_cg( field_offset src, field_offset dest, field_offset temp,
		   Real mass, int prec, imp_ferm_links_t *fn ){
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

    int cgn = mat_invert_cg_field( tsrc, tdest, &qic, mass, fn );

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
    double dtime;

    node0_printf("UML inversion with mass %f\n", mass);
    /* "Precondition" both even and odd sites */
    /* temp <- M_adj * src */

    ks_dirac_adj_op( src, tmp, mass, EVENANDODD, fn );

#if EIGMODE != EIGCG
    /* Skip MILC CPU deflation if using QUDA deflation */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
    if(param.eigen_param.Nvecs > 0 && qic->deflate){

      dtime = - dclock();
      node0_printf("deflating on even sites for mass %g with %d eigenvec\n", mass, param.eigen_param.Nvecs);
      
      deflate(dst, tmp, mass, param.eigen_param.Nvecs, EVEN, 0);

      dtime += dclock();
#ifdef CGTIME
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
#endif
    }
#endif
#endif

    /* dst_e <- (M_adj M)^-1 tmp_e  (even sites only) */
    qic->parity     = EVEN;
#if EIGMODE == EIGCG
    cgn = ks_inc_eigCG_parity(tmp, dst, eigVal, eigVec, &param.eigcgp, qic, mass, fn);
    report_status(qic);
#else
    cgn = ks_congrad_field( tmp, dst, qic, mass, fn );
#endif
    int even_iters = qic->final_iters;

    /* reconstruct odd site solution */
    /* dst_o <-  1/2m (Dslash_oe*dst_e + src_o) */
    dslash_field( dst, ttt, ODD, fn );
    FORODDFIELDSITES_OMP(i,){
      sub_su3_vector( src+i, ttt+i, dst+i);
      scalar_mult_su3_vector( dst+i, 1.0/(2.0*mass), dst+i );
    } END_LOOP_OMP;

#if EIGMODE != EIGCG
    /* Skip MILC CPU deflation if using QUDA deflation */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
    if(param.eigen_param.Nvecs > 0 && qic->deflate){

      dtime = - dclock();
      node0_printf("deflating on odd sites for mass %g with %d eigenvec\n", mass, param.eigen_param.Nvecs);
      
      deflate(dst, tmp, mass, param.eigen_param.Nvecs, ODD, 0);
      
      dtime += dclock();
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
    }
#endif
#endif

    /* Polish off odd sites to correct for possible roundoff error */
    /* dst_o <- (M_adj M)^-1 temp_o  (odd sites only) */
    qic->parity = ODD;
    cgn += ks_congrad_field( tmp, dst, qic, mass, fn );
    qic->final_iters += even_iters;

    //    node0_printf("Entering check_invert_field in mat_invert_uml_field\n");
    //    fflush(stdout);
    //    check_invert_field( dst, src, mass, 1e-6, fn, EVENANDODD);
    destroy_v_field(tmp);
    destroy_v_field(ttt);

    return cgn;
}

#if defined(HAVE_QUDA) && defined(USE_CG_GPU)

/********************************************************************/
/* Multigrid solution of the full Dirac equation for both parities  */
/********************************************************************/

static void *mg_preconditioner = NULL;

int mat_invert_mg_field_gpu(su3_vector *t_src, su3_vector *t_dest, 
			    quark_invert_control *qic,
			    Real mass, imp_ferm_links_t *fn){

#ifdef MULTIGRID
  char myname[] = "mat_invert_mg_field_gpu";
  QudaInvertArgs_t inv_args;
  int i;
  double dtimec = -dclock();
#ifdef CGTIME
  double nflop = 1187;
#endif
  
  /* Initialize qic */
  qic->size_r = 0;
  qic->size_relr = 0;
  qic->final_iters   = 0;
  qic->final_restart = 0;
  qic->converged     = 1;
  qic->final_rsq = 0.;
  qic->final_relrsq = 0.;
  
  /* Compute source norm */
  double source_norm = 0.0;
  FORSOMEFIELDPARITY(i,qic->parity){
    source_norm += (double)magsq_su3vec( &t_src[i] );
  } END_LOOP;
  g_doublesum( &source_norm );
#ifdef CG_DEBUG
  node0_printf("mat_invert_mg_field_gpu: source_norm = %e\n", (double)source_norm);
#endif
  
  /* Provide for trivial solution */
  if(source_norm == 0.0){
    /* Zero the solution and return zero iterations */
    FORSOMEFIELDPARITY(i,qic->parity){
      memset(t_dest + i, 0, sizeof(su3_vector));
    } END_LOOP;
    
    dtimec += dclock();
#ifdef CGTIME
    if(this_node==0){
      printf("CONGRAD5: time = %e (fn_QUDA_MG %s) masses = 1 iters = %d mflops = %e\n",
	     dtimec, prec_label[qic->prec-1], qic->final_iters, 
	     (double)(nflop*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
      fflush(stdout);}
#endif
    
    return 0;
  }
  
  /* Initialize QUDA parameters */
  
  initialize_quda();

  // need to set a dummy value, ignored (for now),
  // MG will eventually support Schur solve
  inv_args.evenodd = QUDA_EVEN_PARITY; 
  /*if(qic->parity == EVEN){
          inv_args.evenodd = QUDA_EVEN_PARITY;
  }else if(qic->parity == ODD){
          inv_args.evenodd = QUDA_ODD_PARITY;
  }else{
    printf("%s: Unrecognised parity\n",myname);
    terminate(2);
  }*/

  inv_args.max_iter = qic->max*qic->nrestart;
#if defined(MAX_MIXED)
  inv_args.mixed_precision = 2;
#elif defined(HALF_MIXED)
  inv_args.mixed_precision = 1;
#else
  inv_args.mixed_precision = 0;
#endif
  
  su3_matrix* fatlink = get_fatlinks(fn);
  su3_matrix* longlink = get_lnglinks(fn);
  const int quda_precision = qic->prec;
  
  double residual, relative_residual;
  int num_iters = 0;
  
  inv_args.naik_epsilon = fn->eps_naik;
  
#if (FERM_ACTION==HISQ)
  inv_args.tadpole = 1.0;
#else
  inv_args.tadpole = u0;
#endif
  
  // for newer versions of QUDA we need to invalidate the gauge field if the links are new
  if ( fn != get_fn_last() || fresh_fn_links(fn) ){
    
    node0_printf("%s: fn, notify: Signal QUDA to refresh links\n", myname);
    cancel_quda_notification(fn);
    set_fn_last(fn);
    /* Hack to cause QUDA to copy new links to the GPU */
    num_iters = -1;

    node0_printf("%s: setting up the MG inverter\n", myname);

    /* Set up the MG inverter when the links change */
    /* FIXME: what do we do if the mgparamfile changes? */
    
    if (mg_preconditioner == NULL ){
      double mg_regen_time = -dclock();
      mg_preconditioner = qudaMultigridCreate(MILC_PRECISION,
                              quda_precision,
                              mass,
                              inv_args,
                              fatlink,
                              longlink,
                              qic->mgparamfile);

      mg_regen_time += dclock();
      node0_printf("%s: MG inverter setup complete. Time = %g\n", myname,
		   mg_regen_time);
    } else {
      node0_printf("%s: MG inverter already set up.  Skipping.\n", myname);
    }
  }

  /* Specify the type of rebuild _if_ a rebuild is required */
  int mg_rebuild_type = 1;
  if (qic->mg_rebuild_type == THINREBUILD) {
    mg_rebuild_type = 0;
  } 

  qudaInvertMG(MILC_PRECISION,
	       quda_precision, 
	       mass,
	       inv_args,
	       qic->resid,
	       qic->relresid,
	       fatlink, 
	       longlink,
	       mg_preconditioner,
	       mg_rebuild_type,
	       t_src, 
	       t_dest,
	       &residual,
	       &relative_residual, 
	       &num_iters);
  
  qic->final_rsq = residual*residual;
  qic->final_relrsq = relative_residual*relative_residual;
  qic->final_iters = num_iters;

  // check for convergence 
  qic->converged = (residual < qic->resid) ? 1 : 0;

  // Cumulative residual. Not used in practice 
  qic->size_r = 0.0;
  qic->size_relr = 0.0; 

  dtimec += dclock();

#ifdef CGTIME
  if(this_node==0){
    printf("CONGRAD5: time = %e (fn_QUDA_MG %s) masses = 1 iters = %d mflops = %e\n",
	   dtimec, prec_label[quda_precision-1], qic->final_iters, 
	   (double)(nflop*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

    //node0_printf("Entering check_invert_field in mat_invert_mg_field_gpu\n");
    //  fflush(stdout);
    //  check_invert_field( t_dest, t_src, mass, 1e-6, fn, EVENANDODD);

  report_status(qic);

  return num_iters;

#else
  node0_printf("mat_invert_mg_field_gpu: ERROR. Multigrid is available only with GPU compilation\n");
  terminate(1);
  return 0;  
#endif
}

void mat_invert_mg_cleanup(void){

#ifdef MULTIGRID
  if(mg_preconditioner != NULL)
    qudaMultigridDestroy(mg_preconditioner);
#endif
}

#endif /* HAVE_QUDA */

/*****************************************************************************/
/* This algorithm solves even and odd sites separately */

/* dst = (M^adj M)^{-1} M^adj src */

int mat_invert_block_cg(su3_vector **src, su3_vector **dst, 
			Real mass, int nsrc, quark_invert_control *qic,
			imp_ferm_links_t *fn ){
  double dtime;
  
  node0_printf("Plain CG inversion from zero initial guess for mass %f\n", mass);

  su3_vector **tmp = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(int is = 0; is < nsrc; is++){
    tmp[is] = create_v_field();
    
    /* "Precondition" both even and odd sites */
    /* tmp <- M_adj * src */
    /* The following operation is done in the prevailing
       precision.  The algorithm needs to be fixed! */
    ks_dirac_adj_op( src[is], tmp[is], mass, EVENANDODD, fn);
  }

  /* Put "exact" low-mode even-site solution in tmp if deflate parameter is true */
  /* Skip MILC CPU deflation if using QUDA deflation */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
  if(param.eigen_param.Nvecs > 0 && qic->deflate){

    dtime = -dclock();
#ifdef CGTIME
    node0_printf("deflating on even sites for mass %g with %d eigenvec\n",
		 mass, param.eigen_param.Nvecs);
#endif

    for(int is = 0; is < nsrc; is++)
      deflate(dst[is], tmp[is], mass, param.eigen_param.Nvecs, EVEN, 0);
    
    dtime += dclock();
#ifdef CGTIME
    node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
#endif
  }
#endif

  /* dst_e <- (M_adj M)^-1 tmp_e  (even sites only) */
  qic->parity = EVEN;
  int cgn = ks_congrad_block_field( nsrc, tmp, dst, qic, mass, fn );
  int even_iters = qic->final_iters;

  /* Deflation on odd sites */
  /* Skip MILC CPU deflation if using QUDA deflation */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
  if(param.eigen_param.Nvecs > 0 && qic->deflate){
    
    dtime = -dclock();
#ifdef CGTIME
    node0_printf("deflating on odd sites for mass %g with %d eigenvec\n", mass, param.eigen_param.Nvecs);
#endif
      
    for(int is = 0; is < nsrc; is++)
      deflate(dst[is], tmp[is], mass, param.eigen_param.Nvecs, ODD, 0);
    
    dtime += dclock();
#ifdef CGTIME
    node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
#endif
  }
#endif

  /* dst_o <- (M_adj M)^-1 tmp_o  (odd sites only) */
  qic->parity = ODD;
  cgn += ks_congrad_block_field( nsrc, tmp, dst, qic, mass, fn );
  qic->final_iters += even_iters;

  for(int is = 0; is < nsrc; is++)
    destroy_v_field(tmp[is]);
  free(tmp);

  return cgn;
}
/*****************************************************************************/
/* This algorithm solves even and odd sites separately with zero initial guess */

/* dst = M_adj * (MM^adj)^{-1} src */

int mat_invert_block_cgz(su3_vector **src, su3_vector **dst, 
			 Real mass, int nsrc, quark_invert_control *qic,
			 imp_ferm_links_t *fn ){
  double dtime;
  
  node0_printf("Plain CG inversion from zero initial guess for mass %f\n", mass);

  /* Create tmps as zero initial guesses */
  su3_vector **tmp = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(int is = 0; is < nsrc; is++)
    tmp[is] = create_v_field();

  /* Put "exact" low-mode even-site solution in tmp if deflate parameter is true */
  /* Skip MILC CPU deflation if using QUDA deflation */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
  if(param.eigen_param.Nvecs > 0 && qic->deflate){
    for(int is = 0; is < nsrc; is++){

      dtime = -dclock();
#ifdef CGTIME
      node0_printf("deflating on even sites for mass %g with %d eigenvec\n",
		   mass, param.eigen_param.Nvecs);
#endif
      
      /* tmp[is] is replaced by the exact low mode solution plus zero high mode */
      deflate(tmp[is], src[is], mass, param.eigen_param.Nvecs, EVEN, 1);
      
      dtime += dclock();
#ifdef CGTIME
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
#endif
    }
  }
#endif

  /* Solve for all modes on even sites using tmp as an initial guess */
  /* tmp_e <- (M_adj M)^-1 src_e  (even sites only) */
  qic->parity = EVEN;
  node0_printf("Calling ks_congrad_block_field\n"); fflush(stdout);
  int cgn = ks_congrad_block_field( nsrc, src, tmp, qic, mass, fn );
  int even_iters = qic->final_iters;

  /* Put "exact" low-mode odd-site solution in tmp if deflate parameter is true */
  /* Skip MILC CPU deflation if using QUDA deflation */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
  if(param.eigen_param.Nvecs > 0 && qic->deflate){
    for(int is = 0; is < nsrc; is++){

      dtime = -dclock();
#ifdef CGTIME
      node0_printf("deflating on odd sites for mass %g with %d eigenvec\n",
		   mass, param.eigen_param.Nvecs);
#endif
      
      deflate(tmp[is], src[is], mass, param.eigen_param.Nvecs, ODD, 1);
      
      dtime += dclock();
#ifdef CGTIME
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
#endif
    }
  }
#endif

  /* Solve for all modes on odd sites using tmp as an initial guess */
  /* dst_o <- (M_adj M)^-1 tmp_o  (odd sites only) */
  qic->parity = ODD;
  cgn += ks_congrad_block_field( nsrc, src, tmp, qic, mass, fn );
  qic->final_iters += even_iters;

  /* Convert solution to normal equation into solution to M dst = src */
  /* dst <- M_adj * tmp */
  /* The following operation is done in the prevailing
     precision.  The algorithm needs to be fixed! */
  for(int is = 0; is < nsrc; is++)
    ks_dirac_adj_op( tmp[is], dst[is], mass, EVENANDODD, fn);

  for(int is = 0; is < nsrc; is++)
    destroy_v_field(tmp[is]);
  free(tmp);

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
  
  su3_vector *ttt = create_v_field();
  double dtime;
  int i;
  
  su3_vector **tmp = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(int is = 0; is < nsrc; is++)
    tmp[is] = create_v_field();
  
  /* "Precondition" both even and odd sites */
  /* tmp <- M_adj * src */
  
  for(int is = 0; is < nsrc; is++){
    ks_dirac_adj_op( src[is], tmp[is], mass, EVENANDODD, fn );

    /* Skip MILC CPU deflation if using QUDA deflation */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
    if(param.eigen_param.Nvecs > 0 && qic->deflate){
      dtime = - dclock();
#ifdef CGTIME
      node0_printf("deflating on even sites for mass %g with %d eigenvec\n",
		   mass, param.eigen_param.Nvecs);
#endif
      
      /* dst[is] is replaced by the exact low mode solution plus zero high mode */
      deflate(dst[is], tmp[is], mass, param.eigen_param.Nvecs, EVEN, 0);
      
      dtime += dclock();
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
    }
#endif
  }
  
  /* dst_e <- (M_adj M)^-1 tmp_e  (even sites only) */
  qic->parity     = EVEN;
  int cgn = ks_congrad_block_field(nsrc, tmp, dst, qic, mass, fn );
  int even_iters = qic->final_iters;
  
  /* reconstruct odd site solution */
  /* dst_o <-  1/2m (Dslash_oe*dst_e + src_o) */
  for(int is = 0; is < nsrc; is++){
    dslash_field( dst[is], ttt, ODD, fn );
    FORODDFIELDSITES_OMP(i,){
      sub_su3_vector( src[is]+i, ttt+i, dst[is]+i);
      scalar_mult_su3_vector( dst[is]+i, 1.0/(2.0*mass), dst[is]+i );
    } END_LOOP_OMP;

    /* Skip MILC CPU deflation if using QUDA deflation */
#if !( defined(USE_CG_GPU) && defined(HAVE_QUDA) )
    if(param.eigen_param.Nvecs > 0 && qic->deflate){
      dtime = - dclock();
      node0_printf("deflating on odd sites for mass %g with %d eigenvec\n",
		   mass, param.eigen_param.Nvecs);
      
      deflate(dst[is], tmp[is], mass, param.eigen_param.Nvecs, ODD, 0);
      
      dtime += dclock();
      node0_printf("Time to deflate %d modes %g\n", param.eigen_param.Nvecs, dtime);
    }
#endif
  }
  
  /* Polish off odd sites to correct for possible roundoff error */
  /* dst_o <- (M_adj M)^-1 temp_o  (odd sites only) */
  qic->parity = ODD;
  cgn += ks_congrad_block_field(nsrc, tmp, dst, qic, mass, fn );
  qic->final_iters += even_iters;
  
  //    check_invert_field( dst, src, mass, 1e-6, fn, EVENANDODD);

  for(int is = 0; is < nsrc; is++)
    destroy_v_field(tmp[is]);
  free(tmp);
  destroy_v_field(ttt);
  
  return cgn;
}

#if defined(HAVE_QUDA) && defined(MULTIGRID)
/*****************************************************************************/
/* This algorithm solves the Dirac equation for both parities using
   staggered multigrid */

int mat_invert_block_mg(su3_vector **src, su3_vector **dst, 
			Real mass, int nsrc, quark_invert_control *qic,
			imp_ferm_links_t *fn){
  
  //int cgn = 0;
  
  /* Temporary until there is multi-rhs support for multigrid */
  //for(int is = 0; is < nsrc; is++) {
  //  cgn += mat_invert_mg_field_gpu(src[is], dst[is], qic, mass, fn );
  //}
  //return cgn;

#ifdef MULTIGRID
  char myname[] = "mat_invert_block_mg";
  QudaInvertArgs_t inv_args;
  int i;
  double dtimec = -dclock();
#ifdef CGTIME
  double nflop = 1187;
#endif

  /* Initialize qic */
  qic->size_r = 0;
  qic->size_relr = 0;
  qic->final_iters   = 0;
  qic->final_restart = 0;
  qic->converged     = 1;
  qic->final_rsq = 0.;
  qic->final_relrsq = 0.;

  
  /* Initialize QUDA parameters */
  initialize_quda();

  // need to set a dummy value, ignored (for now),
  // MG will eventually support Schur solve
  inv_args.evenodd = QUDA_EVEN_PARITY; 
  /*if(qic->parity == EVEN){
          inv_args.evenodd = QUDA_EVEN_PARITY;
  }else if(qic->parity == ODD){
          inv_args.evenodd = QUDA_ODD_PARITY;
  }else{
    printf("%s: Unrecognised parity\n",myname);
    terminate(2);
  }*/

  inv_args.max_iter = qic->max * qic->nrestart;
#if defined(MAX_MIXED)
  inv_args.mixed_precision = 2;
#elif defined(HALF_MIXED)
  inv_args.mixed_precision = 1;
#else
  inv_args.mixed_precision = 0;
#endif

  su3_matrix* fatlink = get_fatlinks(fn);
  su3_matrix* longlink = get_lnglinks(fn);
  const int quda_precision = qic->prec;
  
  double residual, relative_residual;
  int num_iters = 0;
  
  inv_args.naik_epsilon = fn->eps_naik;

#if (FERM_ACTION==HISQ)
  inv_args.tadpole = 1.0;
#else
  inv_args.tadpole = u0;
#endif

  // for newer versions of QUDA we need to invalidate the gauge field if the links are new
  if ( fn != get_fn_last() || fresh_fn_links(fn) ){

    node0_printf("%s: fn, notify: Signal QUDA to refresh links\n", myname);
    cancel_quda_notification(fn);
    set_fn_last(fn);
    /* Hack to cause QUDA to copy new links to the GPU */
    num_iters = -1;

    node0_printf("%s: setting up the MG inverter\n", myname);

    /* Set up the MG inverter when the links change */
    /* FIXME: what do we do if the mgparamfile changes? */

    if (mg_preconditioner == NULL ){
      double mg_regen_time = -dclock();
      mg_preconditioner = qudaMultigridCreate(MILC_PRECISION,
                              quda_precision,
                              mass,
                              inv_args,
                              fatlink,
                              longlink,
                              qic->mgparamfile);

      mg_regen_time += dclock();
      node0_printf("%s: MG inverter setup complete. Time = %g\n", myname, mg_regen_time);
    } else {
      node0_printf("%s: MG inverter already set up.  Skipping.\n", myname);
    }
  }

  /* Specify the type of rebuild _if_ a rebuild is required */
  int mg_rebuild_type = 1;
  if (qic->mg_rebuild_type == THINREBUILD) {
    mg_rebuild_type = 0;
  }

  int cgn = 0;

  /* check norms */
  for (int is = 0; is < nsrc; is++) {

    node0_printf("%s: starting source %d\n", myname, is);

    /* Compute source norm */
    double source_norm = 0.0;
    FORSOMEFIELDPARITY(i,qic->parity){
      source_norm += (double)magsq_su3vec( &src[is][i] );
    } END_LOOP;
    g_doublesum( &source_norm );
#ifdef CG_DEBUG
    node0_printf("%s: source %d source_norm = %e\n", myname, is, (double)source_norm);
#endif

    /* Provide for trivial solution */
    if (source_norm == 0.0) {
      /* Zero the solution and do not update the number of iterations */
      FORSOMEFIELDPARITY(i, qic->parity){
        memset(dst[is] + i, 0, sizeof(su3_vector));
      } END_LOOP;
    }
  }

  qudaInvertMsrcMG(MILC_PRECISION,
        quda_precision, 
        mass,
        inv_args,
        qic->resid,
        qic->relresid,
        fatlink, 
        longlink,
        mg_preconditioner,
        mg_rebuild_type,
        (void**)src,
        (void**)dst,
        &residual,
        &relative_residual, 
        &num_iters,
        nsrc);

  qic->final_rsq = residual*residual;
  qic->final_relrsq = relative_residual*relative_residual;
  qic->final_iters = num_iters;

  // check for convergence 
  qic->converged = (residual < qic->resid) ? 1 : 0;

  // Cumulative residual. Not used in practice 
  qic->size_r = 0.0;
  qic->size_relr = 0.0;

  dtimec += dclock();

#ifdef CGTIME
  if(this_node==0){
    printf("CONGRAD5: time = %e (fn_QUDA %s) masses = 1 srcs = %d iters = %d mflops = %e\n",
           dtimec, prec_label[quda_precision-1], nsrc, qic->final_iters,
           (double)(nflop*nsrc*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);
  }
#endif

    //node0_printf("Entering check_invert_field in mat_invert_mg_field_gpu\n");
    //  fflush(stdout);
    //  check_invert_field( t_dest, t_src, mass, 1e-6, fn, EVENANDODD);

  report_status(qic);

  return num_iters;

#else
  node0_printf("%s: ERROR. Multigrid is available only with GPU compilation\n", myname);
  terminate(1);
  return 0;  
#endif

}
#endif

/*****************************************************************************/
int mat_invert_uml(field_offset src, field_offset dest, field_offset temp,
		   Real mass, int prec, imp_ferm_links_t *fn ){
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

    int cgn = mat_invert_uml_field(tsrc, tdest, &qic, mass, fn);

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
/*****************************************************************************/

int mat_invert_field(su3_vector *src, su3_vector *dst, 
		     quark_invert_control *qic,
		     Real mass, imp_ferm_links_t *fn){

  int cgn = 0;

  switch(qic->inv_type){
  case CGTYPE:

    cgn = mat_invert_cg_field(src, dst, qic, mass, fn );
    break;

  case CGZTYPE:

    cgn = mat_invert_cgz_field(src, dst, qic, mass, fn );
    break;

  case UMLTYPE:

    cgn = mat_invert_uml_field(src, dst, qic, mass, fn );
    break;

  case MGTYPE:
    if(qic->mg_rebuild_type == CGREBUILD){    
      cgn = mat_invert_uml_field(src, dst, qic, mass, fn );
      node0_printf("WARNING: Best practices for inv_type MG is to move forced CG solves to a different set\n");
#ifdef HAVE_QUDA
      /* Force a reload b/c of sloppy link precision changes */
      refresh_fn_links(fn);
#endif
    } else {
      /* inv_type == MGTYPE */
#if defined(USE_CG_GPU) && defined(HAVE_QUDA)
      /* Currently only available through QUDA on GPUs */
      cgn = mat_invert_mg_field_gpu(src, dst, qic, mass, fn );
#else
      node0_printf("mat_invert_field: ERROR. Multigrid is available only with QUDA compilation\n");
      terminate(1);
#endif
    }
    break;

  default:

    node0_printf("mat_invert_field: Bad inv_type %d\n", qic->inv_type);
    terminate(1);
  }
  return cgn;
}

/*****************************************************************************/
/* Generic multi-rhs inversion                                               */
/*****************************************************************************/

int mat_invert_block(su3_vector **src, su3_vector **dst, 
		     Real mass, int nsrc, quark_invert_control *qic,
		     imp_ferm_links_t *fn){
  int cgn = 0;

  switch(qic->inv_type){

  case CGTYPE:
    
    cgn = mat_invert_block_cg(src, dst, mass, nsrc, qic, fn);
    break;

  case CGZTYPE:

    cgn = mat_invert_block_cgz(src, dst, mass, nsrc, qic, fn);
    break;

  case UMLTYPE:

    cgn = mat_invert_block_uml(src, dst, mass, nsrc, qic, fn);
    break;

  case MGTYPE:

    if(qic->mg_rebuild_type == CGREBUILD){    
      cgn = mat_invert_block_uml(src, dst, mass, nsrc, qic, fn);
      node0_printf("WARNING: Best practices for inv_type MG is to move forced CG solves to a different set\n");
#if defined(USE_CG_GPU) && defined(HAVE_QUDA)
      /* Force a reload b/c of sloppy link precision changes */
      refresh_fn_links(fn);
#endif
    } else {
      /* inv_type == MGTYPE */
#if defined(MULTIGRID) && defined(HAVE_QUDA)
      /* Currently only available through QUDA on GPUs */
      cgn = mat_invert_block_mg(src, dst, mass, nsrc, qic, fn);
#else
      node0_printf("mat_invert_block: ERROR. Multigrid is available only with GPU compilation\n");
      terminate(1);
#endif
    }
    break;

  default:

    node0_printf("mat_invert_block: Bad inv_type %d\n", qic->inv_type);
    terminate(1);
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
    ks_dirac_op( src, tmp, mass, parity, fn);

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

/* DEBUG */
static Real 
my_relative_residue(su3_vector *p, su3_vector *q, int parity)
{
  double residue, num, den;
  int i;
  
  residue = 0;
  FORSOMEFIELDPARITY_OMP(i,parity,private(num,den) reduction(+:residue)){
    num = (double)magsq_su3vec( &(p[i]) );
    den = (double)magsq_su3vec( &(q[i]) );
    residue += (den==0) ? 1.0 : (num/den);
  } END_LOOP_OMP;

  g_doublesum(&residue);

  if(parity == EVENANDODD)
    return sqrt(residue/volume);
  else
    return sqrt(2*residue/volume);
}

/*****************************************************************************/
/* FOR TESTING: multiply src by Madj M and check against dest */
void check_invert_field2( su3_vector *src, su3_vector *dest, Real mass,
			  Real tol, imp_ferm_links_t *fn, int parity){
    register int i,k,flag;
    register site *s;
    Real r_diff, i_diff;
    double sum,sum2,dflag,dmaxerr,derr;

    su3_vector *tmp   = create_v_field();
    su3_vector *resid = create_v_field();
    if(tmp==NULL || resid==NULL){
      printf("check_invert_field2(%d): no room for tmp and resid\n",this_node);
      terminate(1);
    }

    /* Compute tmp = -(Madj M) dest */
    ks_dirac_opsq( dest, tmp, mass, parity, fn);

    FORSOMEFIELDPARITY(i,parity){
      add_su3_vector(src+i, tmp+i, resid+i);
    } END_LOOP;

    node0_printf("check_invert_field2: Checking solution with tolerance %e\n",
		 tol);
    sum2=sum=0.0;
    dmaxerr=0;
    flag = 0;

    FORSOMEFIELDPARITY(i,parity){
      sum2 += magsq_su3vec( src+i );
    } END_LOOP;

    g_doublesum( &sum2 );
    double srcnorm = sqrt(sum2);

    FORSOMEFIELDPARITY(i,parity){
      for(k=0;k<3;k++){
	r_diff = resid[i].c[k].real/srcnorm;
	i_diff = resid[i].c[k].imag/srcnorm;
	if( fabs(r_diff) > tol || fabs(i_diff) > tol ){
	  if(flag < 50) /* Don't print too many */
	    printf("site %d color %d  expected ( %.4e , %.4e ) got ( %.4e , %.4e )\n",
		   i,k,
		   src[i].c[k].real, src[i].c[k].imag,
		   -tmp[i].c[k].real, -tmp[i].c[k].imag);
	  flag++;
	}
	derr = r_diff*r_diff + i_diff*i_diff;
	if(derr>dmaxerr)dmaxerr=derr;
	sum += derr;
      }
    } END_LOOP;

    g_doublesum( &sum );
    dflag=flag;
    g_doublesum( &dflag );
    g_doublemax( &dmaxerr );
    double rel_resid = my_relative_residue(resid, dest, parity);

    if(this_node==0){
      printf("Inversion checked, resid = %e rel_resid = %g\n",
	     sqrt(sum), rel_resid);
      printf("Flagged comparisons = %d\n",(int)dflag);
      printf("Max err. = %e\n",sqrt(dmaxerr));
      fflush(stdout);
    }
		 
    destroy_v_field(tmp);
    destroy_v_field(resid);
}

#if 0 /* Haven't been using these */
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

#endif
