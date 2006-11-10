/******* d_mrilu_cl_lean.c - MR-ILU for  clover fermions ****/
/* MIMD version 7 */

/* Modifications:
   1/24/00 adapted from Schroedinger functional version, and
           moved parameters to structures - UMH
   2/12/97 ANSI prototyping U.M.H.
   */

/* Memory stingy version
  "r" overwrites src on even sites
   3/29/00 EVENFIRST is the rule now. CD.
*/

/* Requires qic wilson vector temporaries wv1, wv2, wv3, wv4 */

/* ------------------------------------------------------------
   The matrix to be inverted is written in block even-odd form as


   M = ( R_o     -K D_oe )
       ( -K D_eo  R_e    )

   where R_o and R_e are 1 - K Clov_c/u_0^3 i sigma_mu,nu F_mu,nu
   are the site-diagonal hermitian clover matrices on odd and even
   sites.

   The dslash operators D_oe src and D_eo src are given by

   SUM_dirs ( 
      ( 1 + gamma[dir] ) * U(x,dir) * src(x+dir)
    + ( 1 - gamma[dir] ) * U_adj(x-dir,dir) * src(x-dir)
   )

   with gammas defined in libraries/mb_gamma.c.

   The ILU decomposition results in M = L A U


   M = ( 1            0 ) ( R_o  0   ) ( 1  -K/R_o D_oe )
       ( -K D_eo/R_o  1 ) ( 0    M_e ) ( 0     1        )

   where

         M_e = R_e - K^2 D_eo/R_o D_oe

   acts only on even sites.  Since M_e is not hermitian, it is
   necessary to invert M_e_dag M_e where M_e_dag = M_e^\dagger

 ------------------------------------------------------------ */

/* The source vector is in "src", and the initial guess and answer
   in "dest".  "r" is the residual vector, which is a pointer to src since
   the source is overwritten to save space,
   and "my_mp" is a
   working vector for the conjugate gradient. 
   MinCG = minimum number of iterations.
   MaxCG = maximum number of iterations.
   size_r = desired residual, quit when we reach it.
   (Square root def for residue size_r = sqrt(r*r))
   ILU resides on parity=EVEN so do only even sites
*/

#include "generic_clover_includes.h"

int mrilu_cl(		/* Return value is number of iterations taken */
    field_offset src,	/* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,	/* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic,	/* parameters controlling inversion */
    void *dmp		/* parameters defining the Dirac matrix */
    )
{
  /* Unpack required members of the structures */
  int MinCG = qic->min;		/* minimum number of iterations */
  int MaxCG = qic->max;		/* maximum number of iterations */
  Real RsdCG = qic->resid;	/* desired residual -
                                 normalized as sqrt(r*r)/sqrt(src_e*src_e */
  int flag = qic->start_flag;	/* 0: use a zero initial guess; 1: use dest */
  field_offset tmp = qic->wv1;	/* size of wilson_vector */
  field_offset my_mp = qic->wv2;	/* size of wilson_vector */

  dirac_clover_param *dcp
    = (dirac_clover_param *)dmp;	/* Cast pass-through pointer */

  Real Kappa = dcp->Kappa;	/* hopping */
  Real Clov_c = dcp->Clov_c;	/* Perturbative clover coeff */
  Real U0 = dcp->U0;		/* Tadpole correction to Clov_c */
  /* End of unpacking required members of structures */

  int N_iter;
  register int i;
  register site *s;
  Real size_src;
  double c;
  complex a, ctmp;
  double_complex d;
  register Real MKsq = -Kappa*Kappa;
  register field_offset r;
  Real CKU0 = Kappa*Clov_c/(U0*U0*U0);
#ifdef CGTIME
  double dtime;
#endif
    if(even_sites_on_node!=odd_sites_on_node){
      printf("Need same number of even and odd sites on each node\n");
      terminate(1);
    }

  make_clov(CKU0);

  /* Take the inverse on the odd sublattice */
  make_clovinv(ODD);

  r = src;

  /* MR_ILU: */

  /* Start Inversion */

  /* src = L^(-1)*src */
#ifdef CGTIME
  dtime = -dclock();
#endif

  mult_ldu_site(src, my_mp, ODD);
  dslash_w_site(my_mp, my_mp, PLUS, EVEN);

  /* Normalisation  */
  c = 0.0;
#ifdef SCHROED_FUN
  FOREVENSITES(i,s) if(s->t > 0) {
#else
  FOREVENSITES(i,s) {
#endif
    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,src),
			  (wilson_vector *)F_PT(s,my_mp),
			  Kappa, (wilson_vector *)F_PT(s,src) );
    c += (double)magsq_wvec( (wilson_vector *)F_PT(s,src) );
  }
  g_doublesum(&c);
  size_src = (Real)sqrt((double)c);

/*  if(this_node==0)printf("beginning inversion--size_src=%e\n",
	(double)size_src); */

  /* r overwrites the source (bad for
     dynamical fermions but good for quenched calcs) */
  /* set r = src --- nothing to do */

  /* Initial guess */
  /* set dest = src on odd sites, whatever else you do
     (even if we restart with a nonzero solution vector, the end of the
     subroutine rebuilds the odd component from the even one. The (trivial)
     solution of the odd component of the equation is dest = src, before
     we rotate back to the basis in which  M is not checkerboard-diagonal) */
#ifdef SCHROED_FUN
  FORODDSITES(i,s) if(s->t > 0) {
#else
  FORODDSITES(i,s) {
#endif
    copy_wvec( (wilson_vector *)F_PT(s,src),
	       (wilson_vector *)F_PT(s,dest) );
  }


  /* code if you want to start with dest=0... */
  if(flag == 0) {
    /*  if(this_node==0){printf("dest_0=0\n");fflush(stdout);} */
#ifdef SCHROED_FUN
    FOREVENSITES(i,s) if(s->t > 0) {
#else
    FOREVENSITES(i,s) {
#endif
      clear_wvec( (wilson_vector *)F_PT(s,dest) );
    }
  }
  /* code if you want to start dest with some particular starting value... */
  /* r=src[1]-[L^(-1)*M*U^(-1)]*dest */
  if(flag != 0) {
    /*  if(this_node==0){printf("dest_0  !=0\n");fflush(stdout);} */
    /* we use my_mp temporarily to construct r */
    mult_ldu_site(dest, tmp, EVEN);
    dslash_w_site(dest, my_mp, PLUS, ODD);
    mult_ldu_site(my_mp, tmp, ODD);
    dslash_w_site(tmp, my_mp, PLUS, EVEN);
#ifdef SCHROED_FUN
    FOREVENSITES(i,s) if(s->t > 0) {
#else
    FOREVENSITES(i,s) {
#endif
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,tmp),
			    (wilson_vector *)F_PT(s,my_mp), MKsq,
			    (wilson_vector *)F_PT(s,my_mp) );
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,r),
			    (wilson_vector *)F_PT(s,my_mp), -1.0,
			    (wilson_vector *)F_PT(s,r) );
    }
  }

  c = 0.0;
#ifdef SCHROED_FUN
  FOREVENSITES(i,s) if(s->t > 0) {
#else
  FOREVENSITES(i,s) {
#endif
    c += (double)magsq_wvec( (wilson_vector *)F_PT(s,r) );
  }
  g_doublesum(&c);
  qic->size_r = (Real)sqrt(c)/size_src;
  /* if(this_node==0)    printf("beginning inversion--size_r=%e\n",
    (double)(qic->size_r)); */

  for( N_iter = 0; N_iter < MinCG || (N_iter < MaxCG && RsdCG  < qic->size_r);
      N_iter = N_iter + 1) {

    /*   my_mp = M(u)*r */
    mult_ldu_site(r, tmp, EVEN);
    dslash_w_site(r, my_mp, PLUS, ODD);
    mult_ldu_site(my_mp, tmp, ODD);
    dslash_w_site(tmp, my_mp, PLUS, EVEN);

    /* d = <Mr|r>  */
    /* c = <Mr|Mr>  */
    d = dcmplx((double)0.0,(double)0.0);
    c = 0.0;
#ifdef SCHROED_FUN
    FOREVENSITES(i,s) if(s->t > 0) {
#else
    FOREVENSITES(i,s) {
#endif
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,tmp),
			    (wilson_vector *)F_PT(s,my_mp), MKsq,
			    (wilson_vector *)F_PT(s,my_mp) );
      ctmp = wvec_dot( (wilson_vector *)F_PT(s,my_mp),
		       (wilson_vector *)F_PT(s,r) );
      CSUM(d, ctmp);
      c += (double)magsq_wvec( (wilson_vector *)F_PT(s,my_mp) );
    }
    g_dcomplexsum(&d);
    g_doublesum(&c);
    ctmp=cmplx(c,0.0);
    CDIV(d,ctmp,a);

    /* Note: here we could multiply a with an overrelaxation parameter! */

    /* dest = dest + a*r  */
#ifdef SCHROED_FUN
    FOREVENSITES(i,s) if(s->t > 0) {
#else
    FOREVENSITES(i,s) {
#endif
      c_scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest),
			      (wilson_vector *)F_PT(s,r), &a,
			      (wilson_vector *)F_PT(s,dest) );
    }

    /* r = r - a*mp */
    a.real = -a.real;
    a.imag = -a.imag;
    c = 0.0;
#ifdef SCHROED_FUN
    FOREVENSITES(i,s) if(s->t > 0) {
#else
    FOREVENSITES(i,s) {
#endif
      c_scalar_mult_add_wvec( (wilson_vector *)F_PT(s,r),
			      (wilson_vector *)F_PT(s,my_mp),
			      &a, (wilson_vector *)F_PT(s,r) );
      c += (double)magsq_wvec( (wilson_vector *)F_PT(s,r) );
    }
    g_doublesum(&c);
    qic->size_r = (Real)sqrt(c)/size_src;
    /**if(this_node==0)printf("iteration= %d, residue= %e\n",N_iter,
	    (double)(qic->size_r));**/
  }
#ifdef CGTIME
  dtime = -dclock();
#endif
  if(this_node==0){
    if(N_iter==0)
      printf("MRILU: NO ITERATIONS TAKEN size_r= %.2e\n",qic->size_r);
#ifdef CGTIME
    else
      printf("MRILU: time = %.2e size_r= %.2e iters= %d MF = %.1f\n",
	     dtime,qic->size_r,N_iter,
	     (double)4371*N_iter*even_sites_on_node/(dtime*1e6));
#else
    else
      printf("MRILU: size_r= %.2e iters= %d\n", qic->size_r, N_iter);
#endif
    fflush(stdout);
  }

  /**  if( (qic->size_r) > RsdCG ) {
    if(this_node==0)printf(" MR_ILU_ Not Converged: size_r=%e\n",
	    (double)(qic->size_r));
  } **/

  /* dest = R^(-1)*dest  */
  dslash_w_site(dest, my_mp, PLUS, ODD);
#ifdef SCHROED_FUN
  FORODDSITES(i,s) if(s->t > 0) {
#else
  FORODDSITES(i,s) {
#endif
    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest),
			  (wilson_vector *)F_PT(s,my_mp), Kappa,
			  (wilson_vector *)F_PT(s,my_mp) );
  }
  mult_ldu_site(my_mp, dest, ODD);

  free_clov();
  return(N_iter);
}

