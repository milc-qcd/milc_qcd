/******* d_cgilu_cl_lean.c - CG-ILU for  clover fermions ****/
/* MIMD version 6 */

/* Modifications:
   7/18/01 Uses dslash_w_special - CD
   1/24/00 combined with Schroedinger functional version - UMH
   4/26/98 Moved parameters to structures CD
   8/29/97 ANSI prototyping and added comments C. D.
   Code created by U.M.H.
   */

/* Memory stingy version
   "r" overwrites src on even sites
   "p" overwrites src on odd sites
   3/29/00 EVENFIRST is the rule now. CD.
   */

/* Requires qic wilson vector temporaries wv1 and wv2 */

/* ------------------------------------------------------------
   The matrix to be inverted is written in block even-odd form as


   M = ( R_o     -K D_oe )
       ( -K D_eo  R_e    )

   where R_o and R_e are 1 - K Clov_c/u_0^3 i sigma_mu,nu F_mu,nu
   are the site-diagonal hermitian clover matrices on odd and even
   sites.

   The ILU decomposition results in M = L A U


   M = ( 1            0 ) ( R_o  0   ) ( 1  -K/R_o D_oe )
       ( -K D_eo/R_o  1 ) ( 0    M_e ) ( 0     1        )

   where

         M_e = R_e - K^2 D_eo/R_o D_oe

   acts only on even sites.  Since M_e is not hermitian, it is
   necessary to invert M_e_dag M_e where M_e_dag = M_e^\dagger

 ------------------------------------------------------------ */

#include "generic_clover_includes.h"

/*#define CGTIME*/       /* Uncomment if you want timing info */

int cgilu_cl(            /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    )
{
  /* Unpack required members of the structures */
  int MinCG = qic->min;      /* minimum number of iterations */
  int MaxCG = qic->max;      /* maximum number of iterations */
  Real RsdCG = qic->resid;  /* desired residual - 
				 normalized as sqrt(r*r)/sqrt(src_e*src_e */
  int flag = qic->start_flag;   /* 0: use a zero initial guess; 1: use dest */
  field_offset tmp = qic->wv1;   /* size of wilson_vector */
  field_offset my_mp = qic->wv2;    /* size of wilson_vector */

  dirac_clover_param *dcp 
    = (dirac_clover_param *)dmp; /* Cast pass-through pointer */

  field_offset f_mn = dcp->work_f_mn;  /* size of su3_matrix */
  Real Kappa = dcp->Kappa;     /* hopping */
  Real Clov_c = dcp->Clov_c;   /* Perturbative clover coeff */
  Real U0 = dcp->U0;           /* Tadpole correction to Clov_c */
  /* End of unpacking required members of structures */

  int N_iter;
  register int i;
  register site *s;
  Real size_src;
  double sr, cp, c, d;
  register Real a, b;
  register Real MKsq = -Kappa*Kappa;
  register field_offset r,p;
  Real CKU0 = Kappa*Clov_c/(U0*U0*U0);
#ifdef CGTIME
  double dtime;
#endif
  msg_tag *tage[8],*tago[8];
  int is_startedo, is_startede;
  
  is_startedo = is_startede = 0;

    if(even_sites_on_node!=odd_sites_on_node){
      printf("Need same number of even and odd sites on each node\n");
      terminate(1);
    }
  
  /* Compute R_e and R_o */
  make_clov(CKU0,f_mn);

  /* Invert R_o only, leaving R_e on even sites and 1/R_o on odd sites 
     in "clov" and "clov_diag" */
  make_clovinv(ODD);
  
  /* "r" is the residual vector, which is a pointer to src since the
     source is overwritten to save space.
     "p" and "my_mp" are working vectors for the conjugate gradient.  */
  r = src;
  /* This disgusting trick makes p for each even site actually be
     src on some corresponding odd site */
  p = src + even_sites_on_node*sizeof(site);

  /** if(this_node==0)printf("CGILU: p=%d\n",p); **/
  
  /* CG_ILU: */
  
  /* Start Inversion */
  
#ifdef CGTIME
  dtime = -dclock();
#endif

  /* ---------  src = L^(-1)*src  ------------- */

  /* ( src_o ) = ( 1           0 ) ( src_o )
     ( sec_e )   (-K D_eo/R_o  1 ) ( src_e )

     */
  
  /* mp_o = 1/R_o srce_e */
  mult_ldu(src, my_mp, ODD);
  /* mp_e = D_eo/R_o srce_e */
  dslash_w_special(my_mp, my_mp, PLUS, EVEN, tage, is_startede);
  is_startede = 1;
  
  /* src_e = srce_e + K D_eo/R_o srce_e */
  /* (leaving src_o = src_o)   */
  sr=0.0;
  FOREVENSITESDOMAIN(i,s) {
    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,src), 
			  (wilson_vector *)F_PT(s,my_mp), Kappa, 
			  (wilson_vector *)F_PT(s,src) );
    sr += (double)magsq_wvec( (wilson_vector *)F_PT(s,src) );
  }

  /* --------- size_src = sqrt(src_e*src_e) ---------- */
  g_doublesum(&sr);
  size_src = (Real)sqrt(sr);
  
/*  if(this_node==0)printf("beginning inversion--size_src=%e\n",
			 (double)size_src); */
  
  /* r and p overwrite the source (bad for
     dynamical fermions but good for quenched calcs) */
  
  /* --------- dest_o = src_o ---------- */
  /* Initial guess */
  /* set dest = src on odd sites, whatever else you do
     (even if we restart with a nonzero solution vector, the end of the
     subroutine rebuilds the odd component from the even one. The (trivial)
     solution of the odd component of the equation is dest = src, before
     we rotate back to the basis in which  M is not checkerboard-diagonal) */
  FORODDSITESDOMAIN(i,s) {
    copy_wvec( (wilson_vector *)F_PT(s,src),
	       (wilson_vector *)F_PT(s,dest) );
  }
  
  /* --------- if flag == 0 set dest_e = 0 ---------- */
  /* (then  r = src which is a pointer equivalence so nothing else to do) */
  if(flag == 0) {
/*    if(this_node==0)printf("dest_0=0\n"); */
    FOREVENSITESDOMAIN(i,s) {
      clear_wvec( (wilson_vector *)F_PT(s,dest) );
    }
  }
  /* ---------- otherwise, use the given dest_e --------- */
  /* --------- and compute r_e = src_e[1]-M_e*dest_e --------- */
  if(flag != 0) {
/*    if(this_node==0)    printf("dest_0  !=0\n"); */
    /* tmp_e = R_e dest_e */
    mult_ldu(dest, tmp, EVEN);
    /* mp_o = D_oe dest_e */
    dslash_w_special(dest, my_mp, PLUS, ODD, tago, is_startedo);
    is_startedo = 1;
    /* tmp_o = 1/R_o D_oe dest_e */
    mult_ldu(my_mp, tmp, ODD);
    /* mp_e = D_eo/R_o D_oe dest_e */
    dslash_w_special(tmp, my_mp, PLUS, EVEN, tage, is_startede);
    is_startede = 1;
    /* mp_e = R_e dest_e - K^2 D_eo/R_o D_oe dest_e = M_e dest_e */
    /* r_e = src_e - M_e dest_e */
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,tmp), 
			    (wilson_vector *)F_PT(s,my_mp), MKsq, 
			    (wilson_vector *)F_PT(s,my_mp) );
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,r),
			    (wilson_vector *)F_PT(s,my_mp), -1.0, 
			    (wilson_vector *)F_PT(s,r) );
    }
  }

  /* --------- size_r = sqrt(r * r)/size_src --------- */
  sr=0.0;
  FOREVENSITESDOMAIN(i,s) {
    sr += (double)magsq_wvec( (wilson_vector *)F_PT(s,r) );
  }
  g_doublesum(&sr);
  qic->size_r = (Real)sqrt(sr)/size_src;
  

  /* --------- p_e = M_e_dag*r_e --------- */
  mult_ldu(r, tmp, EVEN);
  dslash_w_special(r, my_mp, MINUS, ODD, tago, is_startedo);
  is_startedo = 1;
  mult_ldu(my_mp, tmp, ODD);
  dslash_w_special(tmp, p, MINUS, EVEN, tage, is_startede);
  is_startede = 1;
  
  /* --------- cp = |p|^2 --------- */
  cp=0.0;
  FOREVENSITESDOMAIN(i,s) {
    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,tmp), 
			  (wilson_vector *)F_PT(s,p), MKsq, 
			  (wilson_vector *)F_PT(s,p) );
    cp += (double)magsq_wvec( (wilson_vector *)(F_PT(s,p)) );
  }
  g_doublesum(&cp);

  /* ---------------------------------------------- */
  /* --------- Beginning of cg iterations --------- */

  for( N_iter = 0; N_iter < MinCG || (N_iter < MaxCG && RsdCG  < qic->size_r); 
       N_iter = N_iter + 1) {

    /* --------- c = cp -------- */
    c=cp;

    /* ---------  mp_e = M_e*p_e --------- */
    /* tmp_e = R_e p_e */
    mult_ldu(p, tmp, EVEN);
    /* mp_o = D_oe p_e */
    dslash_w_special(p, my_mp, PLUS, ODD, tago, is_startedo);
    is_startedo = 1;
    /* tmp_o = 1/R_o D_oe p_e */
    mult_ldu(my_mp, tmp, ODD);
    /* mp_e = D_eo/R_o D_oe p_e */
    dslash_w_special(tmp, my_mp, PLUS, EVEN, tage, is_startede);
    is_startede = 1;
    
    /* mp_e = R_e p_e - K^2 D_eo/R_o D_oe p_e */
    d=0.0;
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,tmp), 
			    (wilson_vector *)F_PT(s,my_mp), MKsq, 
			    (wilson_vector *)F_PT(s,my_mp) );
      d += (double)magsq_wvec( (wilson_vector *)F_PT(s,my_mp) );
    }
    /* --------- a = c/|mp_e|^2 --------- */
    g_doublesum(&d);
    a = (Real)(c/d);
    
    /* --------- dest_e = dest_e + a*p_e --------- */
    /* --------- r_e = r_e - a*mp_e --------- */
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest),
			    (wilson_vector *)F_PT(s,p), a,
			    (wilson_vector *)F_PT(s,dest) );
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,r),
			    (wilson_vector *)F_PT(s,my_mp), -a, 
			    (wilson_vector *)F_PT(s,r) );
    }
    
    /* --------- mp_e M_e_dag*r_e --------- */
    
    mult_ldu(r, tmp, EVEN);
    dslash_w_special(r, my_mp, MINUS, ODD, tago, is_startedo);
    is_startedo = 1;
    mult_ldu(my_mp, tmp, ODD);
    dslash_w_special(tmp, my_mp, MINUS, EVEN, tage, is_startede);
    is_startede = 1;

    cp=0.0;
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,tmp), 
			    (wilson_vector *)F_PT(s,my_mp), MKsq, 
			    (wilson_vector *)F_PT(s,my_mp) );
      cp += (double)magsq_wvec( (wilson_vector *)F_PT(s,my_mp) );
    }
    /* ---------  cp = |mp|^2 --------- */
    g_doublesum(&cp);
    /* ---------  b = |mp|^2/c --------- */
    b = (Real)(cp/c);

    /* ---------  p_e = mp_e + b*p_e --------- */
    sr=0.0;
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,my_mp), 
			    ((wilson_vector *)F_PT(s,p)),
			    b, ((wilson_vector *)F_PT(s,p)) );
      sr += (double)magsq_wvec( ((wilson_vector *)F_PT(s,r)) );
    }
    g_doublesum(&sr);
    
    /* --------- size_r = sqrt(r_e * r_e)/size_src --------- */
    qic->size_r = (Real)sqrt(sr)/size_src;
    /**if(this_node==0)printf("iteration= %d, residue= %e\n",N_iter,
      (double)(qic->size_r));**/
  }
  /* ------------------------------------ */
  /* --------- End of iterations --------- */

#ifdef CGTIME
  dtime += dclock();
#endif
  if(this_node==0){
    if(N_iter==0)
      printf("CGILU: NO ITERATIONS TAKEN size_r= %.2e\n",qic->size_r);
#ifdef CGTIME
    else
      printf("CGILU: time = %.2e size_r= %.2e iters= %d MF = %.1f\n",
	     dtime,qic->size_r,N_iter,
	     (double)8064*N_iter*even_sites_on_node/(dtime*1e6));
#endif
    fflush(stdout);
  }
  
  /** if( (qic->size_r) > RsdCG ) {
    if(this_node==0)printf(" CG_ILU: Not Converged\n");
    } **/
  
  /* --------- dest_o = U^(-1)_oo/R_o dest_o + U^(-1)_oe dest_e --------- */
  /* mp_o = D_oe * dest_e */
  dslash_w_special(dest, my_mp, PLUS, ODD, tago, is_startedo);
  is_startedo = 1;
  /* mp_o = dest_o + K D_oe * dest_e (remember dest_o = original src_o still)*/
  FORODDSITESDOMAIN(i,s) {
    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest), 
			  (wilson_vector *)F_PT(s,my_mp), Kappa, 
			  (wilson_vector *)F_PT(s,my_mp) );
  }
  /* dest_o = 1/R_o dest_o + K/R_o D_oe * dest_e */
  mult_ldu(my_mp, dest, ODD);

  for( i=XUP; i <= TUP; i++) {
    if(is_startedo)cleanup_gather(tago[i]);
    if(is_startedo)cleanup_gather(tago[OPP_DIR(i)]);
    if(is_startede)cleanup_gather(tage[i]);
    if(is_startede)cleanup_gather(tage[OPP_DIR(i)]);
  }
  is_startede = is_startedo = 0;
    
  free_clov();
  return(N_iter);
} /* d_cgilu_cl_lean.c */
