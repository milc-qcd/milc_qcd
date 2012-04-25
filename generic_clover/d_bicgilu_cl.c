/******* d_bicgilu_cl.c - BiCGstab-ILU for  clover fermions ****/
/* MIMD version 7 */

/* Modifications:
   7/18/01 calls dslash_special CD
   1/24/00 combined with Schroedinger functional version - UMH
   4/26/98 Moved parameters to structures CD
   8/29/97 ANSI prototyping and added comments C. D.
   Code created by U.M.H.
   */

/* Memory stingy version
   3/29/00 EVENFIRST is the rule now. CD.
   */

#include "generic_clover_includes.h"

/*#define CG_DEBUG*/

/* The source vector is in "src", and the initial guess and answer in
   "dest". "p" and "mp" are working vectors for the conjugate
   gradient.  
   MaxCG = maximum number of iterations.  

   size_r = desired residual and 
   size_relr = desired relative residual, quit when we reach it.
   (Square def for residue size_r = r*r) 
   */

int bicgilu_cl_field_cpu(    /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    )
{
  /* Unpack required members of the structures */
  int max_restarts = qic->nrestart;  /* Number of restarts */
  int nrestart = 0;
  int restart = qic->max;            /* Restart interval */
  int MaxCG = restart*qic->max;     /* maximum number of iterations */
  Real RsdCG = qic->resid * qic->resid;      /* desired residual - 
				normalized as (r*r)/(src_e*src_e) */
  Real RRsdCG = qic->relresid * qic->relresid;  /* desired relative residual - */
  int flag = qic->start_flag;        /* 0: use a zero initial guess; 
					1: use dest */
  dirac_clover_param *dcp 
    = (dirac_clover_param *)dmp; /* Cast pass-through pointer */

  Real Kappa = dcp->Kappa;     /* hopping */
  Real Clov_c = dcp->Clov_c;   /* Perturbative clover coeff */
  Real U0 = dcp->U0;           /* Tadpole correction to Clov_c */
  /* End of unpacking required members of structures */

  wilson_vector *tmp=NULL,*mp=NULL,*rv=NULL,*sss=NULL,*r=NULL,
    *p=NULL,*ttt=NULL;
  int N_iter;
  register int i;
  register site *s;
  Real size_src, size_src2;
  double rsq, tsq;
  complex ctmp, a, b;
  complex omega, omegam;
  double_complex tdots,rvro,rvv,rvr;
  register Real MKsq = -Kappa*Kappa;
  Real CKU0 = Kappa*Clov_c/(U0*U0*U0);
#ifdef CGTIME
  double dtime;
#endif
  msg_tag *tago[8],*tage[8];
  int is_startedo, is_startede;
  
  is_startedo = is_startede = 0;

  qic->size_r = 0;
  qic->size_relr = 0;
  qic->final_rsq = 0;
  qic->final_relrsq = 0;
  qic->final_iters = 0;
  qic->final_restart = 0;

  //  if(even_sites_on_node!=odd_sites_on_node){
  //    printf("Need same number of even and odd sites on each node\n");
  //    terminate(1);
  //  }
  
  /* Compute R_e and R_o and put in "clov" and "clov_diag" */
  compute_clov(gen_clov, CKU0);

  /* Invert R_o only, leaving R_e on even sites and 1/R_o on odd sites 
     in "clov" and "clov_diag" */
  compute_clovinv(gen_clov, ODD);
  
  /* now we can allocate temporary variables and copy them */
  /* PAD may be used to avoid cache trashing */
#define PAD 0
    
  tmp    = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  mp     = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  rv     = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  sss    = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  r      = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  p      = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  ttt    = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));

  if(tmp == NULL || mp == NULL || mp == NULL || rv == NULL || 
     sss == NULL || r == NULL ){
    printf("bicgilu_cl_field(%d): No room for temporaries\n",this_node);
    terminate(1);
  }

  /* BiCGstab_ILU: */
  
  /* Transform source - result is in r (even) - and dest (odd) */
  
#ifdef CGTIME
  dtime = -dclock();
#endif
  
  /* now we copy src to temporaries */
  FORALLSITES(i,s) {
    r[i] = src[i];
  }

  /* src = L^(-1)*src */
  size_src = ilu_xfm_source(dest, r, mp, Kappa, &is_startede, tage);
  size_src2 = size_src*size_src;

#if 0

  mult_this_ldu_field(gen_clov, r, mp, ODD);
  dslash_w_field_special(mp, mp, PLUS, EVEN, tage, is_startede);
  is_startede = 1;
  
  /* Normalization  */
  rsq = 0.0;
  FOREVENSITESDOMAIN(i,s) {
    scalar_mult_add_wvec( &(r[i]), &(mp[i]), Kappa, &(r[i]) );
    rsq += (double)magsq_wvec( &(r[i]) );
    /* Save transformed source: Overwrite src on even sites */
    copy_wvec( &(r[i]), &(src[i]) );
  }
  g_doublesum(&rsq);
  size_src2 = (Real)rsq;
#endif

  /* Provision for trivial solution */
  if(size_src2 == 0.0){
    clear_wv_field(dest);
    free(tmp); free(mp); free(rv); free(sss); free(r); free(p); free(ttt);
    is_startede = is_startedo = 0;

    cleanup_tmp_links();
    cleanup_dslash_wtemps();
    return 0;
  }


  /* Save transformed source: Overwrite src on even sites */
  FOREVENSITESDOMAIN(i,s) {
    copy_wvec( &(r[i]), &(src[i]) );
  }

#ifdef CG_DEBUG
  node0_printf("beginning inversion--size_src=%e\n",
	       (double)size_src2); fflush(stdout);
#endif
  
  /* Initial guess */
  /* set dest = src on odd sites, whatever else you do
     (even if we restart with a nonzero solution vector, the end of the
     subroutine rebuilds the odd component from the even one. The (trivial)
     solution of the odd component of the equation is dest = src, before
     we rotate back to the basis in which  M is not checkerboard-diagonal) */
  FORODDSITESDOMAIN(i,s) {
    copy_wvec( &(r[i]), &(dest[i]) );
  }

  /* Start BiCG iterations, working on even sites only */
  
  N_iter = 0;
  while(1) {
    if( N_iter % restart == 0 || 
	( ( RsdCG  <= 0 || RsdCG  > qic->size_r   ) &&
	  ( RRsdCG <= 0 || RRsdCG > qic->size_relr ) ) ) {
      
      /* Provision for starting with dest = 0 the first time */
      if(flag == 0) {
#ifdef CG_DEBUG
	node0_printf("dest_0=0\n");fflush(stdout);
#endif
	FOREVENSITESDOMAIN(i,s) {
	  clear_wvec( &(dest[i]) );
	}
	flag = 1;

      } else {
	
	/* Test true residual for convergence */
	/* r=src[1]-[L^(-1)*M*U^(-1)]*dest (even sites) */

#if 0
	mult_this_ldu_field(gen_clov, dest, tmp, EVEN);
	dslash_w_field_special(dest, mp, PLUS, ODD, tago, is_startedo);
	is_startedo = 1;
	mult_this_ldu_field(gen_clov, mp, tmp, ODD);
	dslash_w_field_special(tmp, mp, PLUS, EVEN, tage, is_startede);
	is_startede = 1;
	FOREVENSITESDOMAIN(i,s) {
	  scalar_mult_add_wvec( &(tmp[i]), &(mp[i]), MKsq, &(mp[i]) );
	  scalar_mult_add_wvec( &(src[i]),
				&(mp[i]), -1.0, &(r[i])     );
	}
#endif
	ilu_DRD(dest, mp, tmp, mp, PLUS, tago, &is_startedo,
		  tage, &is_startede);
	FOREVENSITESDOMAIN(i,s) {
	  scalar_mult_add_wvec( &(tmp[i]), &(mp[i]), MKsq, &(mp[i]) );
	  scalar_mult_add_wvec( &(src[i]),
				&(mp[i]), -1.0, &(r[i])     );
	}
      }

      rsq = 0.0;
      FOREVENSITESDOMAIN(i,s) {
	rsq += (double)magsq_wvec( &(r[i]) );
	copy_wvec( &(r[i]), &(p[i]) );
	copy_wvec( &(r[i]), &(rv[i]) );
      }
      g_doublesum(&rsq);
      rvr = dcmplx(rsq,(double)0.0);

      qic->final_rsq    = rsq/size_src2;
      qic->final_relrsq = relative_residue(r, dest, EVEN);
      
#ifdef CG_DEBUG
      node0_printf("start,   true residue= %e, rel residue= %e\n",
		   qic->final_rsq, qic->final_relrsq);
      fflush(stdout);
#endif
      /* Quit when true residual and true relative residual are within
	 tolerance or when we exhaust iterations or restarts */

      if( N_iter >= MaxCG || 
	  nrestart >= max_restarts ||
	  ( ( RsdCG  <= 0 || RsdCG  > qic->final_rsq ) &&
	    ( RRsdCG <= 0 || RRsdCG > qic->final_relrsq ) ) ) break;

      nrestart++;

    }

    /*   mp = M(u)*p */
#if 0
    mult_this_ldu_field(gen_clov, p, tmp, EVEN);
    dslash_w_field_special(p, mp, PLUS, ODD, tago, is_startedo);
    is_startedo = 1;
    mult_this_ldu_field(gen_clov, mp, tmp, ODD);
    dslash_w_field_special(tmp, mp, PLUS, EVEN, tage, is_startede);
    is_startede = 1;
#endif
    ilu_DRD(p, mp, tmp, mp, PLUS, tago, &is_startedo, tage, &is_startede);
    
    /* rvv = <rv|mp> */
    rvv = dcmplx((double)0.0,(double)0.0);
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( &(tmp[i]), &(mp[i]), MKsq, &(mp[i]) );
      ctmp = wvec_dot( &(rv[i]), &(mp[i]) );
      CSUM(rvv,ctmp);
    }
    g_dcomplexsum(&rvv);
    CDIV(rvr,rvv,a);
    
    /* sss = r - a*mp  */
    CMULREAL(a,-1.0,ctmp);
    FOREVENSITESDOMAIN(i,s) {
      c_scalar_mult_add_wvec( &(r[i]), &(mp[i]), &ctmp, &(sss[i]) );
    }
    
    /* ttt = M(u)*sss */
#if 0
    mult_this_ldu_field(gen_clov, sss, tmp, EVEN);
    dslash_w_field_special(sss, sss, PLUS, ODD, tago, is_startedo);
    is_startedo = 1;
    mult_this_ldu_field(gen_clov, sss, tmp, ODD);
    dslash_w_field_special(tmp, ttt, PLUS, EVEN, tage, is_startede);
    is_startede = 1;
#endif
    ilu_DRD(sss, ttt, tmp, sss, PLUS, tago, &is_startedo, tage, &is_startede);

    /* tdots = <ttt|sss>; tsq=|ttt|^2 */
    tdots = dcmplx((double)0.0,(double)0.0);
    tsq = 0.0;
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( &(tmp[i]), &(ttt[i]), MKsq, &(ttt[i]) );
      ctmp = wvec_dot( &(ttt[i]), &(sss[i]) );
      CSUM(tdots, ctmp);
      tsq += (double)magsq_wvec( &(ttt[i]) );
    }
    g_dcomplexsum(&tdots);
    g_doublesum(&tsq);
    
    omega.real = tdots.real/tsq;
    omega.imag = tdots.imag/tsq;
    
    /* dest = dest + omega*sss + a*p */
    /* r = sss - omega*ttt */
    omegam.real = -omega.real;
    omegam.imag = -omega.imag;
    rsq = 0.0;
    rvro = rvr;
    rvr = dcmplx((double)0.0,(double)0.0);
    FOREVENSITESDOMAIN(i,s) {
      c_scalar_mult_add_wvec(&(dest[i]),&(sss[i]), &omega,  &(dest[i]) );
      c_scalar_mult_add_wvec(&(dest[i]),&(p[i]),   &a,      &(dest[i]) );
      c_scalar_mult_add_wvec(&(sss[i]),   &(ttt[i]), &omegam, &(r[i])      );
      ctmp = wvec_dot( &(rv[i]), &(r[i]) );
      CSUM(rvr, ctmp);
      rsq += (double)magsq_wvec( &(r[i]) );
    }
    g_dcomplexsum(&rvr);
    g_doublesum(&rsq);
    
    qic->size_relr = relative_residue(r, dest, EVEN);
    
    CDIV(rvr,rvro,b);
    CMUL(a,b,ctmp);
    CDIV(ctmp,omega,b);
    
    /*   p = r + b*(p - omega*mp)  */
    FOREVENSITESDOMAIN(i,s) {
      c_scalar_mult_add_wvec( &(p[i]),  &(mp[i]), &omegam, &(p[i]) );
      c_scalar_mult_add_wvec( &(r[i]),  &(p[i]),     &b,      &(p[i]) );
    }
    
    qic->size_r = rsq/size_src2;
    
    N_iter++;
#ifdef CG_DEBUG
    node0_printf("iteration= %d, residue= %e, rel residue= %e\n",N_iter,
		 qic->size_r, qic->size_relr);
    fflush(stdout);
#endif
  }

  qic->final_iters = N_iter;
  qic->final_restart = nrestart;

#ifdef CGTIME
  dtime += dclock();
#endif
  if(this_node==0){
    if(N_iter==0)
      printf("BICGILU: NO iterations taken size_r= %.2e rel %.2e\n",
	     qic->final_rsq, qic->final_relrsq);
#ifdef CGTIME
    else
      printf("CGTIME: time = %.2e (bicgilu) size_r= %.2e relr= %.2e iters= %d MF = %.1f\n",
	     dtime,qic->final_rsq,qic->final_relrsq,N_iter,
	     (double)8742*N_iter*even_sites_on_node/(dtime*1e6));
#endif
    fflush(stdout);
  }
  
#ifdef CG_DEBUG
  if( (RsdCG > 0 && qic->final_rsq > RsdCG) &&
      (RRsdCG > 0 && qic->final_relrsq > RRsdCG) )
    {
      node0_printf(" BiCG_ILU: Not Converged: size_r=%.2e rel %.2e \n",
		   qic->final_rsq, qic->final_relrsq);fflush(stdout);
    }
#endif

#if 0
  /* Reconstruct solution on odd sites */

  /* --------- dest_o = U^(-1)_oo/R_o dest_o + U^(-1)_oe dest_e --------- */
  /* --------- dest_o =  R^(-1)*dest_o = dest_o + K/R_o D_oe dest_e ------- */
  /* mp_o = D_oe * dest_e */
  dslash_w_field_special(dest, mp, PLUS, ODD, tago, is_startedo);
  is_startedo = 1;

  /* mp_o = dest_o + K D_oe * dest_e (remember dest_o = original src_o still)*/
  FORODDSITESDOMAIN(i,s) {
    scalar_mult_add_wvec( &(dest[i]), &(mp[i]), Kappa, &(mp[i]) );
  }
  /* dest_o = 1/R_o dest_o + K/R_o D_oe * dest_e */
  mult_this_ldu_field(gen_clov, mp, dest, ODD);
#endif

  ilu_xfm_dest(dest, mp, Kappa, &is_startedo, tago);

  free(tmp); free(mp); free(rv); free(sss); free(r); free(p); free(ttt);

  for( i=XUP; i <= TUP; i++) {
    if(is_startedo)cleanup_gather(tago[i]);
    if(is_startedo)cleanup_gather(tago[OPP_DIR(i)]);
    if(is_startede)cleanup_gather(tage[i]);
    if(is_startede)cleanup_gather(tage[OPP_DIR(i)]);
  }
  is_startede = is_startedo = 0;

  cleanup_tmp_links();
  cleanup_dslash_wtemps();
  return(N_iter);

} /* bicgilu_cl_field */

int bicgilu_cl_site(    /* Return value is number of iterations taken */
    field_offset src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    )
{
  int i;
  site *s;
  wilson_vector *t_src, *t_dest;
  int iters;

  t_src  = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  t_dest = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));

  if(t_src == NULL || t_dest == NULL){
    printf("bicgilu_cl_site(%d): Can't allocate src and dest\n",this_node);
    terminate(1);
  }

  /* copy src and dest to temporary */
  FORALLSITES(i,s) {
    t_src[i] = *(wilson_vector *)F_PT(s,src);
    t_dest[i] = *(wilson_vector *)F_PT(s,dest);
  }

  iters = bicgilu_cl_field(t_src, t_dest, qic, dmp);

  /* copy dest back */
  FORALLSITES(i,s) {
    *(wilson_vector *)F_PT(s,dest) = t_dest[i];
  }

  free(t_dest);
  free(t_src);

  return iters;
} /* bicgilu_cl_site */


/* Plain Wilson variant */

int bicgilu_w_field(     /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the WILSON Dirac matrix */
    )
{
  /* Input parameters are for Dirac Wilson! Need remapping. */
  dirac_clover_param dcp;

  map_dwp_to_dcp(&dcp, (dirac_wilson_param *)dmp);
  return bicgilu_cl_field(src, dest, qic, (void *)&dcp);
}


/* Site-based - kept for backward compatibility */
int bicgilu_w_site(      /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    )
{
  /* Input parameters are for Dirac Wilson! Need remapping. */
  dirac_clover_param dcp;

  map_dwp_to_dcp(&dcp, (dirac_wilson_param *)dmp);
  return bicgilu_cl_site(src, dest, qic, (void *)&dcp);
}
