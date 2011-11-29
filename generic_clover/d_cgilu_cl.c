/******* d_cgilu_cl.c - CG-ILU for  clover fermions ****/
/* MIMD version 7 */

/* Modifications:
   7/18/01 Uses dslash_w_site_special - CD
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

#include "generic_clover_includes.h"

/*#define CG_DEBUG*/

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

int cgilu_cl_field(       /* Return value is number of iterations taken */
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
  int MaxCG = restart*qic->max;      /* maximum number of iterations */
  Real RsdCG = qic->resid * qic->resid;      /* desired residual - 
				normalized as (r*r)/(src_e*src_e) */
  Real RRsdCG = qic->relresid * qic->relresid;  /* desired relative residual - */
  int flag = qic->start_flag;   /* 0: use a zero initial guess; 1: use dest */

  dirac_clover_param *dcp 
    = (dirac_clover_param *)dmp; /* Cast pass-through pointer */

  Real Kappa = dcp->Kappa;     /* hopping */
  Real Clov_c = dcp->Clov_c;   /* Perturbative clover coeff */
  Real U0 = dcp->U0;           /* Tadpole correction to Clov_c */
  /* End of unpacking required members of structures */

  wilson_vector *tmp=NULL,*mp=NULL,*r=NULL,*p=NULL;
  int N_iter;
  register int i;
  register site *s;
  Real size_src, size_src2;
  double sr, cp, c, d;
  register Real a, b;
  register Real MKsq = -Kappa*Kappa;
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
  
  /* Compute R_e and R_o and put in "clov" and "clov_diag" */
  compute_clov(gen_clov, CKU0);

  /* Invert R_o only, leaving R_e on even sites and 1/R_o on odd sites 
     in "clov" and "clov_diag" */
  compute_clovinv(gen_clov, ODD);
  
  /* now we can allocate temporary variables and copy then */
  /* PAD may be used to avoid cache trashing */
#define PAD 0
    
  tmp    = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  mp  = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  r      = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  p      = r  + even_sites_on_node;

  if(tmp == NULL || mp == NULL || r == NULL){
    printf("cgilu_cl_field(%d): No room for temporaries\n", this_node);
    terminate(1);
  }

  /* CG_ILU: */
  
#ifdef CGTIME
  dtime = -dclock();
#endif

  /* now we copy dest and src to temporaries */
  FORALLSITES(i,s) {
    r[i] = src[i];
  }

  size_src = ilu_xfm_source(dest, r, mp, Kappa, &is_startede, tage);
  size_src2 = size_src * size_src;

#if 0
  /* Transform source - result is in r (even) - and dest (odd) */
  
  /* ---------  src = L^(-1)*src  ------------- */

  /* ( src_o ) = ( 1           0 ) ( src_o )
     ( sec_e )   (-K D_eo/R_o  1 ) ( src_e )

     */

  /* mp_o = 1/R_o srce_e */
  mult_this_ldu_field(gen_clov, r, mp, ODD);
  /* mp_e = D_eo/R_o srce_e */
  dslash_w_field_special(mp, mp, PLUS, EVEN, tage, is_startede);
  is_startede = 1;
  
  /* src_e = srce_e + K D_eo/R_o srce_e */
  /* (leaving src_o = src_o)   */
  sr=0.0;
  FOREVENSITESDOMAIN(i,s) {
    scalar_mult_add_wvec( &(r[i]), &(mp[i]), Kappa, &(r[i]) );
    sr += (double)magsq_wvec( &(r[i]) );
    /* Save transformed source: Overwrite src on even sites */
    copy_wvec( &(r[i]), &(src[i]) );
  }

  /* --------- size_src = sqrt(src_e*src_e) ---------- */
  g_doublesum(&sr);
  size_src2 = sr;
#endif

  /* Save transformed source: Overwrite src on even sites */
  FOREVENSITESDOMAIN(i,s) {
    copy_wvec( &(r[i]), &(src[i]) );
  }

#ifdef CG_DEBUG
  node0_printf("beginning inversion--size_src=%e\n",
	       (double)size_src2); fflush(stdout);
#endif
  
  /* --------- dest_o = src_o ---------- */
  /* Initial guess */
  /* set dest = src on odd sites, whatever else you do
     (even if we restart with a nonzero solution vector, the end of the
     subroutine rebuilds the odd component from the even one. The (trivial)
     solution of the odd component of the equation is dest = src, before
     we rotate back to the basis in which  M is not checkerboard-diagonal) */
  FORODDSITESDOMAIN(i,s) {
    copy_wvec( &(r[i]), &(dest[i]) );
  }
  
  /* ---------------------------------------------------- */
  /* ----- Beginning of cg iterations on even sites ----- */
  
  N_iter = 0;
  qic->size_r = 0;
  qic->size_relr = 0;
  qic->final_rsq = 0;
  qic->final_relrsq = 0;

  while(1) {
    if( N_iter % restart == 0 || 
	( ( RsdCG  <= 0 || RsdCG  > qic->size_r   ) &&
	  ( RRsdCG <= 0 || RRsdCG > qic->size_relr) ) ) {
      
      /* --------- if flag == 0 set dest_e = 0 ---------- */
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
	/* --------- and compute r_e = src_e[1]-M_e*dest_e --------- */

#if 0
	/* tmp_e = R_e dest_e */
	mult_this_ldu_field(gen_clov, dest, tmp, EVEN);
	/* mp_o = D_oe dest_e */
	dslash_w_field_special(dest, mp, PLUS, ODD, tago, is_startedo);
	is_startedo = 1;
	/* tmp_o = 1/R_o D_oe dest_e */
	mult_this_ldu_field(gen_clov, mp, tmp, ODD);
	/* mp_e = D_eo/R_o D_oe dest_e */
	dslash_w_field_special(tmp, mp, PLUS, EVEN, tage, is_startede);
	is_startede = 1;
#endif
	ilu_DRD(dest, mp, tmp, mp, PLUS, tago, &is_startedo,
		  tage, &is_startede);

	/* mp_e = R_e dest_e - K^2 D_eo/R_o D_oe dest_e = M_e dest_e */
	/* r_e = src_e - M_e dest_e */
	FOREVENSITESDOMAIN(i,s) {
	  scalar_mult_add_wvec( &(tmp[i]), &(mp[i]), MKsq, &(mp[i]) );
	  scalar_mult_add_wvec( &(src[i]), 
				&(mp[i]), -1.0, &(r[i])     );
	}
      }

      /* --------- size_r = sqrt(r * r)/size_src --------- */
      sr=0.0;
      FOREVENSITESDOMAIN(i,s) {
	sr += (double)magsq_wvec( &(r[i]) );
      }
      g_doublesum(&sr);
      qic->final_rsq = sr/size_src2;
      qic->final_relrsq = relative_residue(r, dest, EVEN);
  
      
#ifdef CG_DEBUG
      node0_printf("start,   true residue= %e, rel residue= %e\n",
		   (double)(qic->final_rsq), (double)(qic->final_relrsq));
      fflush(stdout);
#endif

      /* Quit when true residual and true relative residual are within
	 tolerance or when we exhaust iterations or restarts */

      if( N_iter >= MaxCG || 
	  nrestart >= max_restarts ||
	  ( ( RsdCG  <= 0 || RsdCG  > qic->final_rsq   ) &&
	    ( RRsdCG <= 0 || RRsdCG > qic->final_relrsq) ) ) break;

      nrestart++;

      /* --------- p_e = M_e_dag*r_e --------- */
#if 0
      mult_this_ldu_field(gen_clov, r, tmp, EVEN);
      dslash_w_field_special(r, mp, MINUS, ODD, tago, is_startedo);
      is_startedo = 1;
      mult_this_ldu_field(gen_clov, mp, tmp, ODD);
      dslash_w_field_special(tmp, p, MINUS, EVEN, tage, is_startede);
      is_startede = 1;
#endif
      ilu_DRD(r, p, tmp, mp, MINUS, tago, &is_startedo,
	      tage, &is_startede);

      
      /* --------- cp = |p|^2 --------- */
      cp=0.0;
      FOREVENSITESDOMAIN(i,s) {
	scalar_mult_add_wvec( &(tmp[i]), &(p[i]), MKsq, &(p[i]) );
	cp += (double)magsq_wvec( &(p[i]) );
      }
      g_doublesum(&cp);
    }

    /* --------- c = cp -------- */
    c=cp;
    
    /* ---------  mp_e = M_e*p_e --------- */
#if 0
    /* tmp_e = R_e p_e */
    mult_this_ldu_field(gen_clov, p, tmp, EVEN);
    /* mp_o = D_oe p_e */
    dslash_w_field_special(p, mp, PLUS, ODD, tago, is_startedo);
    is_startedo = 1;
    /* tmp_o = 1/R_o D_oe p_e */
    mult_this_ldu_field(gen_clov, mp, tmp, ODD);
    /* mp_e = D_eo/R_o D_oe p_e */
    dslash_w_field_special(tmp, mp, PLUS, EVEN, tage, is_startede);
    is_startede = 1;
#endif
    ilu_DRD(p, mp, tmp, mp, PLUS, tago, &is_startedo,
	    tage, &is_startede);
    
    /* mp_e = R_e p_e - K^2 D_eo/R_o D_oe p_e */
    d=0.0;
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( &(tmp[i]), &(mp[i]), MKsq, &(mp[i]) );
      d += (double)magsq_wvec( &(mp[i]) );
    }
    /* --------- a = c/|mp_e|^2 --------- */
    g_doublesum(&d);
    a = (Real)(c/d);
    
    /* --------- dest_e = dest_e + a*p_e --------- */
    /* --------- r_e = r_e - a*mp_e --------- */
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( &(dest[i]), &(p[i]), a, &(dest[i]) );
      scalar_mult_add_wvec( &(r[i]), &(mp[i]), -a, &(r[i]) );
    }
    
    /* --------- mp_e M_e_dag*r_e --------- */
    
#if 0
    mult_this_ldu_field(gen_clov, r, tmp, EVEN);
    dslash_w_field_special(r, mp, MINUS, ODD, tago, is_startedo);
    is_startedo = 1;
    mult_this_ldu_field(gen_clov, mp, tmp, ODD);
    dslash_w_field_special(tmp, mp, MINUS, EVEN, tage, is_startede);
    is_startede = 1;
#endif
    ilu_DRD(r, mp, tmp, mp, MINUS, tago, &is_startedo,
	    tage, &is_startede);

    cp=0.0;
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( &(tmp[i]), &(mp[i]), MKsq, &(mp[i]) );
      cp += (double)magsq_wvec( &(mp[i]) );
    }
    /* ---------  cp = |mp|^2 --------- */
    g_doublesum(&cp);
    /* ---------  b = |mp|^2/c --------- */
    b = (Real)(cp/c);

    /* ---------  p_e = mp_e + b*p_e --------- */
    sr=0.0;
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( &(mp[i]), &(p[i]), b, &(p[i]) );
      sr += (double)magsq_wvec( &(r[i]) );
    }
    g_doublesum(&sr);
    
    /* --------- size_r = sqrt(r_e * r_e)/size_src --------- */
    qic->size_r = sr/size_src2;
    qic->size_relr = relative_residue(r, dest, EVEN);

    N_iter++;
#ifdef CG_DEBUG
    node0_printf("iteration= %d, residue= %e, rel residue= %e\n",N_iter,
		 (double)(qic->size_r), (double)(qic->size_relr));
    fflush(stdout);
#endif
  }

  /* ------------------------------------ */
  /* --------- End of iterations --------- */

  qic->final_iters = N_iter;
  qic->final_restart = nrestart;

#ifdef CGTIME
  dtime += dclock();
#endif
  if(this_node==0){
    if(N_iter==0)
      printf("CGILU: NO iterations taken size_r= %.2e rel %.2e\n",
	     qic->final_rsq, qic->final_relrsq);
#ifdef CGTIME
    else
      printf("CGTIME: time = %.2e (cgilu) size_r= %.2e relr= %.2e iters= %d MF = %.1f\n",
	     dtime,qic->final_rsq,qic->final_relrsq,N_iter,
	     (double)8064*N_iter*even_sites_on_node/(dtime*1e6));
#endif
    fflush(stdout);
  }
  
#ifdef CG_DEBUG
  if( (RsdCG > 0 && qic->final_rsq > RsdCG) || 
      (RRsdCG > 0 && qic->final_relrsq > RRsdCG) )
    {
      node0_printf(" BiCG_ILU: Not Converged: size_r=%.2e rel %.2e \n",
		   qic->final_rsq, qic->final_relrsq);fflush(stdout);
    }
#endif

  ilu_xfm_dest(dest, mp, Kappa, &is_startedo, tago);

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

  is_startede = is_startedo = 0;
#endif

  free(tmp); free(mp); free(r);

  for( i=XUP; i <= TUP; i++) {
    if(is_startedo)cleanup_gather(tago[i]);
    if(is_startedo)cleanup_gather(tago[OPP_DIR(i)]);
    if(is_startede)cleanup_gather(tage[i]);
    if(is_startede)cleanup_gather(tage[OPP_DIR(i)]);
  }
  cleanup_tmp_links();
  cleanup_dslash_wtemps();
  return(N_iter);

} /* cgilu_cl_field */

/* Site-based - kept for backward compatibility */

int cgilu_cl_site(       /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
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
    printf("cgilu_cl_site(%d): Can't allocate src and dest\n",this_node);
    terminate(1);
  }

  /* copy src to temporary */
  FORALLSITES(i,s) {
    t_src[i] = *(wilson_vector *)F_PT(s,src);
    t_dest[i] = *(wilson_vector *)F_PT(s,dest);
  }

  iters = cgilu_cl_field(t_src, t_dest, qic, dmp);

  /* copy dest back */
  FORALLSITES(i,s) {
    *(wilson_vector *)F_PT(s,dest) = t_dest[i];
  }

  free(t_dest);
  free(t_src);

  return iters;
} /* cgilu_cl_site */

/* Plain Wilson variants */

int cgilu_w_field(       /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the WILSON Dirac matrix */
    )
{
  /* Input parameters are for Dirac Wilson! Need remapping. */
  dirac_clover_param dcp;

  map_dwp_to_dcp(&dcp, (dirac_wilson_param *)dmp);
  return cgilu_cl_field(src, dest, qic, (void *)&dcp);
}


/* Site-based - kept for backward compatibility */
int cgilu_w_site(        /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    )
{
  /* Input parameters are for Dirac Wilson! Need remapping. */
  dirac_clover_param dcp;

  map_dwp_to_dcp(&dcp, (dirac_wilson_param *)dmp);
  return cgilu_cl_site(src, dest, qic, (void *)&dcp);
}
