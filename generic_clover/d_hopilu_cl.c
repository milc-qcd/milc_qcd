/******* d_hopilu_cl.c - hopping expansion for heavy clover fermions ****/
/* MIMD version 7 */
/* Memory stingy version
   3/29/00 EVENFIRST is the rule now. CD.
   */

/* Modifications:

   4/27/07 Made consistent with other inverters and added Wilson option - CD
   8/02/01 Uses dslash_w_site_special - CD
   4/25/98 Initial version based on UMH's d_cgilu_cl_lean.c CD
   */

#include "generic_clover_includes.h"

/*#define CG_DEBUG*/

int hopilu_cl_field(     /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    )
{
  /* Unpack required members of the structures */
  int MaxHOP = qic->max;      /* maximum number of iterations */
  Real RsdHOP = qic->resid * qic->resid;      /* desired residual - 
				normalized as (r*r)/(src_e*src_e) */
  Real RRsdHOP = qic->relresid * qic->relresid;  /* desired relative residual - */

  dirac_clover_param *dcp 
    = (dirac_clover_param *)dmp; /* Cast pass-through pointer */
  
  Real Kappa = dcp->Kappa;     /* hopping */
  Real Clov_c = dcp->Clov_c;   /* Perturbative clover coeff */
  Real U0 = dcp->U0;           /* Tadpole correction to Clov_c */
  
  /* End of unpacking required members of structures */
  
  wilson_vector *tmp=NULL,*mp=NULL,*r=NULL;

  int N_iter;
  register int i;
  register site *s;
  Real size_src2;
  double sr;
  register Real Ksq = Kappa*Kappa;
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
  
  /* Invert BOTH R_o and R_e in place in "clov" and "clov_diag" */
  compute_clovinv(gen_clov, EVENANDODD);
  
  /* now we can allocate temporary variables and copy them */
  /* PAD may be used to avoid cache trashing */
#define PAD 0
    
  tmp    = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  mp     = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  r      = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));

  if(tmp == NULL || mp == NULL || r == NULL){
    printf("hopilu_cl_field(%d): No room for temporaries\n", this_node);
    terminate(1);
  }

  /* copy src to temporary */
  FORALLSITES(i,s) {
    r[i] = src[i];
  }

  /* HOP_ILU: */
  
#ifdef CGTIME
  dtime = -dclock();
#endif

  /* r <- L^(-1) r */
  ilu_xfm_source(dest, r, mp, Kappa, &is_startede, tage);

  /* --------- dest_o = src_o ---------- */

  /* set dest = src on odd sites, whatever else you do
     (even if we restart with a nonzero solution vector, the end of the
     subroutine rebuilds the odd component from the even one. The (trivial)
     solution of the odd component of the equation is dest = src, before
     we rotate back to the basis in which  M is not checkerboard-diagonal) */
  FORODDSITES(i,s) {
    copy_wvec( &(r[i]), &(dest[i]) ); 
  } 

  /* --------- dest_e <- R_e^(-1) srce_e ---------- */
  mult_this_ldu_field( gen_clov, r, dest, EVEN);
  
  /* ------------ r_e <- dest_e -------------- */
  FOREVENSITESDOMAIN(i,s) {
    copy_wvec( &(dest[i]), &(r[i]) );
  }

  sr=0.0;
  FOREVENSITESDOMAIN(i,s) {
    sr += (double)magsq_wvec( &(r[i]) );
  }

  /* --------- size_src = sqrt(r_e*r_e) ---------- */
  g_doublesum(&sr);
  size_src2 = sr;
  
#ifdef CG_DEBUG
  node0_printf("beginning inversion--size_src=%e\n",
	       (double)size_src2); fflush(stdout);
#endif

  qic->size_r = 0;
  qic->size_relr = 0;
  qic->final_rsq = 0;
  qic->final_relrsq = 0;

  /* --------- Beginning of hop iterations --------- */
  
  for( N_iter = 1; 
       N_iter == 1 ||       /* Force one iteration */
       ( N_iter <= MaxHOP &&
	 ( ( (RsdHOP  > 0) && (RsdHOP  < qic->size_r   ) ) ||
	   ( (RRsdHOP > 0) && (RRsdHOP < qic->size_relr) ) ) );
       N_iter++ )
      {
      
      /* ---- r_e = H r_e =  K^2/R_e D_eo/R_o D_oe r_e ----- */
      /* mp_o = D_oe r_e */
      dslash_w_field_special(r, mp, PLUS, ODD, tago, is_startedo);
      is_startedo = 1;
      /* tmp_o = 1/R_o D_oe r_e */
      mult_this_ldu_field(gen_clov, mp, tmp, ODD);
      /* mp_e = D_eo/R_o D_oe r_e */
      dslash_w_field_special(tmp, mp, PLUS, EVEN, tage, is_startede);
      is_startede = 1;
      /* r_e = 1/R_e D_eo/R_o D_oe r_e */
      mult_this_ldu_field(gen_clov, mp, r, EVEN);
      /* r_e = K^2/R_e D_eo/R_o D_oe r_e */
      sr=0.0;
      FOREVENSITESDOMAIN(i,s) {
	scalar_mult_wvec( &(r[i]), Ksq, &(r[i]) );
	sr += (double)magsq_wvec( &(r[i]) );
      }
      g_doublesum(&sr);
      qic->size_r = sr/size_src2;
      qic->size_relr = relative_residue(r, dest, EVEN);
      qic->final_rsq = qic->size_r;
      qic->final_relrsq = qic->size_relr;

#ifdef CG_DEBUG
      node0_printf("iteration= %d, residue= %e, rel residue= %e\n",N_iter,
		   (double)(qic->size_r), (double)(qic->size_relr));
      fflush(stdout);
#endif
      
      /* --------  dest_e <- dest_e + r_e -------*/
      FOREVENSITESDOMAIN(i,s){
	add_wilson_vector( &(dest[i]), &(r[i]), &(dest[i]));
      }
      
    }

  /* ------------------------------------ */
  /* --------- End of iterations --------- */

  qic->final_iters = N_iter;
  qic->final_restart = 1;

#ifdef CGTIME
  dtime += dclock();
#endif
  if(this_node==0){
    if(N_iter==0)
      printf("HOPILU: NO iterations taken size_r= %.2e rel %.2e\n",
	     qic->final_rsq, qic->final_relrsq);
#ifdef CGTIME
    else
      printf("CGTIME: time = %.2e (hopilu) size_r= %.2e relr= %.2e iters= %d MF = %.1f\n",
	     dtime,qic->final_rsq,qic->final_relrsq,N_iter,
	     (double)3744*N_iter*even_sites_on_node/(dtime*1e6));
#endif
    fflush(stdout);
  }

#ifdef CG_DEBUG
  if( (RsdHOP > 0 && qic->final_rsq > RsdHOP) || 
      (RRsdHOP > 0 && qic->final_relrsq > RRsdHOP) )
    {
      node0_printf(" HOPILU: Not Converged: size_r=%.2e rel %.2e \n",
		   qic->final_rsq, qic->final_relrsq);fflush(stdout);
    }
#endif
  
  /* Reconstruct solution on odd sites */

  ilu_xfm_dest(dest, mp, Kappa, &is_startedo, tago);

  free(tmp); free(mp); free(r);

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

} /* hopilu_cl_field */

int hopilu_cl_site(      /* Return value is number of iterations taken */
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
  }

  iters = hopilu_cl_field(t_src, t_dest, qic, dmp);

  /* copy dest back */
  FORALLSITES(i,s) {
    *(wilson_vector *)F_PT(s,dest) = t_dest[i];
  }

  free(t_dest);
  free(t_src);

  return iters;
} /* hopilu_cl_site */

/* Plain Wilson variant */

int hopilu_w_field(     /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the WILSON Dirac matrix */
    )
{
  /* Input parameters are for Dirac Wilson! Need remapping. */
  dirac_clover_param dcp;

  map_dwp_to_dcp(&dcp, (dirac_wilson_param *)dmp);
  return hopilu_cl_field(src, dest, qic, (void *)&dcp);
}


/* Site-based - kept for backward compatibility */
int hopilu_w_site(      /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    )
{
  /* Input parameters are for Dirac Wilson! Need remapping. */
  dirac_clover_param dcp;

  map_dwp_to_dcp(&dcp, (dirac_wilson_param *)dmp);
  return hopilu_cl_site(src, dest, qic, (void *)&dcp);
}
