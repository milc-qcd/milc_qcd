/******* d_bicgilu_w_lean.c - BiCGstab-ILU for  Wilson fermions ****/
/* MIMD version 6 */

/* Modifications:
   7/17/01 Use dslash_special CD
   8/29/97 ANSI prototyping and added comments C. D.
   Code created by U.M.H.
   */

/* Memory stingy version
   "r" overwrites src on even sites
   "p" overwrites src on odd sites
   3/29/00 EVENFIRST is the rule now. CD.
   */

/* Requires qic wilson vector temporaries wv2, wv3, wv4 */

/* The source vector is in "src", and the initial guess and answer
   in "dest".  "r" is the residual vector, which is a pointer to src since
   the source  is overwritten to save space 
   and "p" and "mmp" are
   working vectors for the conjugate gradient. 
   MaxCG = maximum number of iterations.
   size_r = desired residual, quit when we reach it.
   (Square root def for residue size_r = sqrt(r*r))
   ILU resides on parity=EVEN so do only even sites
   */

#include "generic_wilson_includes.h"


int bicgilu_w(          /* Return value is number of iterations taken */
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

  wilson_vector *mmp,*rv,*sss,*r,*p,*ttt,*t_dest;
  int static first_congrad = 1;

  dirac_wilson_param *dwp 
    = (dirac_wilson_param *)dmp; /* Cast pass-through pointer */
  Real Kappa = dwp->Kappa;     /* hopping */
  /* End of unpacking required members of structures */

  int N_iter;
  register int i;
  register site *s;
  Real size_src;
  double rsq, tsq;
  complex ctmp, a, b;
  complex omega, omegam;
  double_complex tdots,rvro,rvv,rvr;
  register Real MKsq = -Kappa*Kappa;
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

    /* now we can allocate temporary variables and copy then */
    /* PAD may be used to avoid cache trashing */
#define PAD 0
    
    if(first_congrad) {
      mmp    = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
      rv     = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
      sss    = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
      r      = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
      t_dest = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
      p      = r  + even_sites_on_node;
      ttt    = rv + even_sites_on_node;
      first_congrad = 0;
    }

  
  /* BiCGstab_ILU: */
  
  /* Start Inversion */
  
  /* src = L^(-1)*src */
#ifdef CGTIME
  dtime = -dclock();
#endif
  
  /* now we copy dest and src to temporaries */
  FORALLSITES(i,s) {
    t_dest[i] = *(wilson_vector *)F_PT(s,dest);
    r[i] = *(wilson_vector *)F_PT(s,src);
  }

  dslash_on_temp_special(r,mmp,PLUS,EVEN, tage, is_startede);
  is_startede = 1;
  
  /* Normalization  */
  rsq = 0.0;
  FOREVENSITES(i,s) {
    scalar_mult_add_wvec( &r[i], &mmp[i], Kappa, &r[i] );
    rsq += (double)magsq_wvec( &r[i] );
  }
  fflush(stdout);
  g_doublesum(&rsq);
  size_src = (Real)sqrt(rsq);
  
  /** if(this_node==0)printf("beginning inversion--size_src=%e\n",
      (double)size_src); **/
  
  /* r and p overwrite the source (bad for
     dynamical fermions but good for quenched calcs) */
  /* set r = src --- nothing to do */
  
  /* Initial guess */
  /* set dest = src on odd sites, whatever else you do
     (even if we restart with a nonzero solution vector, the end of the
     subroutine rebuilds the odd component from the even one. The (trivial)
     solution of the odd component of the equation is dest = src, before
     we rotate back to the basis in which  M is not checkerboard-diagonal) */
  FORODDSITES(i,s) {
    copy_wvec( &r[i], &t_dest[i] );
  }
  
  
  /* code if you want to start with dest=0... */
  if(flag == 0) {
    /** if(this_node==0)printf("dest_0=0\n"); **/
    FOREVENSITES(i,s) {
      clear_wvec( &t_dest[i] );
    }
  }
  /* code if you want to start dest with some particular starting value... */
  /* r=src[1]-[L^(-1)*M*U^(-1)]*dest */
  if(flag != 0) {
    /** if(this_node==0)printf("dest_0  !=0\n"); **/
    /* we use mmp temporarily to construct r */
    dslash_on_temp_special(t_dest,mmp,PLUS,ODD,tago,is_startedo);
    is_startedo = 1;
    dslash_on_temp_special(mmp,mmp,PLUS,EVEN,tage,is_startede);
    is_startede = 1;
    FOREVENSITES(i,s) {
      scalar_mult_add_wvec( &t_dest[i], &mmp[i], MKsq, &mmp[i] );
      scalar_mult_add_wvec( &r[i], &mmp[i], -1.0, &r[i] );
    }
  }
  
  rsq = 0.0;
  FOREVENSITES(i,s) {
    rsq += (double)magsq_wvec( &r[i] );
    copy_wvec( &r[i], &p[i] );
    copy_wvec( &r[i], &rv[i] );
  }
  g_doublesum(&rsq);
  qic->size_r = (Real)sqrt(rsq)/size_src;
  rvr = dcmplx(rsq,(double)0.0);
  /** if(this_node==0)    printf("beginning inversion--size_r=%e\n",
      (double)(qic->size_r)); **/
  
  for( N_iter = 0; N_iter < MinCG || (N_iter < MaxCG && RsdCG < qic->size_r); 
      N_iter = N_iter + 1) {
    
    /*   mmp = M(u)*p */
    dslash_on_temp_special(p,mmp,PLUS,ODD,tago,is_startedo);
    is_startedo = 1;
    dslash_on_temp_special(mmp,mmp,PLUS,EVEN,tage,is_startede);
    is_startede = 1;
    
    /* rvv = <rv|mmp> */
    rvv = dcmplx((double)0.0,(double)0.0);
    FOREVENSITES(i,s) {
      scalar_mult_add_wvec( &p[i], &mmp[i], MKsq, &mmp[i] );
      ctmp = wvec_dot( &rv[i], &mmp[i] );
      CSUM(rvv,ctmp);
    }
    g_dcomplexsum(&rvv);
    CDIV(rvr,rvv,a);
    
    /* sss = r - a*mmp  */
    CMULREAL(a,-1.0,ctmp);
    FOREVENSITES(i,s) {
      c_scalar_mult_add_wvec( &r[i], &mmp[i], &ctmp, &sss[i] );
    }
    
    /* ttt = M(u)*sss */
    dslash_on_temp_special(sss,sss,PLUS,ODD,tago,is_startedo);
    is_startedo = 1;
    dslash_on_temp_special(sss,ttt,PLUS,EVEN,tage,is_startede);
    is_startede = 1;
    
    /* tdots = <ttt|sss>; tsq=|ttt|^2 */
    tdots = dcmplx((double)0.0,(double)0.0);
    tsq = 0.0;
    FOREVENSITES(i,s) {
      scalar_mult_add_wvec( &sss[i], &ttt[i], MKsq, &ttt[i] );
      ctmp = wvec_dot( (wilson_vector *)F_PT(s,ttt), 
		       (wilson_vector *)F_PT(s,sss) );
      CSUM(tdots, ctmp);
      tsq += (double)magsq_wvec( &ttt[i] );
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
    FOREVENSITES(i,s) {
      c_scalar_mult_add_wvec( &t_dest[i], &sss[i], &omega,  &t_dest[i] );
      c_scalar_mult_add_wvec( &t_dest[i], &p[i],   &a,      &t_dest[i] );
      c_scalar_mult_add_wvec( &sss[i],    &ttt[i], &omegam, &r[i]      );
      ctmp = wvec_dot( &rv[i], &r[i] );
      CSUM(rvr, ctmp);
      rsq += (double)magsq_wvec( &r[i] );
    }
    g_dcomplexsum(&rvr);
    g_doublesum(&rsq);
    
    CDIV(rvr,rvro,b);
    CMUL(a,b,ctmp);
    CDIV(ctmp,omega,b);
    
    /*   p = r + b*(p - omega*mmp)  */
    FOREVENSITES(i,s) {
      c_scalar_mult_add_wvec( &p[i],  &mmp[i], &omegam, &p[i] );
      c_scalar_mult_add_wvec( &r[i],  &p[i],   &b,      &p[i] );
    }
    
    qic->size_r = (Real)sqrt(rsq)/size_src;
    /** if(this_node==0){printf("iteration= %d, residue= %e\n",N_iter,
      (double)(qic->size_r));fflush(stdout);} **/
  }
#ifdef CGTIME
  dtime += dclock();
#endif
  if(this_node==0){
    if(N_iter==0)
      printf("BiCGILU: NO iterations taken size_r= %.2e\n",qic->size_r);
#ifdef CGTIME
    else
      printf("BiCGILU: time = %.2e size_r= %.2e iters= %d MF = %.1f\n",
	     dtime,qic->size_r,N_iter,
	     (double)6342*N_iter*even_sites_on_node/(dtime*1e6));
#endif
    fflush(stdout);
  }
  /** if( (qic->size_r) > RsdCG ) {
    if(this_node==0){printf(" BiCG_ILU_ Not Converged\n");fflush(stdout);}
    } **/
  
  /* dest = R^(-1)*dest  */
  dslash_on_temp_special(t_dest,mmp,PLUS,ODD,tago,is_startedo);
  is_startedo = 1;
  FORODDSITES(i,s) {
    /* Odd site result is placed in dest */
    scalar_mult_add_wvec( &t_dest[i], &mmp[i], Kappa, 
			  (wilson_vector *)F_PT(s,dest) );
  }
  /* Copy back for even sites and free memory */
  FOREVENSITES(i,s) {
    *(wilson_vector *)F_PT(s,dest) = t_dest[i];
  }
  free(mmp); free(rv); free(sss); free(r); free(t_dest);
  first_congrad = 1;

  for( i=XUP; i <= TUP; i++) {
    if(is_startedo)cleanup_gather(tago[i]);
    if(is_startedo)cleanup_gather(tago[OPP_DIR(i)]);
    if(is_startede)cleanup_gather(tage[i]);
    if(is_startede)cleanup_gather(tage[OPP_DIR(i)]);
  }
  is_startede = is_startedo = 0;
    
  return(N_iter);
}

