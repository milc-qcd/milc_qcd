/******* d_cgilu_w_lean.c - CG-ILU for  Wilson fermions ****/
/* MIMD version 6 */

/* Modifications:
   7/18/01 Using dslash_w_special CD
   4/26/98 Moved parameters to structures CD
   8/29/97 ANSI prototyping C. D.
   Code created by U.M.H.
   */

/* double precision accumulations */
/* Memory stingy version
   "r" overwrites src on even sites
   "p" overwrites src on odd sites
   3/29/00 EVENFIRST is the rule now. CD.
   */

/* Requires qic wilson vector temporary wv2 */

/* The source vector is in "src", and the initial guess and answer
   in "dest".  "r" is the residual vector, which is a pointer to src since
   the source  is overwritten to save space 
   and "p" and "mmp" are
   working vectors for the conjugate gradient. 
   MinCG = minimum number of iterations.
   MaxCG = maximum number of iterations.
   size_r = desired residual, quit when we reach it.
   (Square root def for residue size_r = sqrt(r*r))
   ILU resides on parity=EVEN so do only even sites
   */


#include "generic_wilson_includes.h"

/*#define CGTIME*/          /* Uncomment if you want timing info */


int cgilu_w(             /* Return value is number of iterations taken */
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
  field_offset mmp = qic->wv2;    /* size of wilson_vector */

  dirac_wilson_param *dwp 
    = (dirac_wilson_param *)dmp; /* Cast pass-through pointer */
  Real Kappa = dwp->Kappa;     /* hopping */
  /* End of unpacking required members of structures */

  int N_iter;
  register int i;
  register site *s;
  Real size_src;
  double cp,d,dsize_r,dsize_src;
  register Real a, b, c;
  register Real MKsq = -Kappa*Kappa;
  register field_offset r,p;
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
  r = src;
  p = src + even_sites_on_node*sizeof(site);
  /* This disgusting trick makes p for each even site actually be
     src on some corresponding odd site */
  /**if(this_node==0)printf("CGILU: p=%d\n",p);**/
  
  /* CG_ILU: */
  
  /* Start Inversion */
  
  /* src = L^(-1)*src */
#ifdef CGTIME
  dtime = -dclock();
#endif
  
  dslash_w_special(src,mmp,PLUS,EVEN,tage,is_startede);
  is_startede = 1;
  
  /* Normalisation  */
  dsize_src=0.0;
  FOREVENSITES(i,s) {
    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,src), 
			  (wilson_vector *)F_PT(s,mmp),
			 Kappa, (wilson_vector *)F_PT(s,src) );
    dsize_src += (double)magsq_wvec( (wilson_vector *)F_PT(s,src) );
  }
  g_doublesum(&dsize_src);
  size_src = (Real)sqrt(dsize_src);
  
  /**if(this_node==0)printf("beginning inversion--size_src=%e\n",
    (double)size_src);**/
  
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
    copy_wvec( (wilson_vector *)F_PT(s,src),
	      (wilson_vector *)F_PT(s,dest) );
  }
  
  
  /* code if you want to start with dest=0... */
  if(flag == 0) {
    /**if(this_node==0)printf("dest_0=0\n");**/
    FOREVENSITES(i,s) {
      clear_wvec( (wilson_vector *)F_PT(s,dest) );
    }
    dsize_r = 1.00;
    qic->size_r = dsize_r;
    /**if(this_node==0)printf("size_r=%e\n",(double)qic->size_r));**/
    
  }
  /* code if you want to start dest with some particular starting value... */
  /* r=src[1]-[L^(-1)*M*U^(-1)]*dest */
  if(flag != 0) {
    /**if(this_node==0)printf("dest_0  !=0\n");**/
    /* we use mmp temporarily to construct r */
    dslash_w_special(dest,mmp,PLUS,ODD,tago,is_startedo);
    is_startedo = 1;
    dslash_w_special(mmp,mmp,PLUS,EVEN,tage,is_startede);
    is_startede = 1;
    FOREVENSITES(i,s) {
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest),
			   (wilson_vector *)F_PT(s,mmp),MKsq, 
			    (wilson_vector *)F_PT(s,mmp) );
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,r),
			   (wilson_vector *)F_PT(s,mmp),-1.0, 
			    (wilson_vector *)F_PT(s,r) );
    }
    
    dsize_r=0.0;
    FOREVENSITES(i,s) {
      dsize_r += magsq_wvec( (wilson_vector *)F_PT(s,r) );
    }
    g_doublesum(&dsize_r);
    qic->size_r = (Real)sqrt(dsize_r)/size_src;
    /**if(this_node==0)    printf("beginning inversion--size_r=%e\n",
      (double)(qic->size_r));**/
    
  }
  
  /*  p = [L^(-1)*M*U^(-1)]_dag*r  */
  dslash_w_special(r,mmp,MINUS,ODD,tago,is_startedo);
  is_startedo = 1;
  dslash_w_special(mmp,p,MINUS,EVEN,tage,is_startede);
  is_startede = 1;
  
  /* cp = |p|^2 */
  cp=0.0;
  FOREVENSITES(i,s) {
    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,r),
			 (wilson_vector *)F_PT(s,p), MKsq, 
			 (wilson_vector *)F_PT(s,p) );
    cp += (double)magsq_wvec( (wilson_vector *)F_PT(s,p) );
  }
  g_doublesum(&cp);
  
  for( N_iter = 0; N_iter < MinCG || (N_iter < MaxCG && RsdCG < qic->size_r); 
      N_iter = N_iter + 1) {
    
    c=cp;
    /*   mmp = M(u)*p */
    dslash_w_special(p,mmp,PLUS,ODD,tago,is_startedo);
    is_startedo = 1;
    dslash_w_special(mmp,mmp,PLUS,EVEN,tage,is_startede);
    is_startede = 1;
    
    /* d = |mmp|^2  */
    d=0.0;
    FOREVENSITES(i,s) {
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,p),
			   (wilson_vector *)F_PT(s,mmp),MKsq, 
			    (wilson_vector *)F_PT(s,mmp) );
      d += (double)magsq_wvec( (wilson_vector *)F_PT(s,mmp) );
    }
    g_doublesum(&d);
    a = (Real)(c/d);
    
    /* dest = dest + a*p  */
    /* r = r - a*mmp */
    FOREVENSITES(i,s) {
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest),
			   (wilson_vector *)F_PT(s,p), a,
			   (wilson_vector *)F_PT(s,dest) );
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,r),
			   (wilson_vector *)F_PT(s,mmp),-a, 
			    (wilson_vector *)F_PT(s,r) );
    }
    
    /*   mmp = M(u)dag*r  */
    
    dslash_w_special(r,mmp,MINUS,ODD,tago,is_startedo);
    is_startedo = 1;
    dslash_w_special(mmp,mmp,MINUS,EVEN,tage,is_startede);
    is_startede = 1;
    /*   cp = |mmp|^2  */
    cp=0.0;
    FOREVENSITES(i,s) {
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,r),
			   (wilson_vector *)F_PT(s,mmp), MKsq, 
			    (wilson_vector *)F_PT(s,mmp) );
      cp += (double)magsq_wvec( (wilson_vector *)F_PT(s,mmp) );
    }
    g_doublesum(&cp);
    b = (Real)(cp/c);
    /*   p = mmp + b*p  */
    dsize_r=0.0;
    FOREVENSITES(i,s) {
      scalar_mult_add_wvec( (wilson_vector *)F_PT(s,mmp),
			   (wilson_vector *)F_PT(s,p), b, 
			    (wilson_vector *)F_PT(s,p) );
      dsize_r += magsq_wvec( (wilson_vector *)F_PT(s,r) );
    }
    g_doublesum(&dsize_r);
    
    qic->size_r = (Real)sqrt(dsize_r)/size_src;
    /**if(this_node==0)printf("iteration= %d, residue= %e\n",N_iter,
      (double)(qic->size_r));**/
  }
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
	     (double)5616*N_iter*even_sites_on_node/(dtime*1e6));
#endif
    fflush(stdout);
  }
  /**  if( (qic->size_r) > RsdCG ) {
    if(this_node==0)printf(" CG_ILU_ Not Converged\n");
    }**/
  
  /* dest = R^(-1)*dest  */
  dslash_w_special(dest,mmp,PLUS,ODD,tago,is_startedo);
  is_startedo = 1;
  FORODDSITES(i,s) {
    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest), 
			  (wilson_vector *)F_PT(s,mmp),
			 Kappa, (wilson_vector *)F_PT(s,dest) );
  }
  
  for( i=XUP; i <= TUP; i++) {
    if(is_startedo)cleanup_gather(tago[i]);
    if(is_startedo)cleanup_gather(tago[OPP_DIR(i)]);
    if(is_startede)cleanup_gather(tage[i]);
    if(is_startede)cleanup_gather(tage[OPP_DIR(i)]);
  }
  is_startede = is_startedo = 0;
    
  return(N_iter);
}
