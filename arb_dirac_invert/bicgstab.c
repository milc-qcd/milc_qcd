/******* bicgstab.c - biconjugate gradient for SU3/fermions ****/
/* MIMD version 4 */
/* Calls delta0.c for fermion matrix */

/* version of 26 June 96  */

/* This is an implementation of BiCGstab (H. van der Vost, SIAM J. Sci.
   Comp. 13 (1992) 631) used for lattice QCD by A. Frommer et al.,
   WUB 94-10, hep-lat/9404013, Int. J. Mod. Phys., C5 (1994) 1073. */

/* The source vector is in "chi", and the initial guess and answer
   in "psi".  "r" is the residual vector and "rv" the conjugate residual,
   "p", "mp", "tmp", "ttt", "sss" and "v" are working vectors for the
   conjugate gradient.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
*/


#include "arb_dirac_inv_includes.h"
#include <string.h>

int bicgstab( /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (where source is to be created)*/
    field_offset dest,  /* type wilson_vector (answer and initial guess) */
    int MaxCG,          /* maximum number of iterations per restart */
    Real RsdCG,        /* desired residual - 
                           normalized as sqrt(r*r)/sqrt(src_e*src_e */
    Real *size_r,      /* resulting residual */
    int start_flag     /* 0: use a zero initial guess; 1: use dest */
    )
{

  register int i;
  register site *s;
  int iteration; /* counter for iterations */
  double rsq,source_norm,rsqmin,rsqstop;
  msg_tag *tag[8],*tag2[8];

  double dtime,dclock();

  dtime= -dclock();
  iteration = 0;

  /* solve M*psi = chi */
  /* mp <-  M*psi    
   * r,rv,p <- chi - mp
   * rsq = rvr = |r|^2
   * source_norm = |chi|^2
   * rvr = rsq 
   */

  rsq = source_norm = 0.0;
  rsqmin=RsdCG*RsdCG;


    /* code if you want to start with psi=0... 
Kluge to work around Wilson code with point source*/
    if(start_flag == 0) {
        if(this_node==0)printf("psi_0=0\n");
        FORALLSITES(i,s) {
            clear_wvec( ((wilson_vector *)F_PT(s,dest)));
            copy_wvec(((wilson_vector *)F_PT(s,src)), 
			((wilson_vector *)F_PT(s,dest)));
        }
     }   


  delta0( dest, F_OFFSET(mp), PLUS );




  FORALLSITES(i,s){
    sub_wilson_vector( ((wilson_vector *)F_PT(s,src)), &(s->mp), &(s->r) );
    s->rv = s->p = s->r;
    source_norm += (double)magsq_wvec(((wilson_vector *)F_PT(s,src))  );
    rsq += (double)magsq_wvec( &(s->r) );
  }

  g_doublesum( &source_norm );
  g_doublesum( &rsq );

  iteration++ ;	/* iteration counts number of multiplications
		   both by M and M_adjoint */
  total_iters++;

  if(this_node==0){
    printf("bicgstab: source_norm = %e\n",source_norm);
    fflush(stdout);
  } 

  rsqstop = rsqmin * source_norm;
  if( rsq > rsqstop ){
    /* Now have to iterate, not good enough guess! */

    /* solve eq. M psi = chi */

    i = bicgstab_iter( dest ,MaxCG,rsqstop,&rsq,PLUS );
    iteration += i;
    total_iters += i;

    if (this_node==0) {
      printf("bicgstab: iter %d, rsq %g\n",iteration,rsq);
      fflush(stdout);
    } 

  } /* rsq > rsqstop */

  *size_r = rsq;

  return(iteration);

}


/* BICONGRAD */

/* 
 * Subroutine bicgstab_iter does the actual convergence loop
 * plmin = PLUS: solves M * sol = x, MINUS: M_adjoint * sol = x
 */

int bicgstab_iter(
field_offset sol,
int niterr,
double rsqstop,
double *rsq_ptr,
 int plmin)
{
  register int i;
  register site *s;
  int iteration;	/* counter for iterations */
  complex a,b;
  double rsq,ttsq;	/* Sugar's a,b,resid**2,previous resid*2 */
                        /* pkp = cg_p.K.cg_p */
  double dtime,dclock();
  double_complex tdots,rvro,rvv,rvr;
  complex omega,omegam,ctmp;

/*
  complex test;
*/



  /* main loop for inversion - do until convergence or time to restart */
  /* 
   * v <- M * p                  if PLUS, M_adjoint if MINUS
   * rvv <- rv_adjoint.v
   * a <- rvr/rvv
   * sss <- r - a*v
   * ttt <- M * sss              same thing
   * tdots <- ttt.sss
   * ttsq  <- |ttt|^2
   * omega <- tdots/ttsq
   * sol <- sol + omega*sss
   * sol <- sol + a*p
   * r <- sss - omega*ttt
   * rsq <- |r|^2
   * rvro <- rvr
   * rvr <- rv_adjoint.r
   * b <- (rvv/rvro) * (a/omega)
   * p <- p - omega*v
   * p <- r + b*p
   */

  rsq = *rsq_ptr;
  rvr = dcmplx(rsq,(double)0.0);
  iteration = 0;
  do{
    rvv = dcmplx((double)0.0,(double)0.0);


    delta0( F_OFFSET(p), F_OFFSET(v), PLUS );


    FORALLSITES(i,s){
      ctmp =  wvec_dot( &(s->rv), &(s->v));  /* + */
      CSUM(rvv,ctmp);
    }

    g_dcomplexsum(&rvv);
/*
CMUL_J(rvv,rvv,test);
if(test.real<1.e-24){
if(this_node==0) printf("zero divide being patched\n");
rvv.real += 1.e-12;
}
*/
    CDIV(rvr,rvv,a);
    CMULREAL(a,-1.0,ctmp);
    FORALLSITES(i,s){
      c_scalar_mult_add_wvec( &(s->r), &(s->v), &ctmp, &(s->sss) );
    }

    tdots = dcmplx((double)0.0,(double)0.0);
    ttsq = 0.0;
    delta0( F_OFFSET(sss), F_OFFSET(ttt), PLUS );
    FORALLSITES(i,s){
      ctmp = wvec_dot( &(s->ttt), &(s->sss));
      CSUM(tdots, ctmp);
      ttsq += (double)magsq_wvec( &(s->ttt) );


    }
    g_doublesum( &ttsq );
    g_dcomplexsum( &tdots );
    iteration++;

    omega.real = tdots.real/ttsq;
    omega.imag = tdots.imag/ttsq;
    omegam.real = -omega.real;
    omegam.imag = -omega.imag;
    rsq = 0.0;
    rvro = rvr;
    rvr = dcmplx((double)0.0,(double)0.0);
    FORALLSITES(i,s){
      c_scalar_mult_add_wvec( (wilson_vector *)F_PT(s,sol),
               &(s->sss), &omega,(wilson_vector *) F_PT(s,sol) );
      c_scalar_mult_add_wvec( (wilson_vector *)F_PT(s,sol),
               &(s->p),  &a, (wilson_vector *)F_PT(s,sol) );
      c_scalar_mult_add_wvec( &(s->sss), &(s->ttt), &omegam ,&(s->r) );
      ctmp = wvec_dot( &(s->rv),&(s->r));
      CSUM(rvr, ctmp);
      rsq += (double)magsq_wvec( &(s->r) );
    }
    g_doublesum( &rsq );
    g_dcomplexsum( &rvr);
    if(this_node==0 && ((iteration/5)*5==iteration) ) {
      printf("BiCongrad: iter %d, rsq %e, rvr %e, %e\n",
	     iteration,(double)rsq, (double)rvr.real, (double)rvr.imag );
      fflush(stdout);
   } 

    if( rsq <= rsqstop ){
      *rsq_ptr= rsq;
      return(iteration);
    }

/**    dtime += dclock();
    if(this_node==0) {
      printf("BiCongrad: time = %e iters = %d mflops = %e\n",
	     dtime,iteration,(double)(2840.0*volume*iteration/
				      (1.0e6*dtime*numnodes())) );
    } **/
    CDIV(rvr,rvro,b);
    CMUL(a,b,ctmp);
    CDIV(ctmp,omega,b);
    FORALLSITES(i,s){
      c_scalar_mult_add_wvec( &(s->p), &(s->v), &omegam, &(s->p) );
      c_scalar_mult_add_wvec( &(s->r), &(s->p), &b, &(s->p) );
    }

  } while( iteration%niterr != 0);

  *rsq_ptr = rsq;
  return(iteration);

} /* BICONGRAD_ITER */
