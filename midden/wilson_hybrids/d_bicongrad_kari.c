/************************ d_bicongrad_kari.c ************************/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE! */

Received: by higgs.physics.arizona.edu (5.57/Ultrix3.0-C)
	id AA05193; Wed, 25 May 94 15:11:14 MST
Date: Wed, 25 May 94 15:11:14 MST
From: leo@higgs.physics.arizona.edu (Leo Karkkainen)
Message-Id: <9405252211.AA05193@higgs.physics.arizona.edu>
To: doug@klingon.physics.Arizona.EDU

/* Terve Leksa:
 * Loysin koodista viela muutamia virheita; esim. tagin indeksi (0/1)
 * olivat jossain vaarin.  Tein aika ison muutoksen koodiin,
 * pistin sen iteroinnin erilliseen aliohjelmaan, jota sitten
 * kutsutaan kahdesti.  Otin alla olevasta ohjelmasta myos pois
 * sen 'restartin', jos ei suppene niter-iteraatiossa.  Tama
 * kenties pitaa pistaa takaisin, kuitenkaan se ei ole ihan
 * suoraviivaista: jos eka iteraatio ei suppene, niin sita
 * ei voi noin vaan aloittaa alusta, silla arvaus tulee ihan
 * samaksi!  Pitaisi kayttaa toista uuden arvauksen tekemiseen, tai
 * tehda erillinen arvausosa.  Jos taas eka on ok, mutta toka ei
 * suppene, niin silloin tietysti restartataan vain toinen.  Mutta
 * jatin sen nyt pois, jos se ei ole ihan valttamaton.
 *
 * Niin, ja jos alkuarvaus on niin hyva ettei ekaan luuppiin tarvitse
 * menna, niin ei tarvitse toiseenkaan.
 *
 * En ole testannut ohjelmaa -- varmana on virheita!
 *
 * Kari
 */

/*
 * Terve Kari
 * Loysin muuten 'viimeisen' bugin bi_congradista. Nyt se nayttaa
 * suppenevan. Laskeskelin yhdessa kohtaa vaarien vektorien pistetuloa -
 * sen sita saa kun maccimaisesti kopsaa lauseita ja unohtaa sitten vaihtaa
 * argumentit...Koodi on tassa perassa.
 */
/******* bicongrad2.c - biconnugate gradient for SU3/fermions ****/
/* MIMD version 6 */
/* Wilson fermions */
/* Currently only for LU preconditioning */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
   the fermion spinors live on even sites only.  In other words, if
   Dslash_oe is the dslash operator with its source on even sites and
   its result on odd sites, etc.:
   with LU:
   M = 1 - kappa^2 * Dslash_eo * Dslash_oe
*/
#ifdef LU
#define FORMYSITES FOREVENSITES
#else
#define FORMYSITES FORALLSITES
#endif

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "chi", and the initial guess and answer
   in "psi".  "r" is the residual vector and "rv" the conjugate residual,
   "p","mp","ttt","sss" are working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
*/
#include <stdio.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include <lattice.h>
#include "../include/macros.h"
#include "../include/comdefs.h"

void dslash_special();

#ifdef PROTO
int bicongrad_iter(field_offset sol,int niter,double rsqstop,double *rsq_ptr,
		   int plmin,msg_tag *tag[8],msg_tag *tag2[8]);
#endif

int congrad(niter,rsqmin,final_rsq_ptr) 
int niter; Real rsqmin,*final_rsq_ptr;
{
  register int i;
  register site *s;
  int good_enough;
  int iteration; /* counter for iterations */
  complex a;
  double rsq1,rsq2,source_norm,rsqstop;
  msg_tag *tag[8],*tag2[8];
  double dtime,dclock();

  dtime= -dclock();
  iteration = 0;
  good_enough = 0;

  /* mp  <-  M*psi           initial guess for the first inversion
   * ttt <-  M_adjoint*mp    
   * r,rv,p <- chi - ttt
   * rsq1 = rvr = |r|^2
   * source_norm = |chi|^2
   * rvr = rsq 
   */

  rsq1 = rsq2 = source_norm = 0.0;
#ifdef LU
  dslash_special(F_OFFSET(psi),F_OFFSET(psi) ,PLUS,ODD,tag,0 );
  dslash_special(F_OFFSET(psi),F_OFFSET(mp),PLUS,EVEN,tag2,0 );
  FOREVENSITES(i,s){
    scalar_mult_add_wvec( &(s->psi),&(s->mp), -kappa*kappa,&(s->mp) );
  }
  dslash_special(F_OFFSET(mp),F_OFFSET(mp),MINUS,ODD,tag,1 );
  dslash_special(F_OFFSET(mp),F_OFFSET(ttt),MINUS,EVEN,tag2,1 );
  FOREVENSITES(i,s){
    scalar_mult_add_wvec( &(s->mp),&(s->ttt), -kappa*kappa, &(s->ttt) );
    sub_wilson_vector( &(s->chi), &(s->ttt), &(s->r) );
    s->rv = s->p = s->r;
    source_norm += (double)magsq_wvec( &(s->chi) );
    rsq1 += (double)magsq_wvec( &(s->r) );
  }
#endif
  g_doublesum( &source_norm );
  g_doublesum( &rsq1 );

  iteration++ ;	/* iteration counts number of multiplications
		   both by M and M_adjoint */
  total_iters++;

/**  if(this_node==0){
    printf("bicongrad_1: source_norm = %e\n",source_norm);
    fflush(stdout);
  } **/

  rsqstop = rsqmin * source_norm;
  if( rsq1 > rsqstop ){
    /* Now have to iterate, not good enough guess! */

    /* first, solve eq. M_adjoint mp = chi */

    i = bicongrad_iter(F_OFFSET(mp),niter,rsqstop,&rsq1,MINUS,tag,tag2);
    iteration += i;
    total_iters += i;

/**    if (this_node==0) {
      printf("bicongrad_1: iter %d, rsq %g\n",iteration,rsq1);
      fflush(stdout);
    } **/

    /* then, solve M psi = mp */
    /* mp  <-  solution of the first iteration
     * ttt <-  M*psi
     * r,rv,p <- mp - ttt
     * rsq = |r|^2
     * source_norm = |mp|^2
     * rvr = rsq
     */

    rsq2 = source_norm = 0.0;
#ifdef LU
    /* the first part was done above */
    /* dslash_special(F_OFFSET(psi),F_OFFSET(psi) ,PLUS,ODD,tag,1 );*/
    dslash_special(F_OFFSET(psi),F_OFFSET(ttt),PLUS,EVEN,tag2,1 );
    FOREVENSITES(i,s){
      scalar_mult_add_wvec( &(s->psi),&(s->ttt), -kappa*kappa,&(s->ttt));
      sub_wilson_vector( &(s->mp), &(s->ttt), &(s->r) );
      s->rv = s->p = s->r;
      source_norm += (double)magsq_wvec( &(s->mp) );
      rsq2 += (double)magsq_wvec( &(s->r) );
    }
#endif
    g_doublesum( &source_norm );
    g_doublesum( &rsq2 );
    iteration++ ;   /* iteration counts number of multiplications
		       by M_adjoint*M */
    total_iters++;
/**    if(this_node==0) {
      printf("bicongrad_2: source_norm = %e\n",source_norm);
      fflush(stdout);
    } **/
    rsqstop = rsqmin * source_norm;
    if( rsq2 > rsqstop ){
      /* now iterate M psi = mp */

      i = bicongrad_iter(F_OFFSET(psi),niter,rsqstop,&rsq2,PLUS,tag,tag2);
      iteration += i;
      total_iters += i;

/**      if (this_node==0) {
	printf("bicongrad_2: iter %d, rsq %g\n",i,rsq2);
	fflush(stdout);
} **/

    } /* rsq2 > rsqstop2 */

  } /* rsq1 > rsqstop1 */

  for( i=XUP; i <= TUP; i++) {
    cleanup_gather(tag[i]);
    cleanup_gather(tag[OPP_DIR(i)]);
#ifdef LU
    cleanup_gather(tag2[i]);
    cleanup_gather(tag2[OPP_DIR(i)]);
#endif
  }

  *final_rsq_ptr = (Real)(rsq1 + rsq2);   /* ?? */

  return(iteration);
}


/* 
 * Subroutine bicongrad_iter does the actual convergence loop
 * NOTE - this assumes already started gathers, tag and tag2
 * plmin = PLUS: solves M * sol = x, MINUS: M_adjoint * sol = x
 */

int bicongrad_iter(sol,niter,rsqstop,rsq_ptr,plmin,tag,tag2)
  field_offset sol; int niter, plmin; double rsqstop,*rsq_ptr;
  msg_tag *tag[8],*tag2[8];
{
  register int i;
  register site *s;
  int iteration;	/* counter for iterations */
  complex a;
  double rsq,ttsq;	/* Sugar's a,b,resid**2,previous resid*2 */
                        /* pkp = cg_p.K.cg_p */
  double dtime,dclock();
  complex tdots,rvro,rvv,rvr,ctmp;
  complex cc,omega,omegam,b;

  /* main loop for inversion - do until convergence or time to restart */
  /* 
   * v <- M * p                  if PLUS, M_adjoint if MINUS
   * rvv <- rv_adjoint.v
   * a <- rvr/rvv
   * sss <- r - a*v
   * ttt <- M * sss              same thing
   * tdots <- ttt_adjoint.sss
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
  rvr = cmplx((Real)rsq,0.0);
  iteration = 0;
  do{
    rvv = cmplx(0.0,0.0);
#ifdef LU
    dslash_special(F_OFFSET(p),F_OFFSET(p),plmin, ODD,tag, 1);
    dslash_special(F_OFFSET(p),F_OFFSET(v),plmin,EVEN,tag2,1);
    FOREVENSITES(i,s){
      scalar_mult_add_wvec( &(s->p), &(s->v), -kappa*kappa, &(s->v) );
      ctmp =  wvec_dot( &(s->rv), &(s->v));  /* + */
      CSUM(rvv,ctmp);
    }
    g_complexsum(&rvv);
    CDIV(rvr,rvv,a);
    CMULREAL(a,-1.0,ctmp);
    FOREVENSITES(i,s){
      c_scalar_mult_add_wvec( &(s->r), &(s->v), &ctmp, &(s->sss) );
    }
    dslash_special(F_OFFSET(sss),F_OFFSET(sss),plmin,ODD, tag, 1);
    dslash_special(F_OFFSET(sss),F_OFFSET(ttt),plmin,EVEN,tag2,1);
    tdots = cmplx(0.0,0.0);
    ttsq = 0.0;
    FOREVENSITES(i,s){
      scalar_mult_add_wvec(&(s->sss),&(s->ttt), -kappa*kappa, &(s->ttt) );
      ctmp = wvec_dot( &(s->ttt),&(s->sss));
      CSUM(tdots, ctmp);
      ttsq += (double)magsq_wvec( &(s->ttt) );
    }
#endif
    g_doublesum( &ttsq );
    g_complexsum( &tdots );
    iteration++;

    omega.real = tdots.real/ttsq;
    omega.imag = tdots.imag/ttsq;
    omegam.real = -omega.real;
    omegam.imag = -omega.imag;
    rsq = 0.0;
    rvro = rvr;
    rvr = cmplx(0.0,0.0);
    FORMYSITES(i,s){
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
    g_complexsum( &rvr);
/**    if(this_node==0) {
      printf("BiCongrad: iter %d, rsq %e, rvr %e, %e\n",
	     iteration,(double)rsq, (double)rvr.real, (double)rvr.imag );
      fflush(stdout);
   } **/

    if( rsq <= rsqstop ){
      *rsq_ptr= rsq;
      return(iteration);
    }

/**    dtime += dclock();
    if(this_node==0) {
      printf("BICONGRAD: time = %e iters = %d mflops = %e\n",
	     dtime,iteration,(double)(2840.0*volume*iteration/
				      (1.0e6*dtime*numnodes())) );
    } **/
    CDIV(rvr,rvro,b);
    CMUL(a,b,cc);
    CDIV(cc,omega,b);
    FORMYSITES(i,s){
      c_scalar_mult_add_wvec( &(s->p), &(s->v), &omegam, &(s->p) );
      c_scalar_mult_add_wvec( &(s->r), &(s->p), &b, &(s->p) );
    }

  } while( iteration%niter != 0);

  *rsq_ptr = rsq;
  return(iteration);

}

