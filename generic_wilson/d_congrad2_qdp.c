/******* d_congrad2.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* Wilson fermions */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
   the fermion spinors live on even sites only.  In other words, if
   Dslash_oe is the dslash operator with its source on even sites and
   its result on odd sites, etc.:
 
   without LU:
   M = 1 - kappa*( Dslash_eo + DSLASH_oe )
   with LU:
   M = 1 - kappa^2 * Dslash_eo * Dslash_oe
*/
#ifdef LU
#define FORMYSITES FOREVENSITES
#define MYSUBSET QDP_even
#else
#define FORMYSITES FORALLSITES
#define MYSUBSET QDP_all
#endif

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "chi", and the initial guess and answer
   in "phi".  "r" is the residual vector, and "p" and "mp" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
*/

#include "generic_wilson_includes.h"
#include <qdp.h>
#include "../include/generic_qdp.h"

static QDP_ColorMatrix *gaugelink[8];
//static QDP_HalfFermion *dtemp0, *dtemp1[8], *temp1[8], *temp2[8];
static QDP_DiracFermion *temp1[8], *temp2[8], *temp3[8], *temp4[8];
static QDP_DiracFermion *psi, *chi, *cgp, *cgr, *mp, *ttt, *tt1, *tt2, *t1, *t2, *t3;

#define PRESHIFT_LINKS
#define SHIFT_D

static void
setup_cg(void)
{
  static int is_setup=0;
  if(!is_setup) {
    int i;
    is_setup = 1;
    psi = QDP_create_D();
    chi = QDP_create_D();
    cgp = QDP_create_D();
    cgr = QDP_create_D();
    mp = QDP_create_D();
    ttt = QDP_create_D();
    tt1 = QDP_create_D();
    tt2 = QDP_create_D();
    t1 = QDP_create_D();
    t2 = QDP_create_D();
    t3 = QDP_create_D();
    //dtemp0 = QDP_create_H();
    for(i=0; i<4; i++) {
#ifndef PRESHIFT_LINKS
      gaugelink[i] = QDP_create_M();
#endif
    }
    for(i=0; i<8; i++) {
#ifdef PRESHIFT_LINKS
      gaugelink[i] = QDP_create_M();
#endif
      //dtemp1[i] = QDP_create_H();
      //temp1[i] = QDP_create_H();
      //temp2[i] = QDP_create_H();
      temp1[i] = QDP_create_D();
      temp2[i] = QDP_create_D();
      temp3[i] = QDP_create_D();
      temp4[i] = QDP_create_D();
    }
  }
}

#if 0
void
dslash_special_qdp(QDP_DiracFermion *dest, QDP_DiracFermion *src,
		   int sign, QDP_Subset subset, QDP_HalfFermion *temp[]);
#endif
void
dslash_special_qdp(QDP_DiracFermion *dest, QDP_DiracFermion *src,
		   int sign, QDP_Subset subset, QDP_DiracFermion *temp[]);

int
congrad_w(int niter, Real rsqmin, Real *final_rsq_ptr) 
{
  int i;
  int iteration;	/* counter for iterations */
  double source_norm;
  double rsqstop;
  QLA_Real a, b;
  double rsq,oldrsq,pkp;	/* Sugar's a,b,resid**2,previous resid*2 */
				/* pkp = cg_p.K.cg_p */
  QLA_Real mkappa;
  QLA_Real sum;
#ifdef CGTIME
  double dtime;
#endif
#ifdef LU
  mkappa = -kappa*kappa;
#else
  mkappa = -kappa;
#endif

  setup_cg();

  for(i=0; i<4; i++) {
    set_M_from_site(gaugelink[i], F_OFFSET(link[i]),EVENANDODD);
  }
  set_D_from_site(psi, F_OFFSET(psi),EVENANDODD);
  set_D_from_site(chi, F_OFFSET(chi),EVENANDODD);

#ifdef PRESHIFT_LINKS
  {
    QDP_ColorMatrix *tcm;
    tcm = QDP_create_M();
    for(i=0; i<4; i++) {
      QDP_M_eq_sM(tcm, gaugelink[i], QDP_neighbor[i], QDP_backward, QDP_all);
      QDP_M_eq_Ma(gaugelink[i+4], tcm, QDP_all);
    }
    QDP_destroy_M(tcm);
  }
#endif

#ifdef CGTIME
  dtime = -dclock();
#endif

  iteration=0;
 start:
  /* mp <-  M_adjoint*M*psi
     r,p <- chi - mp
     rsq = |r|^2
     source_norm = |chi|^2
  */
  rsq = source_norm = 0.0;

#ifdef LU

  QDP_D_eq_D(cgp, psi, QDP_even);
  dslash_special_qdp(tt1, cgp, 1, QDP_odd, temp1);
  dslash_special_qdp(ttt, tt1, 1, QDP_even, temp2);
  QDP_D_eq_r_times_D_plus_D(ttt, &mkappa, ttt, cgp, QDP_even);

  dslash_special_qdp(tt2, ttt, -1, QDP_odd, temp3);
  dslash_special_qdp(mp, tt2, -1, QDP_even, temp4);
  QDP_D_eq_r_times_D_plus_D(mp, &mkappa, mp, ttt, QDP_even);
  QDP_D_eq_D_minus_D(cgr, chi, mp, QDP_even);
  QDP_D_eq_D(cgp, cgr, QDP_even);

  QDP_r_eq_norm2_D(&sum, chi, QDP_even);
  source_norm = sum;
  QDP_r_eq_norm2_D(&sum, cgr, QDP_even);
  rsq = sum;

#else

  QDP_D_eq_D(cgp, psi, QDP_even);
  dslash_special_qdp(ttt, cgp, 1, QDP_all, temp1);
  QDP_D_eq_r_times_D_plus_D(ttt, &mkappa, ttt, cgp, QDP_all);

  dslash_special_qdp(mp, ttt, -1, QDP_all, temp1);
  QDP_D_eq_r_times_D_plus_D(mp, &mkappa, mp, ttt, QDP_all);

  QDP_D_eq_D_minus_D(cgr, chi, mp, QDP_all);
  QDP_D_eq_D(cgp, cgr, QDP_all);

  QDP_r_eq_norm2_D(&sum, chi, QDP_all);
  source_norm = sum;
  QDP_r_eq_norm2_D(&sum, cgr, QDP_all);
  rsq = sum;

#endif

  iteration++ ;	/* iteration counts number of multiplications
		   by M_adjoint*M */
  total_iters++;
  /**if(this_node==0)printf("congrad2: source_norm = %e\n",source_norm);
     if(this_node==0)printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
     iteration,(double)rsq,(double)pkp,(double)a );**/
  rsqstop = rsqmin * source_norm;
  if( rsq <= rsqstop ){
    *final_rsq_ptr= (Real)rsq;
    return (iteration);
  }

  /* main loop - do until convergence or time to restart */
  /* 
     oldrsq <- rsq
     mp <- M_adjoint*M*p
     pkp <- p.M_adjoint*M.p
     a <- rsq/pkp
     psi <- psi + a*p
     r <- r - a*mp
     rsq <- |r|^2
     b <- rsq/oldrsq
     p <- r + b*p
  */
  do {
    oldrsq = rsq;
#ifdef LU
    dslash_special_qdp(tt1, cgp, 1, QDP_odd, temp1);
    dslash_special_qdp(ttt, tt1, 1, QDP_even, temp2);
    QDP_D_eq_r_times_D_plus_D(ttt, &mkappa, ttt, cgp, QDP_even);

    dslash_special_qdp(tt2, ttt, -1, QDP_odd, temp3);
    dslash_special_qdp(mp, tt2, -1, QDP_even, temp4);
    QDP_D_eq_r_times_D_plus_D(mp, &mkappa, mp, ttt, QDP_even);

    QDP_r_eq_re_D_dot_D(&sum, cgp, mp, QDP_even);
    pkp = sum;
#else
    dslash_special_qdp(ttt, cgp, 1, QDP_all, temp1);
    QDP_D_eq_r_times_D_plus_D(ttt, &mkappa, ttt, cgp, QDP_all);

    dslash_special_qdp(mp, ttt, -1, QDP_all, temp1);
    QDP_D_eq_r_times_D_plus_D(mp, &mkappa, mp, ttt, QDP_all);

    QDP_r_eq_re_D_dot_D(&sum, cgp, mp, QDP_all);
    pkp = sum;
#endif
    iteration++;
    total_iters++;

    a = rsq / pkp;
    QDP_D_peq_r_times_D(psi, &a, cgp, MYSUBSET);
    QDP_D_meq_r_times_D(cgr, &a, mp, MYSUBSET);
    QDP_r_eq_norm2_D(&sum, cgr, MYSUBSET);
    rsq = sum;

    /**if(this_node==0)printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
       iteration,(double)rsq,(double)pkp,(double)a );**/
    if( rsq <= rsqstop ){
      *final_rsq_ptr= (Real)rsq;
#ifdef CGTIME
      dtime += dclock();
      if(this_node==0)
	printf("CONGRAD2: time = %.2e size_r= %.2e iters= %d MF = %.1f\n",
	       dtime,rsq,iteration,
	       (double)6480*iteration*even_sites_on_node/(dtime*1e6));
      //(double)5616*iteration*even_sites_on_node/(dtime*1e6));
#endif
      set_site_from_D(F_OFFSET(psi), psi,EVENANDODD);
      return (iteration);
    }

    b = rsq / oldrsq;
    QDP_D_eq_r_times_D_plus_D(cgp, &b, cgp, cgr, MYSUBSET);

  } while( iteration%niter != 0);

  set_site_from_D(F_OFFSET(psi), psi,EVENANDODD);

  if( iteration < 3*niter ) goto start;
  *final_rsq_ptr= (Real)rsq;
  return(iteration);
}

/************ dslash *************/
#ifdef PRESHIFT_LINKS

/* Special dslash for use by congrad.  Uses restart_gather_site() when
   possible. Last argument is an integer, which will tell if
   gathers have been started.  If is_started=0,use
   start_gather_site, otherwise use restart_gather_site.
   Argument "tag" is a vector of a msg_tag *'s to use for
   the gathers.
   The calling program must clean up the gathers! */
void
dslash_special_qdp(QDP_DiracFermion *dest, QDP_DiracFermion *src,
		   int sign, QDP_Subset subset, QDP_DiracFermion *temp[])
{
  int mu;
  QDP_DiracFermion *vsrc[8];
  QDP_DiracFermion *vdest[8];
  QDP_Shift sh[8];
  QDP_ShiftDir sd[8];
  int dir[8], sgn[8];

#ifndef SHIFT_D
  QDP_Subset othersubset;

  if(subset==QDP_even) othersubset = QDP_odd;
  else if(subset==QDP_odd) othersubset = QDP_even;
  else othersubset = QDP_all;
#endif

  for(mu=0; mu<4; mu++) {
    vsrc[mu] = src;
    vsrc[mu+4] = src;
    vdest[mu] = dest;
    vdest[mu+4] = dest;
    dir[mu] = mu;
    dir[mu+4] = mu;
    sgn[mu] = sign;
    sgn[mu+4] = -sign;
    sh[mu] = QDP_neighbor[mu];
    sh[mu+4] = QDP_neighbor[mu];
    sd[mu] = QDP_forward;
    sd[mu+4] = QDP_backward;
  }

  /* Take Wilson projection for src displaced in up direction, gather
     it to "our site" */

#ifdef SHIFT_D
  QDP_D_veq_sD(temp, vsrc, sh, sd, subset, 8);
#else
  for(mu=0; mu<8; mu++) {
    QDP_H_eq_spproj_D(dtemp1[mu], vsrc[mu], dir[mu], sgn[mu], othersubset);
    QDP_H_eq_sH(temp[mu], dtemp1[mu], sh[mu], sd[mu], subset);
  }
#endif
  //QDP_H_veq_spproj_D(dtemp1, vsrc, dir, sgn, othersubset, 8);
  //QDP_H_veq_sH(temp, dtemp1, sh, sd, subset, 8);

  //QDP_H_veq_spproj_D(dtemp1, vsrc, dir, sgn, othersubset, 4);
  //QDP_H_veq_sH(temp, dtemp1, QDP_neighbor, sd, subset, 4);
  //for(mu=0; mu<4; mu++) {
  //QDP_H_eq_spproj_D(dtemp1[mu], src, mu, sign, othersubset);
  //QDP_H_eq_sH(temp[mu], dtemp1[mu], QDP_neighbor[mu], QDP_forward, subset);
  //}
  //QDP_H_veq_spproj_D(dtemp1+4, vsrc, dir, sgn+4, othersubset, 4);
  //QDP_H_veq_sH(temp+4, dtemp1+4, QDP_neighbor, sd+4, subset, 4);
  //for(mu=0; mu<4; mu++) {
  //QDP_H_eq_spproj_D(dtemp1[mu+4], src, mu, -sign, othersubset);
  //QDP_H_eq_sH(temp[mu+4], dtemp1[mu+4], QDP_neighbor[mu], QDP_backward, subset);
  //}

  /* Set dest to zero */
  /* Take Wilson projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add.
     to dest */

  QDP_D_eq_zero(dest, subset);
  //QDP_D_vpeq_sprecon_M_times_H(vdest, gaugelink, temp, dir, sgn, subset, 8);
  QDP_D_vpeq_wilsonspin_M_times_D(vdest, gaugelink, temp, dir, sgn, subset, 8);
#if 0
  for(mu=0; mu<8; mu++) {
    QDP_D_peq_wilsonspin_M_times_D(dest, gaugelink[mu], temp[mu], dir[mu], sgn[mu], subset);
  }
#endif

#if 0
  {
    QDP_HalfFermion *hf;
    hf = QDP_create_H();
    for(mu=0; mu<8; mu++) {
      QDP_D_peq_wilsonspin_M_times_D(dest, gaugelink[mu], temp[mu], dir[mu], sgn[mu], subset);
    }
    for( ; mu<8; mu++) {
      QDP_H_eq_spproj_D(hf, temp[mu], dir[mu], sgn[mu], subset);
      QDP_D_peq_sprecon_M_times_H(dest, gaugelink[mu], hf, dir[mu], sgn[mu], subset);
    }
    QDP_destroy_H(hf);
  }
#endif

  for(mu=0; mu<8; mu++) {
    QDP_discard_D(temp[mu]);
  }
} /* end of dslash_special_qdp() */

#else

/* Special dslash for use by congrad.  Uses restart_gather_site() when
   possible. Last argument is an integer, which will tell if
   gathers have been started.  If is_started=0,use
   start_gather_site, otherwise use restart_gather_site.
   Argument "tag" is a vector of a msg_tag *'s to use for
   the gathers.
   The calling program must clean up the gathers! */
void
dslash_special_qdp(QDP_DiracFermion *dest, QDP_DiracFermion *src,
		   int sign, QDP_Subset subset, QDP_HalfFermion *temp[])
{
  int mu;
  QDP_DiracFermion *vsrc[4];
  QDP_DiracFermion *vdest[4];
  QDP_ShiftDir fwd[4], bck[4];
  int dir[4], sgn[4], msgn[4];
  QDP_Subset othersubset;

  for(mu=0; mu<4; mu++) {
    vsrc[mu] = src;
    vdest[mu] = dest;
    fwd[mu] = QDP_forward;
    bck[mu] = QDP_backward;
    dir[mu] = mu;
    sgn[mu] = sign;
    msgn[mu] = -sign;
  }

  if(subset==QDP_even) othersubset = QDP_odd;
  else if(subset==QDP_odd) othersubset = QDP_even;
  else othersubset = QDP_all;

  /* Take Wilson projection for src displaced in up direction, gather
     it to "our site" */

  for(mu=0; mu<4; mu++) {
    QDP_H_eq_spproj_D(dtemp1[mu], src, mu, sign, othersubset);
    QDP_H_eq_sH(temp[mu], dtemp1[mu], QDP_neighbor[mu], QDP_forward, subset);
  }

  /* Take Wilson projection for src displaced in down direction,
     multiply it by adjoint link matrix, gather it "up" */

  for(mu=0; mu<4; mu++) {
    QDP_H_eq_spproj_D(dtemp0, src, mu, -sign, othersubset);
    QDP_H_eq_Ma_times_H(dtemp1[4+mu], gaugelink[mu], dtemp0, othersubset);
    QDP_H_eq_sH(temp[4+mu], dtemp1[4+mu], QDP_neighbor[mu], QDP_backward, subset);
  }

  /* Set dest to zero */
  /* Take Wilson projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add.
     to dest */

  QDP_D_eq_zero(dest, subset);

  QDP_D_vpeq_sprecon_M_times_H(vdest, gaugelink, temp, dir, sgn, subset, 4);
  for(mu=0; mu<4; mu++) {
    //QDP_D_peq_sprecon_M_times_H(dest, gaugelink[mu], temp[mu], mu, sign, subset);
    QDP_discard_H(temp[mu]);
  }

  /* Take Wilson projection for src displaced in down direction,
     expand it, and add to dest */

  for(mu=0; mu<4; mu++) {
    QDP_D_peq_sprecon_H(dest, temp[4+mu], mu, -sign, subset);
    QDP_discard_H(temp[4+mu]);
  }

} /* end of dslash_special_qdp() */
#endif
