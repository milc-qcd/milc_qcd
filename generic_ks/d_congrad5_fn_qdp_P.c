/******* d_congrad5_fn_qdp_P.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC standalone Level 2 QDP inverter for FN KS actions */

/* NOTE: This code is actually an include file for d_congrad5_fn_qdp_F.c
   and d_congrad5_fn_qdp_D.c, so any edits should be consistent with this
   purpose. */

/* Entry points (must be redefined to precision-specific names)

   KS_CONGRAD_QDP
   KS_CONGRAD_MILCFIELD2QDP
   KS_CONGRAD_MILC2QDP

*/

/* Kogut-Susskind fermions -- this version for "fat plus Naik" quark
   actions.  */

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* 10/02/01 C. DeTar Consolidated with tmp version */

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "src", and the initial guess and answer
   in "dest".  "resid" is the residual vector, and "cg_p" and "ttt" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq <= rsqmin*source_norm.
	This is different than our old definition of the stopping
	criterion.  To convert an old stopping residual to the new
	one, multiply the old one by sqrt( (2/3)/(8+2*m) )
        This is because the source is obtained from
        a random vector with average squared magnitude 3 on each site.
        Then, on 1/2 the sites, we gather and sum the eight neighboring
        random vectors and add 2*m times the local vector.
            source = M_adjoint*R, on even sites
   reinitialize after niters iterations and try once more.
   parity=EVEN = do only even sites, parity=ODD = do odd sites,
   parity=EVENANDODD = do all sites
*/

#if ( QDP_Precision == 'F' )

#define KS_CONGRAD_QDP       ks_congrad_qdp_F
#define KS_CONGRAD_MILCFIELD2QDP  ks_congrad_milcfield2qdp_F
#define KS_CONGRAD_MILC2QDP  ks_congrad_milc2qdp_F
#define SETUP_DSLASH         setup_dslash_F
#define DSLASH_QDP_FN_SPECIAL2 dslash_qdp_F_fn_special2
#define DSLASH_QDP_FN        dslash_qdp_F_fn
#define BCKLINK              bcklink_F
#define FATLINKS          fatlinks_F
#define LONGLINKS         longlinks_F
#define IMPLINKS          implinks_F

#else

#define KS_CONGRAD_QDP       ks_congrad_qdp_D
#define KS_CONGRAD_MILCFIELD2QDP  ks_congrad_milcfield2qdp_D
#define KS_CONGRAD_MILC2QDP  ks_congrad_milc2qdp_D
#define SETUP_DSLASH         setup_dslash_D
#define DSLASH_QDP_FN_SPECIAL2 dslash_qdp_D_fn_special2
#define DSLASH_QDP_FN        dslash_qdp_D_fn
#define BCKLINK              bcklink_D
#define FATLINKS          fatlinks_D
#define LONGLINKS         longlinks_D
#define IMPLINKS          implinks_D

#endif

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qdp.h"
#include "../include/generic_ks_qdp.h"
#include <lattice_qdp.h>

static int congrad_setup=0;
static QDP_ColorVector *ttt, *tttt, *resid, *cg_p;
static QDP_ColorVector *temp1[16], *temp2[16];
extern QDP_ColorMatrix *BCKLINK[8];
extern QDP_ColorMatrix **FATLINKS, **LONGLINKS, *IMPLINKS[8];

static void
setup_congrad()
{
  int i;

  if(congrad_setup)return;
  SETUP_DSLASH();
  ttt = QDP_create_V();
  tttt = QDP_create_V();
  resid = QDP_create_V();
  cg_p = QDP_create_V();
  for(i=0; i<16; i++) {
    temp1[i] = QDP_create_V();
    temp2[i] = QDP_create_V();
  }
  congrad_setup = 1;
}

int
KS_CONGRAD_QDP(QDP_ColorVector *src, QDP_ColorVector *dest, 
	       quark_invert_control *qic, QLA_Real mass,
	       ferm_links_t *fn)
{
  QLA_Real a,b;		 /* Sugar's a,b,resid**2,last resid*2 */
  QLA_Real rsq,oldrsq,pkp; /* pkp = cg_p.K.cg_p */
  QLA_Real msq_x4;	 /* 4*mass*mass */
  QLA_Real source_norm;	 /* squared magnitude of source vector */
  QLA_Real rsqstop;	 /* stopping residual normalized by source norm */
  QDP_Subset q_parity, q_otherparity;
  int iteration;	 /* counter for iterations */
  int niter = qic->max;
  int nrestart = qic->nrestart;
  QLA_Real rsqmin = qic->resid * qic->resid;
  su3_matrix *t_fatlink;
  su3_matrix *t_longlink;

/* Timing */
  double dtimec;
  double remaptime;
#ifdef CGTIME
  double nflop = 1187;
  if(qic->parity==EVENANDODD) nflop *= 2;
#endif

  setup_congrad();

  if(qic->parity==ODD) {
    q_parity = QDP_odd;
    q_otherparity = QDP_even;
  } else if(qic->parity==EVEN) {
    q_parity = QDP_even;
    q_otherparity = QDP_odd;
  } else {
    q_parity = QDP_all;
    q_otherparity = QDP_all;
  }
  msq_x4 = 4.0*mass*mass;
  iteration = 0;

  t_longlink = fn->fl.lng;
  t_fatlink = fn->fl.fat;

  remaptime = -dclock();
  set4_M_from_field(FATLINKS, t_fatlink, EVENANDODD);
  set4_M_from_field(LONGLINKS, t_longlink, EVENANDODD);

  {
    QDP_ColorMatrix *tcm;
    int i;
    tcm = QDP_create_M();
    for(i=0; i<8; ++i) {
      QDP_M_eq_sM(tcm, IMPLINKS[i], shiftdirs[i], QDP_backward, QDP_all);
      QDP_M_eqm_Ma(BCKLINK[i], tcm, QDP_all);
    }
    QDP_destroy_M(tcm); tcm = NULL;
  }

  remaptime += dclock();
  dtimec = -dclock(); 

  /* initialization process */
  do {
    /* ttt <-  (-1)*M_adjoint*M*dest
       resid,cg_p <- src + ttt
       rsq = |resid|^2
       source_norm = |src|^2
    */

    QDP_V_eq_V(cg_p, dest, q_parity);
    DSLASH_QDP_FN_SPECIAL2(cg_p, tttt, q_otherparity, temp1);
    DSLASH_QDP_FN_SPECIAL2(tttt, ttt, q_parity, temp2);
    iteration++;    /* iteration counts multiplications by M_adjoint*M */

    /* ttt  <- ttt - msq_x4*src	(msq = mass squared) */
    QDP_V_meq_r_times_V(ttt, &msq_x4, dest, q_parity);
    QDP_V_eq_V_plus_V(resid, src, ttt, q_parity);
    QDP_r_eq_norm2_V(&source_norm, src, q_parity);
    QDP_r_eq_norm2_V(&rsq, resid, q_parity);
    oldrsq = rsq;
    QDP_V_eq_zero(cg_p, q_parity);
    rsqstop = rsqmin * source_norm;

    /* main loop - do until convergence or time to restart */
    /*
      oldrsq <- rsq
      ttt <- (-1)*M_adjoint*M*cg_p
      pkp <- (-1)*cg_p.M_adjoint*M.cg_p
      a <- -rsq/pkp
      dest <- dest + a*cg_p
      resid <- resid + a*ttt
      rsq <- |resid|^2
      b <- rsq/oldrsq
      cg_p <- resid + b*cg_p
    */

    while( (rsq>rsqstop) && (iteration%niter!=0) ) {
      b = rsq / oldrsq;
      /* cg_p  <- resid + b*cg_p */
      QDP_V_eq_r_times_V_plus_V(cg_p, &b, cg_p, resid, q_parity);

      oldrsq = rsq;

      DSLASH_QDP_FN_SPECIAL2(cg_p, tttt, q_otherparity, temp1);
      DSLASH_QDP_FN_SPECIAL2(tttt, ttt, q_parity, temp2);
      iteration++;

      /* finish computation of M_adjoint*m*p and p*M_adjoint*m*Kp */
      /* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
      /* pkp  <- cg_p.(ttt - msq*cg_p) */
      QDP_V_meq_r_times_V(ttt, &msq_x4, cg_p, q_parity);
      QDP_r_eq_re_V_dot_V(&pkp, cg_p, ttt, q_parity);

      a = - rsq / pkp;

      /* dest <- dest - a*cg_p */
      /* resid <- resid - a*ttt */
      QDP_V_peq_r_times_V(dest, &a, cg_p, q_parity);
      QDP_V_peq_r_times_V(resid, &a, ttt, q_parity);
      QDP_r_eq_norm2_V(&rsq, resid, q_parity);
    }

  } while( (rsq>rsqstop) && (iteration<nrestart*niter) );

  if( rsq <= rsqstop ) {
    dtimec += dclock();
#ifdef CGTIME
    if(QDP_this_node==0) {
      printf("CONGRAD5: time = %e (fn_qdp %c) masses = 1 iters = %d mflops = %e\n",
	     dtimec, QDP_Precision, iteration,
	     (double)(nflop*volume*iteration/(1.0e6*dtimec*numnodes())) );
#ifdef REMAP
      node0_printf("CGREMAP:  time = %e\n",remaptime);
#endif
      fflush(stdout);
    }
#endif
  } else {
    if(QDP_this_node==0) {
      printf("ks_congrad_qdp: CG not converged after %d iterations, res. = %e wanted %e\n",
	     iteration, rsq, rsqstop);
      fflush(stdout);
    }
  }
  qic->final_rsq = rsq;
  total_iters += iteration;
  return(iteration);
}

/* For field-based src and dest */

int
KS_CONGRAD_MILCFIELD2QDP(su3_vector *f_src, su3_vector *f_dest, 
			 quark_invert_control *qic, Real mass,
			 ferm_links_t *fn )
{
  QLA_Real qmass;
  QDP_ColorVector *src, *dest;
  double remaptime;
  int iteration;
  int parity = qic->parity;

  remaptime = -dclock();
  src = QDP_create_V();
  dest = QDP_create_V();

  set_V_from_field(src, f_src,EVENANDODD);
  set_V_from_field(dest, f_dest,EVENANDODD);

  qmass = (QLA_Real) mass;
  remaptime += dclock();
  iteration = KS_CONGRAD_QDP(src, dest, qic, qmass, fn );
  remaptime -= dclock();

  set_field_from_V(f_dest, dest, parity);

  QDP_destroy_V(dest); dest = NULL;
  QDP_destroy_V(src);  src = NULL;
  remaptime += dclock();

#ifdef CGTIME
#ifdef REMAP
  node0_printf("CGREMAP:  time = %e\n",remaptime);
#endif
#endif

  return(iteration);
}

/* For site-based src and dest */

int
KS_CONGRAD_MILC2QDP(field_offset f_src, field_offset f_dest, 
		    quark_invert_control *qic, Real mass,
		    ferm_links_t *fn )
{
  QLA_Real qmass;
  QDP_ColorVector *src, *dest;
  double remaptime;
  int iteration;
  int parity = qic->parity;

  remaptime = -dclock();
  src = QDP_create_V();
  dest = QDP_create_V();

  set_V_from_site(src, f_src,EVENANDODD);
  set_V_from_site(dest, f_dest,EVENANDODD);

  qmass = (QLA_Real) mass;
  remaptime += dclock();
  iteration = KS_CONGRAD_QDP(src, dest, qic, qmass, fn );
  remaptime -= dclock();

  set_site_from_V(f_dest, dest, parity);

  QDP_destroy_V(dest); dest = NULL;
  QDP_destroy_V(src);  src = NULL;
  remaptime += dclock();

#ifdef CGTIME
#ifdef REMAP
  node0_printf("CGREMAP:  time = %e\n",remaptime);
#endif
#endif

  return(iteration);
}

