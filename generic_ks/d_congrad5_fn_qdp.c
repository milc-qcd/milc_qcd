/******* d_congrad5_fn_qdp.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
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
#include "generic_ks_includes.h"	/* definitions files and prototypes */

static int congrad_setup=0;
static QDP_ColorVector *ttt, *tttt, *resid, *cg_p;
static QDP_ColorVector *temp1[16], *temp2[16];
extern QDP_ColorMatrix *bcklink[8];

void setup_dslash(void);

static void
setup_congrad(void)
{
  int i;

  setup_dslash();
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
ks_congrad_qdp(QDP_ColorVector *src, QDP_ColorVector *dest, QLA_Real mass,
	       int niter, QLA_Real rsqmin, QDP_Subset parity,
	       QLA_Real *final_rsq_ptr)
{
  QLA_Real a,b;		 /* Sugar's a,b,resid**2,last resid*2 */
  QLA_Real rsq,oldrsq,pkp; /* pkp = cg_p.K.cg_p */
  QLA_Real msq_x4;	 /* 4*mass*mass */
  QLA_Real source_norm;	 /* squared magnitude of source vector */
  QLA_Real rsqstop;	 /* stopping residual normalized by source norm */
  QDP_Subset q_parity, q_otherparity;
  int iteration;	 /* counter for iterations */

/* Timing */
#ifdef CGTIME
  double dtimec;
  double nflop;
  nflop = 1187;
  if(parity==QDP_all) nflop *= 2;
#endif

  if(!congrad_setup) setup_congrad();

  if(parity==QDP_odd) {
    q_parity = QDP_odd;
    q_otherparity = QDP_even;
  } else if(parity==QDP_even) {
    q_parity = QDP_even;
    q_otherparity = QDP_odd;
  } else {
    q_parity = QDP_all;
    q_otherparity = QDP_all;
  }
  msq_x4 = 4.0*mass*mass;
  iteration = 0;

  if (!valid_fatlinks) load_fatlinks();
  if (!valid_longlinks) load_longlinks();
  set4_M_from_temp(fatlinks, t_fatlink);
  set4_M_from_temp(longlinks, t_longlink);

  {
    QDP_ColorMatrix *tcm;
    int i;
    tcm = QDP_create_M();
    for(i=0; i<8; ++i) {
      QDP_M_eq_sM(tcm, implinks[i], shiftdirs[i], QDP_backward, QDP_all);
      QDP_M_eqm_Ma(bcklink[i], tcm, QDP_all);
    }
    QDP_destroy_M(tcm);
  }

#ifdef CGTIME
  dtimec = -dclock(); 
#endif

  /* initialization process */
  do {
    /* ttt <-  (-1)*M_adjoint*M*dest
       resid,cg_p <- src + ttt
       rsq = |resid|^2
       source_norm = |src|^2
    */

    QDP_V_eq_V(cg_p, dest, q_parity);
    dslash_qdp_fn_special2(cg_p, tttt, q_otherparity, temp1);
    dslash_qdp_fn_special2(tttt, ttt, q_parity, temp2);
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

      dslash_qdp_fn_special2(cg_p, tttt, q_otherparity, temp1);
      dslash_qdp_fn_special2(tttt, ttt, q_parity, temp2);
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

  } while( (rsq>rsqstop) && (iteration<5*niter) );

  if( rsq <= rsqstop ) {
#ifdef CGTIME
    dtimec += dclock();
    if(QDP_this_node==0) {
      printf("CONGRAD5: time = %e iters = %d mflops = %e\n",
	     dtimec, iteration,
	     (double)(nflop*volume*iteration/(1.0e6*dtimec*numnodes())) );
      fflush(stdout);
    }
#endif
  } else {
    if(QDP_this_node==0) {
      printf("CG not converged after %d iterations, res. = %e wanted %e\n",
	     iteration, rsq, rsqstop);
      fflush(stdout);
    }
  }
  *final_rsq_ptr = rsq;
  total_iters += iteration;
  return(iteration);
}

int
ks_congrad(field_offset f_src, field_offset f_dest, Real mass,
	   int niter, Real rsqmin, int parity, Real *final_rsq_ptr)
{
  QLA_Real qmass, qrsqmin, qfinal_rsq_ptr;
  QDP_ColorVector *src, *dest;
  QDP_Subset q_parity;
  int iteration;

  switch(parity) {
  case(EVEN):  q_parity = QDP_even; break;
  case(ODD):   q_parity = QDP_odd;  break;
  default:     q_parity = QDP_all;  break;
  }

  src = QDP_create_V();
  dest = QDP_create_V();

  set_V_from_field(src, f_src);
  set_V_from_field(dest, f_dest);

  qmass = (QLA_Real) mass;
  qrsqmin = (QLA_Real) rsqmin;
  iteration = ks_congrad_qdp(src, dest, qmass, niter, qrsqmin, q_parity,
			     &qfinal_rsq_ptr);
  *final_rsq_ptr = (Real) qfinal_rsq_ptr;

  set_field_from_V(f_dest, dest);

  QDP_destroy_V(dest);
  QDP_destroy_V(src);

  return(iteration);
}
