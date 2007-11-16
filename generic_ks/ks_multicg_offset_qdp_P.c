/******* ks_multicg_qdp_P.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* Multi-mass CG inverter for staggered fermions */

/* NOTE: This code is actually an include file for ks_multicg_qdp_F.c
   and ks_multicg_qdp_D.c, so any edits should be consistent with this
   purpose. */

/* Based on B. Jegerlehner, hep-lat/9612014.
   See also A. Frommer, S. G\"usken, T. Lippert, B. N\"ockel,"
   K. Schilling, Int. J. Mod. Phys. C6 (1995) 627. */

/* This version is based on d_congrad5_fn.c and d_congrad5_eo.c */

/* For "fat link actions", ie when FN is defined, this version
   assumes connection to nearest neighbor points is stored in fatlinks.
   For actions with a Naik term, it assumes the connection to third
   nearest neighbors is in longlinks. */

#if ( QDP_Precision == 'F' )

#define KS_MULTICG_MASS   ks_multicg_mass_F
#define KS_MULTICG_OFFSET ks_multicg_offset_F
#define DSLASH_QDP_FN_SPECIAL2 dslash_qdp_F_fn_special2
#define SETUP_DSLASH      setup_dslash_F
#define MYREAL            float
#define BCKLINK           bcklink_F
#define FATLINKS          fatlinks_F
#define LONGLINKS         longlinks_F
#define IMPLINKS          implinks_F

#else

#define KS_MULTICG_MASS   ks_multicg_mass_D
#define KS_MULTICG_OFFSET ks_multicg_offset_D
#define DSLASH_QDP_FN_SPECIAL2 dslash_qdp_D_fn_special2
#define SETUP_DSLASH      setup_dslash_D
#define MYREAL            double
#define BCKLINK           bcklink_D
#define FATLINKS          fatlinks_D
#define LONGLINKS         longlinks_D
#define IMPLINKS          implinks_D

#endif

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qdp.h"
#include "../include/generic_ks_qdp.h"
#include <lattice_qdp.h>
#ifndef FN
#error only FN supported
#endif

static int congrad_setup=0;
static QDP_ColorVector *ttt, *tttt, *resid, *cg_p;
static QDP_ColorVector *temp1[16], *temp2[16];
extern QDP_ColorMatrix *BCKLINK[8];
extern QDP_ColorMatrix **FATLINKS, **LONGLINKS, *IMPLINKS[8];

extern void SETUP_DSLASH(void);

static void
setup_congrad(void)
{
  int i;

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

static int
ks_multicg_mass_qdp(	/* Return value is number of iterations taken */
	       QDP_ColorVector *src,	/* source vector (type su3_vector) */
	       QDP_ColorVector **dest,	/* solution vectors */
	       QLA_Real *masses,	/* the masses */
	       int num_masses,	        /* number of masses */
	       int niter,		/* maximal number of CG interations */
	       QLA_Real rsqmin,	        /* desired residue squared */
	       QDP_Subset parity,	/* parity to be worked on */
	       QLA_Real *final_rsq_ptr, /* final residue squared */
	       ferm_links_t *fn           /* Storage for fat and Naik links */
	       )
{
  /* Site su3_vector's resid, cg_p and ttt are used as temporaies */
  double c1, c2, oldrsq;		/* pkp = cg_p.K.cg_p */
  double rsqstop;	/* stopping residual normalized by source norm */
  QLA_Real rsq, pkp;
  QLA_Real source_norm;	/* squared magnitude of source vector */
  QLA_Real msq_xm4;
  double *beta_im1;
  double *zeta_i, *zeta_im1;
  QLA_Real *shifts;
  QLA_Real *beta_i, *alpha, *zeta_ip1;
  QDP_ColorVector **pm;	   /* vectors not involved in gathers */
  QDP_Subset q_parity, q_otherparity;
  int iteration;	   /* counter for iterations */
  int j, j_low;

/* Timing */
#ifdef CGTIME
  double dtimec;
  double nflop;
  dtimec = -dclock(); 
  nflop = 1187;
  if(parity==QDP_all) nflop *= 2;
#endif

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

  shifts = (QLA_Real *)malloc(num_masses*sizeof(QLA_Real));
  zeta_i = (double *)malloc(num_masses*sizeof(double));
  zeta_im1 = (double *)malloc(num_masses*sizeof(double));
  zeta_ip1 = (QLA_Real *)malloc(num_masses*sizeof(QLA_Real));
  beta_i = (QLA_Real *)malloc(num_masses*sizeof(QLA_Real));
  beta_im1 = (double *)malloc(num_masses*sizeof(double));
  alpha = (QLA_Real *)malloc(num_masses*sizeof(QLA_Real));

  j_low = 0;
  for(j=0; j<num_masses; j++) {
    shifts[j] = 4.0*masses[j]*masses[j];
    if (masses[j] < masses[j_low]) {
      j_low = j;
    }
  }

  pm = (QDP_ColorVector **) malloc(num_masses*sizeof(QDP_ColorVector *));
  for(j=0; j<num_masses; j++) if(j!=j_low) {
    pm[j] = QDP_create_V();
    shifts[j] -= shifts[j_low];
  }
  msq_xm4 = -shifts[j_low];

  iteration = 0;

  if(!congrad_setup) setup_congrad();

  set4_M_from_field(FATLINKS, fn->fat, EVENANDODD);
  set4_M_from_field(LONGLINKS, fn->lng, EVENANDODD);
  {
    QDP_ColorMatrix *tcm;
    int i;
    tcm = QDP_create_M();
    for(i=0; i<8; ++i) {
      QDP_M_eq_sM(tcm, IMPLINKS[i], shiftdirs[i], QDP_backward, QDP_all);
      QDP_M_eqm_Ma(BCKLINK[i], tcm, QDP_all);
    }
    QDP_destroy_M(tcm);
    tcm = NULL;
  }

#ifdef CGTIME
  dtimec = -dclock(); 
#endif

  /* initialization process */
  /* start: */

  QDP_r_eq_norm2_V(&source_norm, src, q_parity);
  QDP_V_eq_V(resid, src, q_parity);
  QDP_V_eq_V(cg_p, src, q_parity);
  QDP_V_eq_zero(dest[j_low], q_parity);
  for(j=0; j<num_masses; j++) if(j!=j_low) {
    QDP_V_eq_zero(dest[j], q_parity);
    QDP_V_eq_V(pm[j], resid, q_parity);
  }

  rsq = source_norm;

  iteration++;   /* iteration counts number of multiplications
		    by M_adjoint*M */
  total_iters++;
  rsqstop = rsqmin * source_norm;
  /**node0_printf("congrad: source_norm = %e\n", (double)source_norm);**/

  for(j=0;j<num_masses;j++){
    zeta_im1[j] = zeta_i[j] = 1.0;
    beta_im1[j] = -1.0;
    alpha[j] = 0.0;
  }

  do {
    oldrsq = rsq;
    /* sum of neighbors */

    DSLASH_QDP_FN_SPECIAL2(cg_p, tttt, q_otherparity, temp1);
    DSLASH_QDP_FN_SPECIAL2(tttt, ttt, q_parity, temp2);

    /* finish computation of (-1)*M_adjoint*m*p and (-1)*p*M_adjoint*M*p */
    /* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
    /* pkp  <- cg_p . ttt */
    QDP_V_peq_r_times_V(ttt, &msq_xm4, cg_p, q_parity);
    QDP_r_eq_re_V_dot_V(&pkp, cg_p, ttt, q_parity);

    iteration++;
    total_iters++;

    beta_i[j_low] = - rsq / pkp;

    zeta_ip1[j_low] = 1.0;
    for(j=0;j<num_masses;j++) if(j!=j_low){
      zeta_ip1[j] = zeta_i[j] * zeta_im1[j] * beta_im1[j_low];
      c1 = beta_i[j_low] * alpha[j_low] * (zeta_im1[j]-zeta_i[j]);
      c2 = zeta_im1[j] * beta_im1[j_low] * (1.0+shifts[j]*beta_i[j_low]);
      /*THISBLOWSUP*/
	/** zeta_ip1[j] /= c1 + c2;
	beta_i[j] = beta_i[j_low] * zeta_ip1[j] / zeta_i[j];**/

      /*TRYTHIS*/
      if( c1+c2 != 0.0 )
	zeta_ip1[j] /= c1 + c2; 
      else {
	zeta_ip1[j] = 0.0;
      }
      if( zeta_i[j] != 0.0){
	beta_i[j] = beta_i[j_low] * zeta_ip1[j] / zeta_i[j];
      } else  {
	zeta_ip1[j] = 0.0;
	beta_i[j] = 0.0;
      }	

    }

    /* dest <- dest + beta*cg_p */
    QDP_V_peq_r_times_V(dest[j_low], &beta_i[j_low], cg_p, q_parity);
    for(j=0; j<num_masses; j++) if(j!=j_low) {
      QDP_V_peq_r_times_V(dest[j], &beta_i[j], pm[j], q_parity);
    }

    /* resid <- resid + beta*ttt */
    QDP_V_peq_r_times_V(resid, &beta_i[j_low], ttt, q_parity);
    QDP_r_eq_norm2_V(&rsq, resid, q_parity);

    if( rsq <= rsqstop ) {
      *final_rsq_ptr = (Real)rsq;

      /* Free stuff */
      for(j=0;j<num_masses;j++) if(j!=j_low) QDP_destroy_V(pm[j]);
      free(pm);
      pm = NULL;

      free(zeta_i);   zeta_i = NULL;
      free(zeta_ip1); zeta_ip1 = NULL;
      free(zeta_im1); zeta_im1 = NULL;
      free(beta_i);   beta_i = NULL; 
      free(beta_im1); beta_im1 = NULL;
      free(alpha);    alpha = NULL;
      free(shifts);   shifts = NULL;

#ifdef CGTIME
      dtimec += dclock();
      if(this_node==0) {
	printf("CONGRAD5: time = %e (multicg_qdp %c) iters = %d masses = %d mflops = %e\n",
	       dtimec,QDP_Precision,iteration,num_masses,
	       (double)(nflop*volume*iteration/(1.0e6*dtimec*numnodes())) );
	fflush(stdout);
      }
#endif
      return (iteration);
    }

    alpha[j_low] = rsq / oldrsq;

    for(j=0;j<num_masses;j++) if(j!=j_low){
      /*THISBLOWSUP
      alpha[j] = alpha[j_low] * zeta_ip1[j] * beta_i[j] /
	(zeta_i[j] * beta_i[j_low]);
      */
      /*TRYTHIS*/
      if( zeta_i[j] * beta_i[j_low] != 0.0)
	alpha[j] = alpha[j_low] * zeta_ip1[j] * beta_i[j] /
	  (zeta_i[j] * beta_i[j_low]);
      else {
	alpha[j] = 0.0;
      }
    }

    /* cg_p  <- resid + alpha*cg_p */
    QDP_V_eq_r_times_V_plus_V(cg_p, &alpha[j_low], cg_p, resid, q_parity);
    for(j=0; j<num_masses; j++) if(j!=j_low) {
      QDP_V_eq_r_times_V(ttt, &zeta_ip1[j], resid, q_parity);
      QDP_V_eq_r_times_V_plus_V(pm[j], &alpha[j], pm[j], ttt, q_parity);
    }

    /* scroll the scalars */
    for(j=0;j<num_masses;j++){
      beta_im1[j] = beta_i[j];
      zeta_im1[j] = zeta_i[j];
      zeta_i[j] = zeta_ip1[j];
    }

  } while( iteration < niter );

  node0_printf("ks_multicg_mass_qdp: CG not converged after %d iterations, res. = %e wanted %e\n",
	       iteration, rsq, rsqstop);
  fflush(stdout);

  *final_rsq_ptr = rsq;

  /* Free stuff */
  for(j=0;j<num_masses;j++) if(j!=j_low) QDP_destroy_V(pm[j]);
  free(pm);           pm = NULL;

  free(zeta_i);       zeta_i = NULL;
  free(zeta_ip1);     zeta_ip1 = NULL;
  free(zeta_im1);     zeta_im1 = NULL;
  free(beta_i);       beta_i = NULL;
  free(beta_im1);     beta_im1 = NULL;
  free(alpha);        alpha = NULL;
  free(shifts);       shifts = NULL;

  return(iteration);
}

int
KS_MULTICG_MASS(	/* Return value is number of iterations taken */
	   field_offset f_src,	/* source vector (type su3_vector) */
	   su3_vector **psim,	/* solution vectors */
	   Real *masses,	/* the masses */
	   int num_masses,	/* number of masses */
	   int niter,		/* maximal number of CG interations */
	   Real rsqmin,	/* desired residue squared */
	   int parity,		/* parity to be worked on */
	   Real *final_rsq_ptr,	/* final residue squared */
	   ferm_links_t *fn       /* Storage for fat and Naik links */
	   )
{
  QLA_Real qrsqmin, qfinal_rsq_ptr, *qmasses;
  QDP_ColorVector *src, **dest;
  QDP_Subset q_parity;
  int iteration, i;

  switch(parity) {
  case(EVEN):  q_parity = QDP_even; break;
  case(ODD):   q_parity = QDP_odd;  break;
  default:     q_parity = QDP_all;  break;
  }

  src = QDP_create_V();
  set_V_from_site(src, f_src, EVENANDODD);

  dest = (QDP_ColorVector **) malloc(num_masses*sizeof(QDP_ColorVector *));
  qmasses = (QLA_Real *) malloc(num_masses*sizeof(QLA_Real));
  for(i=0; i<num_masses; i++) {
    dest[i] = QDP_create_V();
    qmasses[i] = (QLA_Real) masses[i];
  }
  qrsqmin = (QLA_Real) rsqmin;

  iteration = ks_multicg_mass_qdp(src, dest, qmasses, num_masses, niter, 
				  qrsqmin, q_parity, &qfinal_rsq_ptr, fn);
  *final_rsq_ptr = (Real) qfinal_rsq_ptr;

  for(i=0; i<num_masses; i++) {
    set_field_from_V(psim[i], dest[i], parity);
    QDP_destroy_V(dest[i]);
  }
  free(qmasses);      qmasses = NULL;
  free(dest);         dest = NULL;
  QDP_destroy_V(src); src = NULL;

  return(iteration);
}

/* Just a wrapper for ks_multicg_mass */
/* Offsets are 4 * mass * mass and must be positive */
int KS_MULTICG_OFFSET(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    ferm_links_t *fn      /* Storage for fat and Naik links */
    )
{
  int i;
  Real *masses;
  int status;
  int num_masses = num_offsets;
  char myname[] = "ks_multicg_offset";

  if( num_offsets==0 )return(0);

  masses = (Real *)malloc(sizeof(Real)*num_offsets);
  if(masses == NULL){
    printf("%s(%d): No room for masses\n",myname,this_node);
    terminate(1);
  }
  for(i = 0; i < num_offsets; i++){
    if(offsets[i] < 0){
      printf("%s(%d): called with negative offset %e\n",
	     myname,this_node,offsets[i]);
      terminate(1);
    }
    masses[i] = sqrt(offsets[i]/4.0);
  }
  status = KS_MULTICG_MASS(src, psim, masses, num_masses, niter, rsqmin,
			   parity, final_rsq_ptr, fn);
  free(masses);
  return status;
}

