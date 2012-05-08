/******* ks_multicg.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* Wrappers for multi-mass CG inverter for staggered fermions

   8/12 C. DeTar added macros for selecting the multicg inverter option
*/


#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/loopend.h"
#include <string.h>

/* Set the KS multicg inverter flavor depending on the macro KS_MULTICG */
/* Defaults to hybrid */

#ifndef KS_MULTICG
#define KS_MULTICG HYBRID
#endif

/* Forward declarations */
static int ks_multicg_fake_field(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* the offsets */
    int num_offsets,	/* number of offsets */
    quark_invert_control qic[], 
    imp_ferm_links_t *fn[]    /* Storage for fermion links */
);

static int ks_multicg_hybrid_field(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* the offsets */
    int num_offsets,	/* number of offsets */
    quark_invert_control qic[],
    imp_ferm_links_t *fn[]    /* Storage for fermion links */
);

// static void ks_multicg_reverse_field(	/* Return value is number of iterations taken */
//     su3_vector *src,	/* source vector (type su3_vector) */
//     su3_vector **psim,	/* solution vectors */
//     ks_param *ksp,	/* the offsets */
//     int num_offsets,	/* number of offsets */
//     quark_invert_control *qic,
//     imp_ferm_links_t *fn      /* Storage for fermion links */
// );
// 
// static void ks_multicg_revhyb_field(	/* Return value is number of iterations taken */
//     su3_vector *src,	/* source vector (type su3_vector) */
//     su3_vector **psim,	/* solution vectors */
//     ks_param *ksp,	/* the offsets */
//     int num_offsets,	/* number of offsets */
//     quark_invert_control *qic,
//     imp_ferm_links_t *fn[]     /* Storage for fermion links */
// );

enum ks_multicg_opt_t {OFFSET, HYBRID, FAKE, REVERSE, REVHYB};


/**********************************************************************/
/*   Wrapper for the multimass inverter with multiple offsets         */
/**********************************************************************/
int ks_multicg_field(   /* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* the offsets */
    int num_offsets,	/* number of offsets */
    quark_invert_control qic[], /* inversion parameters */
    imp_ferm_links_t *fn[]    /* Storage for fermion links */
    )
{

  int iters;

  switch(KS_MULTICG){
    // CD option removed.  It doesn't work for HISQ.
//  case OFFSET:
//    ks_multicg_offset_field( src, psim, ksp, num_offsets, qic, fn);
//    break;
  case HYBRID:
    iters = ks_multicg_hybrid_field( src, psim, ksp, num_offsets, qic, fn);
    break;
  case FAKE:
    iters = ks_multicg_fake_field( src, psim, ksp, num_offsets, qic, fn);
    break;
    // CD option removed.  It doesn't work for HISQ.
//  case REVERSE:
//    ks_multicg_reverse_field( src, psim, ksp, num_offsets, qic, fn);
//    break;
//  case REVHYB:
//    ks_multicg_revhyb_field( src, psim, ksp, num_offsets, qic, fn);
//    break;
  default:
    iters = ks_multicg_hybrid_field( src, psim, ksp, num_offsets, qic, fn);
  }
//  /* Report status information stored in the first entry */
//  report_status(qic+0);
  return iters;
}

/**********************************************************************/
/*   Accessor for string describing the option                        */
/**********************************************************************/
const char *ks_multicg_opt_chr( void )
{
  switch(KS_MULTICG){
  case OFFSET:
    return "OFFSET";
    break;
  case HYBRID:
    return "HYBRID";
    break;
  case FAKE:
    return "FAKE";
    break;
  case REVERSE:
    return "REVERSE";
    break;
  case REVHYB:
    return "REVHYB";
    break;
  default:
    return "HYBRID";
  }
  return NULL;
}

// mock up multicg by repeated calls to ordinary cg
static int ks_multicg_fake_field(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* the offsets */
    int num_offsets,	/* number of offsets */
    quark_invert_control qic[], 
    imp_ferm_links_t *fn_multi[]    /* Storage for fermion links */
    )
{int i,iters=0;

  for(i=0;i<num_offsets;i++){
    ks_congrad_field( src, psim[i], qic+i, 0.5*sqrt(ksp[i].offset), 
		      fn_multi[i] );
    //    report_status(qic+i);
    iters += qic->final_iters;
  }
  return iters;
}

#if 0
static Real *create_offsets_from_ksp(ks_param *ksp, int num_offsets){
  Real *offsets;
  int i;

  offsets = (Real *)malloc(num_offsets*sizeof(Real));
  if(offsets == NULL){
    printf("create_offsets_from_ksp(%d): No room\n",this_node);
    terminate(1);
  }
  
  for(i = 0; i < num_offsets; i++){
    offsets[i] = ksp[i].offset;
  }
  return offsets;
}
#endif

// Do a multimass CG followed by calls to individual CG's
// to finish off.
static int ks_multicg_hybrid_field(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* the offsets */
    int num_offsets,	/* number of offsets */
    quark_invert_control qic[],
    imp_ferm_links_t *fn_multi[]    /* Storage for fermion links */
    )
{
  int i,multi_iters=0,iters=0;

#if defined(HALF_MIXED) && !defined(USE_CG_GPU)
  /* Do multicg in single precision.  (The GPU routine does this automatically for HALF_MIXED) */
  int prec_save = qic[0].prec;
  qic[0].prec = 1;
#endif
  
  /* First we invert as though all masses took the same Naik epsilon */
  multi_iters = iters = 
    ks_multicg_offset_field( src, psim, ksp, num_offsets, qic, fn_multi[0]);
  report_status(qic+0);

#if defined(HALF_MIXED) && !defined(USE_CG_GPU)
  qic[0].prec = prec_save;
#endif

  /* Then we refine using the correct Naik epsilon */
  for(i=0;i<num_offsets;i++){
#ifdef NO_REFINE
    if(fn_multi[i] == fn_multi[0])
      continue;
#endif
    ks_congrad_field_cpu( src, psim[i], qic+i, 0.5*sqrt(ksp[i].offset), fn_multi[i] );
    iters += qic[i].final_iters;
    qic[i].final_iters += multi_iters;
  }

  return iters;
}

#if 0

#include "../include/dslash_ks_redefine.h"

#include "../include/loopend.h"

static su3_vector *ttt,*cg_p;
static su3_vector *resid;
static int first_multicongrad = 1;

/* This choice is currently not supported in QDP or QOP.  Perhaps it
   can be rewritten so it just rearranges the masses and calls the
   standard multimass inverter.  Since it is not supported in QDP or
   QOP, the precision argument "prec" has not effect. */

static void ks_multicg_reverse_field(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* KS parametes, including the offsets */
    int num_offsets,	/* number of offsets */
    quark_invert_control *qic,
    imp_ferm_links_t *fn      /* Storage for fermion links */
    )
{
    char myname[] = "ks_multicg_reverse_field";
    /* Site su3_vector's resid, cg_p and ttt are used as temporaies */
    register int i;
    register site *s;
    int iteration;	/* counter for iterations */
    int num_offsets_now; /* number of offsets still being worked on */
    double c1, c2, rsq, oldrsq, pkp;		/* pkp = cg_p.K.cg_p */
    double source_norm;	/* squared magnitude of source vector */
    double rsqstop;	/* stopping residual normalized by source norm */
    int l_parity=0;	/* parity we are currently doing */
    int l_otherparity=0; /* the other parity */
#ifdef FN
    msg_tag *tags1[16], *tags2[16];	/* tags for gathers to parity and opposite */
#endif
    int special_started;	/* 1 if dslash_special has been called */
    int j, j_low;
    Real *shifts, mass_low, msq_xm4;
    double *zeta_i, *zeta_im1, *zeta_ip1;
    double *beta_i, *beta_im1, *alpha;
    // su3_vector **pm;	/* vectors not involved in gathers */

    // Switch indices
    su3_vector **psim_rev; su3_vector *psim_space;
    su3_vector **pm_rev; su3_vector *pm_space;

    /* Unpack structure */
    /* We don't restart this algorithm, so we adopt the convention of
       taking the product here */
    int niter        = qic->max*qic->nrestart;
    Real rsqmin      = qic->resid * qic->resid;    /* desired squared residual - 
					 normalized as sqrt(r*r)/sqrt(src_e*src_e) */
    int parity       = qic->parity;   /* EVEN, ODD */
    
/* Timing */
#ifdef CGTIME
    double dtimec;
#endif
    double nflop;

    qic->final_iters = 0;
    qic->final_restart = 0;

    //#if FERM_ACTION == HISQ
    //    fn->hl.current_X_set = 0;
    //    restore_fn_links(fn);
    //#endif
    if( num_offsets==0 )return;

    if(fn == NULL){
      printf("%s(%d): Called with NULL fn\n", myname, this_node);
      terminate(1);
    }

    // Switch indices
    psim_rev = (su3_vector **)malloc( sizeof(su3_vector *)*sites_on_node );
    psim_space = (su3_vector *)malloc( sizeof(su3_vector)*sites_on_node*num_offsets );
    pm_rev = (su3_vector **)malloc( sizeof(su3_vector *)*sites_on_node );
    pm_space = (su3_vector *)malloc( sizeof(su3_vector)*sites_on_node*num_offsets );
    if( psim_space == NULL || pm_space == NULL){printf("%s: NO ROOM!\n",myname); exit(0); }
    for( i=0; i<sites_on_node; i++ ){
	psim_rev[i] = &(psim_space[num_offsets*i]);
	pm_rev[i] = &(pm_space[num_offsets*i]);
	for( j=0; j<num_offsets; j++){
	    psim_rev[i][j] = psim[j][i];
	}
    }

/* debug */
#ifdef CGTIME
    dtimec = -dclock(); 
#endif

    nflop = 1205 + 15*num_offsets;
    if(parity==EVENANDODD)nflop *=2;
	
    special_started = 0;
    /* if we want both parities, we will do even first. */
    switch(parity){
	case(EVEN): l_parity=EVEN; l_otherparity=ODD; break;
	case(ODD):  l_parity=ODD; l_otherparity=EVEN; break;
	case(EVENANDODD):  l_parity=EVEN; l_otherparity=ODD; break;
    }

    shifts = (Real *)malloc(num_offsets*sizeof(Real));
    zeta_i = (double *)malloc(num_offsets*sizeof(double));
    zeta_im1 = (double *)malloc(num_offsets*sizeof(double));
    zeta_ip1 = (double *)malloc(num_offsets*sizeof(double));
    beta_i = (double *)malloc(num_offsets*sizeof(double));
    beta_im1 = (double *)malloc(num_offsets*sizeof(double));
    alpha = (double *)malloc(num_offsets*sizeof(double));

    //pm = (su3_vector **)malloc(num_offsets*sizeof(su3_vector *));
    mass_low = 1.0e+20;
    j_low = -1;
    for(j=0;j<num_offsets;j++){
	shifts[j] = ksp[j].offset;
	if (ksp[j].offset < mass_low){
	    mass_low = ksp[j].offset;
	    j_low = j;
	}
    }
    for(j=0;j<num_offsets;j++) if(j!=j_low){
	//pm[j] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
	shifts[j] -= shifts[j_low];
    }
    msq_xm4 = -shifts[j_low];


    iteration = 0;

#define PAD 0
    /* now we can allocate temporary variables and copy then */
    /* PAD may be used to avoid cache thrashing */
    if(first_multicongrad) {
      ttt = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
      cg_p = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
      resid = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
      first_multicongrad = 0;
    }

#ifdef CGTIME
    dtimec = -dclock(); 
#endif



    /* initialization process */
    start:
#ifdef FN
	if(special_started==1) {        /* clean up gathers */
	    cleanup_gathers(tags1, tags2);
	    special_started = 0;
	}
#endif
	num_offsets_now = num_offsets;
	source_norm = 0.0;
	FORSOMEPARITY(i,s,l_parity){
	    source_norm += (double) magsq_su3vec( src+i );
	    su3vec_copy( src+i, &(resid[i]));
	    su3vec_copy(&(resid[i]), &(cg_p[i]));
	    clearvec(&(psim_rev[i][j_low]));
	    for(j=0;j<num_offsets;j++) if(j!=j_low){
		clearvec(&(psim_rev[i][j]));
		su3vec_copy(&(resid[i]), &(pm_rev[i][j]));
	    }
	} END_LOOP;
	g_doublesum( &source_norm );
	rsq = source_norm;

        iteration++ ;  /* iteration counts number of multiplications
                           by M_adjoint*M */
	total_iters++;
	rsqstop = rsqmin * source_norm;
	/**node0_printf("congrad: source_norm = %e\n", (double)source_norm);**/

	for(j=0;j<num_offsets;j++){
	    zeta_im1[j] = zeta_i[j] = 1.0;
	    beta_im1[j] = -1.0;
	    alpha[j] = 0.0;
	}

    do{
	oldrsq = rsq;
	/* sum of neighbors */

#ifdef FN
	if(special_started==0){
	    dslash_fn_field_special( cg_p, ttt, l_otherparity, tags2, 1, 
				     fn );
	    dslash_fn_field_special( ttt, ttt, l_parity, tags1, 1, fn );
	    special_started = 1;
	}
	else {
	    dslash_fn_field_special( cg_p, ttt, l_otherparity, tags2, 0, 
				     fn);
	    dslash_fn_field_special( ttt, ttt, l_parity, tags1, 0, fn );
	}
#else
	dslash_site( F_OFFSET(cg_p), F_OFFSET(ttt), l_otherparity, fn);
	dslash_site( F_OFFSET(ttt), F_OFFSET(ttt), l_parity, fn);
#endif

	/* finish computation of (-1)*M_adjoint*m*p and (-1)*p*M_adjoint*M*p */
	/* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
	/* pkp  <- cg_p . ttt */
	pkp = 0.0;
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(ttt[i]), &(cg_p[i]), msq_xm4, &(ttt[i]) );
	    pkp += (double)su3_rdot( &(cg_p[i]), &(ttt[i]) );
	} END_LOOP;
	g_doublesum( &pkp );
	iteration++;
	total_iters++;

	beta_i[j_low] = -rsq / pkp;

	zeta_ip1[j_low] = 1.0;
	for(j=0;j<num_offsets_now;j++) if(j!=j_low){
	    zeta_ip1[j] = zeta_i[j] * zeta_im1[j] * beta_im1[j_low];
	    c1 = beta_i[j_low] * alpha[j_low] * (zeta_im1[j]-zeta_i[j]);
	    c2 = zeta_im1[j] * beta_im1[j_low] * (1.0+shifts[j]*beta_i[j_low]);
	    /*THISBLOWSUP
	    zeta_ip1[j] /= c1 + c2;
	    beta_i[j] = beta_i[j_low] * zeta_ip1[j] / zeta_i[j];
	    */
	    /*TRYTHIS*/
	    if( c1+c2 != 0.0 )zeta_ip1[j] /= c1 + c2; else zeta_ip1[j] = 0.0;
	    if( zeta_i[j] != 0.0){
		beta_i[j] = beta_i[j_low] * zeta_ip1[j] / zeta_i[j];
	    } else  {
		zeta_ip1[j] = 0.0;
//node0_printf("SETTING A ZERO, j=%d, num_offsets_now=%d\n",j,num_offsets_now);
		//if(j==num_offsets_now-1)node0_printf("REDUCING OFFSETS\n");
		if(j==num_offsets_now-1)num_offsets_now--;
		// don't work any more on finished solutions
		// this only works if largest offsets are last, otherwise
		// just wastes time multiplying by zero
	    }
	}

	/* dest <- dest + beta*cg_p */
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(psim_rev[i][j_low]),
		&(cg_p[i]), (Real)beta_i[j_low], &(psim_rev[i][j_low]));
	    for(j=0;j<num_offsets_now;j++) if(j!=j_low){
		scalar_mult_add_su3_vector( &(psim_rev[i][j]),
		    &(pm_rev[i][j]), (Real)beta_i[j], &(psim_rev[i][j]));
	    }
	} END_LOOP;

	/* resid <- resid + beta*ttt */
	rsq = 0.0;
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(resid[i]), &(ttt[i]),
		(Real)beta_i[j_low], &(resid[i]));
	    rsq += (double)magsq_su3vec( &(resid[i]) );
	} END_LOOP;
	g_doublesum(&rsq);

	if( rsq <= rsqstop ){
	    /* if parity==EVENANDODD, set up to do odd sites and go back */
	    if(parity == EVENANDODD) {
		l_parity=ODD; l_otherparity=EVEN;
		parity=EVEN;	/* so we won't loop endlessly */
		iteration = 0;
		goto start;
	    }

#ifdef FN
	    if(special_started==1) {
		cleanup_gathers(tags1,tags2);
		special_started = 0;
	    }
#endif

	    /* Free stuff */
	    //for(j=0;j<num_offsets;j++) if(j!=j_low) free(pm[j]);
	    //free(pm);

	    free(zeta_i);
	    free(zeta_ip1);
	    free(zeta_im1);
	    free(beta_i);
	    free(beta_im1);
	    free(alpha);
	    free(shifts);

	    for( i=0; i<sites_on_node; i++ ) for( j=0; j<num_offsets; j++){
		 psim[j][i] = psim_rev[i][j];
	    }
	    free(psim_space); free(psim_rev);
	    free(pm_space); free(pm_rev);

#ifdef CGTIME
	    dtimec += dclock();
	    if(this_node==0){
	      printf("CONGRAD5: time = %e (multicg_rev) iters = %d masses = %d mflops = %e\n",
		     dtimec,iteration,num_offsets,
		     (double)(nflop)*volume*
		     iteration/(1.0e6*dtimec*numnodes()));
		fflush(stdout);}
#endif
	    qic->size_r        = (Real)rsq/source_norm;
	    qic->size_relr     = 0.;
	    qic->final_iters   = iteration;
	    qic->final_restart = 0;
	    qic->converged     = 1;
	    return;
        }

	alpha[j_low] = rsq / oldrsq;

	for(j=0;j<num_offsets_now;j++) if(j!=j_low){
	    /*THISBLOWSUP
	    alpha[j] = alpha[j_low] * zeta_ip1[j] * beta_i[j] /
		       (zeta_i[j] * beta_i[j_low]);
	    */
	    /*TRYTHIS*/
	    if( zeta_i[j] * beta_i[j_low] != 0.0)alpha[j] = alpha[j_low] * zeta_ip1[j] * beta_i[j] /
		       (zeta_i[j] * beta_i[j_low]);
	    else alpha[j] = 0.0;
	}

	/* cg_p  <- resid + alpha*cg_p */
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(resid[i]), &(cg_p[i]),
		(Real)alpha[j_low], &(cg_p[i]));
	    for(j=0;j<num_offsets_now;j++) if(j!=j_low){
		scalar_mult_su3_vector( &(resid[i]),
		    (Real)zeta_ip1[j], &(ttt[i]));
		scalar_mult_add_su3_vector( &(ttt[i]), &(pm_rev[i][j]),
		    (Real)alpha[j], &(pm_rev[i][j]));
	    }
	} END_LOOP;

	/* scroll the scalars */
	for(j=0;j<num_offsets_now;j++){
	    beta_im1[j] = beta_i[j];
	    zeta_im1[j] = zeta_i[j];
	    zeta_i[j] = zeta_ip1[j];
	}

    } while( iteration < niter );

    qic->converged = 0;
    node0_printf(
	"%s CG not converged after %d iterations, res. = %e wanted %e\n",
	myname, iteration, rsq, rsqstop);
    fflush(stdout);

    /* if parity==EVENANDODD, set up to do odd sites and go back */
    if(parity == EVENANDODD) {
	l_parity=ODD; l_otherparity=EVEN;
	parity=EVEN;	/* so we won't loop endlessly */
	iteration = 0;
	goto start;
    }

#ifdef FN
    if(special_started==1){	/* clean up gathers */
	cleanup_gathers(tags1, tags2);
	special_started = 0;
    }
#endif

    /* Free stuff */
    //for(j=0;j<num_offsets;j++) if(j!=j_low) free(pm[j]);
    //free(pm);

    free(zeta_i);
    free(zeta_ip1);
    free(zeta_im1);
    free(beta_i);
    free(beta_im1);
    free(alpha);
    free(shifts);

    for( i=0; i<sites_on_node; i++ ) for( j=0; j<num_offsets; j++){
	 psim[j][i] = psim_rev[i][j];
    }
    free(psim_space); free(psim_rev);
    free(pm_space); free(pm_rev);

    qic->size_r        = (Real)rsq/source_norm;
    qic->size_relr     = 0.;
    qic->final_iters   = iteration;
    qic->final_restart = 0;
}

/* This choice is currently not supported in QDP or QOP so the
   precision argument "prec" can affect only the cleanup
   inversions, if at all. */

static void ks_multicg_revhyb_field(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* the offsets */
    int num_offsets,	/* number of offsets */
    quark_invert_control *qic,
    imp_ferm_links_t *fn_multi[]     /* Storage for fermion links */
    )
{
  int i,iters=0,index;
  
  /* First we invert without the Naik epsilon correction*/
  ks_multicg_reverse_field( src, psim, ksp, num_offsets, qic, fn_multi[0]);
  /* Then we polish with the correct Naik epsilon */
  for(i=0;i<num_offsets;i++){
//#if FERM_ACTION == HISQ
//    //    fn->hl.current_X_set = ksp[i].naik_term_epsilon_index;
//    //    restore_fn_links(fn);
//    index = ksp[i].naik_term_epsilon_index;
//#else
//    index = 0;
//#endif
    iters += ks_congrad_field_cpu( src, psim[i], qic, 0.5*sqrt(ksp[i].offset), fn_multi[i] );
    report_status(qic);
  }
}

#endif //if 0

/**********************************************************************/
/*  Wrapper for the multimass inverter with multiple masses          */
/**********************************************************************/

int ks_multicg_mass_field(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* KS parameters, including masses */
    int num_masses,	/* number of masses */
    quark_invert_control qic[],  /* inversion parameters */
    imp_ferm_links_t *fn_multi[]     /* Storage for fat and Naik links */
				)
{
  int i, iters;

  for(i = 0; i < num_masses; i++){
    ksp[i].offset = 4.0*ksp[i].mass*ksp[i].mass;
  }
  
  iters = ks_multicg_field(src, psim, ksp, num_masses, qic, fn_multi);
  
  return iters;
}


/**********************************************************************/
/*  Old-style interface for the multimass inverter with multiple masses          */
/**********************************************************************/

int ks_multicg_mass_site(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *masses,	/* the masses */
    int num_masses,	/* number of masses */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* internal precision for inversion (ignored) */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr, /* final residue squared */
    imp_ferm_links_t *fn_multi[]     /* Storage for fat and Naik links */
    )
{
  su3_vector *in;
  quark_invert_control *qic;
  ks_param *ksp;
  int i, iters;

  /* Set up inversion control structure */
  /* For molecular dynamics they are identical */
  qic = (quark_invert_control *)malloc(num_masses*sizeof(quark_invert_control));
  if(qic == NULL){
    printf("ks_multicg_mass_site: No room for qic\n");
    terminate(1);
  }

  for(i = 0; i < num_masses; i++){
    qic[i].prec = prec;
    qic[i].min = 0;
    qic[i].max = niter;
    qic[i].nrestart = nrestart;
    qic[i].parity = parity;
    qic[i].nsrc = 1;
    qic[i].start_flag = 0;
    qic[i].resid = sqrt(rsqmin);
    qic[i].relresid = 0;
  
  }

  /* Map MILC field to temporary */
  in = create_v_field_from_site_member(src);
  
  /* Set up ksp structure */
  ksp = (ks_param *)malloc(num_masses*sizeof(ks_param));
  if(ksp == NULL){
    printf("ks_multicg_mass_site(%d): No room\n",this_node);
    terminate(1);
  }
  for(i = 0; i < num_masses; i++){
    ksp[i].mass = masses[i];
    ksp[i].offset = 0;
  }
  
  iters = ks_multicg_mass_field(in, psim, ksp, num_masses, qic, fn_multi);
  
  destroy_v_field(in);
  free(ksp);
  
  *final_rsq_ptr = qic[0].final_rsq;
  free(qic);
  return iters;
}


/* Multimass version of M^{-1} src.  */

/* We might also try the UML trick here, and polish with the
   single-mass inverter. */

int mat_invert_multi(
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **dst,	/* solution vectors */
    ks_param *ksp,	/* KS parameters, including masses */
    int num_masses,	/* number of masses */
    quark_invert_control qic[],  /* inversion parameters */
    imp_ferm_links_t *fn_multi[]   /* Storage for fat and Naik links */
    )
{
  int i, tot_iters = 0;
  
  /* Use the single-mass inverter if there is only one mass */
  if(num_masses == 1)
    return mat_invert_uml_field(src, dst[0], qic+0, ksp->mass, fn_multi[0] );
  
  /* Convert masses to offsets for ks_multicg_mass_field */
  for(i = 0; i < num_masses; i++){
    ksp[i].offset = 4.0*ksp[i].mass*ksp[i].mass;
  }

  /* Multimass inversion on separate "parities" to get dst = (M M_adj)^-1 src */

  
  for(i = 0; i < num_masses; i++)
    qic[i].parity = EVEN;
  tot_iters += ks_multicg_mass_field(src, dst, ksp, num_masses, qic, fn_multi);

  for(i = 0; i < num_masses; i++)
    qic[i].parity = ODD;
  tot_iters += ks_multicg_mass_field(src, dst, ksp, num_masses, qic, fn_multi);

  /* Multiply all solutions by Madjoint to get dst = M^-1 * src */
  for(i = 0; i < num_masses; i++){
    ks_dirac_adj_op_inplace( dst[i], ksp[i].mass, EVENANDODD, fn_multi[i]);
  }
  return tot_iters;
}
