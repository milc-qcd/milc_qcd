/******* ks_multicg_rhmc.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* Multi-mass CG inverter for staggered fermions */
/* OBSOLETE!! UNSUPPORTED!! CD 11/2006 */

/* Inverts M^dagger*M+offset for an array of offsets
   For "ordinary" use (not RHMC) offset = 4*mass^2
   Note difference between "offsets", the arguments, and "shifts", which
   are differences between offset and minimimum offset, and which get
   reordered during the calculation.
*/

/* Based on B. Jegerlehner, hep-lat/9612014.
   See also A. Frommer, S. G\"usken, T. Lippert, B. N\"ockel,"
   K. Schilling, Int. J. Mod. Phys. C6 (1995) 627. 

   This version is based on d_congrad5_fn.c and d_congrad5_eo.c 

   For "fat link actions", ie when FN is defined, this version
   assumes connection to nearest neighbor points is stored in fatlink.
   For actions with a Naik term, it assumes the connection to third
   nearest neighbors is in longlink.

*/


#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "../include/dslash_ks_redefine.h"

#include "../include/loopend.h"

static su3_vector *ttt,*cg_p;
static su3_vector *resid;
static int first_multicongrad = 1;


//#ifdef NOTDEFINED
int ks_multicg_reverse(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* desired intermediate precision */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    fn_links_t *fn      /* Storage for fat and Naik links */
    )
{
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
    msg_tag *tags1[16], *tags2[16];	/* tags for gathers to parity and opposite */
    int special_started;	/* 1 if dslash_special has been called */
    int j, j_low;
    Real *shifts, mass_low, msq_xm4;
    double *zeta_i, *zeta_im1, *zeta_ip1;
    double *beta_i, *beta_im1, *alpha;
    // su3_vector **pm;	/* vectors not involved in gathers */

    // Switch indices
    su3_vector **psim_rev; su3_vector *psim_space;
    su3_vector **pm_rev; su3_vector *pm_space;

/* Timing */
#ifdef CGTIME
    double dtimed,dtimec;
#endif
    double nflop;
    if( num_offsets==0 )return(0);

    // Switch indices
    psim_rev = (su3_vector **)malloc( sizeof(su3_vector *)*sites_on_node );
    psim_space = (su3_vector *)malloc( sizeof(su3_vector)*sites_on_node*num_offsets );
    pm_rev = (su3_vector **)malloc( sizeof(su3_vector *)*sites_on_node );
    pm_space = (su3_vector *)malloc( sizeof(su3_vector)*sites_on_node*num_offsets );
    if( psim_space == NULL || pm_space == NULL){printf("NO ROOM!\n"); exit(0); }
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
	shifts[j] = offsets[j];
	if (offsets[j] < mass_low){
	    mass_low = offsets[j];
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
	    source_norm += (double) magsq_su3vec( (su3_vector *)F_PT(s,src) );
	    su3vec_copy((su3_vector *)F_PT(s,src), &(resid[i]));
	    su3vec_copy(&(resid[i]), &(cg_p[i]));
	    clearvec(&(psim_rev[i][j_low]));
	    for(j=0;j<num_offsets;j++) if(j!=j_low){
		clearvec(&(psim_rev[i][j]));
		su3vec_copy(&(resid[i]), &(pm_rev[i][j]));
	    }
	} END_LOOP
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
				      fn);
	    dslash_fn_field_special( ttt, ttt, l_parity, tags1, 1,
				     fn);
	    special_started = 1;
	}
	else {
	    dslash_fn_field_special( cg_p, ttt, l_otherparity, tags2, 0, fn);
	    dslash_fn_field_special( ttt, ttt, l_parity, tags1, 0, fn );
	}
#else
	dslash_site( F_OFFSET(cg_p), F_OFFSET(ttt), l_otherparity, fn );
	dslash_site( F_OFFSET(ttt), F_OFFSET(ttt), l_parity, fn);
#endif

	/* finish computation of (-1)*M_adjoint*m*p and (-1)*p*M_adjoint*M*p */
	/* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
	/* pkp  <- cg_p . ttt */
	pkp = 0.0;
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(ttt[i]), &(cg_p[i]), msq_xm4, &(ttt[i]) );
	    pkp += (double)su3_rdot( &(cg_p[i]), &(ttt[i]) );
	} END_LOOP
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
	} END_LOOP

	/* resid <- resid + beta*ttt */
	rsq = 0.0;
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(resid[i]), &(ttt[i]),
		(Real)beta_i[j_low], &(resid[i]));
	    rsq += (double)magsq_su3vec( &(resid[i]) );
	} END_LOOP
	g_doublesum(&rsq);

	if( rsq <= rsqstop ){
	    /* if parity==EVENANDODD, set up to do odd sites and go back */
	    if(parity == EVENANDODD) {
		l_parity=ODD; l_otherparity=EVEN;
		parity=EVEN;	/* so we won't loop endlessly */
		iteration = 0;
		goto start;
	    }
	    *final_rsq_ptr = (Real)rsq;

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
	      printf("CONGRAD5: time = %e iters = %d mflops = %e\n",
		     dtimec,iteration,
		     (double)(nflop)*volume*
		     iteration/(1.0e6*dtimec*numnodes()));
		fflush(stdout);}
#endif
             return (iteration);
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
	} END_LOOP

	/* scroll the scalars */
	for(j=0;j<num_offsets_now;j++){
	    beta_im1[j] = beta_i[j];
	    zeta_im1[j] = zeta_i[j];
	    zeta_i[j] = zeta_ip1[j];
	}

    } while( iteration < niter );

    node0_printf(
	"CG not converged after %d iterations, res. = %e wanted %e\n",
	iteration, rsq, rsqstop);
    fflush(stdout);

    /* if parity==EVENANDODD, set up to do odd sites and go back */
    if(parity == EVENANDODD) {
	l_parity=ODD; l_otherparity=EVEN;
	parity=EVEN;	/* so we won't loop endlessly */
	iteration = 0;
	goto start;
    }

    *final_rsq_ptr=rsq;

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

    return(iteration);
}
//#endif //NOTDEFINED

int ks_multicg_revhyb(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* desired intermediate precision */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    fn_links_t *fn      /* Storage for fat and Naik links */
    )
{
    int i,j,iters=0; site *s;
    ks_multicg_reverse( src, psim, offsets, num_offsets, niter, 
			rsqmin, prec, parity, final_rsq_ptr, fn);
    for(i=0;i<num_offsets;i++){
      FORSOMEPARITY(j,s,parity){ s->xxx1 = psim[i][j]; } END_LOOP
       iters += ks_congrad( src, F_OFFSET(xxx1), 0.5*sqrt(offsets[i]), 
			    niter/5, rsqmin, prec, 
			    parity, final_rsq_ptr, fn );
      FORSOMEPARITY(j,s,parity){ psim[i][j] = s->xxx1; } END_LOOP
    }
    return(iters);
}
