/******* d_congrad5_fn_two_src2.c  *************************/

/* Two-mass CG inverter for staggered fermions */
/*This program inverts the equations: src1=A*x1
 *                                   src2= (A+shift)*x2
 *by offsetting the soln. vectors: x1 -->x0 + x1(prime)
                                   x2 -->x0 + x2(prime).
By choosing x0 = (src2-src1)/shift, the sources become equal:
                                       common_source = src1 - A*x0 = A*x1(prime)
				       common_source = src1 - A*x0 = ( A + shift )*x2(prime).
These two eqns are then solved using the multimass method
and then the final solns. are:  soln. #1 = x0 + x1(prime)
                                soln. #2 = x0 + x2(prime)               */

/* Based on B. Jegerlehner, hep-lat/9612014.
   See also A. Frommer, S. G\"usken, T. Lippert, B. N\"ockel,"
   K. Schilling, Int. J. Mod. Phys. C6 (1995) 627. */

/* This version is based on d_congrad5_fn.c and d_congrad5_eo.c */

/* For "fat link actions", ie when FN is defined, this version
   assumes connection to nearest neighbor points is stored in fatlink.
   For actions with a Naik term, it assumes the connection to third
   nearest neighbors is in longlink. */

/* Steve Bildstein, 2004 original version
** DT March 2005 clean up, make use of temp vectors the only option
** Also now requires -DFN
*/

#ifndef FN
Until we have EO versions of dslash_field, we require FN
#endif

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/dslash_ks_redefine.h"   /* Actually not used, yet */
#include "../include/loopend.h"

/*#define CGTIME*/

int ks_congrad_two_src(	/* Return value is number of iterations taken */
    field_offset src1,    /* source vector (type su3_vector) */
    field_offset src2,
    field_offset dest1,	/* solution vectors */
    field_offset dest2,
    Real mass1,
    Real mass2,
    int niter,		/* maximal number of CG interations */
    int nrestart,       /* maximal number of CG restarts */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* internal precision for the inversion (ignored) */
    int parity,		/* parity to be worked on */
    Real  *final_rsq_ptr, /* final residue squared */
    imp_ferm_links_t *fn       /* Storage for fermion links */
    )
{
    /* Site su3_vector's resid, cg_p and ttt are used as temporaries */
    register int i;
    register site *s;
    int iteration;       /* counter for iterations */
    double c1, c2, rsq, oldrsq, pkp;            /* pkp = cg_p.K.cg_p */
    double source_norm,source_norm1;	/* squared magnitude of source vector */
    double rsqstop;	/* stopping residual normalized by source norm */
    int l_parity;	/* parity we are currently doing */
    int l_otherparity;	/* the other parity */
    msg_tag *tags1[16], *tags2[16];	/* tags for gathers to parity and opposite */
    int special_started;	/* 1 if dslash_special has been called */
    int j, jud, jstrange;
    double shift, msq_xm4;
    double zeta_i[2], zeta_im1[2], zeta_ip1[2];
    double beta_i[2], beta_im1[2], alpha[2];
    int first;                                                 
    su3_vector *temp;
    su3_vector *destvec1;
    su3_vector *destvec2;
    su3_vector *pm_strange;
    su3_vector *init_guess;
    su3_vector *common_source;
    su3_vector *ttt;
    su3_vector *cg_p;
    su3_vector *resid;

/* Timing */
#ifdef CGTIME
    double dtimed,dtimec;
#endif
    double nflop;

/* debug */
#ifdef CGTIME
    dtimec = -dclock();   
#endif
    first = 0;
    nflop = 1187;	/* THIS LOOKS WRONG - DT */
    if(parity==EVENANDODD)nflop *=2;
	
    special_started = 0;
    /* if we want both parities, we will do even first. */
    switch(parity){
	case(EVEN): l_parity=EVEN; l_otherparity=ODD; break;
	case(ODD):  l_parity=ODD; l_otherparity=EVEN; break;
	case(EVENANDODD):  l_parity=EVEN; l_otherparity=ODD; break;
     }
    jud = 0;
    jstrange = 1;
    temp=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
    destvec1=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
    destvec2=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
    ttt=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
    cg_p=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
    resid=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
    pm_strange = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
    init_guess = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector) );
    common_source = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector) );
    shift = 4.0*( mass2*mass2 - mass1*mass1);
    
    msq_xm4 = -4.0*mass1*mass1;                              
    iteration = 0;

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
	/*This loop calculates init_guess = (phi2 - phi1)/shift = X0 , which
	 * is used to equalize the two sources*/
	
	FORSOMEPARITY(i,s,l_parity){
	    sub_su3_vector( (su3_vector *)F_PT(s,src2),(su3_vector *)F_PT(s,src1),&temp[i] );      
	    scalar_mult_su3_vector(&temp[i],1.0/shift,&init_guess[i] );
       }END_LOOP                                                    
      			                                               
	/*set temp = -D(adj)*D*init_guess. "D" means "D-slash" */              
	    dslash_fn_field_special(init_guess,temp,l_otherparity,tags2,1,fn);
	    dslash_fn_field_special(temp,temp,l_parity,tags1,1,fn);
	    cleanup_gathers(tags1,tags2);
                                                                       
	    source_norm=0.0;    
	    source_norm1=0.0;
/**This part of the code sets common_source = src1  -( D(adj)*D +4*mass1^2) * init_guess.
 * ** or, phi1 and phi2 ----> common_source = phi1 - A*X0  
 * **It also sets resid = cg_p = pm_strange = common_source (an initialization)
 * ** It also initializes dest1=dest2=X0, ie. I am adding (back) the initial
 *    guess, X0, at the beginning.    */

	FORSOMEPARITY(i,s,l_parity) {
		scalar_mult_add_su3_vector( &temp[i],&init_guess[i],msq_xm4,&temp[i] );
		add_su3_vector((su3_vector *)F_PT(s,src1),&temp[i],&common_source[i] );
	        source_norm += (double)magsq_su3vec( &common_source[i] );	
		source_norm1 += (double)magsq_su3vec( (su3_vector *)F_PT(s,src1) );  
		/*pm_strange[i] = cg_p[i] = resid[i] = common_source[i];*/
		su3vec_copy( &common_source[i],&(resid[i]) );
		su3vec_copy(&(resid[i]),&(cg_p[i]) );
                su3vec_copy(&(resid[i]), &pm_strange[i]);
                su3vec_copy(&init_guess[i],&destvec1[i] );
                su3vec_copy(&init_guess[i],&destvec2[i] );
	} END_LOOP

	rsq = source_norm;
        g_doublesum(&rsq);       
        iteration++ ;  /* iteration counts number of multiplications
                           by M_adjoint*M */
	total_iters++;
	rsqstop = (double)rsqmin*source_norm1 ;
	//node0_printf("MCG_STOP_RSQ =  %le ",rsqstop);fflush(stdout);     
        g_doublesum(&rsqstop);                                                 
        /*note- I am using source_norm1.
        * source_norm1 = |phi1|^2
	* source_norm  = |phi1-A*X0|^2
	* Since X = X0 + X(prime)
	* and X0 is constant throughout,
	* it makes sense to calculate
	* X(prime) to the accuracy with
	* which we calculate X in the
	* regular congrad5 inverter*/
	/**node0_printf("congrad: source_norm = %e\n", (double)source_norm);**/    
	for(j=0;j<2;j++){                                                       
	    zeta_im1[j] = zeta_i[j] = 1.0;
	    beta_im1[j] = -1.0;
	    alpha[j] = 0.0;
	}
          
	/* We are now all initialized, as per Krylov paper reference. */
                                              
    do{
	oldrsq = rsq;                                      
	/* sum of neighbors */
/* We now proceed to calculate pkp which is the name for -cg_p*( D(adj)*D + 4*mass1^2)* cg_p   */
	if(special_started==0){
	    dslash_fn_field_special( cg_p, ttt, l_otherparity, tags2, 1, fn );
	    dslash_fn_field_special( ttt, ttt, l_parity, tags1, 1, fn );
	    special_started = 1;
	}
	else {
	    dslash_fn_field_special( cg_p, ttt, l_otherparity, tags2, 0, fn );
	    dslash_fn_field_special( ttt, ttt, l_parity, tags1, 0, fn );
	}

	/* finish computation of (-1)*M_adjoint*m*p and (-1)*p*M_adjoint*M*p */
	/* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
	/* pkp  <- cg_p . ttt */
	pkp = 0.0;
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(ttt[i]), &(cg_p[i]), msq_xm4,
		&(ttt[i]) );
	    pkp += (double)su3_rdot( &(cg_p[i]), &(ttt[i]) );
	} END_LOOP
	g_doublesum( &pkp );
	iteration++;
	total_iters++;

       /* This part of the code calculates beta, zeta, and beta(sigma) and zeta(sigma) */
	
	beta_i[jud] = -rsq / pkp;
	zeta_ip1[jud] = 1.0;
	    zeta_ip1[jstrange] = zeta_i[jstrange] * zeta_im1[jstrange] * beta_im1[jud];
	    c1 = beta_i[jud] * alpha[jud] * (zeta_im1[jstrange]-zeta_i[jstrange]);
	    c2 = zeta_im1[jstrange] * beta_im1[jud] * (1.0 + shift*beta_i[jud]);
	    zeta_ip1[jstrange] /= c1 + c2;
	    beta_i[jstrange] = beta_i[jud] * zeta_ip1[jstrange] / zeta_i[jstrange];
	/* dest <- dest + beta*cg_p */                     
	FORSOMEPARITY(i,s,l_parity){                   
	    scalar_mult_add_su3_vector( &(destvec1[i]),
		&(cg_p[i]), (Real)beta_i[jud], &(destvec1[i]) );
	    scalar_mult_add_su3_vector( &(destvec2[i]),
	        &pm_strange[i], (Real)beta_i[jstrange], &(destvec2[i]) );
	} END_LOOP

	/* resid <- resid + beta*ttt */
	rsq = 0.0;
	
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(resid[i]), &(ttt[i]),       /*ttt = -( D(adj)D+4*mass1^2)*cg_p */
		(Real)beta_i[jud], &(resid[i]));
	    rsq += (double)magsq_su3vec( &(resid[i]) );
	                                                           
	} END_LOOP                                                      
	/**** calculates a stopping resid^2 comparable to that of single-mass congrad*/
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
	    FORALLSITES(i,s){
		su3vec_copy( &(destvec1[i]), (su3_vector *)F_PT(s,dest1) );
		su3vec_copy( &(destvec2[i]), (su3_vector *)F_PT(s,dest2) );
	    }
	    
   if(special_started==1) {
		cleanup_gathers(tags1,tags2);
		special_started = 0;
	    }                                     

#ifdef CGTIME
	    dtimec += dclock();                                               
	    if(this_node==0){printf("CONGRAD5: time = %e (Bildstein 2 source) masses = 1 iters = %d mflops = %e\n",
dtimec,iteration,(double)(nflop*volume*iteration/(1.0e6*dtimec*numnodes())) );
		fflush(stdout);}                                       
#endif  
 
	    free(temp);
	    free(destvec1);
	    free(destvec2);
	    free(pm_strange);
	    free(init_guess);
	    free(common_source);
	    free(ttt);
	    free(cg_p);
	    free(resid);

             return (iteration);  
        }                                                                                            
	alpha[jud] = rsq/oldrsq ;
	alpha[jstrange] = alpha[jud] * zeta_ip1[jstrange] * beta_i[jstrange] /
		       (zeta_i[jstrange] * beta_i[jud]);
	/* cg_p  <- resid + alpha*cg_p */
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(resid[i]), &(cg_p[i]),
		(Real)alpha[jud], &(cg_p[i]));
	    scalar_mult_su3_vector( &(resid[i]),
		(Real)zeta_ip1[jstrange], &(ttt[i]));
	    scalar_mult_add_su3_vector( &(ttt[i]), &pm_strange[i],
		(Real)alpha[jstrange], &pm_strange[i]);
	} END_LOOP

	/* scroll the scalars */
	for(j=0;j<2;j++){
	    beta_im1[j] = beta_i[j];
	    zeta_im1[j] = zeta_i[j];
	    zeta_i[j] = zeta_ip1[j];
	}

    } while( iteration < niter );                   
                                                       
    node0_printf(
	"ks_congrad_two_src: CG not converged after %d iterations, res. = %e wanted %e\n",
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
    FORALLSITES(i,s){
	su3vec_copy( &(destvec1[i]), (su3_vector *)F_PT(s,dest1) );
	su3vec_copy( &(destvec2[i]), (su3_vector *)F_PT(s,dest2) );
    }
    
    if(special_started==1){	/* clean up gathers */
	cleanup_gathers(tags1, tags2);
	special_started = 0;
    }
    free(temp);
    free(destvec1);
    free(destvec2);
    free(pm_strange);
    free(init_guess);
    free(common_source);
    free(ttt);
    free(cg_p);
    free(resid);
return (iteration);  
}                                                        
