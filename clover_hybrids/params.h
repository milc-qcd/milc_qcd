#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */

/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;   /* 1 if it is time to stop */
    /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */
	int iseed;	/* for random numbers */
    /*  REPEATING BLOCK */
	int startflag;  /* what to do for beginning lattice */
	Real beta,kappa; /* gauge coupling, quark hopping parameter */
	int source_start, source_inc, n_sources; /* source time and increment */
	int niter; 	/* maximum number of c.g. iterations */
	Real rsqprop;  /* for deciding on convergence */
	char startfile[MAXFILENAME];
	Real clov_c,u0;
        int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
	int boundary_flag  ;
        int wot_src ;   /*** source for the quark propagator inverter ***/
        int verbose_flag ; /*** flag controlling the amount of debug information to print ***/


	int oper_PION_SOURCE  ;
	int oper_PION2_SOURCE  ;
	int oper_RHO_SOURCE   ; 
	int oper_RHO2_SOURCE ;
	int oper_A1P_SOURCE   ;
	int oper_A1_SOURCE   ;
	int oper_ZEROMP_SOURCE ;
	int oper_ZEROPM_SOURCE ;
	int oper_ZEROPMP_SOURCE ;
	int oper_ZEROPMB_SOURCE ;
	int oper_ZEROMM_SOURCE ; 
	int oper_ZEROMMP_SOURCE ;
	int oper_ONEMP_SOURCE  ;
	int oper_ONEMP2_SOURCE  ;
	int oper_ONEMM_SOURCE  ;
	int oper_ONEPP_SOURCE  ;
	int oper_QQQQ_SOURCE   ;

	/** parameters controlling the smearing of f_mu_nu ***/
	Real  space_simple_weight  ;
	Real  space_norm_factor   ;
	Real  time_simple_weight ;
	Real  time_norm_factor   ;
	int  smearing_level ;

}  params;

#endif /* _PARAMS_H */
