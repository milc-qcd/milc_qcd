#ifndef _PARAMS_H
#define _PARAMS_H

#include "lattice.h"    /* For MAX_KAP */
#include "../include/macros.h"	/* For MAXFILENAME */

/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;	/* 1 if it is time to stop */
    /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;	/* lattice dimensions */
    /*  REPEATING BLOCK */
	int startflag;	/* what to do for beginning lattice */
	int num_kap;	/* number of kappa's */
	Real kappa;	/* quark hopping parameter */
	Real clov_c;	/* clover coefficient */
	Real kap[MAX_KAP];	/* kappa values for multiple propagators */
	Real resid[MAX_KAP];	/* residue for invertion convergence */
	int niter; 	/* maximum number of c.g. iterations */
	int nrestart; 	/* maximum number of c.g. restarts */
	int bc_flag; 	/* gauge boundary condition flag */
	int num_smear;
	Real alpha;	/* APE smearing parameter (Boulder convention) */
	Real ferm_phas[3];	/* fermion phase factors */
	char startfile[MAXFILENAME];
}  params;

#endif	/* _PARAMS_H */
