#ifndef _PARAMS_H
#define _PARAMS_H

#include "lattice.h"  /* For MAX_KAP */
#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/generic_quark_types.h" /* For wilson_quark_source */

/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;   /* 1 if it is time to stop */
    /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;	/* lattice dimensions */
    /*  REPEATING BLOCK */
	int startflag;	/* what to do for beginning lattice */
	int fixflag;    /* whether to gauge fix */
	int saveflag;	/* what to do for saving lattice */
	int startflag_w[MAX_KAP];	/* what to do for beginning wilson vector */
	int saveflag_w[MAX_KAP];	/* what to do for saving wilson vector */
	int num_kap;	/* number of kappa's */
	Real clov_c,u0;	/* clover coefficient, <Tr(U_p)>^{1/4} */
	Real kap[MAX_KAP];	/* kappa values for multiple propagators */
	Real resid[MAX_KAP];	/* residue for invertion convergence */
	wilson_quark_source wqs[MAX_KAP];  /* source parameters */
	int niter; 	/* maximum number of c.g. iterations */
	int nrestart; 	/* maximum number of c.g. restarts */
	char startfile[MAXFILENAME];
	char startfile_w[MAX_KAP][MAXFILENAME];
	char savefile[MAXFILENAME];
	char savefile_w[MAX_KAP][MAXFILENAME];
	char scratchstem_w[MAXFILENAME];
        int scratchflag;
}  params;


#endif /* _PARAMS_H */
