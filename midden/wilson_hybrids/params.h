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
}  params;


#endif /* _PARAMS_H */
