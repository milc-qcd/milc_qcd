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
	int saveflag;   /* what to do with lattice at end */
	Real beta,mass; /* gauge coupling, quark mass */
	Real staple_weight; /* relative weight of staple in fat link */
	int niter; 	/* maximum number of c.g. iterations */
	Real rsqprop;  /* for deciding on convergence */
	int source_start, source_inc, n_sources; /* source time and increment */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
}  params;
#endif /* _PARAMS_H */
