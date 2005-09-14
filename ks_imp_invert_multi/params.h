#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */

/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;   /* 1 if it is time to stop */
    /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */
	int iseed;	/* for random numbers */
	int nflavors1;  /* the number of flavors of first type */
	int nflavors2;  /* the number of flavors of second type */

    /*  REPEATING BLOCK */
	Real beta,mass1,mass2; /* gauge coupling, quark masses */
	Real u0; /* tadpole parameter */
	int niter; 	/* maximum number of c.g. iterations */
	Real rsqmin,rsqprop;  /* for deciding on convergence */
        int source_start, source_inc, n_sources; /* source time and increment */
	int fpi_nmasses;
	Real fpi_mass[MAX_FPI_NMASSES]; 
        params_mminv mminv;
	int startflag;  /* what to do for beginning lattice */
	int fixflag;    /* whether to gauge fix */
	int saveflag;   /* what to do with lattice at end */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
	char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
        int kssaveflag;	/* whether to save ks propagator */
        char kssavefile[MAXFILENAME];
}  params;

#endif /* _PARAMS_H */
