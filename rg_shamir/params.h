#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "defines.h"

/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;   /* 1 if it is time to stop */
    /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */
	int iseed;	/* for random numbers */
	int nflavors;	/* the number of flavors */
    /*  REPEATING BLOCK */
	Real mass;      /*  quark mass */
	Real u0; /* tadpole parameter */
        int nrg; /* RG blocking steps */
	int niter; 	/* maximum number of c.g. iterations */
	Real rsqprop;  /* for deciding on convergence */
	int startflag;  /* what to do for beginning lattice */
	int saveflag;   /* what to do with lattice at end */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
	char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
        char propfile[MAXFILENAME];  /* For blocked quark propagator */
}  params;

#endif /* _PARAMS_H */
