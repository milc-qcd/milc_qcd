#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/generic_quark_types.h"
/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;   /* 1 if it is time to stop */
    /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;	/* lattice dimensions */
#ifdef RANDOM
        int iseed;      /* for random numbers */
#endif
    /*  REPEATING BLOCK */
	int startflag;	/* what to do for beginning lattice */
	int fixflag;    /* whether to gauge fix */
	int saveflag;	/* what to do for saving lattice */
	int startflag_w[MAX_MASSES];	/* what to do for beginning wilson vector */
	int saveflag_w[MAX_MASSES];	/* what to do for saving wilson vector */
	int num_masses;	/* number of masses */
	Real mass[MAX_MASSES];	/* masses values for multiple propagators */
	Real resid[MAX_MASSES];	/* residue for invertion convergence */
	quark_source wqs[MAX_MASSES];  /* source parameters */
	int niter; 	/* maximum number of c.g. iterations */
	int nrestart; 	/* maximum number of c.g. restarts */
	char startfile[MAXFILENAME];
	char startfile_w[MAX_MASSES][MAXFILENAME];
	char savefile[MAXFILENAME];
	char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
	char savefile_w[MAX_MASSES][MAXFILENAME];
	char scratchstem_w[MAXFILENAME];
        int scratchflag;
}  params;

#endif /* _PARAMS_H */
