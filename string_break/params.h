#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "lattice.h"  /* For MAX_LEVEL */

/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;	/* 1 if it is time to stop */
   /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */
	int iseed;	/* for random numbers */
   /*  REPEATING BLOCK */
	int no_smear_level;	/* number of smearing levels (<=MAX_LEVEL) */
	int smear_num[MAX_LEVEL];	/* the number of smearing iterations */
	int off_axis_flag;	/* off-axis Wilson loops or not? */
	int startflag;	/* what to do for beginning lattice */
	int saveflag;	/* what to do with lattice at end */
	Real smear_fac;	/* smearing factor = weight of direct link */
	Real mass;	/* quark mass */
	int r0;		/* size of "fuzzy" source/sink */
	int num_src;	/* number of Gaussian random sources */
	int niter;	/* maximum number of c.g. iterations */
	int nrestart;	/* maximum number of c.g. restarts */
	Real rsqmin;	/* for deciding on convergence */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
	char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
}  params;

#endif /* _PARAMS_H */
