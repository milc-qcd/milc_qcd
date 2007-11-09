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
   /*  REPEATING BLOCK */
	int startflag;  /* what to do for beginning lattice */
	int fixflag;    /* whether to gauge fix */
	int fixflag_ft; /* whether to FFT gauge fix */
	int saveflag;   /* what to do with lattice at end */
	Real beta;	/* gauge coupling */
#ifdef QUARK_PROP
	int num_mass;	/* number of masses */
	Real mass[MAX_NUM_MASS];	/* quark mass */
	Real u0;	/* tadpole parameter */
	int niter;	/* maximum number of c.g. iterations */
	int nrestart;	/* maximum number of c.g. restart */
	Real rsqprop;	/* for deciding on convergence */
	int run_CG_flag[MAX_NUM_MASS];	/* Do inversion or not */
	int ksstartflag[MAX_NUM_MASS];	/* what to do for beginning propagators */
	int kssaveflag[MAX_NUM_MASS];	/* what to do for saving propagators */
	char ksstartfile[MAX_NUM_MASS][MAXFILENAME];
	char kssavefile[MAX_NUM_MASS][MAXFILENAME];
#else
#ifdef IMP_GFIX
	Real u0;	/* tadpole parameter */
#endif
#endif
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
	char stringLFN[MAXFILENAME];	/** ILDG LFN if applicable ***/
}  params;

#endif /* _PARAMS_H */
