#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"	/* For MAXFILENAME */

/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;	/* 1 if it is time to stop */
    /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */
	int iseed;	/* for random numbers */
	int nflavors;	/* the number of flavors */
    /*  REPEATING BLOCK */
	int warms;	/* the number of warmup trajectories */
	int trajecs;	/* the number of real trajectories */
	int steps;	/* number of steps for updating */
	int startflag;  /* what to do for beginning lattice */
	int saveflag;   /* what to do with lattice at end */
	Real beta,mass; /* gauge coupling, quark mass */
	int niter; 	/* maximum number of c.g. iterations */
	int nrestart; 	/* maximum number of c.g. restarts */
	int bc_flag; 	/* gauge boundary condition flag */
	Real rsqmin;	/* for deciding on convergence */
	Real epsilon;	/* time step */
	Real gamma_rv;	/* Reweigh dS/deta coefficient */
	Real ferm_phas[3];	/* fermion phase factors */
	int n_sxw;	/* n for Sexton-Weingarten update */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
	char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
}  params;

#endif	/* _PARAMS_H */
