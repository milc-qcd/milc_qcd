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
	int warms;	/* the number of warmup trajectories */
	int trajecs;	/* the number of real trajectories */
	int steps_over;	/* number of steps for overrelaxed updating */
	int steps_update;	/* number of steps for Metropolis */
	int propinterval;     /* number of trajectories between measurements */
	int startflag;  /* what to do for beginning lattice */
	int fixflag;    /* whether to gauge fix */
	int saveflag;   /* what to do with lattice at end */
	Real beta;	/* gauge coupling */
	Real C;	/* mass/mu parameter */
	Real epsilon;	/* time step */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
	char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
}  params;


#endif /* _PARAMS_H */
