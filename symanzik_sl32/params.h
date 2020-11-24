#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"	/* For MAXFILENAME */

/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;	/* 1 if it is time to stop */
    /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */
	int iseed;	/* for random numbers */
    /*  REPEATING BLOCK */
	int warms;	/* the number of warmup trajectories */
	int trajecs;	/* the number of real trajectories */
	int steps;	/* number of steps for updating */
	int stepsQ;	/* number of steps for qhb */
	int propinterval;     /* number of trajectories between measurements */
	int startflag;  /* what to do for beginning lattice */
	int saveflag;   /* what to do with lattice at end */
#ifndef ANISOTROPY
	Real beta;	/* gauge coupling */
	Real u0;	/* <Tr(U_p)>^{1/4} */
#else 
        short ani_dir; /* direction of anisotropy */
        Real beta[2]; /* beta[0] - 3d-isotropic, beta[1] - anisotropic */
/* The tadpole factor should be consistently defined as Real u0[2]; 
   with the elements u0[0] - 3d-isotropic, u0[1] - anisotropic.
   However, in the following we keep a single tadpole 
   factor u0 since we only aim for anisotropic simulations 
   using u0=1.0. 
*/
        Real u0; 
#endif
	Real epsilon;	/* time step */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
	char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
}  params;

#endif	/* _PARAMS_H */
