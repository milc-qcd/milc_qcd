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
#ifdef ONEMASS
	int nflavors;	/* the number of flavors */
#else
	int nflavors1;	/* the number of flavors of first type */
	int nflavors2;	/* the number of flavors of second type */
#endif
    /*  REPEATING BLOCK */
	int warms;	/* the number of warmup trajectories */
	int trajecs;	/* the number of real trajectories */
	int steps;	/* number of steps for updating */
	int propinterval;     /* number of trajectories between measurements */
	Real beta;      /* gauge coupling */
#ifdef ONEMASS
	Real mass;      /*  quark mass */
#else
	Real mass1,mass2; /*  quark masses */
#endif
	Real u0; /* tadpole parameter */
	int niter; 	/* maximum number of c.g. iterations */
	Real rsqmin,rsqprop;  /* for deciding on convergence */
	Real epsilon;	/* time step */
        char spectrum_request[MAX_SPECTRUM_REQUEST];   /* request list for spectral measurements */
        int source_start, source_inc, n_sources; /* source time and increment */
        int spectrum_multimom_nmasses; 
        Real spectrum_multimom_low_mass;
        Real spectrum_multimom_mass_step;
        int fpi_nmasses;
        Real fpi_mass[MAX_FPI_NMASSES];
	int startflag;  /* what to do for beginning lattice */
	int saveflag;   /* what to do with lattice at end */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
}  params;

#endif /* _PARAMS_H */
