#ifndef _PARAMS_H
#define _PARAMS_H
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
	int propinterval;     /* number of trajectories between measurements */
	int startflag;  /* what to do for beginning lattice */
	int saveflag;   /* what to do with lattice at end */
	Real beta,mass; /* gauge coupling, quark mass */
	int niter; 	/* maximum number of c.g. iterations */
	Real rsqmin,rsqprop;  /* for deciding on convergence */
	Real epsilon;	/* time step */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
	char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
        char ensemble_id[MAXFILENAME];
        int sequence_number;
}  params;

#endif /* _PARAMS_H */
