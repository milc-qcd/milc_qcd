#ifndef _PARAMS_H
#define _PARAMS_H

/* structure for passing simulation parameters to each node */
typedef struct {
        int stopflag;	/* 1 if it is time to stop */
   /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */

   /*  REPEATING BLOCK */
  //	int warms;	/* the number of warmup trajectories */
  //	int trajecs;	/* the number of real trajectories */
  //	int steps;	/* number of steps for updating */
  //	int propinterval;     /* number of trajectories between measurements */
	int startflag;  /* what to do for beginning lattice */
	int saveflag;   /* what to do with lattice at end */
  //	double beta,mass; /* gauge coupling, quark mass */
	int niter; 	/* maximum number of c.g. iterations */
	double rsqmin,rsqprop;  /* for deciding on convergence */
  //	double epsilon;	/* time step */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
        char ensemble_id[MAXFILENAME];
        int sequence_number;

	char startfile_w[MAX_KAP][MAXFILENAME];
	int startflag_w[MAX_KAP];
	char savefile_w[MAX_KAP][MAXFILENAME];
	int saveflag_w[MAX_KAP];

        char smearfile[MAX_KAP][MAXFILENAME];
        char start_ks_prop_file[MAXFILENAME];
        int ks_prop_startflag;
        int num_kap;	/* number of kappa's */
        int num_smear; /* number of smearings */
        Real kap[MAX_KAP];	/* kappa values for multiple propagators */	
        Real d1[MAX_KAP];   /*rotation parameter*/
        int format[MAX_KAP]; /* propagator format. milc =0, fermi=1*/
        int a_format;
}  params;

#endif /* _PARAMS_H */
