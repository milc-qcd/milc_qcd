#ifndef _PARAMS_H
#define _PARAMS_H

/* structure for passing simulation parameters to each node */
typedef struct {
        int stopflag;	/* 1 if it is time to stop */
   /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */
  //int iseed;	/* for random numbers */
  	int nflavors;	/* the number of flavors */
   /*  REPEATING BLOCK */
  //	int warms;	/* the number of warmup trajectories */
  //	int trajecs;	/* the number of real trajectories */
  //	int steps;	/* number of steps for updating */
  //	int propinterval;     /* number of trajectories between measurements */
	int startflag;  /* what to do for beginning lattice */
	int saveflag;   /* what to do with lattice at end */
	
	float rsqmin,rsqprop;  /* for deciding on convergence */
  //	float epsilon;	/* time step */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
        char ensemble_id[MAXFILENAME];
        int sequence_number;

	char startfile_w[MAX_KAP][MAXFILENAME];
	int startflag_w[MAX_KAP];
	char savefile_w[MAX_KAP][MAXFILENAME];
	int saveflag_w[MAX_KAP];

        char a_startfile_w[MAXFILENAME];
        char smearfile[MAX_KAP][MAXFILENAME];
 
        int num_kap;	/* number of kappa's */
        int num_smear; /* number of smearings */
        float kap[MAX_KAP];	/* kappa values for multiple propagators */
        float d1[MAX_KAP];   /*rotation parameter*/
        int format[MAX_KAP]; /* propagator format. milc =0, fermi=1*/	
        float a_d1;   /*rotation parameter*/
        int a_format; 
}  params;

#endif /* _PARAMS_H */
