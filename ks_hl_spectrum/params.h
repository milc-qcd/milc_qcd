#ifndef _PARAMS_H
#define _PARAMS_H

/* structure for passing simulation parameters to each node */
typedef struct {
        int stopflag;	/* 1 if it is time to stop */
   /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */
        char job_id[MAXFILENAME];

   /*  REPEATING BLOCK */
	int startflag;  /* what to do for beginning lattice */
	int saveflag;   /* what to do with lattice at end */
	int niter; 	/* maximum number of c.g. iterations */
	double rsqmin,rsqprop;  /* for deciding on convergence */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
	char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
        char ensemble_id[MAXFILENAME];
        int sequence_number;

	char startfile_w[MAX_KAP][MAXFILENAME];
	int startflag_w[MAX_KAP];
        char src_label_w[MAX_KAP][16];
	char savefile_w[MAX_KAP][MAXFILENAME];
	int saveflag_w[MAX_KAP];
        char savefile_c[MAXFILENAME];
        int saveflag_c;

        char smearfile[MAX_KAP][MAXFILENAME];
        char sink_label[MAX_KAP][16];
        char start_ks_prop_file[MAXFILENAME];
        int ks_prop_startflag;
        Real mass;
        char mass_label[32];
        int num_kap;	/* number of kappa's */
        int log_correlators;
        int num_smear; /* number of smearings */
        Real kap[MAX_KAP];	/* kappa values for multiple propagators */	
        char kap_label[MAX_KAP][32];
        Real d1[MAX_KAP];   /*rotation parameter*/

  /* For the baryon code */
        int  start_ks_strange_flag[MAX_STRANGE];
        char start_ks_strange_file[MAX_STRANGE][MAXFILENAME];
        int  start_ks_light_flag[MAX_LIGHT];
        char start_ks_light_file[MAX_LIGHT][MAXFILENAME];
        Real m_light[MAX_LIGHT];
        Real m_strange[MAX_LIGHT];
        int num_light;
        int num_strange;
}  params;

#endif /* _PARAMS_H */
