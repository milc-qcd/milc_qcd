	complex ploop_fuzz; /*  for Polyakov loop correlations */

/* stuff for baryon density correlations */
	complex ploop; /* For Polyakov loop correlations */
	su3_vector propmat[3];	/* For three source colors */
	su3_vector propmat2[3]; /* for second type of source */
#define NBPRAND 20
#define NBPAVRG 11
	complex bdensum;
	complex bdens[NBPRAND];  /* Temporary for computing two-trace correlations */
	complex avg[NBPAVRG];
	complex gath1[NBPAVRG];
	complex gath2[NBPAVRG];

