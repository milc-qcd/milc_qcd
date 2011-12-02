#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/generic_quark_types.h"	/* For quark_source */

/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;   /* 1 if it is time to stop */
   /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */
#ifdef RANDOM
	int iseed;	/* for random numbers */
#endif
   /*  REPEATING BLOCK */
	int startflag;  /* what to do for beginning lattice */
	int fixflag;    /* whether to gauge fix */
	int saveflag;   /* what to do with lattice at end */
	int startflag_w[MAX_MASSES];	/* what to do for beginning wilson vector */
	int saveflag_w[MAX_MASSES];	/* what to do for saving wilson vector */
        quark_source wqs;  /* source parameters */
#ifdef H0INV
	int startflag_w3[MAX_MASSES];	/* what to do for beginning 3pt wilson vector */
	int saveflag_w3[MAX_MASSES];	/* what to do for saving 3pt wilson vector */
#endif
	int num_masses;	/* number of masses */
	Real mass[MAX_MASSES];	/* masses values for multiple propagators */
	Real resid[MAX_MASSES];	/* residue for invertion convergence */
	Real resid2[MAX_MASSES];	/* residue for invertion convergence */
	int niter; 	/* maximum number of c.g. iterations */
	int nrestart; 	/* maximum number of c.g. restarts */
	char startfile[MAXFILENAME];
	char startfile_w[MAX_MASSES][MAXFILENAME];
#ifdef H0INV
	char startfile_w3[MAX_MASSES][MAXFILENAME];
#endif
	char savefile[MAXFILENAME];
        char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
	char savefile_w[MAX_MASSES][MAXFILENAME];
	char scratchstem_w[MAXFILENAME];
#ifdef H0INV
	char savefile_w3[MAX_MASSES][MAXFILENAME];
	char scratchstem_w3[MAXFILENAME];
#endif
	int scratchflag;
    /* new stuff for i/o... */
	int out_hr0_flag,in_hr0_flag,out_hov_flag,in_hov_flag;
	char in_hov[MAXFILENAME];
	char out_hov[MAXFILENAME];
	char in_hr0[MAXFILENAME];
	char out_hr0[MAXFILENAME];
    /* overlap parameters */
	Real R0;
	Real prec_sign;
	Real zolo_min;
	Real zolo_max;
	Real scalez;
	int maxcg_inner;
	Real resid_inner;
	Real resid_inner_h;
	int nsw;
	int topology;
	/* gauge parameters */
		
#ifdef EIG
	int Nvecs_h0r0;	/* number of eigenvectors for h0(-R0) */
	int Nvecs_h0;	/* number of trial eigenvectors of h0 */
	int Nvecs_hov;	/* number of eigenvectors of H_ov */
        Real eigenvec_qual;
	Real eigenval_tol;	/* Tolerance for the eigenvalue computation */
	Real error_decr;	/* error decrease per Rayleigh minimization */
	Real eigenval_tol_acc;	/* Tolerance for the eigenvalue computation */
	Real error_decr_acc;	/* error decrease per Rayleigh minimization */
	int MaxIter;	/* max  Rayleigh iterations */
	int Maxr0Iter;	/* max  Rayleigh iterations for h(-r0) */
	int Restart;	/* Restart  Rayleigh every so many iterations */
	int Kiters;	/* Kalkreuter iterations */
	double *eigVal_buf;
#endif
	int ndone;	/* Number of (random) sources already done */


        Real beta;     /* gauge coupling */
	Real u0;

#ifdef IMAGISO
  Real delta_iso;
#endif
	int nsmear;
#ifdef DOMAINX
  int cut_x,cut_y,cut_z,cut_t;
#endif

}  params;

#endif /* _PARAMS_H */
