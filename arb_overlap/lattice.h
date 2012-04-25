#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_wilson.h"
#include "../include/generic_quark_types.h"
#include "../include/random.h"   /* For double_prn */
#include "../include/macros.h"   /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#include "params.h"

/* Begin definition of site structure */

#include "../include/su3.h"



/* The lattice is an array of sites.  */
typedef struct {
    /* The first part is standard to all programs */
	/* coordinates of this site */
	short x,y,z,t;
	/* is it even or odd? */
	char parity;
	/* my index in the array */
	int index;
#ifdef SITERAND
	/* The state information for a random number generator */
	double_prn site_prn;
	/* align to double word boundary (kludge for Intel compiler) */
	int space1;
#endif

    /* Now come the physical fields, program dependent */
	/* gauge field */
	su3_matrix link[4];


	su3_matrix tempmat1,tempmat2;
  /* INV for OUTER CG */
#ifdef INV
/* these used in outer inversion */
	wilson_vector chi0;	/* another CG vector */
	wilson_vector pm0;	/* another CG vector */
	wilson_vector mpm0;	/* another CG vector */
#endif
#ifdef RANDOM
	wilson_vector g_rand;	/* to save the random source */
#endif
/* these are used in hdelta0 and in the build routines */
	wilson_vector psi;	/* solution vector */
	wilson_vector chi;	/* source vector */

	wilson_vector r;	/* residue:  */
	/* wilson vector (temporary used in delta0) */
	wilson_vector htmp[2];
	su3_matrix staple;

#ifndef HASENBUSCH
        wilson_vector ransrc1, hmc_src1;
        wilson_vector ransrc2, hmc_src2;
#endif

#ifdef PTP
  complex pmes_prop[MAX_MASSES][4];
  complex print_var;
#endif

#ifdef RHAIR
  complex wave_fn;        /* for the wave function, and its FFT */
#endif


} site;

/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/* The following are global scalars */
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
#ifdef RANDOM
EXTERN	int iseed;		/* random number seed */
#endif
EXTERN	int niter,nrestart;

/* quark masses & associated CG residures */
EXTERN	Real mass[MAX_MASSES],resid[MAX_MASSES],resid_acc[MAX_MASSES];
EXTERN	int num_masses;		/* max number of masses <= MAX_MASSES */

/* parameters of the kernel action */
EXTERN	Real lambda[5],rho[5];
EXTERN Real clover_term;

/* temporary vectors used in delta0.c */
EXTERN    wilson_vector* delta0_tmp1;
EXTERN    half_wilson_vector* delta0_htmp1;
EXTERN    half_wilson_vector* delta0_htmp2;
EXTERN    half_wilson_vector* delta0_htmp3;

/* parameters for the Zolotarev approx. */
EXTERN	Real *shift, *coeff,scalez,resid_inner,resid_inner_h,R0,FSUcoeff0;
EXTERN   Real resid_inner_save;
EXTERN	int Norder,maxcg_inner,kind_of_h0,chirality_flag;
EXTERN  Real prec_sign;
EXTERN    Real zolo_min, zolo_max, zolo_min_save, zolo_max_save;

#ifdef MINN
EXTERN Real resid_inner_run;
EXTERN int do_minn;
#endif
#ifdef MOUT
EXTERN int do_mout;
#endif

/* for the outer CG multi mass*/
#ifdef INV
EXTERN	Real *shift0;
#ifdef MULTI
EXTERN Real *coeff0;
#endif
#endif

#ifdef EIG
/* Eigenvalue related global variables */
EXTERN	int Nvecs_h0r0, Nvecs_h0, Nvecs_hov;	/* number of eigenvectors:
     ``inner'' eigenvals (for h(-R0)), trial h(0), H_ov */
EXTERN	Real eigenval_tol;	/* Tolerance for the eigenvalue computation */
EXTERN	Real error_decr;	/* error decrease per Rayleigh minimization */
EXTERN	Real eigenval_tol_low;	/* Tolerance for the eigenvalue computation */
EXTERN	Real error_decr_low;	/* error decrease per Rayleigh minimization */
EXTERN	Real eigenval_tol_high;	/* Tolerance for the eigenvalue computation */
EXTERN	Real error_decr_high;	/* error decrease per Rayleigh minimization */
EXTERN	int MaxIter;	/* max  Rayleigh iterations */
EXTERN	int Maxr0Iter;	/* max  Rayleigh iterations for h(-r0) */
EXTERN	int Restart;	/* Restart  Rayleigh every so many iterations */
EXTERN	int Kiters;	/* Kalkreuter iterations */
EXTERN  wilson_vector  **eigVec0, **eigVec;
EXTERN  double *eigVal0, *eigVal;
EXTERN  Real eigValcut;
#endif


EXTERN  params par_buf;
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
/* some file names ---- deprectated */
EXTERN	char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;	/* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN	int saveflag;	/* do with lattice: 1=save; */
EXTERN	char startfile_w[MAX_MASSES][MAXFILENAME],
	     savefile_w[MAX_MASSES][MAXFILENAME];
EXTERN	int startflag_w[MAX_MASSES];	/* beginning wilson:
				RELOAD_ASCII, RELOAD_BINARY, RELOAD_PARALLEL */
EXTERN	int saveflag_w[MAX_MASSES];	/* save propagator: SAVE_ASCII,
				SAVE_BINARY, SAVE_PARALLEL */
#ifdef H0INV
EXTERN	char startfile_w3[MAX_MASSES][MAXFILENAME],
	     savefile_w3[MAX_MASSES][MAXFILENAME];
EXTERN	int startflag_w3[MAX_MASSES];	/* beginning wilson:
				RELOAD_ASCII, RELOAD_BINARY, RELOAD_PARALLEL */
EXTERN	int saveflag_w3[MAX_MASSES];	/* save propagator: SAVE_ASCII,
				SAVE_BINARY, SAVE_PARALLEL */
#endif

EXTERN	char in_hov[MAXFILENAME], out_hov[MAXFILENAME];
EXTERN	char in_hr0[MAXFILENAME], out_hr0[MAXFILENAME];
EXTERN	int out_hr0_flag, in_hr0_flag, out_hov_flag, in_hov_flag;
EXTERN	int ndone;	/* Number of (random) sources already done */

EXTERN	int total_iters;
/* boundary-flip switch, made consistent for higher-rep code */
EXTERN  int current_boundary;
EXTERN  int current_boundary_x;
EXTERN  int current_boundary_z;

/* Are these needed?? */
/* Further source description */
#ifdef CONTROL
int n_spins = 4;		/* Number of spins generated */
int source_loc[4] = {0,0,0,0};	/* Source location */
int spins[4] = {0,1,2,3};	/* List of spins generated */
#else
extern int n_spins;
extern int source_loc[4];
extern int spins[4];
#endif


/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN	int this_node;		/* node number of this node */


EXTERN quark_source wqs;
EXTERN quark_source wqstmp;	/* Temporary */
EXTERN dirac_clover_param dcptmp;
 
 
/* information for derivative and link terms */
EXTERN int offset[NLINK][4];
EXTERN int diroff1[8];
EXTERN int signoff1[8];
EXTERN int diroff[8][8];
EXTERN int signoff[8][8];
EXTERN int label[NLINK],off_max;
/* information for gathers (gather-offset) */
EXTERN int goffset[2*NLINK];


EXTERN	gauge_file *startlat_p;
EXTERN	gauge_file *savelat_p;


/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN  double_prn node_prn ;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

EXTERN su3_matrix *ape_links;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 32	/* Number of generic pointers (usually 10) */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN char ** gen_pt[N_POINTERS];

#ifdef NHYP
/* hard-coding of alpha_smear controlled by defines.h */
#ifndef HARD_CODE_SMEAR
EXTERN Real alpha_smear[3];
#endif
/* prototype fields for NHYP links */
EXTERN su3_matrix *gauge_field[4];
EXTERN su3_matrix *gauge_field_thin[4];

EXTERN su3_matrix *hyplink1[4][4]; /* Needed for other stuff, too */
EXTERN su3_matrix *hyplink2[4][4];
EXTERN su3_matrix *Sigma[4];
EXTERN su3_matrix *SigmaH[4];
EXTERN su3_matrix *SigmaH2[4][4];

EXTERN su3_matrix *Staple1[4][4];
EXTERN su3_matrix *Staple2[4][4];
EXTERN su3_matrix *Staple3[4];

EXTERN su3_matrix *LambdaU[4];
EXTERN su3_matrix *Lambda1[4];
EXTERN su3_matrix *Lambda2[4];
/* EXTERN su3_matrix *save_gf[4]; */

/* renamed from tempmat1, for force_nhyp
   also used to compute gauge action */
EXTERN su3_matrix *tempmat_nhyp1;
EXTERN su3_matrix *tempmat_nhyp2;

#endif /* NHYP */



#endif 

#ifdef FIELD
/* field major storage DON't FORGET to MALLOC somewhere */
EXTERN su3_matrix *t_blocked_link2;
EXTERN	triangular *t_clov;
EXTERN	diagonal *t_clov_diag;
#endif






/* parameters of the action */
EXTERN    int   nflavors;
EXTERN    Real beta;
EXTERN    Real u0;


EXTERN    int   nsmear;

EXTERN    long int   ndelta0;

/* variables to keep track of the topology */
EXTERN    int current_topology;
EXTERN    int ct_save;

EXTERN    int repro_cg;
EXTERN    int md_time;
EXTERN    su3_matrix gauge_force_sum[4];

#ifdef MOUT
EXTERN  Real cgcutoff;
#endif

#ifdef PRECON
EXTERN int do_precon, MaxCG_precon;
#endif

EXTERN    Real eigenvec_quality;

#ifdef BUILDMODES
EXTERN int modeflag,Nvecs_k_hold,Nvecs_k_held,use_modes;
EXTERN double *MyeigVal;
EXTERN  half_wilson_vector **MyeigVec;
EXTERN  half_wilson_vector **MyDeigVec;
EXTERN int Nvecs_hov_hold,Nvecs_hov_build;
#endif

#ifdef IMAGISO
EXTERN Real delta_iso;
#endif

#ifdef DOMAINX
EXTERN int cut_x,cut_y,cut_z,cut_t;
#endif

/* _LATTICE_H */
