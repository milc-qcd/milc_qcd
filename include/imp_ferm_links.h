#ifndef _IMP_FERM_LINKS_H
#define _IMP_FERM_LINKS_H

#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/comdefs.h"
#include "../include/generic_quark_types.h"
#include <quark_action.h>  /* Defines FERM_ACTION */

/*-----------------------------------------------------------------*/
/* Define imp_ferm_links_t */

#ifdef HAVE_QOP

#include <qop.h>

#if FERM_ACTION == HISQ

#include "../include/fn_links_qop.h"
#include "../include/hisq_links_qop.h"
#include "../include/fermion_links_qop.h"
#define imp_ferm_links_t fn_links_qop_t
#define create_imp_ferm_links create_fn_links_qop
#define destroy_imp_ferm_links destroy_fn_links_qop


#else /* Not HISQ */

#include "../include/fn_links_qop.h"
#include "../include/fermion_links_qop.h"
#define imp_ferm_links_t fn_links_qop_t
#define create_imp_ferm_links create_fn_links_qop
#define load_imp_ferm_links load_fn_links_qop
#define destroy_imp_ferm_links destroy_fn_links_qop

#endif

#else /* Not HAVE_QOP */

#include "../include/fn_links.h"
#include "../include/eo_links.h"

#if ( FERM_ACTION == HISQ ) || ( FERM_ACTION == FN_TYPE )

#define imp_ferm_links_t fn_links_t
#define create_imp_ferm_links create_fn_links
#define load_imp_ferm_links load_fn_links
#define destroy_imp_ferm_links destroy_fn_links

#else

#define imp_ferm_links_t eo_links_t 
#define create_imp_ferm_links create_eo_links
#define load_imp_ferm_links load_eo_links
#define destroy_imp_ferm_links destroy_eo_links

#endif

#endif

/*-----------------------------------------------------------------*/
/* Prototypes for functions using imp_ferm_links */

int congrad( int niter, int nrestart, Real rsqmin, int parity, Real *rsq );
int ks_congrad( field_offset src, field_offset dest, Real mass,
		int niter, int nrestart, Real rsqmin, int prec, 
		int parity, Real *rsq, imp_ferm_links_t *fn );

int ks_congrad_field( su3_vector *src, su3_vector *dest, 
		      quark_invert_control *qic, Real mass,
		      imp_ferm_links_t *fn);

int ks_congrad_block_field( int nsrc, su3_vector **src, su3_vector **dest, 
			    quark_invert_control *qic, Real mass,
			    imp_ferm_links_t *fn);

int ks_congrad_field_cpu( su3_vector *src, su3_vector *dest, 
			  quark_invert_control *qic, Real mass,
			  imp_ferm_links_t *fn);

int ks_congrad_site( field_offset src, field_offset dest, 
		     quark_invert_control *qic, Real mass,
		     imp_ferm_links_t *fn);

int ks_congrad_parity_cpu( su3_vector *t_src, su3_vector *t_dest, 
			   quark_invert_control *qic, Real mass,
			   imp_ferm_links_t *fn);


#ifdef USE_CG_GPU
#if defined(HAVE_QUDA) || defined(HAVE_GRID) 

#define ks_congrad_parity ks_congrad_parity_gpu
#define ks_congrad_block_parity ks_congrad_block_parity_gpu

#elif USE_CG_QPHIX

#define ks_congrad_parity ks_congrad_parity_qphix
#define ks_congrad_block_parity ks_congrad_block_parity_qphix

#endif

#else

#define ks_congrad_parity ks_congrad_parity_cpu
#define ks_congrad_block_parity ks_congrad_block_parity_cpu

#endif

int ks_congrad_parity( su3_vector *t_src, su3_vector *t_dest, 
		       quark_invert_control *qic, Real mass,
		       imp_ferm_links_t *fn);

int ks_congrad_block_parity( int nsrc, su3_vector **t_src, su3_vector **t_dest, 
			     quark_invert_control *qic, Real mass,
			     imp_ferm_links_t *fn);

int ks_congrad_two_src(	/* Return value is number of iterations taken */
    field_offset src1,    /* source vector (type su3_vector) */
    field_offset src2,
    field_offset dest1,	/* solution vectors */
    field_offset dest2,
    Real mass1,
    Real mass2,
    int niter,		/* maximal number of CG interations */
    int nrestart,       /* maximum number of restarts */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* internal precision for the inversion */
    int parity,		/* parity to be worked on */
    Real  *final_rsq_ptr, 	/* final residue squared */
    imp_ferm_links_t *fn       /* Storage for fermion links */
    );

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[]);
void cleanup_dslash_temps(void);

void dslash_fn_site( field_offset src, field_offset dest, int parity,
		     imp_ferm_links_t *fn);
void dslash_fn_site_special( field_offset src, field_offset dest,
			     int parity, msg_tag **tag, int start,
			     imp_ferm_links_t *fn);
void ddslash_fn_du0_site( field_offset src, field_offset dest, int parity,
			  imp_ferm_links_t *fn, imp_ferm_links_t *fn_dmdu0);

void dslash_fn_field( su3_vector *src, su3_vector *dest, int parity,
		      imp_ferm_links_t *fn);
void dslash_fn_field_special(su3_vector *src, su3_vector *dest,
			     int parity, msg_tag **tag, int start,
			     imp_ferm_links_t *fn);
void ddslash_fn_du0_field( su3_vector *src, su3_vector *dest, int parity,
			   imp_ferm_links_t *fn, imp_ferm_links_t *fn_dmdu0);

void dslash_fn_dir(const su3_vector * const src, su3_vector *dest, int parity,
		   imp_ferm_links_t *fn, int dir, int fb, 
		   Real wtfat, Real wtlong);


void dslash_eo_site( field_offset src, field_offset dest, int parity,
		     imp_ferm_links_t *fn);


/* The following three do not exist yet (3/05 -CD) */
void dslash_eo_site_special( field_offset src, field_offset dest,
			     int parity, msg_tag **tag, int start,
			     imp_ferm_links_t *fn );
void dslash_eo_field( su3_vector *src, su3_vector *dest, int parity,
		      imp_ferm_links_t *fn);
void dslash_eo_field_special( su3_vector *src, su3_vector *dest,
			      int parity, msg_tag **tag, int start,
			      imp_ferm_links_t *fn);

int congrad_ks(            /* Return value is number of iterations taken */
     field_offset src,       /* type su3_vector* (preloaded source) */
     field_offset dest,      /* type su3_vector*  (answer and initial guess) */
     quark_invert_control *qic, /* inverter control */
     void *dmp               /* parameters defining the Dirac matrix */
     );

int ks_invert( /* Return value is number of iterations taken */
    field_offset src,   /* type su3_vector or multi_su3_vector 
			   (preloaded source) */
    field_offset dest,  /* type su3_vector or multi_su3_vector 
			   (answer and initial guess) */
    int (*invert_func)(field_offset src, field_offset dest,
		       quark_invert_control *qic,
		       Real mass, imp_ferm_links_t *fn),
    quark_invert_control *qic, /* inverter control */
    Real mass,
    imp_ferm_links_t *fn
    );

int ks_invert_ksqs( /* Return value is number of iterations taken */
    quark_source *ksqs, /* source parameters */
    int (*source_func_field)(su3_vector *src, 
			      quark_source *ksqs),  /* source function */
    su3_vector *dest,  /* answer and initial guess */
    int (*invert_func)(su3_vector *src, su3_vector *dest,
		       quark_invert_control *qic, Real mass, 
		       imp_ferm_links_t *fn),
    quark_invert_control *qic, /* inverter control */
    Real mass,
    imp_ferm_links_t *fn
		    );

/* in ks_multicg.c */
const char *ks_multicg_opt_chr( void );

int ks_multicg_field(   /* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* KS parameters with offsets defined */
    int num_offsets,	/* number of offsets */
    quark_invert_control qic[], /* inversion parameters */
    imp_ferm_links_t *fn /* Storage for fat and Naik links */
    );

int ks_multicg_offset_field_cpu(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* the offsets */
    int num_offsets,	/* number of offsets */
    quark_invert_control qic[], /* inversion parameters */
    imp_ferm_links_t *fn      /* Storage for fat and Naik links */
    );

int ks_multicg_offset_field_gpu(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* the offsets */
    int num_offsets,	/* number of offsets */
    quark_invert_control qic[], /* inversion parameters */
    imp_ferm_links_t *fn      /* Storage for fat and Naik links */
    );

int ks_multicg_offset_field_qphix(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* the offsets */
    int num_offsets,	/* number of offsets */
    quark_invert_control qic[], /* inversion parameters */
    imp_ferm_links_t *fn      /* Storage for fat and Naik links */
    );

#ifdef USE_CG_GPU
#if defined(HAVE_GRID) || defined(HAVE_QUDA)
#define ks_multicg_offset_field ks_multicg_offset_field_gpu
#elif USE_CG_QPHIX
#define ks_multicg_offset_field ks_multicg_offset_field_qphix
#endif
#else
#define ks_multicg_offset_field ks_multicg_offset_field_cpu
#endif

int ks_multicg_mass_field(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    ks_param *ksp,	/* the KS parameters, including masses */
    int num_masses,	/* number of masses */
    quark_invert_control qic[],  /* inversion parameters */
    imp_ferm_links_t *fn /* Storage for fat and Naik links */
    );


int ks_multicg_mass_site(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *masses,	/* the masses */
    int num_masses,	/* number of masses */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* internal precision for inversion (ignored) */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr, /* final residue squared */
    imp_ferm_links_t *fn /* Storage for fat and Naik links */
);


/* Unlike the procedures above, this one computes dst = M^{-1} src on
   EVENANDODD sites */

int mat_invert_multi(
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **dst,	/* solution vectors */
    ks_param *ksp,	/* KS parameters, including masses */
    int num_masses,	/* number of masses */
    quark_invert_control qic[],  /* inversion parameters */
    imp_ferm_links_t *fn /* Storage for fat and Naik links */
    );

/*
  Return the most recent fermion link field passed to QUDA
  (defined in generic_ks/ks_multicg_offset_gpu.c)
*/
imp_ferm_links_t* get_fn_last();

/*
  Update the fermion link field passed to QUDA
  (defined in generic_ks/ks_multicg_offset_gpu.c)
*/
void set_fn_last(imp_ferm_links_t *fn_last_new);

/* eigen_stuff*.c */

void restore_eigVec(int Nvecs, Real *eigVal, su3_vector **eigVec, int parity,
		    imp_ferm_links_t *fn);
typedef struct {
  int norder ; /* Order of the preconditioning polynomial */
  int which_poly; /* Polynomial selection */
  double minE ; /* Lower end of eigenvalue exclusion window */
  double maxE ; /* Upper end of eigenvalue exclusion window */
  double poly_param_1;
  double poly_param_2;
  double eigmax;
} ks_eigen_poly;

#if defined(HAVE_PRIMME)
#define ks_eigensolve ks_eigensolve_PRIMME
typedef struct {
  int Nvecs ; /* number of eigenvectors */
  int Nvecs_in ; /* number of input starting eigenvectors */
  Real tol ; /* Tolerance for the eigenvalue computation */
  Real error_decr ; /* error decrease per Rayleigh minimization */
  int MaxIter ; /* max  Rayleigh iterations */
  int Restart ; /* Restart  Rayleigh every so many iterations */
  int Kiters ; /* Kalkreuter iterations */
  int parity; 
  ks_eigen_poly poly; /* Preconditioning polynomial */
} ks_eigen_param;
#elif defined(HAVE_ARPACK)
#define ks_eigensolve ks_eigensolve_ARPACK
typedef struct {
  int Nvecs ; /* number of eigenvectors */
  int Nvecs_in ; /* number of input starting eigenvectors */
  Real tol ; /* Tolerance for the eigenvalue computation */
  Real error_decr ; /* error decrease per Rayleigh minimization */
  int MaxIter ; /* max  Rayleigh iterations */
  int nArnoldi ; /* Number of Arnoldi vectors generated in each iteration */
  int parity; 
  ks_eigen_poly poly; /* Preconditioning polynomial */
} ks_eigen_param;
#elif defined(USE_EIG_GPU) && defined(HAVE_GRID)
#define ks_eigensolve ks_eigensolve_Grid
typedef struct {
  int Nvecs; /* number of eigenvectors */
  int Nvecs_in; /* number of input starting eigenvectors */
  Real tol; /* Tolerance for the eigenvalue computation */
  int MaxIter; /* max implicit restarting iterations */
  int Nrestart; /* Lanczos restarts from this number of eigenvalues */
  int Nmax; /* Lanczos iteration stops here and restart from Nrestart */
  int reorth_period; /* Reorthogonalization period */
  ks_eigen_poly poly; /* Preconditioning polynomial */
  char diagAlg[10];
  int parity; 
} ks_eigen_param;
#elif defined(USE_EIG_GPU) && defined(HAVE_QUDA)
#define ks_eigensolve ks_eigensolve_QUDA
typedef struct {
  int Nvecs ; /* number of eigenvectors */
  int Nvecs_in; /* number of input starting eigenvectors */
  Real tol; /* Tolerance for the eigenvalue computation */
  int MaxIter; /* max restarting iterations */
  int Nrestart; /* Lanczos restarts from this number of eigenvalues */
  int Nkr; /* size of the Krylov subspace */
  ks_eigen_poly poly; /* Preconditioning polynomial */
  int blockSize; /* block size for block variant eigensolvers */
  int partfile; /* Whether to save in partfile or not */
  int eigPrec; /* Run the eigensolver in this precision */
  int batchedRotate; /* Size of the rotation space to use for solver */
  int parity;
  double tol_restart; 
} ks_eigen_param;
#elif defined(HAVE_QDP)
#define ks_eigensolve Kalkreuter
typedef struct {
  int Nvecs ; /* number of eigenvectors */
  int Nvecs_in ; /* number of input starting eigenvectors */
  Real tol ; /* Tolerance for the eigenvalue computation */
  Real RelTol ;
  int MaxIter ; /* max  Rayleigh iterations */
  int Restart ; /* Restart  Rayleigh every so many iterations */
  int Kiters ; /* Kalkreuter iterations */
  int parity; 
} ks_eigen_param;
#else
#define ks_eigensolve ks_eigensolve_Kalkreuter_Ritz
typedef struct {
  int Nvecs ; /* number of eigenvectors */
  int Nvecs_in ; /* number of input starting eigenvectors */
  Real tol ; /* Tolerance for the eigenvalue computation */
  Real error_decr ; /* error decrease per Rayleigh minimization */
  int partfile; /* Whether to save in partfile or not */
  int MaxIter ; /* max  Rayleigh iterations */
  int Restart ; /* Restart  Rayleigh every so many iterations */
  int Kiters ; /* Kalkreuter iterations */
  int parity; 
  Real tol_restart;
} ks_eigen_param;
#endif

#define EIGCG 2

int Rayleigh_min(su3_vector *vec, su3_vector **eigVec, Real Tolerance, 
		 Real RelTol, int Nvecs, int MaxIter, int Restart, 
		 ks_eigen_param *eigen_param, imp_ferm_links_t *fn);
int ks_eigensolve_Kalkreuter_Ritz(su3_vector **eigVec, Real *eigVal, 
				  ks_eigen_param *eigen_param, int init );
int ks_eigensolve_PRIMME(su3_vector **eigVec, Real *eigVal,
				  ks_eigen_param *eigen_param, int init );
int Kalkreuter(su3_vector **eigVec, Real *eigVal, ks_eigen_param *eigen_param, int init);
int ks_eigensolve_ARPACK(su3_vector **eigVec, Real *eigVal, 
				  ks_eigen_param *eigen_param, int init );
int ks_eigensolve_Grid( su3_vector ** eigVec, Real * eigVal, ks_eigen_param * eigen_param, int init );
int ks_eigensolve_QUDA( su3_vector ** eigVec, Real * eigVal, ks_eigen_param * eigen_param, int init );
void Matrix_Vec_mult(su3_vector *src, su3_vector *res, ks_eigen_param *eigen_param, 
		     imp_ferm_links_t *fn );
void Precond_Matrix_Vec_mult(su3_vector *src, su3_vector *res, ks_eigen_param *eigen_param, 
			     imp_ferm_links_t *fn );
void cleanup_Matrix();
void measure_chirality(su3_vector *src, double *chirality, int parity);
void print_densities(su3_vector *src, char *tag, int y,int z,int t, 
		     int parity);
void reset_eigenvalues(su3_vector *eigVec[], Real *eigVal,
		       int Nvecs, int parity, imp_ferm_links_t *fn);
void perturb_eigpair(su3_vector *eigVec_new[], Real *eigVal_new,
		     su3_vector *eigVec_old[], Real *eigVal_old,
		     int Nvecs, int parity, imp_ferm_links_t *fn_new,
		     imp_ferm_links_t *fn_old);
void check_eigres(double *resid, su3_vector *eigVec[], Real *eigVal,
		  int Nvecs, int parity, imp_ferm_links_t *fn);
void construct_eigen_other_parity(su3_vector *eigVec[], Real eigVal[], 
				  ks_eigen_param *eigen_param, imp_ferm_links_t *fn);
#ifdef HAVE_QOP

/* fn_links_qop.c */

su3_matrix *get_fatlinks(fn_links_qop_t *fn);
su3_matrix *get_lnglinks(fn_links_qop_t *fn);
su3_matrix *get_fatbacklinks(fn_links_qop_t *fn);
su3_matrix *get_lngbacklinks(fn_links_qop_t *fn);

#else

/* fn_links_milc.c */

su3_matrix *get_fatlinks(fn_links_t *fn);
su3_matrix *get_lnglinks(fn_links_t *fn);
su3_matrix *get_fatbacklinks(fn_links_t *fn);
su3_matrix *get_lngbacklinks(fn_links_t *fn);

#endif

/* fn_links_milc.c only -- for QUDA */
int fresh_fn_links(imp_ferm_links_t *fn);
void refresh_fn_links(imp_ferm_links_t *fn);
void cancel_quda_notification(imp_ferm_links_t *fn);

/* fpi_2.c */
int fpi_2( /* Return value is number of C.G. iterations taken */
  Real *masses,   /* array of masses */
  int nmasses,      /* number of masses */
  Real tol,       /* tolerance for inverter check. */
  imp_ferm_links_t *fn       /* Storage for fat and Naik links */
  );

/* grsource.c */
void grsource(int parity);

/* grsource_imp.c */
void grsource_imp( field_offset dest, Real mass, int parity,
		   imp_ferm_links_t *fn );
void grsource_imp_plus( field_offset dest, field_offset rand, Real mass, int parity,
			imp_ferm_links_t *fn );
void grsource_imp_field( su3_vector *dest, Real mass, int parity,
			 imp_ferm_links_t *fn );
void grsource_plain( field_offset dest, int parity );
void grsource_plain_field( su3_vector *dest, int parity );
void z2rsource_imp( field_offset dest, Real mass, int parity,
		    imp_ferm_links_t *fn );
void z2rsource_plain( field_offset dest, int parity );
void z2rsource_plain_field( su3_vector *dest, int parity );
void checkmul_imp( field_offset src, Real mass,
		   imp_ferm_links_t *fn );

/* inc_eigcg.c */
#define ORTHO_EPS 1e-15
typedef struct {
  int m;             /* Number of vectors kept for the Lanczos part before restart */
  int Nvecs;         /* Number of eigenpairs computed per inversion */
  int Nvecs_curr;    /* Number of eigenpairs currently computed */
  int Nvecs_max;     /* Maximum number of eigenpairs computed in entire incremental eigCG */
  double_complex *H; /* H = -U^+ Dslash^2 U, U: projection onto smaller subspace */
} eigcg_params;
void calc_eigenpairs(Real *eigVal, su3_vector **eigVec, eigcg_params *eigcgp, int parity);
void calc_eigresid(int Nvecs, double *resid, double *norm, Real *eigVal,
		   su3_vector **eigVec, int parity, imp_ferm_links_t *fn);
int ks_eigCG_parity( su3_vector *src, su3_vector *dest, Real *eigVal, su3_vector **eigVec,
		     int m, int Nvecs, quark_invert_control *qic, Real mass, imp_ferm_links_t *fn);
int ks_inc_eigCG_parity( su3_vector *src, su3_vector *dest, Real *eigVal,
			 su3_vector **eigVec, eigcg_params *eigcgp, quark_invert_control *qic,
			 Real mass, imp_ferm_links_t *fn);

/* ks_baryon.c */
int baryon_type_index(char *label);
const char *baryon_type_label(int index);
void ks_baryon_nd(complex *prop[],
		  ks_prop_field *qp0, ks_prop_field *qp1, ks_prop_field *qp2,
		  int num_corr_b, int baryon_type_snk[], int phase[], Real fact[]);

/* ks_meson_mom.c */
void ks_meson_cont_mom(
  complex **prop,           /* prop[m][t] is where result is accumulated */
  su3_vector *src1,         /* quark propagator (to become antiquark) */
  su3_vector *src2,         /* quark propagator */
  int no_q_momenta,         /* no of unique mom/parity values (gt p) */
  int **q_momstore,         /* q_momstore[p] are the momentum components */
  char **q_parity,          /* q_parity[p] the parity of each mom component */
  int no_spin_taste_corr,   /* Number of unique spin-taste assignments (gt g) */
  int num_corr_mom[],       /* number of momentum/parity values (gt k)*/
  int **corr_table,         /* c = corr_table[k] correlator index */
  int p_index[],            /* p = p_index[c] is the momentum index */
  imp_ferm_links_t *fn_src1,/* Needed for some spin-taste ops */
  imp_ferm_links_t *fn_src2,/* Needed for some spin-taste ops */
  int spin_taste_snk[],     /* spin_taste_snk[c] gives the s/t assignment */
  int meson_phase[],        /* meson_phase[c] is the correlator phase */
  Real meson_factor[],      /* meson_factor[c] scales the correlator */
  int num_corr,             /* number of corrs - first index of prop */
  int corr_index[],         /* m = corr_index[c] is the correlator index */
  int r0[]                  /* spatial origin for defining FT phases */
		       );

/* mat_invert.c */
void ks_dirac_op( su3_vector *src, su3_vector *dst, Real mass, 
		  int parity, imp_ferm_links_t *fn);
void ks_dirac_adj_op( su3_vector *src, su3_vector *dst, Real mass,
		      int parity, imp_ferm_links_t *fn);
void ks_dirac_adj_op_inplace( su3_vector *srcdst, Real mass,
			      int parity, imp_ferm_links_t *fn);
int mat_invert_cg( field_offset src, field_offset dest, field_offset temp,
		   Real mass, int prec, imp_ferm_links_t *fn );
int mat_invert_field(su3_vector *src, su3_vector *dst, 
		     quark_invert_control *qic,
		     Real mass, imp_ferm_links_t *fn );
int mat_invert_cg_field(su3_vector *src, su3_vector *dst, 
			quark_invert_control *qic,
			Real mass, imp_ferm_links_t *fn );
int mat_invert_mg_field_gpu(su3_vector *src, su3_vector *dst, 
			    quark_invert_control *qic,
			    Real mass, imp_ferm_links_t *fn );
void mat_invert_mg_cleanup(void);
int mat_invert_uml(field_offset src, field_offset dest, field_offset temp,
		   Real mass, int prec, imp_ferm_links_t *fn );
int mat_invert_uml_field(su3_vector *src, su3_vector *dst, 
			 quark_invert_control *qic,
			 Real mass, imp_ferm_links_t *fn );
int mat_invert_block_cgz(su3_vector **src, su3_vector **dst, 
			 Real mass, int nsrc, quark_invert_control *qic,
			 imp_ferm_links_t *fn);
int mat_invert_block_cg(su3_vector **src, su3_vector **dst, 
			Real mass, int nsrc, quark_invert_control *qic,
			imp_ferm_links_t *fn);
int mat_invert_block_uml(su3_vector **src, su3_vector **dst, 
			 Real mass, int nsrc, quark_invert_control *qic,
			 imp_ferm_links_t *fn);
int mat_invert_block_mg(su3_vector **src, su3_vector **dst, 
			Real mass, int nsrc, quark_invert_control *qic,
			imp_ferm_links_t *fn);
int mat_invert_block(su3_vector **src, su3_vector **dst, 
		     Real mass, int nsrc, quark_invert_control *qic,
		     imp_ferm_links_t *fn);
void check_invert( field_offset src, field_offset dest, Real mass,
		   Real tol, imp_ferm_links_t *fn );
void check_invert_field( su3_vector *src, su3_vector *dest, Real mass,
			 Real tol, imp_ferm_links_t *fn, int parity);
void check_invert_field2( su3_vector *src, su3_vector *dest, Real mass,
			  Real tol, imp_ferm_links_t *fn, int parity);

/* multimass_inverter.c */
#define MAX_MMINV_NMASSES 32
#define MAX_MMINV_SOURCES 16
typedef struct {
  Real masses[MAX_MMINV_NMASSES];
  int nmasses;
  int n_sources;
  int r0[MAX_MMINV_SOURCES][4];
  Real tol;
  Real rsqprop;
} params_mminv;

int multimass_inverter( params_mminv *mminv, imp_ferm_links_t *fn);

/* read_eigen_param.c */
int read_ks_eigen_param(ks_eigen_param *eigen_param, int status, int prompt);

#if 0  /* obsolete */

/* nl_spectrum.c */
int nl_spectrum( Real vmass, field_offset tempvec1, field_offset tempvec2,
		 field_offset tempmat1, field_offset tempmat2,
		 imp_ferm_links_t *fn);

/* spectrum.c */
int spectrum(imp_ferm_links_t *fn);

/* spectrum2.c */
int spectrum2( Real vmass, field_offset temp1, field_offset temp2,
	       imp_ferm_links_t *fn);

/* spectrum_fzw.c */
int spectrum_fzw( Real vmass, field_offset temp1, field_offset temp2,
		  imp_ferm_links_t *fn );

/* spectrum_hybrids.c */
int spectrum_hybrids( Real mass, field_offset temp, Real tol,
		      imp_ferm_links_t *fn);

/* spectrum_mom.c */
int spectrum_mom( Real qmass, Real amass, field_offset temp, Real tol,
		  imp_ferm_links_t *fn);

/* spectrum_multimom.c */
int spectrum_multimom( Real dyn_mass, Real low_mass, Real mass_inc, 
		       int nmasses, Real tol, imp_ferm_links_t *fn);

/* spectrum_nd.c */
int spectrum_nd( Real mass1, Real mass2, Real tol, 
		 imp_ferm_links_t *fn );

/* spectrum_nlpi2.c */
int spectrum_nlpi2( Real qmass, Real amass, field_offset temp, Real tol,
		    imp_ferm_links_t *fn );
void mult_rho0( int fdir,  field_offset src, field_offset dest ) ;
void mult_rhos( int fdir,  field_offset src, field_offset dest ) ;

/* spectrum_singlets */
int spectrum_singlets( Real mass, Real tol, field_offset temp_offset,
		       imp_ferm_links_t *fn );
#endif

/* shift_field.c */

enum shift_dir {
  SHIFT_FORWARD,
  SHIFT_BACKWARD,
  SHIFT_SYMMETRIC
};

void 
shift_field(int dir, enum shift_dir fb, su3_vector *dest, const su3_vector *const src,
	    const su3_matrix *const links, int *refresh_links);

/* spin_taste_ops.c */
#include "../include/flavor_ops.h"
#ifdef NO_GAUGE_FIELD
void 
spin_taste_op_ape_fn( void *fn, int index, int r0[],
		      su3_vector *dest, const su3_vector *const src);
#else
void spin_taste_op_ape_fn(imp_ferm_links_t *fn, int index, int r0[],
			  su3_vector *dest, const su3_vector *const src);
#endif

#endif /* _IMP_FERM_LINKS_H */
