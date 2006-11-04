#ifndef _GENERIC_KS_H
#define _GENERIC_KS_H
/************************ generic_ks.h **********************************
*									*
*  Macros and declarations for generic_ks routines                      *
*  This header is for codes that call generic_ks routines               *
*  MIMD version 7 							*
*									*
*/

#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/generic_quark_types.h"
#include "../include/comdefs.h"
#ifdef HAVE_QDP
#include <qdp.h>
#endif
#ifdef HAVE_QOP
#include <qop.h>
#endif

/* Structure specifying each rotation and reflection of each kind of
	path.  */
#define MAX_PATH_LENGTH 16
typedef struct {
  int dir[MAX_PATH_LENGTH];	/* directions in path */
  int length;		/* length of path */
  Real coeff;	        /* coefficient, including minus sign if backwards */
  Real forwback;	/* +1 if in forward Dslash, -1 if in backward */
} Q_path;

#ifdef HAVE_QDP
typedef struct {
  QLA_Real one_link     ; 
  QLA_Real naik         ;
  QLA_Real three_staple ;
  QLA_Real five_staple  ;
  QLA_Real seven_staple ;
  QLA_Real lepage       ;
} asqtad_path_coeff;
#endif

int congrad( int niter, int nrestart, Real rsqmin, int parity, Real *rsq );
void copy_latvec(field_offset src, field_offset dest, int parity);
void dslash_site( field_offset src, field_offset dest, int parity );
void dslash_site_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );

void checkmul();
void phaseset();
void rephase( int flag );

void prefetch_vector( su3_vector * );
void prefetch_matrix( su3_matrix * );

int ks_congrad_qop(int niter, int nrestart, Real rsqmin, 
		   Real *masses[], int nmass[], 
		   field_offset milc_srcs[], field_offset *milc_sols[],
		   int nsrc, Real* final_rsq_ptr, int milc_parity );

int ks_congrad_qop_site2field(int niter, int nrestart, Real rsqmin, 
			      Real *masses[], int nmass[], 
			      field_offset milc_srcs[], 
			      su3_vector **milc_sols[],
			      int nsrc, Real* final_rsq_ptr, int milc_parity );

int ks_congrad( field_offset src, field_offset dest, Real mass,
		int niter, int nrestart, Real rsqmin, int parity, Real *rsq );

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
    int parity,		/* parity to be worked on */
    Real  *final_rsq_ptr 	/* final residue squared */
    );

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[]);
void cleanup_dslash_temps();

void dslash_fn_site( field_offset src, field_offset dest, int parity );
void dslash_fn_site_special( field_offset src, field_offset dest,
			     int parity, msg_tag **tag, int start );
void ddslash_fn_du0_site( field_offset src, field_offset dest, int parity );

void dslash_fn_field( su3_vector *src, su3_vector *dest, int parity );
void dslash_fn_field_special(su3_vector *src, su3_vector *dest,
			     int parity, msg_tag **tag, int start );
void ddslash_fn_du0_field( su3_vector *src, su3_vector *dest, int parity );

void dslash_eo_site( field_offset src, field_offset dest, int parity );

/* The following three do not exist yet (3/05 -CD) */
void dslash_eo_site_special( field_offset src, field_offset dest,
			     int parity, msg_tag **tag, int start );
void dslash_eo_field( su3_vector *src, su3_vector *dest, int parity );
void dslash_eo_field_special( su3_vector *src, su3_vector *dest,
			      int parity, msg_tag **tag, int start );

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
			void *dmp),
    quark_invert_control *qic, /* inverter control */
    void *dmp                 /* Passthrough Dirac matrix parameters */
    );

/* in ks_multicg.c */
#define OFFSET  0
#define HYBRID  1
#define FAKE    2
#define REVERSE 3
#define REVHYB  4

int ks_multicg_set_opt(char opt_string[]);

int ks_multicg(	        /* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
    );

int ks_multicg_offset(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
    );

int ks_multicg_mass(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *masses,	/* the masses */
    int num_masses,	/* number of masses */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
    );

int ks_multicg_hybrid(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
    );

int ks_multicg_reverse(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *masses,	/* the masses */
    int num_masses,	/* number of masses */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
    );

int ks_multicg_fake(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
    );

int ks_multicg_revhyb(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
    );

/* d_congrad5_fn_qop.c */
#ifdef HAVE_QOP
void create_qop_asqtad_fermion_links( QOP_FermionLinksAsqtad** qop_links );
#endif
void initialize_congrad( void );
void finalize_congrad( void );
int ks_congrad_qop_site2site(int niter, int nrestart, Real rsqmin, 
			     Real *masses[], int nmass[], 
			     field_offset milc_srcs[], 
			     field_offset *milc_sols[],
			     int nsrc, Real* final_rsq_ptr, int milc_parity );

/* dslash_fn_qop_milc.c */
void cleanup_gathers_qop_milc(msg_tag *tags1[], msg_tag *tags2[]);
void cleanup_dslash_qop_milc_temps();
void dslash_fn_qop_milc( su3_matrix *fatlinks, su3_matrix *longlinks,
			 su3_vector *src, su3_vector *dest, int parity );
void dslash_fn_qop_milc_field_special(su3_matrix *fatlinks, 
				      su3_matrix *longlinks,
				      su3_vector *src, su3_vector *dest,
				      int parity, msg_tag **tag, int start );


/* d_congrad_opt.c */

void clear_latvec(field_offset v,int parity);

void scalar_mult_latvec(field_offset src, Real scalar,
			field_offset dest, int parity);
void scalar_mult_add_latvec(field_offset src1, field_offset src2,
			    Real scalar, field_offset dest, int parity);
void scalar2_mult_add_su3_vector(su3_vector *a, Real s1, su3_vector *b, 
				 Real s2, su3_vector *c);
void scalar_mult_add_lathwvec_proj_su3mat(su3_matrix *mom, 
					  half_wilson_vector *back, 
					  half_wilson_vector *forw, 
					  Real coeff[2]);
void scalar2_mult_add_latvec(field_offset src1,Real scalar1,
			     field_offset src2,Real scalar2,
			     field_offset dest,int parity);
/* eigen_stuff.c */
int Rayleigh_min(su3_vector *vec,su3_vector **eigVec,Real Tolerance, 
		 Real RelTol,int Nvecs,int MaxIter,int Restart,int parity);
int Kalkreuter(su3_vector **eigVec, double *eigVal, Real Tolerance, 
	       Real RelTol, int Nvecs, int MaxIter, 
	       int Restart, int iters, int parity) ;

/* f_meas.c */
void f_meas_imp( field_offset phi_off, field_offset xxx_off, Real mass );

/* fpi_2.c */
int fpi_2( /* Return value is number of C.G. iterations taken */
  Real *masses,   /* array of masses */
  int nmasses,      /* number of masses */
  Real tol        /* tolerance for inverter check. */
  );

/* fermion_force_asqtad*.c */
void eo_fermion_force_oneterm( Real eps, Real weight, field_offset x_off );
void eo_fermion_force_twoterms( Real eps, Real weight1, Real weight2,
				field_offset x1_off, field_offset x2_off );

/* fermion_force_fn_multi.c */

int eo_fermion_force_set_opt(char opt_string[]);
void eo_fermion_force_multi( Real eps, Real *residues, su3_vector **xxx, 
			     int nterms );
void eo_fermion_force_asqtad_block( Real eps, Real *residues, su3_vector **xxx, int nterms, int veclength );
void eo_fermion_force_asqtad_multi( Real eps, Real *residues, su3_vector **xxx, int nterms );
void fn_fermion_force_multi( Real eps, Real *residues, su3_vector **multi_x, int nterms );
void fn_fermion_force_multi_reverse( Real eps, Real *residues, su3_vector **multi_x, int nterms );
void fn_fermion_force_multi_june05( Real eps, Real *residues, su3_vector **multi_x, int nterms );

/* fermion_links_fn.c */
void load_fn_links();
void load_fn_links_dmdu0();

/* fermion_links_helpers.c */
void load_longbacklinks(su3_matrix **t_lbl, su3_matrix *t_ll);
void load_fatbacklinks(su3_matrix **t_fbl, su3_matrix *t_fl);
void free_fn_links();
void free_fn_links_dmdu0();

/* ff_opt.c */
void mult_adj_su3_fieldlink_lathwvec( su3_matrix *link,
				      half_wilson_vector **src_pt, 
				      half_wilson_vector *dest);
void mult_su3_sitelink_lathwvec( int dir, 
				 half_wilson_vector **src_pt, 
				 half_wilson_vector *dest);
void scalar_mult_add_lathwvec_proj(anti_hermitmat *mom, 
				   half_wilson_vector *back, 
				   half_wilson_vector *forw, Real coeff[2]);
void scalar_mult_add_lathwvec(half_wilson_vector *dest, 
			      half_wilson_vector *src, Real s[2]);
void mult_su3_fieldlink_lathwvec( su3_matrix *link,
				  half_wilson_vector **src_pt, 
				  half_wilson_vector *dest);
#ifndef VECLENGTH
#define VECLENGTH 1
#endif

typedef struct { su3_vector v[VECLENGTH]; } veclist;
void mult_adj_su3_fieldlink_latveclist( su3_matrix *link,
               veclist **src_pt, veclist *dest, int listlength );
void mult_su3_sitelink_latveclist( int dir, 
	 veclist **src_pt, veclist *dest, int listlength );
void scalar_mult_add_latveclist_proj(anti_hermitmat *mom,
           veclist *back, veclist *forw, Real *coeff, int listlength );
void scalar_mult_add_latveclist( veclist *dest,
            veclist *src, Real *s, int listlength );

/* flavor_ops.c */
void sym_shift(int dir, field_offset src,field_offset dest) ;
void zeta_shift(int n, int *d, field_offset src, field_offset dest ) ;
void eta_shift(int n, int *d, field_offset src, field_offset dest ) ;


void mult_flavor_vector(int mu, field_offset src, field_offset dest ) ;
void mult_flavor_tensor(int mu, int nu, field_offset src, field_offset dest ) ;
void mult_flavor_pseudovector(int mu, field_offset src, field_offset dest ) ;
void mult_flavor_pseudoscalar(field_offset src, field_offset dest ) ;

void mult_spin_vector(int mu, field_offset src, field_offset dest ) ;
void mult_spin_tensor(int mu, int nu, field_offset src, field_offset dest ) ;
void mult_spin_pseudovector(int mu, field_offset src, field_offset dest ) ;
void mult_spin_pseudoscalar(field_offset src, field_offset dest ) ;

/* grsource.c */
void grsource(int parity);

/* grsource_imp.c */
void grsource_imp( field_offset dest, Real mass, int parity);
void grsource_plain( field_offset dest, int parity );

/* jacobi.c */
#include "../include/jacobi.h"

/* mat_invert.c */
int mat_invert_cg( field_offset src, field_offset dest, field_offset temp,
		   Real mass );
int mat_invert_uml(field_offset src, field_offset dest, field_offset temp,
		   Real mass );
void check_invert( field_offset src, field_offset dest, Real mass,
		   Real tol);

/* mu.c and mu_fast.c */
void M_derivatives(field_offset phi_off, field_offset xxx_off, field_offset xxx1_off, Real mass);
void Deriv_O6(field_offset phi_off, field_offset xxx_off, field_offset xxx1_off, Real mass );

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

int multimass_inverter( params_mminv *mminv);

/* nl_spectrum.c */
int nl_spectrum( Real vmass, field_offset tempvec1, field_offset tempvec2,
		 field_offset tempmat1, field_offset tempmat2);

/* path_transport.c */
void path_transport_site( field_offset src, field_offset dest, int parity,
			  int *dir, int length );
void path_transport_field( su3_vector * src, su3_vector * dest, int parity, 
			   int *dir, int length );
void path_transport_hwv_site( field_offset src, field_offset dest, int parity,
			      int *dir, int length );
void path_transport_hwv_field( half_wilson_vector *src, 
			       half_wilson_vector * dest, int parity,
			       int *dir, int length );
void path_transport_connection( su3_matrix * src, su3_matrix * dest, int parity, int *dir, int length );
void link_transport_connection( su3_matrix * src, su3_matrix * dest, su3_matrix * work, int dir );

/* quark_stuff.c */
void make_path_table();
int get_num_q_paths();
Q_path *get_q_paths();
Q_path *get_q_paths_dmdu0();
Real *get_quark_path_coeff();
Real *get_quark_path_coeff_dmdu0();

/* spectrum.c */
int spectrum();

/* spectrum2.c */
int spectrum2( Real vmass, field_offset temp1, field_offset temp2 );

/* spectrum_fzw.c */
int spectrum_fzw( Real vmass, field_offset temp1, field_offset temp2 );

/* spectrum_hybrids.c */
int spectrum_hybrids( Real mass, field_offset temp, Real tol );

/* spectrum_mom.c */
int spectrum_mom( Real qmass, Real amass, field_offset temp, Real tol);

/* spectrum_multimom.c */
int spectrum_multimom( Real dyn_mass, Real low_mass, Real mass_inc, int nmasses, Real tol);

/* spectrum_nd.c */
int spectrum_nd( Real mass1, Real mass2, Real tol );

/* spectrum_nlpi2.c */
int spectrum_nlpi2( Real qmass, Real amass, field_offset temp, Real tol);
void mult_rho0( int fdir,  field_offset src, field_offset dest ) ;
void mult_rhos( int fdir,  field_offset src, field_offset dest ) ;

/* spectrum_singlets */
int spectrum_singlets( Real mass, Real tol, field_offset temp_offset );

#ifdef HAVE_QDP

void dslash_qdp_fn(QDP_ColorVector *src, QDP_ColorVector *dest,
		   QDP_Subset parity);
void dslash_qdp_fn_special(QDP_ColorVector *src, QDP_ColorVector *dest,
			   QDP_Subset parity, QDP_ColorVector *temp[]);
void dslash_qdp_fn_special2(QDP_ColorVector *src, QDP_ColorVector *dest,
			    QDP_Subset parity, QDP_ColorVector *temp[]);
void 
fn_fermion_force_qdp( QDP_ColorMatrix *force[], QDP_ColorMatrix *gf[], 
		      asqtad_path_coeff *coeffs, QDP_ColorVector *in_pt[], 
		      Real eps[], int nsrc );
int ks_congrad_qdp(QDP_ColorVector *src, QDP_ColorVector *dest, QLA_Real mass,
		   int niter, int nrestart, QLA_Real rsqmin, QDP_Subset parity,
		   QLA_Real *final_rsq_ptr);
int ks_multicg_qdp(QDP_ColorVector *src, QDP_ColorVector **dest,
		   QLA_Real *masses, int num_masses, int niter,
		   QLA_Real rsqmin, QDP_Subset parity,
		   QLA_Real *final_rsq_ptr);

#endif /* HAVE_QDP */

#ifdef HAVE_QOP
void load_qop_asqtad_coeffs(QOP_asqtad_coeffs_t *c, Real weight,
			    Real *act_path_coeff);
#endif

#endif /* _GENERIC_KS_H */
