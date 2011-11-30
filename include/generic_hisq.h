OBSOLETE!!
#ifndef _GENERIC_HISQ_H
#define _GENERIC_HISQ_H

/* UNUSED.  - CD */

/************************ generic_hisq.h **********************************
*									*
*  Macros and declarations for generic_hisq routines                      *
*  This header is for codes that call generic_hisq routines               *
*  MIMD version 7 							*
*									*
*/

#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/generic_quark_types.h"
#include "../include/comdefs.h"
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

/* Structure defining the fermion action using paths or optimized
   coefficients */

#if FERM_ACTION ==  HISQ
typedef struct {
  Real *act_path_coeff;    /* For optimized Asqtad action */
  int num_q_paths;         /* For all actions */
  Q_path *q_paths;         /* For all actions */
} ks_component_paths;

typedef struct {
  ks_component_paths p1, p2;
  int umethod;
  Real mass;             /* The mass last used in the p2 coefficients */
} ks_action_paths;

typedef struct {
  // Flags: 1 if the corresponding link field is valid
  // (valid means it corresponds to the current links in the site structure)
  int valid_U_links, valid_V_links, valid_W_links, valid_Y_links, 
    valid_X_links;
  int valid_all_links; // should be 1 if ALL links are valid
  Real valid_Xfat_mass, valid_Xlong_mass;
  // phases...in = 1 if KS and antiperiodic BC signs are absorbed in the links
  int phases_in_U, phases_in_V, phases_in_W, phases_in_Y,
    phases_in_Xfat, phases_in_Xlong;
  su3_matrix *U_link[4]; // original gauge matrices, stored as four fields
  su3_matrix *V_link[4]; // first iteration of fattening
  su3_matrix *Y_unitlink[4]; // unitary projection of V_link, U(3)
  su3_matrix *W_unitlink[4]; // special unitary projection of Y_link, SU(3)
  su3_matrix *X_fatlink[4];
  su3_matrix *X_longlink[4];
} hisq_links_t;

#elif FERM_ACTION == ASQTAD  /* Non-HISQ actions */
typedef struct {
  Real *act_path_coeff;    /* For optimized Asqtad action */
  int num_q_paths;         /* For all actions */
  Q_path *q_paths;         /* For all actions */
} ks_action_paths;

#else
  BOMB THE COMPILIATION
#endif

/* Structure defining the precomputed links for the FN actions */

typedef struct {
  int valid;
  su3_matrix *fat;
  su3_matrix *lng;
  su3_matrix *fatback;
  su3_matrix *lngback;
  ks_action_paths *ap;  /* For EO actions */
#if FERM_ACTION == HISQ
  Real mass;    /* The mass last used in the coefficients */
  hisq_links_t hl;
#endif
#ifdef HAVE_QOP
  int valid_qop_F;
  int valid_qop_D;
  QOP_F3_FermionLinksAsqtad *qop_F_l;
  QOP_D3_FermionLinksAsqtad *qop_D_l;
#endif
} ferm_links_t;

int congrad( int niter, int nrestart, Real rsqmin, int parity, Real *rsq );
void copy_latvec(field_offset src, field_offset dest, int parity);
void dslash_site( field_offset src, field_offset dest, int parity );
void dslash_site_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );
void dslash_field( su3_vector *src, su3_vector *dest, int parity );
void dslash_field_special( su3_vector *src, su3_vector *dest,
    int parity, msg_tag **tag, int start );

void checkmul(void);
void phaseset(void);
void rephase( int flag );

void prefetch_vector( su3_vector * );
void prefetch_matrix( su3_matrix * );

int ks_congrad( field_offset src, field_offset dest, Real mass,
		int niter, int nrestart, Real rsqmin, int prec, 
		int parity, Real *rsq, ferm_links_t *fn );

int ks_congrad_field( su3_vector *src, su3_vector *dest, 
		      quark_invert_control *qic, Real mass,
		      ferm_links_t *fn);

int ks_congrad_site( field_offset src, field_offset dest, 
		     quark_invert_control *qic, Real mass,
		     ferm_links_t *fn);

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
    ferm_links_t *fn       /* Storage for fermion links */
    );

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[]);
void cleanup_dslash_temps(void);

void dslash_fn_site( field_offset src, field_offset dest, int parity,
		     ferm_links_t *fn);
void dslash_fn_site_special( field_offset src, field_offset dest,
			     int parity, msg_tag **tag, int start,
			     ferm_links_t *fn);
void ddslash_fn_du0_site( field_offset src, field_offset dest, int parity,
			  ferm_links_t *fn, ferm_links_t *fn_dmdu0);

void dslash_fn_field( su3_vector *src, su3_vector *dest, int parity,
		      ferm_links_t *fn);
void dslash_fn_field_special(su3_vector *src, su3_vector *dest,
			     int parity, msg_tag **tag, int start,
			     ferm_links_t *fn);
void ddslash_fn_du0_field( su3_vector *src, su3_vector *dest, int parity,
			   ferm_links_t *fn, ferm_links_t *fn_dmdu0);

void dslash_eo_site( field_offset src, field_offset dest, int parity,
		     ferm_links_t *fn);

/* The following three do not exist yet (3/05 -CD) */
void dslash_eo_site_special( field_offset src, field_offset dest,
			     int parity, msg_tag **tag, int start,
			     ferm_links_t *fn );
void dslash_eo_field( su3_vector *src, su3_vector *dest, int parity,
		      ferm_links_t *fn);
void dslash_eo_field_special( su3_vector *src, su3_vector *dest,
			      int parity, msg_tag **tag, int start,
			      ferm_links_t *fn);

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
enum ks_multicg_opt_t {OFFSET, HYBRID, FAKE, REVERSE, REVHYB};
const char *ks_multicg_opt_chr( void );

void ks_multicg_field(	        /* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* desired intermediate precision */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr, /* final residue squared */
    ferm_links_t *fn       /* Storage for fat and Naik links */
    );

int ks_multicg_offset_field(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    quark_invert_control *qic,  /* inversion parameters */
    ferm_links_t *fn      /* Storage for fat and Naik links */
    );

int ks_multicg_mass_field(	/* Return value is number of iterations taken */
    su3_vector *src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *masses,	/* the masses */
    int num_masses,	/* number of masses */
    quark_invert_control *qic,  /* inversion parameters */
    ferm_links_t *fn       /* Storage for fat and Naik links */
    );


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
	       int Restart, int iters, int parity,
	       ferm_links_t *fn );

/* f_meas.c */
void f_meas_imp( field_offset phi_off, field_offset xxx_off, Real mass,
		 ferm_links_t *fn, ferm_links_t *fn_dmdu0);

/* fpi_2.c */
int fpi_2( /* Return value is number of C.G. iterations taken */
  Real *masses,   /* array of masses */
  int nmasses,      /* number of masses */
  Real tol,       /* tolerance for inverter check. */
  ferm_links_t *fn       /* Storage for fat and Naik links */
  );

/* fermion_force_asqtad*.c */
void eo_fermion_force_oneterm( Real eps, Real weight, su3_vector *x_off,
			       int prec, ferm_links_t *fn,
			       ks_action_paths *ap);
void eo_fermion_force_oneterm_site( Real eps, Real weight, field_offset x_off,
				    int prec, ferm_links_t *fn,
				    ks_action_paths *ap);
void eo_fermion_force_twoterms( Real eps, Real weight1, Real weight2,
				su3_vector *x1_off, su3_vector *x2_off,
				int prec, ferm_links_t *fn, 
				ks_action_paths *ap );
void eo_fermion_force_twoterms_site( Real eps, Real weight1, Real weight2,
				     field_offset x1_off, field_offset x2_off,
				     int prec, ferm_links_t *fn, 
				ks_action_paths *ap );
void fermion_force_asqtad_block( Real eps, Real *residues, 
				 su3_vector **xxx, int nterms, int veclength, 
				 int prec, ferm_links_t *fn, 
				 ks_action_paths *ap );
void fermion_force_asqtad_multi( Real eps, Real *residues, 
				 su3_vector **xxx, int nterms, int prec,
				 ferm_links_t *fn, ks_action_paths *ap);

/* fermion_force_fn_multi.c */

enum ks_multiff_opt_t {ASVEC, FNMAT, FNMATREV};
  
const char *ks_multiff_opt_chr( void );

int eo_fermion_force_set_opt(char opt_string[]);
void eo_fermion_force_multi( Real eps, Real *residues, su3_vector **xxx, 
			     int nterms, int prec, ferm_links_t *fn,
			     ks_action_paths *ap );
void fermion_force_asqtad_block( Real eps, Real *residues, 
				 su3_vector **xxx, int nterms, int veclength, 
				 int prec, ferm_links_t *fn,
				 ks_action_paths *ap );
void fermion_force_asqtad_multi( Real eps, Real *residues, su3_vector **xxx, 
				 int nterms, int prec,
				 ferm_links_t *fn, ks_action_paths *ap );
void fermion_force_fn_multi( Real eps, Real *residues, su3_vector **multi_x, 
			     int nterms, int prec, ferm_links_t *fn,
			       ks_action_paths *ap );
void fermion_force_fn_multi_reverse( Real eps, Real *residues, 
				     su3_vector **multi_x, int nterms,
				     ferm_links_t *fn, ks_action_paths *ap);
void fermion_force_fn_multi_june05( Real eps, Real *residues, 
				    su3_vector **multi_x, int nterms,
				    ferm_links_t *fn, ks_action_paths *ap);
void fn_fermion_force_multi_hisq( Real eps, Real *residues, 
				  su3_vector **multi_x, int nterms,
				  ferm_links_t *fn, ks_action_paths *ap );


/* fermion_links_fn.c */
void init_ferm_links(ferm_links_t *fn);
void load_ferm_links(ferm_links_t *fn, ks_action_paths *ap);
void load_ferm_links_dmdu0(ferm_links_t *fn, ks_action_paths *ap);
void invalidate_all_ferm_links(ferm_links_t *fn);

/* fermion_links_fn_load_milc.c */
void load_longlinks(ferm_links_t *fn, ks_component_paths *p);
void load_fatlinks(ferm_links_t *fn, ks_component_paths *p);
void load_longbacklinks(ferm_links_t *fn);
void load_fatbacklinks(ferm_links_t *fn);
void free_fn_links(ferm_links_t *fn);
void free_fn_links_dmdu0(ferm_links_t *fn);
#if FERM_ACTION == HISQ
//void load_fatlinks_hisq( su3_matrix *Src, ks_component_paths *app, 
//			 su3_matrix *Dest );
void load_longlinks_hisq( su3_matrix *Src, ks_component_paths *app, 
			  su3_matrix *Dest );
#endif
void custom_rephase( su3_matrix **internal_links, int flag, int *status_now );

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
#define VECLENGTH 4
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
void grsource_imp( field_offset dest, Real mass, int parity,
		   ferm_links_t *fn );
void grsource_plain( field_offset dest, int parity );
void z2rsource_plain( field_offset dest, int parity );
void checkmul_imp( field_offset src, Real mass,
		   ferm_links_t *fn );

/* fermion_links_fn.c */
void init_ferm_links(ferm_links_t *fn);
void load_ferm_links(ferm_links_t *fn, ks_action_paths *ap);
void load_ferm_links_dmdu0(ferm_links_t *fn, ks_action_paths *ap);
void invalidate_all_ferm_links(ferm_links_t *fn);

/* jacobi.c */
#include "../include/jacobi.h"

/* mat_invert.c */
int mat_invert_cg( field_offset src, field_offset dest, field_offset temp,
		   Real mass, int prec, ferm_links_t *fn );
int mat_invert_uml(field_offset src, field_offset dest, field_offset temp,
		   Real mass, int prec, ferm_links_t *fn );
void check_invert( field_offset src, field_offset dest, Real mass,
		   Real tol, ferm_links_t *fn );

/* mu.c and mu_fast.c */
void M_derivatives(field_offset phi_off, field_offset xxx_off, 
		   field_offset xxx1_off, Real mass,
		   ferm_links_t *fn, ferm_links_t *fn_dmdu0);
void Deriv_O6(field_offset phi_off, field_offset xxx_off, 
	      field_offset xxx1_off, Real mass,
	      ferm_links_t *fn, ferm_links_t *fn_dmdu0);

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

int multimass_inverter( params_mminv *mminv, ferm_links_t *fn);

/* nl_spectrum.c */
int nl_spectrum( Real vmass, field_offset tempvec1, field_offset tempvec2,
		 field_offset tempmat1, field_offset tempmat2,
		 ferm_links_t *fn);

/* path_transport.c */
void link_gather_connection_hisq( su3_matrix *src, su3_matrix *dest,
				  su3_matrix *work, int dir );
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
void path_transport_connection_hisq( su3_matrix * src, su3_matrix **links, su3_matrix * dest,
    int parity, int *dir, int length );
void link_transport_connection_hisq( su3_matrix * src, su3_matrix **links, su3_matrix * dest,
    su3_matrix * work, int dir );

/* quark_stuff.c */
void make_path_table(ks_action_paths *ap, ks_action_paths *ap_dmdu0);
void make_path_tables(ks_action_paths *ap, ks_action_paths *ap_dmdu0,
		      Real mass);
int get_num_q_paths(void);
int get_num_q_paths_1(void);
int get_num_q_paths_2(void);
Q_path *get_q_paths(void);
Q_path *get_q_paths_1(void);
Q_path *get_q_paths_2(void);
int get_num_q_paths(void);
int get_num_q_paths_1(void);
int get_num_q_paths_2(void);
Real *get_quark_path_coeff_1(void);
Real *get_quark_path_coeff_2(void);

/* show_generic_ks_opts.c */
void show_generic_ks_opts( void );

/* spectrum.c */
int spectrum(ferm_links_t *fn);

/* spectrum2.c */
int spectrum2( Real vmass, field_offset temp1, field_offset temp2,
	       ferm_links_t *fn);

/* spectrum_fzw.c */
int spectrum_fzw( Real vmass, field_offset temp1, field_offset temp2,
		  ferm_links_t *fn );

/* spectrum_hybrids.c */
int spectrum_hybrids( Real mass, field_offset temp, Real tol,
		      ferm_links_t *fn);

/* spectrum_mom.c */
int spectrum_mom( Real qmass, Real amass, field_offset temp, Real tol,
		  ferm_links_t *fn);

/* spectrum_multimom.c */
int spectrum_multimom( Real dyn_mass, Real low_mass, Real mass_inc, 
		       int nmasses, Real tol, ferm_links_t *fn);

/* spectrum_nd.c */
int spectrum_nd( Real mass1, Real mass2, Real tol, 
		 ferm_links_t *fn );

/* spectrum_nlpi2.c */
int spectrum_nlpi2( Real qmass, Real amass, field_offset temp, Real tol,
		    ferm_links_t *fn );
void mult_rho0( int fdir,  field_offset src, field_offset dest ) ;
void mult_rhos( int fdir,  field_offset src, field_offset dest ) ;

/* spectrum_singlets */
int spectrum_singlets( Real mass, Real tol, field_offset temp_offset,
		       ferm_links_t *fn );
/* su3_mat_op.c */
void su3_unitarize( su3_matrix *a, su3_matrix *b );
void su3_spec_unitarize( su3_matrix *a, su3_matrix *b, complex *detA );
void su3_unit_der( su3_matrix *u, su3_tensor4 *dwdu, su3_tensor4 *dwdagdu );
void su3_unit_der_spec( su3_matrix *u, su3_matrix *w, su3_matrix *wdag, 
                        su3_tensor4 *dwdu, su3_tensor4 *dwdagdu,
                        su3_tensor4 *dvdu, su3_tensor4 *dvdagdu );
void su3_unitarize_analytic( su3_matrix *V, su3_matrix *W );
void su3_unit_der_analytic( su3_matrix *V, 
			    su3_tensor4 *dwdv, su3_tensor4 *dwdagdv );
void su3_unit_der_rational( su3_matrix *V, su3_tensor4 *dwdv, su3_tensor4 *dwdagdv );
void su3_unitarize_rational( su3_matrix *V, su3_matrix *W );
#endif /* _GENERIC_HISQ_H */
