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

/* Structure specifying each rotation and reflection of each kind of
	path.  */
#define MAX_PATH_LENGTH 16
typedef struct {
  int dir[MAX_PATH_LENGTH];	/* directions in path */
  int length;		/* length of path */
  Real coeff;	/* coefficient, including minus sign if backwards */
#ifdef DM_DU0
  Real coeff2;	/* coefficient for d(Dslash)/d(u0) */
#endif
  Real forwback;	/* +1 if in forward Dslash, -1 if in backward */
} Q_path;

int congrad( int niter, Real rsqmin, int parity, Real *rsq );
void copy_latvec(field_offset src, field_offset dest, int parity);
void dslash_site( field_offset src, field_offset dest, int parity );
void dslash_site_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );
void clear_latvec(field_offset v,int parity);

void scalar_mult_latvec(field_offset src, Real scalar,
			field_offset dest, int parity);
void scalar_mult_add_latvec(field_offset src1, field_offset src2,
			    Real scalar, field_offset dest, int parity);
void scalar2_mult_add_su3_vector(su3_vector *a, Real s1, su3_vector *b, 
				 Real s2, su3_vector *c);

void scalar2_mult_add_latvec(field_offset src1,Real scalar1,
			     field_offset src2,Real scalar2,
			     field_offset dest,int parity);
void checkmul();
void phaseset();
void rephase( int flag );

void prefetch_vector( su3_vector * );
void prefetch_matrix( su3_matrix * );

int ks_congrad( field_offset src, field_offset dest, Real mass,
     int niter, Real rsqmin, int parity, Real *rsq );

int ks_congrad_two_src(	/* Return value is number of iterations taken */
    field_offset src1,    /* source vector (type su3_vector) */
    field_offset src2,
    field_offset dest1,	/* solution vectors */
    field_offset dest2,
    Real mass1,
    Real mass2,
    int niter,		/* maximal number of CG interations */
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

int ks_multicg(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *masses,	/* the masses */
    int num_masses,	/* number of masses */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
    );

/* d_congrad5_fn_qop.c */
void initialize_congrad( void );
void finalize_congrad( void );

void congrad_fn_allocate_qop_fields( Real** qop_fat_links, 
				     Real** qop_long_links, 
				     Real** qop_src, Real** qop_sol );

void congrad_fn_map_milc_to_qop_raw( field_offset milc_src, 
				     field_offset milc_sol,
				     Real* qop_fat_links, 
				     Real* qop_long_links,
				     Real* qop_src, 
				     Real* qop_sol, int milc_parity );


void congrad_fn_map_qop_raw_to_milc( Real* qop_sol, field_offset milc_sol, 
				     int milc_parity );


#ifdef HAVE_QOP
#include <qop.h>
void congrad_fn_set_qop_invert_arg( QOP_invert_arg* qop_invert_arg, Real mass, 
			 int max_iterations, Real min_resid_sq, 
			 int milc_parity );

int ks_congrad_qop( Real* qop_source, Real* qop_solution,
		    Real* qop_fat_links, Real* qop_long_links,
		    QOP_invert_arg* qop_invert_arg, Real* final_rsq_ptr );

#endif

#ifdef HAVE_QDP
/* d_congrad5_fn_qopqdp.c */
void set_M_from_strided_parity_temp(QDP_ColorMatrix *dest, su3_matrix *src,
				    int stride, int parity);
void set_V_from_parity_field(QDP_ColorVector *dest, field_offset src,
			     int parity);
void set_parity_field_from_V(field_offset dest, QDP_ColorVector *src,
			     int parity);
int ks_congrad_qopqdp( QDP_ColorVector *qop_src, QDP_ColorVector *qop_sol,
		       QDP_ColorMatrix *qop_fat_links[4], 
		       QDP_ColorMatrix *qop_long_links[4],
		       QOP_invert_arg* qop_invert_arg, Real* final_rsq_ptr );
#endif

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

/* quark_stuff.c */
void make_path_table();
int get_num_q_paths();
Q_path *get_q_paths();
Real *get_quark_path_coeff();
void eo_fermion_force( Real eps, int nflavors, field_offset x_off );
void eo_fermion_force_3f( Real eps, int nflav1, field_offset x1_off,
	int nflav2, field_offset x2_off  );
void load_longlinks();
void load_fatlinks();
void free_longlinks();
void free_fatlinks();
void path_transport( field_offset src, field_offset dest, int parity,
    int *dir, int length );
void path_transport_hwv( field_offset src, field_offset dest, int parity,
    int *dir, int length );

/* spectrum.c */
int spectrum();

/* spectrum2.c */
int spectrum2( Real vmass, field_offset temp1, field_offset temp2 );

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
int ks_congrad_qdp(QDP_ColorVector *src, QDP_ColorVector *dest, QLA_Real mass,
		   int niter, QLA_Real rsqmin, QDP_Subset parity,
		   QLA_Real *final_rsq_ptr);
int ks_multicg_qdp(QDP_ColorVector *src, QDP_ColorVector **dest,
		   QLA_Real *masses, int num_masses, int niter,
		   QLA_Real rsqmin, QDP_Subset parity,
		   QLA_Real *final_rsq_ptr);

#endif /* HAVE_QDP */
#endif /* _GENERIC_KS_H */
