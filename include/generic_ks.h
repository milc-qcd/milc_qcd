#ifndef _GENERIC_KS_H
#define _GENERIC_KS_H
/************************ generic_ks.h **********************************
*									*
*  Macros and declarations for generic_ks routines                      *
*  This header is for codes that call generic_ks routines               *
*  MIMD version 7 							*
*									*
*/

#include "../include/generic.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/comdefs.h"
#include "../include/info.h"
#include "../include/params_rhmc.h"
#include "../include/imp_ferm_links.h"
#include <stdint.h>

void copy_latvec(field_offset src, field_offset dest, int parity);
void dslash_site( field_offset src, field_offset dest, int parity );
void dslash_site_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );
void dslash_field( su3_vector *src, su3_vector *dest, int parity );
void dslash_field_special( su3_vector *src, su3_vector *dest,
    int parity, msg_tag **tag, int start );

void checkmul(void);

void prefetch_vector( su3_vector * );
void prefetch_matrix( su3_matrix * );

/* charge_utilities.c */
int fill_charge(double charge_table[], int *n, double next_charge);
int index_charge(double charge_table[], int n, double find_charge);
void start_charge(double charge_table[], int *n);

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

/* gauge_action_imp_ks.c */
double imp_gauge_action_ks(void);

/* gauge_force_imp_ks.c */
void imp_gauge_force_ks( Real eps, field_offset mom_off );

/* gauge_measure_ks.c */
void g_measure_ks(void);

/* gauss_smear_ks.c */
void gauss_smear_v_field(su3_vector *src, su3_matrix *t_links,
			 Real width, int iters, int t0);
void laplacian_v_field(su3_vector *src, su3_matrix *t_links, int t0);
void gauss_smear_ks_prop_field(ks_prop_field *src, su3_matrix *t_links,
			       Real width, int iters, int t0);
/* gauss_smear_ks_cpu.c */
void gauss_smear_reuse_2link_cpu( int flag );
void gauss_smear_delete_2link_cpu();
void gauss_smear_v_field_cpu_twolink(su3_vector *src, su3_matrix *t_links,
				     Real width, int iters, int t0);
void klein_gord_field(su3_vector *psi, su3_vector *chi, 
		      su3_matrix *t_links, Real msq, int t0);
void klein_gord_field_twolink(su3_vector *psi, su3_vector *chi, 
			      su3_matrix *t_links, Real msq, int t0);
void gauss_smear_v_field_cpu(su3_vector *src, su3_matrix *t_links,
			     Real width, int iters, int t0);
/* gauss_smear_ks_QUDA.c */
void gauss_smear_v_field_QUDA(su3_vector *src, su3_matrix *t_links,
                              Real width, int iters, int t0);
void gauss_smear_reuse_2link_QUDA( int flag );
void gauss_smear_delete_2link_QUDA();

/* grsource_rhmc.c */
void grsource_imp_rhmc( field_offset dest, params_ratfunc *rf,
			int parity, su3_vector **multi_x, su3_vector *sumvec,
			Real my_rsqmin, int my_niter, int my_prec,
			imp_ferm_links_t *fn, int naik_term_epsilon_index,
			Real naik_term_epsilon);

/* This version takes dest to be a field */
void grsource_imp_rhmc_field( su3_vector *dest, params_ratfunc *rf,
			      int parity, su3_vector **multi_x, su3_vector *sumvec,
			      Real my_rsqmin, int my_niter, int my_prec,
			      imp_ferm_links_t *fn, int naik_term_epsilon_index,
			      Real naik_term_epsilon);

/* ks_ratinv.c */
int ks_rateval(
    su3_vector *dest,   /* answer vector */
    field_offset src,   /* source vector (for a_0 term) */
    su3_vector **psim,  /* solution vectors  from multiroot CG */
    Real *residues,     /* the residues */
    int order,          /* order of approximation */
    int parity          /* parity to be worked on */
    );

int ks_rateval_field(	
    su3_vector *dest,	/* answer vector */
    su3_vector *src,	/* source vector (for a_0 term) */
    su3_vector **psim,	/* solution vectors  from multiroot CG */
    Real *residues,	/* the residues */
    int order,		/* order of approximation */
    int parity		/* parity to be worked on */
			);

int ks_ratinv(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *roots,	/* the roots */
    Real *residues,	/* the residues (if not zero, there are used to scale the shifted target residual) */
    int order,		/* order of rational function approx */
    int my_niter,	/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* desired intermediate precicion */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    imp_ferm_links_t *fn,     /* Fermion links */
    int naik_term_epsilon_index, /* Index of naik term common to this set */
    Real naik_term_epsilon /* Epsilon common to this set */
		);

int ks_ratinv_field(	/* Return value is number of iterations taken */
    su3_vector * src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *roots,	/* the roots */
    Real *residues,     /* residues */
    int order,		/* order of rational function approx */
    int my_niter,	/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* desired intermediate precicion */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    imp_ferm_links_t *fn_const, /* Fermion links */
    int naik_term_epsilon_index, /* Index of naik term common to this set */
    Real naik_term_epsilon /* Epsilon common to this set */
			);

/* naik_epsilon_utilities.c */
int fill_eps_naik(double eps_naik_table[], int *n, double next_eps_naik);
int index_eps_naik(double eps_naik_table[], int n, double find_eps_naik);
void start_eps_naik(double eps_naik_table[], int *n);

/* path_transport.c */
void link_gather_connection_hisq( su3_matrix *src,
				  su3_matrix *dest, su3_matrix *work, int dir );
void path_transport_site( field_offset src, field_offset dest, int parity,
			  int *dir, int length );
void path_transport_field( su3_vector * src, su3_vector * dest, int parity,
			   int *dir, int length );
void path_transport_hwv_site( field_offset src, field_offset dest, int parity,
			      int *dir, int length );
void path_transport_hwv_field( half_wilson_vector *src,
			       half_wilson_vector * dest, int parity,
			       int *dir, int length );
void path_transport_connection( su3_matrix * src, su3_matrix * dest,
    int parity, int *dir, int length );
void link_transport_connection( su3_matrix * src, su3_matrix * dest,
    su3_matrix * work, int dir );
void path_transport_connection_hisq( su3_matrix * src, su3_matrix **links,
    su3_matrix * dest, int parity, int *dir, int length );
void link_transport_connection_hisq( su3_matrix * src, su3_matrix *links,
    su3_matrix * dest, su3_matrix * work, int dir );

/* load_rhmc_params.c */

params_rhmc *load_rhmc_params(char filename[], int n_pseudo);

/* ploop3_ks.c */
complex ploop_ks(void);

/* rephase.c */
void apply_apbc( su3_matrix *links, int r0t );
void phaseset(void);
void rephase( int flag );
void rephase_field_offset( su3_matrix *internal_links, int flag,
			   int *status_now, int r0[] );
void rephase_offset( int flag, int r0[] );

/* reunitarize_ks.c */
void reunitarize_ks(void);

/* show_generic_ks_opts.c */
void show_generic_ks_opts( void );
void show_generic_ks_md_opts( void );

/* su3_mat_op.c */

void unity_su3mat( su3_matrix *dest );
Real su3_norm_frob( su3_matrix *a );
void su3_inverse( su3_matrix *a, su3_matrix *b );
void u3_root_inv( su3_matrix *a, su3_matrix *x, su3_matrix *y);
void u3_unitarize( su3_matrix *a, su3_matrix *b );
void su3_spec_unitarize( su3_matrix *a, su3_matrix *b, complex *detA );
void su3_spec_unitarize_index( su3_matrix *a, su3_matrix *b, complex *detA,
           int index_site, int index_dir );
void u3_unit_der( su3_matrix *u, su3_tensor4 *dwdu, su3_tensor4 *dwdagdu );
void mult_su3_t4m( su3_tensor4 *a, su3_matrix *b, su3_tensor4 *c );
void mult_su3_mt4( su3_matrix *a, su3_tensor4 *b, su3_tensor4 *c );
void add_su3_t4( su3_tensor4 *a, su3_tensor4 *b, su3_tensor4 *c );
void sub_su3_t4( su3_tensor4 *a, su3_tensor4 *b, su3_tensor4 *c );
Real su3_t4_norm_frob( su3_tensor4 *a );
void dumptensor4( su3_tensor4 *a );
void su3_spec_unit_der( su3_matrix *u, su3_tensor4 *dwdu, su3_tensor4 *dwdagdu );
void su3_unit_der_spec( su3_matrix *u, su3_matrix *w, su3_matrix *wdag,
                        su3_tensor4 *dwdu, su3_tensor4 *dwdagdu,
                        su3_tensor4 *dvdu, su3_tensor4 *dvdagdu );
void su3_unit_der_reim( su3_matrix *u,
                        su3_tensor4 *dwdure, su3_tensor4 *dwduim,
                        su3_tensor4 *dwdagdure, su3_tensor4 *dwdagduim );
void su3_unit_der_reim_join(
        su3_tensor4 *dwdure, su3_tensor4 *dwduim,
        su3_tensor4 *dwdagdure, su3_tensor4 *dwdagduim,
        su3_tensor4 *dwdu, su3_tensor4 *dwdagdu );
void u3_unitarize_rational( su3_matrix *V, su3_matrix *W );
void u3_unit_der_rational( su3_matrix *V,
                            su3_tensor4 *dwdv, su3_tensor4 *dwdagdu );
void su3_der_detWY( su3_matrix *y, su3_tensor4 *dwdy, su3_tensor4 *dwdagdy );
void u3_unitarize_analytic( info_t *info, su3_matrix *V, su3_matrix *W);
void u3_unitarize_analytic_index( su3_matrix *V, su3_matrix *W, int index_site, int index_dir );
void u3_unit_der_analytic( info_t *info, su3_matrix *V, su3_tensor4 *dwdv,
			   su3_tensor4 *dwdagdv);
void su3t4_copy( su3_tensor4 *a, su3_tensor4 *b );
int svd3x3(double A[3][3][2], double *sigma, double U[3][3][2],
	   double V[3][3][2], size_t *nflops);

#ifdef GB_BARYON
     /* gb_baryon_mmap.c */
     void set_mmap_src_number(int num);
     void fetch_ksp_from_cache(ks_prop_field **dest,int qknum,int scIdx,int skIdx);
     void fetch_su3v_from_cache(su3_vector **dest,int qknum,int skIdx,int disp,int color);
     void fetch_src_from_cache(su3_vector **dest,int srcnum,int color);
     void fetch_int_from_cache(su3_vector **dest,int stIdx);
     void fetch_spec_from_cache(su3_vector **dest);
     void msync_ksp_from_cache(int qknum,int scIdx,int skIdx);
     void msync_su3v_from_cache(int qknum,int skIdx,int disp,int color);
     void msync_src_from_cache(int srcnum,int color);
     void toss_ksp_from_cache(ks_prop_field **ksp);
     void populate_next_mmap_src(su3_vector *src,int color);
     void create_qk_oct_cache(ks_prop_field **qko,int qknum,int r0[],su3_matrix *links);
     void create_sink_oct_cache(su3_vector **sko,int qknum,int r0[],su3_matrix *links);
     void destroy_qk_oct_cache(int qknum);
     void unmap_qk_oct_cache(int qknum);
     mmap_cache* get_qk_cache_pointer(int qknum);
     void assign_qk_cache_pointer(mmap_cache* qcache, int qknum);
     void destroy_sink_oct_cache(int qknum);
     void copy_qk_oct_cache(int qkcpy, int qkone);
     void create_gb_qk_cache(int numqk);
     void create_gb_src_cache(int numsrc);
     void create_gb_int_cache(int numcur,int nummom);
     void create_gb_spec_cache(int nummom);
     void destroy_gb_qk_cache();
     void destroy_gb_src_cache();
     void destroy_gb_int_cache();
     void destroy_gb_spec_cache();
     void int_map_add_buffer(int stIdx);
     void set_spec_fill(int i,short val);
     short get_spec_fill(int i);
     int int_map_buffer_num(int stIdx);
     int int_map_indices(int ci,int ki,int t,int cidx,int momi);
     int spec_map_indices(int iqkn,int ci,int cj,int ki,int kj,int kk,int disp,int moms,int sc0);

     /* gb_baryon_3pt.c */
     void apply_par_xport_3pt(ks_prop_field *dest, ks_prop_field *src,
                              int n, int dir[], int r0[], short doBW, su3_matrix *links);
     /* gb_baryon_src.c */
     void apply_par_xport_src_v(su3_vector *dest, su3_vector *src,
                                quark_source_sink_op *qss_op, su3_matrix *links);
#endif /* GB_BARYON */

#endif /* _GENERIC_KS_H */
