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

/* gauss_smear_ks.c */
void gauss_smear_v_field(su3_vector *src, su3_matrix *t_links,
			 Real width, int iters, int t0);
void gauss_smear_ks_prop_field(ks_prop_field *src, su3_matrix *t_links,
			       Real width, int iters, int t0);

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

/* rephase.c */
void phaseset(void);
void rephase( int flag );
void rephase_field_offset( su3_matrix *internal_links, int flag, 
			   int *status_now, int r0[] );
void rephase_offset( int flag, int r0[] );

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

#endif /* _GENERIC_KS_H */
