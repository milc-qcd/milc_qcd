#ifndef _FERMION_LINKS_H
#define _FERMION_LINKS_H
/************************ fermion_links.h **********************************
*									*
*  Definitions for fermion links
*  MIMD version 7 							*
*									*
*/

#include <quark_action.h>
#include "../include/imp_ferm_links.h"
#include "../include/ks_action_paths.h"

#ifdef HAVE_QOP
#include "../include/ks_action_coeffs_qop.h"
#include "../include/fermion_links_qop.h"
#else
#include "../include/fermion_links_milc.h"
#endif

/* Define the fermion links "class" for FN-type actions.

   Design considerations:

   The links are HISQ or fn and MILC-style or QOP.  The combinations give
   four choices for now.  QOP links can, in turn, be single or double
   precision. 

   We assume that only one of the four cases is in effect during
   compilation.  The choices are determined by compiler macros.  The
   API uses a generic fermion links type name and generic
   fermion-links procedure names (methods) that hide the choices.  To
   implement this in the C compilation, we provide four separate coded
   versions of the methods (inverters and fermion force routines).  In
   the Makefile we select one of them.  And we typedef the generic
   fermion links type to one of the four specific types.

   With QOP is is possible to have both precisions in use at the same
   time if the application requires it.

   To create the links, it is first necessary to create the link
   coefficients structure, set options, as desired, and then call
   either the HISQ or fn creation utility passing in the
   coefficients and the gauge field.  In both cases a pointer to the
   multipurpose object is returned.


*/

/*********************************************************************/
/* Options */

typedef struct {
  int want_du0;          /* Do we need to calculate the derivative wrto u0? */
  int want_deps;         /* Do we need to calculate the derivative wrto eps? */
  int want_aux;          /* Do we need to keep the auxiliary HISQ
			    links?  (They are needed for the HISQ
			    fermion force, but not if we just
			    computing propagators.) */
  int want_back;         /* Do we want backward links as well? */
} ferm_links_options_t;

/*********************************************************************/
/* The generic fermion links structure                               */

#if FERM_ACTION == HISQ

#ifdef HAVE_QOP
typedef qop_hisq_links_t ferm_links_generic_t;
#else
typedef milc_hisq_links_t ferm_links_generic_t;
#endif

#else

#ifdef HAVE_QOP
typedef qop_asqtad_links_t ferm_links_generic_t;
#else
typedef milc_fm_links_t ferm_links_generic_t;
#endif

#endif

/**********************************************************************/
/* This is the public fermion links type.  Its members are intended to
   be "private", however */
/**********************************************************************/

typedef struct {
  ferm_links_options_t options;
  ferm_links_generic_t *flg;
} fermion_links_t;



/********************************************************************/
/* Fermion force routines */
/********************************************************************/

#include "../include/macros.h"

/* fermion_force_asqtad_qop.c, fermion_force_asqtad.c */
void eo_fermion_force_oneterm( Real eps, Real weight, su3_vector *x_off,
			       int prec, fermion_links_t *fn);
void eo_fermion_force_oneterm_site( Real eps, Real weight, field_offset x_off,
				    int prec, fermion_links_t *fn);
void eo_fermion_force_twoterms( Real eps, Real weight1, Real weight2,
				su3_vector *x1_off, su3_vector *x2_off,
				int prec, fermion_links_t *fn);
void eo_fermion_force_twoterms_site( Real eps, Real weight1, Real weight2,
				     field_offset x1_off, field_offset x2_off,
				     int prec, fermion_links_t *fn);
void fermion_force_asqtad_block( Real eps, Real *residues, 
				 su3_vector **xxx, int nterms, int veclength, 
				 int prec, fermion_links_t *fl );
void
eo_fermion_force_twoterms_field_cpu( Real eps, Real weight1, Real weight2,
                                     half_wilson_vector *temp_x, int prec,
                                     fermion_links_t *fn,
                                     ks_action_paths *ap);
void
eo_fermion_force_twoterms_field_gpu( Real eps, Real weight1, Real weight2,
				     half_wilson_vector *temp_x, int prec,
				     fermion_links_t *fn,
				     ks_action_paths *ap);
#ifdef USE_FF_GPU
#define eo_fermion_force_twoterms_field eo_fermion_force_twoterms_field_gpu
#else
#define eo_fermion_force_twoterms_field eo_fermion_force_twoterms_field_cpu
#endif

void fermion_force_block( Real eps, Real *residues, 
			  su3_vector **xxx, int nterms, int veclength, 
			  int prec, fermion_links_t *fn);


/* fermion_force_asqtad_qop_F.c */

void eo_fermion_force_oneterm_F( Real eps, Real weight, su3_vector *x_off,
				 fermion_links_t *fl);
void eo_fermion_force_twoterms_F( Real eps, Real weight1, Real weight2, 
				  su3_vector *x1_off, su3_vector *x2_off,
				  fermion_links_t *fl);
void fermion_force_multi_F( Real eps, Real *residues, 
			    su3_vector **xxx, int nterms,
			    fermion_links_t *fl );
void fermion_force_block_F( Real eps, Real *residues, 
			    su3_vector **xxx, int nterms, 
			    int veclength, fermion_links_t *fl );

/* fermion_force_hisq_qop_F.c*/

void fermion_force_multi_hisq_F( Real eps, Real *residues, 
				 su3_vector **xxx, int n_orders_naik[],
				 fermion_links_t *fn);

/* fermion_force_asqtad_qop_D.c, fermion_force_hisq_qop_D.c */

void eo_fermion_force_oneterm_D( Real eps, Real weight, su3_vector *x_off,
				  fermion_links_t *fl );
void eo_fermion_force_twoterms_D( Real eps, Real weight1, Real weight2, 
				  su3_vector *x1_off, su3_vector *x2_off,
				  fermion_links_t *fl );
void fermion_force_multi_D( Real eps, Real *residues, 
			    su3_vector **xxx, int nterms,
			    fermion_links_t *fl );
void fermion_force_block_D( Real eps, Real *residues, 
			    su3_vector **xxx, int nterms, 
			    int veclength, fermion_links_t *fl );

/* fermion_force_hisq_qop_D.c*/

void fermion_force_multi_hisq_D( Real eps, Real *residues, 
				 su3_vector **xxx, int n_orders_naik[],
				 fermion_links_t *fn);
/* fermion_force_fn_multi.c */

enum ks_multiff_opt_t {ASVEC, FNMAT, FNMATREV};
  
const char *ks_multiff_opt_chr( void );

void eo_fermion_force_multi( Real eps, Real *residues, su3_vector **xxx, 
			     int nterms, int prec, fermion_links_t *fn);
void fermion_force_fn_multi( Real eps, Real *residues, su3_vector **multi_x, 
			     int nterms, int prec, fermion_links_t *fn);
void fermion_force_fn_multi_reverse( Real eps, Real *residues, 
				     su3_vector **multi_x, int nterms,
				     fermion_links_t *fn);
void show_hisq_force_opts( void );
void show_su3_mat_opts( void );
void show_hisq_links_opts( void );

/* Temporary, until we move completely to field-based gauge links */
/* fermion_links_from_site.c */

fermion_links_t *create_fermion_links_from_site(int prec, int n_naiks, double *eps_naik);
void restore_fermion_links_from_site(fermion_links_t *fl, int prec);

/* fermion_links_milc.c routines dealing with fermion_links_t */

fermion_links_t *create_fermion_links(int precision, int phases_in, su3_matrix *links);
void destroy_fermion_links(fermion_links_t *fl);
imp_ferm_links_t **get_fm_links(fermion_links_t *fl);
imp_ferm_links_t **get_fm_du0_links(fermion_links_t *fl);
ks_action_paths *get_action_paths(fermion_links_t *fl);
void invalidate_fermion_links(fermion_links_t *fl);
void restore_fermion_links(fermion_links_t *fl, int precision, int phases_in, su3_matrix *links);
int valid_fermion_links(fermion_links_t *fl, int precision);

/* fermion_links_eo_milc.c routines dealing with fermion_links_t */

ks_action_paths *get_action_paths_eo(fermion_links_t *fl);



/* fermion_links_hisq_milc.c routines dealing with fermion_links_t */
/* fermion_links_asqtad_qop.c routines dealing with fermion_links_t */

fermion_links_t *create_fermion_links_hisq(int precision, int n_naiks, 
			   double eps_naik[], int phases_in, su3_matrix *links);
void destroy_fermion_links_hisq(fermion_links_t *fl);
void invalidate_fermion_links(fermion_links_t *fl);
void restore_fermion_links_hisq(fermion_links_t *fl, int precision,
				int phases_in, su3_matrix *links);
imp_ferm_links_t **get_fm_links(fermion_links_t *fl);
imp_ferm_links_t **get_fm_du0_links(fermion_links_t *fl);
imp_ferm_links_t *get_fn_deps_links(fermion_links_t *fl);
ks_action_paths_hisq *get_action_paths_hisq(fermion_links_t *fl);
int get_n_naiks_hisq(fermion_links_t *fl);
double *get_eps_naik_hisq(fermion_links_t *fl);
int valid_fermion_links(fermion_links_t *fl, int precision);
char *get_action_parameter_string(fermion_links_t *fl);
hisq_auxiliary_t *get_hisq_auxiliary(fermion_links_t *fl);
#ifdef HAVE_QOP
QOP_asqtad_coeffs_t *get_action_coeffs(fermion_links_t *fl);
QOP_hisq_coeffs_t *get_action_coeffs_hisq(fermion_links_t *fl);
QOP_F3_FermionLinksHisq *get_F_hisq_links(fermion_links_t *fl);
QOP_D3_FermionLinksHisq *get_D_hisq_links(fermion_links_t *fl);
#endif

/* fermion_links.c */

void fermion_links_want_du0(int request);
void fermion_links_want_deps(int request);
void fermion_links_want_aux(int request);
void fermion_links_want_back(int request);

fermion_links_t *create_fermion_links_t(void);
void destroy_fermion_links_t(fermion_links_t *fl);

int fermion_links_get_n_naiks(fermion_links_t *fl);
double *fermion_links_get_eps_naik(fermion_links_t *fl);

/* f_meas.c */
void f_meas_imp_field( int npbp_reps, quark_invert_control *qic, Real mass, 
		       int naik_term_epsilon_index, fermion_links_t *fl);
void f_meas_imp_multi( int n_masses, int npbp_reps, quark_invert_control *qic, 
		       ks_param *my_ksp, fermion_links_t *fl);
/* DEPRECATED: */
void f_meas_imp( int npbp_reps, int prec, 
		 field_offset phi_off, field_offset xxx_off, Real mass,
		 int naik_term_epsilon_index, fermion_links_t *fl);

/* mu.c and mu_fast.c */
void M_derivatives(field_offset phi_off, field_offset xxx_off, 
		   field_offset xxx1_off, Real mass,
		   fermion_links_t *fl);
void Deriv_O6_field( int npbp_reps, quark_invert_control *qic, Real mass, 
		     fermion_links_t *fl, int naik_term_epsilon_index, Real eps);
void Deriv_O6_multi( int n_masses, int npbp_reps, quark_invert_control *qic, 
		     ks_param *ksp, fermion_links_t *fl);
/* DEPRECATED: */
void Deriv_O6(int npbp_reps, int prec, field_offset phi_off, field_offset xxx_off, 
	      field_offset xxx1_off, Real mass,
	      fermion_links_t *fl);

#endif /* _FERMION_LINKS_H */
