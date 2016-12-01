#ifndef _GENERIC_H
#define _GENERIC_H
/************************ generic.h *************************************
*									*
*  Macros and declarations for miscellaneous generic routines           *
*  This header is for codes that call generic routines                  *
*  MIMD version 7 							*
*									*
*/

/* Other generic directory declarations are elsewhere:

   For com_*.c, see comdefs.h
   For io_lat4.c io_ansi.c, io_nonansi.c, io_piofs.c, io_romio.c see io_lat.h
   For io_prop_w.c, see io_wprop.h
*/

#include <stdio.h>
#include "../include/int32type.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/random.h"
#include "../include/file_types.h"
#include "../include/io_lat.h"
#include "../include/generic_quark_types.h"

/* ape_smear.c */

void ape_smear_field(
  su3_matrix *src,       /* Gauge field input unsmeared */
  su3_matrix *dest,      /* Gauge field output smeared */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  );

void ape_smear_field_dir(
  su3_matrix *src,        /* su3_matrix[4] type 
			     input unsmeared links */
  int dir1,               /* link direction to smear */
  su3_matrix *dest,       /* su3_matrix[4] type smeared links */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
			 );

void ape_smear_dir(
  field_offset src,       /* field offset for su3_matrix[4] type 
			     input unsmeared links */
  int dir1,               /* link direction to smear */
  field_offset dest,      /* field offset for su3_matrix type 
			     pointing to a specific direction 
			     output smeared links */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  );

void ape_smear(
  field_offset src,       /* field offset for su3_matrix type 
			     input unsmeared links */
  field_offset dest,      /* field offset for su3_matrix type 
			     output smeared links */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  );

su3_matrix *ape_smear_3D(Real staple_weight, int iters);
su3_matrix *ape_smear_4D(Real staple_weight, int iters);
void destroy_ape_links_3D(su3_matrix *ape_links);
void destroy_ape_links_4D(su3_matrix *ape_links);

/* ax_gauge.c */
void ax_gauge(void);

/* bsd_sum.c */
int32type bsd_sum (char *data,int32type total_bytes);

/* check_unitarity.c */
Real check_unitarity( void );

/* d_linktrsum */
void d_linktrsum(double_complex *linktrsum);

/* d_plaq?.c */
void d_plaquette(double *ss_plaq,double *st_plaq);

/* discretize_wf.c */
void fnal_wavefunction(complex *wf, int stride,
		       int x0, int y0, int z0, int t0, 
		       Real a, char wf_file[]);

/* field_strength.c */
void make_field_strength(
  field_offset link_src,       /* field offset for su3_matrix[4] type 
				  for the source link matrices */
  field_offset field_dest      /* field offset for su3_matrix[6] type
				  for the resulting field strength */
  );

/* field_translation.c */
void shift_gauge(int rshift[]);
void shift_complex(complex *src, int rshift[]);
void shift_su3_vector(su3_vector *src, int rshift[]);
void shift_wilson_vector(wilson_vector *src, int rshift[]);

/* field_utilities.c */
double start_timing(void);
void print_timing(double dtime, char *str);

Real* create_r_field(void);
void clear_r_field(Real *r);
void copy_r_field(Real *dest, Real *src);
void destroy_r_field(Real *r);

complex* create_c_field(void);
void clear_c_field(complex *c);
void copy_c_field(complex *dest, complex *src);
void destroy_c_field(complex *c);

su3_matrix *create_m_field(void);
void clear_m_field(su3_matrix *m);
void copy_m_field(su3_matrix *dest, su3_matrix *src);
void destroy_m_field(su3_matrix *m);

spin_wilson_vector *create_swv_field(void);
void clear_swv_field(spin_wilson_vector *swv);
void copy_swv_field(spin_wilson_vector *dest, spin_wilson_vector *src);
void destroy_swv_field(spin_wilson_vector *swv);

wilson_vector *create_wv_field(void);
void clear_wv_field(wilson_vector *wv);
void copy_wv_field(wilson_vector *dst, wilson_vector *src);
void destroy_wv_field(wilson_vector *wv);

su3_vector *create_v_field(void);
su3_vector *create_v_field_from_site_member(field_offset sv);
void clear_v_field(su3_vector *v);
void copy_v_field(su3_vector *dst, su3_vector *src);
void copy_site_member_from_v_field(field_offset sv, su3_vector *v);
void add_v_fields(su3_vector *vsum, su3_vector *v1, su3_vector *v2);
void destroy_v_field(su3_vector *v);

/* array versions of the above */

Real* create_r_array_field(int n);
void clear_r_array_field(Real *r, int n);
void copy_r_array_field(Real *dest, Real *src, int n);
void destroy_r_array_field(Real *r, int n);

complex* create_c_array_field(int n);
void clear_c_array_field(complex *c, int n);
void copy_c_array_field(complex *dest, complex *src, int n);
void destroy_c_array_field(complex *c, int n);

su3_matrix *create_m_array_field(int n);
void clear_m_array_field(su3_matrix *m, int n);
void copy_m_array_field(su3_matrix *dest, su3_matrix *src, int n);
void destroy_m_array_field(su3_matrix *m, int n);

spin_wilson_vector *create_swv_array_field(int n);
void clear_swv_array_field(spin_wilson_vector *swv, int n);
void copy_swv_array_field(spin_wilson_vector *dest, spin_wilson_vector *src, int n);
void destroy_swv_array_field(spin_wilson_vector *swv, int n);

wilson_vector *create_wv_array_field(int n);
void clear_wv_array_field(wilson_vector *wv, int n);
void copy_wv_array_field(wilson_vector *dst, wilson_vector *src, int n);
void destroy_wv_array_field(wilson_vector *wv, int n);

su3_vector *create_v_array_field(int n);
void clear_v_array_field(su3_vector *v, int n);
void copy_v_array_field(su3_vector *dst, su3_vector *src, int n);
void destroy_v_array_field(su3_vector *v, int n);

ks_prop_field *create_ksp_field(int nc);
ks_prop_field *create_ksp_field_copy(ks_prop_field *k);
void clear_ksp_field(ks_prop_field *ksp);
void copy_ksp_field(ks_prop_field *kspcopy, ks_prop_field *ksp);
void free_ksp_field(ks_prop_field *ksp);
void destroy_ksp_field(ks_prop_field *ksp);

wilson_prop_field * create_wp_field(int nc);
wilson_prop_field * create_wp_field_copy(wilson_prop_field * w);
void clear_wp_field(wilson_prop_field * wp);
void copy_wp_field(wilson_prop_field * wpcopy, wilson_prop_field * wp);
void scalar_mult_add_ksprop_field(ks_prop_field *a, ks_prop_field *b, 
				  Real s, ks_prop_field *c);
void scalar_mult_add_wprop_field(wilson_prop_field *a, wilson_prop_field *b, 
				 Real s, wilson_prop_field *c);
wilson_prop_field *transpose_wp_field(wilson_prop_field * wp);
void free_wp_field(wilson_prop_field * wp);
void rebuild_wp_field(wilson_prop_field * wp);
void destroy_wp_field(wilson_prop_field * wp);

void gauge_field_copy(field_offset src,field_offset dest);
void gauge_field_copy_field_to_site(su3_matrix **src, field_offset dest);
su3_matrix **gauge_field_copy_site_to_field(field_offset src);
void destroy_gauge_field(su3_matrix **f);

void extract_c_from_v(complex *c, su3_vector *v, int color);
void insert_v_from_c(su3_vector *v, complex *c, int color);

void extract_c_from_wv(complex *c, wilson_vector *wv, 
		       int spin, int color);
void insert_wv_from_c(wilson_vector *wv, complex *c, 
		      int spin, int color);
void insert_wv_from_v(wilson_vector *wv, su3_vector *v, int spin);

void extract_v_from_wv(su3_vector *v, wilson_vector *wv, int spin);

void copy_v_from_ksp(su3_vector *v, ks_prop_field *ksp, int color);
void insert_ksp_from_v(ks_prop_field *ksp, su3_vector *v, int color);

void copy_wp_from_wv(wilson_prop_field * wp, wilson_vector *wv, 
		     int color, int spin);
void copy_wv_from_v(wilson_vector *wv, su3_vector *v, int spin);
void copy_wv_from_wp(wilson_vector *wv, wilson_prop_field * wp, 
		     int color, int spin);

void copy_wv_from_swv(wilson_vector *wv, spin_wilson_vector *swv, int spin);
spin_wilson_vector *extract_swv_from_wp(wilson_prop_field * wp, int color);

void copy_wv_from_wprop(wilson_vector *wv, wilson_propagator *wprop, 
			int color, int spin);

/* file_types_milc_usqcd.c */
int ks_prop_milc_to_usqcd(int milc_type);
int ks_prop_usqcd_to_milc(int usqcd_type);
int w_prop_milc_to_usqcd(int milc_type);
int w_prop_usqcd_to_milc(int usqcd_type);


/* gaugefix.c and gaugefix2.c */
void gaugefix( int gauge_dir, Real relax_boost, int max_gauge_iter,
	       Real gauge_fix_tol );

void gaugefix_combo(int gauge_dir,Real relax_boost,int max_gauge_iter,
		    Real gauge_fix_tol, int nvector, 
		    field_offset vector_offset[], int vector_parity[],
		    int nantiherm, field_offset antiherm_offset[], 
		    int antiherm_parity[] );

/* gauge_force_imp.c and gauge_force_symzk1_qop.c */
/* gauge_force_imp.c and gauge_force_symzk1_qop.c */
void imp_gauge_force_cpu( Real eps, field_offset mom_off );
void imp_gauge_force_gpu( Real eps, field_offset mom_off );

#ifdef USE_GF_GPU
#define imp_gauge_force imp_gauge_force_gpu
#else
#define imp_gauge_force imp_gauge_force_cpu
#endif

/* gauge_stuff.c */
double imp_gauge_action(void);
void g_measure(void);
void make_loop_table(void);
void dsdu_qhb_subl(int dir, int subl);
int get_max_length(void);
int get_nloop(void);
int get_nreps(void);
int *get_loop_length(void);
int *get_loop_num(void);
int ***get_loop_table(void);
Real **get_loop_coeff(void);

/* gauge_utilities.c */
su3_matrix * create_G(void);
su3_matrix * create_G_from_site(void);
void copy_G(su3_matrix *dst, su3_matrix *src);
void destroy_G(su3_matrix *t_links);
su3_matrix * create_G2(void);
su3_matrix * create_G2_from_site(void);
void copy_G2(su3_matrix *dst, su3_matrix *src);
void destroy_G2(su3_matrix *t_links);

/* general_staple.c */
void 
compute_gen_staple_field(su3_matrix *staple, int mu, int nu, 
			 su3_matrix *link, int stride,
			 su3_matrix *fatlink, Real coef,
			 su3_matrix *links);

/* glueball_op.c */
void make_glueball_ops(void);
void measure_glueball_ops(void);

/* hvy_pot.c */
void hvy_pot( su3_matrix *links, int max_t, int max_x );

/* io_detect.c */
int get_file_type(char *filename);
int io_detect(char *filename, file_table ft[], int ntypes);
int io_detect_fm(char *filename);
int io_detect_ks_usqcd(char *filename);
int io_detect_w_usqcd(char *filename);

/* io_helpers.c */
gauge_file *save_lattice( int flag, char *filename, char *stringLFN );
gauge_file *reload_lattice( int flag, char *filename);
int ask_corr_file( FILE *fp, int prompt, int *flag, char* filename);
int ask_starting_lattice( FILE *fp, int prompt, int *flag, char *filename );
int ask_ending_lattice( FILE *fp, int prompt, int *flag, char *filename );
int ask_ildg_LFN(FILE *fp, int prompt, int flag, char *stringLFN);
void coldlat(void);
void funnylat(void);
int get_check_tag(FILE *fp, char *tag, char *myname);
int get_f( FILE *fp, int prompt, char *variable_name_string, Real *value );
int get_i( FILE *fp, int prompt, char *variable_name_string, int *value );
char *get_next_tag(FILE *fp, char *tag, char *myname);
int get_vi( FILE *fp, int prompt, char *variable_name_string, 
	    int *value, int nvalues );
int get_vf( FILE *fp, int prompt, char *variable_name_string, 
	    Real *value, int nvalues );
int get_s( FILE *fp, int prompt, char *variable_name_string, char *value );
int get_sn( FILE *fp, int prompt, char *variable_name_string, char *value );
int get_vs( FILE *fp, int prompt, char *tag, char *value[], int nvalues );
int get_prompt( FILE *fp, int *value );

/* io_source_cmplx_fm.c */
void r_source_cmplx_fm_to_site(char *filename, field_offset dest_site,
			       int x0, int y0, int z0, int t0);
void r_source_cmplx_fm_to_field(char *filename, complex *dest_field, int stride,
				int x0, int y0, int z0, int t0);

/* layout_*.c */
int io_node(const int node);
void setup_layout( void );
int node_number(int x,int y,int z,int t);
int node_index(int x,int y,int z,int t);
size_t num_sites(int node);
const int *get_logical_dimensions(void);
const int *get_logical_coordinate(void);
void get_coords(int coords[], int node, int index);

/* make_lattice.c */
void make_lattice(void);
void free_lattice(void);

/* nersc_cksum.c */
u_int32type nersc_cksum( void );

/* make_global_fields.c */
void make_global_fields(void);

/* momentum_twist.c */
void boundary_twist_field(Real bdry_phase[4], int r0[4], int sign, su3_matrix *links);
void boundary_twist_site(Real bdry_phase[4], int r0[4], int sign);
void momentum_twist_site(Real bdry_phase[4], int sign);
void momentum_twist_links(Real bdry_phase[4], int sign, su3_matrix *links);
void rephase_v_field(su3_vector *v, Real bdry_phase[4], int r0[4], int sign);
void rephase_wv_field(wilson_vector *wv, Real bdry_phase[4], int r0[4], int sign);

/* path_product.c */
void path_product( const int *dir, const int length, su3_matrix *tempmat1);
void path_product_field( const int *dir, const int length, 
			 su3_matrix *tempmat1, su3_matrix *links);
void path_product_fields( su3_matrix *Src, const int *dir, 
			  const int length, su3_matrix *tempmat1);
void path_prod_subl(const int *dir, const int length, const int subl,
		    su3_matrix *tempmat1);

/* phases.c */
int decode_phase(char *label);
void mult_c_by_phase(complex *a, complex *b, int ph);


/* plaquette4.c */
void plaquette(Real *ss_plaq,Real *st_plaq);

/* ploop?.c */
complex ploop( void );

/* ploop_staple.c */
complex ploop_staple(Real alpha_fuzz);

/* project_su3_hit.c */
void project_su3(
   su3_matrix *w,         /* input initial guess. output resulting
                             SU(3) matrix */
   su3_matrix *q,         /* starting 3 x 3 complex matrix */
   int Nhit,              /* number of SU(2) hits. 0 for no projection */
   Real tol              /* tolerance for SU(3) projection.
			     If nonzero, treat Nhit as a maximum
			     number of hits.  If zero, treat Nhit
			     as a prescribed number of hits. */ 
   );

/* quark_source.c */
void init_qs(quark_source *qs);
void alloc_cached_c_source(quark_source *qs);
void alloc_cached_v_source(quark_source *qs);
void alloc_cached_wv_source(quark_source *qs);
int ask_starting_source( FILE *fp, int prompt, int *flag, char *filename );
complex *get_cached_c_source(quark_source *qs);
su3_vector *get_cached_v_source(quark_source *qs);
wilson_vector *get_cached_wv_source(quark_source *qs);
void clear_qs(quark_source *qs);
int convert_ksource_to_color(int ksource);
int convert_ksource_to_spin(int ksource);
char *decode_mask(int mask);
int encode_mask(int *mask, char c_mask[]);
void even_and_odd_wall(complex *c, int t0);
void gaussian_source(complex *src, Real r0, 
		     int x0, int y0, int z0, int t0);
int is_complex_source(int source_type);
int is_vector_source(int source_type);
int is_dirac_source(int source_type);

// int v_base_source(su3_vector *src, quark_source *qs);
//int wv_base_source(wilson_vector *src, quark_source *qs);
int v_source_field(su3_vector *src, quark_source *qs);
int v_source_site(field_offset src, quark_source *qs);
int wv_source_field(wilson_vector *src, quark_source *qs);
int wv_source_site(field_offset src, quark_source *qs);
int get_v_quark_source(FILE *fp, int prompt, quark_source *qs);
int get_wv_quark_source(FILE *fp, int prompt, quark_source *qs);
void subset_mask_c(complex *src, int subset, int t0);
void subset_mask_v(su3_vector *src, int subset, int t0);
void subset_mask_wv(wilson_vector *src, int subset, int t0);
void print_source_info(FILE *fp, char prefix[], quark_source *qs);

/* quark_source_io.c */
int choose_usqcd_ks_file_type(int source_type);
int choose_usqcd_w_file_type(int source_type);
void r_source_open(quark_source *qs);
void r_source_close(quark_source *qs);
#ifdef HAVE_QIO
int r_source_cmplx_scidac(QIO_Reader *infile, complex *src,  
			  int x0, int y0, int z0, int t0);
QIO_Reader *r_source_cmplx_scidac_open(char source_file[]);
#endif
int r_source_vector(quark_source *qs);
int r_source_dirac(quark_source *qs);
int w_source_open_ks(quark_source *qs, char *fileinfo);
int w_source_open_dirac(quark_source *qs, char *fileinfo);
void w_source_close(quark_source *qs);
int w_source_ks(su3_vector *src, quark_source *qs);
int w_source_dirac(wilson_vector *src, quark_source *qs);
int w_source_dirac_site(field_offset src, quark_source *qs);
void print_output_quark_source_choices(void);
int parse_output_quark_source_choices(int *flag, int *save_type, 
				      char *descrp, char* savebuf);
int ask_output_quark_source_file( FILE *fp, int prompt, 
				  int *flag, int *source_type,
				  int *t0, char *descrp, char *filename);

/* quark_source_sink_op.c */
void init_qss_op(quark_source_sink_op *qss_op);
void broadcast_quark_source_sink_op_recursive(quark_source_sink_op **qss_op);
quark_source_sink_op *create_qss_op(void);
void destroy_qss_op(quark_source_sink_op *qss_op);
quark_source_sink_op *copy_qss_op_list(quark_source_sink_op *src_qss_op);
void insert_qss_op(quark_source *qs, quark_source_sink_op *qss_op);
void insert_qss_eps_naik_index(int index, quark_source_sink_op *qss_op);
void v_field_op(su3_vector *src, quark_source_sink_op *qss_op, 
		int subset, int t0);
void wv_field_op(wilson_vector *src, quark_source_sink_op *qss_op, 
		 int subset, int t0);
void ksp_sink_op(quark_source_sink_op *qss_op, ks_prop_field *ksp );
void wp_sink_op(quark_source_sink_op *qss_op, wilson_prop_field *wp );
int get_wv_field_op(FILE *fp, int prompt, quark_source_sink_op *qss_op);
int get_v_field_op(FILE *fp, int prompt, quark_source_sink_op *qss_op);
int get_qss_eps_naik(Real *eps_naik, quark_source_sink_op *qss_op);
void print_field_op_info(FILE *fp, char prefix[], 
			 quark_source_sink_op *qss_op);
void print_field_op_info_list(FILE *fp, char prefix[], 
			      quark_source_sink_op *qss_op[], int n);
void set_qss_op_offset(quark_source_sink_op *qss_op, int r0[]);

/* rand_gauge.c */
void rand_gauge(field_offset G);

/* ranmom.c */
void ranmom( void );

/* remap standard I/O */
int remap_stdio_from_args(int argc, char *argv[]);

/* ranstuff.c */
void initialize_prn(double_prn *prn_pt, int seed, int index);
Real myrand(double_prn *prn_pt);
void initialize_site_prn_from_seed(int iseed);

/* restrict_fourier.c */
void setup_restrict_fourier( int *key, int *slice);
void restrict_fourier_site( 
     field_offset src,	 /* src is field to be transformed */
     int size,		 /* Size of field in bytes.  The field must
			    consist of size/sizeof(complex) consecutive
			    complex numbers.  For example, an su3_vector
			    is 3 complex numbers. */
     int isign);	 /* 1 for x -> k, -1 for k -> x */

void restrict_fourier_field( 
     complex *src,       /* src is field to be transformed */
     int size,		 /* Size of field in bytes.  The field must
			    consist of size/sizeof(complex) consecutive
			    complex numbers.  For example, an su3_vector
			    is 3 complex numbers. */
     int isign);	 /* 1 for x -> k, -1 for k -> x */
void cleanup_restrict_fourier(void);

/* reunitarize2.c */
void reunitarize( void );
int reunit_su3(su3_matrix *c);

/* show_generic_opts.c */
void show_generic_opts( void );

/* show_scidac_opts.c */
void show_scidac_opts( void );

/* Do Morninstar-Peardon stout smearing to construct unitary W from
   smeared link V and unsmeared link U */
void stout_smear(su3_matrix *W, su3_matrix *V, su3_matrix *U);


#ifdef QCDOC
void *qcdoc_alloc(size_t nbytes);
void qfree(void *);
#endif

#endif	/* _GENERIC_H */
