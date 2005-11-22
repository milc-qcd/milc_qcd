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

#include "../include/int32type.h"
#include "../include/complex.h"
#include "../include/macros.h"
#include "../include/random.h"
#include "../include/file_types.h"

/* ape_smear.c */
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

/* ax_gauge.c */
void ax_gauge();

/* bsd_sum.c */
int32type bsd_sum (char *data,int32type total_bytes);

/* check_unitarity.c */
Real check_unitarity( void );

/* d_plaq?.c */
void d_plaquette(double *ss_plaq,double *st_plaq);

/* field_strength.c */
void make_field_strength(
  field_offset link_src,       /* field offset for su3_matrix[4] type 
				  for the source link matrices */
  field_offset field_dest      /* field offset for su3_matrix[6] type
				  for the resulting field strength */
  );

/* gaugefix.c and gaugefix2.c */
void gaugefix(int gauge_dir,Real relax_boost,int max_gauge_iter,
	      Real gauge_fix_tol, field_offset diffmat, field_offset sumvec,
	      int nvector, field_offset vector_offset[], int vector_parity[],
	      int nantiherm, field_offset antiherm_offset[], 
	      int antiherm_parity[] );

/* gauge_stuff.c */
double imp_gauge_action();
void imp_gauge_force( Real eps, field_offset mom_off );
void make_loop_table();
void dsdu_qhb_subl(int dir, int subl);

/* glueball_op.c */
void make_glueball_ops();
void measure_glueball_ops();

/* hvy_pot.c */
void hvy_pot( field_offset links );

/* io_detect.c */
int io_detect(char *filename, file_type ft[], int ntypes);

/* layout_*.c */
void setup_layout( void );
int node_number(int x,int y,int z,int t);
int node_index(int x,int y,int z,int t);
size_t num_sites(int node);
int *get_logical_dimensions();
int *get_logical_coordinate();

/* make_lattice.c */
void make_lattice();
void free_lattice();

/* make_global_fields.c */
void make_global_fields();

/* path_product.c */
void path_product( const int *dir, const int length, su3_matrix *tempmat1);
void path_prod_subl(const int *dir, const int length, const int subl,
		    su3_matrix *tempmat1);

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

/* rand_gauge.c */
void rand_gauge(field_offset G);

/* ranmom.c */
void ranmom( void );

/* remap standard I/O */
int remap_stdio_from_args(int argc, char *argv[]);

/* ranstuff.c */
void initialize_prn(double_prn *prn_pt, int seed, int index);
Real myrand(double_prn *prn_pt);

/* restrict_fourier.c */
void setup_restrict_fourier( int *key, int *slice);
void restrict_fourier( 
     field_offset src,	 /* src is field to be transformed */
     field_offset space, /* space is working space, same size as src */
     field_offset space2,/* space2 is working space, same size as src */
                         /* space2 is needed only for non power of 2 */
     int size,		 /* Size of field in bytes.  The field must
			    consist of size/sizeof(complex) consecutive
			    complex numbers.  For example, an su3_vector
			    is 3 complex numbers. */
     int isign);	 /* 1 for x -> k, -1 for k -> x */

/* reunitarize2.c */
void reunitarize( void );
int reunit_su3(su3_matrix *c);

#ifdef QCDOC
void *qcdoc_alloc(size_t nbytes);
void qfree(void *);
#endif

/* map_milc_to_qop.c */
#ifdef HAVE_QOP
#include <qop.h>
QOP_status_t initialize_qop();
su3_matrix **create_raw_G_from_site_links();
void destroy_raw_G(su3_matrix *rawlinks[]);
su3_matrix **create_raw_F_from_site_mom();
void unload_raw_F_to_site_mom(su3_matrix *rawforce[]);
void destroy_raw_F(su3_matrix *rawforce[]);
su3_vector *create_raw_V_from_site(field_offset x);
void destroy_raw_V(su3_vector *rawsu3vec);
#endif

#ifdef HAVE_QDP
#include <qdp.h>

void set_V_from_field(QDP_ColorVector *dest, field_offset src);
void set_H_from_field(QDP_HalfFermion *dest, field_offset src);
void set_D_from_field(QDP_DiracFermion *dest, field_offset src);
void set_M_from_field(QDP_ColorMatrix *dest, field_offset src);

void set_field_from_V(field_offset dest, QDP_ColorVector *src);
void set_field_from_H(field_offset dest, QDP_HalfFermion *src);
void set_field_from_D(field_offset dest, QDP_DiracFermion *src);
void set_field_from_M(field_offset dest, QDP_ColorMatrix *src);

void set_V_from_temp(QDP_ColorVector *dest, su3_vector *src);
void set_H_from_temp(QDP_HalfFermion *dest, half_wilson_vector *src);
void set_M_from_temp(QDP_ColorMatrix *dest, su3_matrix *src);

void set_temp_from_V(su3_vector *dest, QDP_ColorVector *src);
void set_temp_from_M(su3_matrix *dest, QDP_ColorMatrix *src);

void set4_V_from_temp(QDP_ColorVector *dest[], su3_vector *src);
void set4_M_from_temp(QDP_ColorMatrix *dest[], su3_matrix *src);
#endif

#endif	/* _GENERIC_H */
