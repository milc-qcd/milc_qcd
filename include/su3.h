#ifndef _SU3_H
#define _SU3_H

#include "../include/complex.h"
#include "../include/random.h"

/******************************  su3.h **********************************
*									*
*  Defines and subroutine declarations for SU3 simulation		*
*  MIMD version 7 							*
*									*
*/
/* SU(3) */

typedef struct { fcomplex e[3][3]; } fsu3_matrix;
typedef struct { fcomplex c[3]; } fsu3_vector;
typedef struct { 
  fcomplex m01,m02,m12; 
  float m00im,m11im,m22im; 
  float space; } fanti_hermitmat;

typedef struct { dcomplex e[3][3]; } dsu3_matrix;
typedef struct { dcomplex c[3]; } dsu3_vector;
typedef struct { 
  dcomplex m01,m02,m12; 
  double m00im,m11im,m22im; 
  double space; } danti_hermitmat;

#if (PRECISION==1)

#define su3_matrix      fsu3_matrix
#define su3_vector      fsu3_vector
#define anti_hermitmat  fanti_hermitmat

#else

#define su3_matrix      dsu3_matrix
#define su3_vector      dsu3_vector
#define anti_hermitmat  danti_hermitmat

#endif

/* For KS spectroscopy */
typedef struct {
  su3_vector ** v;  /* For the nc vector fields */
  int nc;           /* number of vectors */
} ks_prop_field;

/* Used in HISQ codes */
/* Rank 4 tensor for storing derivatives */
typedef struct { fcomplex t4[3][3][3][3]; } fsu3_tensor4;
typedef struct { dcomplex t4[3][3][3][3]; } dsu3_tensor4;

#if (PRECISION==1)
#define su3_tensor4 fsu3_tensor4
#else
#define su3_tensor4 dsu3_tensor4
#endif

/* SU(2) */
typedef struct { complex e[2][2]; } su2_matrix;

/* Wilson vectors */

/* e.g. 					     */
/* wilson_propagator prop;                           */
/* prop.c[ci].d[si].d[sf].c[cf]                      */
/* ----------------------->    complex               */
/* ----------------->          su3_vector            */
/* ----------->                wilson_vector         */
/* ----->                      spin_wilson_vector    */
/* e.g. 					     */
/* wilson_matrix matr;                               */
/* matr.d[si].c[ci].d[sf].c[cf]                      */
/* ----------------------->    complex               */
/* ----------------->          su3_vector            */
/* ----------->                wilson_vector         */
/* ----->                      color_wilson_vector   */

/* Object with two Dirac and two color indices. A given element
   of a "wilson_propagator" is accessed by
   object.c[color1].d[spin1].d[spin2].c[color2].real , etc.
   As alway, "d" denotes a Dirac index and "c" a color index.
   "1" refers to the source, "2" to the sink.
*/

typedef struct { fsu3_vector d[4]; } fwilson_vector;
typedef struct { fsu3_vector h[2]; } fhalf_wilson_vector;
typedef struct { fwilson_vector c[3]; } fcolor_wilson_vector;
typedef struct { fwilson_vector d[4]; } fspin_wilson_vector;
typedef struct { fcolor_wilson_vector d[4]; } fwilson_matrix;
typedef struct { fspin_wilson_vector c[3]; } fwilson_propagator;

typedef struct { dsu3_vector d[4]; } dwilson_vector;
typedef struct { dsu3_vector h[2]; } dhalf_wilson_vector;
typedef struct { dwilson_vector c[3]; } dcolor_wilson_vector;
typedef struct { dwilson_vector d[4]; } dspin_wilson_vector;
typedef struct { dcolor_wilson_vector d[4]; } dwilson_matrix;
typedef struct { dspin_wilson_vector c[3]; } dwilson_propagator;

#if (PRECISION==1)

#define wilson_vector       fwilson_vector
#define half_wilson_vector  fhalf_wilson_vector
#define color_wilson_vector fcolor_wilson_vector
#define spin_wilson_vector  fspin_wilson_vector
#define wilson_matrix       fwilson_matrix
#define wilson_propagator   fwilson_propagator

#else

#define wilson_vector       dwilson_vector
#define half_wilson_vector  dhalf_wilson_vector
#define color_wilson_vector dcolor_wilson_vector
#define spin_wilson_vector  dspin_wilson_vector
#define wilson_matrix       dwilson_matrix
#define wilson_propagator   dwilson_propagator

#endif

/* For Wilson spectroscopy */
typedef struct {
  spin_wilson_vector ** swv;
  int nc;
} wilson_prop_field;


#define GAMMAFIVE -1    /* some integer which is not a direction */
#define PLUS 1          /* flags for selecting M or M_adjoint */
#define MINUS -1
/* Macros to multiply complex numbers by +-1 and +-i */
#define TIMESPLUSONE(a,b) { (b).real =  (a).real; (b).imag = (a).imag; }
#define TIMESMINUSONE(a,b) { (b).real =  -(a).real; (b).imag = -(a).imag; }
#define TIMESPLUSI(a,b) { (b).real = -(a).imag; (b).imag =  (a).real; }
#define TIMESMINUSI(a,b) { (b).real =  (a).imag; (b).imag = -(a).real; }

/*
* ROUTINES FOR SU(3) MATRIX OPERATIONS
*
* void mult_su3_nn(  su3_matrix *a, su3_matrix *b, su3_matrix *c  )
*	matrix multiply, no adjoints
*	files "m_mat_nn.c", "m_mat_nn.m4"
* void mult_su3_na( su3_matrix *a, su3_matrix *b, su3_matrix *c )
*	matrix multiply, second matrix is adjoint
*	files "m_mat_na.c", "m_mat_na.m4"
* void mult_su3_an( su3_matrix *a, su3_matrix *b, su3_matrix *c )
*	matrix multiply, first matrix is adjoint
*	files "m_mat_an.c", "m_mat_an.m4"
* Real realtrace_su3(  su3_matrix *a, su3_matrix *b )
*	(Re(Tr( A_adjoint*B)) )
*	file "realtr.c"
* complex trace_su3( su3_matrix *a )
*	file "trace_su3.c"
* complex complextrace_su3( su3_matrix *a, su3_matrix *b )
*	(Tr( A_adjoint*B))
*	file "complextr.c"
* complex det_su3( su3_matrix *a )
*	file "det_su3.c"
* void add_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c )
*	file "addmat.c"
* void sub_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c )
*	file "submat.c"
* void scalar_mult_su3_matrix( su3_matrix *a, Real s, su3_matrix *b )
*	file "s_m_mat.c"
* void scalar_mult_add_su3_matrix( su3_matrix *a, su3_matrix *b,
*	Real s, su3_matrix *c)
*	file "s_m_a_mat.c"
* void scalar_mult_sub_su3_matrix( su3_matrix *a, su3_matrix *b,
*	Real s, su3_matrix *c)
*	file "s_m_s_mat.c"
* void c_scalar_mult_su3mat( su3_matrix *src, complex *phase, su3_matrix *dest)
*	file "cs_m_mat.c"
* void c_scalar_mult_add_su3mat( su3_matrix *m1, su3_matrix *m2,
*	complex *phase, su3_matrix *m3)
*	file "cs_m_a_mat.c"
* void c_scalar_mult_sub_su3mat( su3_matrix *m1, su3_matrix *m2,
*	complex *phase, su3_matrix *m3)
*	file "cs_m_s_mat.c"
* void su3_adjoint( su3_matrix *a, su3_matrix *b )
*	file "su3_adjoint.c"
* void make_anti_hermitian( su3_matrix *m3,  anti_hermitmat *ah3 )
*	file "make_ahmat.c"
* void random_anti_hermitian( anti_hermitmat *mat_antihermit, double_prn *prn_pt )
*	(prn_pt passed through to myrand())
*	file "rand_ahmat.c"
* void uncompress_anti_hermitian( anti_hermitmat *mat_anti, su3_matrix *mat )
*	file "uncmp_ahmat.c"
* void compress_anti_hermitian( su3_matrix *mat, anti_hermitmat *mat_anti)
*	file "cmp_ahmat.c"
* void clear_su3mat( su3_matrix *dest );
*       file clear_mat.c
*          dest <- 0.0
* void su3mat_copy( su3_matrix *a, su3_matrix *b )
*	file "su3mat_copy.c"
* void dumpmat( su3_matrix *m )
*       file "dumpmat.c"
* void dumptensor4( su3_tensor4 *m )
*       file "dumpmat.c"
* void eigen_su3_UdU( su3_matrix *U, Real *g0, Real *g1, Real *g2)
*       file "eigen_su3_UdU.c"
*
*
* ROUTINES FOR su3_vector OPERATIONS ( 3 COMPONENT COMPLEX )
*
* void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c )
*	( outer product of A and B)
*	file "su3_proj.c"
* complex su3_dot( su3_vector *a, su3_vector *b )
*	file "su3_dot.c"
* Real su3_rdot( su3_vector *a, su3_vector *b )
*	file "su3_rdot.c", "su3_rdot.m4"
* Real magsq_su3vec( su3_vector *a )
*	file "msq_su3vec.c", "msq_su3vec.m4"
* void su3vec_copy( su3_vector *a, su3_vector *b )
*	file "su3vec_copy.c"
* 
* void mult_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c )
*	 C  <-  A*B
*	file "m_matvec.c", "m_matvec.m4"
* void mult_su3_mat_vec_sum( su3_matrix *a, su3_vector *b, su3_vector *c )
*	 C  <-  C + A*B
*	file "m_matvec_s.c", "m_matvec_s.m4"
* void mult_su3_mat_vec_sum_4dir( su3_matrix *a, su3_vector *b0,
*	su3_vector *b1, su3_vector *b2, su3_vector *b3, su3_vector *c )
*	file "m_mv_s_4dir.c", "m_mv_s_4dir.m4"
*	file "m_mv_s_4di2.m4" is alternate version with pipelined loads.
*	Multiply four su3_vectors by elements of an array of su3_matrices,
*	sum results.
*	C <- A[0]*B0 + A[1]*B1 + A[2]*B2 + A[3]*B3
* void mult_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c )
*	file "m_matvec_ns.c"
* void mult_adj_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c )
*	file "m_amatvec.c", "m_amatvec.m4"
* void mult_adj_su3_mat_vec_4dir( su3_matrix *a, su3_vector *b, su3_vector *c )
*	file "m_amv_4dir.c", "m_amv_4dir.m4"
*	file "m_amv_4di2.m4" is alternate version with pipelined loads.
*	Multiply an su3_vector by adjoints of elements of an array 
*	of su3_matrices, results in an array of su3_vectors.
*	C[i] <- A_adjoint[i]*B, i = 0,1,2,3
* void mult_adj_su3_mat_4vec( su3_matrix *mat, su3_vector *src,
*			     su3_vector *dest0, su3_vector *dest1, 
*			     su3_vector *dest2, su3_vector *dest3  ) ;
*       file m_amv_4vec.c
*	Same as above, but result vectors need not be in an array.
* void mult_adj_su3_mat_vec_sum( su3_matrix *a, su3_vector *b, su3_vector *c )
*	file "m_amatvec_s.c"
* void mult_adj_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c )
*	file "m_amatvec_ns.c"
* void add_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c )
*	file "addvec.c", "addvec.m4"
* void sub_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c )
*	file "subvec.c", "subvec.m4"
* void sub_four_su3_vecs( su3_vector *a, su3_vector *b1, su3_vector *b2,
*   su3_vector *b3, su3_vector *b4 )
*	file "sub4vecs.c", "sub4vecs.m4"
*
* void scalar_mult_su3_vector( su3_vector *a, Real s, su3_vector *c )
*	file "s_m_vec.c"
* void scalar_mult_add_su3_vector( su3_vector *a, su3_vector *b, Real s,
*	su3_vector *c)
*	file "s_m_a_vec.c", "s_m_a_vec.m4"
* void scalar_mult_sum_su3_vector( su3_vector *a, su3_vector *b, Real s )
*	file "s_m_sum_vec.c", "s_m_sum_vec.m4"
* void scalar_mult_sub_su3_vector( su3_vector *a, su3_vector *b, Real s,
*	su3_vector *c )
*	file "s_m_s_vec.c"
* void c_scalar_mult_su3vec( su3_vector *src, complex *phase, su3_vector *dest )
*	file "cs_m_vec.c"
* void c_scalar_mult_add_su3vec( su3_vector *v1, complex *phase, su3_vector *v2)
*	file "cs_m_a_vec.c"
* void c_scalar_mult_sub_su3vec( su3_vector *v1, complex *phase, su3_vector *v2)
*	file "cs_m_s_vec.c"
* void dumpvec( su3_vector *v )
*       file "dumpvec.c"
* void clearvec( su3_vector *v )
*       file "clearvec.c"
* 
* ROUTINES MIXING SU(2) and SU(3)
* 
* void left_su2_hit_n(su2_matrix *u,int p,int q,su3_matrix *link)
*       file "l_su2_hit_n.c"
* void right_su2_hit_a(su2_matrix *u,int p,int q,su3_matrix *link)
*       file "r_su2_hit_a.c"
* void dumpsu2(su2_matrix *u)
*       file "dumpsu2.c"
* void mult_su2_mat_vec_elem_n(su2_matrix *u,complex *x0,complex *x1)
*       file "m_su2_mat_vec_n.c"
* void mult_su2_mat_vec_elem_a(su2_matrix *u,complex *x0,complex *x1);
*       file "m_su2_mat_vec_a.c"
*
* ROUTINES FOR WILSON VECTORS
*
* void mult_mat_wilson_vec( su3_matrix *mat, wilson_vector *src,
*	wilson_vector *dest );
*	file m_mat_wvec.c
*	   dest <- mat*src
* void mult_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
*	half_wilson_vector *dest );
*	file m_mat_hwvec.c
*	   dest <- mat*src
* void mult_adj_mat_wilson_vec( su3_matrix *mat, wilson_vector *src,
*	wilson_vector *dest)
*	file m_amat_wvec.c
*	   dest <- mat_adjoint*src
* void mult_adj_su3_mat_hwvec su3_matrix *mat,
*	half_wilson_vector *src, half_wilson_vector *dest )
*	file m_amat_hwvec.c
*	   dest <- mat_adjoint*src
*
* void add_wilson_vector( wilson_vector *src1, wilson_vector *src2,
*	wilson_vector *dest );
*	file add_wvec.c
*	   dest <- src1+src2
* void sub_wilson_vector( wilson_vector *src1, wilson_vector *src2,
*       wilson_vector *dest );
*	file sub_wvec.c
*	   dest <- src1-src2
*
* void scalar_mult_wvec  wilson_vector *src, Real s, wilson_vector *dest )
*	file s_m_wvec.c
*	   dest <- s*src
* void scalar_mult_hwvec( half_wilson_vector *src, Real s,
*          half_wilson_vector *dest)
*	file s_m_hwvec.c
*	   dest <- s*src
* Real magsq_wvec( wilson_vector *src );
*	file msq_wvec.c
*	   s <- squared magnitude of src
* complex wvec_dot( wilson_vector *src1, wilson_vector *src2 );
*	file wvec_dot.c
*	   c <- dot product of src1 and src2
* complex wvec2_dot( wilson_vector *src1, wilson_vector *src2 );
*	file wvec2_dot.c
*	   c <- dot product of src1 and src2, Used only in Claude's
*	   mrilu.c,  I don't know what the difference is.  DT.
* Real wvec_rdot( wilson_vector *a, wilson_vector *b )
*	wilson_vector *a,*b;
*	file "wvec_rdot.c", "wvec_rdot.m4"
*	   r <- real part of dot product of src1 and src2
* void scalar_mult_add_wvec(  wilson_vector *src1,*src2, Real s,
*           wilson_vector *dest)
*	file s_m_a_wvec.c
*	   dest <- src1 + s*src2
* void scalar_mult_addtm_wvec( wilson_vector *src1,*src2, Real s,
*           wilson_vector *dest)
*	file s_m_atm_wvec.c
*	   dest <- -src1 + s*src2   ("atm"="add to minus")
* void c_scalar_mult_wvec( wilson_vector *v1, complex *phase,
*	    wilson_vector *v2)
*       file "cs_m_wvec.c"
* void c_scalar_mult_add_wvec( wilson_vector *v1, wilson_vector *v2,
*	complex *phase, wilson_vector *v3)
*       file "cs_m_a_wvec.c"
* void c_scalar_mult_add_wvec2( wilson_vector *v1, wilson_vector *v2,
*	complex scalar, wilson_vector *v3)
*       file "cs_m_a_wvec2.c"
*	differs from previous one: value of scalar, not address is arg.
* void wp_shrink( wilson_vector *src, half_wilson_vector *dest,
*	int dir, int sign );
*	file wp_shrink.c , wp_shrink.m4
*	   if(dir = [XYZT]UP) dest <- components of src along eigenvectors
*		of gamma_dir with eigenvalue +1
*	   if(dir = [XYZT]DOWN) dest <- components of src along eigenvectors
*		of gamma_dir with eigenvalue -1
*	   if(sign==MINUS)switch roles of +1 and -1
* void wp_shrink_4dir( wilson_vector *a,  half_wilson_vector *b1,
*	half_wilson_vector *b2, half_wilson_vector *b3,
*	half_wilson_vector *b4,	int sign );
*	file wp_shrink4.c wp_shrink4.m4
*	  Shrink A in X,Y,Z,T directions respectively, results in B1-B4
* void wp_grow(  half_wilson_vector *src, wilson_vector *dest,
*	int dir, int sign );
*	file wp_grow.c , wp_grow.m4
*	   if(dir = [XYZT]UP) dest <- components of src times eigenvectors
*		of gamma_dir with eigenvalue +1
*	   if(dir = [XYZT]DOWN) dest <- components of src times eigenvectors
*		of gamma_dir with eigenvalue -1
*	   if(sign==MINUS)switch roles of +1 and -1
*  	    Note: wp_shrink( +-dir) followed by wp_grow( +-dir) amounts to
*		multiplication by 1+-gamma_dir, or 1-+gamma_dir if sign=MINUS
* void wp_grow_add( half_wilson_vector *src, wilson_vector *dest,
*	int dir, int sign );
*	file wp_grow_a.c , wp_grow_a.m4
*	   wp_grow, and add result to previous contents of dest.
* void grow_add_four_wvecs( wilson_vector *a, half_wilson_vector *b1,
*	half_wilson_vector *b2, half_wilson_vector *b3,
*	half_wilson_vector *b4, int sign, int sum );
*	file grow4wvecs.c grow4wvecs.m4
*         If sum==0
*	  Grow b1-b4 in X,Y,Z,T directions respectively, sum of results in A
*         If sum==1
*	  Grow b1-b4 in X,Y,Z,T directions respectively, add to current A
*
* void mult_by_gamma( wilson_vector *src, wilson_vector *dest, int dir );
*	file mb_gamma.c
*	   dest <- gamma[dir] * src,  dir=[XYZT]UP,GAMMAFIVE
* void mult_by_gamma_left( wilson_matrix *src,  wilson_matrix *dest, int dir );
*	file mb_gamma_l.c
*	   dest <- gamma[dir] * src,  dir=[XYZT]UP,GAMMAFIVE
*	   acts on first index of matrix
* void mult_by_gamma_right( wilson_matrix *src,  wilson_matrix *dest, int dir );
*	file mb_gamma_r.c
*	   dest_ij <- gamma[dir]_ik * src_jk,  dir=[XYZT]UP,GAMMAFIVE
*	   acts on second index of matrix
*
* void mult_swv_by_gamma_l( spin_wilson_vector *src,  spin_wilson_vector *dest, int dir );
*	file mswvb_gamma_l.c
*	   same as mult_by_gamma_left, but acts on a spin_wilson_vector
*
* void mult_swv_by_gamma_r( spin_wilson_vector *src,  spin_wilson_vector *dest, int dir );
*	file mswvb_gamma_r.c
*	   same as mult_by_gamma_right, but acts on a spin_wilson_vector
*
* void su3_projector_w( wilson_vector *a, wilson_vector *b, su3_matrix *c )
*	sum over spins of outer product of A.d[s] and B.d[s]  - a three
*	  by three complex matrix
*	file "su3_proj_w.c"
* void clear_wvec( wilson_vector *dest );
*	file clear_wvec.c
*	   dest <- 0.0
* void copy_wvec( wilson_vector *src, wilson_vector *dest );
*	file copy_wvec.c
*	   dest <- src
* void dump_wilson_vec( wilson_vector *src );
*	file dump_wvec.c
*	   print out a wilson vector
*
* MISCELLANEOUS ROUTINES
*
* Real gaussian_rand_no( double_prn *prn_pt )
*	file "gaussrand.c"
* complex complex_gaussian_rand_no( double_prn *prn_pt )
*	file "gaussrand.c"
* Real z2_rand_no( double_prn *prn_pt );
*	file "z2rand.c"
* void byterevn(int32type w[], int n)
* void byterevn64(int32type w[], int n)
*
*/

Real realtrace_su3(  su3_matrix *a, su3_matrix *b );
complex trace_su3(  su3_matrix *a );
complex complextrace_su3( su3_matrix *a, su3_matrix *b );
complex det_su3( su3_matrix *a );
void sub_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void scalar_mult_su3_matrix( su3_matrix *src, Real scalar, su3_matrix *dest);
void scalar_mult_sub_su3_matrix( su3_matrix *src1, su3_matrix *src2,
	Real scalar, su3_matrix *dest);
void c_scalar_mult_su3mat( su3_matrix *src, complex *scalar,
	su3_matrix *dest);
void c_scalar_mult_add_su3mat( su3_matrix *src1, su3_matrix *src2,
	complex *scalar, su3_matrix *dest);
void c_scalar_mult_sub_su3mat( su3_matrix *src1, su3_matrix *src2,
	complex *scalar, su3_matrix *dest);
void su3_adjoint( su3_matrix *a, su3_matrix *b );
void make_anti_hermitian( su3_matrix *m3, anti_hermitmat *ah3 );
void random_anti_hermitian( anti_hermitmat *mat_antihermit, double_prn *prn_pt );
void uncompress_anti_hermitian( anti_hermitmat *mat_anti, su3_matrix *mat );
void compress_anti_hermitian( su3_matrix *mat, anti_hermitmat *mat_anti);
void clear_su3mat( su3_matrix *dest );
void su3mat_copy( su3_matrix *a, su3_matrix *b );
void dumpmat( su3_matrix *m );
void dumptensor4( su3_tensor4 *m );
void eigen_su3_UdU( su3_matrix *U, Real *g0, Real *g1, Real *g2);

complex su3_dot( su3_vector *a, su3_vector *b );
void su3vec_copy( su3_vector *a, su3_vector *b );
void dumpvec( su3_vector *v );
void clearvec( su3_vector *v );

void mult_su3_mat_vec_sum(  su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_vec_sum( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c );

void scalar_mult_su3_vector(  su3_vector *src, Real scalar, 
	su3_vector *dest);
void scalar_mult_sum_su3_vector( su3_vector *src1, su3_vector *src2,
	Real scalar);
void scalar_mult_sub_su3_vector( su3_vector *src1, su3_vector *src2,
	Real scalar, su3_vector *dest);
void scalar_mult_wvec( wilson_vector *src, Real s, wilson_vector *dest );
void scalar_mult_hwvec( half_wilson_vector *src, Real s, 
    half_wilson_vector *dest );
void scalar_mult_add_wvec( wilson_vector *src1, wilson_vector *src2,
	Real scalar, wilson_vector *dest );
void scalar_mult_addtm_wvec( wilson_vector *src1, wilson_vector *src2,
	Real scalar, wilson_vector *dest );
void c_scalar_mult_wvec(wilson_vector *src1, complex *phase,
	wilson_vector *dest );
void c_scalar_mult_add_wvec(wilson_vector *src1, wilson_vector *src2,
	complex *phase, wilson_vector *dest );
void c_scalar_mult_add_wvec2(wilson_vector *src1, wilson_vector *src2,
	complex s, wilson_vector *dest );
void c_scalar_mult_su3vec( su3_vector *src, complex *phase, su3_vector *dest );
void c_scalar_mult_add_su3vec(su3_vector *v1, complex *phase, su3_vector *v2);
void c_scalar_mult_sub_su3vec(su3_vector *v1, complex *phase, su3_vector *v2);

void left_su2_hit_n(su2_matrix *u,int p,int q,su3_matrix *link);
void right_su2_hit_a(su2_matrix *u,int p,int q,su3_matrix *link);
void dumpsu2(su2_matrix *u);
void mult_su2_mat_vec_elem_n(su2_matrix *u,complex *x0,complex *x1);
void mult_su2_mat_vec_elem_a(su2_matrix *u,complex *x0,complex *x1);

void mult_mat_wilson_vec( su3_matrix *mat, wilson_vector *src,
	wilson_vector *dest );
void mult_adj_mat_wilson_vec( su3_matrix *mat, wilson_vector *src,
	wilson_vector *dest);

void add_wilson_vector( wilson_vector *src1, wilson_vector *src2,
	wilson_vector *dest );
void sub_wilson_vector( wilson_vector *src1, wilson_vector *src2,
       wilson_vector *dest );
Real magsq_wvec( wilson_vector *src );
complex wvec_dot( wilson_vector *src1, wilson_vector *src2 );
complex wvec2_dot( wilson_vector *src1, wilson_vector *src2 );
Real wvec_rdot( wilson_vector *a, wilson_vector *b );

void wp_shrink( wilson_vector *src, half_wilson_vector *dest,
	int dir, int sign );
void wp_grow(  half_wilson_vector *src, wilson_vector *dest,
	int dir, int sign );
void wp_grow_add( half_wilson_vector *src, wilson_vector *dest,
	int dir, int sign );
void mult_by_gamma( wilson_vector *src, wilson_vector *dest, int dir );
void mult_by_gamma_left( wilson_matrix *src,  wilson_matrix *dest, int dir );
void mult_by_gamma_right( wilson_matrix *src,  wilson_matrix *dest, int dir );
void mult_swv_by_gamma_l(spin_wilson_vector *src, spin_wilson_vector *dest, int dir);
void mult_swv_by_gamma_r(spin_wilson_vector *src, spin_wilson_vector *dest, int dir);
void su3_projector_w( wilson_vector *a, wilson_vector *b, su3_matrix *c );
void clear_wvec( wilson_vector *dest );
void copy_wvec( wilson_vector *src, wilson_vector *dest );
void dump_wilson_vec( wilson_vector *src );

Real gaussian_rand_no( double_prn *prn_pt );
complex complex_gaussian_rand_no( double_prn *prn_pt );
Real z2_rand_no( double_prn *prn_pt );
#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);


/********************************************************************/
/* Inline macros                                                    */
/********************************************************************/

/* We have optional SSE and C inline macros for selected library
   routines */

/* All macros are defined and available for selective inlining.  To
   invoke them selectively, define the macro SSE_INLINE and then
   modify the code by replacing the function call with the appropriate
   macro */

/* Inlining can also be implemented globally without requiring any
   changes in the code.  SSE macros are inlined globally by defining
   SSE_GLOBAL_INLINE.  */

/* C macros are similarly inlined selectively with C_INLINE and
   globally with C_GLOBAL_INLINE */

/* SSE macros for selected library functions are available in single
   and double precision */

#if defined SSE_INLINE || defined SSE_GLOBAL_INLINE
#if (PRECISION==1)
#ifdef SSEOPTERON
#include "../sse_opteron/include/inline_sse.h"
#else
#include "../sse/include/inline_sse.h"
#endif
#else
#include "../sse2/include/inline_sse.h"
#endif
#endif

/* C macros are similarly available for selected library functions in
   both precisions */

#if defined C_INLINE || defined C_GLOBAL_INLINE
#include "../libraries/include/inline_C.h"
#endif

/* The following definitions cause the macros to be invoked globally
   if SSE_INLINE and/or C_INLINE is defined.  Note that in each stanza
   the same library routines are repeated, but with either a macro
   definition when available or a prototype definition. */

/* Note that SSE takes precedence over C */

#if defined SSE_GLOBAL_INLINE

/********************************************************************/
/* Our available single-precision SSE macros */

#if (PRECISION == 1)

#define add_su3_vector(...) _inline_sse_add_su3_vector(__VA_ARGS__)
#define mult_su3_nn(...) _inline_sse_mult_su3_nn(__VA_ARGS__)
#define mult_su3_na(...) _inline_sse_mult_su3_na(__VA_ARGS__)
#define mult_su3_an(...) _inline_sse_mult_su3_an(__VA_ARGS__)
#define mult_su3_mat_vec(...) _inline_sse_mult_su3_mat_vec(__VA_ARGS__)
#define mult_adj_su3_mat_vec(...) _inline_sse_mult_adj_su3_mat_vec(__VA_ARGS__)
#define mult_adj_su3_mat_vec_4dir(...) _inline_sse_mult_adj_su3_mat_vec_4dir(__VA_ARGS__)
#define mult_adj_su3_mat_4vec(...) _inline_sse_mult_adj_su3_mat_4vec(__VA_ARGS__)
#define mult_adj_su3_mat_hwvec(...) _inline_sse_mult_adj_su3_mat_hwvec(__VA_ARGS__)
#define mult_su3_mat_hwvec(...) _inline_sse_mult_su3_mat_hwvec(__VA_ARGS__)
#define mult_su3_mat_vec_sum_4dir(...) _inline_sse_mult_su3_mat_vec_sum_4dir(__VA_ARGS__)
#define scalar_mult_add_su3_matrix(a,b,c,d) {Real _temp = c; _inline_sse_scalar_mult_add_su3_matrix(a,b,_temp,d);}
#define scalar_mult_add_su3_vector(a,b,c,d) {Real _temp = c; _inline_sse_scalar_mult_add_su3_vector(a,b,_temp,d);}
#define su3_projector(...) _inline_sse_su3_projector(__VA_ARGS__)
#define sub_four_su3_vecs(...) _inline_sse_sub_four_su3_vecs(__VA_ARGS__)

/********************************************************************/

#else // PRECISION == 2

/********************************************************************/
/* Our available double-precision SSE macros */

#define mult_su3_nn(...) _inline_sse_mult_su3_nn(__VA_ARGS__)
#define mult_su3_na(...) _inline_sse_mult_su3_na(__VA_ARGS__)
#define mult_su3_an(...) _inline_sse_mult_su3_an(__VA_ARGS__)
#define mult_su3_mat_vec(...) _inline_sse_mult_su3_mat_vec(__VA_ARGS__)
#define mult_adj_su3_mat_vec(...) _inline_sse_mult_adj_su3_mat_vec(__VA_ARGS__)
#define mult_adj_su3_mat_vec_4dir(...) _inline_sse_mult_adj_su3_mat_vec_4dir(__VA_ARGS__)
#define mult_adj_su3_mat_4vec(...) _inline_sse_mult_adj_su3_mat_4vec(__VA_ARGS__)
#define mult_adj_su3_mat_hwvec(...) _inline_sse_mult_adj_su3_mat_hwvec(__VA_ARGS__)
#define mult_su3_mat_hwvec(...) _inline_sse_mult_su3_mat_hwvec(__VA_ARGS__)
#define mult_su3_mat_vec_sum_4dir(...) _inline_sse_mult_su3_mat_vec_sum_4dir(__VA_ARGS__)
#define su3_projector(...) _inline_sse_su3_projector(__VA_ARGS__)

/********************************************************************/

#endif // PRECISION
#endif // SSE_GLOBAL_INLINE

#if defined C_GLOBAL_INLINE

/********************************************************************/
/* Our available C-inline macros */

#ifndef add_su3_matrix
#define add_su3_matrix(...) _inline_C_add_su3_matrix(__VA_ARGS__)
#endif

#ifndef add_su3_vector
#define add_su3_vector(...) _inline_C_add_su3_vector(__VA_ARGS__)
#endif

#ifndef grow_add_four_wvecs
#define grow_add_four_wvecs(...)  _inline_C_grow_add_four_wvecs(__VA_ARGS__)
#endif

#ifndef magsq_su3vec
#define magsq_su3vec(...) _inline_C_magsq_su3vec(__VA_ARGS__)
#endif

#ifndef mult_su3_nn
#define mult_su3_nn(...) _inline_C_mult_su3_nn(__VA_ARGS__)
#endif

#ifndef mult_su3_na
#define mult_su3_na(...) _inline_C_mult_su3_na(__VA_ARGS__)
#endif

#ifndef mult_adj_su3_mat_hwvec
#define mult_adj_su3_mat_hwvec(...) _inline_C_mult_adj_su3_mat_hwvec(__VA_ARGS__)
#endif

#ifndef mult_su3_mat_hwvec
#define mult_su3_mat_hwvec(...) _inline_C_mult_su3_mat_hwvec(__VA_ARGS__)
#endif

#ifndef  mult_su3_mat_vec_sum_4dir
#define  mult_su3_mat_vec_sum_4dir(...) _inline_C_mult_su3_mat_vec_sum_4dir(__VA_ARGS__)
#endif

#ifndef scalar_mult_add_hwvec_proj
#define scalar_mult_add_hwvec_proj(...) _inline_C_scalar_mult_add_hwvec_proj(__VA_ARGS__)
#endif

#ifndef scalar_mult_add_su3_matrix
#define scalar_mult_add_su3_matrix(a,b,c,d) {Real _temp = c; _inline_C_scalar_mult_add_su3_matrix(a,b,_temp,d);}
#endif

#ifndef scalar_mult_add_su3_vector
#define scalar_mult_add_su3_vector(a,b,c,d) {Real _temp = c; _inline_C_scalar_mult_add_su3_vector(a,b,_temp,d);}
#endif

#ifndef su3_projector
#define su3_projector(...) _inline_C_su3_projector(__VA_ARGS__)
#endif

#ifndef su3_rdot
#define su3_rdot(...) _inline_C_su3_rdot(__VA_ARGS__)
#endif

#ifndef su3mat_copy
//#define su3mat_copy(...) _inline_C_su3mat_copy(__VA_ARGS__)
#endif

#ifndef sub_su3_matrix
#define sub_su3_matrix(...) _inline_C_sub_su3_matrix(__VA_ARGS__)
#endif

#ifndef sub_su3_vector
#define sub_su3_vector(...) _inline_C_sub_su3_vector(__VA_ARGS__)
#endif

#ifndef sub_four_su3vecs
#define sub_four_su3vecs(...) _inline_C_sub_four_su3vecs(__VA_ARGS__)
#endif

#ifndef wp_shrink_4dir
#define wp_shrink_4dir(...) _inline_C_wp_shrink_4dir(__VA_ARGS__)
#endif

/********************************************************************/

#endif // C_GLOBAL_INLINE

/********************************************************************/
/* Use standard prototypes if macros are not defined */

#ifndef add_su3_matrix
void add_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c );
#endif

#ifndef add_su3_vector
void add_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c );
#endif

#ifndef grow_add_four_wvecs
void grow_add_four_wvecs( wilson_vector *a, half_wilson_vector *b1,
	half_wilson_vector *b2, half_wilson_vector *b3,
	half_wilson_vector *b4, int sign, int sum );
#endif

#ifndef magsq_su3vec
Real magsq_su3vec( su3_vector *a );
#endif

#ifndef mult_su3_nn 
void mult_su3_nn ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
#endif

#ifndef mult_su3_na 
void mult_su3_na ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
#endif

#ifndef mult_su3_an 
void mult_su3_an ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
#endif

#ifndef mult_su3_mat_vec
void mult_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c );
#endif

#ifndef mult_adj_su3_mat_vec
void mult_adj_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c );
#endif

#ifndef mult_adj_su3_mat_vec_4dir
void mult_adj_su3_mat_vec_4dir( su3_matrix *a, su3_vector *b, su3_vector *c );
#endif

#ifndef mult_adj_su3_mat_4vec
void mult_adj_su3_mat_4vec( su3_matrix *mat, su3_vector *src,
			    su3_vector *dest0, su3_vector *dest1, 
			    su3_vector *dest2, su3_vector *dest3  ) ;
#endif

#ifndef mult_adj_su3_mat_hwvec
void mult_adj_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
	half_wilson_vector *dest );
#endif

#ifndef mult_su3_mat_vec_sum_4dir
void mult_su3_mat_vec_sum_4dir( su3_matrix *a, su3_vector *b0,
	su3_vector *b1, su3_vector *b2, su3_vector *b3, su3_vector *c );
#endif

#ifndef mult_su3_mat_hwvec
void mult_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
	half_wilson_vector *dest );
#endif

#ifndef scalar_mult_add_hwvec_proj
void scalar_mult_add_hwvec_proj( su3_matrix * const a, 
				 half_wilson_vector * const b, 
				 half_wilson_vector * const c, 
				 Real * const s, su3_matrix *d );
#endif

#ifndef scalar_mult_add_su3_matrix
void scalar_mult_add_su3_matrix( su3_matrix *src1, su3_matrix *src2,
	Real scalar, su3_matrix *dest);
#endif

#ifndef scalar_mult_add_su3_vector
void scalar_mult_add_su3_vector( su3_vector *src1, su3_vector *src2,
	Real scalar, su3_vector *dest);
#endif

#ifndef su3_projector
void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c );
#endif

#ifndef su3_rdot
Real su3_rdot( su3_vector *a, su3_vector *b );
#endif

#ifndef sub_four_su3_vecs
void sub_four_su3_vecs( su3_vector *a, su3_vector *b1, su3_vector *b2,
	su3_vector *b3, su3_vector *b4 );
#endif

#ifndef sub_su3_matrix
void sub_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c );
#endif

#ifndef sub_su3_vector
void sub_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c );
#endif

#ifndef wp_shrink_4dir
void wp_shrink_4dir( wilson_vector *a,  half_wilson_vector *b1,
	half_wilson_vector *b2, half_wilson_vector *b3,
	half_wilson_vector *b4,	int sign );
#endif

/********************************************************************/

#endif /* _SU3_H */
