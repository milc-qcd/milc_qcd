#ifndef _SU3_H
#define _SU3_H

#include "../include/complex.h"
#include "../include/random.h"

/******************************  su3.h **********************************
*									*
*  Defines and subroutine declarations for SU3 simulation		*
*  MIMD version 6 							*
*									*
*/
/* SU(3) */
typedef struct { fcomplex e[3][3]; } fsu3_matrix;
typedef struct { dcomplex e[3][3]; } dsu3_matrix;
typedef struct { complex e[3][3]; } su3_matrix;
typedef struct { fcomplex c[3]; } fsu3_vector;
typedef struct { dcomplex c[3]; } dsu3_vector;
typedef struct { complex c[3]; } su3_vector;
typedef struct
  { complex m01,m02,m12; Real m00im,m11im,m22im; Real space; } anti_hermitmat;

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

typedef struct { su3_vector d[4]; } wilson_vector;
typedef struct { su3_vector h[2]; } half_wilson_vector;
typedef struct { wilson_vector c[3]; } color_wilson_vector;
typedef struct { wilson_vector d[4]; } spin_wilson_vector;
typedef struct { color_wilson_vector d[4]; } wilson_matrix;
typedef struct { spin_wilson_vector c[3]; } wilson_propagator;

typedef struct { fsu3_vector d[4]; } fwilson_vector;
typedef struct { fsu3_vector h[2]; } fhalf_wilson_vector;
typedef struct { fwilson_vector c[3]; } fcolor_wilson_vector;
typedef struct { fwilson_vector d[4]; } fspin_wilson_vector;
typedef struct { fcolor_wilson_vector d[4]; } fwilson_matrix;
typedef struct { fspin_wilson_vector c[3]; } fwilson_propagator;

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
* void byterevn(int32type w[], int n)
* void byterevn64(int32type w[], int n)
*
*/

Real realtrace_su3(  su3_matrix *a, su3_matrix *b );
complex trace_su3(  su3_matrix *a );
complex complextrace_su3( su3_matrix *a, su3_matrix *b );
complex det_su3( su3_matrix *a );
void add_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c );
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

complex su3_dot( su3_vector *a, su3_vector *b );
Real su3_rdot( su3_vector *a, su3_vector *b );
Real magsq_su3vec( su3_vector *a );
void su3vec_copy( su3_vector *a, su3_vector *b );
void dumpvec( su3_vector *v );
void clearvec( su3_vector *v );

void mult_su3_mat_vec_sum(  su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_vec_sum( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c );

void sub_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c );

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
void wp_shrink_4dir( wilson_vector *a,  half_wilson_vector *b1,
	half_wilson_vector *b2, half_wilson_vector *b3,
	half_wilson_vector *b4,	int sign );
void wp_grow(  half_wilson_vector *src, wilson_vector *dest,
	int dir, int sign );
void wp_grow_add( half_wilson_vector *src, wilson_vector *dest,
	int dir, int sign );
void grow_add_four_wvecs( wilson_vector *a, half_wilson_vector *b1,
	half_wilson_vector *b2, half_wilson_vector *b3,
	half_wilson_vector *b4, int sign, int sum );
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
#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);

/* For inserting Don Holmgren's SSE versions of some library routines */
/* Use only for Gnu C on the P3 or P4 */

/* The header inline_sse.h defines macros that can be used to replace
   subroutine calls by inline assembly code.  Replacement occurs in
   the compilation if the macro SSE_INLINE is also defined. */

/* Define SSE ASM inline macros and, if SSE_INLINE is defined, define
   macros that replace the library calls listed below */

#if defined SSE

#if PRECISION==1
#include "../sse/include/inline_sse.h"
#else
#include "../sse2/include/inline_sse.h"
/* The following routines are not currently available in double precision */
void add_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c );
void scalar_mult_add_su3_vector( su3_vector *src1, su3_vector *src2,
	Real scalar, su3_vector *dest);
void scalar_mult_add_su3_matrix( su3_matrix *src1, su3_matrix *src2,
	Real scalar, su3_matrix *dest);
void sub_four_su3_vecs( su3_vector *a, su3_vector *b1, su3_vector *b2,
	su3_vector *b3, su3_vector *b4 );
#endif

#endif

/* If SSE_INLINE is not defined, the subroutine calls are used and the
   assembly code macros must be invoked individually with their
   otherwise internal names.  See sse/include/inline_sse.h.  This
   allows selective testing */

/* The following inline macros are not currently available in double
   precision so we keep the prototypes for external linkage */

#if ! defined SSE_INLINE

/* The usual case: prototypes for library calls when these are
   externally linked */

void add_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c );
void scalar_mult_add_su3_vector( su3_vector *src1, su3_vector *src2,
	Real scalar, su3_vector *dest);
void scalar_mult_add_su3_matrix( su3_matrix *src1, su3_matrix *src2,
	Real scalar, su3_matrix *dest);
void sub_four_su3_vecs( su3_vector *a, su3_vector *b1, su3_vector *b2,
	su3_vector *b3, su3_vector *b4 );

void mult_su3_nn ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void mult_su3_na ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void mult_su3_an ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void mult_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_su3_mat_vec_sum_4dir( su3_matrix *a, su3_vector *b0,
	su3_vector *b1, su3_vector *b2, su3_vector *b3, su3_vector *c );
void mult_adj_su3_mat_vec_4dir( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_4vec( su3_matrix *mat, su3_vector *src,
			    su3_vector *dest0, su3_vector *dest1, 
			    su3_vector *dest2, su3_vector *dest3  ) ;
void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c );
void mult_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
	half_wilson_vector *dest );
void mult_adj_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
	half_wilson_vector *dest );

#endif

#endif /* _SU3_H */
