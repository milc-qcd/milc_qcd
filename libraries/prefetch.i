# 1 "prefetch.c"
 

# 1 "../include/config.h" 1


 

 


 
 
 

 
 
 



 
 
 

 
 
 
 

 
 
 

 
 

 


 


 


 
   
 



# 3 "prefetch.c" 2







# 1 "prefetch32.c" 1
 













 




# 1 "../include/complex.h" 1



 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 






 







typedef struct {	    
   float real;		    
   float imag;
} complex;
typedef struct {            
   double real;		    
   double imag;
} double_complex;

 
 






 
 


 
complex cmplx(  float x, float y );
complex cadd( complex *a, complex *b );
complex cmul( complex *a, complex *b );
complex csub( complex *a, complex *b );
complex cdiv( complex *a, complex *b );
complex conjg( complex *a );
complex cexp( complex *a );       
complex clog( complex *a );        
complex csqrt( complex *z );        
complex ce_itheta( float theta );    

double_complex dcmplx( double x, double y );
double_complex dcadd( double_complex *a, double_complex *b );
double_complex dcmul( double_complex *a, double_complex *b );
double_complex dcsub( double_complex *a, double_complex *b );
double_complex dcdiv( double_complex *a, double_complex *b );
double_complex dconjg(  double_complex *a );
double_complex dcexp(  double_complex *a ); 
double_complex dclog(  double_complex *a );  
double_complex dcsqrt( double_complex *z );  
double_complex dce_itheta( double theta );   

 
								 

								 

								 


								 

								 


								 


								 

								 


								 


								 



								 


								 


								 


								 

								 

								 

								 

								 


                                                                


                                                                 




# 20 "prefetch32.c" 2

# 1 "../include/su3.h" 1




# 1 "../include/../include/random.h" 1



 

typedef struct {
   
  unsigned long r0,r1,r2,r3,r4,r5,r6;
  unsigned long multiplier,addend,ic_state;
  float scale;
} double_prn;

 

float myrand(double_prn *prn_pt);


# 5 "../include/su3.h" 2


 





 
typedef struct { complex e[3][3]; } su3_matrix;
typedef struct { complex c[3]; } su3_vector;
typedef struct
  { complex m01,m02,m12; float m00im,m11im,m22im; float space; } anti_hermitmat;

 
typedef struct { complex e[2][2]; } su2_matrix;

 

 
 
 
 
 
 
 
 
 
 
 
 
 
 

 






typedef struct { su3_vector d[4]; } wilson_vector;
typedef struct { su3_vector h[2]; } half_wilson_vector;
typedef struct { wilson_vector c[3]; } color_wilson_vector;
typedef struct { wilson_vector d[4]; } spin_wilson_vector;
typedef struct { color_wilson_vector d[4]; } wilson_matrix;
typedef struct { spin_wilson_vector c[3]; } wilson_propagator;




 





 







































































































































































































































































































void mult_su3_nn ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void mult_su3_na ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void mult_su3_an ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
float realtrace_su3(  su3_matrix *a, su3_matrix *b );
complex trace_su3(  su3_matrix *a );
complex complextrace_su3( su3_matrix *a, su3_matrix *b );
complex det_su3( su3_matrix *a );
void add_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void sub_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void scalar_mult_su3_matrix( su3_matrix *src, float scalar, su3_matrix *dest);
void scalar_mult_add_su3_matrix( su3_matrix *src1, su3_matrix *src2,
	float scalar, su3_matrix *dest);
void scalar_mult_sub_su3_matrix( su3_matrix *src1, su3_matrix *src2,
	float scalar, su3_matrix *dest);
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

void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c );
complex su3_dot( su3_vector *a, su3_vector *b );
float su3_rdot( su3_vector *a, su3_vector *b );
float magsq_su3vec( su3_vector *a );
void su3vec_copy( su3_vector *a, su3_vector *b );
void dumpvec( su3_vector *v );
void clearvec( su3_vector *v );

void mult_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_su3_mat_vec_sum(  su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_su3_mat_vec_sum_4dir( su3_matrix *a, su3_vector *b0,
	su3_vector *b1, su3_vector *b2, su3_vector *b3, su3_vector *c );
void mult_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_vec_4dir( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_4vec( su3_matrix *mat, su3_vector *src,
			    su3_vector *dest0, su3_vector *dest1, 
			    su3_vector *dest2, su3_vector *dest3  ) ;
void mult_adj_su3_mat_vec_sum( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c );

void add_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c );
void sub_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c );
void sub_four_su3_vecs( su3_vector *a, su3_vector *b1, su3_vector *b2,
	su3_vector *b3, su3_vector *b4 );

void scalar_mult_su3_vector(  su3_vector *src, float scalar, 
	su3_vector *dest);
void scalar_mult_add_su3_vector( su3_vector *src1, su3_vector *src2,
	float scalar, su3_vector *dest);
void scalar_mult_sum_su3_vector( su3_vector *src1, su3_vector *src2,
	float scalar);
void scalar_mult_sub_su3_vector( su3_vector *src1, su3_vector *src2,
	float scalar, su3_vector *dest);
void scalar_mult_wvec( wilson_vector *src, float s, wilson_vector *dest );
void scalar_mult_hwvec( half_wilson_vector *src, float s, 
    half_wilson_vector *dest );
void scalar_mult_add_wvec( wilson_vector *src1, wilson_vector *src2,
	float scalar, wilson_vector *dest );
void scalar_mult_addtm_wvec( wilson_vector *src1, wilson_vector *src2,
	float scalar, wilson_vector *dest );
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
void mult_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
	half_wilson_vector *dest );
void mult_adj_mat_wilson_vec( su3_matrix *mat, wilson_vector *src,
	wilson_vector *dest);
void mult_adj_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
	half_wilson_vector *dest );

void add_wilson_vector( wilson_vector *src1, wilson_vector *src2,
	wilson_vector *dest );
void sub_wilson_vector( wilson_vector *src1, wilson_vector *src2,
       wilson_vector *dest );
float magsq_wvec( wilson_vector *src );
complex wvec_dot( wilson_vector *src1, wilson_vector *src2 );
complex wvec2_dot( wilson_vector *src1, wilson_vector *src2 );
float wvec_rdot( wilson_vector *a, wilson_vector *b );

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

float gaussian_rand_no( double_prn *prn_pt );
# 1 "../include/../include/int32type.h" 1
 



 










 














typedef int int32type;
typedef unsigned int u_int32type;




# 485 "../include/su3.h" 2

void byterevn(int32type w[], int n);

# 21 "prefetch32.c" 2

# 1 "../include/prefetch.h" 1



 
 
 
 
 


 






# 90 "../include/prefetch.h"


# 115 "../include/prefetch.h"


 
 
 
 
 

void _prefetch_M( su3_matrix * );
void _prefetch_V( su3_vector * );
void _prefetch_W( wilson_vector *);
void _prefetch_H( half_wilson_vector *);
void _prefetch_VV( su3_vector *, su3_vector *);
void _prefetch_VVV( su3_vector *, su3_vector *, su3_vector *);
void _prefetch_VVVV( su3_vector *, su3_vector *, su3_vector *, su3_vector *);
void _prefetch_VVVVV( su3_vector *, su3_vector *, su3_vector *, su3_vector *, su3_vector *);
void _prefetch_WWW( wilson_vector *, wilson_vector *, wilson_vector *);
void _prefetch_WWWW( wilson_vector *, wilson_vector *, wilson_vector *, wilson_vector *);
void _prefetch_WWWWW( wilson_vector *, wilson_vector *, wilson_vector *, 
		      wilson_vector *, wilson_vector *);
void _prefetch_4MVVVV( su3_matrix *, su3_vector *, su3_vector *, su3_vector *, su3_vector *);
void _prefetch_4MWWWW( su3_matrix *, wilson_vector *, wilson_vector *, wilson_vector *, wilson_vector *);
void _prefetch_4MV4V( su3_matrix *, su3_vector *, su3_vector *);
void _prefetch_4MW4W( su3_matrix *, wilson_vector *, wilson_vector *);























# 22 "prefetch32.c" 2


 
 

  






  





 






 







 







 


# 76 "prefetch32.c"

 


# 93 "prefetch32.c"


 
void _prefetch_M( su3_matrix *a0 ){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) +  8 ); dummy = *((float *)( a0 ) + 16 ); ;

}

void _prefetch_V( su3_vector *a0 ){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) + 4  ); ;

}

void _prefetch_W( wilson_vector *a0 ){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) +  8 ); dummy = *((float *)( a0 ) + 16 ); dummy = *((float *)( a0 ) + 22 ); ;

}

void _prefetch_H( half_wilson_vector *a0 ){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) +  8 ); dummy = *((float *)( a0 ) + 10 ); ;

}

void _prefetch_VV( su3_vector *a0, su3_vector *a1){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) + 4  ); ;
  dummy = *((float *)( a1 )      ); dummy = *((float *)( a1 ) + 4  ); ;

}

void _prefetch_VVV( su3_vector *a0, su3_vector *a1, su3_vector *a2){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) + 4  ); ;
  dummy = *((float *)( a1 )      ); dummy = *((float *)( a1 ) + 4  ); ;
  dummy = *((float *)( a2 )      ); dummy = *((float *)( a2 ) + 4  ); ;

}

void _prefetch_VVVV( su3_vector *a0, su3_vector *a1, su3_vector *a2, 
		    su3_vector *a3){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) + 4  ); ;
  dummy = *((float *)( a1 )      ); dummy = *((float *)( a1 ) + 4  ); ;
  dummy = *((float *)( a2 )      ); dummy = *((float *)( a2 ) + 4  ); ;
  dummy = *((float *)( a3 )      ); dummy = *((float *)( a3 ) + 4  ); ;

}
void _prefetch_VVVVV( su3_vector *a0, su3_vector *a1, su3_vector *a2, 
		     su3_vector *a3, su3_vector *a4){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) + 4  ); ;
  dummy = *((float *)( a1 )      ); dummy = *((float *)( a1 ) + 4  ); ;
  dummy = *((float *)( a2 )      ); dummy = *((float *)( a2 ) + 4  ); ;
  dummy = *((float *)( a3 )      ); dummy = *((float *)( a3 ) + 4  ); ;
  dummy = *((float *)( a4 )      ); dummy = *((float *)( a4 ) + 4  ); ;

}

void _prefetch_WWW( wilson_vector *a0, wilson_vector *a1, wilson_vector *a2){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) +  8 ); dummy = *((float *)( a0 ) + 16 ); dummy = *((float *)( a0 ) + 22 ); ;
  dummy = *((float *)( a1 )      ); dummy = *((float *)( a1 ) +  8 ); dummy = *((float *)( a1 ) + 16 ); dummy = *((float *)( a1 ) + 22 ); ;
  dummy = *((float *)( a2 )      ); dummy = *((float *)( a2 ) +  8 ); dummy = *((float *)( a2 ) + 16 ); dummy = *((float *)( a2 ) + 22 ); ;

}

void _prefetch_WWWW( wilson_vector *a0, wilson_vector *a1, 
		    wilson_vector *a2, wilson_vector *a3){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) +  8 ); dummy = *((float *)( a0 ) + 16 ); dummy = *((float *)( a0 ) + 22 ); ;
  dummy = *((float *)( a1 )      ); dummy = *((float *)( a1 ) +  8 ); dummy = *((float *)( a1 ) + 16 ); dummy = *((float *)( a1 ) + 22 ); ;
  dummy = *((float *)( a2 )      ); dummy = *((float *)( a2 ) +  8 ); dummy = *((float *)( a2 ) + 16 ); dummy = *((float *)( a2 ) + 22 ); ;
  dummy = *((float *)( a3 )      ); dummy = *((float *)( a3 ) +  8 ); dummy = *((float *)( a3 ) + 16 ); dummy = *((float *)( a3 ) + 22 ); ;

}

void _prefetch_WWWWW( su3_vector *a0, su3_vector *a1, su3_vector *a2,
		     su3_vector *a3, su3_vector *a4){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) +  8 ); dummy = *((float *)( a0 ) + 16 ); dummy = *((float *)( a0 ) + 22 ); ;
  dummy = *((float *)( a1 )      ); dummy = *((float *)( a1 ) +  8 ); dummy = *((float *)( a1 ) + 16 ); dummy = *((float *)( a1 ) + 22 ); ;
  dummy = *((float *)( a2 )      ); dummy = *((float *)( a2 ) +  8 ); dummy = *((float *)( a2 ) + 16 ); dummy = *((float *)( a2 ) + 22 ); ;
  dummy = *((float *)( a3 )      ); dummy = *((float *)( a3 ) +  8 ); dummy = *((float *)( a3 ) + 16 ); dummy = *((float *)( a3 ) + 22 ); ;
  dummy = *((float *)( a4 )      ); dummy = *((float *)( a4 ) +  8 ); dummy = *((float *)( a4 ) + 16 ); dummy = *((float *)( a4 ) + 22 ); ;

}

void _prefetch_4MVVVV( su3_matrix *a0, su3_vector *a1, su3_vector *a2, 
		      su3_vector *a3, su3_vector *a4){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) +  8 ); dummy = *((float *)( a0 ) + 16 ); dummy = *((float *)( a0 ) + 24 ); dummy = *((float *)( a0 ) + 32 ); dummy = *((float *)( a0 ) + 40 ); dummy = *((float *)( a0 ) + 48 ); dummy = *((float *)( a0 ) + 56 ); dummy = *((float *)( a0 ) + 64 ); dummy = *((float *)( a0 ) + 70 ); ;
  dummy = *((float *)( a1 )      ); dummy = *((float *)( a1 ) + 4  ); ;
  dummy = *((float *)( a2 )      ); dummy = *((float *)( a2 ) + 4  ); ;
  dummy = *((float *)( a3 )      ); dummy = *((float *)( a3 ) + 4  ); ;
  dummy = *((float *)( a4 )      ); dummy = *((float *)( a4 ) + 4  ); ;

}

void _prefetch_4MWWWW( su3_matrix *a0, wilson_vector *a1, wilson_vector *a2, 
		      wilson_vector *a3, wilson_vector *a4){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) +  8 ); dummy = *((float *)( a0 ) + 16 ); dummy = *((float *)( a0 ) + 24 ); dummy = *((float *)( a0 ) + 32 ); dummy = *((float *)( a0 ) + 40 ); dummy = *((float *)( a0 ) + 48 ); dummy = *((float *)( a0 ) + 56 ); dummy = *((float *)( a0 ) + 64 ); dummy = *((float *)( a0 ) + 70 ); ;
  dummy = *((float *)( a1 )      ); dummy = *((float *)( a1 ) +  8 ); dummy = *((float *)( a1 ) + 16 ); dummy = *((float *)( a1 ) + 22 ); ;
  dummy = *((float *)( a2 )      ); dummy = *((float *)( a2 ) +  8 ); dummy = *((float *)( a2 ) + 16 ); dummy = *((float *)( a2 ) + 22 ); ;
  dummy = *((float *)( a3 )      ); dummy = *((float *)( a3 ) +  8 ); dummy = *((float *)( a3 ) + 16 ); dummy = *((float *)( a3 ) + 22 ); ;
  dummy = *((float *)( a4 )      ); dummy = *((float *)( a4 ) +  8 ); dummy = *((float *)( a4 ) + 16 ); dummy = *((float *)( a4 ) + 22 ); ;

}
void _prefetch_4MV4V( su3_matrix *a0, su3_vector *a1, su3_vector *a2){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) +  8 ); dummy = *((float *)( a0 ) + 16 ); dummy = *((float *)( a0 ) + 24 ); dummy = *((float *)( a0 ) + 32 ); dummy = *((float *)( a0 ) + 40 ); dummy = *((float *)( a0 ) + 48 ); dummy = *((float *)( a0 ) + 56 ); dummy = *((float *)( a0 ) + 64 ); dummy = *((float *)( a0 ) + 70 ); ;
  dummy = *((float *)( a1 )      ); dummy = *((float *)( a1 ) + 4  ); ;
  dummy = *((float *)( a2 )      ); dummy = *((float *)( a2 ) +  8 ); dummy = *((float *)( a2 ) + 16 ); dummy = *((float *)( a2 ) + 22 ); ;

}


void _prefetch_4MW4W( su3_matrix *a0, wilson_vector *a1, wilson_vector *a2){
  register float dummy;

  dummy = *((float *)( a0 )      ); dummy = *((float *)( a0 ) +  8 ); dummy = *((float *)( a0 ) + 16 ); dummy = *((float *)( a0 ) + 24 ); dummy = *((float *)( a0 ) + 32 ); dummy = *((float *)( a0 ) + 40 ); dummy = *((float *)( a0 ) + 48 ); dummy = *((float *)( a0 ) + 56 ); dummy = *((float *)( a0 ) + 64 ); dummy = *((float *)( a0 ) + 70 ); ;
  dummy = *((float *)( a1 )      ); dummy = *((float *)( a1 ) +  8 ); dummy = *((float *)( a1 ) + 16 ); dummy = *((float *)( a1 ) + 22 ); ;
  dummy = *((float *)( a2 )      ); dummy = *((float *)( a2 ) +  8 ); dummy = *((float *)( a2 ) + 16 ); dummy = *((float *)( a2 ) + 24 ); dummy = *((float *)( a2 ) + 32 ); dummy = *((float *)( a2 ) + 40 ); dummy = *((float *)( a2 ) + 48 ); dummy = *((float *)( a2 ) + 56 ); dummy = *((float *)( a2 ) + 64 ); dummy = *((float *)( a2 ) + 72 ); dummy = *((float *)( a2 ) + 80 ); dummy = *((float *)( a2 ) + 88 ); dummy = *((float *)( a2 ) + 94 ); ;

}









  
# 10 "prefetch.c" 2




  
