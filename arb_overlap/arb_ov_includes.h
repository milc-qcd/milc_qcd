/****************** arb_overlap_includes.h ******************************/
/*
*  Include files for the arbitrary overlap fermion applications.
*/

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/dirs.h"
#include <string.h>


#ifdef EIG
#include "jacobi.h"
#endif

/* prototypes for functions in high level code */
int setup_p();
int readin(int prompt);

/*functions for kernel action*/

void setup_offset();
void setup_links(int iflag);
void path_stuff(int iflag, int ioffset, int k, int npath[5],
                int dir[5][24][4], int sign[5][24][4], 
		int length[5][24], Real cc[5][24]);
void build_params(Real mass_0);
void delta0(field_offset src,field_offset dest,  int isign);
void hdelta0(field_offset src,field_offset dest, int ikind);

void path(int *dir,int *sign,int length);

void make_clov1();
void f_mu_nu1(field_offset f_mn, int mu, int nu);
void mult_ldu1(field_offset src,field_offset dest,
               field_offset triang,field_offset diag, int parity);
void g5delta0_field(wilson_vector * src, wilson_vector * dest, int sign);
void block_stout(int ns);
void d_plaqstout(double *ss_plaq,double *st_plaq, int step);
void boundary_flip_stout(int sign);
void reunitarize_stout(int step);


#ifdef FIELD
void hdelta0_field(wilson_vector *src,wilson_vector *dest,  int isign);
void delta0_field(wilson_vector *  src,wilson_vector *  dest,  int isign);
void mult_ldu1_field(wilson_vector *src,wilson_vector *dest,
triangular *triang, diagonal *diag, int parity);
int congrad_multi_field( /* Return value is number of iterations taken */
    wilson_vector *src,   /* type wilson_vector (where source is to be created)*/
    wilson_vector **psim, /* elements of the multimass inverter */
    int MaxCG,          /* maximum number of iterations per restart */
    Real RsdCG,        /* desired residual - 
                           normalized as sqrt(r*r)/sqrt(src_e*src_e */
    Real *size_r,      /* resulting residual */
    int start_flag     /* 0: use a zero initial guess; 1: use dest */
    );
#endif

#ifdef PROJ
void wp_shrink_pl( wilson_vector *src, half_wilson_vector *dest,
        int dir, int sign );
void wp_grow_pl(  half_wilson_vector *src, wilson_vector *dest,
        int dir, int sign );
#endif




/*function for overlap*/

void setup_multi();
void setup_inner();

int congrad_multi( /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (where source is to be created)*/
    wilson_vector **psim, /* elements of the multimass inverter */
    int MaxCG,          /* maximum number of iterations per restart */
    Real RsdCG,        /* desired residual - 
                           normalized as sqrt(r*r)/sqrt(src_e*src_e */
    Real *size_r,      /* resulting residual */
    int start_flag     /* 0: use a zero initial guess; 1: use dest */
    );

void step(
	  field_offset src,  /* type wilson_vector (where source is to be created)*/
	  field_offset dest  /* type wilson_vector (answer and initial guess) */
	  );
void step_field( wilson_vector* src,  wilson_vector* dest );

void grsource_dirac(int parity);

void hoverlap(field_offset src,
              field_offset dest);


int cg_outer(
    wilson_vector* src,   /* type wilson_vector (where source is to be created)*/
    wilson_vector* dest,   /* type wilson_vector (where dest is to be created)*/
    Real mass,
    int chirality, /* chirality of source vector */
    Real resid,
    int invflag /* inverter flag */
    );


int congrad_multi_o( /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (where source is to be created)*/
    wilson_vector **psim, /* elements of the multimass inverter */
    int MaxCG,          /* maximum number of iterations per restart */
    Real *RsdCG,        /* desired residual - array (one for each mass)
                           normalized as sqrt(r*r)/sqrt(src*src) */
    Real *size_r,      /* resulting residual */
    int start_flag     /* 0: use a zero initial guess; 1: use dest */
    );


/* linalg_stuff.c */
void copy_Vector(wilson_vector *src, wilson_vector *res ) ; 
void norm2(wilson_vector *vec, double *norm ); 
void dot_product(wilson_vector *vec1, wilson_vector *vec2, 
		 double_complex *dot) ;
void complex_vec_mult_sub(double_complex *cc, wilson_vector *vec1, 
			  wilson_vector *vec2) ;
void complex_vec_mult_add(double_complex *cc, wilson_vector *vec1, 
			  wilson_vector *vec2) ;
void double_vec_mult(double *a, wilson_vector *vec1, 
		     wilson_vector *vec2) ;
void double_vec_mult_sub(double *rr, wilson_vector *vec1,  
			 wilson_vector *vec2) ;
void double_vec_mult_add(double *rr, wilson_vector *vec1,  
			 wilson_vector *vec2) ;
void dax_p_by(double *a, wilson_vector *vec1, double *b, wilson_vector *vec2) ; 
void vec_plus_double_vec_mult(wilson_vector *vec1, double *a, wilson_vector *vec2) ; 

void normalize(wilson_vector *vec) ;
void project_out(wilson_vector *vec, wilson_vector **vector, int Num);
void RotateBasis(wilson_vector **eigVec, Matrix *V);
void constructArray(wilson_vector **eigVec,wilson_vector **MeigVec,
 Matrix *A,  int *converged) ;


/*functions for eigenvalues*/

void build_hr0( int flag, int prec /* 0--restart with stored vectors, 1, fresh start*/ );
int build_h0();
int build_lowest_chi(int *trial_chirality);
int build_hov(int *trial_chirality, int *nzero);


void read_eigen(wilson_vector **eigVec,int Nmodes, int in_flag,
                char *in_file);
void write_eigen(wilson_vector **eigVec,int Nmodes, int out_flag,
                 char *out_file);
void read_eigenval(double *eigVal,int Nmodes, char *in_file);
void write_eigenval(double *eigVal,int Nmodes, char *out_file);
void Matrix_Vec_mult(wilson_vector *src, wilson_vector *res) ;
void cleanup_Matrix() ;
int Rayleigh_min(wilson_vector *vec,wilson_vector **eigVec,Real Tolerance, 
                 Real RelTol,int Nvecs,int MaxIter,int Restart);

#ifdef PRIMME
#define Kalkreuter Kalkreuter_PRIMME
#else
#define Kalkreuter Kalkreuter_Ritz
#endif

int Kalkreuter_PRIMME(wilson_vector **eigVec, double *eigVal, Real Tolerance, 
               Real RelTol, int Nvecs, int MaxIter, 
               int Restart, int iters, int parity) ;
int Kalkreuter_Ritz(wilson_vector **eigVec, double *eigVal, Real Tolerance, 
               Real RelTol, int Nvecs, int MaxIter, 
               int Restart, int iters, int parity) ;
//void dot_product(wilson_vector *vec1, wilson_vector *vec2,
// double_complex *dot) ;





void grsource(int parity);
void boundary_flip(int sign );


int congrad_multi_field( /* Return value is number of iterations taken */
	    wilson_vector  *src,   /* type wilson_vector (where source is to be created)*/
   wilson_vector **psim, /* elements of the multimass inverter */
   int MaxCG,          /* maximum number of iterations per restart */
   Real RsdCG,        /* desired residual -
                                 normalized as sqrt(r*r)/sqrt(src*src) */
   Real *size_r,      /* resulting residual */
   int start_flag     /* 0: use a zero initial guess; 1: use dest */
);


#ifdef INV
int congrad_xxx(
    field_offset src,   /* type wilson_vector (where source is to be created)*/
    Real mass,
    int source_chirality /* chirality sector for inversion (NOT USED HERE, could be used
elsewhere)  */
    );
#endif
/* misc */

#ifdef PTP
void print_var(char *label, int i);
void conv4d(field_offset src,field_offset dest);
#endif

/* NHYP code */
#ifndef NHYP_DEBUG
void compute_fhb( su3_matrix *Q, Real *f, Real b[3][3], int compute_b );
#else
void compute_fhb( su3_matrix *Omega, su3_matrix *Q,
                  Real *f, Real b[3][3], int compute_b );
#endif
void block_nhyp();

void gauge_field_copy(field_offset src,field_offset dest);
void gauge_field_change(field_offset src,field_offset dest);
void grsource_chiral_field(chiral_src* f);


void refresh_links();
void re_setup_inner(double emin, double emax);
Real vectornorm(wilson_vector * w1);
void project_out_eigen(wilson_vector* vec, wilson_vector** evec, int nvec);


#ifdef HYPSTOUT
void  stout_force1();
void  stout_force2();
void  stout_force3();
#endif

Real magsq_hwvec( half_wilson_vector *vec );

void hw_to_w(half_wilson_vector *src1, int ichiral,
	     wilson_vector *dest);
void w_to_hw(wilson_vector *src1, int ichiral,
	     half_wilson_vector *dest);
