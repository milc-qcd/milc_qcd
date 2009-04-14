/* Temporary include file.
   These routined should eventually be included in su3.h
   Started: 071005
   Last modified: 081121 */


#ifndef SU3_MAT_OP_H_

#define SU3_MAT_OP_H_

/* MILC includes */
#include "../include/complex.h"
#include "../include/su3.h"

#if (PRECISION==1)
#define SU3_UNIT_ANALYTIC_EPS 1.0e-6
#define SU3_ROOT_INV_NORM_EPS 1.0e-6
#else
#define SU3_UNIT_ANALYTIC_EPS 1.0e-14
#define SU3_ROOT_INV_NORM_EPS 1.0e-9
#endif

#define SU3_UNIT_DER_EPS 1.0e-6

// If one runs into small eigenvalues when calculating
// the HISQ fermion force, defining this option allows to
// regularize them and preven large spikes in the force
#define HISQ_FORCE_FILTER 5.0e-5


#define SU3_ROOT_INV_MAX_ITER 100
#define MILC_AB_PI 3.1415926535897931159979634685441852
#define MILC_AB_TPI 6.2831853071795862319959269370883703

#define SU3_UNIT_RAT_NTERMS 14
static Real c_l_SU3_UNIT_RAT[ SU3_UNIT_RAT_NTERMS+1 ] = {
         0.0850910,
         2.413330975,
         6.257184884e-01,
         2.925707925e-01,
         1.737405612e-01,
         1.166359792e-01,
         8.372555094e-02,
         6.216038074e-02,
         4.652496186e-02,
         3.423610040e-02,
         2.404754621e-02,
         1.545550091e-02,
         8.436481876e-03,
         3.419245947e-03,
         1.138166539e-03  };
static Real d_l_SU3_UNIT_RAT[ SU3_UNIT_RAT_NTERMS+1 ] = {
         0.0,
         1.361747338e+01,
         3.135687028e+00,
         1.213113539e+00,
         5.596349298e-01,
         2.752627333e-01,
         1.364115846e-01,
         6.543005714e-02,
         2.923946484e-02,
         1.164228894e-02,
         3.887745892e-03,
         9.937321442e-04,
         1.684882417e-04,
         1.585925699e-05,
         5.914114023e-07  };


/* Make 3x3 matrix unity: dest=I */
void unity_su3mat( su3_matrix *dest );


/* Frobenius (or Hilbert-Schmidt) norm of a matrix A:
   norm(A)=sqrt(Tr{A^+A})=sqrt(sum_ij |a_ij|^2) */
Real su3_norm_frob( su3_matrix *a );


/* Inverse of 3x3 complex matrix: B=A^-1
   For 3x3 matrix brute force calculation using
   determinants is about 2.5 times faster than
   LU factorization in "meschach" library */
void su3_inverse( su3_matrix *a, su3_matrix *b );


/* Square root and inverse square root of 
   3x3 complex matrix A: X=A^1/2, Y=A^-1/2
   Denman-Beavers recursion is used
   (Denman, Beavers, Appl. Math. Comp. 2 (1976) 63-94)
   Let X_0=A, Y_0=I, then
   X_{k+1}= 1/2 ( X_k + Y_k^-1 )
   Y_{k+1}= 1/2 ( Y_k + X_k^-1 )
   X converges to A^1/2, Y converges to A^-1/2
   In cases of interest typical recursion converges in 12 steps
   to give Frobenius norm of (X*X-A) and (A*Y*Y-I) 
   less then 10^-8. Therefore on average this routine
   requires 24 3x3 matrix inversions. SLOW! */
void su3_root_inv( su3_matrix *a, su3_matrix *x, su3_matrix *y);


/* Unitarize 3x3 complex matrix:
   B=A*(A^+ A)^-1/2
   B is a U(3) matrix but not an SU(3) matrix(!) */
void su3_unitarize( su3_matrix *a, su3_matrix *b );


/* Special unitarize a unitary matrix:
   B=A/det(A)^1/3
   A is a U(3) matrix and B is an SU(3) matrix.
   This function also returns det(A) */
void su3_spec_unitarize( su3_matrix *a, su3_matrix *b, complex *detA );


/* Special unitarize a unitary matrix:
   B=A/det(A)^1/3
   A is a U(3) matrix and B is an SU(3) matrix.
   This function also returns det(A) */
void su3_spec_unitarize_index( su3_matrix *a, su3_matrix *b, complex *detA, 
           int index_site, int index_dir );


/* Derivative of the unitarized matrix with respect to the original:
   dW/dU and d(W^+)/dU (at fixed U^+ !), where W=U(U^+U)^-1/2
   EXPLAIN indexing convention!!! */
void su3_unit_der( su3_matrix *u, su3_tensor4 *dwdu, su3_tensor4 *dwdagdu );


/* Multiply tensor4 by su3_matrix on the right:
   C_{ijkl}=A_{ijkm}B_{ml} */
void mult_su3_t4m( su3_tensor4 *a, su3_matrix *b, su3_tensor4 *c );


/* Multiply tensor4 by su3_matrix on the left:
   C_{ijkl}=A_{im}B_{mjkl} */
void mult_su3_mt4( su3_matrix *a, su3_tensor4 *b, su3_tensor4 *c );


/* Add two tensor4: C=A+B */
void add_su3_t4( su3_tensor4 *a, su3_tensor4 *b, su3_tensor4 *c );


/* Subtract tensor4: C=A-B */
void sub_su3_t4( su3_tensor4 *a, su3_tensor4 *b, su3_tensor4 *c );


/* Frobenius norm of tensor4 */
Real su3_t4_norm_frob( su3_tensor4 *a );


/* Dump tensor4 on screen */
void dumptensor4( su3_tensor4 *a );


/* Derivative of the special unitarized matrix with respect to the original:
   dW/dU and d(W^+)/dU (at fixed U^+ !), where W=V/det(V)^1/3, V=U(U^+U)^-1/2
   EXPLAIN indexing convention!!! */
void 
  su3_spec_unit_der( su3_matrix *u, su3_tensor4 *dwdu, su3_tensor4 *dwdagdu );


/* Given a derivative of U(3) matrix calculate a 
   derivative of SU(3) matrix
   U -- original matrix,
   W -- U(3) matrix out of U (could be W=U(U^+U)^-1/2 or Re{Tr{WU^+}}=max,
   V=W/det(W)^1/3
   W^-1=W^+ is used here(!)
   EXPLAIN as in Wong, Woloshyn notes!!! */
void su3_unit_der_spec( su3_matrix *u, su3_matrix *w, su3_matrix *wdag, 
                        su3_tensor4 *dwdu, su3_tensor4 *dwdagdu,
                        su3_tensor4 *dvdu, su3_tensor4 *dvdagdu );


/* Derivative of a unitary matrix W with respect X=Re{U} and Y=Im{U}.
   CAUTION! This routine does not fix U^+ as independet variable.
   Insted X, Y are independent variables, so finite difference
   are calculated independently for real and purely imaginary epsilon. */
void su3_unit_der_reim( su3_matrix *u, 
                        su3_tensor4 *dwdure, su3_tensor4 *dwduim, 
                        su3_tensor4 *dwdagdure, su3_tensor4 *dwdagduim );


/* Out of dW/d(Re{U}) and dW/d(Im{U}) construct dW/dU at fixed U^+,
   same for dW^+/dU */
void su3_unit_der_reim_join( 
        su3_tensor4 *dwdure, su3_tensor4 *dwduim, 
        su3_tensor4 *dwdagdure, su3_tensor4 *dwdagduim,
        su3_tensor4 *dwdu, su3_tensor4 *dwdagdu );


/* Unitarization with rational approximation */
void su3_unitarize_rational( su3_matrix *V, su3_matrix *W );


/* Derivative of a unitarized matrix with rational approximation */
void su3_unit_der_rational( su3_matrix *V, 
                            su3_tensor4 *dwdv, su3_tensor4 *dwdagdu );


/* Derivative dW/dY, dW^+/dY, where Y is U(3) matrix,
   W=Y/(detY)^1/3 */
void su3_der_detWY( su3_matrix *y, su3_tensor4 *dwdy, su3_tensor4 *dwdagdy );



/* Analytic unitarization, Hasenfratz, Hoffmann, Schaefer, JHEP05 (2007) 029 */
void su3_unitarize_analytic( su3_matrix *V, su3_matrix *W );



/* Analytic unitarization, Hasenfratz, Hoffmann, Schaefer, JHEP05 (2007) 029 */
void su3_unitarize_analytic_index( su3_matrix *V, su3_matrix *W, int index_site, int index_dir );



/* Analytic derivative of the unitarized matrix with respect to the original:
   dW/dV and d(W^+)/dV, where W=V(V^+V)^-1/2
   EXPLAIN indexing convention!!! */
void su3_unit_der_analytic( su3_matrix *V, 
               su3_tensor4 *dwdv, su3_tensor4 *dwdagdv );



/* copy rank 4 tensor: a -> b */
void su3t4_copy( su3_tensor4 *a, su3_tensor4 *b );



/* **************************************************
   SVD stuff needs to be put into a separate file
   ************************************************** */
/* Input: A -- 3x3 complex matrix,
   Output: sigma[3] -- singular values,
           U,V -- U(3) matrices such, that
           A=U Sigma V^+ */
int svd3x3(double A[3][3][2], double *sigma, double U[3][3][2], double V[3][3][2]);


#endif /* SU3_MAT_OP_H_ */
