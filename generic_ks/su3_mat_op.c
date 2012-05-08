/****************  su3_mat_op.c *************************/
/* MIMD version 7 */
/* Various utilities for 3x3 complex matrix routines
   (arbitrary, not necessarily SU(3) matrices).
   AB 071005 */


#include <stdio.h>
#include <math.h>
#include "generic_ks_includes.h"      /* definitions files and prototypes */
#include "../include/su3_mat_op.h"
#include "../include/info.h"


// this is needed for
// u3_unitarize_analytic and u3_unit_der_analytic
#if(PRECISION==2)
#define U3_UNIT_ANALYTIC_FOLLOW_PREC
#endif

#ifndef HISQ_REUNIT_SVD_REL_ERROR
#define HISQ_REUNIT_SVD_REL_ERROR 1e-8
#endif

#ifndef HISQ_REUNIT_SVD_ABS_ERROR
#define HISQ_REUNIT_SVD_ABS_ERROR 1e-8
#endif

/* Frobenius (or Hilbert-Schmidt) norm of a matrix A:
   norm(A)=sqrt(Tr{A^+A})=sqrt(sum_ij |a_ij|^2) */
Real su3_norm_frob( su3_matrix *a ) {
  Real sum;

  sum = cabs_sq( &(a->e[0][0]) ) + cabs_sq( &(a->e[0][1]) ) 
      + cabs_sq( &(a->e[0][2]) ) + cabs_sq( &(a->e[1][0]) )
      + cabs_sq( &(a->e[1][1]) ) + cabs_sq( &(a->e[1][2]) )
      + cabs_sq( &(a->e[2][0]) ) + cabs_sq( &(a->e[2][1]) )
      + cabs_sq( &(a->e[2][2]) );

  return sqrt( sum );
}


/* Make 3x3 matrix unity: dest=I */
void unity_su3mat( su3_matrix *dest ) {
  register int i,j;
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    dest->e[i][j].real = dest->e[i][j].imag = 0.0;
  }
  dest->e[0][0].real = dest->e[1][1].real = dest->e[2][2].real = 1.0;
}


/* Inverse of 3x3 complex matrix: B=A^-1
   For 3x3 matrix brute force calculation using
   determinants is about 2.5 times faster than
   LU factorization in "meschach" library */
void su3_inverse( su3_matrix *a, su3_matrix *b ) {
  complex det, temp1, temp2, temp3;

  det = det_su3( a );

  if( cabs_sq( &det ) < 1.0e-8 ) {
    fprintf( stdout, "*** WARNING in su3_inverse: determinant is too small\n");
    fprintf( stdout, "det=%18.10g\n", cabs_sq( &det ) );
  }

  /* ** rely on fast macros, avoid function call overhead ** */
  /* b_00 */
  CMUL( a->e[1][1], a->e[2][2], temp1 );
  CMUL( a->e[1][2], a->e[2][1], temp2 );
  CSUB( temp1, temp2, temp3 );
  CDIV( temp3, det, b->e[0][0] );
  /* b_01 */
  CMUL( a->e[0][2], a->e[2][1], temp1 );
  CMUL( a->e[0][1], a->e[2][2], temp2 );
  CSUB( temp1, temp2, temp3 );
  CDIV( temp3, det, b->e[0][1] );
  /* b_02 */
  CMUL( a->e[0][1], a->e[1][2], temp1 );
  CMUL( a->e[0][2], a->e[1][1], temp2 );
  CSUB( temp1, temp2, temp3 );
  CDIV( temp3, det, b->e[0][2] );
  /* b_10 */
  CMUL( a->e[1][2], a->e[2][0], temp1 );
  CMUL( a->e[1][0], a->e[2][2], temp2 );
  CSUB( temp1, temp2, temp3 );
  CDIV( temp3, det, b->e[1][0] );
  /* b_11 */
  CMUL( a->e[0][0], a->e[2][2], temp1 );
  CMUL( a->e[0][2], a->e[2][0], temp2 );
  CSUB( temp1, temp2, temp3 );
  CDIV( temp3, det, b->e[1][1] );
  /* b_12 */
  CMUL( a->e[0][2], a->e[1][0], temp1 );
  CMUL( a->e[0][0], a->e[1][2], temp2 );
  CSUB( temp1, temp2, temp3 );
  CDIV( temp3, det, b->e[1][2] );
  /* b_20 */
  CMUL( a->e[1][0], a->e[2][1], temp1 );
  CMUL( a->e[1][1], a->e[2][0], temp2 );
  CSUB( temp1, temp2, temp3 );
  CDIV( temp3, det, b->e[2][0] );
  /* b_21 */
  CMUL( a->e[0][1], a->e[2][0], temp1 );
  CMUL( a->e[0][0], a->e[2][1], temp2 );
  CSUB( temp1, temp2, temp3 );
  CDIV( temp3, det, b->e[2][1] );
  /* b_22 */
  CMUL( a->e[0][0], a->e[1][1], temp1 );
  CMUL( a->e[0][1], a->e[1][0], temp2 );
  CSUB( temp1, temp2, temp3 );
  CDIV( temp3, det, b->e[2][2] );
}


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
void u3_root_inv( su3_matrix *a, su3_matrix *x, su3_matrix *y) {
  su3_matrix xinv, yinv, x2, y2, diff, Unit;
  Real norm_x, norm_y;
  int iter;

  /* 0th iteration: x=a, y=I */
  su3mat_copy( a, x );
  unity_su3mat( y );
  unity_su3mat( &Unit );

  iter=0;
  do {
    /* get inverses */
    su3_inverse( x, &xinv );
    su3_inverse( y, &yinv );

    /* (k+1) step */
    add_su3_matrix( x, &yinv, &x2 );
    add_su3_matrix( y, &xinv, &y2 );
    scalar_mult_su3_matrix( &x2, 0.5, x );
    scalar_mult_su3_matrix( &y2, 0.5, y );

    /* calculate norm of (x^1/2) * (x^1/2) - x */
    mult_su3_nn( x, x, &x2 );
    sub_su3_matrix( &x2, a, &diff );
    norm_x = su3_norm_frob( &diff );

    /* calculate norm of a * (x^-1/2) * (x^-1/2) - I */
    mult_su3_nn( a, y, &x2 );
    mult_su3_nn( &x2, y, &y2 );
    sub_su3_matrix( &y2, &Unit, &diff );
    norm_y = su3_norm_frob( &diff );

    iter++;

  } while( (iter<U3_ROOT_INV_MAX_ITER) && 
           (norm_x>U3_ROOT_INV_NORM_EPS) &&
           (norm_y>U3_ROOT_INV_NORM_EPS) );

}


/* Unitarize 3x3 complex matrix:
   B=A*(A^+ A)^-1/2
   B is a U(3) matrix but not an SU(3) matrix(!) */
void u3_unitarize( su3_matrix *a, su3_matrix *b ) {
  su3_matrix a2, x, y;

  /* X=A^+ A */
  mult_su3_an( a, a, &a2 );

  /* X=A2^1/2, Y=A2^-1/2 */
  u3_root_inv( &a2, &x, &y );

  /* B=A*Y */
  mult_su3_nn( a, &y, b );
}


/* Special unitarize a unitary matrix:
   B=A/det(A)^1/3
   A is a U(3) matrix and B is an SU(3) matrix.
   This function also returns det(A) */
void su3_spec_unitarize( su3_matrix *a, su3_matrix *b, complex *detA ) {
  complex s;
  Real argdet, r;

  (*detA) = det_su3( a );

  /* take third root, choose argument in [-pi/3, pi/3) branch */
  argdet = carg( detA );
  r = cabs( detA );
  r = exp( log(r)/3 );
  argdet /= 3;

  /* multiply the matrix by inverse root */
  s.real = cos( argdet ) / r;
  s.imag = -sin( argdet ) / r;
  c_scalar_mult_su3mat( a, &s, b );
}


/* Special unitarize a unitary matrix:
   B=A/det(A)^1/3
   A is a U(3) matrix and B is an SU(3) matrix.
   This function also returns det(A) */
void su3_spec_unitarize_index( su3_matrix *a, su3_matrix *b, complex *detA, 
           int index_site, int index_dir ) {
  complex s;
  Real argdet, r;

  (*detA) = det_su3( a );

  /* take third root, choose argument in [-pi/3, pi/3) branch */
  argdet = carg( detA );
  r = cabs( detA );
  r = exp( log(r)/3 );
  argdet /= 3;

  /* multiply the matrix by inverse root */
  s.real = cos( argdet ) / r;
  s.imag = -sin( argdet ) / r;
  c_scalar_mult_su3mat( a, &s, b );
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
   /* store SU(3) matrices */
  if( lattice[index_site].on_step_W[index_dir] < global_current_time_step ) {
    lattice[index_site].on_step_W[index_dir] = global_current_time_step;
    su3mat_copy( &(lattice[index_site].Wlink[index_dir]),
                 &(lattice[index_site].Wlink_previous[index_dir]) );
    su3mat_copy( b, &(lattice[index_site].Wlink[index_dir]) );
  }
#endif /* HISQ_REUNITARIZATION_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
}


/* Derivative of the unitarized matrix with respect to the original:
   dW/dU and d(W^+)/dU (at fixed U^+ !), where W=U(U^+U)^-1/2 */
void u3_unit_der( su3_matrix *u, su3_tensor4 *dwdu, su3_tensor4 *dwdagdu ) {
  su3_matrix up, um, a, b, up12, um12, upu, umu, dw, dwdag;
  int i, j, m, n;
  Real factor;

  factor = 0.5/U3_UNIT_DER_EPS;

  /* loop on components of the original matrix U */
  for( i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      /* temporary storage */
      su3mat_copy( u, &up );
      su3mat_copy( u, &um );

      /* add/subtract epsilon to the real part of [i][j] component */
      (up.e[i][j]).real += U3_UNIT_DER_EPS;
      (um.e[i][j]).real -= U3_UNIT_DER_EPS;

      /* Up12=(U^+ Up)^-1/2 */
      mult_su3_an( u, &up, &a );
      u3_root_inv( &a, &b, &up12 );

      /* Um12=(U^+ Um)^-1/2 */
      mult_su3_an( u, &um, &a );
      u3_root_inv( &a, &b, &um12 );

      /* ** dW/dU ** */
      /* Upu=Up * (U^+ Up)^-1/2 */
      mult_su3_nn( &up, &up12, &upu );

      /* Umu=Um * (U^+ Um)^-1/2 */
      mult_su3_nn( &um, &um12, &umu );

      sub_su3_matrix( &upu, &umu, &dw );

      /* ** dW^+/dU ** */
      mult_su3_na( &up12, u, &upu );
      mult_su3_na( &um12, u, &umu );

      sub_su3_matrix( &upu, &umu, &dwdag );

      /* populate tensor4 */
      for( m=0; m<3; m++) {
        for( n=0; n<3; n++) {
          dwdu->t4[m][i][j][n].real = factor * (dw.e[m][n]).real;
          dwdu->t4[m][i][j][n].imag = factor * (dw.e[m][n]).imag;
          dwdagdu->t4[m][i][j][n].real = factor * (dwdag.e[m][n]).real;
          dwdagdu->t4[m][i][j][n].imag = factor * (dwdag.e[m][n]).imag;
        }
      }
    }
  }
}


/* Derivative of the special unitarized matrix with respect to the original:
   dW/dU and d(W^+)/dU (at fixed U^+ !), where W=V/det(V)^1/3, V=U(U^+U)^-1/2.
   THIS IS LESS STABLE THAN CALCULATING THE DERIVATIVE OF U(3) MATRIX AND
   THEN USING ANALYTIC EXPRESSION FOR DERIVATIVE OF W=V/det(V)^1/3 */
void su3_spec_unit_der( su3_matrix *u, su3_tensor4 *dwdu, su3_tensor4 *dwdagdu ) {
  su3_matrix up, um, a, b, up12, um12, upu, umu, dw, dwdag, upu2, umu2;
  int i, j, m, n;
  Real factor;
  complex det;

  factor = 0.5/U3_UNIT_DER_EPS;

  /* loop on components of the original matrix U */
  for( i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      /* temporary storage */
      su3mat_copy( u, &up );
      su3mat_copy( u, &um );

      /* add/subtract epsilon to the real part of [i][j] component */
      (up.e[i][j]).real += U3_UNIT_DER_EPS;
      (um.e[i][j]).real -= U3_UNIT_DER_EPS;

      /* Up12=(U^+ Up)^-1/2 */
      mult_su3_an( u, &up, &a );
      u3_root_inv( &a, &b, &up12 );

      /* Um12=(U^+ Um)^-1/2 */
      mult_su3_an( u, &um, &a );
      u3_root_inv( &a, &b, &um12 );

      /* ** dW/dU ** */
      /* Upu=Up * (U^+ Up)^-1/2 */
      mult_su3_nn( &up, &up12, &upu2 );

      /* Umu=Um * (U^+ Um)^-1/2 */
      mult_su3_nn( &um, &um12, &umu2 );

      su3_spec_unitarize( &upu2, &upu, &det );
      su3_spec_unitarize( &umu2, &umu, &det );

      sub_su3_matrix( &upu, &umu, &dw );

      /* ** dW^+/dU ** */
      mult_su3_na( &up12, u, &upu2 );
      mult_su3_na( &um12, u, &umu2 );

      su3_spec_unitarize( &upu2, &upu, &det );
      su3_spec_unitarize( &umu2, &umu, &det );

      sub_su3_matrix( &upu, &umu, &dwdag );

      /* populate tensor4 */
      for( m=0; m<3; m++) {
        for( n=0; n<3; n++) {
          dwdu->t4[m][i][j][n].real = factor * (dw.e[m][n]).real;
          dwdu->t4[m][i][j][n].imag = factor * (dw.e[m][n]).imag;
          dwdagdu->t4[m][i][j][n].real = factor * (dwdag.e[m][n]).real;
          dwdagdu->t4[m][i][j][n].imag = factor * (dwdag.e[m][n]).imag;
        }
      }
    }
  }
}


/* Multiply tensor4 by su3_matrix on the right:
   C_{ijkl}=A_{ijkm}B_{ml} */
void mult_su3_t4m( su3_tensor4 *a, su3_matrix *b, su3_tensor4 *c ) {
  int i, j, k, l, m;

  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      for( k=0; k<3; k++) {
        for( l=0; l<3; l++) {
          (c->t4[i][j][k][l]).real = 0.0;
          (c->t4[i][j][k][l]).imag = 0.0;
          for( m=0; m<3; m++) {
            /* real part */
            (c->t4[i][j][k][l]).real +=
               (a->t4[i][j][k][m]).real * (b->e[m][l]).real
             - (a->t4[i][j][k][m]).imag * (b->e[m][l]).imag;
            /* imag part */
            (c->t4[i][j][k][l]).imag +=
               (a->t4[i][j][k][m]).imag * (b->e[m][l]).real
             + (a->t4[i][j][k][m]).real * (b->e[m][l]).imag;
          }
        }
      }
    }
  }
}


/* Multiply tensor4 by su3_matrix on the left:
   C_{ijkl}=A_{im}B_{mjkl} */
void mult_su3_mt4( su3_matrix *a, su3_tensor4 *b, su3_tensor4 *c ) {
  int i, j, k, l, m;

  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      for( k=0; k<3; k++) {
        for( l=0; l<3; l++) {
          (c->t4[i][j][k][l]).real = 0.0;
          (c->t4[i][j][k][l]).imag = 0.0;
          for( m=0; m<3; m++) {
            /* real part */
            (c->t4[i][j][k][l]).real +=
               (a->e[i][m]).real * (b->t4[m][j][k][l]).real
             - (a->e[i][m]).imag * (b->t4[m][j][k][l]).imag;
            /* imag part */
            (c->t4[i][j][k][l]).imag +=
               (a->e[i][m]).imag * (b->t4[m][j][k][l]).real
             + (a->e[i][m]).real * (b->t4[m][j][k][l]).imag;
          }
        }
      }
    }
  }
}


/* Add two tensor4: C=A+B */
void add_su3_t4( su3_tensor4 *a, su3_tensor4 *b, su3_tensor4 *c ) {
  int i, j, k, l;

  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      for( k=0; k<3; k++) {
        for( l=0; l<3; l++) {
          (c->t4[i][j][k][l]).real = 
              (a->t4[i][j][k][l]).real + (b->t4[i][j][k][l]).real;
          (c->t4[i][j][k][l]).imag = 
              (a->t4[i][j][k][l]).imag + (b->t4[i][j][k][l]).imag;
        }
      }
    }
  }
}


/* Subtract tensor4: C=A-B */
void sub_su3_t4( su3_tensor4 *a, su3_tensor4 *b, su3_tensor4 *c ) {
  int i, j, k, l;

  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      for( k=0; k<3; k++) {
        for( l=0; l<3; l++) {
          (c->t4[i][j][k][l]).real = 
              (a->t4[i][j][k][l]).real - (b->t4[i][j][k][l]).real;
          (c->t4[i][j][k][l]).imag = 
              (a->t4[i][j][k][l]).imag - (b->t4[i][j][k][l]).imag;
        }
      }
    }
  }
}


/* Frobenius norm of tensor4 */
Real su3_t4_norm_frob( su3_tensor4 *a ) {
  int i, j, k, l;
  Real sum;

  sum = 0;
  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      for( k=0; k<3; k++) {
        for( l=0; l<3; l++) {
          sum += (a->t4[i][j][k][l]).real * (a->t4[i][j][k][l]).real
               + (a->t4[i][j][k][l]).imag * (a->t4[i][j][k][l]).imag;
        }
      }
    }
  }

  return sqrt( sum );
}


/* Given a derivative of U(3) matrix calculate a 
   derivative of SU(3) matrix
   U -- original matrix,
   W -- U(3) matrix out of U (could be W=U(U^+U)^-1/2 or Re{Tr{WU^+}}=max,
   V=W/det(W)^1/3
   W^-1=W^+ is assumed */
void su3_unit_der_spec( su3_matrix *u, su3_matrix *w, su3_matrix *wdag, 
                        su3_tensor4 *dwdu, su3_tensor4 *dwdagdu,
                        su3_tensor4 *dvdu, su3_tensor4 *dvdagdu ) {
  int i, j, k, l, m;
  su3_tensor4 t1, t2, t3, t4;
  su3_matrix a, b;
  Real factor, argdet, r;
  complex det, cf1, cf2;

  /* calculate W^-1*dW/dU */
  mult_su3_mt4( wdag, dwdu, &t1 );
  /* calculate (W^+)^-1*dW^+/dU */
  mult_su3_mt4( w, dwdagdu, &t2 );

  /* take Tr{W^-1*dW/dU} and Tr{W*dW^+/dU} on outer indeces */
  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      (a.e[i][j]).real = 0.0;
      (a.e[i][j]).imag = 0.0;
      (b.e[i][j]).real = 0.0;
      (b.e[i][j]).imag = 0.0;
      for( m=0; m<3; m++) {
        (a.e[i][j]).real += (t1.t4[m][i][j][m]).real;
        (a.e[i][j]).imag += (t1.t4[m][i][j][m]).imag;
        (b.e[i][j]).real += (t2.t4[m][i][j][m]).real;
        (b.e[i][j]).imag += (t2.t4[m][i][j][m]).imag;
      }
    }
  }

  /* make direct product of Tr{...} and W (and W^+) and multiply by -1/3 */
  factor = -1.0/3;
  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      for( k=0; k<3; k++) {
        for( l=0; l<3; l++) {
          CMUL( w->e[i][l], a.e[j][k], cf1 );
          (t3.t4[i][j][k][l]).real = factor * cf1.real;
          (t3.t4[i][j][k][l]).imag = factor * cf1.imag;
          CMUL( wdag->e[i][l], b.e[j][k], cf2 );
          (t4.t4[i][j][k][l]).real = factor * cf2.real;
          (t4.t4[i][j][k][l]).imag = factor * cf2.imag;
        }
      }
    }
  }

  /* add -1/3Tr{..}W and dW/dU, same for W^+ */
  add_su3_t4( &t3, dwdu, &t1 );
  add_su3_t4( &t4, dwdagdu, &t2 );

  /* get determinant of W */
  det = det_su3( w );
  argdet = carg( &det )/3;
  r = cabs_sq( &det );
  r = exp( log(r)/6 );
  factor = cos( argdet )/r;
  cf1.real = factor;
  cf2.real = factor;
  factor = sin( argdet )/r;
  cf1.imag = -factor;
  cf2.imag = factor;

  /* multiply by det(W)^-1/3, same for W^+ */
  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      for( k=0; k<3; k++) {
        for( l=0; l<3; l++) {
          CMUL( t1.t4[i][j][k][l], cf1, dvdu->t4[i][j][k][l] );
          CMUL( t2.t4[i][j][k][l], cf2, dvdagdu->t4[i][j][k][l] );
        }
      }
    }
  }

}


/* Derivative of a unitary matrix W with respect X=Re{U} and Y=Im{U}.
   CAUTION! This routine does not fix U^+ as independet variable.
   Insted X, Y are independent variables, so finite difference
   are calculated independently for real and purely imaginary epsilon. */
void su3_unit_der_reim( su3_matrix *u, 
                        su3_tensor4 *dwdure, su3_tensor4 *dwduim, 
                        su3_tensor4 *dwdagdure, su3_tensor4 *dwdagduim ) {
  su3_matrix upre, upim, umre, umim, a, b, up12re, up12im, um12re, um12im;
  su3_matrix upure, upuim, umure, umuim, dwre, dwim, dwdagre, dwdagim;
  int i, j, m, n;
  Real factor;

  factor = 0.5/U3_UNIT_DER_EPS;

  /* loop on components of the original matrix U */
  for( i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      /* temporary storage
         THIS CAN BE SPEEDED UP LATER */
      su3mat_copy( u, &upre );
      su3mat_copy( u, &umre );
      su3mat_copy( u, &upim );
      su3mat_copy( u, &umim );

      /* add/subtract epsilon to the real part of [i][j] component */
      (upre.e[i][j]).real += U3_UNIT_DER_EPS;
      (umre.e[i][j]).real -= U3_UNIT_DER_EPS;
      (upim.e[i][j]).imag += U3_UNIT_DER_EPS;
      (umim.e[i][j]).imag -= U3_UNIT_DER_EPS;

      /* Up12=(Up^+ Up)^-1/2 */
      mult_su3_an( &upre, &upre, &a );
      u3_root_inv( &a, &b, &up12re );
      mult_su3_an( &upim, &upim, &a );
      u3_root_inv( &a, &b, &up12im );

      /* Um12=(Um^+ Um)^-1/2 */
      mult_su3_an( &umre, &umre, &a );
      u3_root_inv( &a, &b, &um12re );
      mult_su3_an( &umim, &umim, &a );
      u3_root_inv( &a, &b, &um12im );

      /* ** dW/dU ** */
      /* Upu=Up * (Up^+ Up)^-1/2 */
      mult_su3_nn( &upre, &up12re, &upure );
      mult_su3_nn( &upim, &up12im, &upuim );

      /* Umu=Um * (Um^+ Um)^-1/2 */
      mult_su3_nn( &umre, &um12re, &umure );
      mult_su3_nn( &umim, &um12im, &umuim );

      sub_su3_matrix( &upure, &umure, &dwre );
      sub_su3_matrix( &upuim, &umuim, &dwim );

      /* ** dW^+/dU ** */
      mult_su3_na( &up12re, &upre, &upure );
      mult_su3_na( &up12im, &upim, &upuim );
      mult_su3_na( &um12re, &umre, &umure );
      mult_su3_na( &um12im, &umim, &umuim );

      sub_su3_matrix( &upure, &umure, &dwdagre );
      sub_su3_matrix( &upuim, &umuim, &dwdagim );

      /* populate tensor4 */
      for( m=0; m<3; m++) {
        for( n=0; n<3; n++) {
          dwdure->t4[m][i][j][n].real = factor * (dwre.e[m][n]).real;
          dwdure->t4[m][i][j][n].imag = factor * (dwre.e[m][n]).imag;
          dwduim->t4[m][i][j][n].real = factor * (dwim.e[m][n]).real;
          dwduim->t4[m][i][j][n].imag = factor * (dwim.e[m][n]).imag;
          dwdagdure->t4[m][i][j][n].real = factor * (dwdagre.e[m][n]).real;
          dwdagdure->t4[m][i][j][n].imag = factor * (dwdagre.e[m][n]).imag;
          dwdagduim->t4[m][i][j][n].real = factor * (dwdagim.e[m][n]).real;
          dwdagduim->t4[m][i][j][n].imag = factor * (dwdagim.e[m][n]).imag;
        }
      }
    }
  }
}


/* Out of dW/d(Re{U}) and dW/d(Im{U}) construct dW/dU at fixed U^+,
   same for dW^+/dU */
void su3_unit_der_reim_join( 
        su3_tensor4 *dwdure, su3_tensor4 *dwduim, 
        su3_tensor4 *dwdagdure, su3_tensor4 *dwdagduim,
        su3_tensor4 *dwdu, su3_tensor4 *dwdagdu ) {
  int i, j, k, l;

  /* for a complex function f(z), z=x+iy:
     df/dz=1/2*(df/dx-i*df/dy) */
  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      for( k=0; k<3; k++) {
        for( l=0; l<3; l++) {
          dwdu->t4[i][j][k][l].real = 
            0.5 * ( dwdure->t4[i][j][k][l].real + dwduim->t4[i][j][k][l].imag
//                  + dwdagdure->t4[i][j][k][l].real - dwdagduim->t4[i][j][k][l].imag 
                  );
          dwdu->t4[i][j][k][l].imag = 
            0.5 * ( dwdure->t4[i][j][k][l].imag - dwduim->t4[i][j][k][l].real
//                  + dwdagdure->t4[i][j][k][l].imag - dwdagduim->t4[i][j][k][l].real 
                  );
          dwdagdu->t4[i][j][k][l].real = 
//            0.5 * ( dwdagdure->t4[i][j][k][l].real + dwdagduim->t4[i][j][k][l].imag );
            0.5 * ( dwdagdure->t4[i][j][k][l].real - dwdagduim->t4[i][j][k][l].imag );
          dwdagdu->t4[i][j][k][l].imag = 
//            0.5 * ( dwdagdure->t4[i][j][k][l].imag - dwdagduim->t4[i][j][k][l].real );
            0.5 * ( dwdagdure->t4[i][j][k][l].imag + dwdagduim->t4[i][j][k][l].real );
        }
      }
    }
  }
}



/* Unitarization with rational approximation */
void u3_unitarize_rational( su3_matrix *V, su3_matrix *W ) {
  int l;
  su3_matrix H, T, X, Z;

  /* hermitian matrix: H=V^+V */
  mult_su3_an( V, V, &H );

  clear_su3mat( &X );

  for( l=1; l<=U3_UNIT_RAT_NTERMS; l++) {
    su3mat_copy( &H, &T );

    /* add d_l to H */
    T.e[0][0].real += d_l_U3_UNIT_RAT[l];
    T.e[1][1].real += d_l_U3_UNIT_RAT[l];
    T.e[2][2].real += d_l_U3_UNIT_RAT[l];

    /* calculate inverse of H+d_l */
    su3_inverse( &T, &Z );

    /* add c_l/(H+d_l) */
    scalar_mult_add_su3_matrix( &X, &Z, c_l_U3_UNIT_RAT[l], &T );

    su3mat_copy( &T, &X );
  }

  /* add c_0 */
  X.e[0][0].real += c_l_U3_UNIT_RAT[0];
  X.e[1][1].real += c_l_U3_UNIT_RAT[0];
  X.e[2][2].real += c_l_U3_UNIT_RAT[0];

  mult_su3_nn( V, &X, W );
}



/* Derivative of a unitarized matrix with rational approximation */
void u3_unit_der_rational( su3_matrix *V, su3_tensor4 *dwdv, su3_tensor4 *dwdagdv ) {
  int i, j, l, p, q, r, s;
  su3_matrix Kl[ U3_UNIT_RAT_NTERMS ]; // store intermediate inverted matrices
  su3_matrix Vdag, H, T, X;
  su3_tensor4 B4;
  complex ftmp, ftmp2, ftmp3, ftmp4;

  /* adjoint, needed later */
  su3_adjoint( V, &Vdag );

  /* hermitian matrix: H=V^+V */
  mult_su3_an( V, V, &H );

  for( l=1; l<=U3_UNIT_RAT_NTERMS; l++) {
    su3mat_copy( &H, &T );

    /* add d_l to H */
    T.e[0][0].real += d_l_U3_UNIT_RAT[l];
    T.e[1][1].real += d_l_U3_UNIT_RAT[l];
    T.e[2][2].real += d_l_U3_UNIT_RAT[l];

    /* calculate inverse of H+d_l and store in Kl array */
    su3_inverse( &T, &( Kl[l-1] ) );
  }

  /* zero out tensor4 */
  for( p=0; p<3; p++) {
    for( r=0; r<3; r++) {
      for( s=0; s<3; s++) {
        for( q=0; q<3; q++) {
          (B4.t4[p][r][s][q]).real = 0.0;
          (B4.t4[p][r][s][q]).imag = 0.0;
        }
      }
    }
  }

  clear_su3mat( &X );

  for( l=1; l<=U3_UNIT_RAT_NTERMS; l++) {
    /* assemble tensor4 out of Kls */
    for( p=0; p<3; p++) {
      for( r=0; r<3; r++) {
        for( s=0; s<3; s++) {
          for( q=0; q<3; q++) {
            CMUL( Kl[l-1].e[p][r], Kl[l-1].e[s][q], ftmp );
            (B4.t4[p][r][s][q]).real += c_l_U3_UNIT_RAT[l] * ftmp.real;
            (B4.t4[p][r][s][q]).imag += c_l_U3_UNIT_RAT[l] * ftmp.imag;
          }
        }
      }
    }
    scalar_mult_add_su3_matrix( &X, &Kl[ l-1 ], c_l_U3_UNIT_RAT[ l ], &T );
    su3mat_copy( &T, &X );
  }


  /* multiplication of tensor4 by V and V^+:
     T_{prsq} = V_{pa} T_{absq} V^+_{br},
     TODO: EXPLAIN addition of the first term */
  for( p=0; p<3; p++) {
    for( r=0; r<3; r++) {
      for( s=0; s<3; s++) {
        for( q=0; q<3; q++) {
          if( p==r ) {
            (dwdv->t4[p][r][s][q]).real = (X.e[s][q]).real;
            (dwdv->t4[p][r][s][q]).imag = (X.e[s][q]).imag;
            if( s==q ) {
              (dwdv->t4[p][r][s][q]).real += c_l_U3_UNIT_RAT[0];
            }
          }
          else {
            (dwdv->t4[p][r][s][q]).real = 0.0;
            (dwdv->t4[p][r][s][q]).imag = 0.0;
          }
          (dwdagdv->t4[p][r][s][q]).real = 0.0;
          (dwdagdv->t4[p][r][s][q]).imag = 0.0;
          for( i=0; i<3; i++ ) {
            ftmp.real = 0.0;
            ftmp.imag = 0.0;
            ftmp3.real = 0.0;
            ftmp3.imag = 0.0;
            for( j=0; j<3; j++ ) {
              // multiplications here: left, right for dW/dV
              //                       right, right for dW^+/dV
              CMUL( B4.t4[i][j][s][q], Vdag.e[j][r], ftmp2 );
              CSUM( ftmp, ftmp2 );
              CMUL( B4.t4[p][i][s][j], Vdag.e[j][q], ftmp4 );
              CSUM( ftmp3, ftmp4 );
            }
            CMUL( V->e[p][i], ftmp, ftmp2 );
            (dwdv->t4[p][r][s][q]).real -= ftmp2.real;
            (dwdv->t4[p][r][s][q]).imag -= ftmp2.imag;
            CMUL( Vdag.e[i][r], ftmp3, ftmp4 );
            (dwdagdv->t4[p][r][s][q]).real -= ftmp4.real;
            (dwdagdv->t4[p][r][s][q]).imag -= ftmp4.imag;
          }
        }
      }
    }
  }
}

/* Analytic unitarization, Hasenfratz, Hoffmann, Schaefer, JHEP05 (2007) 029 */
/* Returns 1 if SVD was used */
void u3_unitarize_analytic( info_t *info, su3_matrix *V, su3_matrix *W ) {
// unless the following switch is defined, u3_unitarize_analytic
// uses double precision, if defined it follows the precision
// set in Makefile
#ifdef U3_UNIT_ANALYTIC_FOLLOW_PREC
  su3_matrix Q, Q2, Q3, S1, S2;
  Real c0, c1, c2, S, S3, R, R2, CQ3, RoS, theta, theta3, pi23, denom;
  Real g0, g1, g2, g0sq, g1sq, g2sq, f0, f1, f2, us, vs, ws;
#else /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  double Ve[3][3][2],Qe[3][3][2],Q2e[3][3][2],Q3e[3][3][2],S2e[3][3][2];
  double c0, c1, c2, S, S3, R, RoS, theta, theta3, pi23, denom;
  double g0, g1, g2, g0sq, g1sq, g2sq, f0, f1, f2, us, vs, ws;
#endif /* U3_UNIT_ANALYTIC_FOLLOW_PREC */

  int i, j;
  size_t nflops = 0;

#ifdef HISQ_REUNIT_ALLOW_SVD
  double Qd[3][3][2];
  //  complex cdet;
  double a1re, a1im, a2re, a2im, a3re, a3im, detre, detim, det_check;
  double sigma[3], Uleft[3][3][2], Vright[3][3][2];
  int perform_svd=0;

  /* get determinant for future comparison */
  a1re=((double)V->e[1][1].real)*((double)V->e[2][2].real)-((double)V->e[1][1].imag)*((double)V->e[2][2].imag)
      -((double)V->e[1][2].real)*((double)V->e[2][1].real)+((double)V->e[1][2].imag)*((double)V->e[2][1].imag);
  /* imag part of (U_22 U_33 - U_23 U_32) */
  a1im=((double)V->e[1][1].real)*((double)V->e[2][2].imag)+((double)V->e[1][1].imag)*((double)V->e[2][2].real)
      -((double)V->e[1][2].real)*((double)V->e[2][1].imag)-((double)V->e[1][2].imag)*((double)V->e[2][1].real);
  /* real part of (U_21 U_33 - U_23 U_31) */
  a2re=((double)V->e[1][0].real)*((double)V->e[2][2].real)-((double)V->e[1][0].imag)*((double)V->e[2][2].imag)
      -((double)V->e[1][2].real)*((double)V->e[2][0].real)+((double)V->e[1][2].imag)*((double)V->e[2][0].imag);
  /* imag part of (U_21 U_33 - U_23 U_31) */
  a2im=((double)V->e[1][0].real)*((double)V->e[2][2].imag)+((double)V->e[1][0].imag)*((double)V->e[2][2].real)
      -((double)V->e[1][2].real)*((double)V->e[2][0].imag)-((double)V->e[1][2].imag)*((double)V->e[2][0].real);
  /* real part of (U_21 U_32 - U_22 U_31) */
  a3re=((double)V->e[1][0].real)*((double)V->e[2][1].real)-((double)V->e[1][0].imag)*((double)V->e[2][1].imag)
      -((double)V->e[1][1].real)*((double)V->e[2][0].real)+((double)V->e[1][1].imag)*((double)V->e[2][0].imag);
  /* imag part of (U_21 U_32 - U_22 U_31) */
  a3im=((double)V->e[1][0].real)*((double)V->e[2][1].imag)+((double)V->e[1][0].imag)*((double)V->e[2][1].real)
      -((double)V->e[1][1].real)*((double)V->e[2][0].imag)-((double)V->e[1][1].imag)*((double)V->e[2][0].real);

  /* real part of det */
  detre=((double)V->e[0][0].real)*a1re-((double)V->e[0][0].imag)*a1im
       -((double)V->e[0][1].real)*a2re+((double)V->e[0][1].imag)*a2im
       +((double)V->e[0][2].real)*a3re-((double)V->e[0][2].imag)*a3im;
  /* imag part of det */
  detim=((double)V->e[0][0].imag)*a1re+((double)V->e[0][0].real)*a1im
       -((double)V->e[0][1].imag)*a2re-((double)V->e[0][1].real)*a2im
       +((double)V->e[0][2].imag)*a3re+((double)V->e[0][2].real)*a3im;
  det_check=detre*detre+detim*detim;

  nflops += 67;

//  cdet = det_su3( V );
//  det_check=cdet.real*cdet.real+cdet.imag*cdet.imag;
#endif /* HISQ_REUNIT_ALLOW_SVD */


#ifdef U3_UNIT_ANALYTIC_FOLLOW_PREC
  /* Hermitian matrix: Q=V^+ V */
  mult_su3_an( V, V, &Q );
  nflops += 198;

  /* Q^2 */
  mult_su3_nn( &Q, &Q, &Q2 );
  nflops += 198;

  /* Q^3 */
  mult_su3_nn( &Q2, &Q, &Q3 );
  nflops += 198;

  /* (real) traces */
  c0 = Q.e[0][0].real + Q.e[1][1].real + Q.e[2][2].real;
  c1 = ( Q2.e[0][0].real + Q2.e[1][1].real + Q2.e[2][2].real ) / 2;
  c2 = ( Q3.e[0][0].real + Q3.e[1][1].real + Q3.e[2][2].real ) / 3;
  nflops += 8;

#else /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  // perform matrix multiplications explicitly in double

  // convert matrix V to double
  for(i=0;i<3;i++)for(j=0;j<3;j++) {
    Ve[i][j][0]=V->e[i][j].real;
    Ve[i][j][1]=V->e[i][j].imag;
  }

  /* Hermitian matrix: Q=V^+ V */
  for(j=0;j<3;j++){
    Qe[0][j][0] = 
        Ve[0][0][0]*Ve[0][j][0] + Ve[0][0][1]*Ve[0][j][1]
      + Ve[1][0][0]*Ve[1][j][0] + Ve[1][0][1]*Ve[1][j][1]
      + Ve[2][0][0]*Ve[2][j][0] + Ve[2][0][1]*Ve[2][j][1];
    Qe[0][j][1] = 
        Ve[0][0][0]*Ve[0][j][1] - Ve[0][0][1]*Ve[0][j][0]
      + Ve[1][0][0]*Ve[1][j][1] - Ve[1][0][1]*Ve[1][j][0]
      + Ve[2][0][0]*Ve[2][j][1] - Ve[2][0][1]*Ve[2][j][0];
    Qe[1][j][0] =
        Ve[0][1][0]*Ve[0][j][0] + Ve[0][1][1]*Ve[0][j][1]
      + Ve[1][1][0]*Ve[1][j][0] + Ve[1][1][1]*Ve[1][j][1]
      + Ve[2][1][0]*Ve[2][j][0] + Ve[2][1][1]*Ve[2][j][1];
    Qe[1][j][1] =
        Ve[0][1][0]*Ve[0][j][1] - Ve[0][1][1]*Ve[0][j][0]
      + Ve[1][1][0]*Ve[1][j][1] - Ve[1][1][1]*Ve[1][j][0]
      + Ve[2][1][0]*Ve[2][j][1] - Ve[2][1][1]*Ve[2][j][0];
    Qe[2][j][0] =
        Ve[0][2][0]*Ve[0][j][0] + Ve[0][2][1]*Ve[0][j][1]
      + Ve[1][2][0]*Ve[1][j][0] + Ve[1][2][1]*Ve[1][j][1]
      + Ve[2][2][0]*Ve[2][j][0] + Ve[2][2][1]*Ve[2][j][1];
    Qe[2][j][1] =
        Ve[0][2][0]*Ve[0][j][1] - Ve[0][2][1]*Ve[0][j][0]
      + Ve[1][2][0]*Ve[1][j][1] - Ve[1][2][1]*Ve[1][j][0]
      + Ve[2][2][0]*Ve[2][j][1] - Ve[2][2][1]*Ve[2][j][0];
  }

  nflops += 3*6*11;

  /* Q^2 */
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    Q2e[i][j][0] =Qe[i][0][0]*Qe[0][j][0]-Qe[i][0][1]*Qe[0][j][1];
    Q2e[i][j][1] =Qe[i][0][0]*Qe[0][j][1]+Qe[i][0][1]*Qe[0][j][0];
    Q2e[i][j][0]+=Qe[i][1][0]*Qe[1][j][0]-Qe[i][1][1]*Qe[1][j][1];
    Q2e[i][j][1]+=Qe[i][1][0]*Qe[1][j][1]+Qe[i][1][1]*Qe[1][j][0];
    Q2e[i][j][0]+=Qe[i][2][0]*Qe[2][j][0]-Qe[i][2][1]*Qe[2][j][1];
    Q2e[i][j][1]+=Qe[i][2][0]*Qe[2][j][1]+Qe[i][2][1]*Qe[2][j][0];
  }

  nflops += 9*22;

  /* Q^3 -- WE NEED ONLY DIAGONAL ELEMENTS */
  for(i=0;i<3;i++){
    Q3e[i][i][0]= Q2e[i][0][0]*Qe[0][i][0]-Q2e[i][0][1]*Qe[0][i][1];
    Q3e[i][i][1]= Q2e[i][0][0]*Qe[0][i][1]+Q2e[i][0][1]*Qe[0][i][0];
    Q3e[i][i][0]+=Q2e[i][1][0]*Qe[1][i][0]-Q2e[i][1][1]*Qe[1][i][1];
    Q3e[i][i][1]+=Q2e[i][1][0]*Qe[1][i][1]+Q2e[i][1][1]*Qe[1][i][0];
    Q3e[i][i][0]+=Q2e[i][2][0]*Qe[2][i][0]-Q2e[i][2][1]*Qe[2][i][1];
    Q3e[i][i][1]+=Q2e[i][2][0]*Qe[2][i][1]+Q2e[i][2][1]*Qe[2][i][0];
  }

  nflops += 9*22;

  /* (real) traces */
  c0 = Qe[0][0][0] + Qe[1][1][0] + Qe[2][2][0];
  c1 = ( Q2e[0][0][0] + Q2e[1][1][0] + Q2e[2][2][0] ) / 2;
  c2 = ( Q3e[0][0][0] + Q3e[1][1][0] + Q3e[2][2][0] ) / 3;

  nflops += 8;
#endif /* U3_UNIT_ANALYTIC_FOLLOW_PREC */

  S = c1/3 - c0 * (c0/18);
  if( fabs(S)<U3_UNIT_ANALYTIC_EPS ) {
    /* eigenvalues of Q */
    g0 = c0/3;
    g1 = c0/3;
    g2 = c0/3;

    nflops += 3;
  }
  else {
    R = c2/2 - c0 * (c1/3) + c0 * c0 * (c0/27);
    S = sqrt(S);
    S3 = S*S*S;
    /* treat possible underflow: R/S^3/2>1.0 leads to acos giving NaN */
    RoS = R/S3;

    nflops += 14;

    if( !( fabs(RoS)<1.0 ) ) {
      if( R>0 ) {
        theta = 0.0;
      }
      else {
        theta = MILC_AB_PI;
      }
    } 
    else {
      theta = acos( RoS );
      if(isnan(theta)){
        printf("Hit NaN in u3_unitarize_analytic()\n");
        printf("RoS=%24.18g\n",RoS);
        printf("Matrix V (row-wise):\n");
        for( i=0; i<3; i++ ) {
          for( j=0; j<3; j++ ) {
            printf( "%24.18g %24.18g\n", V->e[i][j].real, V->e[i][j].imag );
          }
        }
        terminate(0);
      }
      nflops += 1;
    }

    /* eigenvalues of Q */
    theta3 = theta/3;
    pi23 = MILC_AB_TPI / 3;
    g0 = c0/3 + 2 * S * cos( theta3 );
    g1 = c0/3 + 2 * S * cos( theta3 + pi23 );
    g2 = c0/3 + 2 * S * cos( theta3 + 2*pi23 );

    nflops += 20;
  }

#ifdef HISQ_REUNIT_ALLOW_SVD
//  if(det_check!=0) {
//    printf("CHECK: determ. of Hermitian      = %28.18g\n",det_check);
//    printf("       as product of eigenvalues = %28.18g\n",g0*g1*g2);
//    printf("       percentage difference     = %f\n",fabs(det_check-g0*g1*g2)/fabs(det_check)*100);
//  }

#ifndef HISQ_REUNIT_SVD_ONLY

  /* conditions to call SVD */
  if(det_check!=0) {
    if( fabs(det_check-g0*g1*g2)/fabs(det_check)>HISQ_REUNIT_SVD_REL_ERROR ) {
      perform_svd=1;
      nflops += 4;
    }
  }
  
  if(det_check<HISQ_REUNIT_SVD_ABS_ERROR) {
    perform_svd=1;
  }
#else /* HISQ_REUNIT_SVD_ONLY */
  /* exclusively use SVD for finding eigenvalues and reunitarization,
     this is slow since Q, Q^2 and Q^3 are calculated anyway;
     this option is for testing: under normal circumstances SVD is
     rarely used, which makes it harder to test, therefore one can
     SVD with this switch */
  perform_svd=1;
#endif /* ifndef HISQ_REUNIT_SVD_ONLY */

  if( 0!=perform_svd ) {

#ifdef U3_UNIT_ANALYTIC_FOLLOW_PREC
  // convert matrix V to double
  for(i=0;i<3;i++)for(j=0;j<3;j++) {
    Qd[i][j][0]=V->e[i][j].real;
    Qd[i][j][1]=V->e[i][j].imag;
  }
#else /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  // convert matrix V to double
  for(i=0;i<3;i++)for(j=0;j<3;j++) {
    Qd[i][j][0]=Ve[i][j][0];
    Qd[i][j][1]=Ve[i][j][1];
  }

#endif /* U3_UNIT_ANALYTIC_FOLLOW_PREC */

  /* call SVD */
  svd3x3(Qd, sigma, Uleft, Vright, &nflops);
  INFO_HISQ_SVD_COUNTER(info)++;
  

#ifndef HISQ_REUNIT_SVD_ONLY
#ifdef HISQ_SVD_VALUES_INFO
  printf("*** Resort to using SVD ***\n");
  printf("*** printing from node %d ***\n",this_node);
  printf("Eigenvalues from cubic equation:\n");
  printf( "  g0 = %28.18f\n", g0 );
  printf( "  g1 = %28.18f\n", g1 );
  printf( "  g2 = %28.18f\n", g2 );
  printf("Eigenvalues from singular value decomposition:\n");
  printf( "  g0 = %28.18f\n", sigma[0]*sigma[0] );
  printf( "  g1 = %28.18f\n", sigma[1]*sigma[1] );
  printf( "  g2 = %28.18f\n", sigma[2]*sigma[2] );
#endif
#endif
  /* construct Uleft*Vright^+ as a reunitarized link W */
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    W->e[i][j].real = Uleft[i][0][0]*Vright[j][0][0]+Uleft[i][0][1]*Vright[j][0][1];
    W->e[i][j].imag =-Uleft[i][0][0]*Vright[j][0][1]+Uleft[i][0][1]*Vright[j][0][0];
    W->e[i][j].real+= Uleft[i][1][0]*Vright[j][1][0]+Uleft[i][1][1]*Vright[j][1][1];
    W->e[i][j].imag+=-Uleft[i][1][0]*Vright[j][1][1]+Uleft[i][1][1]*Vright[j][1][0];
    W->e[i][j].real+= Uleft[i][2][0]*Vright[j][2][0]+Uleft[i][2][1]*Vright[j][2][1];
    W->e[i][j].imag+=-Uleft[i][2][0]*Vright[j][2][1]+Uleft[i][2][1]*Vright[j][2][0];
  }

  nflops += 9*23;

//      for(i=0;i<3;i++)for(j=0;j<3;j++) {
//        printf( "Qd[%d][%d].re=%26.18e  Qd[%d][%d].im=%26.18e\n",
//                i, j, Qd[i][j][0], i, j, Qd[i][j][1] );
//      }





  }
  else {
#endif /* HISQ_REUNIT_ALLOW_SVD */
//TEMP OUTPUT EIGENVALUES OF Q
//  printf( "  g0 = %28.18f\n", g0 );
//  printf( "  g1 = %28.18f\n", g1 );
//  printf( "  g2 = %28.18f\n", g2 );

//TEMP
//printf("Call to eigen_su3_UdU\n\n");
//  eigen_su3_UdU( V, &g0, &g1, &g2);

  /* roots of eigenvalues */
  g0sq = sqrt( g0 );
  g1sq = sqrt( g1 );
  g2sq = sqrt( g2 );

  /* symmetric combinations */
  us = g1sq + g2sq;
  ws = g1sq * g2sq;
  vs = g0sq * us + ws;
  us += g0sq;
  ws *= g0sq;

  if( ws < U3_UNIT_ANALYTIC_EPS ) {
    printf( "WARNING: u3_unitarize_analytic: ws is too small!\n" );
    printf( "  g0 = %28.18f\n", g0 );
    printf( "  g1 = %28.18f\n", g1 );
    printf( "  g2 = %28.18f\n", g2 );
  }

  denom = ws * ( us*vs - ws );

  /* constants in inverse root expression */
  f0 = ( us*vs*vs - ws*(us*us+vs) ) / denom;
  f1 = ( 2*us*vs - ws - us*us*us ) / denom;
  f2 = us / denom;

  nflops += 30;

#ifdef U3_UNIT_ANALYTIC_FOLLOW_PREC
  /* assemble inverse root: Q^-1/2 = f0 + f1*Q + f2*Q^2 */
  scalar_mult_su3_matrix( &Q2, f2, &S1 );
  scalar_mult_add_su3_matrix( &S1, &Q, f1, &S2 );
  S2.e[0][0].real += f0;
  S2.e[1][1].real += f0;
  S2.e[2][2].real += f0;

  /* W = V*S2 */
  mult_su3_nn( V, &S2, W );

  nflops += 18 + 36 + 3;

#else /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  /* assemble inverse root: Q^-1/2 = f0 + f1*Q + f2*Q^2 */
  S2e[0][0][0] = f0 + f1*Qe[0][0][0] + f2*Q2e[0][0][0];
  S2e[0][0][1] =      f1*Qe[0][0][1] + f2*Q2e[0][0][1];
  S2e[0][1][0] =      f1*Qe[0][1][0] + f2*Q2e[0][1][0];
  S2e[0][1][1] =      f1*Qe[0][1][1] + f2*Q2e[0][1][1];
  S2e[0][2][0] =      f1*Qe[0][2][0] + f2*Q2e[0][2][0];
  S2e[0][2][1] =      f1*Qe[0][2][1] + f2*Q2e[0][2][1];
  S2e[1][0][0] =      f1*Qe[1][0][0] + f2*Q2e[1][0][0];
  S2e[1][0][1] =      f1*Qe[1][0][1] + f2*Q2e[1][0][1];
  S2e[1][1][0] = f0 + f1*Qe[1][1][0] + f2*Q2e[1][1][0];
  S2e[1][1][1] =      f1*Qe[1][1][1] + f2*Q2e[1][1][1];
  S2e[1][2][0] =      f1*Qe[1][2][0] + f2*Q2e[1][2][0];
  S2e[1][2][1] =      f1*Qe[1][2][1] + f2*Q2e[1][2][1];
  S2e[2][0][0] =      f1*Qe[2][0][0] + f2*Q2e[2][0][0];
  S2e[2][0][1] =      f1*Qe[2][0][1] + f2*Q2e[2][0][1];
  S2e[2][1][0] =      f1*Qe[2][1][0] + f2*Q2e[2][1][0];
  S2e[2][1][1] =      f1*Qe[2][1][1] + f2*Q2e[2][1][1];
  S2e[2][2][0] = f0 + f1*Qe[2][2][0] + f2*Q2e[2][2][0];
  S2e[2][2][1] =      f1*Qe[2][2][1] + f2*Q2e[2][2][1];

  /* W = V*S2 */
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    W->e[i][j].real =Ve[i][0][0]*S2e[0][j][0]-Ve[i][0][1]*S2e[0][j][1];
    W->e[i][j].imag =Ve[i][0][0]*S2e[0][j][1]+Ve[i][0][1]*S2e[0][j][0];
    W->e[i][j].real+=Ve[i][1][0]*S2e[1][j][0]-Ve[i][1][1]*S2e[1][j][1];
    W->e[i][j].imag+=Ve[i][1][0]*S2e[1][j][1]+Ve[i][1][1]*S2e[1][j][0];
    W->e[i][j].real+=Ve[i][2][0]*S2e[2][j][0]-Ve[i][2][1]*S2e[2][j][1];
    W->e[i][j].imag+=Ve[i][2][0]*S2e[2][j][1]+Ve[i][2][1]*S2e[2][j][0];
  }

  nflops += 57 + 9*22;

#endif /* U3_UNIT_ANALYTIC_FOLLOW_PREC */

#ifdef HISQ_REUNIT_ALLOW_SVD
  } // end of SVD related if
#endif /* HISQ_REUNIT_ALLOW_SVD */

  info->final_flop += nflops;
}


/* Analytic unitarization, Hasenfratz, Hoffmann, Schaefer, JHEP05 (2007) 029 */
void u3_unitarize_analytic_index( su3_matrix *V, su3_matrix *W, int index_site, int index_dir ) {
  su3_matrix Q, Q2, Q3, S1, S2;
  Real c0, c1, c2, S, S3, R, RoS, theta, theta3, pi23, denom;
  Real g0, g1, g2, g0sq, g1sq, g2sq, f0, f1, f2, us, vs, ws;
  int i, j;


  /* Hermitian matrix: Q=V^+ V */
  mult_su3_an( V, V, &Q );

  /* Q^2 */
  mult_su3_nn( &Q, &Q, &Q2 );

  /* Q^3 */
  mult_su3_nn( &Q2, &Q, &Q3 );

  /* (real) traces */
  c0 = Q.e[0][0].real + Q.e[1][1].real + Q.e[2][2].real;
  c1 = ( Q2.e[0][0].real + Q2.e[1][1].real + Q2.e[2][2].real ) / 2;
  c2 = ( Q3.e[0][0].real + Q3.e[1][1].real + Q3.e[2][2].real ) / 3;

  S = c1/3 - c0 * (c0/18);
  if( fabs(S)<U3_UNIT_ANALYTIC_EPS ) {
    /* eigenvalues of Q */
    g0 = c0/3;
    g1 = c0/3;
    g2 = c0/3;
  }
  else {
    R = c2/2 - c0 * (c1/3) + c0 * c0 * (c0/27);
    S = sqrt(S);
    S3 = S*S*S;
    /* treat possible underflow: R/S^3/2>1.0 leads to acos giving NaN */
    RoS = R/S3;
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
   /* DEBUG: set internal variable */
    lattice[index_site].RoS[index_dir]=RoS;
#endif /* HISQ_REUNITARIZATION_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
    if( !( fabs(RoS)<1.0 ) ) {
      if( R>0 ) {
        theta = 0.0;
      }
      else {
        theta = MILC_AB_PI;
      }
    } 
    else {
      theta = acos( RoS );
      if(isnan(theta)){
        printf("Hit NaN in u3_unitarize_analytic()\n");
        printf("RoS=%24.18g\n",RoS);
        printf("Matrix V (row-wise):\n");
        for( i=0; i<3; i++ ) {
          for( j=0; j<3; j++ ) {
            printf( "%24.18g %24.18g\n", V->e[i][j].real, V->e[i][j].imag );
          }
        }
        terminate(0);
      }
    }

    /* eigenvalues of Q */
    theta3 = theta/3;
    pi23 = MILC_AB_TPI / 3;
    g0 = c0/3 + 2 * S * cos( theta3 );
    g1 = c0/3 + 2 * S * cos( theta3 + pi23 );
    g2 = c0/3 + 2 * S * cos( theta3 + 2*pi23 );
  }

  /* roots of eigenvalues */
  g0sq = sqrt( g0 );
  g1sq = sqrt( g1 );
  g2sq = sqrt( g2 );

  /* symmetric combinations */
  us = g1sq + g2sq;
  ws = g1sq * g2sq;
  vs = g0sq * us + ws;
  us += g0sq;
  ws *= g0sq;

  if( ws < U3_UNIT_ANALYTIC_EPS ) {
    printf( "WARNING: u3_unitarize_analytic: ws is too small!\n" );
  }

  denom = ws * ( us*vs - ws );

  /* constants in inverse root expression */
  f0 = ( us*vs*vs - ws*(us*us+vs) ) / denom;
  f1 = ( 2*us*vs - ws - us*us*us ) / denom;
  f2 = us / denom;

  /* assemble inverse root: Q^-1/2 = f0 + f1*Q + f2*Q^2 */
  scalar_mult_su3_matrix( &Q2, f2, &S1 );
  scalar_mult_add_su3_matrix( &S1, &Q, f1, &S2 );
  S2.e[0][0].real += f0;
  S2.e[1][1].real += f0;
  S2.e[2][2].real += f0;

  /* W = V*S2 */
  mult_su3_nn( V, &S2, W );

#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
   /* calculate and store the phase of U(3) matrices */
  complex detW;
  detW = det_su3( W );
  Real argdet = carg( &detW );
  if( lattice[index_site].on_step_Y[index_dir] < global_current_time_step ) {
    lattice[index_site].on_step_Y[index_dir] = global_current_time_step;
    lattice[index_site].phase_Y_previous[index_dir] =
      lattice[index_site].phase_Y[index_dir]; // previous phase in Y
    lattice[index_site].phase_Y[index_dir] = argdet; // current phase in Y
  }
  /* store eigenvalues of (V^+ V) */
  lattice[index_site].gmin[index_dir]=g0;
  lattice[index_site].gmax[index_dir]=g0;
  if(g1<lattice[index_site].gmin[index_dir])
     lattice[index_site].gmin[index_dir]=g1;
  if(g1>lattice[index_site].gmax[index_dir])
     lattice[index_site].gmax[index_dir]=g1;
  if(g2<lattice[index_site].gmin[index_dir])
     lattice[index_site].gmin[index_dir]=g2;
  if(g2>lattice[index_site].gmax[index_dir])
     lattice[index_site].gmax[index_dir]=g2;
  lattice[index_site].denom[index_dir] = denom;
  /* deviation of W from unitarity */
  mult_su3_an( W, W, &Q );
  Q.e[0][0].real -= 1.;
  Q.e[1][1].real -= 1.;
  Q.e[2][2].real -= 1.;
  lattice[index_site].unitW1[index_dir] = su3_norm_frob( &Q );

  /* EXTRA PIECE: if eigenvalue is less than certain cutoff
                  dump V and W matrices */
#define EIGEN_CUTOFF_DEBUG 1e-5
  if( lattice[index_site].gmin[index_dir] < EIGEN_CUTOFF_DEBUG ) {
    printf("*** Eigenvalue less than cutoff=%f occured ***\n",EIGEN_CUTOFF_DEBUG);
    printf("*** printing from node %d ***\n",this_node);
    printf("Matrix V (row-wise):\n");
    for( i=0; i<3; i++ ) {
      for( j=0; j<3; j++ ) {
        printf( "%24.18g %24.18g\n", V->e[i][j].real, V->e[i][j].imag );
      }
    }
    printf("Matrix W (row-wise):\n");
    for( i=0; i<3; i++ ) {
      for( j=0; j<3; j++ ) {
        printf( "%24.18g %24.18g\n", W->e[i][j].real, W->e[i][j].imag );
      }
    }
    printf("All eigenvalues of (V^+V):\n");
    printf( "g0=%24.18g\n", g0 );
    printf( "g1=%24.18g\n", g1 );
    printf( "g2=%24.18g\n", g2 );
  }
#endif /* HISQ_REUNITARIZATION_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
}


/* Analytic derivative of the unitarized matrix with respect to the original:
   dW/dV and d(W^+)/dV, where W=V(V^+V)^-1/2 */
/* Returns 1 if SVD was used */
void u3_unit_der_analytic( info_t *info, su3_matrix *V, su3_tensor4 *dwdv, 
			   su3_tensor4 *dwdagdv ) {
// unless the following switch is defined, u3_unitarize_analytic
// uses double precision, if defined it follows the precision
// set in Makefile
#ifdef U3_UNIT_ANALYTIC_FOLLOW_PREC
  su3_matrix Q, Q2, Q3, S1, S2, W, Q12;
  su3_matrix VVd, VQ, QVd, QQVd, VQQ, VQVd, PVd, RVd, SVd, Vd;
  Real c0, c1, c2, S, S3, RoS, R, theta, theta3, pi23, denom;
  Real g0, g1, g2, g0sq, g1sq, g2sq, f0, f1, f2, us, vs, ws;
  Real u2, u3, u4, u5, u6, u7, u8, v2, v3, v4, v5, v6, w2, w3, w4, w5;
  Real b00, b01, b02, b11, b12, b22, denom3;
  complex der, ctmp, ctmp2;
#else /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  double Ve[3][3][2], Qe[3][3][2], Q2e[3][3][2], Q3e[3][3][2];
  double We[3][3][2], Q12e[3][3][2];
  double VVde[3][3][2], VQe[3][3][2], QVde[3][3][2], QQVde[3][3][2];
  double VQQe[3][3][2], VQVde[3][3][2], PVde[3][3][2], RVde[3][3][2];
  double SVde[3][3][2], Vde[3][3][2];
  double c0, c1, c2, S, S3, RoS, R, theta, theta3, pi23, denom;
  double g0, g1, g2, g0sq, g1sq, g2sq, f0, f1, f2, us, vs, ws;
  double u2, u3, u4, u5, u6, u7, u8, v2, v3, v4, v5, v6, w2, w3, w4, w5;
  double b00, b01, b02, b11, b12, b22, denom3;
  double der_real, der_imag, ctmp_real, ctmp_imag, ctmp2_real, ctmp2_imag;
#endif /* U3_UNIT_ANALYTIC_FOLLOW_PREC */

  int i, j, m, n;
  size_t nflops = 0;

#ifdef HISQ_REUNIT_ALLOW_SVD
  double Qd[3][3][2];
  //  complex cdet;
  double a1re, a1im, a2re, a2im, a3re, a3im, detre, detim, det_check;
  double sigma[3], Uleft[3][3][2], Vright[3][3][2];
  int perform_svd=0;

  /* get determinant for future comparison */
  a1re=((double)V->e[1][1].real)*((double)V->e[2][2].real)-((double)V->e[1][1].imag)*((double)V->e[2][2].imag)
      -((double)V->e[1][2].real)*((double)V->e[2][1].real)+((double)V->e[1][2].imag)*((double)V->e[2][1].imag);
  /* imag part of (U_22 U_33 - U_23 U_32) */
  a1im=((double)V->e[1][1].real)*((double)V->e[2][2].imag)+((double)V->e[1][1].imag)*((double)V->e[2][2].real)
      -((double)V->e[1][2].real)*((double)V->e[2][1].imag)-((double)V->e[1][2].imag)*((double)V->e[2][1].real);
  /* real part of (U_21 U_33 - U_23 U_31) */
  a2re=((double)V->e[1][0].real)*((double)V->e[2][2].real)-((double)V->e[1][0].imag)*((double)V->e[2][2].imag)
      -((double)V->e[1][2].real)*((double)V->e[2][0].real)+((double)V->e[1][2].imag)*((double)V->e[2][0].imag);
  /* imag part of (U_21 U_33 - U_23 U_31) */
  a2im=((double)V->e[1][0].real)*((double)V->e[2][2].imag)+((double)V->e[1][0].imag)*((double)V->e[2][2].real)
      -((double)V->e[1][2].real)*((double)V->e[2][0].imag)-((double)V->e[1][2].imag)*((double)V->e[2][0].real);
  /* real part of (U_21 U_32 - U_22 U_31) */
  a3re=((double)V->e[1][0].real)*((double)V->e[2][1].real)-((double)V->e[1][0].imag)*((double)V->e[2][1].imag)
      -((double)V->e[1][1].real)*((double)V->e[2][0].real)+((double)V->e[1][1].imag)*((double)V->e[2][0].imag);
  /* imag part of (U_21 U_32 - U_22 U_31) */
  a3im=((double)V->e[1][0].real)*((double)V->e[2][1].imag)+((double)V->e[1][0].imag)*((double)V->e[2][1].real)
      -((double)V->e[1][1].real)*((double)V->e[2][0].imag)-((double)V->e[1][1].imag)*((double)V->e[2][0].real);

  /* real part of det */
  detre=((double)V->e[0][0].real)*a1re-((double)V->e[0][0].imag)*a1im
       -((double)V->e[0][1].real)*a2re+((double)V->e[0][1].imag)*a2im
       +((double)V->e[0][2].real)*a3re-((double)V->e[0][2].imag)*a3im;
  /* imag part of det */
  detim=((double)V->e[0][0].imag)*a1re+((double)V->e[0][0].real)*a1im
       -((double)V->e[0][1].imag)*a2re-((double)V->e[0][1].real)*a2im
       +((double)V->e[0][2].imag)*a3re+((double)V->e[0][2].real)*a3im;
  det_check=detre*detre+detim*detim;

  nflops += 67;

//  cdet = det_su3( V );
//  det_check=cdet.real*cdet.real+cdet.imag*cdet.imag;
#endif /* HISQ_REUNIT_ALLOW_SVD */

#ifdef U3_UNIT_ANALYTIC_FOLLOW_PREC
  /* adjoint */
  su3_adjoint( V, &Vd );

  /* Hermitian matrix: Q=V^+ V */
  mult_su3_an( V, V, &Q ); // +9*22 flops
  nflops += 198;

  /* Q^2 */
  mult_su3_nn( &Q, &Q, &Q2 ); // +9*22 flops
  nflops += 198;

  /* Q^3 */
  mult_su3_nn( &Q2, &Q, &Q3 ); // +9*22 flops
  nflops += 198;

  /* (real) traces */
  c0 = Q.e[0][0].real + Q.e[1][1].real + Q.e[2][2].real;
  c1 = ( Q2.e[0][0].real + Q2.e[1][1].real + Q2.e[2][2].real ) / 2;
  c2 = ( Q3.e[0][0].real + Q3.e[1][1].real + Q3.e[2][2].real ) / 3;
  nflops += 8;
  // +8 flops

#else /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  // perform matrix multiplications explicitly in double

  // convert matrix V to double
  for(i=0;i<3;i++)for(j=0;j<3;j++) {
    Ve[i][j][0]=V->e[i][j].real;
    Ve[i][j][1]=V->e[i][j].imag;
  }

  /* adjoint */
  for(i=0;i<3;i++)for(j=0;j<3;j++) {
    Vde[i][j][0]= V->e[j][i].real;
    Vde[i][j][1]=-V->e[j][i].imag;
  }


  /* Hermitian matrix: Q=V^+ V */
  for(j=0;j<3;j++){
    Qe[0][j][0] = 
        Ve[0][0][0]*Ve[0][j][0] + Ve[0][0][1]*Ve[0][j][1]
      + Ve[1][0][0]*Ve[1][j][0] + Ve[1][0][1]*Ve[1][j][1]
      + Ve[2][0][0]*Ve[2][j][0] + Ve[2][0][1]*Ve[2][j][1];
    Qe[0][j][1] = 
        Ve[0][0][0]*Ve[0][j][1] - Ve[0][0][1]*Ve[0][j][0]
      + Ve[1][0][0]*Ve[1][j][1] - Ve[1][0][1]*Ve[1][j][0]
      + Ve[2][0][0]*Ve[2][j][1] - Ve[2][0][1]*Ve[2][j][0];
    Qe[1][j][0] =
        Ve[0][1][0]*Ve[0][j][0] + Ve[0][1][1]*Ve[0][j][1]
      + Ve[1][1][0]*Ve[1][j][0] + Ve[1][1][1]*Ve[1][j][1]
      + Ve[2][1][0]*Ve[2][j][0] + Ve[2][1][1]*Ve[2][j][1];
    Qe[1][j][1] =
        Ve[0][1][0]*Ve[0][j][1] - Ve[0][1][1]*Ve[0][j][0]
      + Ve[1][1][0]*Ve[1][j][1] - Ve[1][1][1]*Ve[1][j][0]
      + Ve[2][1][0]*Ve[2][j][1] - Ve[2][1][1]*Ve[2][j][0];
    Qe[2][j][0] =
        Ve[0][2][0]*Ve[0][j][0] + Ve[0][2][1]*Ve[0][j][1]
      + Ve[1][2][0]*Ve[1][j][0] + Ve[1][2][1]*Ve[1][j][1]
      + Ve[2][2][0]*Ve[2][j][0] + Ve[2][2][1]*Ve[2][j][1];
    Qe[2][j][1] =
        Ve[0][2][0]*Ve[0][j][1] - Ve[0][2][1]*Ve[0][j][0]
      + Ve[1][2][0]*Ve[1][j][1] - Ve[1][2][1]*Ve[1][j][0]
      + Ve[2][2][0]*Ve[2][j][1] - Ve[2][2][1]*Ve[2][j][0];
  }

  nflops += 3*6*11;

  /* Q^2 */
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    Q2e[i][j][0] =Qe[i][0][0]*Qe[0][j][0]-Qe[i][0][1]*Qe[0][j][1];
    Q2e[i][j][1] =Qe[i][0][0]*Qe[0][j][1]+Qe[i][0][1]*Qe[0][j][0];
    Q2e[i][j][0]+=Qe[i][1][0]*Qe[1][j][0]-Qe[i][1][1]*Qe[1][j][1];
    Q2e[i][j][1]+=Qe[i][1][0]*Qe[1][j][1]+Qe[i][1][1]*Qe[1][j][0];
    Q2e[i][j][0]+=Qe[i][2][0]*Qe[2][j][0]-Qe[i][2][1]*Qe[2][j][1];
    Q2e[i][j][1]+=Qe[i][2][0]*Qe[2][j][1]+Qe[i][2][1]*Qe[2][j][0];
  }

  nflops += 9*22;

  /* Q^3 -- WE NEED ONLY DIAGONAL ELEMENTS */
  for(i=0;i<3;i++){
    Q3e[i][i][0]= Q2e[i][0][0]*Qe[0][i][0]-Q2e[i][0][1]*Qe[0][i][1];
    Q3e[i][i][1]= Q2e[i][0][0]*Qe[0][i][1]+Q2e[i][0][1]*Qe[0][i][0];
    Q3e[i][i][0]+=Q2e[i][1][0]*Qe[1][i][0]-Q2e[i][1][1]*Qe[1][i][1];
    Q3e[i][i][1]+=Q2e[i][1][0]*Qe[1][i][1]+Q2e[i][1][1]*Qe[1][i][0];
    Q3e[i][i][0]+=Q2e[i][2][0]*Qe[2][i][0]-Q2e[i][2][1]*Qe[2][i][1];
    Q3e[i][i][1]+=Q2e[i][2][0]*Qe[2][i][1]+Q2e[i][2][1]*Qe[2][i][0];
  }

  nflops += 9*22;

  /* (real) traces */
  c0 = Qe[0][0][0] + Qe[1][1][0] + Qe[2][2][0];
  c1 = ( Q2e[0][0][0] + Q2e[1][1][0] + Q2e[2][2][0] ) / 2;
  c2 = ( Q3e[0][0][0] + Q3e[1][1][0] + Q3e[2][2][0] ) / 3;

  nflops += 8;
#endif /* U3_UNIT_ANALYTIC_FOLLOW_PREC */

  S = c1/3 - c0 * (c0/18); // +4 flops
  if( fabs(S)<U3_UNIT_ANALYTIC_EPS ) {
    /* eigenvalues of Q */
    g0 = c0/3;
    g1 = c0/3;
    g2 = c0/3;
  }
  else {
    R = c2/2 - c0 * (c1/3) + c0 * c0 * (c0/27); // +8 flops
    S = sqrt(S);
    S3 = S*S*S;
    // +6 flops (sqrt counted as 1 flop)
    /* treat possible underflow: R/S^3/2>1.0 leads to acos giving NaN */
    RoS = R/S3; // +1 flops

    nflops += 14;

    if( !( fabs(RoS)<1.0 ) ) {
      if( R>0 ) {
        theta = 0.0;
      }
      else {
        theta = MILC_AB_PI;
      }
    }
    else {
      theta = acos( RoS );
//      if(isnan(theta)){printf("Hit NaN in u3_unit_der_analytic()\n"); terminate(0);}
      if(isnan(theta)){
        printf("Hit NaN in u3_unit_der_analytic()\n");
        printf("RoS=%24.18g\n",RoS);
        printf("Matrix V (row-wise):\n");
        for( i=0; i<3; i++ ) {
          for( j=0; j<3; j++ ) {
            printf( "%24.18g %24.18g\n", V->e[i][j].real, V->e[i][j].imag );
          }
        }
        terminate(0);
      }
      nflops += 1;
    }
    // +1 flops (acos counted as 1 flop here)

    /* eigenvalues of Q */
    theta3 = theta/3;
    pi23 = MILC_AB_TPI / 3;
    g0 = c0/3 + 2 * S * cos( theta3 );
    g1 = c0/3 + 2 * S * cos( theta3 + pi23 );
    g2 = c0/3 + 2 * S * cos( theta3 + 2*pi23 );

    nflops += 20;
    // +20 flops (cos counted as 1 flop here)
  }
//TEMP OUTPUT EIGENVALUES OF Q
//  printf( "  g0 = %28.18f\n", g0 );
//  printf( "  g1 = %28.18f\n", g1 );
//  printf( "  g2 = %28.18f\n", g2 );

#ifdef HISQ_REUNIT_ALLOW_SVD
//  if(det_check!=0) {
//    printf("CHECK: determ. of Hermitian      = %28.18g\n",det_check);
//    printf("       as product of eigenvalues = %28.18g\n",g0*g1*g2);
//    printf("       percentage difference     = %f\n",fabs(det_check-g0*g1*g2)/fabs(det_check)*100);
//  }

#ifndef HISQ_REUNIT_SVD_ONLY
  /* conditions to call SVD */
  if(det_check!=0) {
    if( fabs(det_check-g0*g1*g2)/fabs(det_check)>HISQ_REUNIT_SVD_REL_ERROR ) {
      perform_svd=1;
      nflops += 4;
    }
  }
  if(det_check<HISQ_REUNIT_SVD_ABS_ERROR) {
    perform_svd=1;
  }
#else /* HISQ_REUNIT_SVD_ONLY */
  /* exclusively use SVD for finding eigenvalues and reunitarization,
     this is slow since Q, Q^2 and Q^3 are calculated anyway;
     this option is for testing: under normal circumstances SVD is
     rarely used, which makes it harder to test, therefore one can
     SVD with this switch */
  perform_svd=1;
#endif /* HISQ_REUNIT_SVD_ONLY */

  if( 0!=perform_svd ) {
    
#ifdef U3_UNIT_ANALYTIC_FOLLOW_PREC
    // convert matrix V to double
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
	Qd[i][j][0]=V->e[i][j].real;
	Qd[i][j][1]=V->e[i][j].imag;
      }
#else /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
    // convert matrix V to double
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
	Qd[i][j][0]=Ve[i][j][0];
	Qd[i][j][1]=Ve[i][j][1];
      }
#endif /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
    
    /* call SVD */
    svd3x3(Qd, sigma, Uleft, Vright, &nflops);
    INFO_HISQ_SVD_COUNTER(info)++;
    
#ifndef HISQ_REUNIT_SVD_ONLY
#ifdef HISQ_SVD_VALUES_INFO
    printf("*** Resort to using svd (force) ***\n");
    printf("*** printing from node %d ***\n",this_node);
    printf("Eigenvalues from cubic equation:\n");
    printf( "  g0 = %28.18f\n", g0 );
    printf( "  g1 = %28.18f\n", g1 );
    printf( "  g2 = %28.18f\n", g2 );
    printf("Eigenvalues from singular value decomposition:\n");
    printf( "  g0 = %28.18f\n", sigma[0]*sigma[0] );
    printf( "  g1 = %28.18f\n", sigma[1]*sigma[1] );
    printf( "  g2 = %28.18f\n", sigma[2]*sigma[2] );
#endif
#endif
    
    g0=sigma[0]*sigma[0];
    g1=sigma[1]*sigma[1];
    g2=sigma[2]*sigma[2];

    nflops += 3;
    
  }
#endif /* HISQ_REUNIT_ALLOW_SVD */

#ifdef HISQ_FORCE_FILTER
#ifdef U3_UNIT_ANALYTIC_FOLLOW_PREC
  Real gmin,g_epsilon;
#else /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  double gmin,g_epsilon;
#endif /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  gmin=g0;
  if(g1<gmin) gmin=g1;
  if(g2<gmin) gmin=g2;
  if(gmin<HISQ_FORCE_FILTER) {
    INFO_HISQ_FORCE_FILTER_COUNTER(info)++;
/*    g_epsilon=HISQ_FORCE_FILTER-gmin;
    if(g_epsilon<0) g_epsilon=-g_epsilon;*/
    g_epsilon=HISQ_FORCE_FILTER;
    g0 += g_epsilon;
    g1 += g_epsilon;
    g2 += g_epsilon;
    // +3 flops
#ifdef U3_UNIT_ANALYTIC_FOLLOW_PREC
  // modify also Q and Q2 matrices
  for(i=0;i<3;i++) {
    Q.e[i][i].real+=g_epsilon;
  }
  mult_su3_nn( &Q, &Q, &Q2 );
  nflops += 198;
#else /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  // modify also Q and Q2 matrices
  for(i=0;i<3;i++) {
    Qe[i][i][0]+=g_epsilon;
  }
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    Q2e[i][j][0] =Qe[i][0][0]*Qe[0][j][0]-Qe[i][0][1]*Qe[0][j][1];
    Q2e[i][j][1] =Qe[i][0][0]*Qe[0][j][1]+Qe[i][0][1]*Qe[0][j][0];
    Q2e[i][j][0]+=Qe[i][1][0]*Qe[1][j][0]-Qe[i][1][1]*Qe[1][j][1];
    Q2e[i][j][1]+=Qe[i][1][0]*Qe[1][j][1]+Qe[i][1][1]*Qe[1][j][0];
    Q2e[i][j][0]+=Qe[i][2][0]*Qe[2][j][0]-Qe[i][2][1]*Qe[2][j][1];
    Q2e[i][j][1]+=Qe[i][2][0]*Qe[2][j][1]+Qe[i][2][1]*Qe[2][j][0];
  }

  nflops += 9*22;

#endif /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
    printf("*** Applying eigenvalue filter in the HISQ reunit force ***\n");
    printf("*** printing from node %d ***\n",this_node);
    printf( "g_epsilon=%24.18g\n", g_epsilon );
    printf( "New eigenvalues:\n" );
    printf( "g0=%24.18g\n", g0 );
    printf( "g1=%24.18g\n", g1 );
    printf( "g2=%24.18g\n", g2 );
#endif /* HISQ_REUNITARIZATION_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
  }
#endif /* HISQ_FORCE_FILTER */

  /* roots of eigenvalues */
  g0sq = sqrt( g0 );
  g1sq = sqrt( g1 );
  g2sq = sqrt( g2 );
  // +3 flops (sqrt counted as 1 flop)

  /* symmetric combinations */
  us = g1sq + g2sq;
  ws = g1sq * g2sq;
  vs = g0sq * us + ws;
  us += g0sq; 
  ws *= g0sq;
  // +6 flops

  nflops += 9;

  if( ws < U3_UNIT_ANALYTIC_EPS ) {
    printf( "WARNING: u3_unit_der_analytic: ws is too small!\n" );
  }

  denom = ws * ( us*vs - ws );
  // +3 flops

  /* constants in inverse root expression */
  f0 = ( us*vs*vs - ws*(us*us+vs) ) / denom;
  f1 = ( 2*us*vs - ws - us*us*us ) / denom;
  f2 = us / denom;
  // +15 flops

  nflops += 18;
//TEMP
//  printf( "  f0 = %28.18f\n", f0 );
//  printf( "  f1 = %28.18f\n", f1 );
//  printf( "  f2 = %28.18f\n", f2 );

#ifdef U3_UNIT_ANALYTIC_FOLLOW_PREC
  /* assemble inverse root: Q^-1/2 = f0 + f1*Q + f2*Q^2 */
  scalar_mult_su3_matrix( &Q2, f2, &S1 ); // +18 flops
  scalar_mult_add_su3_matrix( &S1, &Q, f1, &Q12 ); // +36 flops
  Q12.e[0][0].real += f0;
  Q12.e[1][1].real += f0;
  Q12.e[2][2].real += f0;
  // +3 flops
  nflops += 57;
//TEMP
//  dumpmat( &Q12 );

  /* W = V*Q12 */
  mult_su3_nn( V, &Q12, &W ); // +9*22 flops

  nflops += 198;

#else /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  /* assemble inverse root: Q^-1/2 = f0 + f1*Q + f2*Q^2 */
  Q12e[0][0][0] = f0 + f1*Qe[0][0][0] + f2*Q2e[0][0][0];
  Q12e[0][0][1] =      f1*Qe[0][0][1] + f2*Q2e[0][0][1];
  Q12e[0][1][0] =      f1*Qe[0][1][0] + f2*Q2e[0][1][0];
  Q12e[0][1][1] =      f1*Qe[0][1][1] + f2*Q2e[0][1][1];
  Q12e[0][2][0] =      f1*Qe[0][2][0] + f2*Q2e[0][2][0];
  Q12e[0][2][1] =      f1*Qe[0][2][1] + f2*Q2e[0][2][1];
  Q12e[1][0][0] =      f1*Qe[1][0][0] + f2*Q2e[1][0][0];
  Q12e[1][0][1] =      f1*Qe[1][0][1] + f2*Q2e[1][0][1];
  Q12e[1][1][0] = f0 + f1*Qe[1][1][0] + f2*Q2e[1][1][0];
  Q12e[1][1][1] =      f1*Qe[1][1][1] + f2*Q2e[1][1][1];
  Q12e[1][2][0] =      f1*Qe[1][2][0] + f2*Q2e[1][2][0];
  Q12e[1][2][1] =      f1*Qe[1][2][1] + f2*Q2e[1][2][1];
  Q12e[2][0][0] =      f1*Qe[2][0][0] + f2*Q2e[2][0][0];
  Q12e[2][0][1] =      f1*Qe[2][0][1] + f2*Q2e[2][0][1];
  Q12e[2][1][0] =      f1*Qe[2][1][0] + f2*Q2e[2][1][0];
  Q12e[2][1][1] =      f1*Qe[2][1][1] + f2*Q2e[2][1][1];
  Q12e[2][2][0] = f0 + f1*Qe[2][2][0] + f2*Q2e[2][2][0];
  Q12e[2][2][1] =      f1*Qe[2][2][1] + f2*Q2e[2][2][1];

  nflops += 57;
  /* W = V*Q12 */
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    We[i][j][0] =Ve[i][0][0]*Q12e[0][j][0]-Ve[i][0][1]*Q12e[0][j][1];
    We[i][j][1] =Ve[i][0][0]*Q12e[0][j][1]+Ve[i][0][1]*Q12e[0][j][0];
    We[i][j][0]+=Ve[i][1][0]*Q12e[1][j][0]-Ve[i][1][1]*Q12e[1][j][1];
    We[i][j][1]+=Ve[i][1][0]*Q12e[1][j][1]+Ve[i][1][1]*Q12e[1][j][0];
    We[i][j][0]+=Ve[i][2][0]*Q12e[2][j][0]-Ve[i][2][1]*Q12e[2][j][1];
    We[i][j][1]+=Ve[i][2][0]*Q12e[2][j][1]+Ve[i][2][1]*Q12e[2][j][0];
  }

  nflops += 9*22;

#endif /* U3_UNIT_ANALYTIC_FOLLOW_PREC */

  denom3 = 2*denom*denom*denom; // +3 flops

  /* derivatives of coefficients: B_ij=df_i/dc_j */
  u2 = us*us;
  u3 = u2*us;
  u4 = u3*us;
  u5 = u4*us;
  u6 = u5*us;
  u7 = u6*us;
  u8 = u7*us;
  v2 = vs*vs;
  v3 = v2*vs;
  v4 = v3*vs;
  v5 = v4*vs;
  v6 = v5*vs;
  w2 = ws*ws;
  w3 = w2*ws;
  w4 = w3*ws;
  w5 = w4*ws; // +16 flops
  b00 = -w3*u6 +3*u4*( vs*w3 +v4*ws )
        -u3*( v6 +4*w4 +12*v3*w2 ) + u2*( 16*v2*w3 +3*v5*ws )
        -us*( 8*vs*w4 +3*v4*w2) +w5 +v3*w3;
  b00 /= denom3; // + 33 flops
  b01 = -w2*u7 -v2*ws*u6 +u5*( v4 +6*vs*w2 ) -u4*( 5*w3 +v3*ws )
        -u3*( 2*v5 +6*v2*w2 ) +u2*( 10*vs*w3 +6*v4*ws )
        -us*( 3*w4 +6*v3*w2 ) +2*v2*w3;
  b01 /= denom3; // +38 flops
  b02 = w2*u5 +v2*ws*u4 -u3*( v4 +4*vs*w2 )
        +u2*( 4*w3 +3*v3*ws ) -3*v2*w2*us +vs*w3;
  b02 /= denom3; // +22 flops
  b11 = -ws*u8 -v2*u7 +7*vs*ws*u6 +u5*( 4*v3 -5*w2 ) -16*v2*ws*u4
        -u3*( 4*v4 -16*vs*w2 ) -u2*( 3*w3 -12*v3*ws ) -12*v2*w2*us +3*vs*w3;
  b11 /= denom3; // +37 flops
  b12 = ws*u6 +v2*u5 -5*vs*ws*u4 -u3*( 2*v3 -4*w2 ) +6*v2*ws*u2 -6*vs*w2*us+w3;
  b12 /= denom3; // +22 flops
  b22 = -ws*u4 -v2*u3 +3*vs*ws*u2 -3*w2*us;
  b22 /= denom3; // +12 flops

  nflops += 183;

#ifdef U3_UNIT_ANALYTIC_FOLLOW_PREC
  /* ** create several building blocks for derivative ** */
  mult_su3_nn( V, &Q, &VQ );
  mult_su3_na( &Q, V, &QVd );
  mult_su3_na( V, V, &VVd );
  mult_su3_nn( V, &Q2, &VQQ );
  mult_su3_na( &Q2, V, &QQVd );
  mult_su3_na( &VQ, V, &VQVd ); // +6*9*22 flops
  scalar_mult_su3_matrix( &QVd, b01, &S1 ); // +18 flops
  scalar_mult_add_su3_matrix( &S1, &QQVd, b02, &S2 ); // +36 flops
  scalar_mult_add_su3_matrix( &S2, &Vd, b00, &PVd ); // +36 flops
  scalar_mult_su3_matrix( &QVd, b11, &S1 ); // +18 flops
  scalar_mult_add_su3_matrix( &S1, &QQVd, b12, &S2 ); // +36 flops
  scalar_mult_add_su3_matrix( &S2, &Vd, b01, &RVd ); // +36 flops
  scalar_mult_su3_matrix( &QVd, b12, &S1 ); // +18 flops
  scalar_mult_add_su3_matrix( &S1, &QQVd, b22, &S2 ); // +36 flops
  scalar_mult_add_su3_matrix( &S2, &Vd, b02, &SVd ); // +36 flops

  nflops += 1458;
#else /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  // mult_su3_nn( V, &Q, &VQ );
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    VQe[i][j][0] =Ve[i][0][0]*Qe[0][j][0]-Ve[i][0][1]*Qe[0][j][1];
    VQe[i][j][1] =Ve[i][0][0]*Qe[0][j][1]+Ve[i][0][1]*Qe[0][j][0];
    VQe[i][j][0]+=Ve[i][1][0]*Qe[1][j][0]-Ve[i][1][1]*Qe[1][j][1];
    VQe[i][j][1]+=Ve[i][1][0]*Qe[1][j][1]+Ve[i][1][1]*Qe[1][j][0];
    VQe[i][j][0]+=Ve[i][2][0]*Qe[2][j][0]-Ve[i][2][1]*Qe[2][j][1];
    VQe[i][j][1]+=Ve[i][2][0]*Qe[2][j][1]+Ve[i][2][1]*Qe[2][j][0];
  }
  // mult_su3_na( &Q, V, &QVd );
  // EQUIVALENT TO mult_su3_nn( &Q, &Vd, &QVd), Vd=adjoint of V, defined above
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    QVde[i][j][0] =Qe[i][0][0]*Vde[0][j][0]-Qe[i][0][1]*Vde[0][j][1];
    QVde[i][j][1] =Qe[i][0][0]*Vde[0][j][1]+Qe[i][0][1]*Vde[0][j][0];
    QVde[i][j][0]+=Qe[i][1][0]*Vde[1][j][0]-Qe[i][1][1]*Vde[1][j][1];
    QVde[i][j][1]+=Qe[i][1][0]*Vde[1][j][1]+Qe[i][1][1]*Vde[1][j][0];
    QVde[i][j][0]+=Qe[i][2][0]*Vde[2][j][0]-Qe[i][2][1]*Vde[2][j][1];
    QVde[i][j][1]+=Qe[i][2][0]*Vde[2][j][1]+Qe[i][2][1]*Vde[2][j][0];
  }
  // mult_su3_na( V, V, &VVd );
  // EQUIVALENT TO mult_su3_nn( &V, &Vd, &VVd), Vd=adjoint of V, defined above
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    VVde[i][j][0] =Ve[i][0][0]*Vde[0][j][0]-Ve[i][0][1]*Vde[0][j][1];
    VVde[i][j][1] =Ve[i][0][0]*Vde[0][j][1]+Ve[i][0][1]*Vde[0][j][0];
    VVde[i][j][0]+=Ve[i][1][0]*Vde[1][j][0]-Ve[i][1][1]*Vde[1][j][1];
    VVde[i][j][1]+=Ve[i][1][0]*Vde[1][j][1]+Ve[i][1][1]*Vde[1][j][0];
    VVde[i][j][0]+=Ve[i][2][0]*Vde[2][j][0]-Ve[i][2][1]*Vde[2][j][1];
    VVde[i][j][1]+=Ve[i][2][0]*Vde[2][j][1]+Ve[i][2][1]*Vde[2][j][0];
  }
  // mult_su3_nn( V, &Q2, &VQQ );
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    VQQe[i][j][0] =Ve[i][0][0]*Q2e[0][j][0]-Ve[i][0][1]*Q2e[0][j][1];
    VQQe[i][j][1] =Ve[i][0][0]*Q2e[0][j][1]+Ve[i][0][1]*Q2e[0][j][0];
    VQQe[i][j][0]+=Ve[i][1][0]*Q2e[1][j][0]-Ve[i][1][1]*Q2e[1][j][1];
    VQQe[i][j][1]+=Ve[i][1][0]*Q2e[1][j][1]+Ve[i][1][1]*Q2e[1][j][0];
    VQQe[i][j][0]+=Ve[i][2][0]*Q2e[2][j][0]-Ve[i][2][1]*Q2e[2][j][1];
    VQQe[i][j][1]+=Ve[i][2][0]*Q2e[2][j][1]+Ve[i][2][1]*Q2e[2][j][0];
  }
  // mult_su3_na( &Q2, V, &QQVd );
  // EQUIVALENT TO mult_su3_nn( &Q2, &Vd, &QQVd), Vd=adjoint of V, defined above
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    QQVde[i][j][0] =Q2e[i][0][0]*Vde[0][j][0]-Q2e[i][0][1]*Vde[0][j][1];
    QQVde[i][j][1] =Q2e[i][0][0]*Vde[0][j][1]+Q2e[i][0][1]*Vde[0][j][0];
    QQVde[i][j][0]+=Q2e[i][1][0]*Vde[1][j][0]-Q2e[i][1][1]*Vde[1][j][1];
    QQVde[i][j][1]+=Q2e[i][1][0]*Vde[1][j][1]+Q2e[i][1][1]*Vde[1][j][0];
    QQVde[i][j][0]+=Q2e[i][2][0]*Vde[2][j][0]-Q2e[i][2][1]*Vde[2][j][1];
    QQVde[i][j][1]+=Q2e[i][2][0]*Vde[2][j][1]+Q2e[i][2][1]*Vde[2][j][0];
  }
  // mult_su3_na( &VQ, V, &VQVd );
  // EQUIVALENT TO mult_su3_na( &VQ, &Vd, &VQVd), Vd=adjoint of V, defined above
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    VQVde[i][j][0] =VQe[i][0][0]*Vde[0][j][0]-VQe[i][0][1]*Vde[0][j][1];
    VQVde[i][j][1] =VQe[i][0][0]*Vde[0][j][1]+VQe[i][0][1]*Vde[0][j][0];
    VQVde[i][j][0]+=VQe[i][1][0]*Vde[1][j][0]-VQe[i][1][1]*Vde[1][j][1];
    VQVde[i][j][1]+=VQe[i][1][0]*Vde[1][j][1]+VQe[i][1][1]*Vde[1][j][0];
    VQVde[i][j][0]+=VQe[i][2][0]*Vde[2][j][0]-VQe[i][2][1]*Vde[2][j][1];
    VQVde[i][j][1]+=VQe[i][2][0]*Vde[2][j][1]+VQe[i][2][1]*Vde[2][j][0];
  }

  nflops += 9*198;

  // scalar_mult_su3_matrix( &QVd, b01, &S1 );
  // scalar_mult_add_su3_matrix( &S1, &QQVd, b02, &S2 );
  // scalar_mult_add_su3_matrix( &S2, &Vd, b00, &PVd );

  // scalar_mult_su3_matrix( &QVd, b11, &S1 );
  // scalar_mult_add_su3_matrix( &S1, &QQVd, b12, &S2 );
  // scalar_mult_add_su3_matrix( &S2, &Vd, b01, &RVd );

  // scalar_mult_su3_matrix( &QVd, b12, &S1 );
  // scalar_mult_add_su3_matrix( &S1, &QQVd, b22, &S2 );
  // scalar_mult_add_su3_matrix( &S2, &Vd, b02, &SVd );
  for(i=0;i<3;i++)for(j=0;j<3;j++){
    PVde[i][j][0] = b00*Vde[i][j][0] + b01*QVde[i][j][0] + b02*QQVde[i][j][0];
    PVde[i][j][1] = b00*Vde[i][j][1] + b01*QVde[i][j][1] + b02*QQVde[i][j][1];
    RVde[i][j][0] = b01*Vde[i][j][0] + b11*QVde[i][j][0] + b12*QQVde[i][j][0];
    RVde[i][j][1] = b01*Vde[i][j][1] + b11*QVde[i][j][1] + b12*QQVde[i][j][1];
    SVde[i][j][0] = b02*Vde[i][j][0] + b12*QVde[i][j][0] + b22*QQVde[i][j][0];
    SVde[i][j][1] = b02*Vde[i][j][1] + b12*QVde[i][j][1] + b22*QQVde[i][j][1];
  }

  nflops += 9*30;

#endif /* U3_UNIT_ANALYTIC_FOLLOW_PREC */

#ifdef U3_UNIT_ANALYTIC_FOLLOW_PREC
  /* assemble the derivative rank 4 tensor */
  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      for( m=0; m<3; m++) {
        for( n=0; n<3; n++) {
          der.real = 0.0;
          der.imag = 0.0;
          /* dW/dV */
          if( i==m ) {
            der.real += Q12.e[n][j].real;
            der.imag += Q12.e[n][j].imag;
          }
          if( j==n ) {
            der.real += f1 * VVd.e[i][m].real + f2 * VQVd.e[i][m].real;
            der.imag += f1 * VVd.e[i][m].imag + f2 * VQVd.e[i][m].imag;
          } // +8 flops
          CMUL( V->e[i][j], PVd.e[n][m], ctmp );
          CSUM( der, ctmp );
          CMUL( VQ.e[i][j], RVd.e[n][m], ctmp );
          CSUM( der, ctmp );
          CMUL( VQQ.e[i][j], SVd.e[n][m], ctmp );
          CSUM( der, ctmp );
          CMUL( VVd.e[i][m], Q.e[n][j], ctmp );
          der.real += f2 * ctmp.real;
          der.imag += f2 * ctmp.imag;
          dwdv->t4[i][m][n][j].real = der.real;
          dwdv->t4[i][m][n][j].imag = der.imag;
          /* dW^+/dV */
          CMUL( Vd.e[i][j], PVd.e[n][m], der );
          CMUL( QVd.e[i][j], RVd.e[n][m], ctmp );
          CSUM( der, ctmp );
          CMUL( Vd.e[i][m], Vd.e[n][j], ctmp );
          der.real += f1 * ctmp.real;
          der.imag += f1 * ctmp.imag;
          CMUL( QQVd.e[i][j], SVd.e[n][m], ctmp );
          CSUM( der, ctmp );
          CMUL( Vd.e[i][m], QVd.e[n][j], ctmp );
          CMUL( Vd.e[n][j], QVd.e[i][m], ctmp2 );
          der.real += f2 * ( ctmp.real + ctmp2.real );
          der.imag += f2 * ( ctmp.imag + ctmp2.imag );
          dwdagdv->t4[i][m][n][j].real = der.real;
          dwdagdv->t4[i][m][n][j].imag = der.imag;

	  nflops += 84;
        }
      }
    } // +81*(32+50)  flops
  }
#else /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  /* assemble the derivative rank 4 tensor */
  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      for( m=0; m<3; m++) {
        for( n=0; n<3; n++) {
          der_real = 0.0;
          der_imag = 0.0;
          /* dW/dV */
          if( i==m ) {
            der_real += Q12e[n][j][0];
            der_imag += Q12e[n][j][1];
          }
          if( j==n ) {
            der_real += f1 * VVde[i][m][0] + f2 * VQVde[i][m][0];
            der_imag += f1 * VVde[i][m][1] + f2 * VQVde[i][m][1];
          } // +8 flops
          der_real += Ve[i][j][0]*PVde[n][m][0] - Ve[i][j][1]*PVde[n][m][1];
          der_imag += Ve[i][j][0]*PVde[n][m][1] + Ve[i][j][1]*PVde[n][m][0];
          der_real += VQe[i][j][0]*RVde[n][m][0] - VQe[i][j][1]*RVde[n][m][1];
          der_imag += VQe[i][j][0]*RVde[n][m][1] + VQe[i][j][1]*RVde[n][m][0];
          der_real += VQQe[i][j][0]*SVde[n][m][0] - VQQe[i][j][1]*SVde[n][m][1];
          der_imag += VQQe[i][j][0]*SVde[n][m][1] + VQQe[i][j][1]*SVde[n][m][0];
          der_real +=f2*(VVde[i][m][0]*Qe[n][j][0] - VVde[i][m][1]*Qe[n][j][1]);
          der_imag +=f2*(VVde[i][m][0]*Qe[n][j][1] + VVde[i][m][1]*Qe[n][j][0]);
          dwdv->t4[i][m][n][j].real = (Real)der_real;
          dwdv->t4[i][m][n][j].imag = (Real)der_imag;
	  // +32 flops
          /* dW^+/dV */
          der_real  = Vde[i][j][0]*PVde[n][m][0] - Vde[i][j][1]*PVde[n][m][1];
          der_imag  = Vde[i][j][0]*PVde[n][m][1] + Vde[i][j][1]*PVde[n][m][0];
          der_real += QVde[i][j][0]*RVde[n][m][0] - QVde[i][j][1]*RVde[n][m][1];
          der_imag += QVde[i][j][0]*RVde[n][m][1] + QVde[i][j][1]*RVde[n][m][0];
          der_real +=f1*(Vde[i][m][0]*Vde[n][j][0] - Vde[i][m][1]*Vde[n][j][1]);
          der_imag +=f1*(Vde[i][m][0]*Vde[n][j][1] + Vde[i][m][1]*Vde[n][j][0]);
          der_real += QQVde[i][j][0]*SVde[n][m][0]-QQVde[i][j][1]*SVde[n][m][1];
          der_imag += QQVde[i][j][0]*SVde[n][m][1]+QQVde[i][j][1]*SVde[n][m][0];
          ctmp_real = Vde[i][m][0]*QVde[n][j][0] - Vde[i][m][1]*QVde[n][j][1];
          ctmp_imag = Vde[i][m][0]*QVde[n][j][1] + Vde[i][m][1]*QVde[n][j][0];
          ctmp2_real= Vde[n][j][0]*QVde[i][m][0] - Vde[n][j][1]*QVde[i][m][1];
          ctmp2_imag= Vde[n][j][0]*QVde[i][m][1] + Vde[n][j][1]*QVde[i][m][0];
          der_real += f2 * ( ctmp_real + ctmp2_real );
          der_imag += f2 * ( ctmp_imag + ctmp2_imag );
          dwdagdv->t4[i][m][n][j].real = (Real)der_real;
          dwdagdv->t4[i][m][n][j].imag = (Real)der_imag;
	  // +48 flops
	  nflops += 88;
        }
      }
    }
  }
#endif /* U3_UNIT_ANALYTIC_FOLLOW_PREC */
  info->final_flop += nflops;
  
}



/* copy rank 4 tensor: a -> b */
void su3t4_copy( su3_tensor4 *a, su3_tensor4 *b ) {
  int i, j, k, l;

  for( i=0; i<3; i++ ) {
    for( j=0; j<3; j++ ) {
      for( k=0; k<3; k++ ) {
        for( l=0; l<3; l++ ) {
          b->t4[i][j][k][l].real = a->t4[i][j][k][l].real;
          b->t4[i][j][k][l].imag = a->t4[i][j][k][l].imag;
        }
      }
    }
  }
}




/* **************************************************
   SVD stuff needs to be put into a separate file
   ************************************************** */
/* Singular value decomposition (SVD) for
   3x3 complex matrix
   A.Bazavov, Feb 20 2009


   Algorithm sketch:

   SVD is performed in two steps:
     1) 3x3 complex matrix is reduced to real bidiagonal form
        with Householder transformations,
        3x3 matrix requires three left and two right such
        transformations
     2) bidiagonal matrix has the form
        [ b00 b01   0 ]
        [   0 b11 b12 ]
        [   0   0 b22 ]
        It is iteratively diagonalized with QR algorithm with shifts
        (it constructs Given rotations).
        There are many special cases (such as b00==0, b11==0, etc.)
        that are handled separately. If b01==0 then an auxiliary
        routine svd2x2bidiag is used to decompose the lower 2x2 block,
        if b12==0 the same is done for upper block.
        svd2x2bidiag is a separate routine because there are special
        cases for 2x2 superdiagonal matrix that need to be handled.

   This routine needs to be stable for singular matrices. Therefore, most
   of the operations are done to avoid underflow/overflow, for example,
   if norm=sqrt(a^2+b^2) then the calculation proceeds as:
     min=min(a,b), max=max(a,b)
     norm=max*sqrt(1+(min/max)^2)
   and so on.
*/

#include <stdio.h>
#include <math.h>


/* debugging define: prints a lot(!) */
/*#define SVD3x3_DEBUG*/

/* define precision for chopping small values */
#define SVD3x3_PREC 5e-16

/* defines that allow to remap input arrays easily,
   in the routine internal arrays in double precision are used */
#define A00re A[0][0][0]
#define A00im A[0][0][1]
#define A01re A[0][1][0]
#define A01im A[0][1][1]
#define A02re A[0][2][0]
#define A02im A[0][2][1]
#define A10re A[1][0][0]
#define A10im A[1][0][1]
#define A11re A[1][1][0]
#define A11im A[1][1][1]
#define A12re A[1][2][0]
#define A12im A[1][2][1]
#define A20re A[2][0][0]
#define A20im A[2][0][1]
#define A21re A[2][1][0]
#define A21im A[2][1][1]
#define A22re A[2][2][0]
#define A22im A[2][2][1]
#define U00re U[0][0][0]
#define U00im U[0][0][1]
#define U01re U[0][1][0]
#define U01im U[0][1][1]
#define U02re U[0][2][0]
#define U02im U[0][2][1]
#define U10re U[1][0][0]
#define U10im U[1][0][1]
#define U11re U[1][1][0]
#define U11im U[1][1][1]
#define U12re U[1][2][0]
#define U12im U[1][2][1]
#define U20re U[2][0][0]
#define U20im U[2][0][1]
#define U21re U[2][1][0]
#define U21im U[2][1][1]
#define U22re U[2][2][0]
#define U22im U[2][2][1]
#define V00re V[0][0][0]
#define V00im V[0][0][1]
#define V01re V[0][1][0]
#define V01im V[0][1][1]
#define V02re V[0][2][0]
#define V02im V[0][2][1]
#define V10re V[1][0][0]
#define V10im V[1][0][1]
#define V11re V[1][1][0]
#define V11im V[1][1][1]
#define V12re V[1][2][0]
#define V12im V[1][2][1]
#define V20re V[2][0][0]
#define V20im V[2][0][1]
#define V21re V[2][1][0]
#define V21im V[2][1][1]
#define V22re V[2][2][0]
#define V22im V[2][2][1]
#define b00 P[0][0][0]
#define b01 P[0][1][0]
#define b02 P[0][2][0]
#define b10 P[1][0][0]
#define b11 P[1][1][0]
#define b12 P[1][2][0]
#define b20 P[2][0][0]
#define b21 P[2][1][0]
#define b22 P[2][2][0]


/* forward declaration */
int svd2x2bidiag(double *a00, double *a01, double *a11, 
		 double U2[2][2], double V2[2][2], size_t *nflops);


/* Input: A -- 3x3 complex matrix,
   Output: sigma[3] -- singular values,
           U,V -- U(3) matrices such, that
           A=U Sigma V^+ */
int svd3x3(double A[3][3][2], double *sigma, double U[3][3][2], 
	   double V[3][3][2], size_t *nf) {
  double Ad[3][3][2], P[3][3][2], Q[3][3][2];
  double U1[3][3][2], U2[3][3][2], U3[3][3][2], V1[3][3][2], V2[3][3][2];
  double UO3[3][3], VO3[3][3], v[3][2];
  double UO2[2][2], VO2[2][2];
  register double a, b, c, d, factor, norm, min, max, taure, tauim, beta;
  register double m11, m12, m22, dm, lambdamax, cosphi, sinphi, tanphi, cotphi;
  register int i, iter;
  size_t nflops = 0;

  /* format of external matrices A, U and V can be arbitrary,
     therefore this routine uses defines (above) to access them
     and never reads A, U and V directly */

  /* original matrix can be in single precision,
     so copy it into double */
  Ad[0][0][0]=(double)A00re; Ad[0][0][1]=(double)A00im;
  Ad[0][1][0]=(double)A01re; Ad[0][1][1]=(double)A01im;
  Ad[0][2][0]=(double)A02re; Ad[0][2][1]=(double)A02im;
  Ad[1][0][0]=(double)A10re; Ad[1][0][1]=(double)A10im;
  Ad[1][1][0]=(double)A11re; Ad[1][1][1]=(double)A11im;
  Ad[1][2][0]=(double)A12re; Ad[1][2][1]=(double)A12im;
  Ad[2][0][0]=(double)A20re; Ad[2][0][1]=(double)A20im;
  Ad[2][1][0]=(double)A21re; Ad[2][1][1]=(double)A21im;
  Ad[2][2][0]=(double)A22re; Ad[2][2][1]=(double)A22im;


  i=0;

  /* *** Step 1: build first left reflector v,
                 calculate first left rotation U1,
                 apply to the original matrix A *** */
  /* calculate norm of ( A[10] )
                       ( A[20] ) vector
     with minimal loss of accuracy (similar to BLAS) */
  c = 1.;
  factor = fabs( Ad[1][0][0] );
  a = fabs( Ad[1][0][1] );
  if( a!=0 ) {
    if( factor < a ) {
      c = 1 + (factor/a)*(factor/a);
      factor = a;
    }
    else {
      c = 1 + (a/factor)*(a/factor);
    }
  }
  a = fabs( Ad[2][0][0] );
  if( a!=0 ) {
    if( factor < a ) {
      c = 1 + c*(factor/a)*(factor/a);
      factor = a;
    }
    else {
      c += (a/factor)*(a/factor);
    }
  }
  a = fabs( Ad[2][0][1] );
  if( a!=0 ) {
    if( factor < a ) {
      c = 1 + c*(factor/a)*(factor/a);
      factor = a;
    }
    else {
      c += (a/factor)*(a/factor);
    }
  }
  norm = factor*sqrt(c);

  nflops += 15;

  if( norm==0 && Ad[0][0][1]==0 ) { /* no rotation needed */
#ifdef SVD3x3_DEBUG
    printf("Step 1: no rotation needed\n");
#endif /* SVD3x3_DEBUG */
    U1[0][0][0]=1.; U1[0][0][1]=0.;
    U1[0][1][0]=0.; U1[0][1][1]=0.;
    U1[0][2][0]=0.; U1[0][2][1]=0.;
    U1[1][0][0]=0.; U1[1][0][1]=0.;
    U1[1][1][0]=1.; U1[1][1][1]=0.;
    U1[1][2][0]=0.; U1[1][2][1]=0.;
    U1[2][0][0]=0.; U1[2][0][1]=0.;
    U1[2][1][0]=0.; U1[2][1][1]=0.;
    U1[2][2][0]=1.; U1[2][2][1]=0.;
    P[0][0][0]=Ad[0][0][0]; P[0][0][1]=Ad[0][0][1];
    P[1][0][0]=Ad[1][0][0]; P[1][0][1]=Ad[1][0][1];
    P[2][0][0]=Ad[2][0][0]; P[2][0][1]=Ad[2][0][1];
    P[0][1][0]=Ad[0][1][0]; P[0][1][1]=Ad[0][1][1];
    P[1][1][0]=Ad[1][1][0]; P[1][1][1]=Ad[1][1][1];
    P[2][1][0]=Ad[2][1][0]; P[2][1][1]=Ad[2][1][1];
    P[0][2][0]=Ad[0][2][0]; P[0][2][1]=Ad[0][2][1];
    P[1][2][0]=Ad[1][2][0]; P[1][2][1]=Ad[1][2][1];
    P[2][2][0]=Ad[2][2][0]; P[2][2][1]=Ad[2][2][1];
  }
  else {


    /* get the norm of full first column of A matrix */
    c=1.;
    factor = norm;
    a = fabs( Ad[0][0][0] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + (factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    a = fabs( Ad[0][0][1] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + c*(factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    beta = factor*sqrt(c); /* norm of first column */
    if( Ad[0][0][0]>0 ) {
      beta = -beta;
    }

#ifdef SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* SVD3x3_DEBUG */


    /* a=Re(A_00-beta), b=Im(A_00-beta) */
    a=Ad[0][0][0]-beta; b=Ad[0][0][1];
    /* norm=sqrt(a^2+b^2) */
    c=fabs(a); d=fabs(b);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      norm = max;
    }
    else {
      c = min/max;
      norm = max*sqrt(1+c*c);
    }
    /* c=a/norm, d=b/norm */
    c=a/norm; d=b/norm;


    /* construct reflector (vector "v" for Householder transformation)
       v_0=1 */
    v[0][0]=1.; v[0][1]=0.;
    /* v_1=A_10/(A_00-beta)=A_10/(a+ib)=(A_10*(a-ib))/norm^2=(A_10/norm)*((a-ib)/norm)
          =(A_10/norm)*(c-id)=|a=Re(A_10)/norm,b=Im(A_10)/norm|=(a+ib)*(c-id)
          =(a*c+b*d)+i(b*c-a*d) */
    a=Ad[1][0][0]/norm; b=Ad[1][0][1]/norm;
    v[1][0]=a*c+b*d;
    v[1][1]=b*c-a*d;
    /* v_2=A_20/(A_00-beta)=A_20/(a+ib)=(A_20*(a-ib))/norm^2=(A_20/norm)*((a-ib)/norm)
          =(A_20/norm)*(c-id)=|a=Re(A_20)/norm,b=Im(A_20)/norm|=(a+ib)*(c-id)
          =(a*c+b*d)+i(b*c-a*d) */
    a=Ad[2][0][0]/norm; b=Ad[2][0][1]/norm;
    v[2][0]=a*c+b*d;
    v[2][1]=b*c-a*d;
#ifdef SVD3x3_DEBUG
for(i=0;i<3;i++) {
  printf("v[%d].re=%28.18e  v[%d].im=%28.18e\n",i,v[i][0],i,v[i][1]);
}
#endif /* SVD3x3_DEBUG */

    /* calcualate tau (coefficient for reflector) */
    taure=(beta-Ad[0][0][0])/beta;
    tauim=(Ad[0][0][1])/beta;


    /* assemble left unitary matrix U1=I-tau^+*v*v^+ (store in U1[3][3][2])
       U1_00=A_00/beta */
    U1[0][0][0]=(Ad[0][0][0])/beta;
    U1[0][0][1]=(Ad[0][0][1])/beta;
    /* U1_10=A_10/beta */
    U1[1][0][0]=(Ad[1][0][0])/beta;
    U1[1][0][1]=(Ad[1][0][1])/beta;
    /* U1_20=A_20/beta */
    U1[2][0][0]=(Ad[2][0][0])/beta;
    U1[2][0][1]=(Ad[2][0][1])/beta;
    /* U1_01=-tau^+*v_1^+=-(tau*v_1)^+ */
    U1[0][1][0]=-(taure*v[1][0]-tauim*v[1][1]);
    U1[0][1][1]=taure*v[1][1]+tauim*v[1][0];
    /* U1_11=1-tau^+*v_1*v_1^+ */
    a=v[1][0]*v[1][0]+v[1][1]*v[1][1];
    U1[1][1][0]=1-taure*a;
    U1[1][1][1]=tauim*a;
    /* U1_21=-tau^+*v_2*v_1^+ */
      /* v_2*v_1^+ */
      a=v[2][0]*v[1][0]+v[2][1]*v[1][1];
      b=-v[2][0]*v[1][1]+v[2][1]*v[1][0];
    U1[2][1][0]=-(taure*a+tauim*b);
    U1[2][1][1]=-(taure*b-tauim*a);
    /* U1_02=-tau^+*v_2^+=-(tau*v_2)^+ */
    U1[0][2][0]=-(taure*v[2][0]-tauim*v[2][1]);
    U1[0][2][1]=taure*v[2][1]+tauim*v[2][0];
    /* U1_12=-tau^+*v_1*v_2^+ */
      /* v_1*v_2^+ */
      a=v[1][0]*v[2][0]+v[1][1]*v[2][1];
      b=-v[1][0]*v[2][1]+v[1][1]*v[2][0];
    U1[1][2][0]=-(taure*a+tauim*b);
    U1[1][2][1]=-(taure*b-tauim*a);
    /* U1_22=1-tau^+*v_2*v_2^+ */
    a=v[2][0]*v[2][0]+v[2][1]*v[2][1];
    U1[2][2][0]=1-taure*a;
    U1[2][2][1]=tauim*a;

    nflops += 91;


    /* apply the transformation to A matrix and store the result in P
       P=U^+A */
    P[0][0][0]=beta;
    P[0][0][1]=0;
    P[1][0][0]=0;
    P[1][0][1]=0;
    P[2][0][0]=0;
    P[2][0][1]=0;
    /* P_01=U1_00^+*A_01+U1_10^+*A_11+U1_20^+*A_21 */
    P[0][1][0]=U1[0][0][0]*Ad[0][1][0]+U1[0][0][1]*Ad[0][1][1]
              +U1[1][0][0]*Ad[1][1][0]+U1[1][0][1]*Ad[1][1][1]
              +U1[2][0][0]*Ad[2][1][0]+U1[2][0][1]*Ad[2][1][1];
    P[0][1][1]=U1[0][0][0]*Ad[0][1][1]-U1[0][0][1]*Ad[0][1][0]
              +U1[1][0][0]*Ad[1][1][1]-U1[1][0][1]*Ad[1][1][0]
              +U1[2][0][0]*Ad[2][1][1]-U1[2][0][1]*Ad[2][1][0];
    /* P_02=U1_00^+*A_02+U1_10^+*A_12+U1_20^+*A_22 */
    P[0][2][0]=U1[0][0][0]*Ad[0][2][0]+U1[0][0][1]*Ad[0][2][1]
              +U1[1][0][0]*Ad[1][2][0]+U1[1][0][1]*Ad[1][2][1]
              +U1[2][0][0]*Ad[2][2][0]+U1[2][0][1]*Ad[2][2][1];
    P[0][2][1]=U1[0][0][0]*Ad[0][2][1]-U1[0][0][1]*Ad[0][2][0]
              +U1[1][0][0]*Ad[1][2][1]-U1[1][0][1]*Ad[1][2][0]
              +U1[2][0][0]*Ad[2][2][1]-U1[2][0][1]*Ad[2][2][0];
    /* P_11=U1_01^+*A_01+U1_11^+*A_11+U1_21^+*A_21 */
    P[1][1][0]=U1[0][1][0]*Ad[0][1][0]+U1[0][1][1]*Ad[0][1][1]
              +U1[1][1][0]*Ad[1][1][0]+U1[1][1][1]*Ad[1][1][1]
              +U1[2][1][0]*Ad[2][1][0]+U1[2][1][1]*Ad[2][1][1];
    P[1][1][1]=U1[0][1][0]*Ad[0][1][1]-U1[0][1][1]*Ad[0][1][0]
              +U1[1][1][0]*Ad[1][1][1]-U1[1][1][1]*Ad[1][1][0]
              +U1[2][1][0]*Ad[2][1][1]-U1[2][1][1]*Ad[2][1][0];
    /* P_12=U1_01^+*A_02+U1_11^+*A_12+U1_21^+*A_22 */
    P[1][2][0]=U1[0][1][0]*Ad[0][2][0]+U1[0][1][1]*Ad[0][2][1]
              +U1[1][1][0]*Ad[1][2][0]+U1[1][1][1]*Ad[1][2][1]
              +U1[2][1][0]*Ad[2][2][0]+U1[2][1][1]*Ad[2][2][1];
    P[1][2][1]=U1[0][1][0]*Ad[0][2][1]-U1[0][1][1]*Ad[0][2][0]
              +U1[1][1][0]*Ad[1][2][1]-U1[1][1][1]*Ad[1][2][0]
              +U1[2][1][0]*Ad[2][2][1]-U1[2][1][1]*Ad[2][2][0];
    /* P_21=U1_02^+*A_01+U1_12^+*A_11+U1_22^+*A_21 */
    P[2][1][0]=U1[0][2][0]*Ad[0][1][0]+U1[0][2][1]*Ad[0][1][1]
              +U1[1][2][0]*Ad[1][1][0]+U1[1][2][1]*Ad[1][1][1]
              +U1[2][2][0]*Ad[2][1][0]+U1[2][2][1]*Ad[2][1][1];
    P[2][1][1]=U1[0][2][0]*Ad[0][1][1]-U1[0][2][1]*Ad[0][1][0]
              +U1[1][2][0]*Ad[1][1][1]-U1[1][2][1]*Ad[1][1][0]
              +U1[2][2][0]*Ad[2][1][1]-U1[2][2][1]*Ad[2][1][0];
    /* P_22=U1_02^+*A_02+U1_12^+*A_12+U1_22^+*A_22 */
    P[2][2][0]=U1[0][2][0]*Ad[0][2][0]+U1[0][2][1]*Ad[0][2][1]
              +U1[1][2][0]*Ad[1][2][0]+U1[1][2][1]*Ad[1][2][1]
              +U1[2][2][0]*Ad[2][2][0]+U1[2][2][1]*Ad[2][2][1];
    P[2][2][1]=U1[0][2][0]*Ad[0][2][1]-U1[0][2][1]*Ad[0][2][0]
              +U1[1][2][0]*Ad[1][2][1]-U1[1][2][1]*Ad[1][2][0]
              +U1[2][2][0]*Ad[2][2][1]-U1[2][2][1]*Ad[2][2][0];

    nflops += 9*12;

  }
#ifdef SVD3x3_DEBUG
  {
    int j;
    printf("Left unitary matrix U1:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
	printf( "U1[%d][%d].re=%26.18e  U1[%d][%d].im=%26.18e\n",
		i, j, U1[i][j][0], i, j, U1[i][j][1] );
      }
  }
#endif /* SVD3x3_DEBUG */



  /* *** Step 2: build first right reflector v,
                 calculate first right rotation V1,
                 apply to the matrix P from step 1 *** */
  /* calculate norm of ( P[02] )
     with minimal loss of accuracy */
  a=fabs( P[0][2][0] ); b=fabs( P[0][2][1] );
  /* norm=sqrt(a^2+b^2) */
  if( a>b ) {
    max=a; min=b;
  }
  else {
    max=b; min=a;
  }
  if( min==0 ) {
    norm = max;
  }
  else {
    c = min/max;
    norm = max*sqrt(1+c*c);
  }

  if( norm==0 && P[0][1][1]==0 ) { /* no rotation needed */
#ifdef SVD3x3_DEBUG
    printf("Step 2: no rotation needed\n");
#endif /* SVD3x3_DEBUG */
    V1[0][0][0]=1.; V1[0][0][1]=0.;
    V1[0][1][0]=0.; V1[0][1][1]=0.;
    V1[0][2][0]=0.; V1[0][2][1]=0.;
    V1[1][0][0]=0.; V1[1][0][1]=0.;
    V1[1][1][0]=1.; V1[1][1][1]=0.;
    V1[1][2][0]=0.; V1[1][2][1]=0.;
    V1[2][0][0]=0.; V1[2][0][1]=0.;
    V1[2][1][0]=0.; V1[2][1][1]=0.;
    V1[2][2][0]=1.; V1[2][2][1]=0.;
    Q[0][0][0]=P[0][0][0]; Q[0][0][1]=P[0][0][1];
    Q[1][0][0]=P[1][0][0]; Q[1][0][1]=P[1][0][1];
    Q[2][0][0]=P[2][0][0]; Q[2][0][1]=P[2][0][1];
    Q[0][1][0]=P[0][1][0]; Q[0][1][1]=P[0][1][1];
    Q[1][1][0]=P[1][1][0]; Q[1][1][1]=P[1][1][1];
    Q[2][1][0]=P[2][1][0]; Q[2][1][1]=P[2][1][1];
    Q[0][2][0]=P[0][2][0]; Q[0][2][1]=P[0][2][1];
    Q[1][2][0]=P[1][2][0]; Q[1][2][1]=P[1][2][1];
    Q[2][2][0]=P[2][2][0]; Q[2][2][1]=P[2][2][1];
  }
  else {
    /* get the norm of (P_01 P_02) row vector */
    c=1.;
    factor = norm;
    a = fabs( P[0][1][0] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + (factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    a = fabs( P[0][1][1] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + c*(factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    beta = factor*sqrt(c); /* norm of (P_01 P_02) row vector */
    if( P[0][1][0]>0 ) {
      beta = -beta;
    }

#ifdef SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* SVD3x3_DEBUG */


    /* a=Re(P_01^+-beta), b=Im(P_01^+-beta) */
    a=P[0][1][0]-beta; b=-P[0][1][1];
    /* norm=sqrt(a^2+b^2) */
    c=fabs(a); d=fabs(b);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      norm = max;
    }
    else {
      c = min/max;
      norm = max*sqrt(1+c*c);
    }
    /* c=a/norm, d=b/norm */
    c=a/norm; d=b/norm;


    /* construct reflector (vector "v" for Householder transformation) */
    /* v_0=0 */
    v[0][0]=0.; v[0][1]=0.;
    /* v_1=1 */
    v[1][0]=1.; v[1][1]=0.;
    /* v_2=P_02^+/(P_01^+-beta)=P_02^+/(a+ib)=(P_02^+*(a-ib))/norm^2=(P_02^+/norm)*((a-ib)/norm)
          =(P_02^+/norm)*(c-id)=|a=Re(P_02^+)/norm,b=Im(P_02^+)/norm|=(a+ib)*(c-id)
          =(a*c+b*d)+i(b*c-a*d) */
    a=P[0][2][0]/norm; b=-P[0][2][1]/norm;
    v[2][0]=a*c+b*d;
    v[2][1]=b*c-a*d;

    nflops += 27;

#ifdef SVD3x3_DEBUG
for(i=0;i<3;i++) {
  printf("v[%d].re=%28.18e  v[%d].im=%28.18e\n",i,v[i][0],i,v[i][1]);
}
#endif /* SVD3x3_DEBUG */

    /* calcualate tau (coefficient for reflector) */
    taure=(beta-P[0][1][0])/beta;
    tauim=-P[0][1][1]/beta;

    /* assemble right unitary matrix V1=I-tau^+*v*v^+ (store in V1[3][3][2]) */
    V1[0][0][0]=1.;
    V1[0][0][1]=0.;
    V1[1][0][0]=0.;
    V1[1][0][1]=0.;
    V1[2][0][0]=0.;
    V1[2][0][1]=0.;
    V1[0][1][0]=0.;
    V1[0][1][1]=0.;
    V1[0][2][0]=0.;
    V1[0][2][1]=0.;
    /* V1_11=P_01^+/beta */
    V1[1][1][0]=P[0][1][0]/beta;
    V1[1][1][1]=-P[0][1][1]/beta;
    /* V1_21=P_02^+/beta */
    V1[2][1][0]=P[0][2][0]/beta;
    V1[2][1][1]=-P[0][2][1]/beta;
    /* V1_12=-tau^+*v_2^+=-(tau*v_2)^+ */
    V1[1][2][0]=-(taure*v[2][0]-tauim*v[2][1]);
    V1[1][2][1]=taure*v[2][1]+tauim*v[2][0];
    /* V1_22=1-tau^+*v_2*v_2^+ */
    a=v[2][0]*v[2][0]+v[2][1]*v[2][1];
    V1[2][2][0]=1-taure*a;
    V1[2][2][1]=tauim*a;


    /* apply the transformation to P matrix and store the result in Q
       Q=PV */
    Q[0][0][0]=P[0][0][0];
    Q[0][0][1]=0.;
    Q[1][0][0]=0.;
    Q[1][0][1]=0.;
    Q[2][0][0]=0.;
    Q[2][0][1]=0.;
    Q[0][1][0]=beta;
    Q[0][1][1]=0.;
    Q[0][2][0]=0.;
    Q[0][2][1]=0.;
    /* Q_11=P_11*V1_11+P_12*V_21 */
    Q[1][1][0]=P[1][1][0]*V1[1][1][0]-P[1][1][1]*V1[1][1][1]
              +P[1][2][0]*V1[2][1][0]-P[1][2][1]*V1[2][1][1];
    Q[1][1][1]=P[1][1][0]*V1[1][1][1]+P[1][1][1]*V1[1][1][0]
              +P[1][2][0]*V1[2][1][1]+P[1][2][1]*V1[2][1][0];
    /* Q_12=P_11*V1_12+P_12*V_22 */
    Q[1][2][0]=P[1][1][0]*V1[1][2][0]-P[1][1][1]*V1[1][2][1]
              +P[1][2][0]*V1[2][2][0]-P[1][2][1]*V1[2][2][1];
    Q[1][2][1]=P[1][1][0]*V1[1][2][1]+P[1][1][1]*V1[1][2][0]
              +P[1][2][0]*V1[2][2][1]+P[1][2][1]*V1[2][2][0];
    /* Q_21=P_21*V1_11+P_22*V_21 */
    Q[2][1][0]=P[2][1][0]*V1[1][1][0]-P[2][1][1]*V1[1][1][1]
              +P[2][2][0]*V1[2][1][0]-P[2][2][1]*V1[2][1][1];
    Q[2][1][1]=P[2][1][0]*V1[1][1][1]+P[2][1][1]*V1[1][1][0]
              +P[2][2][0]*V1[2][1][1]+P[2][2][1]*V1[2][1][0];
    /* Q_22=P_21*V1_12+P_22*V_22 */
    Q[2][2][0]=P[2][1][0]*V1[1][2][0]-P[2][1][1]*V1[1][2][1]
              +P[2][2][0]*V1[2][2][0]-P[2][2][1]*V1[2][2][1];
    Q[2][2][1]=P[2][1][0]*V1[1][2][1]+P[2][1][1]*V1[1][2][0]
              +P[2][2][0]*V1[2][2][1]+P[2][2][1]*V1[2][2][0];

    nflops += 15 + 7*8;
  }
#ifdef SVD3x3_DEBUG
printf("Right unitary matrix V1:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
      printf( "V1[%d][%d].re=%26.18e  V1[%d][%d].im=%26.18e\n",
              i, j, V1[i][j][0], i, j, V1[i][j][1] );
    }
#endif /* SVD3x3_DEBUG */



  /* *** Step 3: build second left reflector v,
                 calculate second left rotation U2,
                 apply to the matrix Q *** */
  /* calculate norm of ( Q[21] )
     with minimal loss of accuracy (similar to BLAS) */
  c=fabs(Q[2][1][0]); d=fabs(Q[2][1][1]);
  if( c>d ) {
    max=c; min=d;
  }
  else {
    max=d; min=c;
  }
  if( min==0 ) {
    norm = max;
  }
  else {
    c = min/max;
    norm = max*sqrt(1+c*c);
  }

  if( norm==0 && Q[1][1][1]==0 ) { /* no rotation needed */
#ifdef SVD3x3_DEBUG
    printf("Step 3: no rotation needed\n");
#endif /* SVD3x3_DEBUG */
    U2[0][0][0]=1.; U2[0][0][1]=0.;
    U2[0][1][0]=0.; U2[0][1][1]=0.;
    U2[0][2][0]=0.; U2[0][2][1]=0.;
    U2[1][0][0]=0.; U2[1][0][1]=0.;
    U2[1][1][0]=1.; U2[1][1][1]=0.;
    U2[1][2][0]=0.; U2[1][2][1]=0.;
    U2[2][0][0]=0.; U2[2][0][1]=0.;
    U2[2][1][0]=0.; U2[2][1][1]=0.;
    U2[2][2][0]=1.; U2[2][2][1]=0.;
    P[0][0][0]=Q[0][0][0]; P[0][0][1]=Q[0][0][1];
    P[1][0][0]=Q[1][0][0]; P[1][0][1]=Q[1][0][1];
    P[2][0][0]=Q[2][0][0]; P[2][0][1]=Q[2][0][1];
    P[0][1][0]=Q[0][1][0]; P[0][1][1]=Q[0][1][1];
    P[1][1][0]=Q[1][1][0]; P[1][1][1]=Q[1][1][1];
    P[2][1][0]=Q[2][1][0]; P[2][1][1]=Q[2][1][1];
    P[0][2][0]=Q[0][2][0]; P[0][2][1]=Q[0][2][1];
    P[1][2][0]=Q[1][2][0]; P[1][2][1]=Q[1][2][1];
    P[2][2][0]=Q[2][2][0]; P[2][2][1]=Q[2][2][1];
  }
  else {
    /* get the norm of (Q_11 Q_21) column vector */
    c=1.;
    factor = norm;
    a = fabs( Q[1][1][0] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + (factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    a = fabs( Q[1][1][1] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + c*(factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    beta = factor*sqrt(c); /* norm of (Q_11 Q_21) column vector */
    if( Q[1][1][0]>0 ) {
      beta = -beta;
    }

#ifdef SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* SVD3x3_DEBUG */


    /* a=Re(Q_11-beta), b=Im(Q_11-beta) */
    a=Q[1][1][0]-beta; b=Q[1][1][1];
    /* norm=sqrt(a^2+b^2) */
    c=fabs(a); d=fabs(b);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      norm = max;
    }
    else {
      c = min/max;
      norm = max*sqrt(1+c*c);
    }
    /* c=a/norm, d=b/norm */
    c=a/norm; d=b/norm;

    /* construct reflector (vector "v" for Householder transformation) */
    /* v_0=0 */
    v[0][0]=0.; v[0][1]=0.;
    /* v_1=1 */
    v[1][0]=1.; v[1][1]=0.;
    /* v_2=Q_21/(Q_11-beta)=Q_21/(a+ib)=(Q_21*(a-ib))/norm^2=(Q_21/norm)*((a-ib)/norm)
          =(Q_21/norm)*(c-id)=|a=Re(Q_21)/norm,b=Im(Q_21)/norm|=(a+ib)*(c-id)
          =(a*c+b*d)+i(b*c-a*d) */
    a=Q[2][1][0]/norm; b=Q[2][1][1]/norm;
    v[2][0]=a*c+b*d;
    v[2][1]=b*c-a*d;

    nflops += 27;
#ifdef SVD3x3_DEBUG
for(i=0;i<3;i++) {
  printf("v[%d].re=%28.18e  v[%d].im=%28.18e\n",i,v[i][0],i,v[i][1]);
}
#endif /* SVD3x3_DEBUG */


    /* calcualate tau (coefficient for reflector) */
    taure=(beta-Q[1][1][0])/beta;
    tauim=Q[1][1][1]/beta;


    /* assemble right unitary matrix U2=I-tau^+*v*v^+ (store in U2[3][3][2]) */
    U2[0][0][0]=1.;
    U2[0][0][1]=0.;
    U2[1][0][0]=0.;
    U2[1][0][1]=0.;
    U2[2][0][0]=0.;
    U2[2][0][1]=0.;
    U2[0][1][0]=0.;
    U2[0][1][1]=0.;
    U2[0][2][0]=0.;
    U2[0][2][1]=0.;
    /* U2_11=Q_11/beta */
    U2[1][1][0]=Q[1][1][0]/beta;
    U2[1][1][1]=Q[1][1][1]/beta;
    /* U2_21=Q_21/beta */
    U2[2][1][0]=Q[2][1][0]/beta;
    U2[2][1][1]=Q[2][1][1]/beta;
    /* U2_12=-tau^+*v_2^+=-(tau*v_2)^+ */
    U2[1][2][0]=-(taure*v[2][0]-tauim*v[2][1]);
    U2[1][2][1]=taure*v[2][1]+tauim*v[2][0];
    /* U2_22=1-tau^+*v_2*v_2^+ */
    a=v[2][0]*v[2][0]+v[2][1]*v[2][1];
    U2[2][2][0]=1-taure*a;
    U2[2][2][1]=tauim*a;
#ifdef SVD3x3_DEBUG
printf("Left unitary matrix U2:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
      printf( "U2[%d][%d].re=%26.18e  U2[%d][%d].im=%26.18e\n",
              i, j, U2[i][j][0], i, j, U2[i][j][1] );
    }
#endif /* SVD3x3_DEBUG */


    /* apply the transformation to Q matrix and store the result in P
       P=U^+Q */
    P[0][0][0]=Q[0][0][0];
    P[0][0][1]=0.;
    P[1][0][0]=0.;
    P[1][0][1]=0.;
    P[2][0][0]=0.;
    P[2][0][1]=0.;
    P[0][1][0]=Q[0][1][0];
    P[0][1][1]=0.;
    P[0][2][0]=0.;
    P[0][2][1]=0.;
    P[1][1][0]=beta;
    P[1][1][1]=0.;
    P[2][1][0]=0.;
    P[2][1][1]=0.;
    /* P_12=U2_11^+*Q_12+U2_21^+*Q_22 */
    P[1][2][0]=U2[1][1][0]*Q[1][2][0]+U2[1][1][1]*Q[1][2][1]
              +U2[2][1][0]*Q[2][2][0]+U2[2][1][1]*Q[2][2][1];
    P[1][2][1]=U2[1][1][0]*Q[1][2][1]-U2[1][1][1]*Q[1][2][0]
              +U2[2][1][0]*Q[2][2][1]-U2[2][1][1]*Q[2][2][0];
    /* P_22=U2_12^+*Q_12+U2_22^+*Q_22 */
    P[2][2][0]=U2[1][2][0]*Q[1][2][0]+U2[1][2][1]*Q[1][2][1]
              +U2[2][2][0]*Q[2][2][0]+U2[2][2][1]*Q[2][2][1];
    P[2][2][1]=U2[1][2][0]*Q[1][2][1]-U2[1][2][1]*Q[1][2][0]
              +U2[2][2][0]*Q[2][2][1]-U2[2][2][1]*Q[2][2][0];

    nflops += 15 + 7*8;

  }



  /* *** Step 4: build second right reflector v,
                 calculate second right rotation V2,
                 apply to the matrix P *** */
  if( P[1][2][1]==0 ) { /* no rotation needed */
#ifdef SVD3x3_DEBUG
    printf("Step 4: no rotation needed\n");
#endif /* SVD3x3_DEBUG */
    V2[0][0][0]=1.; V2[0][0][1]=0.;
    V2[0][1][0]=0.; V2[0][1][1]=0.;
    V2[0][2][0]=0.; V2[0][2][1]=0.;
    V2[1][0][0]=0.; V2[1][0][1]=0.;
    V2[1][1][0]=1.; V2[1][1][1]=0.;
    V2[1][2][0]=0.; V2[1][2][1]=0.;
    V2[2][0][0]=0.; V2[2][0][1]=0.;
    V2[2][1][0]=0.; V2[2][1][1]=0.;
    V2[2][2][0]=1.; V2[2][2][1]=0.;
    Q[0][0][0]=P[0][0][0]; Q[0][0][1]=P[0][0][1];
    Q[1][0][0]=P[1][0][0]; Q[1][0][1]=P[1][0][1];
    Q[2][0][0]=P[2][0][0]; Q[2][0][1]=P[2][0][1];
    Q[0][1][0]=P[0][1][0]; Q[0][1][1]=P[0][1][1];
    Q[1][1][0]=P[1][1][0]; Q[1][1][1]=P[1][1][1];
    Q[2][1][0]=P[2][1][0]; Q[2][1][1]=P[2][1][1];
    Q[0][2][0]=P[0][2][0]; Q[0][2][1]=P[0][2][1];
    Q[1][2][0]=P[1][2][0]; Q[1][2][1]=P[1][2][1];
    Q[2][2][0]=P[2][2][0]; Q[2][2][1]=P[2][2][1];
  }
  else {
    /* calculate norm of ( P[12] ) */
    c=fabs(P[1][2][0]); d=fabs(P[1][2][1]);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      beta = max;
    }
    else {
      c = min/max;
      beta = max*sqrt(1+c*c);
    }

    if( P[1][2][0]>0 ) {
      beta = -beta;
    }

#ifdef SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* SVD3x3_DEBUG */

    /* assemble right unitary matrix V1=I-tau^+*v*v^+ (store in V1[3][3][2]) */
    V2[0][0][0]=1.;
    V2[0][0][1]=0.;
    V2[1][0][0]=0.;
    V2[1][0][1]=0.;
    V2[2][0][0]=0.;
    V2[2][0][1]=0.;
    V2[0][1][0]=0.;
    V2[0][1][1]=0.;
    V2[0][2][0]=0.;
    V2[0][2][1]=0.;
    V2[1][1][0]=1.;
    V2[1][1][1]=0.;
    V2[2][1][0]=0.;
    V2[2][1][1]=0.;
    V2[1][2][0]=0.;
    V2[1][2][1]=0.;
    /* V2_22=1-tau^+*v_2*v_2^+=1-tau^+ */
    V2[2][2][0]=P[1][2][0]/beta;
    V2[2][2][1]=-P[1][2][1]/beta;
#ifdef SVD3x3_DEBUG
printf("Right unitary matrix V2:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
      printf( "V2[%d][%d].re=%26.18e  V2[%d][%d].im=%26.18e\n",
              i, j, V2[i][j][0], i, j, V2[i][j][1] );
    }
#endif /* SVD3x3_DEBUG */


    /* apply the transformation to P matrix and store the result in Q
       Q=PV */
    Q[0][0][0]=P[0][0][0];
    Q[0][0][1]=0.;
    Q[1][0][0]=0.;
    Q[1][0][1]=0.;
    Q[2][0][0]=0.;
    Q[2][0][1]=0.;
    Q[0][1][0]=P[0][1][0];
    Q[0][1][1]=0.;
    Q[0][2][0]=0.;
    Q[0][2][1]=0.;
    Q[1][1][0]=P[1][1][0];
    Q[1][1][1]=0.;
    Q[1][2][0]=beta;
    Q[1][2][1]=0.;
    Q[2][1][0]=0.;
    Q[2][1][1]=0.;
    /* Q_22=P_22*V2_22 */
    Q[2][2][0]=P[2][2][0]*V2[2][2][0]-P[2][2][1]*V2[2][2][1];
    Q[2][2][1]=P[2][2][0]*V2[2][2][1]+P[2][2][1]*V2[2][2][0];

    nflops += 12;
  }



  /* *** Step 5: build third left reflector v,
                 calculate third left rotation U3,
                 apply to the matrix P *** */
  if( Q[2][2][1]==0 ) { /* no rotation needed */
#ifdef SVD3x3_DEBUG
    printf("Step 5: no rotation needed\n");
#endif /* SVD3x3_DEBUG */
    U3[0][0][0]=1.; U3[0][0][1]=0.;
    U3[0][1][0]=0.; U3[0][1][1]=0.;
    U3[0][2][0]=0.; U3[0][2][1]=0.;
    U3[1][0][0]=0.; U3[1][0][1]=0.;
    U3[1][1][0]=1.; U3[1][1][1]=0.;
    U3[1][2][0]=0.; U3[1][2][1]=0.;
    U3[2][0][0]=0.; U3[2][0][1]=0.;
    U3[2][1][0]=0.; U3[2][1][1]=0.;
    U3[2][2][0]=1.; U3[2][2][1]=0.;
    P[0][0][0]=Q[0][0][0]; P[0][0][1]=Q[0][0][1];
    P[1][0][0]=Q[1][0][0]; P[1][0][1]=Q[1][0][1];
    P[2][0][0]=Q[2][0][0]; P[2][0][1]=Q[2][0][1];
    P[0][1][0]=Q[0][1][0]; P[0][1][1]=Q[0][1][1];
    P[1][1][0]=Q[1][1][0]; P[1][1][1]=Q[1][1][1];
    P[2][1][0]=Q[2][1][0]; P[2][1][1]=Q[2][1][1];
    P[0][2][0]=Q[0][2][0]; P[0][2][1]=Q[0][2][1];
    P[1][2][0]=Q[1][2][0]; P[1][2][1]=Q[1][2][1];
    P[2][2][0]=Q[2][2][0]; P[2][2][1]=Q[2][2][1];
  }
  else {
    /* calculate norm of ( Q[22] ) */
    c=fabs(Q[2][2][0]); d=fabs(Q[2][2][1]);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      beta = max;
    }
    else {
      c = min/max;
      beta = max*sqrt(1+c*c);
    }

    if( Q[2][2][0]>0 ) {
      beta = -beta;
    }

#ifdef SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* SVD3x3_DEBUG */

    /* assemble left unitary matrix U3=I-tau^+*v*v^+ (store in U3[3][3][2]) */
    U3[0][0][0]=1.;
    U3[0][0][1]=0.;
    U3[1][0][0]=0.;
    U3[1][0][1]=0.;
    U3[2][0][0]=0.;
    U3[2][0][1]=0.;
    U3[0][1][0]=0.;
    U3[0][1][1]=0.;
    U3[0][2][0]=0.;
    U3[0][2][1]=0.;
    U3[1][1][0]=1.;
    U3[1][1][1]=0.;
    U3[2][1][0]=0.;
    U3[2][1][1]=0.;
    U3[1][2][0]=0.;
    U3[1][2][1]=0.;
    /* U3_22=1-tau^+*v_2*v_2^+=1-tau^+ */
    U3[2][2][0]=Q[2][2][0]/beta;
    U3[2][2][1]=Q[2][2][1]/beta;
#ifdef SVD3x3_DEBUG
printf("Left unitary matrix U3:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
      printf( "U3[%d][%d].re=%26.18e  U3[%d][%d].im=%26.18e\n",
              i, j, U3[i][j][0], i, j, U3[i][j][1] );
    }
#endif /* SVD3x3_DEBUG */


    /* apply the transformation to Q matrix and store the result in P
       P=U^+Q */
    P[0][0][0]=Q[0][0][0];
    P[0][0][1]=0.;
    P[1][0][0]=0.;
    P[1][0][1]=0.;
    P[2][0][0]=0.;
    P[2][0][1]=0.;
    P[0][1][0]=Q[0][1][0];
    P[0][1][1]=0.;
    P[0][2][0]=0.;
    P[0][2][1]=0.;
    P[1][1][0]=Q[1][1][0];
    P[1][1][1]=0.;
    P[1][2][0]=Q[1][2][0];
    P[1][2][1]=0.;
    P[2][1][0]=0.;
    P[2][1][1]=0.;
    P[2][2][0]=beta;
    P[2][2][1]=0.;

    nflops += 6;

  }




  /* *** This part starts with a bidiagonal matrix and uses
         QR algorithm with shifts to eliminate the superdiagonal *** */
  /* prepare left and right real orthogonal matrices that
     accumulate Givens rotations from QR algorithm */
  UO3[0][0]=1.; UO3[0][1]=0.; UO3[0][2]=0.;
  UO3[1][0]=0.; UO3[1][1]=1.; UO3[1][2]=0.;
  UO3[2][0]=0.; UO3[2][1]=0.; UO3[2][2]=1.;
  VO3[0][0]=1.; VO3[0][1]=0.; VO3[0][2]=0.;
  VO3[1][0]=0.; VO3[1][1]=1.; VO3[1][2]=0.;
  VO3[2][0]=0.; VO3[2][1]=0.; VO3[2][2]=1.;

  iter=0;

#ifdef SVD3x3_DEBUG
printf( "QR iteration: %d\n", iter );
printf( "%+20.16e %+20.16e %+20.16e\n", b00, b01, b02 );
printf( "%+20.16e %+20.16e %+20.16e\n", b10, b11, b12 );
printf( "%+20.16e %+20.16e %+20.16e\n", b20, b21, b22 );
#endif /* SVD3x3_DEBUG */

  do {

    iter++;
    if(iter>300) return 1;

    /* chop small superdiagonal elements */
    if( fabs(b01) < SVD3x3_PREC*(fabs(b00)+fabs(b11)) ) {
      b01=0;
    }
    if( fabs(b12) < SVD3x3_PREC*(fabs(b00)+fabs(b22)) ) {
      b12=0;
    }

    nflops += 4;

    /* Cases:
       b01=b12=0 -- matrix is already diagonalized,
       b01=0 -- need to work with 2x2 lower block,
       b12=0 -- need to work with 2x2 upper block,
       else -- normal iteration */
    if( !(b01==0 && b12==0) ) {
      if( b01==0 ) {
#ifdef SVD3x3_DEBUG
	printf( "Entering case b01==0\n" );
#endif /* SVD3x3_DEBUG */
        /* need to diagonalize 2x2 lower block */
	svd2x2bidiag( &b11, &b12, &b22, UO2, VO2, &nflops );

        /* multiply left UO3 matrix */
        for(i=0;i<3;i++) {
          a=UO3[i][1]; b=UO3[i][2];
          UO3[i][1]=a*UO2[0][0]+b*UO2[1][0];
          UO3[i][2]=a*UO2[0][1]+b*UO2[1][1];
        }
        /* multiply right VO3 matrix */
        for(i=0;i<3;i++) {
          a=VO3[i][1]; b=VO3[i][2];
          VO3[i][1]=a*VO2[0][0]+b*VO2[1][0];
          VO3[i][2]=a*VO2[0][1]+b*VO2[1][1];
        }

	nflops += 36;

      }
      else {
        if( b12==0 ) {
#ifdef SVD3x3_DEBUG
	  printf( "Entering case b12==0\n" );
#endif /* SVD3x3_DEBUG */
          /* need to diagonalize 2x2 upper block */
	  svd2x2bidiag( &b00, &b01, &b11, UO2, VO2, &nflops );

          /* multiply left UO3 matrix */
          for(i=0;i<3;i++) {
            a=UO3[i][0]; b=UO3[i][1];
            UO3[i][0]=a*UO2[0][0]+b*UO2[1][0];
            UO3[i][1]=a*UO2[0][1]+b*UO2[1][1];
          }
          /* multiply right VO3 matrix */
          for(i=0;i<3;i++) {
            a=VO3[i][0]; b=VO3[i][1];
            VO3[i][0]=a*VO2[0][0]+b*VO2[1][0];
            VO3[i][1]=a*VO2[0][1]+b*VO2[1][1];
          }

	  nflops += 36;
        }
        else {
          /* normal 3x3 iteration */

          /* QR shift does not work if there are zeros
             on the diagonal, therefore first check
             for special cases: b00==0 or b11==0 or b22==0 */

          if( b00==0 ) {
#ifdef SVD3x3_DEBUG
printf( "Entering case b00==0\n" );
#endif /* SVD3x3_DEBUG */
            /* b01 can be rotated away to create b02,
               and then b02 can be rotated away
               (both are left rotations) */
            if( fabs(b01)>fabs(b11) ) {
              cotphi=b11/b01;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b01/b11;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][0]; b=UO3[i][1];
              UO3[i][0]=a*cosphi-b*sinphi;
              UO3[i][1]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this generates b02 */
            b11=b01*sinphi+b11*cosphi;
            b02=-b12*sinphi;
            b12=b12*cosphi;
            b01=0.;
            if( fabs(b02)>fabs(b22) ) {
              cotphi=b22/b02;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b02/b22;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][0]; b=UO3[i][2];
              UO3[i][0]=a*cosphi-b*sinphi;
              UO3[i][2]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix */
            b22=b02*sinphi+b22*cosphi;
            b02=0.;

	    nflops += 56;
          }
          else if( b11==0 ) {
#ifdef SVD3x3_DEBUG
printf( "Entering case b11==0\n" );
#endif /* SVD3x3_DEBUG */
            /* b12 is rotated away with left rotation */
            if( fabs(b12)>fabs(b22) ) {
              cotphi=b22/b12;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b12/b22;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][1]; b=UO3[i][2];
              UO3[i][1]=a*cosphi-b*sinphi;
              UO3[i][2]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix */
            b22=b12*sinphi+b22*cosphi;
            b12=0.;

	    nflops += 27;
          }
          else if( b22==0 ) {
#ifdef SVD3x3_DEBUG
printf( "Entering case b22==0\n" );
#endif /* SVD3x3_DEBUG */
            /* b12 is rotated away and b02 appears,
               then b02 is rotated away, both are
               right rotations */
            if( fabs(b12)>fabs(b11) ) {
              cotphi=b11/b12;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b12/b11;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply right VO3 matrix */
            for(i=0;i<3;i++) {
              a=VO3[i][1]; b=VO3[i][2];
              VO3[i][1]= a*cosphi+b*sinphi;
              VO3[i][2]=-a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix */
            b02=-b01*sinphi;
            b01=b01*cosphi;
            b11=b11*cosphi+b12*sinphi;
            b12=0.;
            /* second rotation removes b02 */
            if( fabs(b02)>fabs(b00) ) {
              cotphi=b00/b02;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b02/b00;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply right VO3 matrix */
            for(i=0;i<3;i++) {
              a=VO3[i][0]; b=VO3[i][2];
              VO3[i][0]= a*cosphi+b*sinphi;
              VO3[i][2]=-a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix */
            b00=b00*cosphi+b02*sinphi;
            b02=0.;
	    
	    nflops += 64;
          }
          else {
            /* full iteration with QR shift */
#ifdef SVD3x3_DEBUG
printf( "Entering case of normal QR iteration\n" );
#endif /* SVD3x3_DEBUG */

            /* find max eigenvalue of bottom 2x2 minor */
            m11=b11*b11+b01*b01;
            m22=b22*b22+b12*b12;
            m12=b11*b12;
            dm=(m11-m22)/2;

            /* safely calculate sqrt */
            c=fabs(dm); d=fabs(m12);
            if( c>d ) {
              max=c; min=d;
            }
            else {
              max=d; min=c;
            }
            if( min==0 ) {
              norm = max;
            }
            else {
              c = min/max;
              norm = max*sqrt(1+c*c);
            }

            if( dm>=0 ) {
              lambdamax=m22-(m12*m12)/(dm+norm);
            }
            else {
              lambdamax=m22+(m12*m12)/(norm-dm);
            }

            /* calculate first Givens rotation (on the right) */
            a=b00*b00-lambdamax;
            b=b00*b01;
            if( 0==b ) {
              cosphi=1.;
              sinphi=0.;
            }
            else {
              if( fabs(b)>fabs(a) ) {
                cotphi=-a/b;
                sinphi=1./sqrt(1+cotphi*cotphi);
                cosphi=cotphi*sinphi;
              }
              else {
                tanphi=-b/a;
                cosphi=1./sqrt(1+tanphi*tanphi);
                sinphi=tanphi*cosphi;
              }
	      nflops += 7;
            }
            /* multiply right VO3 matrix */
            for(i=0;i<3;i++) {
              a=VO3[i][0]; b=VO3[i][1];
              VO3[i][0]=a*cosphi-b*sinphi;
              VO3[i][1]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this generate b10 */
            a=b00; b=b01;
            b00=a*cosphi-b*sinphi;
            b01=a*sinphi+b*cosphi;
            b10=-b11*sinphi;
            b11=b11*cosphi; 

            /* calculate second Givens rotation (on the left) */
            if(0==b10) {
              cosphi=1.;
              sinphi=0.;
            }
            else {
              if( fabs(b10)>fabs(b00) ) {
                cotphi=-b00/b10;
                sinphi=1/sqrt(1+cotphi*cotphi);
                cosphi=cotphi*sinphi;
              }
              else {
                tanphi=-b10/b00;
                cosphi=1/sqrt(1+tanphi*tanphi);
                sinphi=tanphi*cosphi;
              }

	      nflops += 7;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][0]; b=UO3[i][1];
              UO3[i][0]= a*cosphi-b*sinphi;
              UO3[i][1]= a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this generates b02 */
            b00=b00*cosphi-b10*sinphi;
            a=b01; b=b11;
            b01=a*cosphi-b*sinphi;
            b11=a*sinphi+b*cosphi;
            b02=-b12*sinphi;
            b12=b12*cosphi;
            b10=0.;

            /* calculate third Givens rotation (on the right) */
            if(0==b02) {
              cosphi=1.;
              sinphi=0.;
            }
            else {
              if( fabs(b02)>fabs(b01) ) {
                cotphi=-b01/b02;
                sinphi=1/sqrt(1+cotphi*cotphi);
                cosphi=cotphi*sinphi;
              }
              else {
                tanphi=-b02/b01;
                cosphi=1/sqrt(1+tanphi*tanphi);
                sinphi=tanphi*cosphi;
              }

	      nflops += 7;
            }
            /* multiply right VO3 matrix */
            for(i=0;i<3;i++) {
              a=VO3[i][1]; b=VO3[i][2];
              VO3[i][1]=a*cosphi-b*sinphi;
              VO3[i][2]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this generates b21 */
            b01=b01*cosphi-b02*sinphi;
            a=b11; b=b12;
            b11=a*cosphi-b*sinphi;
            b12=a*sinphi+b*cosphi;
            b21=-b22*sinphi;
            b22=b22*cosphi;
            b02=0.;

            /* calculate fourth Givens rotation (on the left) */
            if(0==b21) {
              cosphi=1.;
              sinphi=0.;
            }
            else {
              if( fabs(b21)>fabs(b11) ) {
                cotphi=-b11/b21;
                sinphi=1/sqrt(1+cotphi*cotphi);
                cosphi=cotphi*sinphi;
              }
              else {
                tanphi=-b21/b11;
                cosphi=1/sqrt(1+tanphi*tanphi);
                sinphi=tanphi*cosphi;
              }

	      nflops += 7;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][1]; b=UO3[i][2];
              UO3[i][1]= a*cosphi-b*sinphi;
              UO3[i][2]= a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this eliminates b21 */
            b11=b11*cosphi-b21*sinphi;
            a=b12; b=b22;
            b12=a*cosphi-b*sinphi;
            b22=a*sinphi+b*cosphi;
            b21=0.;

	    nflops += 127;
          }
        } /* end of normal 3x3 iteration */
      }
    }
#ifdef SVD3x3_DEBUG
printf( "QR iteration: %d\n", iter );
printf( "%+20.16e %+20.16e %+20.16e\n", b00, b01, b02 );
printf( "%+20.16e %+20.16e %+20.16e\n", b10, b11, b12 );
printf( "%+20.16e %+20.16e %+20.16e\n", b20, b21, b22 );
#endif /* SVD3x3_DEBUG */
  }
  while( b01!=0 || b12!=0 );


  /* make singular values positive */
  if(b00<0) {
    b00=-b00;
    VO3[0][0]=-VO3[0][0];
    VO3[1][0]=-VO3[1][0];
    VO3[2][0]=-VO3[2][0];
  }
  if(b11<0) {
    b11=-b11;
    VO3[0][1]=-VO3[0][1];
    VO3[1][1]=-VO3[1][1];
    VO3[2][1]=-VO3[2][1];
  }
  if(b22<0) {
    b22=-b22;
    VO3[0][2]=-VO3[0][2];
    VO3[1][2]=-VO3[1][2];
    VO3[2][2]=-VO3[2][2];
  }



  /* Q=U1*U2 (U2 is block diagonal with U2_00=1) */
  Q[0][0][0]=U1[0][0][0]; Q[0][0][1]=U1[0][0][1];
  Q[1][0][0]=U1[1][0][0]; Q[1][0][1]=U1[1][0][1];
  Q[2][0][0]=U1[2][0][0]; Q[2][0][1]=U1[2][0][1];
  /* Q_01=U1_01*U2_11+U1_02*U2_21 */
  Q[0][1][0]=U1[0][1][0]*U2[1][1][0]-U1[0][1][1]*U2[1][1][1]
            +U1[0][2][0]*U2[2][1][0]-U1[0][2][1]*U2[2][1][1];
  Q[0][1][1]=U1[0][1][0]*U2[1][1][1]+U1[0][1][1]*U2[1][1][0]
            +U1[0][2][0]*U2[2][1][1]+U1[0][2][1]*U2[2][1][0];
  /* Q_02=U1_01*U2_12+U1_02*U2_22 */
  Q[0][2][0]=U1[0][1][0]*U2[1][2][0]-U1[0][1][1]*U2[1][2][1]
            +U1[0][2][0]*U2[2][2][0]-U1[0][2][1]*U2[2][2][1];
  Q[0][2][1]=U1[0][1][0]*U2[1][2][1]+U1[0][1][1]*U2[1][2][0]
            +U1[0][2][0]*U2[2][2][1]+U1[0][2][1]*U2[2][2][0];
  /* Q_11=U1_11*U2_11+U1_12*U2_21 */
  Q[1][1][0]=U1[1][1][0]*U2[1][1][0]-U1[1][1][1]*U2[1][1][1]
            +U1[1][2][0]*U2[2][1][0]-U1[1][2][1]*U2[2][1][1];
  Q[1][1][1]=U1[1][1][0]*U2[1][1][1]+U1[1][1][1]*U2[1][1][0]
            +U1[1][2][0]*U2[2][1][1]+U1[1][2][1]*U2[2][1][0];
  /* Q_12=U1_11*U2_12+U1_12*U2_22 */
  Q[1][2][0]=U1[1][1][0]*U2[1][2][0]-U1[1][1][1]*U2[1][2][1]
            +U1[1][2][0]*U2[2][2][0]-U1[1][2][1]*U2[2][2][1];
  Q[1][2][1]=U1[1][1][0]*U2[1][2][1]+U1[1][1][1]*U2[1][2][0]
            +U1[1][2][0]*U2[2][2][1]+U1[1][2][1]*U2[2][2][0];
  /* Q_21=U1_21*U2_11+U1_22*U2_21 */
  Q[2][1][0]=U1[2][1][0]*U2[1][1][0]-U1[2][1][1]*U2[1][1][1]
            +U1[2][2][0]*U2[2][1][0]-U1[2][2][1]*U2[2][1][1];
  Q[2][1][1]=U1[2][1][0]*U2[1][1][1]+U1[2][1][1]*U2[1][1][0]
            +U1[2][2][0]*U2[2][1][1]+U1[2][2][1]*U2[2][1][0];
  /* Q_22=U1_21*U2_12+U1_22*U2_22 */
  Q[2][2][0]=U1[2][1][0]*U2[1][2][0]-U1[2][1][1]*U2[1][2][1]
            +U1[2][2][0]*U2[2][2][0]-U1[2][2][1]*U2[2][2][1];
  Q[2][2][1]=U1[2][1][0]*U2[1][2][1]+U1[2][1][1]*U2[1][2][0]
            +U1[2][2][0]*U2[2][2][1]+U1[2][2][1]*U2[2][2][0];

  /* Q=Q*U3 (U3 is block diagonal with U3_00=1, U3_11=1)
     (this changes only third column of Q */
  a=Q[0][2][0]*U3[2][2][0]-Q[0][2][1]*U3[2][2][1];
  b=Q[0][2][0]*U3[2][2][1]+Q[0][2][1]*U3[2][2][0];
  Q[0][2][0]=a; Q[0][2][1]=b;
  a=Q[1][2][0]*U3[2][2][0]-Q[1][2][1]*U3[2][2][1];
  b=Q[1][2][0]*U3[2][2][1]+Q[1][2][1]*U3[2][2][0];
  Q[1][2][0]=a; Q[1][2][1]=b;
  a=Q[2][2][0]*U3[2][2][0]-Q[2][2][1]*U3[2][2][1];
  b=Q[2][2][0]*U3[2][2][1]+Q[2][2][1]*U3[2][2][0];
  Q[2][2][0]=a; Q[2][2][1]=b;

  nflops += 102;

  /* final U=Q*UO3
     (unitary times orthogonal that accumulated Givens rotations) */
  U00re=Q[0][0][0]*UO3[0][0]+Q[0][1][0]*UO3[1][0]+Q[0][2][0]*UO3[2][0];
  U00im=Q[0][0][1]*UO3[0][0]+Q[0][1][1]*UO3[1][0]+Q[0][2][1]*UO3[2][0];
  U01re=Q[0][0][0]*UO3[0][1]+Q[0][1][0]*UO3[1][1]+Q[0][2][0]*UO3[2][1];
  U01im=Q[0][0][1]*UO3[0][1]+Q[0][1][1]*UO3[1][1]+Q[0][2][1]*UO3[2][1];
  U02re=Q[0][0][0]*UO3[0][2]+Q[0][1][0]*UO3[1][2]+Q[0][2][0]*UO3[2][2];
  U02im=Q[0][0][1]*UO3[0][2]+Q[0][1][1]*UO3[1][2]+Q[0][2][1]*UO3[2][2];
  U10re=Q[1][0][0]*UO3[0][0]+Q[1][1][0]*UO3[1][0]+Q[1][2][0]*UO3[2][0];
  U10im=Q[1][0][1]*UO3[0][0]+Q[1][1][1]*UO3[1][0]+Q[1][2][1]*UO3[2][0];
  U11re=Q[1][0][0]*UO3[0][1]+Q[1][1][0]*UO3[1][1]+Q[1][2][0]*UO3[2][1];
  U11im=Q[1][0][1]*UO3[0][1]+Q[1][1][1]*UO3[1][1]+Q[1][2][1]*UO3[2][1];
  U12re=Q[1][0][0]*UO3[0][2]+Q[1][1][0]*UO3[1][2]+Q[1][2][0]*UO3[2][2];
  U12im=Q[1][0][1]*UO3[0][2]+Q[1][1][1]*UO3[1][2]+Q[1][2][1]*UO3[2][2];
  U20re=Q[2][0][0]*UO3[0][0]+Q[2][1][0]*UO3[1][0]+Q[2][2][0]*UO3[2][0];
  U20im=Q[2][0][1]*UO3[0][0]+Q[2][1][1]*UO3[1][0]+Q[2][2][1]*UO3[2][0];
  U21re=Q[2][0][0]*UO3[0][1]+Q[2][1][0]*UO3[1][1]+Q[2][2][0]*UO3[2][1];
  U21im=Q[2][0][1]*UO3[0][1]+Q[2][1][1]*UO3[1][1]+Q[2][2][1]*UO3[2][1];
  U22re=Q[2][0][0]*UO3[0][2]+Q[2][1][0]*UO3[1][2]+Q[2][2][0]*UO3[2][2];
  U22im=Q[2][0][1]*UO3[0][2]+Q[2][1][1]*UO3[1][2]+Q[2][2][1]*UO3[2][2];

  nflops += 90;

  /* Q=V1*V2 (V1 is block diagonal with V2_11=1,
              V2 is block diagonal with V2_11=1, V2_22=1) */
  Q[0][0][0]=V1[0][0][0]; Q[0][0][1]=V1[0][0][1];
  Q[1][0][0]=V1[1][0][0]; Q[1][0][1]=V1[1][0][1];
  Q[2][0][0]=V1[2][0][0]; Q[2][0][1]=V1[2][0][1];
  Q[0][1][0]=V1[0][1][0]; Q[0][1][1]=V1[0][1][1];
  Q[0][2][0]=V1[0][2][0]; Q[0][2][1]=V1[0][2][1];
  Q[1][1][0]=V1[1][1][0]; Q[1][1][1]=V1[1][1][1];
  Q[2][1][0]=V1[2][1][0]; Q[2][1][1]=V1[2][1][1];
  Q[1][2][0]=V1[1][2][0]*V2[2][2][0]-V1[1][2][1]*V2[2][2][1];
  Q[1][2][1]=V1[1][2][0]*V2[2][2][1]+V1[1][2][1]*V2[2][2][0];
  Q[2][2][0]=V1[2][2][0]*V2[2][2][0]-V1[2][2][1]*V2[2][2][1];
  Q[2][2][1]=V1[2][2][0]*V2[2][2][1]+V1[2][2][1]*V2[2][2][0];

  /* final V=Q*VO3
     (unitary times orthogonal that accumulated Givens rotations) */
  V00re=Q[0][0][0]*VO3[0][0]+Q[0][1][0]*VO3[1][0]+Q[0][2][0]*VO3[2][0];
  V00im=Q[0][0][1]*VO3[0][0]+Q[0][1][1]*VO3[1][0]+Q[0][2][1]*VO3[2][0];
  V01re=Q[0][0][0]*VO3[0][1]+Q[0][1][0]*VO3[1][1]+Q[0][2][0]*VO3[2][1];
  V01im=Q[0][0][1]*VO3[0][1]+Q[0][1][1]*VO3[1][1]+Q[0][2][1]*VO3[2][1];
  V02re=Q[0][0][0]*VO3[0][2]+Q[0][1][0]*VO3[1][2]+Q[0][2][0]*VO3[2][2];
  V02im=Q[0][0][1]*VO3[0][2]+Q[0][1][1]*VO3[1][2]+Q[0][2][1]*VO3[2][2];
  V10re=Q[1][0][0]*VO3[0][0]+Q[1][1][0]*VO3[1][0]+Q[1][2][0]*VO3[2][0];
  V10im=Q[1][0][1]*VO3[0][0]+Q[1][1][1]*VO3[1][0]+Q[1][2][1]*VO3[2][0];
  V11re=Q[1][0][0]*VO3[0][1]+Q[1][1][0]*VO3[1][1]+Q[1][2][0]*VO3[2][1];
  V11im=Q[1][0][1]*VO3[0][1]+Q[1][1][1]*VO3[1][1]+Q[1][2][1]*VO3[2][1];
  V12re=Q[1][0][0]*VO3[0][2]+Q[1][1][0]*VO3[1][2]+Q[1][2][0]*VO3[2][2];
  V12im=Q[1][0][1]*VO3[0][2]+Q[1][1][1]*VO3[1][2]+Q[1][2][1]*VO3[2][2];
  V20re=Q[2][0][0]*VO3[0][0]+Q[2][1][0]*VO3[1][0]+Q[2][2][0]*VO3[2][0];
  V20im=Q[2][0][1]*VO3[0][0]+Q[2][1][1]*VO3[1][0]+Q[2][2][1]*VO3[2][0];
  V21re=Q[2][0][0]*VO3[0][1]+Q[2][1][0]*VO3[1][1]+Q[2][2][0]*VO3[2][1];
  V21im=Q[2][0][1]*VO3[0][1]+Q[2][1][1]*VO3[1][1]+Q[2][2][1]*VO3[2][1];
  V22re=Q[2][0][0]*VO3[0][2]+Q[2][1][0]*VO3[1][2]+Q[2][2][0]*VO3[2][2];
  V22im=Q[2][0][1]*VO3[0][2]+Q[2][1][1]*VO3[1][2]+Q[2][2][1]*VO3[2][2];

  nflops += 102;

  /* singular values */
  sigma[0]=b00; sigma[1]=b11; sigma[2]=b22;

  *nf += nflops;

  return 0;
}


/* SVD of 2x2 real matrix brought to the form:
    [ a00 a01]
    [   0 a11]
   This routine eliminates off-diagonal element, handling special cases */
int svd2x2bidiag(double *a00, double *a01, double *a11, 
		 double U2[2][2], double V2[2][2], size_t *nf) {
  register double sinphi, cosphi, tanphi, cotphi;
  register double a, b, min, max, abs00, abs01, abs11;
  register double lna01a11, lna00, ln_num, tau, t;
  register double P00, P01, P10, P11;
  register int isign;
  size_t nflops = 0;

  U2[0][0]=1.; U2[0][1]=0.;
  U2[1][0]=0.; U2[1][1]=1.;
  V2[0][0]=1.; V2[0][1]=0.;
  V2[1][0]=0.; V2[1][1]=1.;

  if( *a00==0 ) {
    if( *a11==0 ) {
      cosphi=1.;
      sinphi=0.;
    }
    else {
      if( fabs(*a11)>fabs(*a01) ) {
        cotphi=-(*a01)/(*a11);
        sinphi=1/sqrt(1+cotphi*cotphi);
        cosphi=cotphi*sinphi;
      }
      else {
        tanphi=-(*a11)/(*a01);
        cosphi=1/sqrt(1+tanphi*tanphi);
        sinphi=tanphi*cosphi;
      }
      nflops += 6;
    }
    /* multiply matrix A */
    (*a00)=cosphi*(*a01)-sinphi*(*a11);
    (*a01)=0.; (*a11)=0.;
    /* exchange columns in matrix V */
    V2[0][0]=0.; V2[0][1]=1.;
    V2[1][0]=1.; V2[1][1]=0.;
    /* U is just Givens rotation */
    U2[0][0]= cosphi; U2[0][1]= sinphi;
    U2[1][0]=-sinphi; U2[1][1]= cosphi;

    nflops += 3;
  }
  else if( *a11==0 ) {
    if( *a01==0 ) {
      cosphi=1.;
      sinphi=0.;
    }
    else {
      if( fabs(*a01)>fabs(*a00) ) {
        cotphi=-(*a00)/(*a01);
        sinphi=1/sqrt(1+cotphi*cotphi);
        cosphi=cotphi*sinphi;
      }
      else {
        tanphi=-(*a01)/(*a00);
        cosphi=1/sqrt(1+tanphi*tanphi);
        sinphi=tanphi*cosphi;
      }
      nflops += 7;
    }
    /* multiply matrix A */
    (*a00)=cosphi*(*a00)-sinphi*(*a01);
    (*a01)=0.; (*a11)=0.;
    /* V is just Givens rotation */
    V2[0][0]= cosphi; V2[0][1]= sinphi;
    V2[1][0]=-sinphi; V2[1][1]= cosphi;
    nflops += 3;
  }
  else if( *a01==0 ){ /* nothing to be done */
    ;
  }
  else {
    /* need to calculate ( a11^2+a01^2-a00^2 )/( 2*a00*a01 )
       avoiding overflow/underflow,
       use logarithmic coding */
    abs01=fabs(*a01); abs11=fabs(*a11);
    if(abs01>abs11) {
      min=abs11; max=abs01;
    }
    else {
      min=abs01; max=abs11;
    }
    a=min/max;
    lna01a11=2*log(max)+log(1+a*a);

    abs00=fabs(*a00);
    lna00=2*log(abs00);
    if( lna01a11>lna00 ) {
      /* subtract smaller from larger, overall "+" */
      isign=1;
      ln_num=lna01a11+log(1.-exp(lna00-lna01a11));
    }
    else {
      /* subtract larger from smaller, need to change order, overall "-" */
      isign=-1;
      ln_num=lna00+log(1.-exp(lna01a11-lna00));
    }
    a=ln_num-log(2)-log(abs00)-log(abs01);
    tau=exp(a);
    tau*=isign;
    if(*a00<0.)
      {
        tau*=-1;
      }
    if(*a01<0.)
      {
        tau*=-1;
      }

    /* calculate b=sqrt(1+tau^2) */
    a=fabs(tau);
    if( a>1. ) {
      max=a; min=1.;
    }
    else {
      max=1.; min=a;
    }
    if( min==0 ) {
      b = max;
    }
    else {
      a = min/max;
      b = max*sqrt(1+a*a);
    }
    if(tau>=0.) {
      t = 1.0/(tau + b);
    }
    else {
      t = 1.0/(tau - b);
    }

    /* calculate b=sqrt(1+t^2) */
    a=fabs(t);
    if( a>1. ) {
      max=a; min=1.;
    }
    else {
      max=1.; min=a;
    }
    if( min==0 ) {
      b = max;
    }
    else {
      a = min/max;
      b = max*sqrt(1+a*a);
    }
    cosphi=1./b;
    sinphi=t*cosphi;

    /* transform matrix A so it has othogonal columns */
    P00= cosphi*(*a00)-sinphi*(*a01);
    P10=-sinphi*(*a11);
    P01= sinphi*(*a00)+cosphi*(*a01);
    P11= cosphi*(*a11);

    /* prepare V  */
    V2[0][0]= cosphi; V2[0][1]= sinphi;
    V2[1][0]=-sinphi; V2[1][1]= cosphi;

    /* make column with the largest norm first column */
    if( sqrt(P00*P00+P10*P10)<sqrt(P01*P01+P11*P11) ) {
      a=P00; P00=P01; P01=a;
      a=P10; P10=P11; P11=a;
      /* exchange columns in matrix V2 */
      a=V2[0][0]; V2[0][0]=V2[0][1]; V2[0][1]=a;
      a=V2[1][0]; V2[1][0]=V2[1][1]; V2[1][1]=a;
    }

    /* calculate left Givens rotation and diagonalize */
    if( P10==0 ) {
      cosphi=1.;
      sinphi=0.;
    }
    else {
      if( fabs(P10)>fabs(P00) ) {
        cotphi=-P00/P10;
        sinphi=1/sqrt(1+cotphi*cotphi);
        cosphi=cotphi*sinphi;
      }
      else {
        tanphi=-P10/P00;
        cosphi=1/sqrt(1+tanphi*tanphi);
        sinphi=tanphi*cosphi;
      }
      nflops += 7;
    }
    *a00=P00*cosphi-P10*sinphi;
    *a01=0.;
    *a11=P01*sinphi+P11*cosphi;

    /* U is just Givens rotation */
    U2[0][0]= cosphi; U2[0][1]= sinphi;
    U2[1][0]=-sinphi; U2[1][1]= cosphi;

    nflops += 56;
  }

  *nf += nflops;

  return 0;
}


void show_su3_mat_opts( void ){

#ifdef MILC_GLOBAL_DEBUG
    node0_printf("MILC_GLOBAL_DEBUG ***********************\n");
#endif

#ifdef HISQ_REUNITARIZATION_DEBUG
  node0_printf("HISQ_REUNITARIZATION_DEBUG is ON\n");
#endif
}
