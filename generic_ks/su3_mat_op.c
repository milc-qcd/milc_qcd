/* Temporary file with additional 3x3 complex matrix routines
   (arbitrary, not necessarily SU(3) matrices).
   AB 071005 */


#include <stdio.h>
#include <math.h>
#include "generic_ks_includes.h"      /* definitions files and prototypes */
#include "../include/su3_mat_op.h"



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
void su3_root_inv( su3_matrix *a, su3_matrix *x, su3_matrix *y) {
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

  } while( (iter<SU3_ROOT_INV_MAX_ITER) && 
           (norm_x>SU3_ROOT_INV_NORM_EPS) &&
           (norm_y>SU3_ROOT_INV_NORM_EPS) );

}


/* Unitarize 3x3 complex matrix:
   B=A*(A^+ A)^-1/2
   B is a U(3) matrix but not an SU(3) matrix(!) */
void su3_unitarize( su3_matrix *a, su3_matrix *b ) {
  su3_matrix a2, x, y;

  /* X=A^+ A */
  mult_su3_an( a, a, &a2 );

  /* X=A2^1/2, Y=A2^-1/2 */
  su3_root_inv( &a2, &x, &y );

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
#ifdef HISQ_REUNITARIZATION_DEBUG
   /* store SU(3) matrices */
  if( lattice[index_site].on_step_W[index_dir] < global_current_time_step ) {
    lattice[index_site].on_step_W[index_dir] = global_current_time_step;
    su3mat_copy( &(lattice[index_site].Wlink[index_dir]),
                 &(lattice[index_site].Wlink_previous[index_dir]) );
    su3mat_copy( b, &(lattice[index_site].Wlink[index_dir]) );
  }
#endif /* HISQ_REUNITARIZATION_DEBUG */
}


/* Derivative of the unitarized matrix with respect to the original:
   dW/dU and d(W^+)/dU (at fixed U^+ !), where W=U(U^+U)^-1/2 */
void su3_unit_der( su3_matrix *u, su3_tensor4 *dwdu, su3_tensor4 *dwdagdu ) {
  su3_matrix up, um, a, b, up12, um12, upu, umu, dw, dwdag;
  int i, j, m, n;
  Real factor;

  factor = 0.5/SU3_UNIT_DER_EPS;

  /* loop on components of the original matrix U */
  for( i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      /* temporary storage */
      su3mat_copy( u, &up );
      su3mat_copy( u, &um );

      /* add/subtract epsilon to the real part of [i][j] component */
      (up.e[i][j]).real += SU3_UNIT_DER_EPS;
      (um.e[i][j]).real -= SU3_UNIT_DER_EPS;

      /* Up12=(U^+ Up)^-1/2 */
      mult_su3_an( u, &up, &a );
      su3_root_inv( &a, &b, &up12 );

      /* Um12=(U^+ Um)^-1/2 */
      mult_su3_an( u, &um, &a );
      su3_root_inv( &a, &b, &um12 );

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

  factor = 0.5/SU3_UNIT_DER_EPS;

  /* loop on components of the original matrix U */
  for( i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      /* temporary storage */
      su3mat_copy( u, &up );
      su3mat_copy( u, &um );

      /* add/subtract epsilon to the real part of [i][j] component */
      (up.e[i][j]).real += SU3_UNIT_DER_EPS;
      (um.e[i][j]).real -= SU3_UNIT_DER_EPS;

      /* Up12=(U^+ Up)^-1/2 */
      mult_su3_an( u, &up, &a );
      su3_root_inv( &a, &b, &up12 );

      /* Um12=(U^+ Um)^-1/2 */
      mult_su3_an( u, &um, &a );
      su3_root_inv( &a, &b, &um12 );

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

  factor = 0.5/SU3_UNIT_DER_EPS;

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
      (upre.e[i][j]).real += SU3_UNIT_DER_EPS;
      (umre.e[i][j]).real -= SU3_UNIT_DER_EPS;
      (upim.e[i][j]).imag += SU3_UNIT_DER_EPS;
      (umim.e[i][j]).imag -= SU3_UNIT_DER_EPS;

      /* Up12=(Up^+ Up)^-1/2 */
      mult_su3_an( &upre, &upre, &a );
      su3_root_inv( &a, &b, &up12re );
      mult_su3_an( &upim, &upim, &a );
      su3_root_inv( &a, &b, &up12im );

      /* Um12=(Um^+ Um)^-1/2 */
      mult_su3_an( &umre, &umre, &a );
      su3_root_inv( &a, &b, &um12re );
      mult_su3_an( &umim, &umim, &a );
      su3_root_inv( &a, &b, &um12im );

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
void su3_unitarize_rational( su3_matrix *V, su3_matrix *W ) {
  int l;
  su3_matrix H, T, X, Z;

  /* hermitian matrix: H=V^+V */
  mult_su3_an( V, V, &H );

  clear_su3mat( &X );

  for( l=1; l<=SU3_UNIT_RAT_NTERMS; l++) {
    su3mat_copy( &H, &T );

    /* add d_l to H */
    T.e[0][0].real += d_l_SU3_UNIT_RAT[l];
    T.e[1][1].real += d_l_SU3_UNIT_RAT[l];
    T.e[2][2].real += d_l_SU3_UNIT_RAT[l];

    /* calculate inverse of H+d_l */
    su3_inverse( &T, &Z );

    /* add c_l/(H+d_l) */
    scalar_mult_add_su3_matrix( &X, &Z, c_l_SU3_UNIT_RAT[l], &T );

    su3mat_copy( &T, &X );
  }

  /* add c_0 */
  X.e[0][0].real += c_l_SU3_UNIT_RAT[0];
  X.e[1][1].real += c_l_SU3_UNIT_RAT[0];
  X.e[2][2].real += c_l_SU3_UNIT_RAT[0];

  mult_su3_nn( V, &X, W );
}



/* Derivative of a unitarized matrix with rational approximation */
void su3_unit_der_rational( su3_matrix *V, su3_tensor4 *dwdv, su3_tensor4 *dwdagdv ) {
  int i, j, l, p, q, r, s;
  su3_matrix Kl[ SU3_UNIT_RAT_NTERMS ]; // store intermediate inverted matrices
  su3_matrix Vdag, H, T, X;
  su3_tensor4 A4, B4;
  complex ftmp, ftmp2, ftmp3, ftmp4;

  /* adjoint, needed later */
  su3_adjoint( V, &Vdag );

  /* hermitian matrix: H=V^+V */
  mult_su3_an( V, V, &H );

  for( l=1; l<=SU3_UNIT_RAT_NTERMS; l++) {
    su3mat_copy( &H, &T );

    /* add d_l to H */
    T.e[0][0].real += d_l_SU3_UNIT_RAT[l];
    T.e[1][1].real += d_l_SU3_UNIT_RAT[l];
    T.e[2][2].real += d_l_SU3_UNIT_RAT[l];

    /* calculate inverse of H+d_l and store in Kl array */
    su3_inverse( &T, &( Kl[l-1] ) );
  }

  /* zero out tensor4 */
  for( p=0; p<3; p++) {
    for( r=0; r<3; r++) {
      for( s=0; s<3; s++) {
        for( q=0; q<3; q++) {
          (A4.t4[p][r][s][q]).real = 0.0;
          (A4.t4[p][r][s][q]).imag = 0.0;
          (B4.t4[p][r][s][q]).real = 0.0;
          (B4.t4[p][r][s][q]).imag = 0.0;
        }
      }
    }
  }

  clear_su3mat( &X );

  for( l=1; l<=SU3_UNIT_RAT_NTERMS; l++) {
    /* assemble tensor4 out of Kls */
    for( p=0; p<3; p++) {
      for( r=0; r<3; r++) {
        for( s=0; s<3; s++) {
          for( q=0; q<3; q++) {
            CMUL( Kl[l-1].e[p][r], Kl[l-1].e[s][q], ftmp );
            (B4.t4[p][r][s][q]).real += c_l_SU3_UNIT_RAT[l] * ftmp.real;
            (B4.t4[p][r][s][q]).imag += c_l_SU3_UNIT_RAT[l] * ftmp.imag;
          }
        }
      }
    }
    scalar_mult_add_su3_matrix( &X, &Kl[ l-1 ], c_l_SU3_UNIT_RAT[ l ], &T );
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
              (dwdv->t4[p][r][s][q]).real += c_l_SU3_UNIT_RAT[0];
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
void su3_unitarize_analytic( su3_matrix *V, su3_matrix *W ) {
  su3_matrix Q, Q2, Q3, S1, S2;
  Real c0, c1, c2, S, S3, R, R2, CQ3, RoS, theta, theta3, pi23, denom;
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
  if( fabs(S)<SU3_UNIT_ANALYTIC_EPS ) {
    /* eigenvalues of Q */
    g0 = c0/3;
    g1 = c0/3;
    g2 = c0/3;
  }
  else {
    R = c2/2 - c0 * (c1/3) + c0 * c0 * (c0/27);
    R2 = R*R;
    CQ3 = S*S*S;
    S = sqrt(S);
    S3 = S*S*S;
    /* treat possible underflow: R/S^3/2>1.0 leads to acos giving NaN */
    RoS = R/S3;
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
        printf("Hit NaN in su3_unitarize_analytic()\n");
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

  if( ws < SU3_UNIT_ANALYTIC_EPS ) {
    printf( "WARNING: su3_unitarize_analytic: ws is too small!\n" );
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

}


/* Analytic unitarization, Hasenfratz, Hoffmann, Schaefer, JHEP05 (2007) 029 */
void su3_unitarize_analytic_index( su3_matrix *V, su3_matrix *W, int index_site, int index_dir ) {
  su3_matrix Q, Q2, Q3, S1, S2;
  Real c0, c1, c2, S, S3, R, R2, CQ3, RoS, theta, theta3, pi23, denom;
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
  if( fabs(S)<SU3_UNIT_ANALYTIC_EPS ) {
    /* eigenvalues of Q */
    g0 = c0/3;
    g1 = c0/3;
    g2 = c0/3;
  }
  else {
    R = c2/2 - c0 * (c1/3) + c0 * c0 * (c0/27);
    R2 = R*R;
    CQ3 = S*S*S;
    S = sqrt(S);
    S3 = S*S*S;
    /* treat possible underflow: R/S^3/2>1.0 leads to acos giving NaN */
    RoS = R/S3;
#ifdef HISQ_REUNITARIZATION_DEBUG
   /* DEBUG: set internal variable */
    lattice[index_site].RoS[index_dir]=RoS;
#endif /* HISQ_REUNITARIZATION_DEBUG */
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
        printf("Hit NaN in su3_unitarize_analytic()\n");
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

  if( ws < SU3_UNIT_ANALYTIC_EPS ) {
    printf( "WARNING: su3_unitarize_analytic: ws is too small!\n" );
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
#endif /* HISQ_REUNITARIZATION_DEBUG */
}


/* Analytic derivative of the unitarized matrix with respect to the original:
   dW/dV and d(W^+)/dV, where W=V(V^+V)^-1/2 */
void su3_unit_der_analytic( su3_matrix *V, 
               su3_tensor4 *dwdv, su3_tensor4 *dwdagdv ) {
  su3_matrix Q, Q2, Q3, S1, S2, W, Q12;
  Real c0, c1, c2, S, S3, RoS, R, R2, CQ3, theta, theta3, pi23, denom;
  Real g0, g1, g2, g0sq, g1sq, g2sq, f0, f1, f2, us, vs, ws;
  Real u2, u3, u4, u5, u6, u7, u8, v2, v3, v4, v5, v6, w2, w3, w4, w5;
  Real b00, b01, b02, b11, b12, b22, denom3;
  complex der, ctmp, ctmp2;
  su3_matrix VVd, VQ, QVd, QQVd, VQQ, VQVd, PVd, RVd, SVd, Vd;
  int i, j, m, n;


  /* adjoint */
  su3_adjoint( V, &Vd );

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
  if( fabs(S)<SU3_UNIT_ANALYTIC_EPS ) {
    /* eigenvalues of Q */
    g0 = c0/3;
    g1 = c0/3;
    g2 = c0/3;
  }
  else {
    R = c2/2 - c0 * (c1/3) + c0 * c0 * (c0/27);
    R2 = R*R;
    CQ3 = S*S*S;
    S = sqrt(S);
    S3 = S*S*S;
    /* treat possible underflow: R/S^3/2>1.0 leads to acos giving NaN */
    RoS = R/S3;
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
//      if(isnan(theta)){printf("Hit NaN in su3_unit_der_analytic()\n"); terminate(0);}
      if(isnan(theta)){
        printf("Hit NaN in su3_unit_der_analytic()\n");
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

  if( ws < SU3_UNIT_ANALYTIC_EPS ) {
    printf( "WARNING: su3_unit_der_analytic: ws is too small!\n" );
  }

  denom = ws * ( us*vs - ws );

  /* constants in inverse root expression */
  f0 = ( us*vs*vs - ws*(us*us+vs) ) / denom;
  f1 = ( 2*us*vs - ws - us*us*us ) / denom;
  f2 = us / denom;

  /* assemble inverse root: Q^-1/2 = f0 + f1*Q + f2*Q^2 */
  scalar_mult_su3_matrix( &Q2, f2, &S1 );
  scalar_mult_add_su3_matrix( &S1, &Q, f1, &Q12 );
  Q12.e[0][0].real += f0;
  Q12.e[1][1].real += f0;
  Q12.e[2][2].real += f0;

  /* W = V*Q12 */
  mult_su3_nn( V, &Q12, &W );

  denom3 = 2*denom*denom*denom;

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
  w5 = w4*ws;
  b00 = -w3*u6 +3*u4*( vs*w3 +v4*ws )
        -u3*( v6 +4*w4 +12*v3*w2 ) + u2*( 16*v2*w3 +3*v5*ws )
        -us*( 8*vs*w4 +3*v4*w2) +w5 +v3*w3;
  b00 /= denom3;
  b01 = -w2*u7 -v2*ws*u6 +u5*( v4 +6*vs*w2 ) -u4*( 5*w3 +v3*ws )
        -u3*( 2*v5 +6*v2*w2 ) +u2*( 10*vs*w3 +6*v4*ws )
        -us*( 3*w4 +6*v3*w2 ) +2*v2*w3;
  b01 /= denom3;
  b02 = w2*u5 +v2*ws*u4 -u3*( v4 +4*vs*w2 )
        +u2*( 4*w3 +3*v3*ws ) -3*v2*w2*us +vs*w3;
  b02 /= denom3;
  b11 = -ws*u8 -v2*u7 +7*vs*ws*u6 +u5*( 4*v3 -5*w2 ) -16*v2*ws*u4
        -u3*( 4*v4 -16*vs*w2 ) -u2*( 3*w3 -12*v3*ws ) -12*v2*w2*us +3*vs*w3;
  b11 /= denom3;
  b12 = ws*u6 +v2*u5 -5*vs*ws*u4 -u3*( 2*v3 -4*w2 ) +6*v2*ws*u2 -6*vs*w2*us+w3;
  b12 /= denom3;
  b22 = -ws*u4 -v2*u3 +3*vs*ws*u2 -3*w2*us;
  b22 /= denom3;

  /* ** create several building blocks for derivative ** */
  mult_su3_nn( V, &Q, &VQ );
  mult_su3_na( &Q, V, &QVd );
  mult_su3_na( V, V, &VVd );
  mult_su3_nn( V, &Q2, &VQQ );
  mult_su3_na( &Q2, V, &QQVd );
  mult_su3_na( &VQ, V, &VQVd );
  scalar_mult_su3_matrix( &QVd, b01, &S1 );
  scalar_mult_add_su3_matrix( &S1, &QQVd, b02, &S2 );
  scalar_mult_add_su3_matrix( &S2, &Vd, b00, &PVd );
  scalar_mult_su3_matrix( &QVd, b11, &S1 );
  scalar_mult_add_su3_matrix( &S1, &QQVd, b12, &S2 );
  scalar_mult_add_su3_matrix( &S2, &Vd, b01, &RVd );
  scalar_mult_su3_matrix( &QVd, b12, &S1 );
  scalar_mult_add_su3_matrix( &S1, &QQVd, b22, &S2 );
  scalar_mult_add_su3_matrix( &S2, &Vd, b02, &SVd );

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
          }
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
        }
      }
    }
  }
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
