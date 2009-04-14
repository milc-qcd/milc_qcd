/*
    Test matrix and tensor4 operations from su3_mat_op.c
    Good for finding reunitarization, derivative issues, etc.
    AB, 1/25/8 - started
*/

#define CONTROL

#include "ks_imp_includes.h"    /* definitions files and prototypes */
#include "../include/su3_mat_op.h"
#include <quark_action.h>


#define MAX_AMPLITUDE 10 /* amplitude of the random matrix */
#define MAX_REPEAT 10000 /* number of repetitions */


int main(void) {
  su3_matrix Unit, A, B, C, D;
  int irpt, i, j;
  Real norm, g0, g1, g2;
  complex det;
#ifdef TEST_SU3_UNITARIZATION_DERIVATIVE
  su3_tensor4 dydv, dydagdv;
  int k, l;
#ifdef HISQ_FORCE_FILTER
  printf("HISQ_FORCE_FILTER=%f\n",HISQ_FORCE_FILTER);
#endif
#endif /* TEST_SU3_UNITARIZATION_DERIVATIVE */

  /* initialize random number generator */
  srand( 0 );

  /* unit matrix */
  unity_su3mat( &Unit );

#if ( PRECISION==1 )
  printf( "Testing in SINGLE precision\n" );
#else
  printf( "Testing in DOUBLE precision\n" );
#endif /* PRECISION */


  /* MAIN REPETITION LOOP */
  for( irpt=0; irpt<MAX_REPEAT; irpt++ ) {
    /* randomize matrix */
    for( i=0; i<3; i++) {
      for( j=0; j<3; j++) {
        A.e[i][j].real
            = MAX_AMPLITUDE*( 2 * ( ((double)rand())/RAND_MAX ) - 1 );
        A.e[i][j].imag 
            = MAX_AMPLITUDE*( 2 * ( ((double)rand())/RAND_MAX ) - 1 );
      }
    }
    printf( "*** Iteration %d\n", irpt );
    printf( "Random 3x3 complex matrix A:\n" );
    dumpmat( &A );
    det=det_su3( &A );
    printf( "det(A): (%18.10g,%18.10g)\n\n", det.real, det.imag );

#ifdef TEST_SU3_INVERSION
    su3_inverse( &A, &B );
    printf( "Inverted matrix B=A^-1:\n" );
    dumpmat( &B );

    /* calculate (A*A^-1 - I) */
    mult_su3_nn( &A, &B, &C );
    sub_su3_matrix( &C, &Unit, &B );

    /* calculate Frobenius norm */
    norm = su3_norm_frob( &B );
    printf( "Frobenius norm of (A*B-I): %18.10g\n\n", norm );
#endif /* TEST_SU3_INVERSION */


#ifdef TEST_SU3_UNITARIZATION
#if ( UNITARIZATION_METHOD==UNITARIZE_ROOT )
    printf("Unitarization method = UNITARIZE_ROOT\n");
    su3_unitarize( &A, &B );
    printf( "Analytic reunitarization: B=A*(A^+A)^-1/2:\n" );
    dumpmat( &B );
    det = det_su3( &B );
    printf( "det(B): (%18.10g,%18.10g)\n\n", det.real, det.imag );
    mult_su3_an( &B, &B, &C );
    sub_su3_matrix( &C, &Unit, &D );
    norm = su3_norm_frob( &D );
    printf( "Deviation from unity: %18.10g\n\n", norm );
    eigen_su3_UdU( &A, &g0, &g1, &g2);
    printf( "Eigenvalues of A^+A: %18.10g%18.10g%18.10g\n\n", g0, g1, g2 );

#elif ( UNITARIZATION_METHOD==UNITARIZE_RATIONAL )
    printf("Unitarization method = UNITARIZE_RATIONAL\n");
    su3_unitarize_rational( &A, &B );
    printf( "Analytic reunitarization: B=A*(A^+A)^-1/2:\n" );
    dumpmat( &B );
    det = det_su3( &B );
    printf( "det(B): (%18.10g,%18.10g)\n\n", det.real, det.imag );
    mult_su3_an( &B, &B, &C );
    sub_su3_matrix( &C, &Unit, &D );
    norm = su3_norm_frob( &D );
    printf( "Deviation from unity: %18.10g\n\n", norm );
    eigen_su3_UdU( &A, &g0, &g1, &g2);
    printf( "Eigenvalues of A^+A: %18.10g%18.10g%18.10g\n\n", g0, g1, g2 );

#elif ( UNITARIZATION_METHOD==UNITARIZE_ANALYTIC )
    printf("Unitarization method = UNITARIZE_ANALYTIC\n");
    su3_unitarize_analytic( &A, &B );
    printf( "Analytic reunitarization: B=A*(A^+A)^-1/2:\n" );
    dumpmat( &B );
    det = det_su3( &B );
    printf( "det(B): (%18.10g,%18.10g)\n\n", det.real, det.imag );
    mult_su3_an( &B, &B, &C );
    sub_su3_matrix( &C, &Unit, &D );
    norm = su3_norm_frob( &D );
    printf( "Deviation from unity: %18.10g\n\n", norm );
    eigen_su3_UdU( &A, &g0, &g1, &g2);
    printf( "Eigenvalues of A^+A: %18.10g%18.10g%18.10g\n\n", g0, g1, g2 );
#else
    printf("Unknown unitarization method\n"); terminate(0);
#endif

#endif /* TEST_SU3_UNITARIZATION */


#ifdef TEST_SU3_UNITARIZATION_DERIVATIVE
#if ( UNITARIZATION_METHOD==UNITARIZE_ROOT )
    printf("Unitarization method = UNITARIZE_ROOT\n");
    su3_unit_der( &A, &dydv, &dydagdv );
#elif ( UNITARIZATION_METHOD==UNITARIZE_RATIONAL )
    printf("Unitarization method = UNITARIZE_RATIONAL\n");
    su3_unit_der_rational( &A, &dydv, &dydagdv );
#elif ( UNITARIZATION_METHOD==UNITARIZE_ANALYTIC )
    printf("Unitarization method = UNITARIZE_ANALYTIC\n");
    su3_unit_der_analytic( &A, &dydv, &dydagdv );
#else
    printf("Unknown unitarization method\n"); terminate(0);
#endif
    // print t4's
    printf( "dW/dA (4-rank tensor):\n" );
    for(i=0;i<3;i++) {
      for(j=0;j<3;j++) {
        for(k=0;k<3;k++) {
          for(l=0;l<3;l++) {
            printf("%d%d%d%dre=%24.18f  %d%d%d%dim=%24.18f\n",
                   i,j,k,l,dydv.t4[i][j][k][l].real,
                   i,j,k,l,dydv.t4[i][j][k][l].imag);
          }
        }
      }
    }
    norm = su3_t4_norm_frob( &dydv );
    printf("dW/dA norm: %28.18f\n\n",norm);

    printf( "dW^+/dA (4-rank tensor):\n" );
    for(i=0;i<3;i++) {
      for(j=0;j<3;j++) {
        for(k=0;k<3;k++) {
          for(l=0;l<3;l++) {
            printf("%d%d%d%dre=%24.18f  %d%d%d%dim=%24.18f\n",
                   i,j,k,l,dydagdv.t4[i][j][k][l].real,
                   i,j,k,l,dydagdv.t4[i][j][k][l].imag);
          }
        }
      }
    }
    norm = su3_t4_norm_frob( &dydagdv );
    printf("dW^+/dA norm: %28.18f\n\n",norm);
#endif /* TEST_SU3_UNITARIZATION_DERIVATIVE */


  } /* end of MAIN REPETITION LOOP */

  return 0;
}

