/*
    Test matrix and tensor4 operations from su3_mat_op.c
    Good for finding reunitarization, derivative issues, etc.
    AB, 1/25/8 - started
*/

#define CONTROL

#include "ks_imp_includes.h"    /* definitions files and prototypes */
#include "../generic_ks/su3_mat_op.h"


#define MAX_AMPLITUDE 10 /* amplitude of the random matrix */
#define MAX_REPEAT 10000 /* number of repetitions */



int main(void) {
  su3_matrix Unit, A, B, B2, B3, C;
  int irpt, i, j;
  Real norm;
  complex det;

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
    /* unitarize: aproximate reunitarization using recursion for ^-1/2 */
    su3_unitarize( &A, &B );
    printf( "Approximate reunitarization: B=A*(A^+A)^-1/2:\n" );
    dumpmat( &B );
    det = det_su3( &B );
    printf( "det(B): (%18.10g,%18.10g)\n\n", det.real, det.imag );

    /* unitarize_rational: approximate reunitarization 
         using raional function approximation for ^-1/2 */
    su3_unitarize_rational( &A, &B2 );
    printf( "Approximate reunitarization with rational function: B2=A*(A^+A)^-1/2:\n" );
    dumpmat( &B2 );
    det = det_su3( &B2 );
    printf( "det(B2): (%18.10g,%18.10g)\n\n", det.real, det.imag );

    /* unitarize_analytic: exact reunitarization */
    su3_unitarize_analytic( &A, &B3 );
    printf( "Analytic reunitarization: B3=A*(A^+A)^-1/2:\n" );
    dumpmat( &B3 );
    det = det_su3( &B3 );
    printf( "det(B3): (%18.10g,%18.10g)\n\n", det.real, det.imag );

#endif /* TEST_SU3_UNITARIZATION */

  } /* end of MAIN REPETITION LOOP */

  return 0;
}

