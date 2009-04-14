/******************  eigen_su3_UdU.c  (in su3.a) ************************
*									*
* void eigen_su3_UdU( su3_matrix *U, Real *g0, Real *g1, Real *g2)	*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include <math.h>

/* Return: g0, g1, g2 are eigenvalues of U^+ U matrix */
void eigen_su3_UdU( su3_matrix *U, Real *g0, Real *g1, Real *g2) {
  su3_matrix H,Q,Q2,Q3;
  Real c0,c1,c2,c03,c2abs,c2max,eps,S,theta,theta3,pi23;

  /* construct Hermitian matrix H=U^+ U */
  mult_su3_an( U, U, &H );

  /* make traceless Hermitian Q=H-1/3TrH */
  c0=(trace_su3( &H )).real;
  c03=c0/3;
  su3mat_copy( &H, &Q );
  Q.e[0][0].real -= c03;
  Q.e[1][1].real -= c03;
  Q.e[2][2].real -= c03;

  /* Q^2 and Q^3 and traces */
  mult_su3_nn( &Q, &Q, &Q2 );
  mult_su3_nn( &Q2, &Q, &Q3 );
  c1 = ( Q2.e[0][0].real + Q2.e[1][1].real + Q2.e[2][2].real ) / 2;
  c2 = ( Q3.e[0][0].real + Q3.e[1][1].real + Q3.e[2][2].real ) / 3;

  c2abs = fabs( c2 );
  c2max = 2*pow( c1/3., 1.5);

  eps = (c2max - c2abs)/c2max;

  if( eps < 0 ) { // only due to rounding errors, means theta=0
    theta=0;
  }
  else {
    if ( eps < 1.0e-3 ) { // c2 is close to the maximum, use series expansion
      theta = sqrt(2*eps)*( 1 + ( (1./12.) + ( (3./160.) + ( (5./896.) 
            + ( (35./18432.) + (63./90112.)*eps ) *eps) *eps) *eps) *eps);
    }
    else { // normal case
      theta = acos( c2abs/c2max );
    }
  }

  S = sqrt( c1/3. );
  theta3 = theta/3;
  pi23 = 6.2831853071795862319959269370883703 / 3;
  *g0 = c03 + 2*S*cos( theta3 );
  *g1 = c03 + 2*S*cos( theta3 + pi23 );
  *g2 = c03 + 2*S*cos( theta3 + 2*pi23 );

}

