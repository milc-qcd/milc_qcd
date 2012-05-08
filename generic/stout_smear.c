/*********************** stout_smear.c **********************************/
/* MIMD version 7 */

/* Adapted from the Chroma code */
/* Original idea M. Peardon and C. Morningstar Phys.Rev.D69:054501 (2004) */
/* Algorithmic improvements to handle 0/0 cases by Robert Edwards */

#include "generic_includes.h"

/*--------------------------------------------------------------------*/
/* Construct Q from smeared link V and unsmeared link U */
/* Quick and dirty code - can be optimized for SU(3) */

static void 
get_Q_from_VUadj(su3_matrix *Q, su3_matrix *V, su3_matrix *U){

  complex x;
  complex tr;
  su3_matrix Om;

  x = cmplx(0, 0.5);  /* i/2 */

  /* Om = V U^adj */
  mult_su3_na( V, U, &Om );

  /* Q = i/2(Om^adj - Om) */
  su3_adjoint( &Om, Q );
  sub_su3_matrix( Q, &Om, &Om );
  c_scalar_mult_su3mat( &Om, &x, Q );
  /* Q = Q - Tr Q/3 */
  tr = trace_su3( Q );
  CDIVREAL(tr, 3., tr);
  CSUB(Q->e[0][0],tr,Q->e[0][0]);
  CSUB(Q->e[1][1],tr,Q->e[1][1]);
  CSUB(Q->e[2][2],tr,Q->e[2][2]);
}

/*--------------------------------------------------------------------*/
/* Get the coefficients f of the expansion of exp(iQ)                   */
static void 
get_fs_from_Qs( complex f[3], su3_matrix *Q, su3_matrix *QQ){

  su3_matrix QQQ;
  Real trQQQ, trQQ, c0, c1;
  Real c0abs, c0max, theta;
  Real eps, sqtwo = sqrt(2.);
  Real u, w, u_sq, w_sq, xi0;
  Real cosu, sinu, cosw, sin2u, cos2u, ucosu, usinu, ucos2u, usin2u;
  Real denom, subexp1, subexp2, subexp3, subexp;

  mult_su3_nn ( Q, QQ, &QQQ );

  trQQ  = (trace_su3( QQ )).real;
  trQQQ = (trace_su3( &QQQ )).real;

  c0 = 1./3. * trQQQ;
  c1 = 1./2. * trQQ;

  if( c1 < 4.0e-3  )
    { // RGE: set to 4.0e-3 (CM uses this value). I ran into nans with 1.0e-4
      // =======================================================================
      //
      // Corner Case 1: if c1 < 1.0e-4 this implies c0max ~ 3x10^-7
      //    and in this case the division c0/c0max in arccos c0/c0max can be undefined
      //    and produce NaN's

      // In this case what we can do is get the f-s a different way. We go back to basics:
      //
      // We solve (using maple) the matrix equations using the eigenvalues
      //
      //  [ 1, q_1, q_1^2 ] [ f_0 ]       [ exp( iq_1 ) ]
      //  [ 1, q_2, q_2^2 ] [ f_1 ]   =   [ exp( iq_2 ) ]
      //  [ 1, q_3, q_3^2 ] [ f_2 ]       [ exp( iq_3 ) ]
      //
      // with q_1 = 2 u w, q_2 = -u + w, q_3 = - u - w
      //
      // with u and w defined as  u = sqrt( c_1/ 3 ) cos (theta/3)
      //                     and  w = sqrt( c_1 ) sin (theta/3)
      //                          theta = arccos ( c0 / c0max )
      // leaving c0max as a symbol.
      //
      //  we then expand the resulting f_i as a series around c0 = 0 and c1 = 0
      //  and then substitute in c0max = 2 ( c_1/ 3)^(3/2)
      //
      //  we then convert the results to polynomials and take the real and imaginary parts:
      //  we get at the end of the day (to low order)

      //                  1    2
      //   f0[re] := 1 - --- c0  + h.o.t
      //                 720     
      //
      //               1       1           1        2
      //   f0[im] := - - c0 + --- c0 c1 - ---- c0 c1   + h.o.t
      //               6      120         5040       
      //
      //
      //             1        1            1        2
      //   f1[re] := -- c0 - --- c0 c1 + ----- c0 c1  +  h.o.t
      //             24      360         13440        f
      //
      //                 1       1    2    1     3    1     2
      //   f1[im] := 1 - - c1 + --- c1  - ---- c1  - ---- c0   + h.o.t
      //                 6      120       5040       5040
      //
      //               1   1        1    2     1     3     1     2
      //   f2[re] := - - + -- c1 - --- c1  + ----- c1  + ----- c0  + h.o.t
      //               2   24      720       40320       40320   
      //
      //              1        1              1        2
      //   f2[im] := --- c0 - ---- c0 c1 + ------ c0 c1  + h.o.t
      //             120      2520         120960

      //  We then express these using Horner's rule for more stable evaluation.
      //
      //  =====================================================================

      f[0].real = 1. - c0*c0/720.;
      f[0].imag = -(c0/6.)*(1. - (c1/20.)*(1. - (c1/42.))) ;

      f[1].real =  c0/24.*(1. - c1/15.*(1. - 3.*c1/112.)) ;
      f[1].imag =  1.-c1/6.*(1. - c1/20.*(1. - c1/42.)) - c0*c0/5040. ;

      f[2].real = 0.5*(-1. + c1/12.*(1. - c1/30.*(1. - c1/56.)) + c0*c0/20160.);
      f[2].imag = 0.5*(c0/60.*(1. - c1/21.*(1. - c1/48.)));

    }
  else
    {
      // =======================================================================
      // Normal case: Do as per Morningstar-Peardon paper
      // =======================================================================

      c0abs = fabs( c0 );
      c0max = 2*pow( c1/3., 1.5);

      // =======================================================================
      // Now work out theta. In the paper the case where c0 -> c0max even when c1 is reasonable
      // Has never been considered, even though it can arise and can cause the arccos function
      // to fail
      // Here we handle it with series expansion
      // =======================================================================
      eps = (c0max - c0abs)/c0max;

      if( eps < 0 ) {
        // =====================================================================
        // Corner Case 2: Handle case when c0abs is bigger than c0max.
        // This can happen only when there is a rounding error in the ratio, and that the
        // ratio is really 1. This implies theta = 0 which we'll just set.
        // =====================================================================
        theta = 0;
      }
      else if ( eps < 1.0e-3 ) {
        // =====================================================================
        // Corner Case 3: c0->c0max even though c1 may be actually quite reasonable.
        // The ratio |c0|/c0max -> 1 but is still less than one, so that a
        // series expansion is possible.
        // SERIES of acos(1-epsilon): Good to O(eps^6) or with this cutoff to O(10^{-18}) Computed with Maple.
        //  BTW: 1-epsilon = 1 - (c0max-c0abs)/c0max = 1-(1 - c0abs/c0max) = +c0abs/c0max
        //
        // ======================================================================
        theta = sqtwo*sqrt(eps)*( 1 + ( (1./12.) + ( (3./160.) + ( (5./896.) + ( (35./18432.) + (63./90112.)*eps ) *eps) *eps) *eps) *eps);
      }
      else {
        theta = acos( c0abs/c0max );
      }

      u = sqrt(c1/3.)*cos(theta/3.);
      w = sqrt(c1)*sin(theta/3.);

      u_sq = u*u;
      w_sq = w*w;

      if( fabs(w) < 0.05 ) {
        xi0 = 1. - (1./6.)*w_sq*( 1. - (1./20.)*w_sq*( 1. - (1./42.)*w_sq ) );
      }
      else {
        xi0 = sin(w)/w;
      }

      cosu = cos(u);
      sinu = sin(u);
      cosw = cos(w);
      sin2u = sin(2*u);
      cos2u = cos(2*u);
      ucosu = u*cosu;
      usinu = u*sinu;
      ucos2u = u*cos2u;
      usin2u = u*sin2u;

      denom = 9.*u_sq - w_sq;

      subexp1 = u_sq - w_sq;
      subexp2 = 8*u_sq*cosw;
      subexp3 = (3*u_sq + w_sq)*xi0;

      f[0].real = ( (subexp1)*cos2u + cosu*subexp2 + 2*usinu*subexp3 ) / denom ;
      f[0].imag = ( (subexp1)*sin2u - sinu*subexp2 + 2*ucosu*subexp3 ) / denom ;

      subexp = (3*u_sq -w_sq)*xi0;

      f[1].real = (2*(ucos2u - ucosu*cosw)+subexp*sinu)/denom;
      f[1].imag = (2*(usin2u + usinu*cosw)+subexp*cosu)/denom;

      subexp=3*xi0;

      f[2].real = (cos2u - cosu*cosw -usinu*subexp) /denom ;
      f[2].imag = (sin2u + sinu*cosw -ucosu*subexp) /denom ;

      if( c0 < 0 ) {

        // f[0] = conj(f[0]);
        f[0].imag *= -1;
       
        //f[1] = -conj(f[1]);
        f[1].real *= -1;
       
        //f[2] = conj(f[2]);
        f[2].imag *= -1;

      }
    } // End of if( corner_caseP ) else {}
}

void exp_iQ( su3_matrix *T, su3_matrix *Q )
{

  complex f[3];
  su3_matrix QQ;

  mult_su3_nn( Q, Q, &QQ );
  get_fs_from_Qs( f, Q, &QQ );

  /*   f[0] + f[1]*Q + f[2]*QQ */

  clear_su3mat( T );
  T->e[0][0] = f[0];
  T->e[1][1] = f[0];
  T->e[2][2] = f[0];

  c_scalar_mult_add_su3mat( T, Q, &f[1], T );
  c_scalar_mult_add_su3mat( T, &QQ, &f[2], T );
}

/* Do Morninstar-Peardon stout smearing to construct unitary W from
   smeared link V and unsmeared link U */

void stout_smear(su3_matrix *W, su3_matrix *V, su3_matrix *U)
{

  su3_matrix Q, T;

  get_Q_from_VUadj( &Q, V, U);

  /* T = exp(iQ) */
  exp_iQ( &T, &Q );

  /* W = exp(iQ) U */
  mult_su3_nn( &T, U, W);
}
