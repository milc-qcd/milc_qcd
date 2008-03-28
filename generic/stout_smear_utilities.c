/*********************** stout_smear_utilities.c ***************************/
/* MIMD version 7 */

/* Adapted for MILC from the Chroma code */
/* Original idea M. Peardon and C. Morningstar Phys.Rev.D69:054501 (2004) */
/* Algorithmic improvements to handle 0/0 cases by Robert Edwards */
/* 2/2008 CD - corrected some coding types */

#include "stout_smear.h"

/*--------------------------------------------------------------------*/
/* Construct Q from smeared link V and unsmeared link U */
/* Quick and dirty code - can be optimized for SU(3) */

/* a = traceless-hermitian part of b.  b and a may be equivalent. */
static void
traceless_hermitian_su3(su3_matrix *a, su3_matrix *b)
{
  complex t;
  su3_matrix c;
  
  su3_adjoint( b, &c );
  add_su3_matrix( &c, b, a );
  
  t = trace_su3( a );
  
  CDIVREAL(t, 3., t);
  CSUB(a->e[0][0],t,a->e[0][0]);
  CSUB(a->e[1][1],t,a->e[1][1]);
  CSUB(a->e[2][2],t,a->e[2][2]);
  
  scalar_mult_su3_matrix( a, 0.5, a );
}

/*--------------------------------------------------------------------*/

/*
* return Tr( A*B )   						*
*/

static complex 
complextrace_su3_nn( su3_matrix *a, su3_matrix *b ) {
  register int i,j;
  register Real sumr, sumi;
  complex sum;
  for(sumr=0.0,sumi=0.0,i=0;i<3;i++)for(j=0;j<3;j++){
      sumr+= a->e[i][j].real*b->e[j][i].real - a->e[i][j].imag*b->e[j][i].imag;
      sumi+= a->e[i][j].real*b->e[j][i].imag + a->e[i][j].imag*b->e[j][i].real;
    }
  sum.real= sumr; sum.imag=sumi; 
  return sum;
}

/*--------------------------------------------------------------------*/

static void 
get_Q_from_VUadj(su3_matrix *Q, su3_matrix *V, su3_matrix *U){

  complex minusI;
  complex tr;
  su3_matrix Omega;

  minusI = cmplx(0, -1);  /* -i */

  /* Omega = V U^adj */
  mult_su3_na( V, U, &Omega );

  /* Q = traceless hermitian part of [-i Omega] */
  c_scalar_mult_su3mat( &Omega, &minusI, &Omega );
  traceless_hermitian_su3( Q, &Omega );
}

/*--------------------------------------------------------------------*/
/* Get the coefficients f of the expansion of exp(iQ)                
   and if do_bs = true
   get the coefficients b1 and b2 needed for the expansion of d exp(iQ)/dt */

static void 
get_fs_and_bs_from_Qs( complex f[3], complex b1[3], complex b2[3], 
		       su3_matrix *Q, su3_matrix *QQ, int do_bs )
{

  su3_matrix QQQ;
  Real trQQQ, trQQ, c0, c1;
  Real c0abs, c0max, theta;
  Real eps, sqtwo = sqrt(2.);
  Real u, w, u_sq, w_sq, xi0, xi1;
  Real cosu, sinu, cosw, sinw, sin2u, cos2u, ucosu, usinu, ucos2u, usin2u;
  Real denom;

  Real r_1_re[3], r_1_im[3], r_2_re[3], r_2_im[3];
  Real b_denom;

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
      //              1        1            1        2 
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
      //  to get the b-s we use the fact that
      //                                      b2_i = d f_i / d c0
      //                                 and  b1_i = d f_i / d c1
      //
      //  where the derivatives are partial derivatives
      //
      //  And we just differentiate the polynomials above (keeping the same level
      //  of truncation) and reexpress that as Horner's rule
      // 
      //  This clearly also handles the case of a unit gauge as no c1, u etc appears in the 
      //  denominator and the arccos is never taken. In this case, we have the results in 
      //  the raw c0, c1 form and we don't need to flip signs and take complex conjugates.
      //
      //  (not CD) I checked the expressions below by taking the difference between the Horner forms
      //  below from the expanded forms (and their derivatives) above and checking for the
      //  differences to be zero. At this point in time maple seems happy.
      //  ==================================================================
          
      f[0].real = 1. - c0*c0/720.;
      f[0].imag = -(c0/6.)*(1. - (c1/20.)*(1. - (c1/42.))) ;
      
      f[1].real =  c0/24.*(1. - c1/15.*(1. - 3.*c1/112.)) ;
      f[1].imag =  1.-c1/6.*(1. - c1/20.*(1. - c1/42.)) - c0*c0/5040. ;
      
      f[2].real = 0.5*(-1. + c1/12.*(1. - c1/30.*(1. - c1/56.)) + c0*c0/20160.);
      f[2].imag = 0.5*(c0/60.*(1. - c1/21.*(1. - c1/48.)));
      
      if( do_bs ) {
	//  partial f0/ partial c0
	b2[0].real = -c0/360.;
	b2[0].imag =  -(1./6.)*(1.-(c1/20.)*(1.-c1/42.));
        
	// partial f0 / partial c1
	//
	b1[0].real = 0;
	b1[0].imag = (c0/120.)*(1.-c1/21.);
        
	// partial f1 / partial c0
	//
	b2[1].real = (1./24.)*(1.-c1/15.*(1.-3.*c1/112.));
	b2[1].imag = -c0/2520.;
	
        
	// partial f1 / partial c1
	b1[1].real = -c0/360.*(1. - 3.*c1/56. );
	b1[1].imag = -1./6.*(1.-c1/10.*(1.-c1/28.));
        
	// partial f2/ partial c0
	b2[2].real = 0.5*c0/10080.;
	b2[2].imag = 0.5*(  1./60.*(1.-c1/21.*(1.-c1/48.)) );
        
	// partial f2/ partial c1
	b1[2].real = 0.5*(  1./12.*(1.-(2.*c1/30.)*(1.-3.*c1/112.)) ); 
	b1[2].imag = 0.5*( -c0/1260.*(1.-c1/24.) );
        
      } // do_bs
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
	// 
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
      
      if( do_bs) {
	
	if( fabs(w) < 0.05 ) { 
	  xi1 = -( 1./3. - (1./30.)*w_sq*( 1. - (1./28.)*w_sq*( 1. - (1./54.)*w_sq ) ) );
	}
	else { 
	  xi1 = cos(w)/w_sq - sin(w)/(w_sq*w);
	}
      }

      cosu = cos(u);
      sinu = sin(u);
      cosw = cos(w);
      sinw = sin(w);
      sin2u = sin(2*u);
      cos2u = cos(2*u);
      ucosu = u*cosu;
      usinu = u*sinu;
      ucos2u = u*cos2u;
      usin2u = u*sin2u;
      
      denom = 9.*u_sq - w_sq;

      {
	Real subexp1, subexp2, subexp3;

	subexp1 = u_sq - w_sq;
	subexp2 = 8*u_sq*cosw;
	subexp3 = (3*u_sq + w_sq)*xi0;
	
	f[0].real = ( (subexp1)*cos2u + cosu*subexp2 + 2*usinu*subexp3 ) / denom ;
	f[0].imag = ( (subexp1)*sin2u - sinu*subexp2 + 2*ucosu*subexp3 ) / denom ;

      }
      {
	Real subexp;
	
	subexp = (3*u_sq -w_sq)*xi0;
	
	f[1].real = (2*(ucos2u - ucosu*cosw)+subexp*sinu)/denom;
	f[1].imag = (2*(usin2u + usinu*cosw)+subexp*cosu)/denom;
      }
      {
	Real subexp;

	subexp=3*xi0;
      
	f[2].real = (cos2u - cosu*cosw -usinu*subexp) /denom ;
	f[2].imag = (sin2u + sinu*cosw -ucosu*subexp) /denom ;
      }

      if( do_bs )
	{
	  {
	      Real subexp1, subexp2, subexp3;
	      //          r_1[0]=Double(2)*cmplx(u, u_sq-w_sq)*exp2iu
	      //          + 2.0*expmiu*( cmplx(8.0*u*cosw, -4.0*u_sq*cosw)
	      //              + cmplx(u*(3.0*u_sq+w_sq),9.0*u_sq+w_sq)*xi0 );
	      
	      subexp1 = u_sq - w_sq;
	      subexp2 = 8.*cosw + (3.*u_sq + w_sq)*xi0 ;
	      subexp3 = 4.*u_sq*cosw - (9.*u_sq + w_sq)*xi0 ;
	      
	      r_1_re[0] = 2.*(ucos2u - sin2u *(subexp1)+ucosu*( subexp2 )- sinu*( subexp3 ) );
	      r_1_im[0] = 2.*(usin2u + cos2u *(subexp1)-usinu*( subexp2 )- cosu*( subexp3 ) );
	      
	  }
	  {
	      Real subexp1, subexp2;

	      // r_1[1]=cmplx(2.0, 4.0*u)*exp2iu + expmiu*cmplx(-2.0*cosw-(w_sq-3.0*u_sq)*xi0,2.0*u*cosw+6.0*u*xi0);
	      
	      subexp1 = cosw+3.*xi0;
	      subexp2 = 2.*cosw + xi0*(w_sq - 3.*u_sq);
	      
	      r_1_re[1] = 2.*((cos2u - 2.*usin2u) + usinu*subexp1) - cosu*subexp2;
	      r_1_im[1] = 2.*((sin2u + 2.*ucos2u) + ucosu*subexp1) + sinu*subexp2;
          }
	  {
	    Real subexp;
	    // r_1[2]=2.0*timesI(exp2iu)  +expmiu*cmplx(-3.0*u*xi0, cosw-3*xi0);
	    
	    subexp = cosw - 3.*xi0;
	    r_1_re[2] = -2.*sin2u -3.*ucosu*xi0 + sinu*subexp;
	    r_1_im[2] = 2.*cos2u  +3.*usinu*xi0 + cosu*subexp;
	  }
          
	  {
	    Real subexp;
	    //r_2[0]=-2.0*exp2iu + 2*cmplx(0,u)*expmiu*cmplx(cosw+xi0+3*u_sq*xi1,
	    //                                                 4*u*xi0);
	    
	    subexp = cosw + xi0 + 3.*u_sq*xi1;
	    r_2_re[0] = -2.*(cos2u + u*( 4.*ucosu*xi0 - sinu*subexp) );
	    r_2_im[0] = -2.*(sin2u - u*( 4.*usinu*xi0 + cosu*subexp) );
	  }
	  {
	    Real subexp;
          
	    // r_2[1]= expmiu*cmplx(cosw+xi0-3.0*u_sq*xi1, 2.0*u*xi0);
	    // r_2[1] = timesMinusI(r_2[1]);
	    
	    subexp =  cosw + xi0 - 3.*u_sq*xi1;
	    r_2_re[1] =  2.*ucosu*xi0 - sinu*subexp;
	    r_2_im[1] = -2.*usinu*xi0 - cosu*subexp;
	    
	  }
	  {
	    Real subexp;
	    //r_2[2]=expmiu*cmplx(xi0, -3.0*u*xi1);
	    
	    subexp = 3.*xi1;
            
	    r_2_re[2] =    cosu*xi0 - usinu*subexp ;
	    r_2_im[2] = -( sinu*xi0 + ucosu*subexp ) ;
	  }
          
	  b_denom=2.*denom*denom;
          
	  {
	    Real subexp1, subexp2, subexp3;
	    int j;

	    subexp1 = 2.*u;
	    subexp2 = 3.*u_sq - w_sq;
	    subexp3 = 2.*(15.*u_sq + w_sq);
	    
	    for(j=0; j < 3; j++) { 
	      
	      b1[j].real=( subexp1*r_1_re[j] + subexp2*r_2_re[j] - subexp3*f[j].real )/b_denom;
	      b1[j].imag=( subexp1*r_1_im[j] + subexp2*r_2_im[j] - subexp3*f[j].imag )/b_denom;
	    }
	  }
	  {
	    Real subexp1, subexp2;
	    int j;
	    
	    subexp1 = 3.*u;
	    subexp2 = 24.*u;
	    
	    for(j=0; j < 3; j++) { 
	      b2[j].real=( r_1_re[j] - subexp1*r_2_re[j] - subexp2 * f[j].real )/b_denom;
	      b2[j].imag=( r_1_im[j] - subexp1*r_2_im[j] - subexp2 * f[j].imag )/b_denom;
	    }
	  }
	  
	  // Now flip the coefficients of the b-s
	  if( c0 < 0 ) 
	    {
	      //b1_site[0] = conj(b1_site[0]);
	      b1[0].imag *= -1;
	      
	      //b1_site[1] = -conj(b1_site[1]);
	      b1[1].real *= -1;
	      
	      //b1_site[2] = conj(b1_site[2]);
	      b1[2].imag *= -1;
	      
	      //b2_site[0] = -conj(b2_site[0]);
	      b2[0].real *= -1;
	      
	      //b2_site[1] = conj(b2_site[1]);
	      b2[1].imag *= -1;
	      
	      //b2_site[2] = -conj(b2_site[2]);
	      b2[2].real *= -1;
	    }
	} // end of if (do_bs)
      
      // Now when everything is done flip signs of the b-s (can't do this before
      // as the unflipped f-s are needed to find the b-s
      
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

static void
quadr_comb( su3_matrix *T, su3_matrix *Q, su3_matrix *QQ, complex f[3])
{

  /*   T = f[0] + f[1]*Q + f[2]*QQ */

  clear_su3mat( T );
  T->e[0][0] = f[0];
  T->e[1][1] = f[0];
  T->e[2][2] = f[0];

  c_scalar_mult_add_su3mat( T, Q, &f[1], T );
  c_scalar_mult_add_su3mat( T, QQ, &f[2], T );
}

static void 
exp_iQ( su3_matrix *T, su3_matrix *Q )
{

  complex f[3];
  su3_matrix QQ;
  complex b1[3], b2[3];
  int do_bs = 0;

  mult_su3_nn( Q, Q, &QQ );
  get_fs_and_bs_from_Qs( f, b1, b2, Q, &QQ, do_bs);

#if 0
  {
    printf("f = (%.10f, %.10f) (%.10f, %.10f) (%.10f, %.10f)\n",
	   f[0].real,f[0].imag,f[1].real,f[1].imag,f[2].real,f[2].imag);
    get_fs_and_bs_from_Qs ( f, b1, b2, Q, &QQ, 1);
    printf("f = (%.10f, %.10f) (%.10f, %.10f) (%.10f, %.10f)\n",
	   f[0].real,f[0].imag,f[1].real,f[1].imag,f[2].real,f[2].imag);
    printf("b1 = (%.10f, %.10f) (%.10f, %.10f) (%.10f, %.10f)\n",
	   b1[0].real,b1[0].imag,b1[1].real,b1[1].imag,b1[2].real,b1[2].imag);
    printf("b2 = (%.10f, %.10f) (%.10f, %.10f) (%.10f, %.10f)\n",
	   b2[0].real,b2[0].imag,b2[1].real,b2[1].imag,b2[2].real,b2[2].imag);
  }
#endif

  quadr_comb( T, Q, &QQ, f);

}

/*--------------------------------------------------------------------*/

/* Do Morningstar-Peardon stout smearing to construct unitary W from
   the smeared link V and the unsmeared link U.

   Smearing applies to any scheme and not just the APE smearing
   described in MP.  For the general smearing case we still have 
 
      V = k*U + C 

   where C is the analog of the sum of APE terms in MP, U is the
   unsmeared link and k is any constant.  Note that Q in the first
   stout smearing step is the traceless hermitian part of -i V Uadj,
   which is the same as the traceless hermitian part of -i C Uadj, as
   used in MP. */

void stout_smear(su3_matrix *W, su3_matrix *V, su3_matrix *U)
{
  
  su3_matrix Q, tmp;

  get_Q_from_VUadj( &Q, V, U);

  /* tmp = exp(iQ) */
  exp_iQ( &tmp, &Q );

  /* W = exp(iQ) U */
  mult_su3_nn( &tmp, U, W);
}

/*--------------------------------------------------------------------*/

/* Compute force_diag and ILambda from the smeared link V, unsmeared
   link U and the force force_W for the W link (traceless antihermitian).
   These terms apply to any smearing V, not just APE as in MP.  The
   force_diag term is the first one in MP Eq (75).  The ILambdas are
   the insertions in the remaining terms.

   Once force_diag and ILambda are computed for all links, the
   subsidiary force is then computed by adding force_diag to the sum
   of all path contributions. The path contributions are obtained by
   walking the nontrivial paths in the smearing V, replacing each link
   on the path with ILambda, one at a time.  CHECK SIGNS WHEN YOU DO THIS! */

void stout_force_terms(su3_matrix *force_diag, su3_matrix *ILambda, 
		       su3_matrix *V, su3_matrix *U, su3_matrix *force_W)
{
  su3_matrix Q, QQ, USigp;
  su3_matrix Gamma;
  complex f[3], b1[3], b2[3];
  complex plusI = cmplx(0.,1.);
  int do_bs = 1;
  

  get_Q_from_VUadj( &Q, V, U);
  mult_su3_nn( &Q, &Q, &QQ );

  get_fs_and_bs_from_Qs( f, b1, b2, &Q, &QQ, do_bs);

  {
    int i, j;
    printf("Result Q\n");
    
    for(i = 0; i < 3; i++){
      for(j = 0; j < 3; j++)
	printf("%f + %f*I, ",Q.e[i][j].real, Q.e[i][j].imag);
      printf("\n");
    }
  }

  {
    su3_matrix tmp;

    /* tmp = exp(iQ) */
    quadr_comb( &tmp, &Q, &QQ, f);

    {
      int i, j;
      printf("Result exp(iQ)\n");
      
      for(i = 0; i < 3; i++){
	for(j = 0; j < 3; j++)
	  printf("%f + %f*I, ",tmp.e[i][j].real, tmp.e[i][j].imag);
	printf("\n");
      }
    }

    /* Note force_W is U' Sigma' = exp(iQ) U Sigma' in MP notation */
    /* force_diag = exp(-iQ) * force_W * exp(iQ) = first term in MP (75) */
    mult_su3_an( &tmp, force_W, &USigp );
    mult_su3_nn( &USigp, &tmp, force_diag );
  }

  /* Construction of Gamma from MP Eq (74) */

  {
    complex tr[3];
    su3_matrix B1, B2;

    /* B1 = b1[0] + b1[1]*Q + b1[2]*Q^2 */
    /* B2 = b2[0] + b2[1]*Q + b2[2]*Q^2 */
    quadr_comb( &B1, &Q, &QQ, b1 );
    quadr_comb( &B2, &Q, &QQ, b2 );

    /* Gamma = tr(U * Sigma' * B1) Q + tr(U * Sigma' * B2) Q^2 */
    tr[0] = cmplx(0,0);
    tr[1] = complextrace_su3_nn( &USigp, &B1);
    tr[2] = complextrace_su3_nn( &USigp, &B2);
    
    quadr_comb( &Gamma, &Q, &QQ, tr );
  }
  {
    su3_matrix tmp;

    /* Gamma += f_1 * U * Sigma' */
    c_scalar_mult_add_su3mat( &Gamma, &USigp, &f[1], &Gamma);
    
    /* Gamma += f_2 * Q * U * Sigma' */
    mult_su3_nn( &Q, &USigp, &tmp );
    c_scalar_mult_add_su3mat( &Gamma, &tmp, &f[2], &Gamma);
    
    /* Gamma += f_2 * U * Sigma' * Q */
    mult_su3_nn( &USigp, &Q, &tmp );
    c_scalar_mult_add_su3mat( &Gamma, &tmp, &f[2], &Gamma);
  }

  /* Lambda is the traceless-hermitian part of I*Gamma */

  traceless_hermitian_su3( ILambda, &Gamma);

  /* Multiply by I to make traceless antihermitian */

  c_scalar_mult_su3mat( ILambda, &plusI, ILambda );

} /* stout_force_terms */

int main(){
  su3_matrix FW, U, V, W, Fdiag, ILambda;

  int i,j;
  
  printf("Enter V\n");

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      scanf("%lf %lf",&V.e[i][j].real, &V.e[i][j].imag);

  printf("Enter U\n");

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      scanf("%lf %lf",&U.e[i][j].real, &U.e[i][j].imag);

  printf("Enter FW\n");

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      scanf("%lf %lf",&FW.e[i][j].real, &FW.e[i][j].imag);

  stout_smear( &W, &V, &U);

  printf("Result W\n");

  for(i = 0; i < 3; i++){
    for(j = 0; j < 3; j++)
      printf("(%f, %f) ",W.e[i][j].real, W.e[i][j].imag);
    printf("\n");
  }

  stout_force_terms(&Fdiag, &ILambda, &V, &U, &FW);

  printf("Result Fdiag\n");

  for(i = 0; i < 3; i++){
    for(j = 0; j < 3; j++)
      printf("(%f, %f) ",Fdiag.e[i][j].real, Fdiag.e[i][j].imag);
    printf("\n");
  }

  printf("Result ILambda\n");

  for(i = 0; i < 3; i++){
    for(j = 0; j < 3; j++)
      printf("(%f, %f) ",ILambda.e[i][j].real, ILambda.e[i][j].imag);
    printf("\n");
  }

  
}
