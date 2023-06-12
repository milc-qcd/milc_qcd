/************************* integrate.c **********************/
/* Integrators and associated utilities for the Wilson Flow */

#include "wilson_flow_includes.h"

/* Highest level function that integrates the flow, called from control.c */
void
run_gradient_flow() {

  /* RK integration variables */
  int i;
  double flowtime;

  /* Wilson flow output variables */
  double Et_C, Es_C, Et_W, Es_W, Et_S, Es_S, charge;
  double old_value=0, new_value=0;
  double der_value=0;

  /* Print flow output column labels */
  node0_printf("#LABEL time Clover_t Clover_s Plaq_t Plaq_s Rect_t Rect_s charge\n");
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  node0_printf("#LABEL2 time stepsize distance local_tol/distance\n");
#endif
  fflush(stdout);

  /* Calculate and print initial flow output */
  fmunu_fmunu(&Et_C, &Es_C, &charge);
  gauge_action_w_s( &Et_W, &Es_W, &Et_S, &Es_S );
#if (MILC_PRECISION==1)
  node0_printf("GFLOW: %g %g %g %g %g %g %g %g\n", 0.0, Et_C, Es_C, Et_W, Es_W, Et_S, Es_S, charge);
#else
  node0_printf("GFLOW: %g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n", 0.0, Et_C, Es_C, Et_W, Es_W, Et_S, Es_S, charge);
#endif
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
#if (MILC_PRECISION==1)
    node0_printf("ADAPT: %g %g %g %g\n", 0.0, stepsize, 0.0, 0.0 );
#else
    node0_printf("ADAPT: %.16g %.16g %.16g %.16g\n", 0.0, stepsize, 0.0, 0.0 );
#endif
#endif
  fflush(stdout);

#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  steps_rejected = 0; // count rejected steps in adaptive schemes
#endif
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  is_first_step = 1; // need to know the first step for FSAL
  // set the permutation array for FSAL, this saves copying
  // K[3] to K[0] after each step
  indK[0] = 0; indK[1] = 1; indK[2] = 2; indK[3] = 3;
#endif
  is_final_step = 0;
  flowtime = 0;
  i = 0;
  /* Loop over the flow time */
  while( stoptime==AUTO_STOPTIME || ( flowtime<stoptime && is_final_step==0 ) ) {
    /* Adjust last time step to fit exactly stoptime */
    if( stepsize>stoptime-flowtime && stoptime!=AUTO_STOPTIME ) {
      stepsize = stoptime-flowtime;
      is_final_step = 1;
    }
//      printf("%g\n", stepsize);

    /* Perform one flow step (most of the computation is here) */
    flow_step();

    flowtime += stepsize;
    i++;

    /* Calculate and print current flow output */
    fmunu_fmunu(&Et_C, &Es_C, &charge);
    gauge_action_w_s( &Et_W, &Es_W, &Et_S, &Es_S );
#if (MILC_PRECISION==1)
    node0_printf("GFLOW: %g %g %g %g %g %g %g %g\n", flowtime, Et_C, Es_C, Et_W, Es_W, Et_S, Es_S, charge);
#else
    node0_printf("GFLOW: %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n", flowtime, Et_C, Es_C, Et_W, Es_W, Et_S, Es_S, charge);
#endif
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
#if (MILC_PRECISION==1)
    node0_printf("ADAPT: %g %g %g %g\n", flowtime, stepsize, dist, local_tol/dist );
#else
    node0_printf("ADAPT: %g %.16g %.16g %.16g\n", flowtime, stepsize, dist, local_tol/dist );
#endif
#endif
    fflush(stdout);

    /* Automatic determination of stoptime:                         */
    /*  t^2 E > 0.45 and d/dt { t^2 E } > 0.35                      */
    /*  Bounds need to be adjusted with scale determination cutoff  */
    if( stoptime==AUTO_STOPTIME ) {

      old_value = new_value;
      new_value = flowtime*flowtime*(Et_C+Es_C);
      der_value = flowtime*(new_value-old_value)/stepsize;

      if( new_value > 0.45 && der_value > 0.35 )
        break;
    } /* end: auto stoptime */

#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
    if( is_final_step==0 ) {
      // adjust step size for the next step except if it is final
      stepsize = stepsize * SAFETY * pow( local_tol/dist, 1/3. );
    }
#endif

  } /* end: flowtime loop */

  /* Save and print the number of steps */
  total_steps = i;
  node0_printf("Number of steps = %i\n", total_steps);
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  node0_printf("Number of rejected steps = %i\n", steps_rejected);
#endif
  fflush(stdout);

}




/* Resets all entries in an anti-hermition matrix to zero */
void
clear_anti_hermitian( anti_hermitmat *dest )
{
  dest->m01.real = 0;
  dest->m01.imag = 0;
  dest->m02.real = 0;
  dest->m02.imag = 0;
  dest->m12.real = 0;
  dest->m12.imag = 0;
  dest->m00im = 0;
  dest->m11im = 0;
  dest->m22im = 0;
}

/* Sets the provided matrix to the identity */
void
set_identity( su3_matrix *dest )
{
  /* Set diagonals to 1 */
  dest->e[0][0].real = 1.;
  dest->e[0][0].imag = 0.;
  dest->e[1][1].real = 1.;
  dest->e[1][1].imag = 0.;
  dest->e[2][2].real = 1.;
  dest->e[2][2].imag = 0.;

  /* Set off diagonals to 0 */
  dest->e[0][1].real = 0.;
  dest->e[0][1].imag = 0.;
  dest->e[0][2].real = 0.;
  dest->e[0][2].imag = 0.;
  dest->e[1][0].real = 0.;
  dest->e[1][0].imag = 0.;
  dest->e[1][2].real = 0.;
  dest->e[1][2].imag = 0.;
  dest->e[2][0].real = 0.;
  dest->e[2][0].imag = 0.;
  dest->e[2][1].real = 0.;
  dest->e[2][1].imag = 0.;
}

/* Multiplies an antihermitian matrix with an su3 matrix, in that order */
void
mult_ah_su3_nn( anti_hermitmat *a, su3_matrix *b, su3_matrix *dest )
{
  int j;
  register Real a0r,a0i,a1r,a1i,a2r,a2i;
  register Real b0r,b0i,b1r,b1i,b2r,b2i;

  for(j=0;j<3;j++){

    a0i=a->m00im; //a0r is 0
    b0r=b->e[0][j].real; b0i=b->e[0][j].imag;
    a1r=a->m01.real; a1i=a->m01.imag;
    b1r=b->e[1][j].real; b1i=b->e[1][j].imag;
    a2r=a->m02.real; a2i=a->m02.imag;
    b2r=b->e[2][j].real; b2i=b->e[2][j].imag;

    dest->e[0][j].real = a1r*b1r - a0i*b0i - a1i*b1i + a2r*b2r - a2i*b2i;
    dest->e[0][j].imag = a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

    a0r=a->m01.real; a0i=a->m01.imag; //negate a0r
    b0r=b->e[0][j].real; b0i=b->e[0][j].imag;
    a1i=a->m11im; //a1r is 0
    b1r=b->e[1][j].real; b1i=b->e[1][j].imag;
    a2r=a->m12.real; a2i=a->m12.imag;
    b2r=b->e[2][j].real; b2i=b->e[2][j].imag;

    dest->e[1][j].real = a2r*b2r - a0i*b0i - a1i*b1i - a0r*b0r - a2i*b2i;
    dest->e[1][j].imag = a0i*b0r - a0r*b0i + a1i*b1r + a2r*b2i + a2i*b2r;

    a0r=a->m02.real; a0i=a->m02.imag; //negate a0r
    b0r=b->e[0][j].real; b0i=b->e[0][j].imag;
    a1r=a->m12.real; a1i=a->m12.imag; //negate a1r
    b1r=b->e[1][j].real; b1i=b->e[1][j].imag;
    a2i=a->m22im; //a2r is 0
    b2r=b->e[2][j].real; b2i=b->e[2][j].imag;

    dest->e[2][j].real = -1*a0r*b0r - a0i*b0i - a1r*b1r - a1i*b1i - a2i*b2i;
    dest->e[2][j].imag = a0i*b0r - a0r*b0i - a1r*b1i + a1i*b1r + a2i*b2r;
  }
}

/* Multiplies an antihermitian matrix by a real scalar */
void
scalar_mult_ah( anti_hermitmat *a, Real c, anti_hermitmat *dest )
{
  dest->m01.real = c*a->m01.real;
  dest->m01.imag = c*a->m01.imag;
  dest->m02.real = c*a->m02.real;
  dest->m02.imag = c*a->m02.imag;
  dest->m12.real = c*a->m12.real;
  dest->m12.imag = c*a->m12.imag;
  dest->m00im = c*a->m00im;
  dest->m11im = c*a->m11im;
  dest->m22im = c*a->m22im;
}

/* Adds a matrix times a real scalar to a matrix (all antihermitian) */
void
scalar_mult_add_ah( anti_hermitmat *a, anti_hermitmat *b, Real c,
                    anti_hermitmat *dest )
{
  dest->m01.real =  a->m01.real + c*b->m01.real;
  dest->m01.imag =  a->m01.imag + c*b->m01.imag;
  dest->m02.real =  a->m02.real + c*b->m02.real;
  dest->m02.imag =  a->m02.imag + c*b->m02.imag;
  dest->m12.real =  a->m12.real + c*b->m12.real;
  dest->m12.imag =  a->m12.imag + c*b->m12.imag;
  dest->m00im =  a->m00im + c*b->m00im;
  dest->m11im =  a->m11im + c*b->m11im;
  dest->m22im =  a->m22im + c*b->m22im;
}

/* Computes the traceless, antihermition projection of a matrix */
/*  B = (A-A.dag)/2 - Tr{ (A-A.dag)/2 }                         */
void
anti_hermitian_traceless_proj( su3_matrix *a, anti_hermitmat *dest )
{
  /* Reused for diagonal elements and trace */
  register Real a00im, a11im, a22im;
  register Real tr3;

  /* Compute off-diagonal elements (no trace needed) */
  dest->m01.real = (a->e[0][1].real - a->e[1][0].real)*0.5;
  dest->m01.imag = (a->e[0][1].imag + a->e[1][0].imag)*0.5;
  dest->m02.real = (a->e[0][2].real - a->e[2][0].real)*0.5;
  dest->m02.imag = (a->e[0][2].imag + a->e[2][0].imag)*0.5;
  dest->m12.real = (a->e[1][2].real - a->e[2][1].real)*0.5;
  dest->m12.imag = (a->e[1][2].imag + a->e[2][1].imag)*0.5;

  /* Compute 1/3 of the trace of the antihermitian projection */
  a00im = a->e[0][0].imag;
  a11im = a->e[1][1].imag;
  a22im = a->e[2][2].imag;
  tr3 = (a00im + a11im + a22im)*0.33333333333333333;
  dest->m00im = a00im - tr3;
  dest->m11im = a11im - tr3;
  dest->m22im = a22im - tr3;
}

/* Approximates the exponential of an anti-hermition matrix */
/*  Uses a taylor series expansion about 0, to order n      */
void
exp_anti_hermitian( anti_hermitmat *a, su3_matrix *dest, int n )
{
  register int order;
  su3_matrix identity, temp1, temp2;

  /* Initialize the identity and exponential approx */
  set_identity(&identity);
  set_identity(&temp1);

  /* Loop over expansion of exponential starting at the end  */
  /*  exp(a) = I + a/1*(I + a/2*(I + a/3*(...(I + a/n)...))) */
  for(order=n; order>0; order--) {
    mult_ah_su3_nn(a, &temp1, &temp2);
    scalar_mult_add_su3_matrix(&identity, &temp2, 1./order, &temp1);
  }

  /* Copy the result into destination */
  su3mat_copy(&temp1, dest);
}


/* copy antihermitian matrix: a->b */
void ahmat_copy( anti_hermitmat *a, anti_hermitmat *b ) {
  b->m01 = a->m01;
  b->m02 = a->m02;
  b->m12 = a->m12;
  b->m00im = a->m00im;
  b->m11im = a->m11im;
  b->m22im = a->m22im;
}

//#define USE_SLOW_COMMUTATOR_AH
#ifdef USE_SLOW_COMMUTATOR_AH
/* commutator of anti-Hermitian matrices,
   this is a slow version that relies on uncompressing the matrices
   and generic multiplication */
void
commutator_ah( anti_hermitmat *a, anti_hermitmat *b, anti_hermitmat *c ) {

  su3_matrix temp1, temp2, temp3, temp4;

  uncompress_anti_hermitian( a, &temp1 );
  uncompress_anti_hermitian( b, &temp2 );

  mult_su3_nn( &temp1, &temp2, &temp3 );
  mult_su3_nn( &temp2, &temp1, &temp4 );
  sub_su3_matrix( &temp3, &temp4, &temp1 );
  compress_anti_hermitian( &temp1, c );
}
#else
/* commutator of anti-Hermitian matrices,
   direct calculation */
void
commutator_ah( anti_hermitmat *a, anti_hermitmat *b, anti_hermitmat *c ) {

  Real temp01r, temp02r, temp12r;
  Real temp01i, temp02i, temp12i;

  temp01r  = b->m00im*a->m01.imag-a->m00im*b->m01.imag;
  temp01r += b->m01.imag*a->m11im-a->m01.imag*b->m11im;
  temp01r += b->m02.real*a->m12.real-a->m02.real*b->m12.real;
  temp01r += b->m02.imag*a->m12.imag-a->m02.imag*b->m12.imag;

  temp02r  = b->m00im*a->m02.imag-a->m00im*b->m02.imag;
  temp02r += a->m01.real*b->m12.real-b->m01.real*a->m12.real;
  temp02r += b->m02.imag*a->m22im-a->m02.imag*b->m22im;
  temp02r += b->m01.imag*a->m12.imag-a->m01.imag*b->m12.imag;

  temp12r  = b->m11im*a->m12.imag-a->m11im*b->m12.imag;
  temp12r += b->m01.real*a->m02.real-a->m01.real*b->m02.real;
  temp12r += b->m12.imag*a->m22im-a->m12.imag*b->m22im;
  temp12r += b->m01.imag*a->m02.imag-a->m01.imag*b->m02.imag;

  temp01i  = a->m00im*b->m01.real-b->m00im*a->m01.real;
  temp01i += a->m01.real*b->m11im-b->m01.real*a->m11im;
  temp01i += b->m02.imag*a->m12.real-a->m02.imag*b->m12.real;
  temp01i += a->m02.real*b->m12.imag-b->m02.real*a->m12.imag;

  temp02i  = a->m00im*b->m02.real-b->m00im*a->m02.real;
  temp02i += a->m02.real*b->m22im-b->m02.real*a->m22im;
  temp02i += a->m01.imag*b->m12.real-b->m01.imag*a->m12.real;
  temp02i += a->m01.real*b->m12.imag-b->m01.real*a->m12.imag;

  temp12i  = a->m11im*b->m12.real-b->m11im*a->m12.real;
  temp12i += a->m12.real*b->m22im-b->m12.real*a->m22im;
  temp12i += a->m01.imag*b->m02.real-b->m01.imag*a->m02.real;
  temp12i += b->m01.real*a->m02.imag-a->m01.real*b->m02.imag;


  c->m00im  = b->m01.imag*a->m01.real-a->m01.imag*b->m01.real;
  c->m00im += b->m02.imag*a->m02.real-a->m02.imag*b->m02.real;
  c->m00im *= 2;
  c->m11im  = a->m01.imag*b->m01.real-b->m01.imag*a->m01.real;
  c->m11im += b->m12.imag*a->m12.real-a->m12.imag*b->m12.real;
  c->m11im *= 2;
  c->m22im  = a->m02.imag*b->m02.real-b->m02.imag*a->m02.real;
  c->m22im += a->m12.imag*b->m12.real-b->m12.imag*a->m12.real;
  c->m22im *= 2;
  c->m01.real = temp01r;
  c->m01.imag = temp01i;
  c->m02.real = temp02r;
  c->m02.imag = temp02i;
  c->m12.real = temp12r;
  c->m12.imag = temp12i;
}
#endif

/* inverse derivative of the matrix exponential,
   required for generic RKMK methods */
void
dexpinv( anti_hermitmat *u, anti_hermitmat *v, int q, anti_hermitmat *d ) {
  // Bernoulli numbers normalized with k!, i.e. this array is B_k/k!
  Real BernoulliK[11] = { 1, -1/2., 1/12., 0, -1/720., 0, 1/30240., 0, -1/1209600., 0, 1/47900160. };

  anti_hermitmat w;
  int register k;

  ahmat_copy( v, &w );
  ahmat_copy( v, d );
  for( k=1; k<q; k++ ) {
    commutator_ah( u, &w, &w );
    if( BernoulliK[k]==0 ) continue;
    scalar_mult_add_ah( d, &w, BernoulliK[k], d );
  }
}


//#define USE_STRICT_DISTANCE
#ifdef USE_STRICT_DISTANCE
/* distance between SU(3) matrices:
   maximum difference element-wise,
   real and imaginary parts are treated separately,
   this is the strictest criterium */
Real
su3mat_distance( su3_matrix *a, su3_matrix *b ) {

  Real dmax = 0, temp;
  int register i, j;

  for( i=0; i<3; i++ ) {
    for( j=0; j<3; j++ ) {
      temp = fabs(a->e[i][j].real-b->e[i][j].real);
      if( dmax<temp ) dmax = temp;
      temp = fabs(a->e[i][j].imag-b->e[i][j].imag);
      if( dmax<temp ) dmax = temp;
    }
  }
  return dmax;
}
#else
/* distance between SU(3) matrices:
   normalized root of the average element-wise
   distance squared -- milder and smoother criteria,
   defined in Ramos, Fritzsch, 1301.4388 */
Real
su3mat_distance( su3_matrix *a, su3_matrix *b ) {

  Real temp = 0, re, im;
  int register i, j;

  for( i=0; i<3; i++ ) {
    for( j=0; j<3; j++ ) {
      re = a->e[i][j].real-b->e[i][j].real;
      im = a->e[i][j].imag-b->e[i][j].imag;
      temp += re*re + im*im;
    }
  }
  // NOTE: the normalization in 1301.4388
  // is 1/N_c^2 outside the square root
  temp = sqrt(temp) / 9;
  return temp;
}
#endif


#if GF_INTEGRATOR==INTEGRATOR_LUSCHER || GF_INTEGRATOR==INTEGRATOR_CK \
 || GF_INTEGRATOR==INTEGRATOR_BBB || GF_INTEGRATOR==INTEGRATOR_CF3
/* A single step for a 2N-storage Runge-Kutta scheme
 * where the right hand side of the flow equation is evaluated
 * and the fields are updated
 *  A: accumulating matrix over all smear steps in a single time step
 *  S: staple (action dependent); recalculated before each smear step
 *  U: gauge links
 *  cA,cB: constants
 *    Calculating a single smear step is done by
 *    (this follows usual convention on 2N-storage RK schemes):
 *    A = cA*A + proj(S*U) -> update accumulation matrix
 *    U = exp(cB*A)*U   -> update gauge links
 */
void
integrate_RK_2N_one_stage( Real cA, Real cB )
{
  register int dir, i;
  register site *s;

  /* Temporary matrix holders */
  anti_hermitmat *Acur, tempA1;
  su3_matrix *U, tempS1, tempS2;

  /* Calculate the new staple */
  staple();

  FORALLUPDIR(dir)
    FORALLSITES(i, s) {
      /* Retrieve the current link and accumulation matrix */
      U = &(s->link[dir]);
      Acur = &(s->accumulate[dir]);

      /* Update the accumulation matrix A = cA*A + proj(U*S) */
      mult_su3_na( U, &(s->staple[dir]), &tempS1 );
      anti_hermitian_traceless_proj( &tempS1, &tempA1 );
      scalar_mult_add_ah( &tempA1, Acur, cA, Acur );

      /* Update the links U = exp(cB*A)*U */
      scalar_mult_ah( Acur, cB, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, exp_order );
      mult_su3_nn( &tempS1, U, &tempS2 );
      su3mat_copy( &tempS2, U );
  }
}

/* one step of a low-storage Runge-Kutta scheme,
   this includes Luscher, arXiv:1006.4518 [hep-lat]
   or any other 2N-storage scheme */
void
integrate_RK_2N()
{
  register int dir, i;
  register site *s;

  /* Clear the accumulation matrix */
  FORALLSITES(i, s)
    FORALLUPDIR(dir)
      clear_anti_hermitian(&(s->accumulate[dir]));

  /* Infinitesimal stout smearing */
  for( i=0; i<N_stages; i++ ) {
    /* be careful with stepsize: Nathan's convention on the staple
       is such that stepsize should be taken negative */
    integrate_RK_2N_one_stage( A_2N[i], -B_2N[i]*stepsize );
  }
}
#elif GF_INTEGRATOR==INTEGRATOR_RKMK4 || GF_INTEGRATOR==INTEGRATOR_RKMK5 || GF_INTEGRATOR==INTEGRATOR_RKMK8
/* generic integrator of Runge-Kutta-Munthe-Kaas type,
   requires dexpinv evaluation at each step (nested commutators) */
void
integrate_RKMK_generic() {

  register int dir, i, i_rk, j_rk;
  register site *s;
  su3_matrix tempS1, tempS2;
  anti_hermitmat tempA1, *Atemp;

  // store the initial state of the gauge field,
  // it is used at every stage in this integrator format
  FORALLUPDIR(dir)
    FORALLSITES(i, s)
      su3mat_copy( &(s->link[dir]), &(s->link0[dir]) );

  // loop over RK stages
  for( i_rk=0; i_rk<N_stages; i_rk++ ) {
    if( i_rk!=0 ) {
      FORALLUPDIR(dir)
        FORALLSITES(i, s) {

          Atemp = &(s->accumulate[dir]);
          clear_anti_hermitian( Atemp );
          for( j_rk=0; j_rk<i_rk; j_rk++ ) {
            // accumulate a_i1*K1 + a_i2*K2 + ...
            scalar_mult_add_ah( Atemp, &(s->K[j_rk][dir]), a_RK[i_rk][j_rk], Atemp );
          }
          // update the link
          scalar_mult_ah( Atemp, -stepsize, Atemp );
          exp_anti_hermitian( Atemp, &tempS1, exp_order );
          mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
      }
    }

    // get the right hand side of the flow equation from the staple
    // and store in s->K[i_rk]
    staple();

    FORALLUPDIR(dir)
      FORALLSITES(i, s) {

        mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
        anti_hermitian_traceless_proj( &tempS1, &(s->K[i_rk][dir]) );
        if( i_rk!=0 )
          dexpinv( &(s->accumulate[dir]), &(s->K[i_rk][dir]), p_order, &(s->K[i_rk][dir]));
    }
  }
  // final RK stage
  FORALLUPDIR(dir)
    FORALLSITES(i, s) {
      clear_anti_hermitian( &tempA1 );
      for( i_rk=0; i_rk<N_stages; i_rk++ ) {
        // accumulate b_1*K1 + b_2*K2 + ...
        scalar_mult_add_ah( &tempA1, &(s->K[i_rk][dir]), b_RK[i_rk], &tempA1 );
      }
      // update the link
      scalar_mult_ah( &tempA1, -stepsize, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, exp_order );
      mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
  }
}
#elif GF_INTEGRATOR==INTEGRATOR_RKMK3
/* Third order Runge-Kutta-Munthe-Kaas type,
   requires a single commutator at the last stage */
void
integrate_RKMK3() {

  register int dir, i, i_rk, j_rk;
  register site *s;
  su3_matrix tempS1, tempS2;
  anti_hermitmat tempA1, tempA2, *Atemp;

  // store the initial state of the gauge field,
  // it is used at every stage in this integrator format
  FORALLUPDIR(dir)
    FORALLSITES(i, s)
      su3mat_copy( &(s->link[dir]), &(s->link0[dir]) );

  // loop over RK stages
  for( i_rk=0; i_rk<N_stages; i_rk++ ) {
    if( i_rk!=0 ) {
      FORALLUPDIR(dir)
        FORALLSITES(i, s) {

          Atemp = &(s->accumulate[dir]);
          clear_anti_hermitian( Atemp );
          for( j_rk=0; j_rk<i_rk; j_rk++ ) {
            // accumulate a_i1*K1 + a_i2*K2 + ...
            scalar_mult_add_ah( Atemp, &(s->K[j_rk][dir]), a_RK[i_rk][j_rk], Atemp );
          }
          // update the link
          scalar_mult_ah( Atemp, -stepsize, Atemp );
          exp_anti_hermitian( Atemp, &tempS1, exp_order );
          mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
      }
    }

    // get the right hand side of the flow equation from the staple
    // and store in s->K[i_rk]
    staple();

    FORALLUPDIR(dir)
      FORALLSITES(i, s) {

        mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
        anti_hermitian_traceless_proj( &tempS1, &(s->K[i_rk][dir]) );
    }
  }
  // final RK stage
  FORALLUPDIR(dir)
    FORALLSITES(i, s) {
      clear_anti_hermitian( &tempA1 );
      for( i_rk=0; i_rk<N_stages; i_rk++ ) {
        // accumulate b_1*K1 + b_2*K2 + ...
        scalar_mult_add_ah( &tempA1, &(s->K[i_rk][dir]), b_RK[i_rk], &tempA1 );
      }
      // the only commutator in this scheme is at the last stage
      commutator_ah( &tempA1, &(s->K[0][dir]), &tempA2 );
      scalar_mult_add_ah( &tempA1, &tempA2, -stepsize/6, &tempA1 );
      // update the link
      scalar_mult_ah( &tempA1, -stepsize, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, exp_order );
      mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
  }
}
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3
void
integrate_adapt_RK_2N_one_stage( Real cA, Real cB, int istep )
{
  register int dir, i;
  register site *s;

  /* Temporary matrix holders */
  anti_hermitmat *Acur, tempA1;
  su3_matrix *U, tempS1, tempS2;

  /* Calculate the new staple */
  staple();

  FORALLUPDIR(dir)
    FORALLSITES(i, s) {
      /* Retrieve the current link and accumulation matrix */
      U = &(s->link[dir]);
      Acur = &(s->accumulate[dir]);

      /* Update the accumulation matrix A = cA*A + proj(U*S) */
      mult_su3_na( U, &(s->staple[dir]), &tempS1 );
      anti_hermitian_traceless_proj( &tempS1, &(s->K[istep][dir]) );
      scalar_mult_add_ah( &(s->K[istep][dir]), Acur, cA, Acur );

      /* Update the links U = exp(cB*A)*U */
      scalar_mult_ah( Acur, cB, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, exp_order );
      mult_su3_nn( &tempS1, U, &tempS2 );
      su3mat_copy( &tempS2, U );
  }
}

/* Adaptive scheme based on Luscher's in 2N-storage format */
void
integrate_adapt_RK_2N()
{
  register int dir, i;
  register site *s;
  int is_repeat = 1;
  anti_hermitmat tempA1;
  su3_matrix tempS1, tempS2;
  Real temp;

  FORALLSITES(i, s)
    FORALLUPDIR(dir) {
      /* Clear the accumulation matrix */
      clear_anti_hermitian(&(s->accumulate[dir]));
      /* Store the initial state of the gauge field */
      su3mat_copy( &(s->link[dir]), &(s->link0[dir]) );
  }

  do {
    /* Make one RK step */
    for( i=0; i<N_stages; i++ ) {
      /* be careful with stepsize: Nathan's convention on the staple
         is such that stepsize should be taken negative */
      integrate_adapt_RK_2N_one_stage( A_2N[i], -B_2N[i]*stepsize, i );
    }

    dist = 0;
    FORALLUPDIR(dir)
      FORALLSITES(i, s) {
        /* Construct lower order approximation */
        scalar_mult_ah( &(s->K[0][dir]), Lambda[0], &tempA1 );
        scalar_mult_add_ah( &tempA1, &(s->K[1][dir]), Lambda[1], &tempA1 );
        scalar_mult_add_ah( &tempA1, &(s->K[2][dir]), Lambda[2], &tempA1 );
        // NOTE: Nathan Brown's convention: stepsize is negative
        scalar_mult_ah( &tempA1, -stepsize, &tempA1 );
        exp_anti_hermitian( &tempA1, &tempS1, exp_order );
        mult_su3_nn( &tempS1, &(s->link0[dir]), &tempS2 );
        /* Calculate distance between the two approximations */
        temp = su3mat_distance( &(s->link[dir]), &tempS2 );
        /* Find the maximum over the local volume */
        if( dist<temp ) dist = temp;
    }
    /* Get the global maximum distance */
    g_floatmax( &dist );
    /* Check if tolerance is exceeded, redo the step
       except if it is final */
    if( dist>local_tol && is_final_step==0 ) {
      // adjust step size
      stepsize = stepsize * SAFETY * pow( local_tol/dist, 1/3. );
      // record failed step
      steps_rejected++;
      // copy over the original state of the gauge field
      FORALLSITES(i, s)
        FORALLUPDIR(dir)
          su3mat_copy( &(s->link0[dir]), &(s->link[dir]) );
    }
    else {
      is_repeat = 0;
    }
  } while( is_repeat==1 );
}
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
/* Bogacki-Shampine adaptive 3(2) embedded pair */
void
integrate_adapt_bs() {

  register int dir, i, i_rk, j_rk;
  register site *s;
  su3_matrix tempS1, tempS2;
  anti_hermitmat tempA1, tempA2, *Atemp;
  Real temp;
  int is_repeat = 1;

  // store the initial state of the gauge field,
  // it is used at every stage in this integrator format
  // and if the step gets rejected
  FORALLUPDIR(dir)
    FORALLSITES(i, s)
      su3mat_copy( &(s->link[dir]), &(s->link0[dir]) );

  // get the first force evaluation on the very first step of integration
  if( is_first_step==1 ) {
    staple();
    FORALLUPDIR(dir)
      FORALLSITES(i, s) {
        mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
        anti_hermitian_traceless_proj( &tempS1, &(s->K[indK[0]][dir]) );
    }
    is_first_step = 0;
  }

  do {
    // loop over RK stages, skip 0 due to FSAL
    for( i_rk=1; i_rk<N_stages; i_rk++ ) {
        FORALLUPDIR(dir)
          FORALLSITES(i, s) {

            Atemp = &(s->accumulate[dir]);
            clear_anti_hermitian( Atemp );
            for( j_rk=0; j_rk<i_rk; j_rk++ ) {
              // accumulate a_i1*K1 + a_i2*K2 + ...
              scalar_mult_add_ah( Atemp, &(s->K[indK[j_rk]][dir]), a_RK[i_rk][j_rk], Atemp );
            }
            // update the link
            scalar_mult_ah( Atemp, -stepsize, Atemp );
            exp_anti_hermitian( Atemp, &tempS1, exp_order );
            mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
        }
        // get the right hand side of the flow equation from the staple
        // and store in s->K[i_rk]
        // NOTE: here FSAL property is used, so the force in K[0] is
        //       already filled from the previous step
        staple();

        FORALLUPDIR(dir)
          FORALLSITES(i, s) {
            mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
            anti_hermitian_traceless_proj( &tempS1, &(s->K[indK[i_rk]][dir]) );
        }
    }
    // final RK stage that gives fourth-order local approximation
    FORALLUPDIR(dir)
      FORALLSITES(i, s) {
        clear_anti_hermitian( &tempA1 );
        for( i_rk=0; i_rk<N_stages; i_rk++ ) {
          // accumulate b_1*K1 + b_2*K2 + ...
          scalar_mult_add_ah( &tempA1, &(s->K[indK[i_rk]][dir]), b_RK[i_rk], &tempA1 );
        }
        // the only commutator in this scheme is at the last stage
        commutator_ah( &tempA1, &(s->K[indK[0]][dir]), &tempA2 );
        scalar_mult_add_ah( &tempA1, &tempA2, -stepsize/6, &tempA1 );
        // update the link
        scalar_mult_ah( &tempA1, -stepsize, &tempA1 );
        exp_anti_hermitian( &tempA1, &tempS1, exp_order );
        mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
    }
    // additional stage
    staple();
    dist = 0;
    FORALLUPDIR(dir)
      FORALLSITES(i, s) {
        mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
        anti_hermitian_traceless_proj( &tempS1, &(s->K[indK[3]][dir]) );
        clear_anti_hermitian( &tempA1 );
        for( i_rk=0; i_rk<4; i_rk++ ) {
          // accumulate b'_1*K1 + b'_2*K2 + b'_3*K3 + b'_4*K4
          // NOTE: b' coefficients are stored as a_RK[3][0], a_RK[3][1], etc.
          scalar_mult_add_ah( &tempA1, &(s->K[indK[i_rk]][dir]), a_RK[3][i_rk], &tempA1 );
        }
        // get the lower (third) order estimate
        scalar_mult_ah( &tempA1, -stepsize, &tempA1 );
        exp_anti_hermitian( &tempA1, &tempS1, exp_order );
        mult_su3_nn( &tempS1, &(s->link0[dir]), &tempS2 );
        /* Calculate distance between the two approximations */
        temp = su3mat_distance( &(s->link[dir]), &tempS2 );
        /* Find the maximum over the local volume */
        if( dist<temp ) dist = temp;
    }
    /* Get the global maximum distance */
    g_floatmax( &dist );
    /* Check if tolerance is exceeded, redo the step
       except if it is final */
    if( dist>local_tol && is_final_step==0 ) {
      // adjust step size
      stepsize = stepsize * SAFETY * pow( local_tol/dist, 1/3. );
      // record failed step
      steps_rejected++;
      // copy over the original state of the gauge field
      FORALLSITES(i, s)
        FORALLUPDIR(dir)
          su3mat_copy( &(s->link0[dir]), &(s->link[dir]) );
    }
    else {
      is_repeat = 0;
    }
  } while( is_repeat==1 );
  // permute indices to read the force on the next step
  i = indK[0];
  indK[0] = indK[3];
  indK[3] = i;
}
#endif



#ifdef DEBUG_FIELDS
#define REPACK_TO_DOUBLE
#ifdef REPACK_TO_DOUBLE
#define MATRIX_TYPE dsu3_matrix
#else
#define MATRIX_TYPE fsu3_matrix
#endif

void repack_site( su3_matrix *a, MATRIX_TYPE *b ) {
  int dir,i,j;

  for(dir = 0; dir < 4; dir++){
    for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
      b[dir].e[i][j].real = a[dir].e[i][j].real;
      b[dir].e[i][j].imag = a[dir].e[i][j].imag;
    }
  }
}

void dump_double_lattice() {

  int i, x, y, z, t, nxa, nya, nza, nta;
  const int *nsq;
  char filename[1024];
  MATRIX_TYPE *tbuf = NULL;
  char myname[] = "dump_double_lattice";
  FILE *fp;

  // local sizes
  nsq = get_logical_dimensions();
  nxa = nx/nsq[0];
  nya = ny/nsq[1];
  nza = nz/nsq[2];
  nta = nt/nsq[3];

  // make a node file name
  sprintf( filename, "fields_on_node.%04d", this_node );

  tbuf = (MATRIX_TYPE *)malloc(nxa*4*sizeof(MATRIX_TYPE));
  if(tbuf == NULL){
    printf("%s(%d): No room for tbuf\n",myname,this_node);
    terminate(1);
  }

  fp = fopen( filename, "wb" );

  // loop over fields and store
  for( t=0; t<nta; t++) for( z=0; z<nza; z++) for( y=0; y<nya; y++ ) {
    for( x=0; x<nxa; x++) {
      i = node_index( x, y, z, t );
      repack_site( &lattice[i].link[0], &tbuf[4*x] );
    }

    if( (int)fwrite( (void*)tbuf, sizeof(MATRIX_TYPE), 4*nxa, fp ) != 4*nxa )
    {
      printf("dump_double_lattice: Node %d gauge configuration write error\n",
             this_node);
      fflush(stdout);
      terminate(1);
    }

  }

  fclose( fp );

}
#endif
