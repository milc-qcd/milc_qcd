/************************* integrate.c **********************/
/* Integrators and associated utilities for the Wilson Flow */

#include "wilson_flow_includes.h"

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

/* Computes a single smearing steps (three smear steps per time step)
 *  A: accumulating matrix over all smear steps in a single time step
 *  S: staple (action dependent); recalculated before each smear step
 *  U: gauge links
 *  c1,c2: constants
 *    Calculating a single smear step is done by:
 *    A += c1*proj(S*U) -> update accumulation matrix
 *    U = exp(c2*A)*U   -> update gauge links
 */
void
stout_smear_step( Real c1, Real c2 )
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

      /* Update the accumulation matrix A += c1*proj(U*S) */
      mult_su3_na( U, &(s->staple[dir]), &tempS1 ); 
      anti_hermitian_traceless_proj( &tempS1, &tempA1 );
      scalar_mult_add_ah( Acur, &tempA1, c1, Acur );

      /* Update the links U = exp(c2*A)*U */
      scalar_mult_ah( Acur, c2, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, exp_order );
      mult_su3_nn( &tempS1, U, &tempS2 );
      su3mat_copy( &tempS2, U );
  }
}

/* Luscher's integration routine (essentially Runga-Kutta) */
/*  Outlined in 'arXiv:1006.4518 [hep-lat]'                */
void 
stout_step_rk()
{
  register int dir, i;
  register site *s;

  /* Clear the accumulation matrix */
  FORALLSITES(i, s)
    FORALLUPDIR(dir)
      clear_anti_hermitian(&(s->accumulate[dir]));

  /* Infinitesimal stout smearing */
  stout_smear_step( 17./36.*stepsize, -9./17. );
  stout_smear_step( -8./9.*stepsize, 1. );
  stout_smear_step( 3./4.*stepsize, -1. );
}
