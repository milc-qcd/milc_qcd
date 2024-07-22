/************************* flow_helpers.c **********************/
/* Utilities associated with integrating the Wilson Flow 
   or otherwise generic utilities used across multiple modules */

#include "wilson_flow_includes.h"
#include <string.h>


void print_observables( char *TAG, double flowtime, 
                        double *Et_WS, double *Es_WS, 
                        double *Et_C, double *Es_C, 
                        double *charge ) {

  #if (MILC_PRECISION==1)
    #if (REPORT == OLDREPORT) 
      node0_printf("%s %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n", 
        TAG, flowtime, Et_C[0], Es_C[0], Et_WS[0], Es_WS[0], Et_WS[1], Es_WS[1], 
        charge[0]);
    #else // REPORT == NEWREPORT || REPORT == NO_REPORT
      node0_printf("%s %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n", 
        TAG, flowtime, Et_C[0], Es_C[0], Et_C[1], Es_C[1], Et_WS[0], Es_WS[0], Et_WS[1], Es_WS[1], 
        charge[0], charge[1], charge[2], charge[3]);
    #endif
  #else
    #if (REPORT == OLDREPORT) 
      node0_printf("%s %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n", 
        TAG, flowtime, Et_C[0], Es_C[0], Et_WS[0], Es_WS[0], Et_WS[1], Es_WS[1], 
        charge[0]);
    #else // REPORT == NEWREPORT || REPORT == NO_REPORT
      node0_printf("%s %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n", 
        TAG, flowtime, Et_C[0], Es_C[0], Et_C[1], Es_C[1], Et_WS[0], Es_WS[0], Et_WS[1], Es_WS[1], 
        charge[0], charge[1], charge[2], charge[3]);
    #endif
  #endif
}

void clear_field ( su3_matrix **this, size_t ncomp ) {

  if ( this == NULL ) {
    node0_printf("Tried to clear unallocated field\n"); 
    terminate(1);
  }
  if ( this[0] == NULL ) {
    node0_printf("Tried to clear unallocated field\n"); 
    terminate(1);
  }
  memset( this[0], '\0', ncomp * sizeof(su3_matrix) * sites_on_node );
}

void destroy_field ( su3_matrix ***this ) {

  if ( *this == NULL ) {
    node0_printf("Tried to clear unallocated field\n"); 
    terminate(1);
  }
  if ( **this == NULL ) {
    node0_printf("Tried to free unallocated field\n"); 
    terminate(1);
  }

  free( **this );
  free( *this );
  *this = NULL;
}

su3_matrix ** new_field( size_t ncomp ) {

  register int i;
  register site *s;
  
  su3_matrix **this = (su3_matrix **)malloc( ncomp * sizeof(su3_matrix*) );
  if(this == NULL) {
    printf( "new_field: can't malloc this\n" );
    fflush(stdout); terminate(1);
  }

  this[0] = (su3_matrix *)malloc( ncomp * sizeof(su3_matrix) * sites_on_node );
  if(this[0] == NULL) {
    printf( "new_field: can't malloc this[0]\n" );
    fflush(stdout); terminate(1);
  }
  memset( this[0], '\0', ncomp * sizeof(su3_matrix) * sites_on_node );
  for ( i = 1; i < ncomp; i++ ) 
    this[i] = this[0] + i * sites_on_node;

  return ( this );
}

su3_matrix ** new_links_from_site( int region_flag, size_t nlink ) {

  register int i, dir;
  register site *s;
  if( nlink > 4 ) {
    node0_printf( "new_links_from_site: can't create than four links from site; instead %ld\n", nlink );
    fflush(stdout); terminate(1);
  }
  su3_matrix **this = new_field( nlink );

  for (dir = XUP; dir < nlink; dir++ )
  FORALLSITES(i, s) 
  IF_BLOCKED(s, block_stride)       
  IF_REGION(s, region_flag) 
    su3mat_copy( &(s->link[dir]), &(this[dir][i]) );

  return ( this );
}


/* Clear an anti-hermition field */
void clear_anti_hermitian_field ( anti_hermitmat **this, size_t ncomp ) {

  if ( this == NULL ) {
    node0_printf("Tried to clear unallocated anti_hermitian field\n"); 
    terminate(1);
  }
  if ( this[0] == NULL ) {
    node0_printf("Tried to clear unallocated anti_hermitian field\n"); 
    terminate(1);
  }
  memset( this[0], '\0', ncomp * sizeof(anti_hermitmat) * sites_on_node );
}

/* Destroy an anti-hermition field */
void destroy_anti_hermitian_field ( anti_hermitmat ***this ) {

  if ( *this == NULL ) {
    node0_printf("Tried to clear unallocated anti_hermitian field\n"); 
    terminate(1);
  }
  if ( **this == NULL ) {
    node0_printf("Tried to free unallocated anti_hermitian field\n"); 
    terminate(1);
  }

  free( **this );
  free( *this );
  *this = NULL;
}

/* Create a new, empty anti-hermitian field */
anti_hermitmat ** new_anti_hermitian_field( size_t ncomp ) {

  register int i;
  register site *s;
  anti_hermitmat **this = (anti_hermitmat **)malloc( ncomp * sizeof(anti_hermitmat*) );
  if(this == NULL) {
    printf( "new_anti_hermitian_field: can't malloc this\n" );
    fflush(stdout); terminate(1);
  }

  this[0] = (anti_hermitmat *)malloc( ncomp * sizeof(anti_hermitmat) * sites_on_node );
  if(this[0] == NULL) {
    printf( "new_anti_hermitian_field: can't malloc this[0]\n" );
    fflush(stdout); terminate(1);
  }
  memset( this[0], '\0', ncomp * sizeof(anti_hermitmat) * sites_on_node );
  for ( i = 1; i < ncomp; i++ ) 
    this[i] = this[0] + i * sites_on_node;

  return ( this );
}

/* Destroy an anti-hermition twodim_field */
void destroy_anti_hermitian_twodim_field ( anti_hermitmat ****this ) {

  if ( *this == NULL ) {
    node0_printf("Tried to clear unallocated anti_hermitian twodim_field\n"); 
    terminate(1);
  }
  if ( **this == NULL ) {
    node0_printf("Tried to free unallocated anti_hermitian twodim_field\n"); 
    terminate(1);
  }
  if ( ***this == NULL ) {
    node0_printf("Tried to free unallocated anti_hermitian twodim_field\n"); 
    terminate(1);
  }

  free( ***this );
  free( **this );
  free( *this );
  *this = NULL;
}

/* Create a new, empty anti-hermitian twodim_field */
anti_hermitmat *** new_anti_hermitian_twodim_field( size_t dim1, size_t dim2 ) {

  register int i1,i2;
  register site *s;
  anti_hermitmat ***this = (anti_hermitmat ***)malloc( dim1 * sizeof(anti_hermitmat**) );
  if(this == NULL) {
    printf( "new_anti_hermitian_field: can't malloc this\n" );
    fflush(stdout); terminate(1);
  }
  this[0] = (anti_hermitmat **)malloc( dim1 * sizeof(anti_hermitmat*) * dim2 );
  if(this[0] == NULL) {
    printf( "new_anti_hermitian_field: can't malloc this[0]\n" );
    fflush(stdout); terminate(1);
  }
  this[0][0] = (anti_hermitmat *)malloc( dim1 * dim2 * sizeof(anti_hermitmat) * sites_on_node );
  if(this[0][0] == NULL) {
    printf( "new_anti_hermitian_field: can't malloc this[0]\n" );
    fflush(stdout); terminate(1);
  }
  memset( this[0][0], '\0', dim1 * dim2 * sizeof(anti_hermitmat) * sites_on_node );
  for ( i1 = 0; i1 < dim1; i1++ ) {
    this[i1] = this[0] + i1 * dim2;
    for ( i2 = 0; i2 < dim2; i2++ ) 
      this[i1][i2] = this[0][0] + ( i2 + i1 * dim2 ) * sites_on_node;
  }

  return ( this );
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