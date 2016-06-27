/***************** rcorr.c *****************************************/

/* Measure the correlator from the given densities                   */

/* MIMD version 7 */

/* 04/03/15 C. DeTar (after topo_rcorr.c) */

#include "rcorr_includes.h"

static void 
add_real_array_field( complex *dest, complex *src, int count ){
  int i, j;

  FORALLFIELDSITES(i){
    for(j = 0; j < count; j++){
      dest[count*i+j].real += src[count*i+j].real; /* Here we use only the real part */
    }
  }
}

static void 
mulreal_c_field( complex *c, double x, int count )
{
  int i, j;
  FORALLFIELDSITES(i){
    for(j = 0; j < count; j++){
      CMULREAL(c[count*i+j], x, c[count*i+j]);
    }
  }
}

static void 
sum_c_array_field( complex *dest, complex *src, int count ){
  int i, j;

  FORALLFIELDSITES(i){
    for(j = 0; j < count; j++){
      CSUM(dest[count*i+j], src[count*i+j]);
    }
  }
}

static void 
add_corr( complex *dest, complex *src, int count )
{
  int i, j;

  FORALLFIELDSITES(i){
    for(j = 0; j < count; j++)
      dest[i].real += 
	(src[i*count+j].real*src[i*count+j].real + src[i*count+j].imag*src[i*count+j].imag)/volume; 
    
    dest[i].imag = 0;
  }
}

static void 
sub_corr( complex *dest, complex *src, int count )
{
  int i, j;

  FORALLFIELDSITES(i){
    for(j = 0; j < count; j++)
      dest[i].real -= (src[i*count+j].real*src[i*count+j].real + 
		       src[i*count+j].imag*src[i*count+j].imag)/volume; 
    dest[i].imag = 0;
  }
}

static void 
copy_mul_c2r( Real *dest, complex *src, double x )
{
  int i;
  FORALLFIELDSITES(i){
    dest[i] = src[i].real * x;
  }
}

/******************************************************************************/
/* Input fields are have NMU values per site (one for each current component) */

Real *
rcorr(complex *qin_sloppy, int nrand_sloppy, complex *qin_diff, int nrand_diff){
  //  complex *qtmp;
  complex *out2;
  Real *q;
  int jrand;
  int i;

  /* qtmp contains the real data for all the current components for all sites */
  //  qtmp = create_c_array_field(NMU);
  /* out2 contains the real dot products of the currents for all sites */
  out2 = create_c_field();
//  if(qtmp == NULL || out2 == NULL){
//    fprintf(stderr,"No room for qtmp/out/out2\n");
//    terminate(1);
//  }
  
//  /* Compute sum of squares */
//  for (jrand = 0; jrand < nrand; jrand++){
//    clear_c_array_field(qtmp, NMU);
//    add_real_array_field(qtmp, qin[jrand], NMU);
//    
//    /* The forward transform is done separately for each current component */
//    restrict_fourier_field((complex *)qtmp, NMU*sizeof(complex), FORWARDS);
//    
//    /* Take the dot product of the current with itself and accumulate */
//    sub_corr(out2, qtmp, NMU);
//  }
  
//  /* Compute the square of the sum */
//  /* Accumulate sums in qtmp */
//  clear_c_array_field(qtmp, NMU);
//  for (jrand = 0; jrand < nrand; jrand++){
//    add_real_array_field(qtmp, qin[jrand], NMU);
//  }

  /* Average the accumulated results */
  mulreal_c_field(qin_sloppy, 1./((double) nrand_sloppy), NMU);
  mulreal_c_field(qin_diff, 1./((double) nrand_diff), NMU);

  /* Correct the sloppy result by adding the difference between precise and sloppy*/
  sum_c_array_field( qin_sloppy, qin_diff, NMU);
  
  
  /* The forward transform is done separately for each current component */
  restrict_fourier_field((complex *)qin_sloppy, NMU*sizeof(complex), FORWARDS);

  /* Square and sum over components to get the correlator */
  add_corr(out2, qin_sloppy, NMU);
  
  /* For a consistency check. */
  double qtot = 0.;
  FORALLFIELDSITES(i){
    qtot += out2[i].real;
  }
  g_doublesum(&qtot);
  //  qtot /= (nrand*(nrand-1));
  node0_printf("qtot = %g\n",qtot);

  /* Backward FFT */
  restrict_fourier_field(out2, sizeof(complex), BACKWARDS);
  
//  if(q == NULL ){
//    fprintf(stderr,"No room for q\n");
//    terminate(1);
//  }
  
  /* Normalize the result and copy the real part */
  //  copy_mul_c2r(q, out2, 1./(nrand*(nrand-1)) );
  q = create_r_field();
  copy_mul_c2r(q, out2, 1.);
  destroy_c_field(out2);

  return q;
} /* rcorr.c */


