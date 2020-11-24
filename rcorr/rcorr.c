/***************** rcorr.c *****************************************/

/* Measure the correlator from the given densities                   */

/* MIMD version 7 */

/* 04/03/15 C. DeTar (after topo_rcorr.c) */

#include "rcorr_includes.h"

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
dot_corr( complex *dest, complex *src, int count )
{
  int i, j;

  FORALLFIELDSITES(i){
    dest[i].real = 0.;
    /* Sum only the spatial components -- hence count-1 */
    for(j = 0; j < count-1; j++)
      dest[i].real += 
	(src[i*count+j].real*src[i*count+j].real + src[i*count+j].imag*src[i*count+j].imag)/volume; 
    
    dest[i].imag = 0.;
  }
}

static void 
sum_field_c2r( Real *dest, complex *src )
{
  int i;
  FORALLFIELDSITES(i){
    dest[i] += src[i].real;
  }
}

static void 
sum_sq_field_c2r( Real *dest, complex *src )
{
  int i;
  FORALLFIELDSITES(i){
    dest[i] += src[i].real*src[i].real;
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

static void
mean_var_r_field(Real *q, Real *q2, int n)
{
  int i;
  Real d;

  if(n > 1){
    FORALLFIELDSITES(i){
      q2[i] = q2[i]/n;
      q[i] = q[i]/n;
      d = q2[i] - q[i]*q[i];
      if(d < 0.)d = 0.;
      q2[i] = d/(n-1);
    }
  } else {
    clear_r_field(q2);
  }
}

/******************************************************************************/
/* Input fields are have NMU values per site (one for each current component)
   There is one such field per random number group */
/* Output field has one real value per site, one such field for each
   blocking size */

void
rcorr(Real *qblock[], Real *q2block[], 
      complex *qin_sloppy[], int nrand_sloppy, 
      complex *qin_diff[], int nrand_diff,
      int nblock, int block_size[]){
  complex *qtmp;
  int jrand;

  /* Average qin_diff, the differnece between precise and sloppy. */
  /* Result in qcorr */
  complex *qcorr = create_c_array_field(NMU);

  /* Add up the results in qin_diff */
  
  for(jrand = 0; jrand < nrand_diff; jrand++)
    sum_c_array_field(qcorr, qin_diff[jrand], NMU);

  /* Average the accumulated results */
  if(nrand_diff > 0)
    mulreal_c_field(qcorr, 1./((double) nrand_diff), NMU);

  /* Correct the sloppy results for each random block by adding the
     difference between precise and sloppy*/
  for(jrand = 0; jrand < nrand_sloppy; jrand++)
    sum_c_array_field( qin_sloppy[jrand], qcorr, NMU);

  destroy_c_array_field(qcorr, NMU);
  
  /* Create blocked averages of the corrected result, 
     and compute the correlation within each block */
  complex *out = create_c_field();
  qtmp = create_c_array_field(NMU);
  for(int ib = 0; ib < nblock; ib++){

    clear_r_field(qblock[ib]);
    clear_r_field(q2block[ib]);
    int bs = block_size[ib];   /* Size of one block */
    int nsamp = 0;  /* Number of blocks of the given size */

    for(jrand = 0; jrand < nrand_sloppy; jrand += bs){

      /* Compute average of current density for this block */
      clear_c_array_field(qtmp, NMU);
      for(int kb = 0; kb < bs; kb++)
	sum_c_array_field(qtmp, qin_sloppy[jrand+kb], NMU);
      mulreal_c_field(qtmp, 1./((double) bs), NMU);
  
      /* The forward FT is done separately for each current component */
      restrict_fourier_field((complex *)qtmp, NMU*sizeof(complex), FORWARDS);
      
      /* Square and sum over components to get the correlator */
      dot_corr(out, qtmp, NMU);

      /* Consistency check */
      double qtot = 0.;
      int i;
      FORALLFIELDSITES(i){
	qtot += out[i].real;
      }
      g_doublesum(&qtot);
      node0_printf("qtot[%d][%d] = %g\n",ib, jrand, qtot);
      
      /* The backward FT */
      restrict_fourier_field(out, sizeof(complex), BACKWARDS);

      /* Accumulate the result for this block size */
      sum_field_c2r(qblock[ib], out);

      /* Accumulate the squares of the results */
      sum_sq_field_c2r(q2block[ib], out);

      nsamp++;
    }

    /* Compute the mean and the variance of the mean */

    mean_var_r_field(qblock[ib], q2block[ib], nsamp);
  }

  destroy_c_field(qtmp);
  destroy_c_field(out);

} /* rcorr.c */
