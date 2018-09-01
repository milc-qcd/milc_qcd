/***************** linearlsq.c ***************************************/

/* Linear least squares */

/* MIMD version 7 */

/* 06/27/16 C. DeTar */

#include "rcorr_includes.h"

double linearlsq(double *m, double *sdm, double *b, double *sdb,
		 double x[], double y[], double sd[], int n)
{
  double sd2inv;
  double sumx,sumx2,sum1,sumxy,sumy,sumy2,denom;
  double chisq;
  int df, i;

  *m = 0.;
  *b = 0.;
  *sdm = 0.;
  *sdb = 0.;
  
  if(n == 0)return 0.;

  if(n == 1){
    *b = y[0];
    return 0.;
  }
  /* Initialize sums */
  sumx = sumx2 = sum1 = sumxy = sumy = sumy2 = 0.;

  for(i = 0; i < n; i++){
    sd2inv = 1./(sd[i]*sd[i]);
    /* Prepare to multiply, which is more efficient than dividing */
    sum1 += sd2inv;
    sumx += x[i]*sd2inv;
    sumx2 += x[i]*x[i]*sd2inv;
    sumxy += x[i]*y[i]*sd2inv;
    sumy += y[i]*sd2inv;
    sumy2 += y[i]*y[i]*sd2inv;
  }

  /* Calculate best slope, intercept, errors, chisquare, and df */
  denom = sumx2*sum1 - sumx*sumx;
  *m = (sumxy*sum1 - sumx*sumy)/denom;
  *b = (sumx2*sumy - sumx*sumxy)/denom;
  *sdm = sqrt(sum1/denom);
  *sdb = sqrt(sumx2/denom);
  chisq = sumy2 - *m*sumxy - *b*sumy;
  df = n - 2;

  return chisq;
}
