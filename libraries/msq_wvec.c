 /********************  msq_wvec.c  (in su3.a) ********************
*
*Real msq_wvec(wilson_vector *vec)
*  squared magnitude of a Wilson vector
* 
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
Real magsq_wvec( wilson_vector *vec ){
  register int i;
  register Real sum;
  sum=0.0;
  for(i=0;i<4;i++)sum += magsq_su3vec( &(vec->d[i]) );
  return(sum);

#else /* Fast version */
Real magsq_wvec( wilson_vector *vec ){

#ifdef NATIVEDOUBLE
  register double ar,ai,sum;
#else
  register Real ar,ai,sum;
#endif

  ar=vec->d[0].c[0].real; ai=vec->d[0].c[0].imag;
  sum = ar*ar + ai*ai;
  ar=vec->d[0].c[1].real; ai=vec->d[0].c[1].imag;
  sum += ar*ar + ai*ai;
  ar=vec->d[0].c[2].real; ai=vec->d[0].c[2].imag;
  sum += ar*ar + ai*ai;

  ar=vec->d[1].c[0].real; ai=vec->d[1].c[0].imag;
  sum += ar*ar + ai*ai;
  ar=vec->d[1].c[1].real; ai=vec->d[1].c[1].imag;
  sum += ar*ar + ai*ai;
  ar=vec->d[1].c[2].real; ai=vec->d[1].c[2].imag;
  sum += ar*ar + ai*ai;

  ar=vec->d[2].c[0].real; ai=vec->d[2].c[0].imag;
  sum += ar*ar + ai*ai;
  ar=vec->d[2].c[1].real; ai=vec->d[2].c[1].imag;
  sum += ar*ar + ai*ai;
  ar=vec->d[2].c[2].real; ai=vec->d[2].c[2].imag;
  sum += ar*ar + ai*ai;

  ar=vec->d[3].c[0].real; ai=vec->d[3].c[0].imag;
  sum += ar*ar + ai*ai;
  ar=vec->d[3].c[1].real; ai=vec->d[3].c[1].imag;
  sum += ar*ar + ai*ai;
  ar=vec->d[3].c[2].real; ai=vec->d[3].c[2].imag;
  sum += ar*ar + ai*ai;

  return((Real)sum);
#endif
}
