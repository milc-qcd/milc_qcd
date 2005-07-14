/****************** msq_su3vec.c  (in su3.a) ****************************/
/* MIMD version 7 */
/*									*
* Real magsq_su3vec( su3_vector *a )					*
* return squared magnitude of an SU3 vector
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
Real magsq_su3vec( su3_vector *a ){
register Real sum;
register int i;
    for(i=0,sum=0.0;i<3;i++)sum += a->c[i].real*a->c[i].real
	+ a->c[i].imag*a->c[i].imag;
    return(sum);
}

#else
#ifdef NATIVEDOUBLE /* IBM RS6000 version */
Real magsq_su3vec(su3_vector *a){

  register double ar,ai,sum;

  ar=a->c[0].real; ai=a->c[0].imag;
  sum = ar*ar + ai*ai;

  ar=a->c[1].real; ai=a->c[1].imag;
  sum += ar*ar + ai*ai;

  ar=a->c[2].real; ai=a->c[2].imag;
  sum += ar*ar + ai*ai;

  return((Real)sum);
}
#else
Real magsq_su3vec( su3_vector *a ){
register Real temp,sum;
    sum=0.0;
    temp = a->c[0].real*a->c[0].real; sum += temp;
    temp = a->c[0].imag*a->c[0].imag; sum += temp;
    temp = a->c[1].real*a->c[1].real; sum += temp;
    temp = a->c[1].imag*a->c[1].imag; sum += temp;
    temp = a->c[2].real*a->c[2].real; sum += temp;
    temp = a->c[2].imag*a->c[2].imag; sum += temp;
    return(sum);
}
#endif  /* End of "#ifdef NATIVEDOUBLE" */
#endif /* end ifdef FAST */
