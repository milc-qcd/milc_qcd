/******************  complextr.c  (in su3.a) ****************************
*									*
* complex complextrace_su3( su3_matrix *a,*b)				*
* return Tr( A_adjoint*B )   						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex complextrace_su3( su3_matrix *a, su3_matrix *b ) {
register int i,j;
register Real sumr, sumi;
complex sum;
    for(sumr=0.0,sumi=0.0,i=0;i<3;i++)for(j=0;j<3;j++){
     sumr+= a->e[i][j].real*b->e[i][j].real + a->e[i][j].imag*b->e[i][j].imag;
     sumi+= a->e[i][j].real*b->e[i][j].imag - a->e[i][j].imag*b->e[i][j].real;
    }
    sum.real= sumr; sum.imag=sumi; 
    return(sum);
}
