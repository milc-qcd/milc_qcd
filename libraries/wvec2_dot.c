/****************** wvec2_dot.c  (in su3.a) ***************************/
/* MIMD version 7 */
/*									*
* complex wvec2_dot( wilson_vector *a, wilson_vector *b )		*
* return dot product of two wilson_vectors = a-dagger times b		*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex wvec2_dot( wilson_vector *a, wilson_vector *b ){
    complex temp;
    wilson_vector c;
    register int i,j;

    temp.real = wvec_rdot(a,b);    	

    for(i=0;i<4;i++){
    for(j=0;j<3;j++){
        c.d[i].c[j].real = -(a->d[i].c[j].imag);
        c.d[i].c[j].imag = a->d[i].c[j].real;
    }
    }

    temp.imag = wvec_rdot(&c,b);

    return(temp);
}
