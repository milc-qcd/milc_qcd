/********************  clear_wvec.c  (in su3.a) ********************
*
*void clear_wilson_vector( wilson_vector *dest )
*  clear a Wilson vector
* dest  <-  zero_vector
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clear_wvec( wilson_vector *dest ){
register int i,j;
    for(i=0;i<4;i++)for(j=0;j<3;j++){
	dest->d[i].c[j].real = dest->d[i].c[j].imag = 0.0;
    }
}
