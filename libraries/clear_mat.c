/********************  clear_mat.c  (in su3.a) ********************
*
*void clear_su3mat( su3_matrix *dest )
*  clear an SU3 matrix
* dest  <-  zero_matrix
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clear_su3mat( su3_matrix *dest ){
register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	dest->e[i][j].real = dest->e[i][j].imag = 0.0;
    }
}
