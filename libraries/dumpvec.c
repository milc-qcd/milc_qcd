/*******************  dumpvec.c  (in su3.a) *****************************
*									*
*  void dumpvec( su3_vector *vec )					*
*  print out a 3 element complex vector					*
*/
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dumpvec( su3_vector *v ){
int j;
    for(j=0;j<3;j++)printf("(%.2e,%.2e)\t",
	v->c[j].real,v->c[j].imag);
    printf("\n");
}
