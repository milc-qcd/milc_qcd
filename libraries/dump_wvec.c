/****************  dump_wvec.c  (in su3.a) ***********************
*									*
*  void dump_wvec( wilson_vector *v )					*
*  Print out a Wilson vector 						*
*/
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dump_wvec( wilson_vector *v ){
register int i,j;
    for(i=0;i<4;i++){
        for(j=0;j<3;j++)printf("(%.2e,%.2e)\t",
	    v->d[i].c[j].real,v->d[i].c[j].imag);
        printf("\n");
    }
    printf("\n");
}
