/******************  dumpmat.c  (in su3.a) ******************************
*									*
*  void dumpmat( su3_matrix *mat )					*
*  print out a 3x3 complex matrix					*
*/
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dumpmat( su3_matrix *m ){
int i,j;
    for(i=0;i<3;i++){
	for(j=0;j<3;j++)printf("(%.2e,%.2e)\t",
	    m->e[i][j].real,m->e[i][j].imag);
	printf("\n");
    }
    printf("\n");
}
