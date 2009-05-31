/******************  dumptensor4.c  (in su3.a) **************************
*									*
*  void dumptensor4( su3_matrix *mat )					*
*  print out a 3x3x3x3 tensor
*/
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dumptensor4( su3_tensor4 *a ) {
  int i, j, k, l;

  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      for( k=0; k<3; k++) {
        for( l=0; l<3; l++) {
          fprintf( stdout, " Re[%d][%d][%d][%d]=%18.10g    Im[%d][%d][%d][%d]=%18.10g\n",
                   i, j, k, l, (a->t4[i][j][k][l]).real, 
                   i, j, k, l, (a->t4[i][j][k][l]).imag );
        }
      }
    }
  }
}


