/*******************  clearvec.c  (in su3.a) *****************************
*									*
*  void clearvec( su3_vector *vec )					*
*  clear a 3 element complex vector					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clearvec( su3_vector *v ){
     v->c[0].real = v->c[0].imag = 0.0;
     v->c[1].real = v->c[1].imag = 0.0;
     v->c[2].real = v->c[2].imag = 0.0;
}
