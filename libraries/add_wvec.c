/********************  add_wvec.c  (in su3.a) ********************
*
*void add_wilson_vector(wilson_vector *src1,*src2,*dest)
*  add two Wilson vectors
* dest  <-  src1 + src2
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void add_wilson_vector( wilson_vector *src1, wilson_vector *src2,
       wilson_vector *dest ){
   register int i;
   for(i=0;i<4;i++)add_su3_vector( &(src1->d[i]), &(src2->d[i]), &(dest->d[i]));
}
