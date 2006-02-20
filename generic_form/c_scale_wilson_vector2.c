/**************************** c_scale_wilson_vector2.c ****************************/
/* MIMD version 7 */
/*
 *  Scale a su3_vectore by a complex vector
 *
 *  m ----> m * scale
 */


#include "../include/complex.h"
#include "../include/su3.h"


void c_scale_wilson_vector2(wilson_vector *m , complex *scale)
{
  int i , j  ;
  complex z ; 

  for( i = 0 ; i < 3 ; ++i)
    for(j=0 ; j < 4 ; ++j)
    {
      CMUL(m->d[j].c[i] ,(*scale), z);
      m->d[j].c[i]  =  z ;
    }

}
