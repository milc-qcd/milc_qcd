/**************************** 3submat.c ********************************/
/* MIMD version 6 */
/*									*
* subtract su3 matrices:  D = D - A - B - C 
*
*  This routine has not been written for speed.
*
*/

#include "../include/complex.h"
#include "../include/su3.h"

/*** a -= B ****/
#define CDEC(a,b) { (a).real -= (b).real; (a).imag -= (b).imag; }



/* subtract su3 matrices */
void sub3_su3_matrix(su3_matrix *a, su3_matrix*b, su3_matrix*c, su3_matrix *d) 
{
  register int i,j;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      CDEC( d->e[i][j], a->e[i][j] );
      CDEC( d->e[i][j], b->e[i][j] );
      CDEC( d->e[i][j], c->e[i][j] );
    }

}
