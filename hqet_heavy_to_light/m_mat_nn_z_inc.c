/*******************  m_mat_nn_z_inc.c  (in su3.a) ****************************
*									*
* void mult_su3_nn( su3_matrix *a,*b,*c )				*
* matrix multiply, no adjoints 						*
* C  = C + z*A*B							*
*
* where  A,B and C are su3 matrices; and z is a complex number
*
* I have removed the MACRO option for the handcoded
* version
*/
#include "../include/complex.h"
#include "../include/su3.h"

void mult_su3_nn_z_inc(complex z,su3_matrix *a,su3_matrix *b,su3_matrix *c)
{
  register int i,j,k;
  register complex x,y,yy;

    for(i=0;i<3;i++)for(j=0;j<3;j++)
    {
      x.real=x.imag=0.0;
      for(k=0;k<3;k++){
	CMUL( a->e[i][k] , b->e[k][j] , y );
	CSUM( x , y );
      }
      CMUL(x ,z, yy);
      CSUM(c->e[i][j] , yy);
    }

/*    CMUL(a,b,c)       c = a * b					      */
/*    CSUM(a,b)         a += b						      */



}

