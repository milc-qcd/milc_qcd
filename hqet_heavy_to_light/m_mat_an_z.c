/******************  m_mat_an_z.c  (in su3.a) *****************************
*									*
* void mult_su3_an_z( su3_matrix *a,*b,*c )				*
* matrix multiply, first matrix is adjoint 				*
* C <-  A_adjoint*B*z							*
*
* whers A,B, and C are su3 matrices.
* z is a scalar
*
* I have taken out the hand coded part of this subroutine.
*
*/
#include "../include/complex.h"
#include "../include/su3.h"


void mult_su3_an_z(complex z,su3_matrix *a,su3_matrix *b,su3_matrix *c )
{

  register int i,j,k;
  register complex x,y;
    for(i=0;i<3;i++)for(j=0;j<3;j++)
    {
      x.real=x.imag=0.0;
      for(k=0;k<3;k++)
      {
	CMULJ_( a->e[k][i] , b->e[k][j], y );
	CSUM( x , y );
	}
      CMUL(x , z , c->e[i][j]);
    }

/*    CMUL(a,b,c)       c = a * b					      */
/*    CMULJ_(a,b,c)     c = conjg(a) * b				      */
/*    CSUM(a,b)         a += b						      */

}








