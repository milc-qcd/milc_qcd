/**************  m_su2_mat_vec_n.c (in su3.a) **********************
*									*
*  su2 matrix times vector                              		*
*/

#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mult_su2_mat_vec_elem_n(su2_matrix *u,complex *x0,complex *x1)
{
  /* Multiplies the complex column spinor (x0, x1) by the SU(2) matrix u */
  /* and puts the result in (x0,x1).  */
  /* Thus x <- u * x          */
  /* C. DeTar 3 Oct 1990 */
  
  complex z0, z1, t0, t1;

  t0 = *x0; t1 = *x1;

  CMUL(u->e[0][0], t0, z0);
  CMUL(u->e[0][1], t1, z1);
  CADD(z0, z1, *x0);
  CMUL(u->e[1][0], t0, z0);
  CMUL(u->e[1][1], t1, z1);
  CADD(z0, z1, *x1);

} /* m_su2_mat_vec_elem_n.c */

