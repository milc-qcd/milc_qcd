/**************  m_su2_mat_vec_a.c (in su3.a) **********************
*									*
*  adjoint su2 matrix times vector                             		*
*/

#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mult_su2_mat_vec_elem_a(su2_matrix *u,complex *x0,complex *x1)
{
  /* Multiplies the complex row spinor (x0, x1) by the adjoint of the */
  /* SU(2) matrix u and puts the result in (x0,x1).  */
  /* Thus x <-  x * u-adj       */
  /* C. DeTar 3 Oct 1990 */
  
  complex z0, z1, t0, t1;

  t0 = *x0; t1 = *x1;

  CMUL_J(t0, u->e[0][0], z0);
  CMUL_J(t1, u->e[0][1], z1);
  CADD(z0, z1, *x0);
  CMUL_J(t0, u->e[1][0], z0);
  CMUL_J(t1, u->e[1][1], z1);
  CADD(z0, z1, *x1);

} /* m_su2_mat_vec_a.c */

