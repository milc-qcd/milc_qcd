/******************  m_mat_wvec.c  (in su3.a) ********************
*									*
*void mult_mat_wilson_vec(su3_matrix *mat, wilson_vector *src,*dest)	*
*  multiply a Wilson vector by a matrix					*
* dest  <-  mat*src							*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mult_mat_wilson_vec(  su3_matrix *mat, wilson_vector *src,
	wilson_vector *dest ){
    register int i;
    for(i=0;i<4;i++)mult_su3_mat_vec(mat, &(src->d[i]), &(dest->d[i]) );
}
