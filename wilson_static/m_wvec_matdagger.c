/**************** m_wvec_matdagger.c *************************/
/* MIMD version 7 */
/*
 *  out_vec = in_vec * (m)^{dagger}
 *
 *  where out_vec, in_vec are wilson su3 vectors
 *  and m is an su3 matrix
 */

#include "../include/complex.h"
#include "../include/su3.h"

void mult_wilson_vec_matdag(wilson_vector *out_vec,
			    wilson_vector *in_vec, su3_matrix *m)
{
  int ispin,i,k ;
  complex z ;

  for(ispin = 0 ; ispin < 4 ; ++ ispin)
    for(i=0; i < 3 ;++i)
    {
      out_vec->d[ispin].c[i].real = 0.0  ;
      out_vec->d[ispin].c[i].imag = 0.0  ;

      for(k=0 ; k < 3 ; ++k)
      {

	CMUL_J(in_vec->d[ispin].c[k],m->e[i][k] , z);

	out_vec->d[ispin].c[i].real += z.real  ;
	out_vec->d[ispin].c[i].imag += z.imag  ;

      }

    }
		



}
