/**************************** cs_wv.c ****************************/
/* MIMD version 6 */
/*
*  $Header: /lqcdproj/detar/cvsroot/milc_qcd/wilson_static/cs_wv.c,v 1.1 2005/02/23 00:06:10 detar Exp $
*
*  Overwrite the a wilson vector with the wilso vector
*  times a complex scalar
*
*  m --> phase * m
*
*/
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_wilsonvec(wilson_vector *m, complex *phase)
{
  register int ic, ispin;
  complex z ;

  for(ic=0;ic<3;ic++)
    for(ispin=0;ispin<4; ispin++)
    {
      z = cmul(&m->d[ispin].c[ic],phase) ;
      m->d[ispin].c[ic] =  z;
    }

}


