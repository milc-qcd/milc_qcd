/************* msw_gamma_r.c (in su3.a) **************************/

#include "generic_form_includes.h"

/* 
  Multiply a "Wilson matrix" (spin_wilson_vector) by a gamma matrix
  acting on the column index
  (This is the second index, or equivalently, multiplication on the right)
  usage:  mb_gamma_r( src, dest, dir)
	spin_wilson_vector *src,*dest;
	int dir;    dir = XUP, YUP, ZUP, TUP or GAMMAFIVE
*/

void mult_sw_by_gamma_r(spin_wilson_vector * src,spin_wilson_vector * dest, int dir)
{
  register int c2,s1,s2,s;	/* column indices, color and spin */

  if(gamma_initialized==0)make_gammas(gamma_mat);

  /* For compatibility */
  if(dir == GAMMAFIVE)dir = G5;
  if(dir >= MAXGAMMA)
    {
      printf("mult_sw_by_gamma_r: Illegal gamma index %d\n",dir);
      exit(1);
    }

  for(s=0;s<4;s++){
    s2 = gamma_mat[dir].row[s].column;
    switch (gamma_mat[dir].row[s].phase){
    case 0:
      for(s1=0;s1<4;s1++)for(c2=0;c2<3;c2++){
	dest->d[s1].d[s2].c[c2] = src->d[s1].d[s].c[c2];}
      break;
    case 1:
      for(s1=0;s1<4;s1++)for(c2=0;c2<3;c2++){
	TIMESPLUSI( src->d[s1].d[s].c[c2], dest->d[s1].d[s2].c[c2] );}
      break;
    case 2:
      for(s1=0;s1<4;s1++)for(c2=0;c2<3;c2++){
	TIMESMINUSONE( src->d[s1].d[s].c[c2], dest->d[s1].d[s2].c[c2] );}
      break;
    case 3:
      for(s1=0;s1<4;s1++)for(c2=0;c2<3;c2++){
	TIMESMINUSI( src->d[s1].d[s].c[c2], dest->d[s1].d[s2].c[c2] );}
    }
  }
}
