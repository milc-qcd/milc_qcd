/************* msw_gamma_l.c (in su3.a) **************************/
/*
  Multiply a "spin Wilson vector" by a gamma matrix
  acting on the row index
  (This is the first index, or equivalently, multiplication on the left)
  usage:  mult_sw_by_gamma_l( src, dest, dir)
	spin_wilson_vector *src,*dest;
	int dir;    dir = any of the gamma matrix types in gammatypes.h


  Originally from MILC su3.a

  Modifications:
  UMH - modified for spin Wilson vector
  4/29/97 C. DeTar generalized to any gamma matrix.
*/

#include "generic_form_includes.h"

void mult_sw_by_gamma_l(spin_wilson_vector * src,spin_wilson_vector * dest, int dir)
{
  register int c2,s1,s2,s;	/* column indices, color and spin */

  if(gamma_initialized==0)make_gammas(gamma_mat);

  /* For compatibility */
  if(dir == GAMMAFIVE)dir = G5;

  if(dir >= MAXGAMMA)
    {
      printf("mult_sw_by_gamma_l: Illegal gamma index %d\n",dir);
      exit(1);
    }

  for(s1=0;s1<4;s1++){
    s = gamma_mat[dir].row[s1].column;
    switch (gamma_mat[dir].row[s1].phase){
    case 0:
      for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	dest->d[s1].d[s2].c[c2] = src->d[s].d[s2].c[c2];}
      break;
    case 1:
      for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	TIMESPLUSI( src->d[s].d[s2].c[c2], dest->d[s1].d[s2].c[c2] );}
      break;
    case 2:
      for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	TIMESMINUSONE( src->d[s].d[s2].c[c2], dest->d[s1].d[s2].c[c2] );}
      break;
    case 3:
      for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	TIMESMINUSI( src->d[s].d[s2].c[c2], dest->d[s1].d[s2].c[c2] );}
    }
  }
}
